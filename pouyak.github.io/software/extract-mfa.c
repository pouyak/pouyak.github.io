/* Written by Pouya Kheradpour (pouyak@gmail.com)
 * Copyright (c) 2005-2015 Massachusetts Institute of Technology
 * Released under the MIT license
 * https://pouya.kheradpour.com/software/LICENSE.mit */
/* version 20160416 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pk.h"
#include "rfile.h"

#define MAXSEQNAMELEN		300

/* indexing for alignment */
#define SEQidx(sp, loc)		((sp)*alignment->alloc+(loc))
#define SEQNAMEidx(sp)		((sp)*(MAXSEQNAMELEN+1))

struct alignment
{
	int first_file;
	int last_file;
	char* seqs;
	char* names;
	int length;
	int start;
	int end;
	int alloc;
	int num_sp;
};

/* reads the alignment in filename updating alignment, returns 1 if error */
int read_alignment(char* filename, unsigned long long int sp_mask, struct alignment* alignment, int offset);

/* loads the alignment corresponding to chr, start, end. returns 0 if error or not found */
int load_alignment(struct rfile* files, int num_files, long long int sp_mask, char *chr, int start, int end, struct alignment* alignment, int search_start);

/* number of non-zero bits in a integer */
unsigned char bitstr_num_nonzero (long long int str);

int main (int argc, char** argv)
{
	struct rfile *files;
	int num_files;
	struct alignment* alignment = (struct alignment*) malloc(sizeof(struct alignment));
	char *line = NULL;
	size_t line_alloc = 0;
	int line_len;
	long long int sp_mask = -1;
	int buffer_alloc = 100000;
	char *buffer = (char*) malloc ((buffer_alloc)*sizeof(char));
	unsigned char *print_pos = (unsigned char*) malloc (buffer_alloc * sizeof(unsigned char));
	unsigned char skip_first = 0;

	int i;

	if (argc < 3)
	{
		printf("USAGE: %s [OPTIONS] <NumOrg> <AlignFile> < <Locs>\n", argv[0]);
		printf("    -k         Species mask to use. NOTE: this /includes/ the target species (-1 for all) [default: -1]. \n");
		printf("    NumOrg:    Number of organisms in each alignment file (optional when -k is not -1)\n");
		printf("    AlignFile: File containing list of all the alignment files (chr start stop filename).\n");
		printf("    Locs:      File containing locations (name chr start end [strand])\n");
		return 1;
	}

	for (i=1; i<(argc-2); i++)
		if (STREQ(argv[i], ""))
			continue; 
		else if (STREQ(argv[i], "-k"))
			sp_mask = atoll(argv[++i]); 
		else
			ASSERTE(0, "Unrecognized command line argument!\n");


	if (sp_mask == -1)
	{
		int num_sp = atoi(argv[argc-2]);
		ASSERTE(num_sp > 0, "The number of species should be greater than 0\n");
		sp_mask = (((unsigned long long int) 1) << num_sp) - 1;
	}

	if ((sp_mask & 1) == 0)
	{
		/* we should always read in the first species, so we can figure out where to cut the
		 * alignment. keep track here so we know to skip it when printing out */
		skip_first = 1;
		sp_mask = sp_mask | 1;
	}


	num_files = read_rfile(argv[argc-1], &files);

	/* setup alignment */
	alignment->num_sp = bitstr_num_nonzero(sp_mask); /* number of species in the alignment... */
	alignment->alloc = 1000000; 
	alignment->seqs = (char*) malloc (alignment->alloc*alignment->num_sp*sizeof(char));
	alignment->names = (char*) malloc ((MAXSEQNAMELEN+1)*alignment->num_sp*sizeof(char));
	alignment->first_file = alignment->last_file = -1;

	/* read in regions and output alignments */
	while ((line_len = getline(&line, &line_alloc, stdin)) != -1)
	{
		int cut_start, cut_end, sp;
		int skip;
		char strand;
		int num_match;
		unsigned char flip_strand = 0;
		int search_start = 0;

		char cur_chr[MAXCHRNAMELEN+1];
		int cur_start, cur_end;

		if (line[line_len - 1] == '\n')
			line[--line_len] = '\0';

		num_match = sscanf(line, "%*s %s %u %u %c", cur_chr, &cur_start, &cur_end, &strand);
		ASSERTE(cur_start <= cur_end, "%s has start follow end.\n", line);

		/* cut out everything except the name of the line */
		line[strcspn(line, " \t")] = '\0'; 

		if (num_match == 4 && strand == '-')
			flip_strand = 1;

		/* load the alignment for this region */
		while (load_alignment(files, num_files, sp_mask, cur_chr, cur_start, cur_end, alignment, search_start))
		{
			int overlap_start, overlap_end, cut_len;

			/* update where to start looking next time */
			search_start = alignment->last_file + 1;

			/* extract overlapping region from this alignment */
			overlap_start = MAX(alignment->start, cur_start);
			overlap_end = MIN(alignment->end, cur_end);
			
			/* find the start of the cut */
			skip = overlap_start - alignment->start;
			for (cut_start=0; (skip > 0) || alignment->seqs[SEQidx(0,cut_start)]=='-'; cut_start++)
				if (alignment->seqs[SEQidx(0,cut_start)]!='-')
					skip--;

			skip = overlap_end - overlap_start + 1;
			for (cut_end = cut_start; skip > 0; cut_end++)
				if (alignment->seqs[SEQidx(0,cut_end)]!='-')
					skip--;
			cut_end--;

			cut_len = cut_end - cut_start + 1;

			/* take the cut region from each alignment and print it*/

			/* start with the header */
			if (num_match == 4)
				printf("%s %s %u %u %c\n", line, cur_chr, overlap_start, overlap_end, strand);
			else
				printf("%s %s %u %u\n", line, cur_chr, overlap_start, overlap_end);

			/* reallocate space as necessary for buffer/print_pos... keep these the same size out of convenience */
			if ((cut_len + 1) > buffer_alloc)
			{
				while ((cut_len + 1) > buffer_alloc)
					buffer_alloc *= 2;
				buffer = (char*) realloc (buffer, buffer_alloc * sizeof(char));
				print_pos = (unsigned char*) realloc (print_pos, buffer_alloc * sizeof(unsigned char));
			}

			/* indicate bases to keep (must have non-gap in at least one of the printed species */
			memset(print_pos, 0, sizeof(unsigned char)*cut_len);
			for (sp=skip_first; sp<alignment->num_sp; sp++)
			{
				int st = SEQidx(sp,cut_start);

				for (i=0; i<cut_len; i++)
					if (alignment->seqs[st + i] != '-')
						print_pos[i] = 1;
			}

			/* output alignment */
			for (sp=skip_first; sp<alignment->num_sp; sp++)
			{
				int st = SEQidx(sp,cut_start);
				int j = 0;

				if (flip_strand)
				{
					for (i=cut_len-1; i>=0; i--)
						if (print_pos[i])
							switch (alignment->seqs[st + i])
							{
								case 'a': buffer[j++] = 't'; break;
								case 't': buffer[j++] = 'a'; break;
								case 'g': buffer[j++] = 'c'; break;
								case 'c': buffer[j++] = 'g'; break;
								case 'A': buffer[j++] = 'T'; break;
								case 'T': buffer[j++] = 'A'; break;
								case 'G': buffer[j++] = 'C'; break;
								case 'C': buffer[j++] = 'G'; break;
								default:  buffer[j++] = alignment->seqs[st + i];
							}
				}
				else
					for (i=0; i<cut_len; i++)
						if (print_pos[i])
							buffer[j++] = alignment->seqs[st + i];

				buffer[j] = '\0';

				printf("%s\n%s\n", alignment->names+SEQNAMEidx(sp), buffer);
			}

			printf("\n");
		}
	}

	for (i=0; i<num_files; i++)
		free(files[i].name);

	free(files);
	free(buffer);
	free(alignment->names);
	free(alignment->seqs);
	free(alignment);
	free(print_pos);
	return 0;
}

int load_alignment(struct rfile* files, int num_files, long long int sp_mask, char *chr, int start, int end, struct alignment* alignment, int search_start)
{
	int f1, f2, i, bases;

	for (f1=search_start; f1<num_files; f1++)
		if (STREQ(files[f1].chr, chr) && OVERLAPS(files[f1].start,files[f1].end,start,end))
			break;

	if (f1 >= num_files)
		return 0;

	/* see if additional loaded files are necessary must be perfectly adjacent 
	 * (requires that input files be sorted) */
	for (f2=f1+1; 
		f2 < num_files && STREQ(files[f2].chr, chr) && (files[f2-1].end + 1) == files[f2].start && OVERLAPS(files[f2].start,files[f2].end,start,end);
		f2++);
	f2--;

	/* no need to reload already loaded file (if we have /more/ loaded at the end it shouldn't
	 * really slow down queries, so that is ok) */
	if (f1 == alignment->first_file && f2 <= alignment->last_file)
		return 1;

	if (f1 == alignment->first_file && f2 > alignment->last_file)
		/* if a prefix of the specific files are loaded, go ahead and continue to load */
		i = alignment->last_file + 1;
	else
	{
		/* otherwise, start over */
		alignment->first_file = alignment->last_file = -1;
		alignment->length = 0;
		i = f1;
	}

	for (; i<=f2; i++)
		if (read_alignment(files[i].name, sp_mask, alignment, alignment->length) == 1)
			return 0;

	/* make sure that the number of aligned bases is correct */
	for (i=0, bases=0; i<alignment->length; i++)
		if (alignment->seqs[SEQidx(0,i)] != '-')
			bases++;

	if (bases != (files[f2].end - files[f1].start + 1))
		return 0;

	alignment->start = files[f1].start;
	alignment->end = files[f2].end;
	alignment->first_file = f1;
	alignment->last_file = f2;

	return 1;
}

int read_alignment(char* filename, unsigned long long int sp_mask, struct alignment* alignment, int offset)
{
	int loc=0, sp;
	int fail=0;
	char* line = NULL;
	FILE *fp = gzpopen(filename);
	size_t line_alloc = 0;
	int line_len = 0;

	ASSERTE(fp != NULL, "Cannot open: %s\n", filename);

	/* shift over one so we can shift back over on the first > 
	 * we gain an extra bit because we use unsigned long long int here */
	sp_mask = sp_mask << 1;

	/* read in all the sequences with the -'s remaining */
	for(sp=-1; (line_len = getline(&line, &line_alloc, fp)) != -1;)
	{
		if (line[0] == '>')
		{
			if (loc != 0)
				alignment->length = loc;

			sp_mask = sp_mask >> 1;

			/* if we have read in all the sequences, short-circuit */
			if (sp_mask == 0)
				break;
		}

		if (!(sp_mask & 1))
			continue;

		if (line[0] == '>')
		{
			sp++;

			/* fix up the name a little bit... */
			line[strcspn(line, "\n |")] = '\0';
			strncpy(alignment->names+SEQNAMEidx(sp), line, MAXSEQNAMELEN);
			alignment->names[SEQNAMEidx(sp) + MAXSEQNAMELEN] = '\0';
			alignment->seqs[SEQidx(sp,offset)] = '\0';
			loc = offset;
		}
		else
		{
			/* throw away trailing \n */
			if (line[line_len-1] == '\n')
				line[--line_len] = '\0';

			while ((line_len + loc + 1) > alignment->alloc)
			{
				/* double space assigned to alignment */
				int i;
				alignment->seqs = (char*) realloc (alignment->seqs, 2*alignment->alloc*alignment->num_sp*sizeof(char));
				for (i=alignment->num_sp-1; i>0; i--)
					memcpy(alignment->seqs + SEQidx(i*2, 0), alignment->seqs + SEQidx(i, 0), alignment->alloc);
				alignment->alloc *= 2;
			}

			/* notice the \0 is also copied */
			memcpy(alignment->seqs + SEQidx(sp,loc), line, line_len+1);

			/* loc always corresponds to the location of the \0 */
			loc += line_len;
		}
	}

	if (loc != 0)
		alignment->length = loc;
	
	for (sp=0; sp < alignment->num_sp; sp++)
	{
		if (alignment->seqs[SEQidx(sp,offset)] == '\0')
		{
			/* empty alignment means put all gaps... */
			memset(alignment->seqs + SEQidx(sp,offset), '-', alignment->length*sizeof(char)-offset);
			alignment->seqs[SEQidx(sp,alignment->length)] = '\0';
		}
		else if (strlen(alignment->seqs + SEQidx(sp,0)) != ((unsigned) alignment->length))
			goto failure;
	}

end:
	pclose(fp);
	free(line);
	return fail;

failure:
	fprintf(stderr, "Alignment file malformed %s\n", filename);
	fail=1;
	goto end; 
}


unsigned char bitstr_num_nonzero (long long int str)
{
	int len=0;
	for (; str!=0; str>>=1)
		len += str & 1;
	return len;
}


