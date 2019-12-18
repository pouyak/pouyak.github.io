/* Written by Pouya Kheradpour (pouyak@gmail.com)
 * Copyright (c) 2005-2015 Massachusetts Institute of Technology
 * Released under the MIT license
 * https://pouya.kheradpour.com/software/LICENSE.mit */
/* version 20160416 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "pk.h"
#include "pk-bl.h"

#define MAXMOTLEN			200
#define NUMCHAR				16 /* including degeneracies */

/* subtracted from cutoffs so that we get values that equal cut exactly */
#define EPSILON		1e-6

/* indexing for alignment.index_corr */
#define ICidx(sp,loc)		((num_sp-1) * (loc) + ((sp)-1))

/* indexing for alignment.seqs */
#define SEQidx(sp, loc)		(num_sp*(loc) + (sp))

#define BITON(bstr, bit)		(((bstr) & (1 << ((bit)-1))) != 0)

/* used to index PWMs... notice, because of base order we can use a trick to do reverse lookups */
#define RFMidx(pos,chr,len,rev) ((rev) ? RMidx(pos,chr,len) : Midx(pos,chr))
#define RMidx(pos,chr,len) ((4*len) - (4*(pos) + (chr)) - 1)
#define Midx(pos,chr)	(4*(pos) + (chr))

/* indexing for match_positions; note: species 0 is not done */
#define MPidx(sp, dir)		(2*((sp)-1) + (dir))

/* order of bases is (LSB to MSB): A C G T */
const unsigned char alpha[26] = {1,14,2,13,0,0,4,11,0,0,12,0,3,15,0,0,0,5,6,8,0,7,9,0,10,0};
const unsigned char alpha_rev[26] = {8,7,4,11,0,0,2,13,0,0,3,0,12,15,0,0,0,10,6,1,0,14,9,0,5,0};
const unsigned char alpha_num_on[26] = {1,3,1,3,0,0,1,3,0,0,2,0,2,4,0,0,0,2,2,1,0,3,2,0,2,0};

const char *bin2char = "nACMGRSVTWYHKDBN";
const char *bin2char_rev = "nTGKCYSBAWRDMHVN";

/* pwm representation of motifs; 4 floats per position */
typedef float pwmrep[MAXMOTLEN*4];

/* string containing bit representation of motif for fwd and rev directions */
typedef unsigned char bitrep[2][MAXMOTLEN+1];

/* structure holding motif... 
 * fwd and rev are both arrays where each position is a bitstring representing the
 * admissible bases */
struct motif
{
	char name[MAXMOTLEN+1];
	float blcut;
	int window;
	unsigned char allow_flip;
	unsigned char len;

	/* consensus sequence representations of the motif */
	/* bit representation of motif, 0 for fwd, 1 for rev */
	int ns;
	bitrep *s; 

	/* pwm representations of the motif */
	int np; 
	pwmrep *p;
};

/* stores a list of motifs */
struct motif_list
{
	unsigned int num;
	unsigned int* vals;
};

struct alignment
{
	int* lens;

	/* indexed using SEQidx(species, location) */
	unsigned char* seqs; /* ACGTN -> 1, 2, 4, 8, 15 */

	/*  indicates seqs[SEQidx(i,index_corr[ICidx(i,j)])] is aligned to seqs[SEQidx(0,j)] */
	/*  note: if j has a gap in the given base, then the NEXT ungapped position is used */
	int* index_corr;
};

unsigned char* parse_consm_sp_str(char* consm_sp_str, unsigned char num_sp);
unsigned int bitstr_length (unsigned int motif_bitstr);
unsigned int create_bitstr_motifs(unsigned int motif_bitstr, struct motif** motifs_ptr, int def_window, unsigned char def_allow_flip, float def_blcut);
struct motif_list *create_motif_hashes(struct motif* motifs, unsigned int num_motifs, int hash_len);
unsigned int read_motifs(char* motif_file, struct motif** motifs_ptr, int def_window, unsigned char def_allow_flip, float def_blcut);
int find_ali_hash(struct alignment* alignment, int num_sp, int loc, int hash_len);
int find_bitstr_num(struct alignment* alignment, int num_sp, unsigned int motif_bitstr, unsigned char motif_bitstr_len, unsigned char strand, int loc);
int read_alignment(struct alignment **alignment_ptr, char ***names_ptr, unsigned char flip_strand);
unsigned char *cons_match(struct alignment* alignment, int num_sp, int window, unsigned char unique_wmatch, struct motif* motif, int loc, int* match_positions, float* match_scores);

/* prints matching motifs */
void print_motif_match(struct alignment* alignment, char **names, int num_sp, struct consm* consm, unsigned char unique_wmatch, unsigned char ali_scan_strands, int loc, char* chr, unsigned int start, struct motif* motif, long long int match_mask, int print_numflank, unsigned char print_match_score);

int main (int argc, char** argv)
{
	struct motif* motifs;
	char* consm_sp_str = NULL;
	char* bl_file = NULL;
	char* blp_file = NULL;
	char* motif_file = NULL;
	float def_blcut = 0.0;
	struct consm consm = {NULL, NULL, NULL, NULL, NULL};
	unsigned int num_motifs = 0;
	int num_sp = 0;
	int def_window = 0;
	int scan_strands = 3;
	unsigned char def_allow_flip = 0;
	unsigned int motif_bitstr = 0;
	unsigned char motif_bitstr_len = 0;
	unsigned char unique_wmatch = 1;
	long long int match_mask = -1;
	struct motif_list *motif_hash_lookup = NULL;
	int print_numflank = -1;
	unsigned char print_match_score = 0;
	unsigned char pop_mode = 0;
	int hash_len = 8;

	char *line = NULL;
	size_t line_alloc = 0;
	int line_len;

	unsigned int i;

	/* parse command line arguments */
	for (i=1; i<((unsigned int) argc-1); i++)
		if (STREQ(argv[i], ""))
			continue; 
		else if (STREQ(argv[i], "-m"))
			motif_file = argv[++i];
		else if (STREQ(argv[i], "-M"))
			motif_bitstr = (1 << atoi(argv[++i]))-1;
		else if (STREQ(argv[i], "-r"))
			motif_bitstr = atoi(argv[++i]);
		else if (STREQ(argv[i], "-n"))
			num_sp = atoi(argv[++i]);
		else if (STREQ(argv[i], "-w"))
			def_window = atoi(argv[++i]);
		else if (STREQ(argv[i], "-h"))
			hash_len = atoi(argv[++i]);
		else if (STREQ(argv[i], "-c"))
			consm_sp_str = argv[++i];
		else if (STREQ(argv[i], "-b"))
			bl_file = argv[++i];
		else if (STREQ(argv[i], "-bp"))
			blp_file = argv[++i];
		else if (STREQ(argv[i], "-l"))
			def_blcut = atof(argv[++i]) != 0;
		else if (STREQ(argv[i], "-s"))
			def_allow_flip = atoi(argv[++i]) != 0;
		else if (STREQ(argv[i], "-t"))
			scan_strands = atoi(argv[++i]);
		else if (STREQ(argv[i], "-k"))
			match_mask = atoll(argv[++i]);
		else if (STREQ(argv[i], "-u"))
			unique_wmatch = atoi(argv[++i]) != 0;
		else if (STREQ(argv[i], "-v"))
			print_numflank = atoi(argv[++i]);
		else if (STREQ(argv[i], "-V"))
			print_match_score = atoi(argv[++i]) != 0;
		else
			ASSERTE(0, "Unrecognized command line argument!\n");

	if (num_sp == -1)
	{
		def_window = -1;
		def_allow_flip = 0;
		print_numflank = -1;

		pop_mode = 1;
	}

	ASSERTE(scan_strands == 1 || scan_strands == 2 || scan_strands == 3, "-t must be either 1, 2 or 3\n");

	/* we should be able to get around this in the future */
	ASSERTE(num_sp < 64, "Cannot have more than 63 species\n");

	if (num_sp == 0 || (motif_file == NULL && motif_bitstr == 0))
	{
		printf("USAGE: %s [OPTIONS] < ALIGNFILE\n", argv[0]);
		printf("    -n    Number of species; -1 for population output mode [required]\n");
		printf("    -m    Motif file (format: Motifs[; sep] [name window flip blcut]; string of 'X's indicate PWM to follow)) [either -m, -M or -r required]\n");
		printf("    -M    Scan for all non-degenerate n-mers [either -m, -M or -r required]\n");
		printf("    -r    Scan for all non-degenerate k-mers specified by bitstring mask [either -m, -M or -r required]\n");
		printf("    -w    Default match window; a window of -1 indicates that gaps should not be moved to the next aligned base [default: 0]\n");
		printf("    -s    Default whether a match is considered conserved if it occurs on opposite of the aligned strand [default: 0]\n");
		printf("    -t    Strands of the alignment to scan (1 for fwd, 2 for rev, 3 for both) [default: 3]\n");
		printf("    -c    Comma separated list of integer value representing species conservation is required for match [default: all]\n");
		printf("    -b    Branch length file. If supplied then a branch length is added to the end of each match\n");
		printf("    -bp   Branch length file in parent tree format\n");
		printf("    -l    Default branch length cutoff which is considered a match (used only if -c not supplied, -b/-bp required) [default: 0]\n");
		printf("    -u    When using windows/flips, allow each instance in another genome to count only once [default: 1]\n");
		printf("    -v    Print the location/sequence of matches in the genome (# indicates number of flanking bases to include; -1 to disable, -2 for just locations) [default: -1]\n");
		printf("    -V    Print the PWM match score: (score - cutoff) / (max_score - cutoff) (for consensus, always prints 1; not well defined for compound motifs) [default: 0]\n");
		printf("    -k    Value to AND against match (effectively, restricts to species of indicated) (-1 for none) [default: -1]\n");
		printf("    -h    Hash length to use for initial motif matching (higher number faster but more memory) [default: 8]\n");
		return 1;
	}

	ASSERTE(!pop_mode || (bl_file == NULL && blp_file == NULL && consm_sp_str == NULL && match_mask == -1), "When in population mode, cannot specify tree or match mask\n");

	if (motif_bitstr > 0)
	{
		num_motifs = create_bitstr_motifs(motif_bitstr, &motifs, def_window, def_allow_flip, def_blcut);
		motif_bitstr_len = bitstr_length(motif_bitstr);
	}
	else
	{
		num_motifs = read_motifs(motif_file, &motifs, def_window, def_allow_flip, def_blcut);
		motif_hash_lookup = create_motif_hashes(motifs, num_motifs, hash_len);
	}

	/* read in appropriate branch length files */
	if (bl_file != NULL)
		consm.bl = read_bl_file(bl_file, num_sp);
	else if (blp_file != NULL)
		read_blp_file(blp_file, num_sp, &consm);

	/* get array corresponding to species needed for a "match" */
	if (consm_sp_str != NULL)
		consm.sp = parse_consm_sp_str(consm_sp_str, num_sp);

	if (match_mask == -1)
		match_mask = (1 << (num_sp-1)) - 1;

	/* read in files and perform motif matching */
	while ((line_len = getline(&line, &line_alloc, stdin)) != -1)
	{
		char chr[MAXCHRNAMELEN+1];
		unsigned int start;
		int loc;
		unsigned char ali_scan_strands = scan_strands;
		unsigned char num_match;
		unsigned char strand;
		char c_strand;
		struct alignment *alignment;
		char **names = NULL;

		if (line[line_len - 1] == '\n')
			line[--line_len] = '\0';

		num_match = sscanf(line, "%*s %s %u %*d %c", chr, &start, &c_strand);

		strand = (num_match == 3 && c_strand == '-') ? 2 : 1;

		/*NOTE: scan_strand corresponds to the strand of the alignment file; however, the alignment file
		 * is always flipped so that its on the plus strand; thus ali_scan_strand specifies which *genome*
		 * strands to scan */

		/* with align_file, the strand is not always 1, but we want the ali_scan_strand to 
		 * always be the one we read in so we set it accordingly */
		if (strand == 2 && scan_strands != 3)
			ali_scan_strands = 3 - scan_strands;

		/* get the alignment; notice, if strand == 2, we flip the alignment to get it back to plus strand */
		if (pop_mode)
			num_sp = read_alignment(&alignment, &names, strand == 2);
		else
			ASSERTE(num_sp == read_alignment(&alignment, NULL, strand == 2), "Number of species in alignment inconsistent!\n");

		/* print all conserved matches for each motif */
		if (motif_bitstr == 0)
			for (loc=0; loc<alignment->lens[0]; loc++)
			{
				unsigned int m, m_hash = find_ali_hash(alignment, num_sp, loc, hash_len);
				
				for (m=0; m<motif_hash_lookup[m_hash].num; m++)
					print_motif_match(alignment, names, num_sp, &consm, unique_wmatch, ali_scan_strands, loc, chr, start, &motifs[motif_hash_lookup[m_hash].vals[m]], match_mask, print_numflank, print_match_score);
			}
		else if (ali_scan_strands == 3)
			for (loc=0; loc<alignment->lens[0]; loc++)
			{
				int f_num = find_bitstr_num(alignment, num_sp, motif_bitstr, motif_bitstr_len, 1, loc);
				int r_num = find_bitstr_num(alignment, num_sp, motif_bitstr, motif_bitstr_len, 2, loc);

				if (f_num != -1)
					print_motif_match(alignment, names, num_sp, &consm, unique_wmatch, 3, loc, chr, start, &motifs[f_num], match_mask, print_numflank, print_match_score);

				/* do not print the same motif twice */
				if (r_num != -1 && r_num != f_num)
					print_motif_match(alignment, names, num_sp, &consm, unique_wmatch, 3, loc, chr, start, &motifs[r_num], match_mask, print_numflank, print_match_score);
			}
		else
			for (loc=0; loc<alignment->lens[0]; loc++)
			{
				int f_num = find_bitstr_num(alignment, num_sp, motif_bitstr, motif_bitstr_len, ali_scan_strands, loc);
				if (f_num != -1)
					print_motif_match(alignment, names, num_sp, &consm, unique_wmatch, ali_scan_strands, loc, chr, start, &motifs[f_num], match_mask, print_numflank, print_match_score);
			}

		free(alignment->lens);
		free(alignment->seqs);
		free(alignment->index_corr);
		free(alignment);
	}

	/* free up memory for cleanliness */
	free(line);
	if (motif_hash_lookup != NULL)
	{
		unsigned int num_hashes = 1 << (2 * hash_len);
		for (i=0; i<num_hashes; i++)
			free(motif_hash_lookup[i].vals);
		free(motif_hash_lookup);
	}
	free(consm.bl);
	free(consm.sp);
	free(consm.mat);
	free(consm.tlist);
	free(consm.len);
	for (i=0; i<num_motifs; i++)
	{
		free(motifs[i].s);
		free(motifs[i].p);
	}
	free(motifs);
	return 0;
}

/* prints loc+i using conv string; when i is outside of range of len, prints out lowercase */
/* NOTE: x's are printed when flanking is outside of range of alignment */
#define GET_SEQ_HELPER(ch, conv, sp, loc, i, len) do {									\
	if (((loc)+(i))>=0 && ((loc)+(i)) < alignment->lens[sp])							\
	{																					\
		if ((i)<0 || (i)>=(len))														\
			(ch) =  tolower(conv[alignment->seqs[SEQidx((sp),(loc)+(i))]]);				\
		else																			\
			(ch) =  conv[alignment->seqs[SEQidx((sp),(loc)+(i))]];						\
	}																					\
	else																				\
		(ch) = 'x';																		\
} while (0)

/* prints sequence len bases from loc from species in indicated dir +/- flank */
#define GET_SEQ(str, sp, loc, dir, len, flank) do {										\
	int i;																				\
	for (i=0; i<((len)+2*flank); i++)													\
		if ((dir)==0)																	\
			GET_SEQ_HELPER((str)[i], bin2char, sp, loc, i-flank, len);					\
		else																			\
			GET_SEQ_HELPER((str)[i], bin2char_rev, sp, loc, flank+(len)-1-i, len);		\
	(str)[i] = '\0';																	\
} while (0)																				

struct strsp {
	char* str;
	int sp;
};

int strsp_cmp (const void *a, const void *b)
{
	return strcmp(((struct strsp*) a)->str, ((struct strsp*) b)->str);
}

void print_motif_match(struct alignment* alignment, char **names, int num_sp, struct consm* consm, unsigned char unique_wmatch, unsigned char ali_scan_strands, int loc, char* chr, unsigned int start, struct motif* motif, long long int match_mask, int print_numflank, unsigned char print_match_score)
{
	long long int match_bitstr[2];
	unsigned char strand;

	int* match_positions = (print_numflank != -1) ? ((int*) malloc((num_sp-1) * 2 * sizeof(int))) : NULL;
	float* match_scores = print_match_score ? ((float*) malloc(2 * sizeof(float))) : NULL;

	/* NOTE: this uses a unified cons_match that no longer directly uses match_mask (which
	 * may slow down some executions), and also returns the match array, rather than
	 * match bitstrings (which we recreate below) */
	unsigned char *match = cons_match(alignment, num_sp, motif->window, unique_wmatch * (motif->allow_flip + 1), motif, loc, match_positions, match_scores);

	{
		/* create match_bitstr for backwards compatibility */
		int i, d;

		for (d=0; d<2; d++)
		{
			match_bitstr[d] = 0;
			for (i=1; i<num_sp; i++)
				match_bitstr[d] = match_bitstr[d] | (((match[i] & (1 << d)) != 0) << (i-1));
			match_bitstr[d] = match_bitstr[d] & match_mask;
		}
	}

	for (strand = 1; strand <= 2; strand++)
		/* make sure to only print matches to the desired alignment strand if the match was found */
		if ((match[0] & ali_scan_strands & strand))
		{
			long long int str_match_bitstr = (motif->allow_flip ? (match_bitstr[0] | match_bitstr[1]) : match_bitstr[strand-1]);

			float bl; 

			if (consm->sp != NULL && !consm->sp[str_match_bitstr])
				continue;
			
			bl = get_bl(num_sp, consm, str_match_bitstr);

			if (bl < -0.5 || bl >= (motif->blcut-EPSILON))
			{
				printf("%s %s %d %d %c", motif->name, chr, start + loc, start + loc + motif->len - 1, (strand == 1) ? '+' : '-');

				if (names == NULL)
				{
					printf(" %Ld", str_match_bitstr);

					if (bl > -0.5)
						printf(" %f", bl);

					if (print_numflank >= 0)
					{
						int sp;
						char *temp_seq = (char *) malloc(sizeof(char) * (motif->len + 2 * print_numflank + 1));

						/* there must be a match on the same strand; print it */
						GET_SEQ(temp_seq, 0, loc, strand-1, motif->len, print_numflank);
						printf(" 0:+:%s", temp_seq);

						for (sp=1; sp<num_sp; sp++)
						{
							printf(";");
							if (match[sp] & (motif->allow_flip ? 3 : strand))
							{
								int sp_loc = alignment->index_corr[ICidx(sp,loc)]; 
								/* if flips are not allowed, print only matches on the same strand otherwise, 
								 * print the closest one (breaking ties by going with same strand */
								if (!motif->allow_flip || !(match[sp] & (3-strand)) || ((match[sp] & strand) && (abs(match_positions[MPidx(sp,strand-1)]-sp_loc) <= abs(match_positions[MPidx(sp,2-strand)]-sp_loc))))
								{
									GET_SEQ(temp_seq, sp, match_positions[MPidx(sp,strand-1)], strand-1, motif->len, print_numflank);
									printf("%d:+:%s", match_positions[MPidx(sp,strand-1)]-sp_loc, temp_seq);
								}
								else
								{
									GET_SEQ(temp_seq, sp, match_positions[MPidx(sp,2-strand)], 2-strand, motif->len, print_numflank);
									printf("%d:-:%s", match_positions[MPidx(sp,2-strand)]-sp_loc, temp_seq);
								}
							}
							else if (motif->window == -2)
							{
								/* if window == -2, print the value regardless */
								GET_SEQ(temp_seq, sp, alignment->index_corr[ICidx(sp,loc)], strand-1, motif->len, print_numflank);
								printf("%d:+:%s", 0, temp_seq);
							}
						}

						free(temp_seq);
					}
					else if (print_numflank == -2)
					{
						int sp; 

						for (sp=1; sp<num_sp; sp++)
						{
							printf("%c", sp == 1 ? ' ' : ';');
							if (match[sp] & (motif->allow_flip ? 3 : strand))
							{
								int sp_loc = alignment->index_corr[ICidx(sp,loc)]; 
								/* if flips are not allowed, print only matches on the same strand otherwise, 
								 * print the closest one (breaking ties by going with same strand */
								if (!motif->allow_flip || !(match[sp] & (3-strand)) || ((match[sp] & strand) && (abs(match_positions[MPidx(sp,strand-1)]-sp_loc) <= abs(match_positions[MPidx(sp,2-strand)]-sp_loc))))
									printf("%d", match_positions[MPidx(sp,strand-1)]-sp_loc);
								else
									printf("%d", match_positions[MPidx(sp,2-strand)]-sp_loc);
							}
						}
					}
				}

				if (print_match_score)
					printf(" %f", match_scores[strand-1]);

				if (names != NULL)
				{
					/* code specific to population genetics analysis */
					int sp; 
					char *temp_seq = (char *) malloc(sizeof(char) * (motif->len + 1) * num_sp);
					struct strsp *seqs = (struct strsp*) malloc(sizeof(struct strsp) * (num_sp-1));

					/* output format is:
					 * REF_SEQ|ALT_SEQ1:MOTIF MATCH?:IND1;IND2;IND3...|ALT_SEQ2:...
					 */
					GET_SEQ(temp_seq, 0, loc, strand-1, motif->len, 0);
					printf(" %s", temp_seq);

					/* get the sequences for each of the alternate chromosomes */
					for (sp=1; sp<num_sp; sp++)
					{
						seqs[sp-1].str = temp_seq + sp * (motif->len + 1);
						seqs[sp-1].sp = sp;
						
						/* make sure this position has something aligned to it (i.e. not a gap), because we use window < 0 */
						if ((loc+1) < alignment->lens[0] && alignment->index_corr[ICidx(sp,loc)] == alignment->index_corr[ICidx(sp,loc+1)])
						{
							seqs[sp-1].str[0] = '-';
							seqs[sp-1].str[1] = '\0';
						}
						else
							GET_SEQ(seqs[sp-1].str, sp, alignment->index_corr[ICidx(sp,loc)], strand-1, motif->len, 0);
					}

					qsort(seqs, num_sp-1, sizeof(struct strsp), strsp_cmp);

					/* print out list of individuals with changes to their match sequence */
					for (sp=1; sp<num_sp; sp++)
						if (STREQ(seqs[sp-1].str, temp_seq))
							ASSERTE(match[seqs[sp-1].sp] & strand, "Inconsistency in sequence matching! %s %d\n", chr, start + loc);
						else if (sp > 1 && STREQ(seqs[sp-1].str, seqs[sp-2].str))
						{
							ASSERTE((match[seqs[sp-1].sp] & strand) == (match[seqs[sp-2].sp] & strand), "Inconsistency in sequence matching! %s %d\n", chr, start + loc);
							printf(";%s", names[seqs[sp-1].sp]);
						}
						else
							printf("|%s:%c:%s", seqs[sp-1].str, (match[seqs[sp-1].sp] & strand) ? '1' : '0', names[seqs[sp-1].sp]);

					free(temp_seq);
					free(seqs);
				}

				printf("\n");
			}
		}

	free(match);
	free(match_positions);
	free(match_scores);
}

#undef GET_SEQ
#undef GET_SEQ_HELPER

/* parses the comma separated list of species required for a cons match and produces   
 * a 0/1 array where given the bitstring set of species with a match indicates whether
 * the match is conserved ; this is deprecated */
unsigned char* parse_consm_sp_str(char* consm_sp_str, unsigned char num_sp)
{
	char* last_start;
	unsigned int num_sp_comb = 1 << (num_sp-1);
	unsigned char* consm_sp = (unsigned char*) calloc(num_sp_comb, sizeof(unsigned char));
	
	/* go through consm_sp_str looking for comma's */
	for (last_start=consm_sp_str; (*consm_sp_str) != '\0'; consm_sp_str++)
		if ((*consm_sp_str) == ',')
		{
			/* as long as we do not have an empty range, convert the last number */
			if (last_start != consm_sp_str)
				consm_sp[atoi(last_start)] = 1;
			last_start = consm_sp_str+1;
		}
		else
			ASSERTE((*consm_sp_str) >= '0' && (*consm_sp_str) <= '9', "Malformed species conservation; only numbers and commas permitted.\n");

	/* read a final number if the string does not end with a comma */
	if ((*last_start) != '\0') 
		consm_sp[atoi(last_start)] = 1;

	return consm_sp;
}

unsigned int bitstr_num_nonzero (unsigned int motif_bitstr)
{
	int len=0;
	for (; motif_bitstr!=0; motif_bitstr>>=1)
		len += motif_bitstr & 1;
	return len;
}

/* returns the length of a bitstring floor(log2(motif_bitstr)) */
unsigned int bitstr_length (unsigned int motif_bitstr)
{
	int len=0;
	for (; motif_bitstr!=0; motif_bitstr>>=1, len++);
	return len;
}

/* creates all motifs specified by the bitstring in a systematic way */
unsigned int create_bitstr_motifs(unsigned int motif_bitstr, struct motif** motifs_ptr, int def_window, unsigned char def_allow_flip, float def_blcut)
{
	unsigned int bitstr_numnz = bitstr_num_nonzero(motif_bitstr);
	unsigned int bitstr_len = bitstr_length(motif_bitstr);
	unsigned int num_motifs = 1 << (2 * bitstr_numnz), i, j, k; 
	struct motif* motifs = (struct motif*) malloc(num_motifs * sizeof(struct motif));

	/* 0123 -> (ACGT)1248 */
	const unsigned char base2seq[4] = {1,2,4,8};
	const char base2char[4] = "ACGT";

	/* NOTE: the hash is "stored" in the opposite direction as the motif hash (here it is stored left to right) */
	for (i=0; i<num_motifs; i++)
	{
		motifs[i].len = bitstr_len;
		motifs[i].s = (bitrep*) malloc(sizeof(bitrep));
		motifs[i].name[bitstr_len] = '\0';
		motifs[i].ns = 1;

		/* convert each base to a bitstring */
		for (j=0, k=0; j<motifs[i].len; j++)
		{
			if (motif_bitstr & (1 << (bitstr_len - j - 1)))
			{
				unsigned char base = (i >> (2*(bitstr_numnz - k - 1)))%4;
				motifs[i].name[j] = base2char[base];
				motifs[i].s[0][0][j] = base2seq[base];
				motifs[i].s[0][1][motifs[i].len-j-1] = base2seq[3-base];
				k++;
			}
			else
			{
				motifs[i].name[j] = 'N';
				motifs[i].s[0][0][j] = 15;
				motifs[i].s[0][1][motifs[i].len-j-1] = 15;
			}
		}

		motifs[i].blcut = def_blcut;
		motifs[i].window = def_window;
		motifs[i].allow_flip = def_allow_flip;
	}

	*motifs_ptr = motifs;

	return num_motifs;
}


/* checks to see if i is a power of 2... used to see if we need to duplicate size of arrays below */
unsigned char is_power2 (unsigned int i)
{
	while (i > 1)
	{
		if (i % 2)
			return 0;
		i = i / 2;
	}

	return 1;
}

void add_val_to_hash (struct motif_list *motif_hash_lookup, unsigned int hash, unsigned int num)
{
	/* make sure this value is not already in the hash table */
	/* duplicates can occur from opposite strand or from multiple submotifs */
	if (motif_hash_lookup[hash].num != 0 && motif_hash_lookup[hash].vals[motif_hash_lookup[hash].num-1] == num)
		return;

	/* allocate any additional memory required */
	if (motif_hash_lookup[hash].num == 0)
		motif_hash_lookup[hash].vals = (unsigned int*) malloc(sizeof(unsigned int));
	else if (is_power2(motif_hash_lookup[hash].num))
		motif_hash_lookup[hash].vals = (unsigned int*) realloc(motif_hash_lookup[hash].vals, 2*motif_hash_lookup[hash].num*sizeof(unsigned int));

	/* add to hash table */
	motif_hash_lookup[hash].vals[motif_hash_lookup[hash].num++] = num;
}

/*NOTE: these two functions build up the hash starting at the last position hashed and going to the first! */
/* recursive function to add ms to motif_hash_lookup */
void add_to_hash (struct motif_list *motif_hash_lookup, unsigned char hash_len, unsigned int pre_hash, unsigned char* ms, unsigned int m_num, unsigned char m_len)
{
	unsigned char j;

	if (hash_len == 0)
		add_val_to_hash (motif_hash_lookup, pre_hash, m_num);
	else 
		for (j=0; j<4; j++)
			if (hash_len > m_len || BITON(ms[hash_len-1],1+j))
				add_to_hash (motif_hash_lookup, hash_len-1, pre_hash | (j << ((hash_len-1)*2)), ms, m_num, m_len);
}

/* recursive function to add ms to motif_hash_lookup (for PWMs) */
void add_to_hash_pwm (struct motif_list *motif_hash_lookup, unsigned char hash_len, unsigned int pre_hash, float *pwm, unsigned char dir, unsigned int m_num, int m_len, float score)
{
	/* because all the PWM values are <= 0 and starts at -cut, if its ever lower than 0 we give up */
	if (score < 0)
		return;

	if (hash_len == 0)
		add_val_to_hash (motif_hash_lookup, pre_hash, m_num);
	else
	{
		unsigned char j;
		if (hash_len > m_len)
			for (j=0; j<4; j++)
				add_to_hash_pwm (motif_hash_lookup, hash_len-1, pre_hash | (j << ((hash_len-1)*2)), pwm, dir, m_num, m_len, score);
		else
			for (j=0; j<4; j++)
				add_to_hash_pwm (motif_hash_lookup, hash_len-1, pre_hash | (j << ((hash_len-1)*2)), pwm, dir, m_num, m_len, score + pwm[RFMidx(hash_len-1,j,m_len,dir)]);
	}
}

/* actually create the "hash" table... actually a list of all motifs that start (on either
 * strand) with motif represented by the index */
struct motif_list *create_motif_hashes(struct motif* motifs, unsigned int num_motifs, int hash_len)
{
	unsigned int num_hashes = 1 << (2 * hash_len);
	struct motif_list *motif_hash_lookup = (struct motif_list*) calloc(num_hashes, sizeof(struct motif_list));
	unsigned int m; 
	int n;
	unsigned char d;

	for (m=0; m<num_motifs; m++)
	{
		for (n=0; n<motifs[m].ns; n++)
			for (d=0; d<2; d++)
				add_to_hash (motif_hash_lookup, hash_len, 0, motifs[m].s[n][d], m, motifs[m].len);

		for (n=0; n<motifs[m].np; n++)
			for (d=0; d<2; d++)
				add_to_hash_pwm (motif_hash_lookup, hash_len, 0, motifs[m].p[n], d, m, motifs[m].len, 1);
	}

	/* free up extra memory allocated */
	for (m=0; m<num_hashes; m++)
		if (motif_hash_lookup[m].vals  != NULL)
			motif_hash_lookup[m].vals = (unsigned int*) realloc(motif_hash_lookup[m].vals, motif_hash_lookup[m].num * sizeof(unsigned int));

	return motif_hash_lookup;
}

/* reads in the motifs from a file and populates motifs array */
unsigned int read_motifs(char* motif_file, struct motif** motifs_ptr, int def_window, unsigned char def_allow_flip, float def_blcut)
{
	FILE* fp = fopen(motif_file, "r");
	struct motif* motifs;
	char *line = NULL;
	unsigned int i, num_reserved = 1 << 5;
	size_t line_alloc = 0;

	ASSERTE(fp != NULL, "Cannot open: %s\n", motif_file);

	/* count the number of lines in the file */
	motifs = (struct motif*) malloc(num_reserved * sizeof(struct motif));

	for (i=0; ; i++)
	{
		char *str;
		int tlen, n, num_scan, j, line_len;

		if ((line_len = getline(&line, &line_alloc, fp)) == -1)
			break;

		if (line[line_len - 1] == '\n')
			line[--line_len] = '\0';

		if (i == num_reserved)
		{
			num_reserved *= 2; 
			motifs = (struct motif*) realloc(motifs, num_reserved * sizeof(struct motif));
		}

		num_scan = sscanf(line, "%*s %s %d %hhu %f", motifs[i].name, &motifs[i].window, &motifs[i].allow_flip, &motifs[i].blcut);

		/* strip away the rest of the line after a space */
		line[strcspn(line, " \t")] = '\0';

		/* fill in all missing values */
		if (num_scan < 1)
			strcpy(motifs[i].name, line);
		if (num_scan < 2)
			motifs[i].window = def_window;
		if (num_scan < 3)
			motifs[i].allow_flip = def_allow_flip;
		if (num_scan < 4)
			motifs[i].blcut = def_blcut;

		/* count the number of submotifs */
		tlen = strlen(line);
		motifs[i].ns = (line[0] != 'X') ? 1 : 0;
		motifs[i].np = (line[0] == 'X') ? 1 : 0;
		for (j=0; j<tlen; j++)
			if (line[j] == ';')
			{
				if (line[j+1] == 'X')
					motifs[i].np++;
				else
					motifs[i].ns++;
				line[j] = '\0';
			}

		motifs[i].s = (bitrep*) malloc(sizeof(bitrep)*motifs[i].ns);
		motifs[i].p = (pwmrep*) malloc(sizeof(pwmrep)*motifs[i].np);

		/* get length of motif (use length of first submotif; check for consistency later) */
		motifs[i].len = strlen(line);
		ASSERTE(motifs[i].len <= MAXMOTLEN, "Motif %s too long!\n", line);

		/* convert each base to a bitstring */
		for (n=0, str = line; n<motifs[i].ns; n++, str += motifs[i].len+1)
		{
			ASSERTE(strlen(str) == motifs[i].len, "Motif %s has submotifs of inconsistent lengths!\n", motifs[i].name);

			/* skip over the PWMs */
			while (str[0] == 'X')
			{
				str += motifs[i].len+1;
				ASSERTE(strlen(str) == motifs[i].len, "Motif %s has submotifs of inconsistent lengths!\n", motifs[i].name);
			}

			for (j=0; j<motifs[i].len; j++)
			{
				int c_idx = toupper(str[j]) - 'A';
				ASSERTE((c_idx >= 0) && (c_idx < 26) && (alpha_num_on[c_idx]>0), "Motif %s contains an illegal character!\n", str);
				motifs[i].s[n][0][j] = alpha[c_idx];
				motifs[i].s[n][1][motifs[i].len-j-1] = alpha_rev[c_idx];
			}
		}

		/* read in the PWMs */
		for (n=0; n<motifs[i].np; n++)
		{
			float pwm_cut;
			unsigned int k;
			ASSERTE(fscanf(fp, "%f\n", &pwm_cut) == 1, "Error reading PWM cutoff of Motif %s!\n", motifs[i].name);

			for (j=0; j<motifs[i].len; j++)
			{
				float max_val; 
				ASSERTE(fscanf(fp, "%*s %f %f %f %f\n", &motifs[i].p[n][Midx(j,0)], &motifs[i].p[n][Midx(j,1)], &motifs[i].p[n][Midx(j,2)], &motifs[i].p[n][Midx(j,3)]) == 4, "Error reading PWM of Motif %s!\n", motifs[i].name);

				/* normalize so that max value in PWM row is 0 so that we can short circuit searches when total score falls below the cutoff */
				max_val = motifs[i].p[n][Midx(j,0)];
				for (k=1; k<4; k++)
					max_val = MAX(max_val, motifs[i].p[n][Midx(j,k)]);
				pwm_cut -= max_val; 
				for (k=0; k<4; k++)
					motifs[i].p[n][Midx(j,k)] -= max_val;
			}

			/* normalize so that the cutoff is -1, that way we don't have to store it */
			if (pwm_cut > -(EPSILON*motifs[i].len))
			{
				/* if the cutoff is exactly 0, then we need a perfect match to pass;
				 * using - 1e100for all non-zero values in the matrix */

				/* if the cutoff is greater than zero, then set the whole matrix to
				 * -1e100 */
				for (j=0; j<motifs[i].len; j++)
					for (k=0; k<4; k++)
						if (motifs[i].p[n][Midx(j,k)] < 0 || pwm_cut > (EPSILON*motifs[i].len))
							motifs[i].p[n][Midx(j,k)] = -1e100;
			}
			else
			{
				/* otherwise divide by |pwm_cut| = -pwm_cut */
				for (j=0; j<motifs[i].len; j++)
					for (k=0; k<4; k++)
						motifs[i].p[n][Midx(j,k)] = -motifs[i].p[n][Midx(j,k)] / pwm_cut;
			}
		}
	}

	motifs = (struct motif*) realloc(motifs, i * sizeof(struct motif));

	*motifs_ptr = motifs;

	free(line);
	fclose(fp);
	return i;
}

/* like above, but always does length hash_len, always does plus strand and does not return -1 (instead, returns /some/ possible hash) */
/* we want others to also be 0 because if there is an N in the alignment, we want to permit that to match an N in a motif
 * to maintain consistency with if a hash were not used (remember, returning an incorrect hash when no match is possible
 * is O.K. because a later step will prevent a motif from being printed) */
int find_ali_hash(struct alignment* alignment, int num_sp, int loc, int hash_len)
{
	unsigned int num=0;
	unsigned char i;

	/* ACGT(1248)->0123 (but other also to 0) */
	const signed char seq2base[16] = {0,0,1,0,2,0,0,0,3,0,0,0,0,0,0,0};

	/* notice: this hash is "reverse" from the hash used for find_bitstr_num */
	for (i=0; i<hash_len && (loc+i) < alignment->lens[0]; i++)
		num = num | (seq2base[alignment->seqs[SEQidx(0,i+loc)]] << (i*2));

	return num;		
}

/* returns the motif corresponding to the nmer starting at loc */
int find_bitstr_num(struct alignment* alignment, int num_sp, unsigned int motif_bitstr, unsigned char motif_bitstr_len, unsigned char strand, int loc)
{
	unsigned int num = 0;
	unsigned char i;

	/* ACGT(1248)->0123 (other is -1) */
	const signed char seq2base[16] = {-1,0,1,-1,2,-1,-1,-1,3,-1,-1,-1,-1,-1,-1,-1};

	if ((loc+motif_bitstr_len) > alignment->lens[0])
		return -1;

	if (strand == 1)
	{
		for (i=0; i<motif_bitstr_len; i++)
			if (motif_bitstr & (1 << (motif_bitstr_len - i - 1)))
			{
				signed char n = seq2base[alignment->seqs[SEQidx(0,i+loc)]];
				if (n == -1)
					return -1;
				
				num = (num << 2) | n;
			}
	}
	else
	{
		for (i=0; i<motif_bitstr_len; i++)
			if (motif_bitstr & (1 << (motif_bitstr_len - i - 1)))
			{
				signed char n = seq2base[alignment->seqs[SEQidx(0,(motif_bitstr_len - i - 1)+loc)]];
				if (n == -1)
					return -1;
				
				num = (num << 2) | (3-n);
			}
	}

	return num;		
}

#define ASidx(sp,loc)	((sp)*(align_len+1)+loc)

/* returns 1 for failure */
struct alignment* format_alignment(char* align_seqs, int align_len, int num_sp)
{	
	int sp, sp2;
	int loc;

	struct alignment *alignment = (struct alignment *) malloc(sizeof(alignment));

	alignment->seqs = (unsigned char*) malloc(num_sp*align_len*sizeof(unsigned char));
	alignment->index_corr = (int*) malloc((num_sp-1)*align_len*sizeof(int));
	alignment->lens = (int*) malloc(num_sp*sizeof(int));

	for (sp=0; sp<num_sp; sp++)
		alignment->lens[sp] = 0;

	/* go through the alignments */
	for (loc=0; loc<align_len; loc++)
		/* for each species, see if there is a non-gap */
		for (sp=0; sp<num_sp; sp++)
			if (align_seqs[ASidx(sp,loc)] != '-')
			{
				int c_idx = toupper(align_seqs[ASidx(sp,loc)]) - 'A';

				ASSERTE((c_idx >= 0) && (c_idx < 26) && (alpha_num_on[c_idx]==1 || alpha_num_on[c_idx]==4), "Block of alignment malformed: invalid character\n");

				alignment->seqs[SEQidx(sp,alignment->lens[sp])] = alpha[c_idx];

				if (sp == 0)
					/* find the corresponding locations in the other alignments */
					for (sp2=1; sp2<num_sp; sp2++)
						alignment->index_corr[ICidx(sp2,alignment->lens[0])] = alignment->lens[sp2];

				alignment->lens[sp]++;
			}

	return alignment;
}

#define GETLINELEN(str,len) do { (str) = NULL; size_t temp; (len) = getline(&(str), &temp, stdin); (str)[--len] = '\0'; } while (0)
#define GETLINE(str) do { int temp2; GETLINELEN(str, temp2); } while (0)

/* loads alignment file with variable number of people as produced by apply-var */
/* unlike read_alignment in motif-match, returns num_sp */
/* TODO: use static variables to reduce the number of mallocs */
int read_alignment(struct alignment **alignment_ptr, char ***names_ptr, unsigned char flip_strand)
{
	unsigned char sp;
	int align_len;

	char *align_seqs = NULL;
	char **names = NULL;
	int num_sp_alloc = 16;
	int num_sp; 

	/* get name for first sequence */
	ASSERTE(fgetc(stdin) == '>', "Block of alignment malformed: missing >\n");
	if (names_ptr != NULL)
	{
		names = (char **) malloc(sizeof(char*)); 
		GETLINE(names[0]);
	}
	else
		while (fgetc(stdin) != '\n');

	/* get first line */
	GETLINELEN(align_seqs, align_len);

	for (num_sp=1; ; num_sp++)
	{
		/* check for and discard >; also deals with blank separating line */
		if (fgetc(stdin) != '>')
			break;

		if (num_sp == 1 || num_sp == num_sp_alloc)
		{
			num_sp_alloc *= 2;
			align_seqs = (char *) realloc(align_seqs, (num_sp_alloc * (align_len + 1) + 3)* sizeof(char));

			if (names_ptr != NULL)
				names = (char **) realloc(names, num_sp_alloc * sizeof(char*)); 
		}

		/* get the rest of the > line */
		if (names_ptr != NULL)
			GETLINE(names[num_sp]);
		else
			while (fgetc(stdin) != '\n');

		fgets(align_seqs + ASidx(num_sp,0), align_len + 3, stdin);
		ASSERTE(((unsigned int) align_len) == (strlen(align_seqs + ASidx(num_sp,0))-1), "Length of lines do not match.\n");
		align_seqs[ASidx(num_sp,align_len)] = '\0'; /* remove newline */
	}

	/* reverse complement, if needed */
	if (flip_strand)
		for (sp=0; sp<num_sp; sp++)
			rev_comp(align_seqs + ASidx(sp,0), align_len);

	*alignment_ptr = format_alignment(align_seqs, align_len, num_sp);

	free(align_seqs);
	if (names_ptr != NULL)
		*names_ptr = names;
	return num_sp;
}

#undef GETLINE
#undef GETLINELEN
#undef ASidx

/* true iff B \subset D */
#define BMATCH(D,B)			(((D) | (B)) == (D))

/* notice that this short circuits */
#define DOIFMATCHES_HELP(expr_s, expr_p, dir, sp, loc, ns, np) do {														\
	unsigned char i;																									\
	int n;																												\
	for (n=0; n<(ns); n++)																								\
	{																													\
		for(i=0; i<motif->len && BMATCH(motif->s[n][dir][i], alignment->seqs[SEQidx(sp,loc+i)]); i++); 					\
		if (motif->len == i) {(expr_s); break;}																			\
	}																													\
	if (n == (ns))																										\
		for (n=0; n<(np); n++)																							\
		{																												\
			float score = 1;																							\
			for(i=0; i<motif->len; i++)																					\
			{																											\
				if (seq2base[alignment->seqs[SEQidx(sp,loc+i)]] != -1)													\
					score += motif->p[n][RFMidx(i,seq2base[alignment->seqs[SEQidx(sp,loc+i)]],motif->len,dir)];			\
				else																									\
					break;																								\
				if (score < 0)																							\
					break;																								\
			}																											\
			if (i == motif->len) {(expr_p); break;}																		\
		}																												\
} while (0)

#define DOIFMATCHES(expr, dir, sp, loc, ns, np) DOIFMATCHES_HELP(expr, expr, dir, sp, loc, ns, np)

/* sets var to whether or not motif matches in direction dir for species sp at loc 
 * requires variable i to exist in scope for usage */
#define SETIFMATCHES(var, dir, sp, loc, ns, np) DOIFMATCHES((var) = 1, dir, sp, loc, ns, np)

/* if var == -1, fills in the value indicated by loc if there is a match at loc in the main species */
#define FILLINMATCHLOC(var, dir, sp, loc, ns, np) if ((var) == -1) DOIFMATCHES((var) = (loc), dir, sp, loc, ns, np)

/* this large helper function is created so that the code is replicated for different constant values (e.g. ns/np) and
 * the appropriate optimizations can be done for these */
#define CONS_MATCH_HELPER(fill_positions, ns, np)	do {																			\
	if ((loc + motif->len) <= alignment->lens[0])																					\
		for (d=0; d<=1; d++)																										\
			DOIFMATCHES(match[0] = match[0] | (1 << d), d, 0, loc, ns, np);															\
	/* return if no match */																										\
	if (!match[0])																													\
		return match;																												\
	/* either the fwd or rev strand match... we have to scan the other sequences */													\
	for (sp=1; sp<num_sp; sp++)																										\
	{																																\
		/* loc_sp is the corresponding location in the genome under consideration */												\
		int l, search_start[2], search_end[2], loc_sp = alignment->index_corr[ICidx(sp,loc)];										\
		/* if loc_sp is the same as loc+1's loc_sp then loc is approx aligned to a gap in sp (see caveat below) */					\
		if (!allow_gap_match && (loc+1) < alignment->lens[0] && loc_sp == alignment->index_corr[ICidx(sp,loc+1)])     				\
			continue;																												\
		/* these indicate where we should start/stop scanning for motifs in the genome (for either strand)							\
		 * the values for the two strands may differ when unique_wmatch == 1 */														\
		search_start[0] = search_start[1] = (loc_sp < window) ? 0 : (loc_sp - window);												\
		/* NOTE: there is a possibility these will be negative if alignment->lens[sp] < motif->len */								\
		search_end[0] = search_end[1] = MIN(loc_sp + window, alignment->lens[sp] - motif->len);										\
		/* correct search_end/search_start when excluding regions that are covered by another instance of the same motif */			\
		if (unique_wmatch > 0)																										\
		{																															\
			/* expand region searched after loc for the same motif in the main genome */											\
			for(;  (next_match[0] == -1 || next_match[1] == -1) 																	\
				&& (loc + next_searched + motif->len) <= alignment->lens[0] 														\
				&& (alignment->index_corr[ICidx(sp,loc+next_searched)] - loc_sp) <= 2*window										\
				; next_searched++)																									\
				for (d=0; d<=1; d++)																								\
					FILLINMATCHLOC(next_match[d], d, 0, loc+next_searched, ns, np);													\
			/* expand region searched before loc */																					\
			for(;  (prev_match[0] == -1 || prev_match[1] == -1) 																	\
				&& (loc - prev_searched) >= 0 																						\
				&& (loc_sp - alignment->index_corr[ICidx(sp,loc-prev_searched)]) <= 2*window										\
				; prev_searched++)																									\
				for (d=0; d<=1; d++)																								\
					FILLINMATCHLOC(prev_match[d], d, 0, loc-prev_searched, ns, np);													\
			/* update search_start/search_end to avoid regions covered by other motifs */											\
			for (d=0; d<=1; d++)																									\
			{																														\
				/* break ties by assigning match to motif further along in alignment */												\
				if (prev_match[d] != -1)																							\
					search_start[d] = MAX(search_start[d], (loc_sp + alignment->index_corr[ICidx(sp,prev_match[d])] + 1)/2);		\
				if (next_match[d] != -1)																							\
					search_end[d] = MIN(search_end[d], (loc_sp + alignment->index_corr[ICidx(sp,next_match[d])] + 1)/2 - 1);		\
			}																														\
			/* when unique_wmatch == 2, a match in the main species on either strand can "take" a match within the window */		\
			if (unique_wmatch == 2)																									\
			{																														\
				search_start[0] = search_start[1] = MAX(search_start[0], search_start[1]);											\
				search_end[0] = search_end[1] = MIN(search_end[0], search_end[1]);													\
			}																														\
		}																															\
		/* finally, do all the searching in the other genome */																		\
		for (d=0; d<=1; d++)																										\
		{																															\
			unsigned char mt = 0;																									\
			for (l=search_start[d]; l<=search_end[d] && !mt; l++)																	\
				SETIFMATCHES(mt, d, sp, l, ns, np);																					\
			match[sp] = match[sp] | (mt << d);																						\
																																	\
			/* if fill_loc=1, then find the /closest/ match in the genome; this can be greatly optimized */							\
			if (fill_positions)																										\
			{																														\
				match_positions[MPidx(sp, d)] = -1;																					\
				if (mt)																												\
					for (l=search_start[d]; l<=search_end[d]; l++)																	\
						if ((match_positions[MPidx(sp, d)] == -1) || (abs(l-loc_sp) < abs(match_positions[MPidx(sp, d)]-loc_sp)))	\
							DOIFMATCHES(match_positions[MPidx(sp, d)] = l, d, sp, l, ns, np);										\
			}																														\
		}																															\
	}																																\
} while (0)

unsigned char *cons_match(struct alignment* alignment, int num_sp, int window, unsigned char unique_wmatch, struct motif* motif, int loc, int* match_positions, float* match_scores)
{
	/* NOTE: in this function, for all the 2 element arrays, index 0 corresponds to forward strand and 1 corresponds to the reverse strand */

	int sp, next_searched = 1, prev_searched = 1, next_match[2] = {-1,-1}, prev_match[2] = {-1,-1};
	unsigned char d;
	const signed char seq2base[16] = {-1,0,1,-1,2,-1,-1,-1,3,-1,-1,-1,-1,-1,-1,-1};
	unsigned char allow_gap_match = 1;

	/* 0 for none, 1 for fwd, 2 for rev, 3 for both */
	unsigned char *match = (unsigned char *) calloc(num_sp, sizeof(unsigned char));

	if (window < 0)
	{
		/* note: we will still allow gaps if the next base in the alignment has no corresponding base in the current alignment; e.g. A- to -A 
		 *       this is a limitation of how the data is currently stored in the alignment (it looks identical to A to A) */
		allow_gap_match = 0;
		window = 0;
	}
	
#define SWITCH_ON_NP(fill_positions, ns) do {									\
		switch (motif->np) { 													\
			case 0: CONS_MATCH_HELPER(fill_positions, ns, 0); break; 			\
			case 1: CONS_MATCH_HELPER(fill_positions, ns, 1); break;			\
			default: CONS_MATCH_HELPER(fill_positions, ns, motif->np); break;	\
	}																			\
} while (0)

#define SWITCH_ON_NS(fill_positions) do {										\
	switch(motif->ns) {															\
		case 0: SWITCH_ON_NP(fill_positions, 0); break;							\
		case 1:	SWITCH_ON_NP(fill_positions, 1); break;							\
		default: SWITCH_ON_NP(fill_positions, motif->ns); break;				\
	}																			\
} while (0)

	if (match_positions == NULL)
		SWITCH_ON_NS(0);
	else
		SWITCH_ON_NS(1);

	if (match_scores != NULL)
	{
		/* just do this the easy way for now... rematch to the target genome */
		int d;
		for (d=0; d<2; d++)
			if (match[0] & (1 << d))
				DOIFMATCHES_HELP(match_scores[d] = 1, match_scores[d] = score, d, 0, loc, motif->ns, motif->np);
	}

	return match;
}

#undef CONS_MATCH_HELPER
#undef SWITCH_ON_NP
#undef SWITCH_ON_NS
#undef FILLINMATCHLOC
#undef SETIFMATCHES
#undef DOIFMATCHES
#undef DOIFMATCHES_HELP
#undef BMATCH

