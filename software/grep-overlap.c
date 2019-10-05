/* Written by Pouya Kheradpour (pouyak@gmail.com)
 * Copyright (c) 2005-2015 Massachusetts Institute of Technology
 * Released under the MIT license
 * http://compbio.mit.edu/pouyak/software/LICENSE.mit */
/* version 20140430 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pk.h"

#define FILL_REGION(reg, fmt) if (sscanf((reg)->line, fmt, (reg)->chr, &(reg)->start, &(reg)->end, &(reg)->str) == 3 || ((reg)->str != '+' && (reg)->str != '-')) (reg)->str = 0

#define MATCH_MODE_NORMAL				0
#define MATCH_MODE_NORMAL_PLUS_NO		1
#define MATCH_MODE_ONLYMATCHES			2
#define MATCH_MODE_ONLYNONMATCHES		3
#define MATCH_MODE_PREFIX				4
#define MATCH_MODE_NORMAL_PLUS_ALL		5

char* fill_format (int num_skip); 
int col_match_parse (int **match_cols_ptr, char *str);
int check_col_match (char* s1, char* s2, int num_col_match, int* match_cols);

struct region
{
	char *line;
	char chr[MAXCHRNAMELEN+1];
	long long int start;
	long long int end;
	char str;
	struct region* next;
};

int main (int argc, char** argv)
{
	FILE *fp1, *fp2;

	/* this could have probably been implemented somewhat more quickly with some sort of
	 * sliding array to reduce the number of required mallocs... there are some subtleties 
	 * that make this non-trivial, and it looks like this is already very fast, but 
	 * something to consider */

	/* the first one on this list is used as the new item we read in */
	struct region* rlist = (struct region*) malloc (sizeof(struct region));
	size_t rlist_line_alloc = 0;

	/* we use and keep track of this rlist_end ptr so that our output is guaranteed to be
	 * first sorted by fp1 then fp2 */
	struct region* rlist_end = rlist;

	/* indicates with fp2 is empty */
	unsigned char fp2_empty = 0;

	/* use these to ensure that our input files are both sorted */
	char fp1_last_chr[MAXCHRNAMELEN+1] = "";
	long long int fp1_last_start = 0;
	char fp2_last_chr[MAXCHRNAMELEN+1] = "";
	long long int fp2_last_start = 0;
	int fp1_skip_col = 1;
	int fp2_skip_col = 1;

	char *fp1_fmt, *fp2_fmt;

	int i, padding = 0;

	int match_mode = MATCH_MODE_NORMAL;
	char sep_char = '|';
	char out_fmt[100];

	unsigned char match_r1_olen = 0;
	unsigned char match_r2_olen = 0;
	unsigned char require_str = 0; /* 2 for opposite strand */
	unsigned char start1_zeroidx = 0;
	unsigned char start2_zeroidx = 0;

	int num_col_match = 0;
	int *match_cols = NULL;

	int line_read;

	rlist->next = NULL;
	rlist->line = NULL;

	if (argc<3)
	{
		printf("USAGE: %s [options] <F1> <F2> > <OVERLAPS>\n", argv[0]);
		printf("Input format (F1 and F2): \n");
		printf("   [Skip Columns]Chr Start End [strand+-] [addl text] (whitespace separated)\n");
		printf("   Lines must be sorted first by Chr, then by Start\n");
		printf("   - may be used to indicate usage of stdin\n");
		printf("Output format: \n");
		printf("   L1|L2|len(R1)|len(R2)|OLen|Chr|OStart|OEnd\n");
		printf("   L1, L2 are the lines from F1 and F2, respectively\n");
		printf("   len(R1), len(R2) are the corresponding lengths\n");
		printf("   OLen, Chr, OStart, OEnd indicate the length and coordinates of the overlap (including the distance padding to S2)\n");
		printf("Options: \n");
		printf("-s char: separating character [default: |]\n");
		printf("-o: just show matching lines in F1\n");
		printf("-v: just show non-matching lines in F1\n");
		printf("-t/-T: require same/opposite strand match (note: if no strand is indicated, neither will match; -t0/-t1/-t2 or aliases for default/-t/-T)\n");
		printf("-p: print all lines in F1 with 0/1 prefix indicating not match/match\n");
		printf("-w num: consider regions that are w distance apart as overlapping\n");
		printf("-c num, -c1 num, -c2 num: number of columns to skip for before the Chr for all input, F1 and F2, respectively [default: 1]\n");
		printf("   a value of -1 specifies format of: Chr Start End [skip] [skip] [strand+-]\n");
		printf("-z, -z1, -z2: zero indexed starts for all, F1, F2, respectively\n");
		printf("-r1, -r2: require that OLen equal len(R1) or len(R2), respectively\n");
		printf("-m list: comma separated list of columns that must match exactly; negative for must not match\n");
		printf("-e1/-E1: print out lines in F1 that do not match anything/all lines with len(R2) as 0\n");

		return 1;
	}

	ASSERTE(!STREQ(argv[argc-2], "-") ||  !STREQ(argv[argc-1], "-"), "Both input file streams cannot be stdin.\n");

	fp1 = gzpopen(argv[argc-2]);
	fp2 = gzpopen(argv[argc-1]);

	for (i=1; i<argc-2; i++)
		if (STREQ(argv[i], ""))
			continue; 
		else if (STREQ(argv[i], "-o"))
			match_mode = MATCH_MODE_ONLYMATCHES;
		else if (STREQ(argv[i], "-v"))
			match_mode = MATCH_MODE_ONLYNONMATCHES;
		else if (STREQ(argv[i], "-p"))
			match_mode = MATCH_MODE_PREFIX;
		else if (STREQ(argv[i], "-e1"))
			match_mode = MATCH_MODE_NORMAL_PLUS_NO;
		else if (STREQ(argv[i], "-E1"))
			match_mode = MATCH_MODE_NORMAL_PLUS_ALL;
		else if (STREQ(argv[i], "-w"))
			padding = atoi(argv[++i]);
		else if (STREQ(argv[i], "-s"))
			sep_char = argv[++i][0];
		else if (STREQ(argv[i], "-c"))
			fp1_skip_col = fp2_skip_col = atoi(argv[++i]);
		else if (STREQ(argv[i], "-c1"))
			fp1_skip_col = atoi(argv[++i]);
		else if (STREQ(argv[i], "-c2"))
			fp2_skip_col = atoi(argv[++i]);
		else if (STREQ(argv[i], "-m"))
			num_col_match = col_match_parse(&match_cols, argv[++i]);
		else if (STREQ(argv[i], "-z"))
			start1_zeroidx = start2_zeroidx = 1; 
		else if (STREQ(argv[i], "-z1"))
			start1_zeroidx = 1; 
		else if (STREQ(argv[i], "-z2"))
			start2_zeroidx = 1; 
		else if (STREQ(argv[i], "-r1"))
			match_r1_olen = 1; 
		else if (STREQ(argv[i], "-r2"))
			match_r2_olen = 1; 
		else if (STREQ(argv[i], "-t0"))
			require_str = 0;
		else if (STREQ(argv[i], "-t") || STREQ(argv[i], "-t1"))
			require_str = 1; 
		else if (STREQ(argv[i], "-T") || STREQ(argv[i], "-t2"))
			require_str = 2; 
		else
			ASSERTE(0, "Invalid command line argument.\n");

	ASSERTE(padding >= 0, "-w argument must be >= 0.\n");

	/* fill out output format.. notice we allocate plenty of space for out_fmt above */
	if (match_mode == MATCH_MODE_NORMAL || match_mode == MATCH_MODE_NORMAL_PLUS_NO || match_mode == MATCH_MODE_NORMAL_PLUS_ALL)
	{
		strcpy(out_fmt, "%s %s %Ld %Ld %Ld %s %Ld %Ld\n");

		for (i=0; out_fmt[i] != '\n'; i++)
			if (out_fmt[i] == ' ')
				out_fmt[i] = sep_char;
	}
	else if (match_mode == MATCH_MODE_PREFIX)
	{
		strcpy(out_fmt, "%d %s\n");
		out_fmt[2] = sep_char;
	}
	else
		strcpy(out_fmt, "%s\n");

	/* fill in the formats for input*/
	fp1_fmt = fill_format(fp1_skip_col);
	fp2_fmt = fill_format(fp2_skip_col);

	while ((line_read = getline(&rlist->line, &rlist_line_alloc, fp1)) != -1)
	{
		/* used for -v and -o modes of operation */
		int has_match = 0;

		/* used for iterating through rlist */
		struct region* rlist_ptr = NULL;

		rlist->line[--line_read] = '\0';

		FILL_REGION(rlist, fp1_fmt);

		if (start1_zeroidx)
			rlist->start++;

		/* make sure that the data point makes sense */
		ASSERTE(rlist->start <= rlist->end, "%s:%s has start follow end.\n", argv[argc-2], rlist->line);

		/* check to make sure fp1 is sorted */
		ASSERTE(ISINORDER(fp1_last_chr, fp1_last_start, rlist->chr, rlist->start), "%s:%s is out of order.\n", argv[argc-2], rlist->line);
		strcpy(fp1_last_chr, rlist->chr);
		fp1_last_start = rlist->start;

		/* scan through the list, eliminate unneeded entries */
		for (rlist_ptr=rlist; rlist_ptr->next != NULL; )
			if (ISINORDER(rlist_ptr->next->chr, rlist_ptr->next->end+1, rlist->chr, rlist->start))
			{
				/* this item will never be overlapped with... delete it */
				struct region* temp = rlist_ptr->next->next;

				if (rlist_ptr->next == rlist_end)
					rlist_end = rlist_ptr;

				free(rlist_ptr->next->line);
				free(rlist_ptr->next);
				rlist_ptr->next = temp;
			}
			else if (ISINORDER(rlist->chr, rlist->start, rlist_ptr->next->chr, rlist_ptr->next->start))
				/* a short circuit to stop after the stack item's start is after the current start; 
				 * we won't be deleting any more */
				break;
			else
				rlist_ptr=rlist_ptr->next;

		/* if fp2 is empty and rlist->next is empty, we are done (unless we are doing something with non-matches) */
		if (fp2_empty && rlist->next == NULL && match_mode != MATCH_MODE_ONLYNONMATCHES && match_mode != MATCH_MODE_PREFIX && match_mode != MATCH_MODE_NORMAL_PLUS_NO && match_mode != MATCH_MODE_NORMAL_PLUS_ALL)
			break;

		/* read in new regions until a region that starts after the end of this region *
		 * is read in */
		while (
				/* if this file is empty, we cannot scan for any more new items */
				!fp2_empty && (
				
				/* if the list is empty... read in more */
				rlist->next == NULL || 
				
				/* or if what we read in from F1 ends after our stack (i.e. the
				 * last stack item starts before one after the item we just read in)
				 * we want them out of order so we can ensure we have nothing left
				 * in fp2 that will overlap with rlist */
				ISINORDER(rlist_end->chr, rlist_end->start, rlist->chr, rlist->end)
				)
		)
		{
			size_t temp_alloc = 0;
			rlist_end->next = (struct region*) malloc (sizeof(struct region));

			rlist_end->next->line = NULL;
			rlist_end->next->next = NULL;
			
			if ((line_read = getline (&rlist_end->next->line, &temp_alloc, fp2)) == -1)
			{
				fp2_empty = 1;
				free(rlist_end->next->line);
				free(rlist_end->next);
				rlist_end->next = NULL;
			}
			else
			{
				rlist_end->next->line[--line_read] = '\0';

				FILL_REGION(rlist_end->next, fp2_fmt);

				if (start2_zeroidx)
					rlist_end->next->start++;

				/* make sure that the data point makes sense */
				ASSERTE(rlist_end->next->start <= rlist_end->next->end, "%s:%s has start follow end.\n", argv[argc-1], rlist_end->next->line);

				/* check to make sure fp2 is sorted */
				ASSERTE(ISINORDER(fp2_last_chr, fp2_last_start, rlist_end->next->chr, rlist_end->next->start), "%s:%s is out of order.\n", argv[argc-1], rlist_end->next->line);
				strcpy(fp2_last_chr, rlist_end->next->chr);
				fp2_last_start = rlist_end->next->start;

				/* add padding for the padding */
				rlist_end->next->start -= padding;
				rlist_end->next->end += padding;

				/* throw out if this entry is already behind us */
				if (ISINORDER(rlist_end->next->chr, rlist_end->next->end+1, rlist->chr, rlist->start))
				{
					free(rlist_end->next->line);
					free(rlist_end->next);
					rlist_end->next = NULL;
				}
				else
					rlist_end = rlist_end->next;
			}
		}

		if (match_mode == MATCH_MODE_NORMAL_PLUS_ALL)
			printf(out_fmt, rlist->line, "", rlist->end - rlist->start + 1, (long long int) 0, (long long int) 0, rlist->chr, (long long int) 0, (long long int) 0);

		/* scan through the list and print overlaps */
		for (rlist_ptr=rlist; rlist_ptr->next != NULL; rlist_ptr=rlist_ptr->next)
			if (ISINORDER(rlist->chr, rlist->end+1, rlist_ptr->next->chr, rlist_ptr->next->start))
				/* a short circuit to stop after the stack is past the end of the current end; also
				 * needed to make sure we are in the right chromosome! */
				break;
			else
			{
				long long int overlap_size = OVERLAP(rlist_ptr->next->start, rlist_ptr->next->end, rlist->start, rlist->end);
				if (overlap_size > 0
					&& (!match_r1_olen || overlap_size == (rlist->end - rlist->start + 1))
					&& (!match_r2_olen || overlap_size == (rlist_ptr->next->end - rlist_ptr->next->start + 1 - 2 * padding))
					&& (require_str==0 || (rlist->str && rlist_ptr->next->str && ((require_str == 1) == (rlist->str == rlist_ptr->next->str))))
					&& (num_col_match==0 || check_col_match(rlist->line, rlist_ptr->next->line, num_col_match, match_cols))
				)
				{
					has_match = 1;

					if (match_mode == MATCH_MODE_NORMAL || match_mode == MATCH_MODE_NORMAL_PLUS_NO || match_mode == MATCH_MODE_NORMAL_PLUS_ALL)
						printf(out_fmt, rlist->line, rlist_ptr->next->line, rlist->end - rlist->start + 1, rlist_ptr->next->end - rlist_ptr->next->start + 1 - 2 * padding, overlap_size, rlist->chr, MAX(rlist_ptr->next->start,rlist->start), MIN(rlist_ptr->next->end,rlist->end));
					else
						/* we can short circuit in match_mode 1, 2 and 3 */
						break;
				}
			}

		if ((match_mode == MATCH_MODE_ONLYMATCHES && has_match) || (match_mode == MATCH_MODE_ONLYNONMATCHES && !has_match))
			printf(out_fmt, rlist->line);
		else if (match_mode == MATCH_MODE_PREFIX)
			printf(out_fmt, has_match ? 1 : 0, rlist->line);
		else if (match_mode == MATCH_MODE_NORMAL_PLUS_NO && !has_match)
			printf(out_fmt, rlist->line, "", rlist->end - rlist->start + 1, (long long int) 0, (long long int) 0, rlist->chr, (long long int) 0, (long long int) 0);
	}

	free(fp1_fmt);
	free(fp2_fmt);
	free(match_cols);

	pclose(fp1);
	pclose(fp2);
	while (rlist != NULL)
	{
		rlist_end = rlist;
		rlist = rlist->next;
		free(rlist_end->line);
		free(rlist_end);
	}

	return 0;
}

/* when not doing strand, the strand may be left unfilled.. we don't care about that */
char* fill_format (int num_skip)
{
	char* fmt; 
	int i;
	
	if (num_skip == -1)
	{
		/* bed format... adding 1 to start has to happen elsewhere */
		fmt = (char*) malloc(sizeof(char) * 22);
		strcpy(fmt, "%s %Ld %Ld %*s %*s %c");
	}
	else
	{
		fmt = (char*) malloc(sizeof(char)*(4*num_skip+14));
		for (i=0; i<(4*num_skip); i+=4)
			strcpy(fmt+i, "%*s ");
		strcpy(fmt+4*num_skip, "%s %Ld %Ld %c");
	}
	
	return fmt;
}


/* require columns indicated in match_cols match between s1 and s2
 * note: if column does not exist, then check_col_match returns false */
int check_col_match (char* s1, char* s2, int num_col_match, int* match_cols)
{
	int i, c, cn;

	for (c=1, cn=0; cn<num_col_match; cn++)
	{
#define ISWHITE(x) ((x) == ' ' || (x) == '\t' || (x) == '\0')
		for (;c < abs(match_cols[cn]); c++)
		{
			for (; !ISWHITE(*s1); s1++);
			for (; !ISWHITE(*s2); s2++);

			if (*s1 == '\0' || *s2 == '\0')
				return 0;

			s1++;
			s2++;
		}

		for (i=0; !ISWHITE(s1[i]) && !ISWHITE(s2[i]) && s1[i] == s2[i]; i++);

		if ((ISWHITE(s1[i]) && ISWHITE(s2[i])) != (match_cols[cn] > 0))
			return 0;
#undef ISWHITE
	}

	return 1;
}


int abs_int_cmp (const void *a, const void *b)
{
	if (abs(*((int*) a)) > abs(*((int*) b)))
		return 1;
	if (abs(*((int*) a)) < abs(*((int*) b)))
		return -1;

	return 0;
}

int col_match_parse (int **match_cols_ptr, char *str)
{
	char *temp = strdup(str);
	int n = 0;

	free (*match_cols_ptr);
	*match_cols_ptr = NULL;

	if (temp[0] != '\0')
	{
		int i, j;
		n = 1;
		for (i=0; temp[i] != '\0'; i++)
			if (temp[i] == ',')
				n++;
		*match_cols_ptr = (int *) malloc(sizeof(int) * n);

		for (i=0, j=0, n=0; ; i++)
			if (temp[i] == ',' || temp[i] == '\0')
			{
				if (temp[i] == '\0')
				{
					(*match_cols_ptr)[n++] = atoi(temp + j);
					break;
				}
				else
				{
					temp[i] = '\0';
					(*match_cols_ptr)[n++] = atoi(temp + j);
				}
				j = i + 1;
			}
		qsort(*match_cols_ptr, n, sizeof(int), abs_int_cmp); 
	}
	
	free(temp);
	return n;
}

