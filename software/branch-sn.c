/* Written by Pouya Kheradpour (pouyak@gmail.com)
 * Copyright (c) 2005-2015 Massachusetts Institute of Technology
 * Released under the MIT license
 * http://compbio.mit.edu/pouyak/software/LICENSE.mit */
/* version 20160730 */

/* Accompanying code for Kheradpour, Stark, et al. (2007):
 * Reliable prediction of regulator targets using 12 Drosophila genomes
 *
 * Software not to be redistributed without permission from authors
 * Support contact email: manoli at mit.edu
 *
 * PURPOSE:
 *   1) calculates branch length score (BLS) and its confidence
 *   2) identifies potential control motifs
 *   3) counts the number of instances at each confidence level
 *
 * COMPILATION:
 *   gcc -Wall -O3 -m32 -lm -o branch-sn branch-sn.c
 *
 * SAMPLE USAGE:
 *
 * REQUIRED FILES (. in formats indicates spacer column):
 *
 * <INPUTFILE> with lines of format:
 *    motif . . . . species_bitstring
 *    where species_bitstring is, for e.g. 5 if both the first and
 *    third informant species species match, but not the second
 *    these represent each match to each motif and potential
 *    control motif in the desired regions
 *
 * <MOTIFFILE> with one motif (no whitespace) per line (should not include control motifs)
 *
 * <BLFILE> with lines of format:
 *    species_bitstring . . branch_length
 *
 * <NUMSP> is the number of species includingthe target species
 *    thus, 2^(<NUMSP>-1)-1 is the greatest species_bitstring
 *
 * First, compute <CONTORLFILE> (control motifs for each motif):
 *    branch-sn -b <BLFILE> -n <NUMSP> -e <MOTIFFILE> -pc 1 <INPUTFILE> | grep -w -f <MOTIFFILE> > <CONTROLFILE>
 *    NOTE: control motifs have the same composition and are within +/-20% of the number of instances
 *    NOTE: we do not permit one of the non-control motifs to be control motifs
 *    NOTE: this only works if the motif in <INPUTFILE> is a consensus sequence
 *
 * Next, use these and compute the branch length to confidence mapping:
 *    branch-sn -f <CONTROLFILE> -b <BLFILE> -n <NUMSP> <INPUTFILE> > <OUTFILE>
 *
 * <OUTFILE> has the format:
 *    column 1:      motif
 *    columns 2-103: for branch lengths MAXBL*(0.00, 0.01, ..., 0.99, 1.00):
 *                   confidence;#instances;branch_length
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <search.h>

#include "pk.h"
#include "pk-bl.h"

#define MAXMOTLEN			200
#define ALPHSIZE			26
#define MAXWSTRLEN			3000

/* indexing for motif->bin_cnts (if w == 0 then can just directly index with bin)*/
#define BCidx(w,bin)		((num_bins) * (w) + (bin))

struct motif
{
	char str[MAXMOTLEN+1];
	unsigned char base_cnts[ALPHSIZE];
	unsigned int *bin_cnts;
	unsigned int total_cnt;
	unsigned char c_ex; /* indicates excluded from being a control */
};

/* compares motifs on the basis of base_cnts (i.e. motifs with the same composition will
 * have return value 0) */
int motifs_cmp (const void *a, const void *b)
{
	return memcmp(((struct motif*)a)->base_cnts, ((struct motif*)b)->base_cnts, sizeof(unsigned char)*ALPHSIZE);
}

/* read in all the motifs along with their counts at each of the blcuts
 * final motifs table is sorted such that motifs with the same composition are together */
int read_motifs (char* fn, struct motif **motifs_ptr, unsigned char num_bins, unsigned char num_sp, long long int sp_mask, struct consm *consm);

/* same as read_motifs except uses output from motif_match directly without use of mm2branch-sn */
int read_motifs_mm (char* fn, struct motif **motifs_ptr, unsigned char num_bins, unsigned char num_sp, long long int sp_mask, struct consm *consm, int num_windows, unsigned int *windows, unsigned char motif_match_input);

/* !!!NOTE!!!: must make sure that n_c and nc != 0 */
double count_to_confid (double corr_z, unsigned char corr_real, unsigned int n, unsigned int nc, unsigned int n_c, unsigned int nc_c);

void print_motif_line (char* mstr, unsigned int*  m_bins, unsigned int* c_bins, float sp_mask_bl, unsigned char num_bins, double corr_z, unsigned char corr_real, FILE* ofp);

unsigned char spbit_to_bin(long long int spbit, unsigned char num_sp, long long int sp_mask, float sp_mask_bl, struct consm *consm, unsigned char num_bins);

int main (int argc, char** argv)
{
	unsigned char num_sp = 0;
	int num_motifs = 0;
	unsigned char num_bins = 101;
	struct motif* motifs;
	char* bl_file = NULL;
	char* blp_file = NULL;
	char* c_file = NULL;
	char* ex_file = NULL;
	float tol_diff = 0.2;
	long long int sp_mask = -1;
	int print_controls = 0;
	unsigned char motif_match_input = 1;
	double corr_z = 0;
	struct consm consm = {NULL, NULL, NULL, NULL, NULL};
	float sp_mask_bl;
	unsigned int* windows = NULL;
	int num_windows = 1;
	char* outfile_fmt = NULL;
	unsigned char corr_real = 0;

	long int i;
	int j, w;

	for (i=1; i<argc-1; i++)
		if (STREQ(argv[i], ""))
			continue;
		else if (STREQ(argv[i], "-n"))
			num_sp = atoi(argv[++i]);
		else if (STREQ(argv[i], "-m"))
		{
			num_bins = atoi(argv[++i]);
			ASSERTE(num_bins == atoi(argv[i]), "Number of bins must be less than %d.\n", 1 << (8*sizeof(unsigned char)));
		}
		else if (STREQ(argv[i], "-k"))
			sp_mask = atoi(argv[++i]);
		else if (STREQ(argv[i], "-t"))
			tol_diff = atof(argv[++i]);
		else if (STREQ(argv[i], "-e"))
			ex_file = argv[++i];
		else if (STREQ(argv[i], "-f"))
			c_file = argv[++i];
		else if (STREQ(argv[i], "-b"))
			bl_file = argv[++i];
		else if (STREQ(argv[i], "-bp"))
			blp_file = argv[++i];
		else if (STREQ(argv[i], "-o"))
			outfile_fmt = argv[++i];
		else if (STREQ(argv[i], "-w"))
		{
			++i;
			int len = strlen(argv[i]);

			num_windows = 1;
			for (j=0; j<len; j++)
				if (argv[i][j] == ' ')
					num_windows++;

			if (windows)
				free(windows);

			if (len == 0)
				windows = NULL;
			else
			{
				windows = (unsigned int*) malloc(sizeof(unsigned int)*num_windows);
				num_windows = 1;
				windows[0] = atoi(argv[i]);
				for (j=0; j<len; j++)
					if (argv[i][j] == ' ')
						windows[num_windows++] = atoi(argv[i] + j + 1);
			}
		}
		else if (STREQ(argv[i], "-z"))
			corr_z = atof(argv[++i]);
		else if (STREQ(argv[i], "-C"))
			corr_real = atoi(argv[++i]) != 0;
		else if (STREQ(argv[i], "-mm"))
			motif_match_input = atoi(argv[++i]);
		else if (STREQ(argv[i], "-pc"))
			print_controls = atoi(argv[++i]) != 0;
		else
			ASSERTE(0, "Invalid command line argument.\n");

	if (num_sp == 0 || (bl_file == NULL && blp_file == NULL))
	{
		printf("USAGE: %s <-n int> <-b file> [-t frac] [-n num] [-f file] [-e file] [-pc 0/1] [-mm 0/1/2] [-z float] [-C 0/1] <input file>\n", argv[0]);
		printf("    -n    Number of species (0 <= #Bls < 2^(NS - 1)) [required]\n");
		printf("    -f    List of controls for each motif (overrides -e, -pc, -t)\n");
		printf("    -e    List of motifs to exclude from controls\n");
		printf("    -k    Species mask (-1 for all) [default: -1]\n");
		printf("    -b    Branch length file [either -b or -bp required]\n");
		printf("    -bp   Branch length file in parent tree format [either -b or -bp required].\n");
		printf("    -m    Number of bins to use for output [default: 101]\n");
		printf("    -t    Tolerance permitted between motifs and matched controls [default: 0.2]\n");
		printf("    -z    Z-value to use to create confidence intervals for ratios when generating confidence value [default: 0]\n");
		printf("    -C    Use real motif correction when computing confidence value [default: 0]\n");
		printf("    -o    If more than one window is included, output file format for sprintfu [required with more than one -w window]\n");
		printf("    -mm   0: Input from mm2branch-sn, 1: input from motif-match, 2: motif and species columns from motif-match, 3: input from motif-match -v >= 0, 4: same as 3 but with only motif and species columns [default: 1]\n");
		printf("    -w    Window to use (-mm must be 1-4) (space separated list permitted) \n");
		printf("    -pc   Prints control motifs chosen for each motif [default: 0]\n");
		return 1;
	}

	ASSERTE(num_windows == 1 || outfile_fmt != NULL, "-o is required when -w has more than one value\n");
	ASSERTE(windows != NULL || outfile_fmt == NULL, "-w is required with -o\n");
	ASSERTE(windows != NULL || motif_match_input <= 2, "-w is required with -mm > 2\n");

	if (sp_mask < 0)
		sp_mask = (1 << (num_sp-1)) - 1;

	/* read in appropriate branch length files */
	if (bl_file != NULL)
		consm.bl = read_bl_file(bl_file, num_sp);
	else if (blp_file != NULL)
		read_blp_file(blp_file, num_sp, &consm);

	sp_mask_bl = get_bl (num_sp, &consm, sp_mask);

	/* create motifs array */
	if (!motif_match_input)
		num_motifs = read_motifs(argv[argc-1], &motifs, num_bins, num_sp, sp_mask, &consm);
	else
		num_motifs = read_motifs_mm(argv[argc-1], &motifs, num_bins, num_sp, sp_mask, &consm, num_windows, windows, motif_match_input);

	/* create hash table for motifs for quick lookup */
	hcreate (num_motifs*2);
	for (i=0; i<num_motifs; i++)
	{
		ENTRY e;
		e.key = motifs[i].str;
		e.data = (void*) i;
		hsearch(e,ENTER);
	}

	for (w=0; w<num_windows; w++)
	{
		FILE* ofp;
		if (outfile_fmt == NULL)
			ofp = stdout;
		else
		{
			char outfile[MAXFILENAMELEN+1];
			sprintf(outfile,outfile_fmt,windows[w]);
			ofp = fopen(outfile, "w");
		}

		if (c_file != NULL)
		{
			/* read controls explicitly from c_file */
			FILE* fp = fopen(c_file, "r");
			ENTRY e, *ep;
			size_t len = 0;
			char *line = NULL;
			ssize_t read;

			/* read in a line from the file */
			while ((read = getline(&line, &len, fp)) != -1)
			{
				/* get rid of \n */
				if (line[read-1] == '\n')
					line[read-1] = '\0';

				/* look up the location of the main motif from the hash file */
				e.key = strtok(line, "\t ");
				ep = hsearch(e, FIND);

				if (ep != NULL)
				{
					/* if found, begin */
					struct motif* m = &motifs[(long int) ep->data];
					unsigned int *c_bins = (unsigned int*) calloc(num_bins,sizeof(unsigned int));

					if (print_controls)
						fprintf(ofp, "%s", m->str);

					e.key = strtok(NULL, "\t ");
					while (e.key != NULL)
					{
						/* if found, update c_bins */
						ep = hsearch(e, FIND);
						if (ep != NULL)
						{
							for (i=0; i<num_bins; i++)
								c_bins[i] += motifs[(long int) ep->data].bin_cnts[BCidx(w,i)];
							if (print_controls && m->total_cnt*(1-tol_diff) <= motifs[(long int) ep->data].total_cnt && m->total_cnt*(1+tol_diff) >= motifs[(long int) ep->data].total_cnt)
								fprintf(ofp, "\t%s", motifs[(long int) ep->data].str);
						}
						e.key = strtok(NULL, "\t ");

					}

					/* print the line corresponding to this motif */
					if (print_controls)
						fprintf(ofp, "\n");
					else
						print_motif_line(m->str, m->bin_cnts + num_bins*w, c_bins, sp_mask_bl, num_bins, corr_z, corr_real, ofp);
					free(c_bins);
				}
			}

			fclose(fp);
			if (line)
				free(line);
		}
		else
		{
			/* same_end is actually one AFTER last with same composition */
			int same_start = 0, same_end = 0;

			/* exclude motifs */
			if (ex_file != NULL)
			{
				ENTRY e, *ep;
				FILE* fp = fopen (ex_file, "r");
				char ex_motif[MAXMOTLEN+1];
				e.key = ex_motif;

				while (fscanf(fp, "%s", e.key) != EOF)
				{
					ep = hsearch(e, FIND);
					if (ep != NULL)
						motifs[(long int) ep->data].c_ex = 1;
				}
				fclose(fp);
			}

			/* scan through motifs and find matching instances */
			for (i=0; i<num_motifs; i++)
			{
				/* same_start and same_end correspond to the range of motifs with the same
				 * composition as the current one */
				/* same_end is actually one AFTER last with same composition */
				if (i >= same_end)
				{
					same_start = i;

					/* correct same_end */
					for (same_end=i+1; same_end<num_motifs && motifs_cmp(&motifs[i], &motifs[same_end]) == 0; same_end++);
				}

				if (!print_controls)
				{
					int k;
					/* print confidence at each bl cut */
					unsigned int *c_bins = (unsigned int*) calloc(num_bins,sizeof(unsigned int));

					for (j=same_start; j<same_end; j++)
						if (j!=i && !motifs[j].c_ex
								 && motifs[i].total_cnt*(1-tol_diff) <= motifs[j].total_cnt
								 && motifs[i].total_cnt*(1+tol_diff) >= motifs[j].total_cnt
						)
							for (k=0; k<num_bins; k++)
								c_bins[k] += motifs[j].bin_cnts[BCidx(w,k)];
					print_motif_line(motifs[i].str, motifs[i].bin_cnts + num_bins * w, c_bins, sp_mask_bl, num_bins, corr_z, corr_real, ofp);

					free(c_bins);
				}
				else
				{
					fprintf(ofp,"%s", motifs[i].str);
					for (j=same_start; j<same_end; j++)
						if (j!=i && !motifs[j].c_ex
								 && motifs[i].total_cnt*(1-tol_diff) <= motifs[j].total_cnt
								 && motifs[i].total_cnt*(1+tol_diff) >= motifs[j].total_cnt
						)
							fprintf(ofp,"\t%s", motifs[j].str);
					fprintf(ofp,"\n");
				}
			}
		}

		fclose(ofp);
	}

	/* free up used memory */
	for (i=0; i<num_motifs; i++)
		free(motifs[i].bin_cnts);
	free(windows);
	free(motifs);

	free(consm.mat);
	free(consm.tlist);
	free(consm.len);
	free(consm.bl);
	return 0;
}

unsigned int* cumul_array (unsigned int* s, int n)
{
	unsigned int* t = (unsigned int*) malloc (n * sizeof(unsigned int));

	int i;
	t[n-1] = s[n-1];
	for (i=n-1; i>0; i--)
		t[i-1] = t[i] + s[i-1];

	return t;
}

void print_motif_line (char* mstr, unsigned int*  m_bins, unsigned int* c_bins, float sp_mask_bl, unsigned char num_bins, double corr_z, unsigned char corr_real, FILE* ofp)
{
	unsigned char j;

	unsigned int* m_bins_c = cumul_array(m_bins, num_bins);
	unsigned int* c_bins_c = cumul_array(c_bins, num_bins);

	fprintf(ofp, "%s", mstr);
	for (j=0; j<num_bins; j++)
	{
		fprintf(ofp, "\t");

		if (c_bins_c[0] > 0 && m_bins_c[j] > 0)
			fprintf(ofp, "%lf", count_to_confid(corr_z, corr_real, m_bins_c[0], m_bins_c[j], c_bins_c[0], c_bins_c[j]));
		else
			/* use -1 confidence to indicate unknown */
			fprintf(ofp, "-");

		fprintf(ofp, ";%d;%f;%d", m_bins_c[j], ((float)j*sp_mask_bl)/(num_bins-1), c_bins_c[j]);
	}
	fprintf(ofp, "\n");

	free(m_bins_c);
	free(c_bins_c);
}

double binconfid(double Z, unsigned int n, unsigned int na)
{
	if (n == 0) /* NOTE: this is just to prevent overflow... actually binconfid is undefined! */
		return 0.0;

	double p = na / ((double) n);

	return (p+((Z*Z)/(2.0*n))+Z*sqrt(((p*(1.0-p))/n)+((Z*Z)/(4.0*n*n))))/(1.0+Z*Z/n);
}

double count_to_confid (double corr_z, unsigned char corr_real, unsigned int n, unsigned int nc, unsigned int n_c, unsigned int nc_c)
{
	/* NOTE: we no longer attempt return negative confidences!
	 * this is important because if r_c == 1, then divide by zero if corr_real */

	double r = binconfid(-corr_z, n, nc); /* correct cons rate of motif down */
	double r_c = binconfid(corr_z, n_c, nc_c); /* correct cons rate of controls up */

	if (r_c >= r)
		return 0.0;

	if (corr_real)
		return 1.0 - ((r_c/r) * ((1.0 - r)/(1.0 - r_c)));
	else
		return 1.0 - (r_c/r);
}

int read_motifs (char* fn, struct motif **motifs_ptr, unsigned char num_bins, unsigned char num_sp, long long int sp_mask, struct consm *consm)
{
	int num_motifs = 0;
	int size_reserved = 1 << 5;
	struct motif* motifs = (struct motif*) malloc(size_reserved*sizeof(struct motif));
	long long int num_bl = 1 << (num_sp-1);
	float sp_mask_bl = get_bl (num_sp, consm, sp_mask);

	FILE* fp = STREQ(fn, "-") ? stdin : fopen(fn, "r");
	ASSERTE(fp != NULL, "File could not be opened!\n");

	while (fscanf(fp, "%s", motifs[num_motifs].str) != EOF)
	{
		struct motif* m;
		long long int i;

		num_motifs++;

		if ((size_reserved-1) < num_motifs)
		{
			size_reserved *= 2;
			motifs = (struct motif*) realloc(motifs,size_reserved*sizeof(struct motif));
		}

		m = &motifs[num_motifs-1];

		m->c_ex = 0;
		memset(m->base_cnts, 0, sizeof(unsigned char)*ALPHSIZE);
		m->bin_cnts = (unsigned int*) calloc(num_bins, sizeof(unsigned int));
		m->total_cnt = 0;

		for (i=0; m->str[i]; i++)
		{
			/* note; with -f all this base composition stuff is ignored allowing for motifs
			 * to have names that include non-alpha characters */
			int c = toupper(m->str[i]) - 'A';
			if (c >= 0 && c < 26)
				m->base_cnts[c]++;
		}

		for (i=0; i<num_bl; i++)
		{
			unsigned int c;
			fscanf(fp, "%u", &c);
			m->bin_cnts[spbit_to_bin(i, num_sp, sp_mask, sp_mask_bl, consm, num_bins)] += c;
			m->total_cnt += c;
		}
	}

	/* sort motifs so that motifs with the same composition are put together */
	qsort (motifs, num_motifs, sizeof(struct motif), motifs_cmp);

	*motifs_ptr = motifs;

	fclose(fp);

	return num_motifs;
}

struct motif_idx
{
	char str[MAXMOTLEN+1];
	int idx;
};

int motif_idx_cmp (const void *a, const void *b)
{
	return strcmp(((struct motif_idx*)a)->str, ((struct motif_idx*)b)->str);
}

int read_motifs_mm (char* fn, struct motif **motifs_ptr, unsigned char num_bins, unsigned char num_sp, long long int sp_mask, struct consm *consm, int num_windows, unsigned int *windows, unsigned char motif_match_input)
{
	int num_motifs = 0, w;
	int size_reserved = 1 << 5;
	struct motif* motifs = (struct motif*) malloc(size_reserved*sizeof(struct motif));
	struct motif_idx m, *mp, **mpp;
	void* motif_idx_tree = NULL;
	long long int *bl;
	float sp_mask_bl = get_bl (num_sp, consm, sp_mask);

	FILE* fp = STREQ(fn, "-") ? stdin : fopen(fn, "r");
	ASSERTE(fp != NULL, "File could not be opened!\n");

	char *in_fmt = windows==NULL ? (motif_match_input==1 ? "%s %*s %*s %*s %*s %d%*[^\n]\n" : "%s %d\n") : (motif_match_input==1 ? "%s %*s %*s %*s %*s %*d %s%*[^\n]\n" : (motif_match_input==2 ? "%s %s\n" : (motif_match_input==3 ? "%s %*s %*s %*s %*s %*d %*[^;];%s%*[^\n]\n" : "%s %*[^;];%s\n")));

	bl = (long long int*) malloc(sizeof(long long int)*num_windows);

	while (1)
	{
		if (windows == NULL)
		{
			if (fscanf(fp, in_fmt, m.str, bl) == EOF)
				break;
		}
		else
		{
			char wstr[MAXWSTRLEN+1], *wstr_ptr;
			int sp = 0;
			if (fscanf(fp, in_fmt, m.str, wstr) == EOF)
				break;

			memset(bl, 0, sizeof(long long int)*num_windows);
			for (wstr_ptr=wstr, sp=0; *wstr_ptr != '\0'; sp++)
			{
				if (*wstr_ptr != ';')
				{
					unsigned int wind = abs(atoi(wstr_ptr));
					for (w=0; w<num_windows; w++)
						if (wind <= windows[w])
							bl[w] |= 1 << sp;
				}

				/* notice this works for when the motif-match parameter -v is -2 or >= 0, because it skips the extra :.:.* */
				for (; *wstr_ptr != ';' && *wstr_ptr != '\0'; wstr_ptr++);
				if (*wstr_ptr == ';')
					wstr_ptr++;
			}
		}

		/* lookup in tree */
		mpp = tfind((void *) (&m), &motif_idx_tree, motif_idx_cmp);

		/* if not found, setup the motif */
		if (mpp == NULL)
		{
			int i;

			/* add to array */
			num_motifs++;
			if ((size_reserved-1) < num_motifs)
			{
				size_reserved *= 2;
				motifs = (struct motif*) realloc(motifs,size_reserved*sizeof(struct motif));
			}

			strcpy(motifs[num_motifs-1].str, m.str);
			memset(motifs[num_motifs-1].base_cnts, 0, sizeof(unsigned char)*ALPHSIZE);
			motifs[num_motifs-1].c_ex = 0;
			motifs[num_motifs-1].bin_cnts = (unsigned int*) calloc(num_bins*num_windows, sizeof(unsigned int));
			motifs[num_motifs-1].total_cnt = 0;

			/* figure out base_cnts */
			for (i=0; motifs[num_motifs-1].str[i]; i++)
			{
				/* note; with -f all this base composition stuff is ignored allowing for motifs
				 * to have names that include non-alpha characters */
				int c = toupper(motifs[num_motifs-1].str[i]) - 'A';
				if (c >= 0 && c < 26)
					motifs[num_motifs-1].base_cnts[c]++;
			}

			/* add to tree */
			mp = (struct motif_idx*) malloc(sizeof(struct motif_idx));
			strcpy(mp->str, m.str);
			mp->idx = num_motifs-1;
			mpp = tsearch((void *) mp, &motif_idx_tree, motif_idx_cmp);
		}

		/* NOTE: we do not check bl[w], so we cannot accept motif-match files with more species */
		for (w=0; w<num_windows; w++)
			motifs[(*mpp)->idx].bin_cnts[BCidx(w,spbit_to_bin(bl[w], num_sp, sp_mask, sp_mask_bl, consm, num_bins))]++;
		motifs[(*mpp)->idx].total_cnt++;
	}

	/* NOTE: the tree is not deleted even though it could be now */

	/* sort motifs so that motifs with the same composition are put together */
	qsort (motifs, num_motifs, sizeof(struct motif), motifs_cmp);

	*motifs_ptr = motifs;

	fclose(fp);
	free(bl);

	return num_motifs;
}

unsigned char spbit_to_bin(long long int spbit, unsigned char num_sp, long long int sp_mask, float sp_mask_bl, struct consm *consm, unsigned char num_bins)
{
	return (get_bl (num_sp, consm, spbit & sp_mask) / sp_mask_bl) * (num_bins-1);
}

