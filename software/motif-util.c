/* Written by Pouya Kheradpour (pouyak@gmail.com)
 * Copyright (c) 2005-2015 Massachusetts Institute of Technology
 * Released under the MIT license
 * http://compbio.mit.edu/pouyak/software/LICENSE.mit */
/* version 20140120 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "pk.h"

#define MAXMOTLEN	200

/* subtracted from cut so that we get values that equal cut exactly */
#define EPSILON		1e-10

#define RFMidx(pos,chr,len,rev) ((rev) ? RMidx(pos,chr,len) : Midx(pos,chr))

#define RMidx(pos,chr,len) ((4*len) - (4*(pos) + (chr)) - 1)
#define Midx(pos,chr)	(4*(pos) + (chr))

#define ISSET(num,bit)	(((num) >> (bit)) & 1)

/* for converting between log-odds and frequency formats */
#define F2LG(freq, bk) (log(((freq) + (bk) * lg_pseudo) / (1.0 + lg_pseudo) / (bk)) / log(2.0))
#define LG2F(lg, bk) ((pow(2.0, (lg)) * (1 + lg_pseudo) - lg_pseudo) * (bk))

struct motif
{
	int len;
	int n;
	double pwm[MAXMOTLEN*4];
	int cluster;
	int parents[2];
	double clust_corr;
	char* name;
	double cut;
};

struct id_score
{
	int id;
	double score;
};

struct motif_comp
{
	int m1;
	int m2;

	/* indicates m2 is reversed */
	int rev;

	/* how many unpaired bases are there at the beginning of m1
	 * negative for unpaired bases of m2 */
	int offset;

	double comp;
};

/* compares in descending order */
int id_score_cmp (const void *a, const void *b)
{
	if (((struct id_score*)a)->score < ((struct id_score*)b)->score)
		return 1;
	else if (((struct id_score*)a)->score < ((struct id_score*)b)->score)
		return -1;
	return 0;
}

/* order of bases is (LSB to MSB): A C G T */
const unsigned char alpha[26] = {1,14,2,13,0,0,4,11,0,0,12,0,3,15,0,0,0,5,6,8,0,7,9,0,10,0};
const unsigned char alpha_num_on[26] = {1,3,1,3,0,0,1,3,0,0,2,0,2,4,0,0,0,2,2,1,0,3,2,0,2,0};

/* returns the number of motifs; allocates enough memory for all motifs
 * plus all the possible merges */
int read_consensus (char* fn, struct motif **motifs_ptr, char* input_fmt);

/* like read_consensus except reads matrix file */
int read_matrix (char* fn, struct motif **motifs_ptr);

/* like read_consensus except reads log-odds format motifs */
int read_logodds (char* fn, struct motif **motifs_ptr);

/* like read_consensus except reads alignment file */
int read_alignment (char* fn, struct motif **motifs_ptr);

/* like read_consensus, except read directly from input string (space separated list) */
int read_consensus_string (char* motif_str, struct motif **motifs_ptr);

/* fills in a motif_comp for motifs m1 and m2 and returns 1 if the corresponding best comp
 * is greater than cut; otherwise, returns 0*/
int comp_motifs(struct motif_comp* mc, struct motif* motifs, int m1, int m2, int strand_allow, double cut);

/* merge the motifs m1, m2 into m */
void weighted_motif_merge(struct motif* m, struct motif* m1, struct motif* m2, int rev, int offset);
void or_motif_merge(struct motif* m, struct motif* m1, struct motif* m2, int rev, int offset);

/* writes sequence of motif to str (and returns a pointer to the first non-N of str)
 * end of string is one after last non-N is made to \0 (all N's are either N or \0) */
char* motif_to_str(char* str, struct motif* motif, unsigned char trim_motif, unsigned char original_name);

/* outputs a motif and its parents */
void print_motif(struct motif* motifs, int m1);

/* prints out the palindrome score for each motif */
void run_motif_palindrome(struct motif* motifs, int num_motifs);

/* prints out the information content for each motif */
void run_motif_infoc(struct motif* motifs, int num_motifs);

/* if motifs2==NULL -> prints out mean difference for motifs in motifs1, otherwise, prints out mean difference between motifs in motifs1 and motifs2 */ 
void run_motif_meandiff(struct motif* motifs1, int num_motifs1, struct motif* motifs2, int num_motifs2, int strand_allow);

/* runs clustering of motifs and prints out output; changes motifs and num_motifs so requires pointer */
void run_motif_clust(struct motif** motifs_ptr, int *num_motifs_ptr, int strand_allow, double cut);

/* runs motif annotation and prints output */
void run_motif_anno(struct motif* motifs, int num_motifs, struct motif* motifs_anno, int num_motifs_anno, int strand_allow, double cut);

/* returns the best comparison between m1 and m2 and fills best_rev_ptr and best_offset_ptr
 * if they are not NULL
 * strand allow: 1 for fwd, 2 for rev, 3 for both */
double pearson_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr);

/* allows for up to motif_comp_param gaps; with motif_comp_param=0 same as pearson_motif_comp */
double peargap_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr);

/* computes intersection divided by union of matching sequences at best offset */
double fracmatch_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr);

/* returns the best comparison between m1 and m2 and fills best_rev_ptr and best_offset_ptr
 * if they are not NULL
 * strand allow: 1 for fwd, 2 for rev, 3 for both */
double wprod_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr);

/* returns the minimum number of positions that are different (negative) */
double ndiff_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr);

/* returns the maximum number of positions that are the same but not N */
double nsame_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr);

/* function to convert pwm column -> consensus character using manolis' correction strategy */
char correction_pwm_to_char(double *pwm_col);

/* function to convert pwm column -> consensus character where any amount of weight on a base makes it count */
char epsilon_pwm_to_char(double *pwm_col);

/* global variable passed to motif_to_str by many functions */
static unsigned char trim_motifs = 1;

/* global variable passed to motif_to_str by many functions */
static unsigned char original_name = 1;

/* global variable that indicates that in clustering the number should be printed rather than the consensus */
static unsigned char number_output = 0;

/* global variable that indicates that in clustering the output should be in tree format with () */
static unsigned char tree_output = 1;

/* global variable that indicates that in clustering the cluster correlation should be outputted */
static unsigned char correlation_output = 0;

/* global variable that indicates what output type should be used when clustering should be used for the final motif 
 * 0: motif_consensus name; 1: matrix output; 2: log-odds output */
static unsigned char output_mode = 0;

/* global variable that indicates that in all the output motifs should be reverse complimented */
static unsigned char reverse_output = 0;

/* global variable that indicates the maximum offset we allow when comparing; similar to strand_allow */
static int max_offset = MAXMOTLEN;

/* global variables with parameters for going to/from log-odds matrix */
static double lg_bg[4] = {0.25, 0.25, 0.25, 0.25};
static double lg_pseudo = 0.001;

/* global variable indicating which motif comparison function to use throughout */
typedef double(*motif_comp_func_ptr)(struct motif*, struct motif*, int, int*, double*);
motif_comp_func_ptr motif_comp_func = &pearson_motif_comp;

/* global variable indicating which function to use to convert a column of an PWM to a consensus character */
typedef char(*pwm_to_char_func_ptr)(double *);
pwm_to_char_func_ptr pwm_to_char_func = &correction_pwm_to_char;

typedef void(*motif_merge_func_ptr)(struct motif*, struct motif*, struct motif*, int, int);
motif_merge_func_ptr motif_merge_func = &weighted_motif_merge;

/* global variable; argument for wprod_motif_comp and peargap_motif_comp */
static int motif_comp_param = 0;

int main (int argc, char** argv)
{
	double cut = INFINITY;
	int strand_allow = 3;
	int num_motifs = -1;
	int num_motifs_anno = 0;

	char* mode = "clust-anno";

	struct motif *motifs;
	struct motif *motifs_anno = NULL;

	int i;

	/* all arguments take at least 1 parameter... we stop looking at parameters after argc-1 */
	for (i=1; i<argc-1; i++)
		if (STREQ(argv[i], ""))
			continue; 
		else if (STREQ(argv[i], "-c"))
			cut = atof(argv[++i]) - EPSILON;
		else if (STREQ(argv[i], "-r"))
			strand_allow = atoi(argv[++i]) * 2 | 1;
		else if (STREQ(argv[i], "-on"))
			number_output = atoi(argv[++i]) != 0;
		else if (STREQ(argv[i], "-ot"))
			tree_output = atoi(argv[++i]) != 0;
		else if (STREQ(argv[i], "-oc"))
			correlation_output = atoi(argv[++i]) != 0;
		else if (STREQ(argv[i], "-om"))
		{
			i++;
			if (STREQ(argv[i], "consensus") || STREQ(argv[i], "c") || STREQ(argv[i], "0"))
				output_mode = 0;
			else if (STREQ(argv[i], "matrix") || STREQ(argv[i], "m") || STREQ(argv[i], "1"))
				output_mode = 1;
			else if (STREQ(argv[i], "logodds") || STREQ(argv[i], "l") || STREQ(argv[i], "2"))
				output_mode = 2;
		}
		else if (STREQ(argv[i], "-or"))
			reverse_output = atoi(argv[++i]) != 0;
		else if (STREQ(argv[i], "-oo"))
			original_name = atoi(argv[++i]) != 0;
		else if (STREQ(argv[i], "-tm"))
			trim_motifs = atoi(argv[++i]) != 0;
		else if (argv[i][0] == '-' && (argv[i][1] == 'i' || argv[i][1] == 'a'))
		{
			int *num_ptr;
			struct motif** motifs_ptr;

			if (argv[i][1] == 'i')
			{
				num_ptr = &num_motifs;
				motifs_ptr = &motifs;
			}
			else
			{
				num_ptr = &num_motifs_anno;
				motifs_ptr = &motifs_anno;
			}

			if (STREQ(argv[i]+2, "") || STREQ(argv[i] + 2, "s"))
				*num_ptr = read_consensus(argv[i+1], motifs_ptr, "%s");
			else if (STREQ(argv[i]+2, "1"))
				*num_ptr = read_consensus(argv[i+1], motifs_ptr, "%s%*[^\n]\n");
			else if (STREQ(argv[i]+2, "c"))
				*num_ptr = read_consensus_string(argv[i+1], motifs_ptr);
			else if (STREQ(argv[i]+2, "n"))
				*num_ptr = read_consensus(argv[i+1], motifs_ptr, "%s %s%*[^\n]\n");
			else if (STREQ(argv[i]+2, "m"))
				*num_ptr = read_matrix(argv[i+1], motifs_ptr);
			else if (STREQ(argv[i]+2, "a"))
				*num_ptr = read_alignment(argv[i+1], motifs_ptr);
			else if (STREQ(argv[i]+2, "l"))
				*num_ptr = read_logodds(argv[i+1], motifs_ptr);
			else
				ASSERTE(0, "Unrecognized input motif type.\n");

			++i;
		}
		else if (STREQ(argv[i], "-m"))
			mode =  argv[++i];
		else if (STREQ(argv[i], "-f"))
		{
			i++;
			if (STREQ(argv[i], "wprod"))
				motif_comp_func = &wprod_motif_comp;
			else if (STREQ(argv[i], "ndiff"))
				motif_comp_func = &ndiff_motif_comp;
			else if (STREQ(argv[i], "nsame"))
				motif_comp_func = &nsame_motif_comp;
			else if (STREQ(argv[i], "pearson"))
				motif_comp_func = &pearson_motif_comp;
			else if (STREQ(argv[i], "peargap"))
				motif_comp_func = &peargap_motif_comp;
			else if (STREQ(argv[i], "fracmatch"))
				motif_comp_func = &fracmatch_motif_comp;
			else
				ASSERTE(0, "Unrecognized comparison function.\n");
		}
		else if (STREQ(argv[i], "-fp"))
		{
			i++;
			if (STREQ(argv[i], "epsilon"))
				pwm_to_char_func = &epsilon_pwm_to_char;
			else if (STREQ(argv[i], "correction"))
				pwm_to_char_func = &correction_pwm_to_char;
			else
				ASSERTE(0, "Unrecognized PWM->Consensus function.\n");
		}
		else if (STREQ(argv[i], "-fm"))
		{
			i++;
			if (STREQ(argv[i], "weighted"))
				motif_merge_func = &weighted_motif_merge;
			else if (STREQ(argv[i], "or"))
				motif_merge_func = &or_motif_merge;
			else
				ASSERTE(0, "Unrecognized motif merge function.\n");
		}
		else if (STREQ(argv[i], "-p"))
			motif_comp_param =  atoi(argv[++i]);
		else if (STREQ(argv[i], "-po"))
			max_offset =  atoi(argv[++i]);
		else if (STREQ(argv[i], "-lp"))
			lg_pseudo =  atof(argv[++i]);
		else if (STREQ(argv[i], "-la"))
			lg_bg[0] =  atof(argv[++i]);
		else if (STREQ(argv[i], "-lc"))
			lg_bg[1] =  atof(argv[++i]);
		else if (STREQ(argv[i], "-lg"))
			lg_bg[2] =  atof(argv[++i]);
		else if (STREQ(argv[i], "-lt"))
			lg_bg[3] =  atof(argv[++i]);
		else
			ASSERTE(0, "Unrecognized command-line parameter.\n");
	ASSERTE(i == argc, "Unrecognized command-line parameter.\n");

	if (num_motifs == -1)
	{
		printf("USAGE: %s [OPTIONS] \n", argv[0]);
		printf("    -i$   Input file\n");
		printf("    -a$   Annotation file\n");
		printf("      $   Input format (s: list format; m: matrix format; a: alignment format; 1: only first column; n: motif then name; c: list on command line; l: log-odds) [default: s] \n");
		printf("    -c    Cutoff for merges and annotations (inf uses less memory for clust) [default: inf]\n");
		printf("    -f    Comparison function (pearson/wprod/ndiff/nsame/peargap/fracmatch) [default: pearson]\n");
		printf("    -fp   PWM->Consensus conversion function (correction/epsilon) [default: correction]\n");
		printf("    -fm   PWM Merging function (weighted/or) [default: weighted]\n");
		printf("    -p    Parameter for comparison functions (window for wprod, max number of gaps for peargap) [default: 0]\n");
		printf("    -po   Maximum offset permitted by comparison functions\n");
		printf("    -r    Allow reversal of motifs for clustering and annotation [default: 1]\n");
		printf("    -m    Indicates mode to run in (clust-anno/palindrome/meandiff/infoc/anno/clust) [default: clust-anno]\n");
		printf("    -om   Output cluster mode (consensus/matrix/logodds) [default: consensus]\n");
		printf("    -on   Output cluster members in number format rather than consensus [default: 0]\n");
		printf("    -oo   When possible, output original motif names [default: 1]\n");
		printf("    -oc   Output cluster with the correlation [default: 0]\n");
		printf("    -or   Output motifs as reverse compliments [default: 0]\n");
		printf("    -ot   Output clusters as a tree [default: 1]\n");
		printf("    -tm   Trim motifs, removing flanking Ns [default: 1]\n");
		printf("    -l$   Log odds parameters; $=p,a,c,g,t for pseudocount and base frequencies [default: 0.001, 0.25, 0.25, 0.25, 0.25]\n");
		printf("          can be specified multiple times for different values for -ao/-io or output (must precede -io/-ao)\n");
		return 1;
	}

	{
		double s = lg_bg[0] + lg_bg[1] + lg_bg[2] + lg_bg[3];
		for (i=0; i<4; i++)
			lg_bg[i] = lg_bg[i] / s;
	}

	if (STREQ(mode, "palindrome"))
		run_motif_palindrome(motifs, num_motifs);
	else if (STREQ(mode, "infoc"))
		run_motif_infoc(motifs, num_motifs);
	else if (STREQ(mode, "meandiff"))
		run_motif_meandiff(motifs, num_motifs, motifs_anno, num_motifs_anno, strand_allow);
	else if ((STREQ(mode, "clust-anno") && motifs_anno == NULL) || STREQ(mode, "clust"))
		/* clust forces clustering, even if there is an annotation file */
		run_motif_clust(&motifs, &num_motifs, strand_allow, cut);
	else if ((STREQ(mode, "clust-anno") && motifs_anno != NULL) || STREQ(mode, "anno"))
	{
		/* if motifs_anno is null, just do against self */
		if (motifs_anno == NULL)
			run_motif_anno(motifs, num_motifs, motifs, num_motifs, strand_allow, cut);
		else
			run_motif_anno(motifs, num_motifs, motifs_anno, num_motifs_anno, strand_allow, cut);
	}
	else
		ASSERTE(0, "Mode not recognized!\n");

	for (i=0; i<num_motifs_anno; i++)
		if (motifs_anno[i].name != NULL)
			free(motifs_anno[i].name);
	for (i=0; i<num_motifs; i++)
		if (motifs[i].name != NULL)
			free(motifs[i].name);
	free(motifs_anno);
	free(motifs);
	return 0;
}

/* runs motif annotation and prints output */
void run_motif_anno(struct motif* motifs, int num_motifs, struct motif* motifs_anno, int num_motifs_anno, int strand_allow, double cut)
{
	int i, j;

	struct id_score* matched_motifs = (struct id_score*) malloc(sizeof(struct id_score)*num_motifs_anno);

	for (i=0; i<num_motifs; i++)
	{
		int num_matched_motifs = 0;
		char out_motif[MAXMOTLEN+1];
		
		for (j=0; j<num_motifs_anno; j++)
		{
			double comp = (*motif_comp_func)(&motifs[i], &motifs_anno[j], strand_allow, NULL, NULL);

			if (comp >= cut)
			{
				matched_motifs[num_matched_motifs].id = j;
				matched_motifs[num_matched_motifs].score = comp;
				num_matched_motifs++;
			}
		}


		qsort (matched_motifs, num_matched_motifs, sizeof(struct id_score), id_score_cmp);

		printf("%s\t", motif_to_str(out_motif, &motifs[i], trim_motifs, original_name));

		for (j=0; j<num_matched_motifs; j++)
			printf("%s%d",  j ? ";" : "", matched_motifs[j].id);

		printf("\t");

		for (j=0; j<num_matched_motifs; j++)
			printf("%s%f",  j ? ";" : "", matched_motifs[j].score);

		printf("\t");

		for (j=0; j<num_matched_motifs; j++)
			printf("%s%s",  j ? ";" : "", motif_to_str(out_motif, &motifs_anno[matched_motifs[j].id], trim_motifs, original_name));

		printf("\n");

	}

	free(matched_motifs);
}

/* if motifs2==NULL -> prints out mean difference for motifs in motifs1, otherwise, prints out mean difference between motifs in motifs1 and motifs2 */ 
void run_motif_meandiff(struct motif* motifs1, int num_motifs1, struct motif* motifs2, int num_motifs2, int strand_allow)
{
	double sum = 0;
	int i, j;

	if (motifs2 == NULL)
	{
		for (i=0; i<num_motifs1-1; i++)
			for (j=i+1; j<num_motifs1; j++)
				sum += (*motif_comp_func)(&motifs1[i], &motifs1[j], strand_allow, NULL, NULL);

		if (num_motifs1 > 1)
			printf("%lf\n", sum / (((num_motifs1-1) * num_motifs1) / 2));
	}
	else
	{
		for (i=0; i<num_motifs1; i++)
			for (j=0; j<num_motifs2; j++)
				sum += (*motif_comp_func)(&motifs1[i], &motifs2[j], strand_allow, NULL, NULL);

		if (num_motifs1 > 0 && num_motifs2 > 0)
			printf("%lf\n", sum / (num_motifs1 * num_motifs2));
	}
}

/* prints out the palindrome score for each motif */
void run_motif_palindrome(struct motif* motifs, int num_motifs)
{
	int i;
	char out_motif[MAXMOTLEN+1];

	for (i=0; i<num_motifs; i++)
		printf("%s\t%f\n", motif_to_str(out_motif, &motifs[i], trim_motifs, original_name), (*motif_comp_func)(&motifs[i], &motifs[i], 2, NULL, NULL));
}

/* prints out the total information content of the motif */
void run_motif_infoc(struct motif* motifs, int num_motifs)
{
	int i, j;
	char out_motif[MAXMOTLEN+1];

	for (i=0; i<num_motifs; i++)
	{
		double s = 0;

		for (j=0; j<(4*motifs[i].len); j++)
			if (motifs[i].pwm[j] > 0)
				s += motifs[i].pwm[j] * log(motifs[i].pwm[j]);

		printf("%s\t%f\n", motif_to_str(out_motif, &motifs[i], trim_motifs, original_name), 2 * motifs[i].len + s / log(2.0));
	}
}

/* runs clustering of motifs and prints out output; also updates motifs_ptr and num_motifs_ptr */
void run_motif_clust(struct motif** motifs_ptr, int *num_motifs_ptr, int strand_allow, double cut)
{
	int i;

	struct motif* motifs = *motifs_ptr;
	int num_motifs = *num_motifs_ptr;

	if (isinf(cut) != 1)
	{
		struct motif_comp *motif_comps = (struct motif_comp*) malloc(sizeof(struct motif_comp)*(num_motifs*(num_motifs-1))/2);
		int num_motif_comps = 0;
		int m1, m2;

		/* allocate extra space for the additional motifs */
		motifs = (struct motif*) realloc(motifs, (2*num_motifs-1)*sizeof(struct motif));

		/* get initial motif comparison numbers */
		for (m1=0; m1<num_motifs-1; m1++)
			for (m2=m1+1; m2<num_motifs; m2++)
				if (comp_motifs(&motif_comps[num_motif_comps], motifs, m1, m2, strand_allow, cut))
					num_motif_comps++;

		/* perform motif merges */
		while (num_motif_comps > 0)
		{
			int maxc = 0;

			/* find best merge (there must be one because motif_comps only stores >= cut) */
			for (i=1; i<num_motif_comps; i++)
				if (motif_comps[i].comp > motif_comps[maxc].comp)
					maxc = i;

			/* perform the desired merge */
			m1 = motif_comps[maxc].m1;
			m2 = motif_comps[maxc].m2;

			(*motif_merge_func)(&motifs[num_motifs++], &motifs[m1], &motifs[m2], motif_comps[maxc].rev, motif_comps[maxc].offset); 

			/* setup cluster stats on motifs... these don't set themselves */
			motifs[num_motifs-1].cluster = motifs[m1].cluster = motifs[m2].cluster = num_motifs-1;
			motifs[num_motifs-1].parents[0] = m1;
			motifs[num_motifs-1].parents[1] = m2;
			motifs[num_motifs-1].clust_corr = motif_comps[maxc].comp;
			motifs[num_motifs-1].name = NULL;
			motifs[num_motifs-1].cut = (motifs[m1].cut == motifs[m2].cut) ? motifs[m1].cut : NAN;

			/* remove the motif_comps involving m1 or m2 */
#define ANY_COMP_MATCH(c,a,b)	((c).m1==(a)||(c).m2==(a)||(c).m1==(b)||(c).m2==(b))
			for (i=0; i<num_motif_comps; )
				if (ANY_COMP_MATCH(motif_comps[num_motif_comps-1],m1,m2))
					num_motif_comps--;
				else if (ANY_COMP_MATCH(motif_comps[i],m1,m2))
					motif_comps[i] = motif_comps[--num_motif_comps];
				else
					i++;
#undef ANY_COMP_MATCH

			/* add in new motif_comps involving the new motif */
			for (i=0; i<num_motifs-1; i++)
				if ((motifs[i].cluster == i) && comp_motifs(&motif_comps[num_motif_comps], motifs, i, num_motifs-1, strand_allow, cut))
					 num_motif_comps++;
		}

		free(motif_comps);
	}

	/* print out final motifs, sorted by the order the first motif in the cluster occurred in the original file */
	{
		unsigned char *seen = (unsigned char*) calloc(num_motifs, sizeof(unsigned char));

		for (i=0; i<num_motifs; i++)
			/* only look at original motifs... everything after that will be a cluster that will 
			 * necessarily have a smaller motif it is part of */
			if (motifs[i].parents[0] == -1)
			{
				/* find the corresponding motif cluster */
				int j=i;
				while (motifs[j].cluster != j)
					j = motifs[j].cluster;
				
				if (!seen[j])
					print_motif(motifs,j);
				seen[j] = 1;
			}

		free(seen);
	}

	*num_motifs_ptr = num_motifs;
	*motifs_ptr = (struct motif*) realloc(motifs, num_motifs*sizeof(struct motif));
}

void consensus_to_motif (struct motif *m, char* motif_str)
{
	int i, j, c_idx;

	ASSERTE(m->len <= MAXMOTLEN, "Input motif too long!\n");

	m->n = 1;
	m->parents[0] = m->parents[1] = -1;
	m->clust_corr = 1;
	m->name = NULL;
	m->cut = F2LG(0.1, 0.25); /* smallest difference is between freq=0 (neg) and freq=0.25 (0) */
	
	for (i=0; i<m->len; i++)
	{
		if (motif_str[i] == '.')
			c_idx = 'N' - 'A';
		else
		{
			c_idx = toupper(motif_str[i])-'A';
			ASSERTE(c_idx >= 0 && c_idx < 26, "Unrecognized character in motif %s!\n", motif_str);
		}

		for (j=0; j<4; j++)
			m->pwm[Midx(i,j)] = ((double)ISSET(alpha[c_idx],j))/alpha_num_on[c_idx];

		/* consensus motifs must be perfect matches */
		m->cut += F2LG(1.0/alpha_num_on[c_idx], 0.25);
	}
}

/* returns the number of motifs and puts them in motifs_ptr */
int read_consensus (char* fn, struct motif **motifs_ptr, char* input_fmt)
{
	int num_motifs = 0;
	char motif_str[MAXMOTLEN+1];
	char motif_name[MAXMOTLEN+1] = "";
	int size_reserved = 1 << 5;
	int num_read = 0;
	
	struct motif* motifs = (struct motif*) malloc(size_reserved*sizeof(struct motif));

	FILE* fp = gzpopen(fn);

	while ((num_read = fscanf(fp, input_fmt, motif_str, motif_name)) != EOF)
	{
		if (num_read > 1)
			ASSERTE(strlen(motif_name) <= MAXMOTLEN, "Input motif name too long!\n");

		num_motifs++;
		
		if (size_reserved < num_motifs)
		{
			size_reserved *= 2;
			motifs = (struct motif*) realloc(motifs,size_reserved*sizeof(struct motif));
		}
		
		motifs[num_motifs-1].len = strlen(motif_str);
		motifs[num_motifs-1].cluster = num_motifs-1;
		consensus_to_motif(&motifs[num_motifs-1], motif_str);

		if (num_read != 1)
			motifs[num_motifs-1].name = strdup(motif_name);
	}

	*motifs_ptr = (struct motif*) realloc(motifs, num_motifs*sizeof(struct motif));

	pclose(fp);

	return num_motifs;
}

/* list of motifs specified in str (space separated) */
int read_consensus_string (char* motif_str, struct motif **motifs_ptr)
{
	int num_motifs = 0;
	int size_reserved = 1 << 5;

	struct motif* motifs = (struct motif*) malloc(size_reserved*sizeof(struct motif));
	
	while (*motif_str != '\0') {
		int i;

		num_motifs++;
		
		if (size_reserved < num_motifs)
		{
			size_reserved *= 2;
			motifs = (struct motif*) realloc(motifs,size_reserved*sizeof(struct motif));
		}

		for (i=0; motif_str[i] != '\0' && motif_str[i] != ' ' && motif_str[i] != '\n'; i++);
		motifs[num_motifs-1].len = i;
		motifs[num_motifs-1].cluster = num_motifs-1;
		consensus_to_motif(&motifs[num_motifs-1], motif_str);
		motif_str += i + (motif_str[i] != '\0');
	}

	*motifs_ptr = (struct motif*) realloc(motifs, num_motifs*sizeof(struct motif));

	return num_motifs;
}

/* like read_consensus, but matrices */
int read_matrix (char* fn, struct motif **motifs_ptr)
{
	int num_motifs = 0, i;
	int size_reserved = 1 << 5;
	double sum;
	char *line = NULL;
	size_t len = 0;

	struct motif* motifs = (struct motif*) malloc(size_reserved*sizeof(struct motif));
	struct motif* m = NULL;

	FILE* fp = gzpopen(fn);

	while (getline(&line, &len, fp) != -1)
	{
		if (line[0] == '>')
		{
			num_motifs++;
			
			if (size_reserved < num_motifs)
			{
				size_reserved *= 2;
				motifs = (struct motif*) realloc(motifs,size_reserved*sizeof(struct motif));
			}
			
			m = &motifs[num_motifs-1];

			m->len = 0;
			m->n = 1;
			m->cluster = num_motifs-1;
			m->parents[0] = m->parents[1] = -1;
			m->cut = NAN;

			for (i=0; line[i] != '\0' && line[i] != ' ' && line[i] != '\t' && line[i] != '\n'; i++);

			/* the next field (if available) is the cutoff */
			if (line[i] != '\n')
			{
				char* end = NULL;
				m->cut = strtod(line+i+1, &end);
				
				if (!(end == NULL || *end == ' ' || *end == '\t' || *end == '\n') || end < (line + i + 2))
					m->cut = NAN;
			}

			line[i]='\0';

			/* copy in name (ignore the >) */
			m->name = strdup(line+1);

			continue;
		}

		ASSERTE(m->len < MAXMOTLEN, "Input motif too long!\n");

		sscanf(line, "%*s %lf %lf %lf %lf", &m->pwm[Midx(m->len,0)], &m->pwm[Midx(m->len,1)], &m->pwm[Midx(m->len,2)], &m->pwm[Midx(m->len,3)]);

		for (i=0; i<4; i++)
			if (m->pwm[Midx(m->len,i)] < 0)
				m->pwm[Midx(m->len,i)] = 0;

		sum = m->pwm[Midx(m->len,0)] + m->pwm[Midx(m->len,1)] + m->pwm[Midx(m->len,2)] + m->pwm[Midx(m->len,3)];
		for (i=0; i<4; i++)
			if (sum <= 0)
				m->pwm[Midx(m->len,i)] = 0.25;
			else
				m->pwm[Midx(m->len,i)] /= sum;

		m->len++;
	}

	*motifs_ptr = (struct motif*) realloc(motifs, num_motifs*sizeof(struct motif));

	pclose(fp);
	free(line);

	return num_motifs;
}

int read_logodds (char* fn, struct motif **motifs_ptr)
{
	int num_motifs = 0, i;
	int size_reserved = 1 << 5;
	double sum, cfact;
	char *line = NULL;
	size_t len = 0;

	struct motif* motifs = (struct motif*) malloc(size_reserved*sizeof(struct motif));
	struct motif* m = NULL;

	FILE* fp = gzpopen(fn);

	int len_left = 0;

	while (getline(&line, &len, fp) != -1)
	{
		int p;
		if (len_left == 0)
		{
			ASSERTE(line[0] == 'X', "Invalid log-odds input format!\n");

			num_motifs++;
			
			/* must have enough room for all motifs + all clusters */
			if (size_reserved < num_motifs)
			{
				size_reserved *= 2;
				motifs = (struct motif*) realloc(motifs,size_reserved*sizeof(struct motif));
			}
			
			m = &motifs[num_motifs-1];

			/* get the length by the number of X's */
			for (m->len=0; line[m->len] == 'X'; m->len++);
			ASSERTE(m->len < MAXMOTLEN, "Input motif too long!\n");

			/* copy in name */
			for (i=m->len+1; line[i] != '\0' && line[i] != ' ' && line[i] != '\t' && line[i] != '\n'; i++);
			line[i]='\0';
			m->name = strdup(line+m->len+1);

			m->n = 1;
			m->cluster = num_motifs-1;
			m->parents[0] = m->parents[1] = -1;

			/* get the cutoff */
			ASSERTE(getline(&line, &len, fp) != -1, "Missing log-odds threshold!\n");
			sscanf(line, "%lf", &m->cut);

			/* number of positions we have to read in */
			len_left = m->len;

			continue;
		}

		p = m->len - len_left;

		/* convert to matrix format */
		sscanf(line, "%*s %lf %lf %lf %lf", &m->pwm[Midx(p,0)], &m->pwm[Midx(p,1)], &m->pwm[Midx(p,2)], &m->pwm[Midx(p,3)]);

		/* note: we correct for shifts of the distribution, but not taking a product of the whole matrix
		 * (e.g. if it isn't in log(2)) */

		/* first get sum of initially converted values */
		sum = 0;
		for (i=0; i<4; i++)
			sum += LG2F(m->pwm[Midx(p,i)], lg_bg[i]);

		/* now compute correction factor so that the final values add to 1 */
		cfact = log((sum + lg_pseudo)/(1 + lg_pseudo))/log(2);

		/* we have to also adjust m->cut to account for changes across the row */
		m->cut -= cfact;

		/* now recompute values with correction factor considered */
		for (i=0; i<4; i++)
			m->pwm[Midx(p,i)] = LG2F(m->pwm[Midx(p,i)] - cfact, lg_bg[i]);

		/* ensure values are a probability distribution */
		for (i=0; i<4; i++)
			if (m->pwm[Midx(p,i)] < 0)
				m->pwm[Midx(p,i)] = 0;

		sum = m->pwm[Midx(p,0)] + m->pwm[Midx(p,1)] + m->pwm[Midx(p,2)] + m->pwm[Midx(p,3)];
		for (i=0; i<4; i++)
			if (sum <= 0)
				m->pwm[Midx(p,i)] = 0.25;
			else
				m->pwm[Midx(p,i)] /= sum;

		len_left--;
	}

	ASSERTE(len_left == 0, "Invalid log-odds input format!\n");

	*motifs_ptr = (struct motif*) realloc(motifs, num_motifs*sizeof(struct motif));

	pclose(fp);
	free(line);

	return num_motifs;
}


/* returns the number of motifs; allocates enough memory for all motifs
 * plus all the possible merges (which is another n-1) */
int read_alignment (char* fn, struct motif **motifs_ptr)
{
	int num_motifs = 0, i;
	int size_reserved = 1 << 5;
	int num_str = 0;
	char *line = NULL;
	size_t len = 0;

	struct motif* motifs = (struct motif*) malloc(size_reserved*sizeof(struct motif));
	struct motif* m = NULL;

	int base_count[MAXMOTLEN*4];

	FILE* fp = gzpopen(fn);

	while (getline(&line, &len, fp) != -1)
	{
		int line_len = strlen(line);
		if (line[line_len-1] == '\n')
			line[--line_len] = '\0';

		/* either file has just one motif, or motifs are separated by > lines or 
		 * by blank lines */
		if (num_motifs == 0 || line[0] == '>' || line[0] == '\0')
		{
			if (num_motifs > 0)
			{
				if (num_str == 0)
					num_motifs--;
				else
					for (i=0; i<4*m->len; i++)
						m->pwm[i] = ((double)base_count[i])/num_str;
			}

			num_motifs++;
			
			if (size_reserved < num_motifs)
			{
				size_reserved *= 2;
				motifs = (struct motif*) realloc(motifs,size_reserved*sizeof(struct motif));
			}
			
			m = &motifs[num_motifs-1];

			m->len = 0;
			m->n = 1;
			m->cluster = num_motifs-1;
			m->parents[0] = m->parents[1] = -1;

			if (line[0] == '>')
			{
				for (i=0; line[i] != '\0' && line[i] != ' ' && line[i] != '\t' && line[i] != '\n'; i++);
				line[i]='\0';
				m->name = (char *) malloc(sizeof(char) * i);
				strcpy(m->name, line+1);
			}
			else
				m->name = NULL;

			
			num_str = 0;

			memset(base_count, 0, sizeof(int)*MAXMOTLEN*4);

			if (line[0] == '>' || line[0] == '\0')
				continue;
		}

		if (m->len == 0)
			m->len = line_len;

		ASSERTE(line_len == m->len, "Length of input sequences does not match!\n");

		ASSERTE(m->len < MAXMOTLEN, "Input motif too long!\n");

		for (i=0; i<m->len; i++)
			switch(line[i])
			{
				case 'A': case 'a':
					base_count[Midx(i,0)]++;
					break;
				case 'C': case 'c':
					base_count[Midx(i,1)]++;
					break;
				case 'G': case 'g':
					base_count[Midx(i,2)]++;
					break;
				case 'T': case 't':
					base_count[Midx(i,3)]++;
					break;
				default:
					ASSERTE(0, "Invalid character\n");
			}

		num_str++;
	}

	if (num_motifs > 0)
	{
		if (num_str == 0)
			num_motifs--;
		else
			for (i=0; i<4*m->len; i++)
				m->pwm[i] = ((double)base_count[i])/num_str;
	}

	*motifs_ptr = (struct motif*) realloc(motifs, num_motifs*sizeof(struct motif));

	pclose(fp);
	free(line);

	return num_motifs;
}


double compute_corr_div(struct motif* m)
{
	int i, j;
	double sum_sq = 0;
	
	for (i=0; i<m->len; i++)
		for (j=0; j<4; j++)
			sum_sq += m->pwm[Midx(i,j)] * m->pwm[Midx(i,j)];

	return sqrt(sum_sq - (m->len/4.0));
}

#define TLEN(offset) (MAX(m2->len+(offset), m1->len) - MIN(offset, 0))

/* set the best_ values; break ties by taking the shorter */
#define SET_IF_BETTER(conv, offset, rev)  do {                      \
	if ((conv) > best_conv ||                                       \
	  (((conv) == best_conv) && TLEN(offset) < TLEN(best_offset)))  \
	{                                                               \
		best_conv = (conv);                                         \
		best_offset = (offset);                                     \
		best_rev = (rev);                                           \
	}} while (0)

/* in these functions: 
 * start and end refer to the start (inclusive) and end (exclusive) position in m1
 * offset is how much m2 starts ahead of m1 */
double pearson_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr)
{
	double div = compute_corr_div(m1) * compute_corr_div(m2);

	double best_conv = -2.0*div;

	int offset, best_offset=0, best_rev=0, rev;

	if ((div*div) < (EPSILON*EPSILON))
		return 0.0;

	for (offset=MAX(1-m2->len, -max_offset); offset<MIN(m1->len, max_offset); offset++)
	{
		/* start and end of overlapping region */
		int start = MAX(offset, 0);
		int end = MIN(m2->len+offset, m1->len);
		
		for (rev=0; rev<=1; rev++)
			if ((strand_allow>>rev)&1)
			{
				int i, j;
				double conv = 0;

				for (i=start; i<end; i++)
					for (j=0; j<4; j++)
						conv += m1->pwm[Midx(i,j)] * m2->pwm[RFMidx(i-offset,j,m2->len,rev)];

				SET_IF_BETTER(conv-(end-start)/4.0, offset, rev);
			}
	}

	if (best_rev_ptr != NULL)
		*best_rev_ptr = best_rev;

	if (best_offset_ptr != NULL)
		*best_offset_ptr = best_offset;

	return best_conv/div;
}

/* allows for up to motif_comp_param gaps; with motif_comp_param=0 same as pearson_motif_comp */
double peargap_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr)
{

	/* i/j - 0-indexed position in m1/m2 of length I/J
	 *
	 * uses dynamic programming with the following recursion to 
	 * fill in best_conv:
	 *    V[ i, j, gap ] = max( V[i-1, j-1, gap] + SCORE[i,j],
	 *                          V[i-1, j,   gap-1], 
	 *                          V[i,   j-1, gap-1]
	 *                          V[i,   j, gap-1])
	 * where:
	 *    V[i, -1, 0] = V[-1, j, 0] = 0 
	 * and the answer is:
	 *    max(V[I-1, *, gap], V[*, J-1, gap])
	 *
	 * because gaps are allowed at the beginning or end for free
	 */

#define Vidx(i, j, g) ((g) * ((1+m1->len) * (1+m2->len)) + (i+1) * (1+m2->len) + (j+1))
#define SET_MAX(a,b) do { if ((b) > (a)) (a) = (b); } while (0)

	double div = compute_corr_div(m1) * compute_corr_div(m2);
	double best_conv = -(m1->len+m2->len);
	int rev, i, j, g, k;

	ASSERTE(best_rev_ptr == NULL && best_offset_ptr == NULL, "peargap comparison function cannot be used for clustering.\n");

	if ((div*div) < (EPSILON*EPSILON))
		return 0.0;

	for (rev=0; rev<=1; rev++)
		if ((strand_allow>>rev)&1)
		{
			double *V = (double *) calloc ((1+m1->len) * (1+m2->len) * (motif_comp_param+1), sizeof(double));
			
			for (i=0; i<m1->len; i++)
				for (j=0; j<m2->len; j++)
					for (g=0; g<=motif_comp_param; g++)
					{
						/* NOTE: j indexes the position in the second sequence here, not the base like in pearson_motif_comp */
						V[Vidx(i,j,g)] = V[Vidx(i-1, j-1, g)] - 1.0/4.0;
						
						for (k=0; k<4; k++)
							V[Vidx(i,j,g)] += m1->pwm[Midx(i,k)] * m2->pwm[RFMidx(j,k,m2->len,rev)];

						if (g > 0)
						{
							SET_MAX(V[Vidx(i,j,g)], V[Vidx(i-1, j, g-1)]);
							SET_MAX(V[Vidx(i,j,g)], V[Vidx(i, j-1, g-1)]);
							SET_MAX(V[Vidx(i,j,g)], V[Vidx(i, j, g-1)]);
						}
					}

			for (i=0; i<m1->len; i++)
				SET_MAX(best_conv, V[Vidx(i,m2->len-1, motif_comp_param)]);

			for (j=0; j<m2->len; j++)
				SET_MAX(best_conv, V[Vidx(m1->len-1,j, motif_comp_param)]);

			free(V);
		}

	return best_conv/div;
}

#define DO_IF_ALL_COLS(cmd, expr) do { for (j=0; j<4; j++) if (!(expr)) break; if (j==4) cmd;} while (0)
#define GTE(val) (((val) >= EPSILON) ? 1 : 0)

/* if one is exactly a subset of the other (or exactly the same), returns 1 otherwise, returns -(# differences) */
double ndiff_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr)
{
	int best_conv = -(m1->len+m2->len);
	int offset, best_offset=0, best_rev=0, rev;

	for (offset=MAX(1-m2->len, -max_offset); offset<MIN(m1->len, max_offset); offset++)
	{
		int i, j, num;
		/* start and end of overlapping region */
		int start = MAX(offset,0);
		int end = MIN(m2->len+offset, m1->len);
		
		/* for the non-overlapping regions: get length for each motif and number of Ns */
		int num_m1_nonover_N = 0, num_m2_nonover_N = 0, num_m1_nonover = 0, num_m2_nonover = 0; 

		int final_motif_len = TLEN(offset);

		for (i=0; i<m1->len; i++)
			if (i < start || i >= end)
			{
				num_m1_nonover++;
				DO_IF_ALL_COLS(num_m1_nonover_N++, GTE(m1->pwm[Midx(i,j)]));
			}

		for (i=0; i<m2->len; i++)
			if (i < -offset || i >= m1->len-offset)
			{
				num_m2_nonover++;
				DO_IF_ALL_COLS(num_m2_nonover_N++, GTE(m2->pwm[Midx(i,j)]));
			}

		for (rev=0; rev<=1; rev++)
			if ((strand_allow>>rev)&1)
			{
				/* check if m1 is a superset (i.e. matches every sequence) of m2 */
				num=0; 
				for (i=start; i<end; i++)
					DO_IF_ALL_COLS(num++, GTE(m1->pwm[Midx(i,j)]) >= GTE(m2->pwm[RFMidx(i-offset,j,m2->len,rev)]));

				if ((num + num_m2_nonover + num_m1_nonover_N) == final_motif_len)
				{
					best_conv = 1;
					best_offset = offset;
					best_rev = rev;
					goto FoundBest;
				}

				/* check if m2 is a superset of m1 */
				num=0; 
				for (i=start; i<end; i++)
					DO_IF_ALL_COLS(num++, GTE(m1->pwm[Midx(i,j)]) <= GTE(m2->pwm[RFMidx(i-offset,j,m2->len,rev)]));

				if ((num + num_m1_nonover + num_m2_nonover_N) == final_motif_len)
				{
					best_conv = 1;
					best_offset = offset;
					best_rev = rev;
					goto FoundBest;
				}

				/* compute the negative of the total number of non-matching positions */
				num = -final_motif_len + num_m1_nonover_N + num_m2_nonover_N;
				for (i=start; i<end; i++)
					DO_IF_ALL_COLS(num++, GTE(m1->pwm[Midx(i,j)]) == GTE(m2->pwm[RFMidx(i-offset,j,m2->len,rev)]));

				SET_IF_BETTER(num, offset, rev);
			}
	}

FoundBest:
	if (best_rev_ptr != NULL)
		*best_rev_ptr = best_rev;

	if (best_offset_ptr != NULL)
		*best_offset_ptr = best_offset;

	return (double) best_conv;
}

/* outputs the number of positions that are the same (pattern of >= EPSILON the same), but not N */
double nsame_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr)
{
	int best_conv = -(m1->len+m2->len);
	int offset, best_offset=0, best_rev=0, rev;

	for (offset=MAX(1-m2->len, -max_offset); offset<MIN(m1->len, max_offset); offset++)
	{
		int i, j, num;

		/* start and end of overlapping region */
		int start = MAX(offset,0);
		int end = MIN(m2->len+offset, m1->len);
		
		for (rev=0; rev<=1; rev++)
			if ((strand_allow>>rev)&1)
			{
				/* compute the negative of the total number of non-matching positions */
				num = 0;
				for (i=start; i<end; i++)
					/* check that the position is not an N */
					for (j=0; j<4; j++) 
						if (!GTE(m1->pwm[Midx(i,j)])) 
						{
							/* add 1 to num if the column has the same pattern of non-1 */
							DO_IF_ALL_COLS(num++, GTE(m1->pwm[Midx(i,j)]) == GTE(m2->pwm[RFMidx(i-offset,j,m2->len,rev)]));
							break;
						}

				SET_IF_BETTER(num, offset, rev);
			}
	}

	if (best_rev_ptr != NULL)
		*best_rev_ptr = best_rev;

	if (best_offset_ptr != NULL)
		*best_offset_ptr = best_offset;

	return (double) best_conv;
}


/* also uses global variable motif_comp_param to indicate the window */
double wprod_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr)
{
	double best_conv = -1;

	int offset, best_offset=0, best_rev=0, rev;

	for (offset=MAX(1-m2->len, -max_offset); offset<MIN(m1->len, max_offset); offset++)
	{
		/* start and end of overlapping region */
		int start = MAX(offset,0);
		int end = MIN(m2->len+offset, m1->len);

		/* start of the window */
		int wstart;

		for (wstart=start; wstart < end; wstart++)
		{
			int wend = MIN(wstart+motif_comp_param, end);

			for (rev=0; rev<=1; rev++)
				if ((strand_allow>>rev)&1)
				{
					int i, j;
					double conv = 0;

					for (i=wstart; i<wend; i++)
						for (j=0; j<4; j++)
							conv += m1->pwm[Midx(i,j)] * m2->pwm[RFMidx(i-offset,j,m2->len,rev)];

					SET_IF_BETTER(conv, offset, rev);
				}
		}
	}

	if (best_rev_ptr != NULL)
		*best_rev_ptr = best_rev;

	if (best_offset_ptr != NULL)
		*best_offset_ptr = best_offset;

	return best_conv/motif_comp_param;
}

/* converts a matrix motif into a log-odds (space must be allocated) and enforces that
 * the maximum value for any row is 0 (returning the corrected cut) */
double matrix_to_norm_logodds(double *lg, struct motif* m)
{
	int i, j;

	double cut = m->cut;
	for (i=0; i<m->len; i++)
	{
		double max = -INFINITY;
		for (j=0; j<4; j++)
		{
			lg[Midx(i,j)] = F2LG(m->pwm[Midx(i,j)], lg_bg[j]);

			if (lg[Midx(i,j)] > max)
				max = lg[Midx(i,j)];
		}

		cut -= max;
		for (j=0; j<4; j++)
			lg[Midx(i,j)] -= max;
	}

	return cut;
}

/* three unsigned long long int values --
 * number of matching for both, only motif1, only motif2 */
struct match_counts 
{
	unsigned long long int both;
	unsigned long long int m1;
	unsigned long long int m2;
};

struct match_counts count_matches(double *g1, double *g2, int l1, int l2, double c1, double c2, int p1, int p2, int rev)
{
	struct match_counts out = {0,0,0};
	int j;

	if (p1 >= l1 && p2 >= l2)
	{
		if (c1 <= 0 && c2 <= 0)
			out.both = 1;
		else if (c1 <= 0)
			out.m1 = 1;
		else if (c2 <= 0)
			out.m2 = 1;
	}
	else if (c1 <= 0 || c2 <= 0)
		for (j=0; j<4; j++)
		{
			struct match_counts new = count_matches(g1, g2, l1, l2, 
				(p1 >= 0 && p1 < l1) ? (c1 - g1[Midx(p1, j)])            : c1, 
				(p2 >= 0 && p2 < l2) ? (c2 - g2[RFMidx(p2, j, l2, rev)]) : c2, 
				p1+1, p2+1, rev);
			out.both += new.both;
			out.m1 += new.m1;
			out.m2 += new.m2;
		}

	return out;
}

/* computes intersection divided by union of matching sequences at best offset 
 * ignores background frequencies */
double fracmatch_motif_comp(struct motif* m1, struct motif* m2, int strand_allow, int *best_rev_ptr, double *best_offset_ptr)
{
	double best_conv = 0;
	int offset, best_offset=0, best_rev=0, rev;

	/* convert matrices into log-odds */
	double *g1 = (double*) alloca(sizeof(double) * m1->len * 4);
	double *g2 = (double*) alloca(sizeof(double) * m2->len * 4);
	double c1 = matrix_to_norm_logodds(g1, m1);
	double c2 = matrix_to_norm_logodds(g2, m2);

	if (isnan(c1) || isnan(c2))
		return NAN;
	
	/* maximum score may occur when completely non-overlapping */
	for (offset=MAX(-m2->len, -max_offset); offset<MIN(m1->len+1, max_offset); offset++)
	{
		for (rev=0; rev<=1; rev++)
			if ((strand_allow>>rev)&1)
			{
				struct match_counts out;

				if (offset >= 0)
					out = count_matches(g1, g2, m1->len, m2->len, c1, c2, 0, -offset, rev);
				else
					out = count_matches(g1, g2, m1->len, m2->len, c1, c2, offset, 0, rev);

				SET_IF_BETTER(((double) out.both) / ((double) (out.both + out.m1 + out.m2)), offset, rev);
			}
	}

	if (best_rev_ptr != NULL)
		*best_rev_ptr = best_rev;

	if (best_offset_ptr != NULL)
		*best_offset_ptr = best_offset;

	return best_conv;
}

/* fills in a motif_comp for motifs m1 and m2 and returns 1 if the corresponding best comp
 * is greater than cut; otherwise, returns 0*/
int comp_motifs(struct motif_comp* mc, struct motif* motifs, int m1, int m2, int strand_allow, double cut)
{
	int rev;
	double offset;
	double comp = (*motif_comp_func)(&motifs[m1], &motifs[m2], strand_allow, &rev, &offset);
	
	if (comp >= cut)
	{
		mc->m1 = m1;
		mc->m2 = m2;
		mc->rev = rev;
		mc->offset = offset;
		mc->comp = comp;
		return 1;
	}

	return 0;
}

/* if position exists in M, returns corresponding PWM val, otherwise, returns 0.25 */
#define InRange(M,pos) 				((pos)>=0 && (pos)<(M)->len)
#define SafePWMVal(M, pos, chr, rev) (InRange(M,pos) ? M->pwm[RFMidx(pos,chr,(M)->len,rev)] : 0.25)
	
/* merge the motifs m1, m2 into m as a weighted average of the two motifs */
void weighted_motif_merge(struct motif* m, struct motif* m1, struct motif* m2, int rev, int offset)
{
	int i, j;

	int off1 = MIN(offset,0);
	int off2 = off1-offset;
	
	m->len = MAX(m2->len+offset, m1->len) - off1;

	ASSERTE(m->len <= MAXMOTLEN, "A motif too long was constructed!\n");

	m->n = m1->n + m2->n;

	for (i=0; i<m->len; i++)
		for (j=0; j<4; j++)
			m->pwm[Midx(i,j)] = (SafePWMVal(m1,off1+i,j,0)   * m1->n 
			                   + SafePWMVal(m2,off2+i,j,rev) * m2->n) / m->n;
}

/* merge using >= logic... also trim motif */
void or_motif_merge(struct motif* m, struct motif* m1, struct motif* m2, int rev, int offset)
{
	/* start by running weighted_motif_merge */
	weighted_motif_merge(m,m1,m2,rev,offset);

	/* now take each of the columns and equalize them out */
	int i, j, num_match, num_skip_start = 0, num_skip_end = 0, at_front = 1;

	for (i=0; i<m->len; i++)
	{
		num_match = 0;
		for (j=0; j<4; j++)
			if (m->pwm[Midx(i,j)] >= EPSILON)
				num_match++;

		if (num_match != 4)
			at_front = num_skip_end = 0;
		else 
		{
			if (at_front)
				num_skip_start++;
			num_skip_end++;
		}

		for (j=0; j<4; j++)
			if (m->pwm[Midx(i,j)] >= EPSILON)
				m->pwm[Midx(i,j)] = 1.0/num_match;
	}

	m->len -= num_skip_start + num_skip_end; 

	if (num_skip_start > 0)
		memmove(m->pwm, m->pwm + 4*num_skip_start, m->len * 4 * sizeof(m->pwm[0]));
}


void print_known_parents(struct motif* motifs, int m1)
{
	if (motifs[m1].parents[0] == -1)
	{
		if (number_output)
			printf("%d", m1);
		else
		{
			char out_motif[MAXMOTLEN+1];
			printf("%s", motif_to_str(out_motif, &motifs[m1], trim_motifs, original_name));
		}
	}
	else
	{
		if (tree_output)
			printf("(");
		print_known_parents(motifs, motifs[m1].parents[0]);
		printf(tree_output ? "," : ";");
		print_known_parents(motifs, motifs[m1].parents[1]);
		if (tree_output)
			printf(")");

		if (correlation_output)
			printf(":%lf", motifs[m1].clust_corr);
	}
}

char epsilon_pwm_to_char(double *pwm_col)
{
	char *bin2char = "nACMGRSVTWYHKDBN";
	int bitrep = 0; 
	int j;

	for (j=0; j<4; j++)
		if (pwm_col[j] >= EPSILON)
			bitrep |= 1 << j;

	return bin2char[bitrep]; 
}

char correction_pwm_to_char(double *pwm_col)
{
	/* correction for degeneracy */
	const double correction[4] = {0.5, 2.0/3, 5.0/6, 1};

	char best_char = 'N';
	double best_score = 1e-10;
	int c;

	for (c=0; c<26; c++)
		if (alpha_num_on[c] > 0)
		{
			double score = -correction[alpha_num_on[c]-1];

			int j;
			for (j=0; j<4; j++)
				score += ISSET(alpha[c],j)*pwm_col[j];

			if (score > best_score)
			{
				best_score = score;
				best_char = 'A' + c;
			}
		}
	
	return best_char; 
}

/* writes sequence of motif to str (and returns a pointer to the first non-N of str)
 * one after last non-N is made into \0 
 * if all the characters are N, return value points to last N */
char* motif_to_str(char* str, struct motif* motif, unsigned char trim_motif, unsigned char original_name)
{
	/* use these to strip away N's from start/end */
	int first_non_n = -1;
	int last_non_n = 0;
	int i, j;

	if (original_name && motif->name != NULL)
		return motif->name;
	
	for (i=0; i<motif->len; i++)
	{
		double pos[4];
		for (j=0; j<4; j++)
			pos[j] = motif->pwm[RFMidx(i,j,motif->len,reverse_output)];
		str[i] = (*pwm_to_char_func)(pos);

		if (str[i] != 'N')
		{
			if (first_non_n == -1)
				first_non_n = i;

			last_non_n = i;
		}
	}

	if (!trim_motif)
	{
		str[motif->len] = '\0';
		return str;
	}

	if (first_non_n == -1)
		last_non_n = first_non_n = motif->len-1;

	str[last_non_n+1] = '\0';

	return str + first_non_n;
}

/* outputs a motif; only used in clustering mode */
void print_motif(struct motif* motifs, int m1)
{
	char out_motif[MAXMOTLEN+1];
	char* out_motif_n = motif_to_str(out_motif, &motifs[m1], trim_motifs, 0);

	if (output_mode == 2)
	{
		/* log odds output */
		int i, j;
		for (i=0; i<motifs[m1].len; i++)
			printf("X");
		printf("\t");
		print_known_parents(motifs, m1);
		printf("\n%lf\n", motifs[m1].cut);

		for (i=0; i<motifs[m1].len; i++)
		{
			printf("%c", out_motif[i]=='\0' ? 'N' : out_motif[i]);
			for (j=0; j<4; j++)
				printf(" % f", F2LG(motifs[m1].pwm[RFMidx(i,j,motifs[m1].len,reverse_output)],  lg_bg[j]));
			printf("\n");
		}
	}
	else
	{
		/* either standard matrix or consensus output */
		printf("%s%s\t", output_mode ? ">" : "", out_motif_n);
		print_known_parents(motifs, m1);
		printf("\n");
		if (output_mode == 1)
		{
			int i, j;

			for (i=0; i<motifs[m1].len; i++)
			{
				printf("%c", out_motif[i]=='\0' ? 'N' : out_motif[i]);
				for (j=0; j<4; j++)
					printf(" %f", motifs[m1].pwm[RFMidx(i,j,motifs[m1].len,reverse_output)]);
				printf("\n");
			}
		}
	}
}

