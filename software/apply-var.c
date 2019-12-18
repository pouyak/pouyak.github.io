/* Written by Pouya Kheradpour (pouyak@gmail.com)
 * Copyright (c) 2005-2015 Massachusetts Institute of Technology
 * Released under the MIT license
 * https://pouya.kheradpour.com/software/LICENSE.mit */
/* version 20140120 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <search.h>

#include "pk.h"

#define Aidx(c, p) ((c)*(max_seq_len+1) + (p))

/* makes an individual index het (must be positive);
 * notice the inverse is the same */
#define HET(x) (-(x)-1)

/* removes het if negative */
#define NOHET(x) (((x) < 0) ? HET(x) : (x))

void rev_comp(char* seq, unsigned int len);

struct mod
{
	char chr[MAXCHRNAMELEN+1];
	long long int start;
	long long int end;

	/* a negative ind indicates het for the -ind-1 individual */
	int *inds;
	int num_inds;
	struct mod* next;
	char *seq;
};

struct chrom
{
	/* a -1 indicates + -- all other chromosomes */
	int *inds;
	int num_inds;

	struct mod** mods_ptr;
	int num_mods;
	int mods_alloc;
};

void update_mods(FILE* fp, struct mod* mlist, char** line_ptr, size_t *line_alloc_ptr, char ***ind_names_ptr, void **ind_tree_ptr, int *num_inds_ptr);

int make_chroms (struct mod* mlist, struct chrom** chroms_ptr, int* max_seq_len_ptr, int num_inds);

int main (int argc, char** argv)
{
	FILE *mod_fp = NULL;
	FILE *seq_fp = NULL;
	int i;
	size_t line_alloc = 0, seq_alloc = 0;
	char *line = NULL, *seq = NULL;
	struct chrom* chroms = NULL;
	char *align = NULL;
	int align_alloc = 0; /* in terms of maximum sequence length */

	int num_inds = 0;
	void *ind_tree = NULL;
	char **ind_names = NULL;

	unsigned char error_mode = 1;

	/* very similar to grep-overlap's linked list data structure 
	 * we use the first item in the list to hold the region of the sequence we have read
	 * in */
	struct mod *mlist = (struct mod*) calloc(1, sizeof(struct mod));

	for (i=1; i<(argc-1); i++)
		if (STREQ(argv[i], ""))
			continue; 
		else if (STREQ(argv[i], "-v"))
			mod_fp =  gzpopen(argv[++i]);
		else if (STREQ(argv[i], "-s"))
			seq_fp =  gzpopen(argv[++i]);
		else if (STREQ(argv[i], "-e"))
			error_mode =  atoi(argv[++i]);
		else
			ASSERTE(0, "Unrecognized command line argument!\n");

	if (mod_fp == NULL)
	{
		printf("USAGE: %s [OPTIONS]\n", argv[0]);
		printf("    -v    variation file (id chr start stop replace bases individuals) prefix hets with ? [required]\n");
		printf("    -s    sequence file (from extract-mfa -k 1) [default: -]\n");
		printf("    -e    error mode; 1 to fail on errors [default: 1]\n");
		return 1;
	}

	if (seq_fp == NULL)
		seq_fp =  gzpopen("-");

	/* read in sequence */
	while (getline(&line, &line_alloc, seq_fp) != -1)
	{
		unsigned char flip = 0;
		int num_read;
		int seq_len, max_seq_len;

		int num_chroms;

		printf("%s", line);
		
		{
			unsigned char strand;
			num_read = sscanf(line, "%*s %s %Ld %Ld %c", mlist->chr, &mlist->start, &mlist->end, &strand);
			if (num_read == 4 && strand == '-')
				flip = 1;
		}

		getline(&line, &line_alloc, seq_fp);

		max_seq_len = seq_len = getline(&seq, &seq_alloc, seq_fp) - 1;

		/* get rid of the \n */
		seq[seq_len] = '\0';
		
		ASSERTE(seq_len == (mlist->end - mlist->start +1), "%s %Ld sequence length does not match region start and end\n", mlist->chr, mlist->start);

		/* read and discard the separating blank line */
		ASSERTE(getline(&line, &line_alloc, seq_fp) == 1, "Improperly formatted sequence input\n");

		/* update the list of mods appropriately */
		update_mods(mod_fp, mlist, &line, &line_alloc, &ind_names, &ind_tree, &num_inds);

		/* get list of chromosomes to print out (one for each individual with
		 * some modification, grouping those individuals with an identical set */
		num_chroms = make_chroms(mlist, &chroms, &max_seq_len, num_inds);

		/* apply modifications */
		if (flip)
			/* if the input is reverse complemented, we reverse complement the input
			 * and then again the output */
			rev_comp(seq, seq_len);

		/* allocate space for the alignment */
		if (((max_seq_len+1) * num_chroms) > align_alloc)
		{
			align_alloc = (max_seq_len+1) * num_chroms;
			align = (char *) realloc(align, sizeof(char) * align_alloc);
		}

		/* and copy in the reference (which has number 0); make sure to also copy \0*/
		for (i=0; seq[i] != '\0'; i++)
			align[i] = tolower(seq[i]);
		align[i] = '\0';

		/* apply modifications; note: chromosome 0 is always reference and must have no mods */
		for (i=1; i<num_chroms; i++)
		{
			int mi;
			
			/* make a copy of the current reference (including the \0) */
			memcpy(align + Aidx(i, 0), align, seq_len+1);

			/* make the appropriate modifications */
			for (mi=0; mi<chroms[i].num_mods; mi++)
			{
				struct mod* m = chroms[i].mods_ptr[mi];

				int k;

				int alt_length = strlen(m->seq);

				/* the first alternate base is usually the one we start with, unless
				 * some modified bases are not within our region and thus are assigned
				 * to unused bases here (one for each base) */
				int alt_start = MIN(alt_length, MAX(0, mlist->start - m->start));

				/* this is one past the position in the alt sequence that we align
				 * to a base in the reference */
				int alta_end = MIN(alt_length, MIN(mlist->end,m->end)-m->start + 1);

				/* alta_end <= i < alt_end defines the positions to place immediately 
				 * after the last reference aligned base... note: if the reference
				 * bases go beyond the end of the sequence, then we discard all bases
				 * and alt_start = alta_end; otherwise alta_end == seqlen */
				int alt_end = m->end <= mlist->end ? alt_length : alta_end;

				/* bases to set as - after the ones filled in with alt */
				int del_bases = MAX(0, MIN(mlist->end,m->end) - MAX(mlist->start,m->start) + 1 - alt_length + alt_start);

				/* get the first position in the alignment to modify */
				int cb = 0;
				for (k=mlist->start; ;cb++)
					if (align[cb] != '-')
					{
						k++;
						if (k>m->start)
							break;
					}

				/* modify the bases that are properly aligned */
				for (k=alt_start; k<alta_end; cb++)
					if (align[cb] != '-')
					{
						/* when error_mode != 0 -- do not allow a double modification unless this doesn't actually change the underlying base 
						 * some vcf files are just broken and will have these double mutations */
						ASSERTE(error_mode == 0 || align[cb] == align[Aidx(i, cb)] || (align[Aidx(i, cb)] == '-' && toupper(align[cb]) == m->seq[k]), "Double modification found: %s %Ld\n", m->chr, m->start);

						/* NOTE: when there is an double modification, we default to deletion (if it occurred) or the last modification 
						 *       because modifications can now be ignored, this can lead to identical sequences being outputted! 
						 * otherwise, do not change bases where we match the reference */
						if (toupper(align[cb]) != m->seq[k] && align[Aidx(i, cb)] != '-')
							align[Aidx(i, cb)] = m->seq[k];
						k++;
					}

				/* place any other bases possible to positions that are already gapped
				 * (e.g. due to another chrom) */
				for (; k<alt_end; cb++)
					if (align[cb] == '-')
						align[Aidx(i, cb)] = m->seq[k++];
					else
						break;

				/* a gap must be produced in the alignment for the additional bases, right before
				 * the position that is currently cb */
				if (k < alt_end)
				{
					int ii;
					int gap = alt_end - k;

					for (ii=0; ii<=i; ii++)
					{
						/* make sure to also copy the \0 */
						memmove(align + Aidx(ii, cb) + gap, align + Aidx(ii, cb), seq_len - cb + 1);
						
						/* and fill in -'s */
						memset(align + Aidx(ii, cb), '-', sizeof(char)*gap);
					}

					seq_len += gap;

					/* fill in the final bases */
					for (; k<alt_end; cb++)
						align[Aidx(i, cb)] = m->seq[k++];
				}

				/* delete any additional bases necessary */
				for (k=0; k < del_bases; cb++)
					if (align[cb] != '-')
					{
						/* we shouldn't be deleting a base that has already been changed when error_mode != 0*/
						ASSERTE(error_mode == 0 || align[cb] == align[Aidx(i, cb)], "Double modification found: %s %Ld\n", m->chr, m->start);
						
						/* when error_mode is 0 and this base has already been modified, we default to removal */
						k++;
						align[Aidx(i, cb)] = '-';
					}
			}
		}
		
		/* output resulting sequences */
		for (i=0; i<num_chroms; i++)
		{
			int j;

			printf(">%s", chroms[i].inds[0] == -1 ? "+" : ind_names[chroms[i].inds[0]]);
			for (j=1; j<chroms[i].num_inds; j++)
				printf(";%s", ind_names[chroms[i].inds[j]]);

			printf("\n");

			if (flip)
				rev_comp(align + Aidx(i,0), seq_len);

			printf("%s\n", align + Aidx(i,0));
		}

		printf("\n");

		/* clean up chroms */
		for (i=0; i<num_chroms; i++)
		{
			free(chroms[i].inds);
			free(chroms[i].mods_ptr);
		}
	}

	free(align);
	free(chroms);
	free(seq);
	free(line);
	while (mlist != NULL)
	{
		struct mod* temp = mlist->next;
		free(mlist->seq);
		free(mlist->inds);
		free(mlist);
		mlist = temp;
	}
	for (i=0; i<num_inds; i++)
		free(ind_names[i]);
	free(ind_names);
	tdestroy(ind_tree, free);
	pclose(seq_fp);
	pclose(mod_fp);
	return 0;
}

int int_cmp (const void *a, const void *b)
{
	if ((*((int*) a)) > (*((int*) b)))
		return 1;
	if ((*((int*) a)) < (*((int*) b)))
		return -1;

	return 0;
}

int ind_cmp (const void *a, const void *b)
{
	if (NOHET(*((int*) a)) > NOHET(*((int*) b)))
		return 1;
	if (NOHET(*((int*) a)) < NOHET(*((int*) b)))
		return -1;

	/* break ties by putting het first */
	return int_cmp(a, b);
}

/* order by smallest int */
int chrom_cmp (const void *a, const void *b)
{
	return int_cmp(((struct chrom*)a)->inds, ((struct chrom*)b)->inds);
}

int make_chroms (struct mod* mlist, struct chrom** chroms_ptr, int* max_seq_len_ptr, int num_inds)
{
	int num_chroms = 1; /* first one is already accounted for (see below) */
	struct chrom *chroms = *chroms_ptr;
	
	static int chroms_alloc = 32;

	struct mod *mlist_ptr;

	int i;

	if (chroms == NULL)
		chroms = (struct chrom *) malloc(sizeof(struct chrom) * chroms_alloc);

	/* the first chromosome always has no modifications (but may include individuals) */
	chroms[0].inds = (int *) malloc(sizeof(int)*(num_inds+1));
	chroms[0].num_inds = 1;
	chroms[0].num_mods = 0;
	chroms[0].mods_ptr = NULL;
	chroms[0].inds[0] = -1;

	/* identify all the overlapping mods */
	for (mlist_ptr=mlist; mlist_ptr->next != NULL; mlist_ptr=mlist_ptr->next)
		if (ISINORDER(mlist->chr, mlist->end+1, mlist_ptr->next->chr, mlist_ptr->next->start))
			/* a short circuit to stop after the stack is past the end of the current end; also
			 * needed to make sure we are in the right chromosome! */
			break;
		else if (OVERLAP(mlist_ptr->next->start, mlist_ptr->next->end, mlist->start, mlist->end) > 0)
		{
			struct mod* m = mlist_ptr->next;

			/* indicates whether individuals that are part of this mod have been used in another chromosome */
			unsigned char *ind_used = (unsigned char *) calloc (m->num_inds, sizeof(unsigned char));

			/* list of used list for each individual old chromosome */
			int *used_list = (int *) calloc (m->num_inds, sizeof(int));

			/* indicates the new modification is a het */
			unsigned char *used_list_het = (unsigned char *) calloc (m->num_inds, sizeof(unsigned char));

			/* update max_seq_len_ptr if the new sequence will be longer than the old sequence */
			*max_seq_len_ptr += MAX(0, ((int) strlen(m->seq)) - (m->end - m->start + 1));

			/* look at each old chromosome... we go in reverse order to avoid looking at the new ones
			 * and so we end up on the wt chromosome */
			for (i=num_chroms-1; i >= 0; i--)
			{
				int j, k;
				int num_match = 0;
				int nc = i;

				/* we require a new chromosome if there are any hets that match a chromosome or
				 * if we are considering i==0 and there are any unused inds */
				int override = 0;

				/* NOTE: find overlaps in chroms[i].inds and m->inds; note both lists are sorted by NOHET value*/
				for (j=0, k=0; j<chroms[i].num_inds; j++)
					for (; k<m->num_inds && chroms[i].inds[j] >= NOHET(m->inds[k]); k++)
						if (chroms[i].inds[j] == NOHET(m->inds[k]))
						{
							ind_used[k] = 1;
							used_list[num_match] = j;
							used_list_het[num_match] = m->inds[k] < 0;

							if (used_list_het[num_match])
								override = 1;

							num_match++;
						}

				/* we do something special when we are on chromosome 0 (which is the catch all chromosome 
				 * notice: used_list stores something fundamentally different now */
				if (i==0)
					for (j=0; j<m->num_inds; j++)
						if (!ind_used[j])
						{
							override = 1;
							break;
						}

				if (!override && num_match == 0)
					/* if no individuals overlap, we don't have to do anything for this 
					 * chromosome... move on*/
					continue;
				else if (override || num_match < chroms[i].num_inds)
				{
					/* if some but not all individuals overlap, then we create a new chromosome */
					nc = num_chroms++;

					if (num_chroms >=  chroms_alloc)
					{
						chroms_alloc *= 2;
						chroms = (struct chrom*) realloc(chroms, sizeof(struct chrom) * chroms_alloc);
					}
					
					/* we can never add new values to an old chrom (except 0) */
					chroms[nc].inds = (int *) malloc(sizeof(int)*(i==0 ? num_inds : num_match));
					chroms[nc].num_inds = 0;

					/* go through and remove the inds from the old chromosome and add to the new one
					 * do this backwards to avoid changes in the indexes in used_list */
					for (j=num_match-1; j>=0; j--)
					{
						chroms[nc].inds[chroms[nc].num_inds++] = chroms[i].inds[used_list[j]];

						/* hets are not removed from the old chromosome */
						if (!used_list_het[j])
							chroms[i].inds[used_list[j]] = chroms[i].inds[--chroms[i].num_inds];
					}

					if (i == 0)
						/* just add the new chromosomes here */
						for (j=0; j<m->num_inds; j++)
							if (!ind_used[j])
							{
								/* hets are added to both the new and old chromosomes */
								if (m->inds[j] < 0)
									chroms[i].inds[chroms[i].num_inds++] = NOHET(m->inds[j]);
								chroms[nc].inds[chroms[nc].num_inds++] = NOHET(m->inds[j]);
							}

					chroms[nc].mods_alloc = 4;
					chroms[nc].mods_ptr = NULL;
					
					/* have to copy old ones */
					chroms[nc].num_mods = chroms[i].num_mods;
				}

				chroms[nc].num_mods++;

				while (chroms[nc].mods_alloc < chroms[nc].num_mods)
					chroms[nc].mods_alloc *= 2;
				
				chroms[nc].mods_ptr = (struct mod**) realloc(chroms[nc].mods_ptr, sizeof(struct mod*) * chroms[nc].mods_alloc);
				
				/* resort the individuals; TODO: make this unnecessary by ensuring they are sorted on input */
				qsort(chroms[i].inds, chroms[i].num_inds, sizeof(int), int_cmp);

				if (i != nc)
				{
					memcpy(chroms[nc].mods_ptr, chroms[i].mods_ptr, sizeof(struct mod *) * chroms[i].num_mods);
					qsort(chroms[nc].inds, chroms[nc].num_inds, sizeof(int), int_cmp);
				}

				chroms[nc].mods_ptr[chroms[nc].num_mods-1] = m;
			}

			free(ind_used);
			free(used_list);
			free(used_list_het);
		}

	/* sort the chroms (individuals are guaranteed to be sorted) */
	qsort(chroms, num_chroms, sizeof(struct chrom), chrom_cmp);

	*chroms_ptr = chroms;

	return num_chroms;
}

struct ind_node 
{
	char *name;
	int idx;
};

int ind_node_cmp (const void *a, const void *b)
{
	return strcmp(((struct ind_node*)a)->name, ((struct ind_node*)b)->name);
}

void update_mods(FILE* fp, struct mod* mlist, char** line_ptr, size_t *line_alloc_ptr, char ***ind_names_ptr, void **ind_tree_ptr, int *num_inds_ptr)
{
	static unsigned int file_empty = 0;
	static struct mod* mlist_end = NULL;
	struct mod* mlist_ptr;
	static int ind_names_alloc = 0;

	if (mlist_end == NULL)
		mlist_end = mlist;

	/* scan through the list and eliminate unneeded entries */
	for (mlist_ptr=mlist; mlist_ptr->next != NULL; )
		if (ISINORDER(mlist_ptr->next->chr, mlist_ptr->next->end+1, mlist->chr, mlist->start))
			{
				/* this item will never be overlapped with... delete it */
				struct mod* temp = mlist_ptr->next->next;
				
				if (mlist_ptr->next == mlist_end)
					mlist_end = mlist_ptr;

				free(mlist_ptr->next->seq);
				free(mlist_ptr->next->inds);
				free(mlist_ptr->next);
				mlist_ptr->next = temp;
			}
			else if (ISINORDER(mlist->chr, mlist->start, mlist_ptr->next->chr, mlist_ptr->next->start))
				/* a short circuit to stop after the stack item's start is after the current start; 
				 * we won't be deleting any more */
				break;
			else
				mlist_ptr=mlist_ptr->next;

	/* if the file is empty we just return now */
	if (file_empty)
		return;
	
	/* read in new regions until a region that starts after the end of this region *
	 * is read in */
	while (
			/* if the list is empty... read in more */
			mlist->next == NULL || 
			
			/* or if the current region ends after our stack (i.e. the
			 * last stack item starts before one after the item we just read in)
			 * we want them out of order so we can ensure we have nothing left
			 * in fp2 that will overlap with mlist */
			ISINORDER(mlist_end->chr, mlist_end->start, mlist->chr, mlist->end)
	)
	{
		if (getline(line_ptr, line_alloc_ptr, fp) == -1)
		{
			file_empty = 1;
			return;
		}
		else
		{
			int i;

			char *lt = *line_ptr;
			mlist_end->next = (struct mod*) malloc (sizeof(struct mod));

			mlist_end->next->next = NULL;

			/* read in the new modification */
			strsep(&lt, "\t"); /* discard the first line */
			strcpy(mlist_end->next->chr, strsep(&lt, "\t"));
			mlist_end->next->start = atoll(strsep(&lt, "\t"));
			mlist_end->next->end = atoll(strsep(&lt, "\t"));

			/* throw out if this entry is already behind us */
			if (ISINORDER(mlist_end->next->chr, mlist_end->next->end+1, mlist->chr, mlist->start))
			{
				free(mlist_end->next);
				mlist_end->next = NULL;
				continue;
			}

			mlist_end->next->seq = strdup(strsep(&lt, "\t"));
			for (i=0; mlist_end->next->seq[i] != '\0'; i++)
				mlist_end->next->seq[i] = toupper(mlist_end->next->seq[i]);

			/* read the semi-colon separated list of individuals */
			mlist_end->next->num_inds = 0;
			if (lt[0] != '\0' && lt[0] != '\t')
			{
				int j, k;
				mlist_end->next->num_inds++;

				for (j=0; lt[j] != '\0' && lt[j] != '\t'; j++)
					if (lt[j] == ';')
						mlist_end->next->num_inds++;
				mlist_end->next->inds = (int *) malloc(sizeof(int)*mlist_end->next->num_inds);
				
				for (j=0; j<mlist_end->next->num_inds; j++)
				{
					/* we convert the input strings to numbers and store those */
					struct ind_node m, **mpp;
					int is_het = 0;

					m.name = strsep(&lt, ";\t\n");

					/* hets are represented with a starting ? in the name and stored as negative values in inds */
					if (m.name[0] == '?')
					{
						is_het = 1;
						m.name++;
					}

					mpp = tfind((void*) &m, ind_tree_ptr, ind_node_cmp);

					if (mpp == NULL)
					{
						struct ind_node *mp = (struct ind_node*) malloc(sizeof(struct ind_node));

						if (*num_inds_ptr >= ind_names_alloc)
						{
							if (ind_names_alloc == 0)
								ind_names_alloc = 32;
							else
								ind_names_alloc *= 2;

							*ind_names_ptr = (char**) realloc(*ind_names_ptr, sizeof(char *) * ind_names_alloc);
						}

						(*ind_names_ptr)[*num_inds_ptr] = strdup(m.name);
						mp->name = (*ind_names_ptr)[*num_inds_ptr];
						mp->idx = *num_inds_ptr;
						mpp = tsearch((void *) mp, ind_tree_ptr, ind_node_cmp);
						(*num_inds_ptr)++;
					}

					mlist_end->next->inds[j] = is_het ? HET((*mpp)->idx) : (*mpp)->idx;
				}

				/* sort and unique ind list */
				qsort(mlist_end->next->inds, mlist_end->next->num_inds, sizeof(int), ind_cmp);

				/* NOTE: if list mixes het and non-het for same individual, het is maintained */
				for (j=0, k=0; j<mlist_end->next->num_inds; j++)
					if (j==0 || NOHET(mlist_end->next->inds[j]) != NOHET(mlist_end->next->inds[j-1]))
					  mlist_end->next->inds[k++] = mlist_end->next->inds[j];
				mlist_end->next->num_inds = k;
			}

			/* make sure that the data point makes sense */
			ASSERTE(mlist_end->next->start <= mlist_end->next->end, "%s has start follow end.\n", *line_ptr);

			mlist_end = mlist_end->next;
		}
	}
}

