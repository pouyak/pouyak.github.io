/* Written by Pouya Kheradpour (pouyak@gmail.com)
 * Copyright (c) 2005-2015 Massachusetts Institute of Technology
 * Released under the MIT license
 * http://compbio.mit.edu/pouyak/software/LICENSE.mit */
/* version 20160416 */

#ifndef PK_BL_H
#define PK_BL_H

/* indexing for consm->mat. note similarity to ICidx */
#define CMidx(sp,edge)		((num_sp-1) * (edge) + ((sp)-1))

/* structure holding data for determining acceptance and conservation level of motifs */
struct consm
{
	/* allowed species combinations */
	unsigned char* sp;

	/* branch length of species combinations */
	float* bl;

	/* data used to compute bl in its absence, -1 indicates root */
	unsigned char* mat;
	float* len;
	unsigned int* tlist; /* array used to store list of matching species */
};

/* reads a branch length file where the first column is the species bistring and the
 * forth column is the total branch length */
float* read_bl_file(char* bl_file, unsigned char num_sp)
{
	unsigned int num_sp_comb = 1 << (num_sp-1);
	float* consm_bl = (float*) calloc(num_sp_comb, sizeof(float));

	unsigned int index;
	float bl;

	FILE* fp = fopen(bl_file, "r");
	ASSERTE(fp != NULL, "Cannot open: %s\n", bl_file);

	while (fscanf(fp, "%d %*s %*s %f\n", &index, &bl) != EOF)
		consm_bl[index]=bl;
	fclose(fp);

	return consm_bl;
}

/* reads a branch length file in parent tree format where the second column indicates the parent node (numbered from 0) and the third column indicates the length of the edge to that node */
void read_blp_file(char* bl_file, unsigned char num_sp, struct consm *consm)
{
	int i;
	unsigned char sp;

	int num_edges = 2*num_sp - 2;
	FILE* fp = fopen(bl_file, "r");
	ASSERTE(fp != NULL, "Cannot open: %s\n", bl_file);

	int* prt = (int*) malloc(sizeof(int)*num_edges);
	consm->mat = (unsigned char*) calloc((unsigned int) num_edges*(num_sp-1), sizeof(unsigned char));
	consm->len = (float*) malloc(sizeof(float)*num_edges);

	/* not used here... allocated here once */
	consm->tlist = (unsigned int*) malloc(sizeof(unsigned int)*(num_sp-1));

	for (i=0; i<num_edges; i++)
		fscanf(fp, "%*s %d %f\n", &prt[i], &consm->len[i]);
	
	/* for each path from target to other species, indicate which edges are needed */
	for (sp=1; sp<num_sp; sp++)
	{
		/* go from 0 to root, marking everything as 1 */
		for (i=0; i != -1; i = prt[i])
			consm->mat[CMidx(sp,i)] = 1;

		/* go from sp to root, switching 1<->0  */
		for (i=sp; i != -1; i = prt[i])
			consm->mat[CMidx(sp,i)] = 1 - consm->mat[CMidx(sp,i)];
	}
	fclose(fp);
	free(prt);
}

float get_bl (unsigned char num_sp, struct consm *consm, long long int match_bitstr)
{
	if (consm->bl == NULL && consm->mat == NULL)
		return -1; 
	else if (consm->bl != NULL)
		return consm->bl[match_bitstr];
	else
	{
		float bl = 0;
		int num_edges = 2*num_sp - 2;
		long long int i;
		int num_sp_match = 0, j;

		/* create list of matching species */
		for (i=1; match_bitstr; i++, match_bitstr >>= 1)
			if (match_bitstr & 1)
				consm->tlist[num_sp_match++] = i;

		if (num_sp_match > 0)
			for (i=0; i<num_edges; i++)
				/* go in reverse order because later species more likely to have more edges */
				for (j=num_sp_match-1; j>=0; j--)
					if (consm->mat[CMidx(consm->tlist[j],i)])
					{
						bl += consm->len[i];

						/* short circuit to avoid double counting */
						break;
					}
		return bl;
	}
}

#endif /* PK_BL_H */

