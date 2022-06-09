/* Written by Pouya Kheradpour (pouyak@gmail.com)
 * Copyright (c) 2005-2015 Massachusetts Institute of Technology
 * Released under the MIT license
 * https://pouya.kheradpour.com/software/LICENSE.mit */
/* version 20160416 */

#ifndef PK_RFILE_H
#define PK_RFILE_H

struct rfile
{
	char *name;
	char chr[MAXCHRNAMELEN];
	int start;
	int end;
};

/* compare function for files (first by chromosome, then by start) */
int rfile_cmp (const void *a, const void *b)
{
	int chr_cmp = strcmp(((struct rfile*) a)->chr, ((struct rfile*) b)->chr);

	if (chr_cmp != 0)
		return chr_cmp;

	if (((struct rfile*) a)->start < ((struct rfile*) b)->start)
		return -1;

	if (((struct rfile*) a)->start > ((struct rfile*) b)->start)
		return 1;

	return 0;
}

/* reads in all the file names, returns the number found */
int read_rfile(char* input_file, struct rfile** files_ptr)
{
	int num_files, num_files_alloc = 10;
	struct rfile *files = (struct rfile*) malloc(num_files_alloc * sizeof(struct rfile));
	FILE* fp = fopen(input_file, "r");
	ASSERTE(fp != NULL, "Cannot open: %s\n", input_file);

	for (num_files=0; fscanf(fp, "%s %d %d %as", files[num_files].chr, &files[num_files].start, &files[num_files].end, &files[num_files].name) != EOF; num_files++)
	{
		if ((num_files+1) >= num_files_alloc)
		{
			num_files_alloc *= 2;
			files = (struct rfile*) realloc(files, num_files_alloc * sizeof(struct rfile));
		}
	}

	fclose(fp);

	/* clean up unused space */
	files = (struct rfile*) realloc(files, num_files * sizeof(struct rfile));

	/* sort files by chromosome then start */
	qsort (files, num_files, sizeof(struct rfile), rfile_cmp);

	*files_ptr = files;
	return num_files;
}

#endif /* PK_RFILE_H */

