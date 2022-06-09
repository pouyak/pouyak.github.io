/* Written by Pouya Kheradpour (pouyak@gmail.com)
 * Copyright (c) 2005-2015 Massachusetts Institute of Technology
 * Released under the MIT license
 * https://pouya.kheradpour.com/software/LICENSE.mit */
/* version 20140430 */

#ifndef PK_H
#define PK_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define STRINGIFY(x)            #x
#define STR(x)                  STRINGIFY(x)

#define CONCAT_XY2(x,y)         x##y
#define CONCAT_XY(x,y)          CONCAT_XY2(x,y)

#define STREQ(x,y)              (strcmp(x,y)==0)

#define UNUSED_VAR(VAR)         (void) VAR

#define ASSERTE(b,...)          do {if (!(b)) {fprintf(stderr, __VA_ARGS__); exit(1);}} while(0)

#define MAX(a,b)                (((a)>(b)) ? (a) : (b))
#define MIN(a,b)                (((a)<(b)) ? (a) : (b))
#define MAXCHRNAMELEN           200
#define MAXFILENAMELEN          500

#define OVERLAPS(start1,end1,start2,end2)   (MIN(end1,end2) >= MAX(start1,start2))
#define OVERLAP(start1,end1,start2,end2)    (MIN(end1,end2) - MAX(start1,start2) + 1)
#define ISINORDER(chr1, pos1, chr2, pos2)   (strcmp(chr1, chr2) < 0 || (strcmp(chr1, chr2) == 0 && (pos1) <= (pos2)))

/* open the file in a gzcat */
FILE* gzpopen(char* fn)
{
	FILE* fp;

	if (STREQ(fn, "-"))
		fp = popen("gunzip -cf</dev/stdin", "r");
	else
	{
		char *f = (char*) malloc(sizeof(char) * (strlen(fn) + 14));
		ASSERTE(access(fn, R_OK) != -1, "Cannot read: %s\n", fn);
		strcpy(f, "gunzip -cf<'");
		strcat(f, fn);
		strcat(f, "'");
		fp = popen(f, "r");
		free(f);
	}

	ASSERTE(fp != NULL, "Cannot open: %s\n", fn);

	return fp;
}

#define REV(dest,src)	do {switch (src)  \
{  \
	case 'a': (dest) = 't';	break;  \
	case 't': (dest) = 'a';	break;  \
	case 'g': (dest) = 'c';	break;  \
	case 'c': (dest) = 'g';	break;  \
	case 'A': (dest) = 'T';	break;  \
	case 'T': (dest) = 'A';	break;  \
	case 'G': (dest) = 'C';	break;  \
	case 'C': (dest) = 'G';	break;  \
	default:  (dest) = src;  \
} } while(0)

/* reverse complements seq in line */
void rev_comp(char* seq, unsigned int len)
{
	unsigned int i;
	char temp;
	for (i=0; i<len/2; i++)
	{
		temp = seq[i];

		REV(seq[i], seq[len-i-1]);
		REV(seq[len-i-1], temp);
	}

	if (len%2 == 1)
		REV(seq[len/2], seq[len/2]);
}

/* generates random number from [0,n) */
int randn(int n)
{
	/* rand produces [0,RAND_MAX] */

	/* apparently current implementations of rand have equally "random" 
	 * higher and lower bits so we can just use modulus */

	/* compute floor((RAND_MAX+1) / n) * n - 1 but avoid potential overflow 
	 * which would happen if the number is a power of 2 */
	int rm = (((RAND_MAX - n) + 1) / n) * n + (n - 1);
	int r;

	ASSERTE((n-1) <= RAND_MAX, "randn: %d out of range!", n);

	while ((r = rand()) > rm);

	return r % n;
}

#endif /* PK_H */

