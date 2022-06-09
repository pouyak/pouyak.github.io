/* Written by Pouya Kheradpour (pouyak@gmail.com)
 * Copyright (c) 2005-2015 Massachusetts Institute of Technology
 * Released under the MIT license
 * https://pouya.kheradpour.com/software/LICENSE.mit */
/* version 20140120 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pk.h"

/********* INPUT MUST BE SORTED BY CHR THEN START!! ***********/
/* e.g. first pipe through sort -k2,2 -k3,3n ***/
/* input file must be of the form: type chr start stop */
/* command line arguments should be in order list of type precedence
 * (best precedence first) */

#define MAXLINELEN		1000

struct end_list {
	unsigned int type;
	char chr[MAXCHRNAMELEN+1];
	unsigned int loc;
	struct end_list* next;
};

/* the print function doesn't actually print until a region comes in that is unlike
 * the previous, this way it can merge similar regions */
unsigned int last_region_type;
char last_region_chr[MAXCHRNAMELEN+1] = "";
unsigned int last_region_start;
unsigned int last_region_end;

void print_region(char **type_list, unsigned int num_type, unsigned int type, char* chr, unsigned int end);
unsigned int type(char **type_list, unsigned int num_type, char *str);

void print_region(char **type_list, unsigned int num_type, unsigned int type, char* chr, unsigned int end)
{
	/* empty regions aren't kept... */
	if (end == last_region_end)
		return;

	if (STREQ(chr, last_region_chr) && type == last_region_type)
		/* merge the two regions, if possible */
		last_region_end = end;
	else
	{
		/* print the previous region and start a new one */
		if (last_region_type != num_type)
			printf("%s\t%s\t%d\t%d\n", type_list[last_region_type], last_region_chr, last_region_start, last_region_end);

		last_region_type = type;
		strcpy(last_region_chr, chr);
		last_region_start = last_region_end+1;
		last_region_end = end;
	}
}

unsigned int type(char **type_list, unsigned int num_type, char *str)
{
	unsigned int i;
	for (i=0; i<num_type; i++)
		if (STREQ(str, type_list[i]))
			return i;
	return num_type;
}

int main (int argc, char** argv)
{
	char **type_list = argv + 1;
	unsigned int num_type = argc-1;
	unsigned int *type_counts = (unsigned int*) calloc(num_type, sizeof(unsigned int));

	char last_chr[MAXCHRNAMELEN+1] = "";
	char new_chr[MAXCHRNAMELEN+1];
	char new_line[MAXLINELEN+1];
	char new_type_name[MAXLINELEN+1];
	unsigned int last_start = 0;
	unsigned int new_start;
	unsigned int new_end; 
	unsigned int new_type;
	unsigned int best_cur_type = num_type;
	last_region_type = num_type;

	unsigned char at_eof = 0;

	/* first one is a dummy to simplify code a bit */
	struct end_list* ends = (struct end_list*) malloc(sizeof(struct end_list));
	struct end_list* ends_ptr;
	ends->next = NULL;
	
	while (!at_eof)
	{
		at_eof = scanf("%" STR(MAXLINELEN) "[^\n]\n", new_line) == EOF;
		
		if (!at_eof)
		{
			sscanf(new_line, "%" STR(MAXLINELEN) "s %" STR(MAXCHRNAMELEN) "s %u %u", new_type_name, new_chr,  &new_start, &new_end);
			ASSERTE(ISINORDER(last_chr, last_start, new_chr, new_start), "%s is out of order.\n", new_line);

			strcpy(last_chr, new_chr);
			last_start = new_start;
		}

		/* first go through and print out any regions that ended before this start */
		/* notice we keep a pointer to the previous in the list to facilitate easy removal */
		for (ends_ptr=ends; ends_ptr->next!=NULL && (at_eof || ISINORDER(ends_ptr->next->chr, ends_ptr->next->loc, new_chr, new_start-1));)
		{
			struct end_list* cur = ends_ptr->next;
			/* print region to the end */
			print_region(type_list, num_type, best_cur_type, cur->chr, cur->loc);

			/***** remove the region *****/
			/* reduce the count of the region of this type by one */
			type_counts[cur->type]--;

			/* update best_cur_type, as needed */
			for (;best_cur_type < num_type && type_counts[best_cur_type] == 0; best_cur_type++);

			/* remove the end from the list */
			ends_ptr->next = ends_ptr->next->next; 
			free(cur);
		}

		/**** add the new region to the list and end the previous region ****/
		if (!at_eof)
		{
			struct end_list* new_end_item = (struct end_list*) malloc(sizeof(struct end_list));

			new_type = type(type_list, num_type, new_type_name);

			/* not in type list... just ignore this new item */
			if (new_type == num_type)
				continue;

			/* print out region ending before this start */
			print_region(type_list, num_type, best_cur_type, new_chr, new_start-1);

			/* add this start to type list */
			type_counts[new_type]++;
			if (new_type < best_cur_type)
				best_cur_type = new_type;

			/* add end to end list */
			strcpy(new_end_item->chr, new_chr);
			new_end_item->type = new_type;
			new_end_item->loc = new_end;
			for (ends_ptr=ends; ends_ptr->next!=NULL && ISINORDER(ends_ptr->next->chr, ends_ptr->next->loc, new_chr, new_end); ends_ptr = ends_ptr->next);
			new_end_item->next = ends_ptr->next;
			ends_ptr->next = new_end_item;
		}
	}

	/* flush the last region */
	print_region(type_list, num_type, num_type, "", last_region_end+1);
	free(ends);
	free(type_counts);

	return 0;
}

