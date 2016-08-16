#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "outputs.h"

outputlist_t *read_outputlist(char *filename)
{
	FILE * file;
	
	int counter;
	double redshift;
	char line[128];
	
    outputlist_t *thisOutputlist;
    
	if(file_exist(filename) == 1)
	{
		file = fopen(filename, "rt");
		counter = 0;
        while(fgets(line, 128, file) != NULL) counter++;
        
        thisOutputlist = initOutputlist(counter);
        
        fseek(file, 0, SEEK_SET);
        
		while(fgets(line, 128, file) != NULL)
		{
		      /* get a line, up to 80 chars from fr.  done if NULL */
		      sscanf (line, "%le", &redshift);

              thisOutputlist->redshift[counter] = redshift;
              
		      counter ++;
		      if(counter > thisOutputlist->num)
		      {
				fprintf(stderr, "Error in reading outputlist (inputs.c).\n");
				exit(EXIT_FAILURE);
		      }
		}
		assert(counter == thisOutputlist->num);
		fclose(file);  /* close the file prior to exiting the routine */  
		return thisOutputlist;
	}else{
		return  NULL;
	}
}

outputlist_t *initOutputlist(int num)
{
	outputlist_t *thisOutputlist;
	
	thisOutputlist = malloc(sizeof(outputlist_t));
	if(thisOutputlist == NULL)
	{
		fprintf(stderr, "thisOutputlist in initOutputlist (inputs.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	thisOutputlist->num = num;
    
    thisOutputlist->redshift = malloc(sizeof(double)*num);
	
	return thisOutputlist;
}

void deallocateOutputlist(outputlist_t *thisOutputlist)
{
	if(thisOutputlist != NULL) 
    {
        if(thisOutputlist->redshift != NULL) free(thisOutputlist->redshift);
        
        free(thisOutputlist);
    }
}
