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

#include "utils.h"
#include "inputs.h"

inputlist_t *read_inputlist(char *filename)
{
	FILE * file;
	
	int counter;
	double redshift;
    int dens_snap, xray_snap, lya_snap;
    double lumX, alphaX, nuX_min;
    double lumLya, alphaLya, nuLya_min;
	char line[128];
    
    inputlist_t *thisInputlist;
	
	if(file_exist(filename) == 1)
	{
		file = fopen(filename, "rt");
		counter = 0;
        while(fgets(line, 128, file) != NULL) counter++;
        
        thisInputlist = initInputlist(counter);
        
        fseek(file, 0, SEEK_SET);
        
        counter = 0;
		while(fgets(line, 128, file) != NULL)
		{
		      /* get a line, up to 80 chars from fr.  done if NULL */
		      sscanf (line, "%le\t%d\t%d\t%le\t%le\t%le\t%d\t%le\t%le\t%le", &redshift, &dens_snap, &xray_snap, &lumX, &alphaX, &nuX_min, &lya_snap, &lumLya, &alphaLya, &nuLya_min);

              thisInputlist->redshift[counter] = redshift;
              thisInputlist->dens_snap[counter] = dens_snap;
              thisInputlist->xray_snap[counter] = xray_snap;
              thisInputlist->lya_snap[counter] = lya_snap;
              thisInputlist->lumX[counter] = lumX;
              thisInputlist->alphaX[counter] = alphaX;
              thisInputlist->nuX_min[counter] = nuX_min;
              thisInputlist->lumLya[counter] = lumLya;
              thisInputlist->alphaLya[counter] = alphaLya;
              thisInputlist->nuLya_min[counter] = nuLya_min;
              
		      counter ++;
		      if(counter > thisInputlist->num)
		      {
				fprintf(stderr, "Error in reading inputlist (inputs.c).\n");
				exit(EXIT_FAILURE);
		      }
		}
		assert(counter == thisInputlist->num);
		fclose(file);  /* close the file prior to exiting the routine */  
		return thisInputlist;
	}else{
		return  NULL;
	}
}

inputlist_t *initInputlist(int num)
{
	inputlist_t *thisInputlist;
	
	thisInputlist = malloc(sizeof(inputlist_t));
	if(thisInputlist == NULL)
	{
		fprintf(stderr, "thisInputlist in initInputlist (inputs.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	thisInputlist->num = num;
    
    thisInputlist->redshift = malloc(sizeof(double)*num);
    thisInputlist->dens_snap = malloc(sizeof(int)*num);
    thisInputlist->xray_snap = malloc(sizeof(int)*num);
    thisInputlist->lya_snap = malloc(sizeof(int)*num);
    thisInputlist->lumX = malloc(sizeof(double)*num);
    thisInputlist->alphaX = malloc(sizeof(double)*num);
    thisInputlist->nuX_min = malloc(sizeof(double)*num);
    thisInputlist->lumLya = malloc(sizeof(double)*num);
    thisInputlist->alphaLya = malloc(sizeof(double)*num);
    thisInputlist->nuLya_min = malloc(sizeof(double)*num);
	
	return thisInputlist;
}

void deallocateInputlist(inputlist_t *thisInputlist)
{
	if(thisInputlist != NULL) 
    {
        if(thisInputlist->redshift != NULL) free(thisInputlist->redshift);
        if(thisInputlist->dens_snap != NULL) free(thisInputlist->dens_snap);
        if(thisInputlist->xray_snap != NULL) free(thisInputlist->xray_snap);
        if(thisInputlist->lya_snap != NULL) free(thisInputlist->lya_snap);
        if(thisInputlist->lumX != NULL) free(thisInputlist->lumX);
        if(thisInputlist->alphaX != NULL) free(thisInputlist->alphaX);
        if(thisInputlist->nuX_min != NULL) free(thisInputlist->nuX_min);
        if(thisInputlist->lumLya != NULL) free(thisInputlist->lumLya);
        if(thisInputlist->alphaLya != NULL) free(thisInputlist->alphaLya);
        if(thisInputlist->nuLya_min != NULL) free(thisInputlist->nuLya_min);

        free(thisInputlist);
    }
}
