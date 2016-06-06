#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>	//included because M_PI is not defined in <math.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "fftw_array_tools.h"

/*-------------------------------------------------------------------------------------*/
/* Reading files to fftw arrays */
/*-------------------------------------------------------------------------------------*/


#ifdef __MPI
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, char *filename)
#endif
{
	float *tmparray;
		
	tmparray = (float*)malloc(sizeof(float)*local_n0*nbins*nbins);
	if(tmparray == NULL)
	{
		fprintf(stderr, "tmparray in read_grid (grid.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
#ifdef __MPI
	int success;
	int resultlen;
	char msg[MPI_MAX_ERROR_STRING];

	MPI_File mpifile;
	MPI_Offset offset;
	MPI_Status status;
	
	offset = (local_0_start*nbins*nbins*sizeof(float));
	
	success = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &mpifile);
	if(success != MPI_SUCCESS)
	{
		MPI_Error_string(success, msg, &resultlen);
		fprintf(stderr, "MPI_File_open(): %s\n", msg);
		exit(-1);
	}
	MPI_File_read_at_all(mpifile,offset,tmparray, local_n0*nbins*nbins,MPI_FLOAT,&status);
	MPI_File_close(&mpifile);
#else
	FILE *fp;
	fp = fopen(filename, "rb");
	fread(tmparray, sizeof(float), nbins*nbins*nbins, fp);
	fclose(fp);
#endif
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				toThisArray[i*nbins*nbins+j*nbins+k] = (double)tmparray[i*nbins*nbins+j*nbins+k]+0.*I;
			}
		}
	}
	free(tmparray);
}


#ifdef __MPI
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, char *filename)
#endif
{
	double *tmparray;
	
	tmparray = (double*)malloc(sizeof(double)*local_n0*nbins*nbins);
	if(tmparray == NULL)
	{
		fprintf(stderr, "tmparray in read_grid (grid.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
#ifdef __MPI
	int success;
	int resultlen;
	char msg[MPI_MAX_ERROR_STRING];

	MPI_File mpifile;
	MPI_Offset offset;
	MPI_Status status;
	
	offset = (local_0_start*nbins*nbins*sizeof(double));
	
	success = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &mpifile);
	if(success != MPI_SUCCESS)
	{
		MPI_Error_string(success, msg, &resultlen);
		fprintf(stderr, "MPI_File_open(): %s\n", msg);
		exit(-1);
	}
	MPI_File_read_at_all(mpifile,offset,tmparray, local_n0*nbins*nbins,MPI_DOUBLE,&status);
	MPI_File_close(&mpifile);
#else
	FILE *fp;
	fp = fopen(filename, "rb");
	fread(tmparray, sizeof(double), nbins*nbins*nbins, fp);
	fclose(fp);
#endif
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				toThisArray[i*nbins*nbins+j*nbins+k] = tmparray[i*nbins*nbins+j*nbins+k]+0.*I;
			}
		}
	}
	free(tmparray);
}

/*-------------------------------------------------------------------------------------*/
/* Writing fftw arrays to file */
/*-------------------------------------------------------------------------------------*/

#ifdef __MPI
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, char *filename)
#endif
{
	float *tmparray;
	
	tmparray = (float*)malloc(sizeof(float)*local_n0*nbins*nbins);
	if(tmparray == NULL)
	{
		fprintf(stderr, "tmparray in write_grid_to_file_float (grid.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				tmparray[i*nbins*nbins+j*nbins+k] = (float)creal(thisArray[i*nbins*nbins+j*nbins+k]);
			}
		}
	}
		
#ifdef __MPI
	MPI_File mpifile;
	MPI_Offset offset;
	MPI_Status status;
	
	offset = (local_0_start*nbins*nbins*sizeof(float));

	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
	MPI_File_write_at_all(mpifile, offset, tmparray, local_n0*nbins*nbins, MPI_FLOAT, &status);
	MPI_File_close(&mpifile);
#else
	FILE * fp;
	
	fp = fopen(filename, "wb");
	fwrite(tmparray, sizeof(float), nbins*nbins*nbins, fp);
	fclose(fp);
#endif
	
	free(tmparray);
}

#ifdef __MPI
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename)
#else
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, char *filename)
#endif
{
	double *tmparray;
	
	tmparray = (double*)malloc(sizeof(double)*local_n0*nbins*nbins);
	if(tmparray == NULL)
	{
		fprintf(stderr, "tmparray in write_grid_to_file_double (grid.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				tmparray[i*nbins*nbins+j*nbins+k] = (double)creal(thisArray[i*nbins*nbins+j*nbins+k]);
			}
		}
	}
	
#ifdef __MPI
	MPI_File mpifile;
	MPI_Offset offset;
	MPI_Status status;
	
	offset = (local_0_start*nbins*nbins*sizeof(double));
	
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
	MPI_File_write_at_all(mpifile, offset, tmparray, local_n0*nbins*nbins, MPI_DOUBLE, &status);
	MPI_File_close(&mpifile);
#else
	FILE * fp;
	
	fp = fopen(filename, "wb");
	fwrite(tmparray, sizeof(double), nbins*nbins*nbins, fp);
	fclose(fp);
#endif
	
	free(tmparray);
}