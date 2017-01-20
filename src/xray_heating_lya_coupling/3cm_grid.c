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
#include "3cm_grid.h"

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR 21CM GRID */
/*-------------------------------------------------------------------------------------*/

grid_3cm_t *init3cmgrid()
{
	grid_3cm_t *this3cmGrid;
	
	this3cmGrid = malloc(sizeof(grid_3cm_t));
	
	this3cmGrid->nbins = 0;
	this3cmGrid->box_size = 0.;
	
	this3cmGrid->z = 0.;
	this3cmGrid->z_old = 0.;
	
	this3cmGrid->temp_old = 0.;
	
	this3cmGrid->temp = 0.;
	this3cmGrid->Ts_inv = NULL;
	this3cmGrid->Tb = NULL;
	
	this3cmGrid->local_n0 = 0;
	this3cmGrid->local_0_start = 0;
	
	return this3cmGrid;
}

void deallocate_3cmgrid(grid_3cm_t *this3cmGrid)
{
	if(this3cmGrid->Ts_inv != NULL) fftw_free(this3cmGrid->Ts_inv);
	if(this3cmGrid->Tb != NULL) fftw_free(this3cmGrid->Tb);

	free(this3cmGrid);
}

grid_3cm_t *allocate_3cmgrid(int nbins, float box_size)
{
#ifdef __MPI
	ptrdiff_t alloc_local, local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	grid_3cm_t *this3cmGrid;

	this3cmGrid = init3cmgrid();
	
	this3cmGrid->nbins = nbins;
	this3cmGrid->box_size = box_size;
	
	this3cmGrid->local_n0 = nbins;
	this3cmGrid->local_0_start = 0;
    local_n0 = nbins;
#ifdef __MPI	
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	
	this3cmGrid->local_n0 = local_n0;
	this3cmGrid->local_0_start = local_0_start;
	
	this3cmGrid->Ts_inv = fftw_alloc_complex(alloc_local);
	this3cmGrid->Tb = fftw_alloc_complex(alloc_local);
#else
	this3cmGrid->Ts_inv = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	this3cmGrid->Tb = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
#endif
	
	initialize_3cmgrid(this3cmGrid->Ts_inv, nbins, local_n0, 0.);
	initialize_3cmgrid(this3cmGrid->Tb, nbins, local_n0, 0.);

	return this3cmGrid;
}

void initialize_3cmgrid(fftw_complex *thisArray, int nbins, int local_n0, double value)
{
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				thisArray[i*nbins*nbins+j*nbins+k] = value +0.*I;
			}
		}
	}
}

void set_temperature_3cmgrid(grid_3cm_t *this3cmGrid, double value)
{
    this3cmGrid->temp = value;
}

void write_Tb_3cm_field_file(grid_3cm_t *this3cmGrid, char *filename)
{
	int nbins = this3cmGrid->nbins;
	int local_n0 = this3cmGrid->local_n0;
	
#ifdef __MPI
	int local_0_start = this3cmGrid->local_0_start;
	
	write_grid_to_file_float(this3cmGrid->Tb, nbins, local_n0, local_0_start, filename);
#else
	write_grid_to_file_float(this3cmGrid->Tb, nbins, local_n0, filename);
#endif
}

void write_Ts_3cm_field_file(grid_3cm_t *this3cmGrid, char *filename)
{
	int nbins = this3cmGrid->nbins;
	int local_n0 = this3cmGrid->local_n0;
	
#ifdef __MPI
	int local_0_start = this3cmGrid->local_0_start;
	
	write_grid_to_file_float(this3cmGrid->Ts_inv, nbins, local_n0, local_0_start, filename);
#else
	write_grid_to_file_float(this3cmGrid->Ts_inv, nbins, local_n0, filename);
#endif
}
