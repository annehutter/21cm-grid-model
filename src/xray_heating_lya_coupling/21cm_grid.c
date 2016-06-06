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
#include "21cm_grid.h"

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR 21CM GRID */
/*-------------------------------------------------------------------------------------*/

grid_21cm_t *init21cmgrid()
{
	grid_21cm_t *this21cmGrid;
	
	this21cmGrid = malloc(sizeof(grid_21cm_t));
	
	this21cmGrid->nbins = 0;
	this21cmGrid->box_size = 0.;
	
	this21cmGrid->temp = NULL;
	this21cmGrid->dens = NULL;
	this21cmGrid->Xe = NULL;
	this21cmGrid->XHI = NULL;
	this21cmGrid->XHeI = NULL;
	this21cmGrid->XHeII = NULL;
	this21cmGrid->Ts_inv = NULL;
	this21cmGrid->Tb = NULL;
	
	this21cmGrid->local_n0 = 0;
	this21cmGrid->local_0_start = 0;
	
	return this21cmGrid;
}

void deallocate_21cmgrid(grid_21cm_t *this21cmGrid)
{
	if(this21cmGrid->temp == NULL) fftw_free(this21cmGrid->temp);
	if(this21cmGrid->dens == NULL) fftw_free(this21cmGrid->dens);
	if(this21cmGrid->Xe == NULL) fftw_free(this21cmGrid->Xe);
	if(this21cmGrid->XHI == NULL) fftw_free(this21cmGrid->XHI);
	if(this21cmGrid->XHeI == NULL) fftw_free(this21cmGrid->XHeI);
	if(this21cmGrid->XHeII == NULL) fftw_free(this21cmGrid->XHeII);
	if(this21cmGrid->Ts_inv == NULL) fftw_free(this21cmGrid->Ts_inv);
	if(this21cmGrid->Tb == NULL) fftw_free(this21cmGrid->Tb);

	free(this21cmGrid);
}

grid_21cm_t *allocate_21cmgrid(int nbins, float box_size)
{
#ifdef __MPI
	ptrdiff_t alloc_local, local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	grid_21cm_t *this21cmGrid;

	this21cmGrid = init21cmgrid();
	
	this21cmGrid->nbins = nbins;
	this21cmGrid->box_size = box_size;
	
	this21cmGrid->local_n0 = nbins;
	this21cmGrid->local_0_start = 0;
#ifdef __MPI	
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	
	this21cmGrid->local_n0 = local_n0;
	this21cmGrid->local_0_start = local_0_start;
	
	this21cmGrid->temp = fftw_alloc_complex(alloc_local);
	this21cmGrid->dens = fftw_alloc_complex(alloc_local);
	this21cmGrid->Xe = fftw_alloc_complex(alloc_local);
	this21cmGrid->XHI = fftw_alloc_complex(alloc_local);
	this21cmGrid->XHeI = fftw_alloc_complex(alloc_local);
	this21cmGrid->XHeII = fftw_alloc_complex(alloc_local);
	this21cmGrid->Ts_inv = fftw_alloc_complex(alloc_local);
	this21cmGrid->Tb = fftw_alloc_complex(alloc_local);
#else
	this21cmGrid->temp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	this21cmGrid->dens = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	this21cmGrid->Xe = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	this21cmGrid->XHI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	this21cmGrid->XHeI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	this21cmGrid->XHeII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	this21cmGrid->Ts_inv = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	this21cmGrid->Tb = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
#endif
	
	initialize_21cmgrid(this21cmGrid->temp, nbins, local_n0, 0.);
	initialize_21cmgrid(this21cmGrid->dens, nbins, local_n0, 1.);
	initialize_21cmgrid(this21cmGrid->Xe, nbins, local_n0, 0.);
	initialize_21cmgrid(this21cmGrid->XHI, nbins, local_n0, 1.);
	initialize_21cmgrid(this21cmGrid->XHeI, nbins, local_n0, 1.);
	initialize_21cmgrid(this21cmGrid->XHeII, nbins, local_n0, 0.);
	initialize_21cmgrid(this21cmGrid->Ts_inv, nbins, local_n0, 0.);
	initialize_21cmgrid(this21cmGrid->Tb, nbins, local_n0, 0.);

	return this21cmGrid;
}

void initialize_21cmgrid(fftw_complex *thisArray, int nbins, int local_n0, double value)
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


void read_density_21cmgrid(grid_21cm_t *this21cmGrid, char *filename, int double_precision)
{
#ifdef __MPI
	ptrdiff_t local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	int nbins;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	
	if(double_precision == 1)
	{
#ifdef __MPI
	local_0_start = thisGrid->local_0_start;
	read_grid_doubleprecision(this21cmGrid->dens, nbins, local_n0, local_0_start, filename);
#else
	read_grid_doubleprecision(this21cmGrid->dens, nbins, local_n0, filename);
#endif
	}else{
#ifdef __MPI
	local_0_start = thisGrid->local_0_start;
	read_grid(this21cmGrid->dens, nbins, local_n0, local_0_start, filename);
#else
	read_grid(this21cmGrid->dens, nbins, local_n0, filename);
#endif
	}
}

double get_mean_Xe_21cmgrid(grid_21cm_t *this21cmGrid)
{
	int nbins = this21cmGrid->nbins;
	int local_n0 = this21cmGrid->local_n0;
	double result;
	
	double sum = 0.;
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				sum += this21cmGrid->Xe[i*nbins*nbins+j*nbins+k];
			}
		}
	}
#ifdef __MPI
	MPI_Allreduce(&sum, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	result /= (nbins*nbins*nbins);
	
	return result;
}
