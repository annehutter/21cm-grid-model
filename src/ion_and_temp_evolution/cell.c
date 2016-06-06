#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "cell.h"

cell_t *initCell()
{
	cell_t *thisCell;
	
	thisCell = malloc(sizeof(cell_t));
	if(thisCell == NULL)
	{
		fprintf(stderr, "thisCell in initCell (cell.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	thisCell->XHI = 1.;
	thisCell->XHII = 1.e-10;
	thisCell->XHeI = 1.;
	thisCell->XHeII = 1.e-10;
	thisCell->XHeIII = 1.e-10;
	
	thisCell->temp = 0.;
	
	thisCell->heatHI = 0.;
	thisCell->heatHeI = 0.;
	thisCell->heatHeII = 0.;
	
	thisCell->photIonHI = 0.;
	thisCell->photIonHeI = 0.;
	thisCell->photIonHeII = 0.;
	
	thisCell->dens = 1.;
	
	return thisCell;
}

void deallocate_cell(cell_t *thisCell)
{
	free(thisCell);
}

cell_t *initCell_temp(double temp, double dens)
{
	cell_t *thisCell;
	
	thisCell = initCell();
	thisCell->temp = temp;
	thisCell->dens = dens;
	
	return thisCell;
}

void update_X_cell(cell_t *thisCell, double XHII, double XHeII, double XHeIII)
{
	thisCell->XHI = 1.-XHII;
	thisCell->XHII = XHII;
	thisCell->XHeI = 1.-XHeII-XHeIII;
	thisCell->XHeII = XHeII;
	thisCell->XHeIII = XHeIII;
}

void update_photIon_cell(cell_t *thisCell, double photIonHI, double photIonHeI, double photIonHeII)
{
	thisCell->photIonHI = photIonHI;
	thisCell->photIonHeI = photIonHeI;
	thisCell->photIonHeII = photIonHeII;
}

void update_dens_cell(cell_t *thisCell, double dens)
{
	thisCell->dens = dens;
}

void update_temp_cell(cell_t *thisCell, double temp)
{
	thisCell->temp = temp;
}

void update_heat_cell(cell_t *thisCell, double heatHI, double heatHeI, double heatHeII)
{
	thisCell->heatHI = heatHI;
	thisCell->heatHeI = heatHeI;
	thisCell->heatHeII = heatHeII;
}