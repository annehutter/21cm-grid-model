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
#include "recomb_rates.h"
#include "evolution_loop.h"
#include "solve_ionfraction.h"
#include "solve_temperature.h"
#include "evolution_loop.h"

double TCMB(double z)
{
	return 2.73*(1.+z);
}

int main (int argc, char *argv[])
{
	double zstart = 300.;
	double zend = 10.;
	double dz = 0.01;
  
	cell_t *cell;
	cell = initCell_temp(TCMB(zstart), 1.);
	
	recomb_t * recomb_rates;
	recomb_rates = calcRecRate();
	
	int dim_matrix = 5;
	
	
	
	for(int i=0; i<1; i++)
	{		
		for(double z = zstart; z>zend; z = z-dz)
		{
			if(z<15.) update_photIon_cell(cell, 3.e-19/(1.+z)*(1.+zstart), 1.e-19/(1.+z)*(1.+zstart), 0.);
			printf("photHI = %e\n", cell->photIonHI);

			calc_step(recomb_rates, cell, dim_matrix, z, dz);
			printf("z = %e:\tT = %e\t T = T0/(1+z)^2 = %e\t XHII = %e\t XHeII = %e\t XHeIII = %e\n", z, cell->temp, TCMB(zstart)/((1.+zstart)*(1.+zstart))*((1.+z-dz)*(1.+z-dz)), cell->XHII, cell->XHeII, cell->XHeIII);
		}
	}
	  
	deallocate_recomb(recomb_rates);
	deallocate_cell(cell);
}