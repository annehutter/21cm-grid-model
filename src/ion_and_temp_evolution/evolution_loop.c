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

#define CUB(X) ((X)*(X)*(X))

void do_ion_integration_step(int dim_matrix, double X_old[], double dt, double X[], recomb_t *thisRecombRates, cell_t *thisCell, int N)
{
	for (int i=0; i<N; i++)
	{
		ModifiedPatankarRK(dim_matrix, X_old, dt, X, thisRecombRates, thisCell);
		if(X[1]/X_old[1] > 2. && N<2000)
		{
			do_ion_integration_step(dim_matrix, X_old, dt*0.5, X, thisRecombRates, thisCell, 2*N);
		}
		initializeVector(dim_matrix, X_old, X);
	}
}

void do_temp_integration_step(double *temp_old, double dt_sub, double dt, double n, double dn, cell_t *thisCell, int N)
{
	double temp;
	for (int i=0; i<N; i++)
	{
		temp = PatankarRKnoncon(*temp_old, dt_sub, dt, n, dn, thisCell);
// 		printf("temp = %e\t temp_old = %e\n", temp, *temp_old);
		if(temp/(*temp_old) > 2.)
		{
			do_temp_integration_step(temp_old, dt_sub*0.5, dt, n, dn, thisCell, N*2);
		}
		*temp_old = temp;
	}
}

void calc_step(recomb_t *thisRecombRates, cell_t *thisCell, int dim_matrix, double z, double dz)
{
	double X_old[dim_matrix], X_old_local[dim_matrix];
	double X[dim_matrix];
	
	double temp_old;
	
	double dens;
	
	double n, n_old, dn;
	
	double dt;
	
	double Y = 0.24;
	
	dt = dt_from_dz(z, dz);
	
	X_old[0] = thisCell->XHI;
	X_old[1] = thisCell->XHII;
	X_old[2] = thisCell->XHeI;
	X_old[3] = thisCell->XHeII;
	X_old[4] = thisCell->XHeIII;
	temp_old = thisCell->temp;
	dens = thisCell->dens;
	
	dens = dens*1.9e-7;
	
// 	printf("dens = %e\t temp = %e\n", dens, temp_old);
	
	initializeVector(dim_matrix, X_old_local, X_old);
	
	do_ion_integration_step(dim_matrix, X_old, dt, X, thisRecombRates, thisCell, 1);
	
	n_old = dens*CUB(1.+z+dz)*((1.-Y)*(X_old_local[0] + 2.*X_old_local[1]) + Y*(X_old_local[2] + 2.*X_old_local[3] + 3.*X_old_local[4]));
	n = dens*CUB(1.+z)*((1.-Y)*(X[0] + 2.*X[1]) + Y*(X[2] + 2.*X[3] + 3.*X[4]));
	dn = n_old - n;
	
// 	printf("n_old = %e\t n = %e\t dn = %e\n", n_old, n, dn);
	
	do_temp_integration_step(&temp_old, dt, dt, n_old, dn, thisCell, 1);
	
	thisCell->XHI = X[0];
	thisCell->XHII = X[1];
	thisCell->XHeI = X[2];
	thisCell->XHeII = X[3];
	thisCell->XHeIII = X[4];
	thisCell->temp = temp_old;
}

double dt_from_dz(double z, double dz)
{
	double ti = t_at_z(z);
	double tf = t_at_z(z-dz);
	
	return tf-ti;
}

double t_at_z(double z)
{
	const double H0 = 2.3e-18;
	const double omega_m = 0.27;
	
	return 2./(3.*H0*omega_m)*pow(1.+z,-1.5);
}