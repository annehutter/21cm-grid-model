#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <gsl/gsl_integration.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "phys_const.h"
#include "redshift_tools.h"

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

/*-------------------------------------------------------------------------------------*/
/* Hubble constant redshift z----------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/

double z_hubble(double h, double omega_m, double omega_l, double z)
{
	double Hubble = h*100; //in km s^-1 Mpc^-1
	return Hubble*sqrt(omega_m*pow(1.+z,3)+omega_l);
}

/*-------------------------------------------------------------------------------------*/
/* time between two redshifts for a flat Universe -------------------------------------*/
/*-------------------------------------------------------------------------------------*/

double z_time_from_redshiftinterval_flatuniverse(double h, double omega_m, double omega_l, double zmin, double zmax)
{
	double prefactor = 2./(3*h*100.*km_cm/Mpc_cm*sqrt(omega_l));
	double tmp = sqrt(omega_l/omega_m);
	
	return prefactor*(asinh(tmp*pow(1.+zmin, -1.5)) - asinh(tmp*pow(1.+zmax, -1.5)));
}

/*-------------------------------------------------------------------------------------*/
/* computes 1+z', whereas z' is redshift that has a comoving distance x from redshift z */
/*-------------------------------------------------------------------------------------*/

double z_distance_to_redshift(double z, double x, double h, double omega_m)
{
	return pow(1./sqrt(1.+z)-h*100.*km_cm/Mpc_cm*sqrt(omega_m)/(2.*clight_cm)*x ,-2);
}

/*-------------------------------------------------------------------------------------*/
/* global number densities ------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/

double z_calc_numdensity_per_mp_at_redshift(double h, double omega_b, double z)
{
	return 3.*SQR(h*100.*km_cm/Mpc_cm)/(8.*M_PI*G)/mp_g*omega_b*CUB(1.+z);
}

double z_calc_numdensity_H_at_redshift(double h, double omega_b, double z, double Y)
{
	return 3.*SQR(h*100.*km_cm/Mpc_cm)/(8.*M_PI*G)/mp_g*omega_b*CUB(1.+z)*(1.-Y);
}

double z_calc_numdensity_He_at_redshift(double h, double omega_b, double z, double Y)
{
	return 3.*SQR(h*100.*km_cm/Mpc_cm)/(8.*M_PI*G)/mp_g*omega_b*CUB(1.+z)*Y*0.25;
}