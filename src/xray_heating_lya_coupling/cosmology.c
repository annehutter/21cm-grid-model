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

#include "phys_const.h"
#include "cosmology.h"

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

cosmology_t *initCosmology()
{
	cosmology_t *thisCosmology;
	
	thisCosmology = malloc(sizeof(cosmology_t));
	
	thisCosmology->h = 0.;
	thisCosmology->omega_m = 0.;
	thisCosmology->omega_l = 0.;
	thisCosmology->omega_b = 0.;
	
	thisCosmology->z = 0.;
	thisCosmology->Hubble_z = 0.;
	
	thisCosmology->Y = 0.;
	thisCosmology->nH_z = 0.;
	thisCosmology->nHe_z = 0.;
    
    thisCosmology->f3He = 0.;
	
	return thisCosmology;
}

void deallocate_cosmology(cosmology_t *thisCosmology)
{
	free(thisCosmology);
}

cosmology_t *allocate_cosmology(double h, double omega_m, double omega_l, double omega_b, double z, double Y)
{
	cosmology_t *thisCosmology;
	
	thisCosmology = initCosmology();
	
	thisCosmology->h = h;
	thisCosmology->omega_m = omega_m;
	thisCosmology->omega_l = omega_l;
	thisCosmology->omega_b = omega_b;
	
	thisCosmology->z = z;
	thisCosmology->Hubble_z = calc_Hubble_z(h, omega_m, omega_l, z);
	
	thisCosmology->Y = Y;
	thisCosmology->nH_z = calc_nH_z(h, omega_b, z, Y);
	thisCosmology->nHe_z = calc_nHe_z(h, omega_b, z, Y);
    
    thisCosmology->f3He = 1.e-5;
	
	return thisCosmology;
}

void cosmology_update_z(cosmology_t *thisCosmology, double z)
{
	double h = thisCosmology->h;
	double omega_m = thisCosmology->omega_m;
	double omega_l = thisCosmology->omega_l;
	double omega_b = thisCosmology->omega_b;
	double Y = thisCosmology->Y;
	
	thisCosmology->z = z;
	thisCosmology->Hubble_z = calc_Hubble_z(h, omega_m, omega_l, z);
	thisCosmology->nH_z = calc_nH_z(h, omega_b, z, Y);
	thisCosmology->nHe_z = calc_nHe_z(h, omega_b, z, Y);
}

/*-------------------------------------------------------------------------------------*/
/* QUANTITIES AT REDSHIFT Z */
/*-------------------------------------------------------------------------------------*/

double calc_Hubble_z(double h, double omega_m, double omega_l, double z)
{
	return h*100.*km_cm/Mpc_cm*sqrt(omega_m*CUB(1.+z)+omega_l);
}

double calc_nH_z(double h, double omega_b, double z, double Y)		//number density per proton / H-atom
{
	return 3.*SQR(h*100.*km_cm/Mpc_cm)/(8.*M_PI*G)/mp_g*omega_b*CUB(1.+z)*(1.-Y);
}

double calc_nHe_z(double h, double omega_b, double z, double Y)		//number density per He-atom
{
	return 3.*SQR(h*100.*km_cm/Mpc_cm)/(8.*M_PI*G)/mp_g*omega_b*CUB(1.+z)*Y*0.25;
}
