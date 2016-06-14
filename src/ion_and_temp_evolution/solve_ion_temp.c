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

#include "../xray_heating_lya_coupling/phys_const.h"
#include "../xray_heating_lya_coupling/chem_const.h"
#include "../xray_heating_lya_coupling/cosmology.h"
#include "../xray_heating_lya_coupling/collisional_coupling.h"
#include "../xray_heating_lya_coupling/xray_heating.h"
#include "../xray_heating_lya_coupling/wouthuysen_effect.h"
#include "../xray_heating_lya_coupling/21cm_grid.h"
#include "../xray_heating_lya_coupling/spin_temperature.h"
#include "../xray_heating_lya_coupling/compute_grid.h"

#include "recomb_rates.h"
#include "solve_ion_temp.h"

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

double update_ion(grid_21cm_t *this21cmGrid, recomb_t *thisRecombRates, int index, double z, double z_old, double gamma, double n, double Hubble)
{
	double temp = this21cmGrid->temp[index];
	
	double dens = n*this21cmGrid->dens[index];
	
	double Xe_old = this21cmGrid->Xe[index];
	
	double recomb = thisRecombRates->recHII[temp_recomb_index(thisRecombRates, temp)];
	
	double z_old_CUB = sqrt(CUB(1.+z_old));
	
	double exponent = sqrt(4.*dens*recomb*gamma/SQR(Hubble) + 2.25);
	double a = pow(1.+z, exponent);
	double b = pow(1.+z_old, exponent);
	
	double a_plus_b = a+b;
	double b_minus_a = b-a;
	
// 	printf("Xe_old = %e\t Xe = %e\tgamma = %e\n", Xe_old, Xe_old/sqrt(CUB(1.+z))*z_old_CUB * ( -1.5*b_minus_a + a_plus_b*exponent + 2.*gamma*b_minus_a/(Hubble*Xe_old*z_old_CUB) ) / ( 1.5*b_minus_a + a_plus_b*exponent + 2.*dens*recomb*b_minus_a*Xe_old*z_old_CUB/Hubble ), gamma);
	return Xe_old/sqrt(CUB(1.+z))*z_old_CUB * ( -1.5*b_minus_a + a_plus_b*exponent + 2.*gamma*b_minus_a/(Hubble*Xe_old*z_old_CUB) ) / ( 1.5*b_minus_a + a_plus_b*exponent + 2.*dens*recomb*b_minus_a*Xe_old*z_old_CUB/Hubble );
}

double update_temp(grid_21cm_t *this21cmGrid, int index, double sqrt_cub_z, double sqrt_cub_z_old, double factor, double xray_heating, double Hubble)
{
	double temp_old = this21cmGrid->temp_old[index];
	
	double dens = this21cmGrid->dens[index];
	double dens_old = this21cmGrid->dens_old[index];
	
	double Xe = this21cmGrid->Xe[index];
	double Xe_old = this21cmGrid->Xe_old[index];
	
// 	printf("temp_old = %e\t temp = %e\t heating = %e\t %e\t%e\n", temp_old, temp_old*SQR(1.+dens)*(1.+Xe)/(SQR(1.+dens_old)*(1.+Xe_old)) + 4.*xray_heating/(21.*boltzman_cgs*Hubble)*(sqrt_cub_z-sqrt_cub_z_old*factor), xray_heating, temp_old*SQR(1.+dens)*(1.+Xe)/(SQR(1.+dens_old)*(1.+Xe_old)), 4.*xray_heating/(21.*boltzman_cgs*Hubble)*(sqrt_cub_z-sqrt_cub_z_old*factor));
	
	return temp_old*SQR(1.+dens)*(1.+Xe)/(SQR(1.+dens_old)*(1.+Xe_old)) + 4.*xray_heating/(21.*boltzman_cgs*Hubble)*(sqrt_cub_z-sqrt_cub_z_old*factor);
}

double update_dens(grid_21cm_t *this21cmGrid, int index, double D)
{
	double dens_old = this21cmGrid->dens_old[index];
	
	return dens_old*D;
}

void overwrite_old_values_21cmgrid(grid_21cm_t *this21cmGrid)
{
  	int nbins = this21cmGrid->nbins;
	int local_n0 = this21cmGrid->local_n0;
	int index;
	
	this21cmGrid->z_old = this21cmGrid->z;
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				index = i*nbins*nbins+j*nbins+k;
				this21cmGrid->Xe_old[index] = this21cmGrid->Xe[index];
				this21cmGrid->dens_old[index] = this21cmGrid->dens[index];
				this21cmGrid->temp_old[index] = this21cmGrid->temp[index];
			}
		}
	}
}

void compute_temp_ion_grid(grid_21cm_t *this21cmGrid, recomb_t *thisRecombRates, cosmology_t *thisCosmology, xray_grid_t *thisXray_grid)
{
	int nbins = this21cmGrid->nbins;
	int local_n0 = this21cmGrid->local_n0;
	int index;
	
	double z = this21cmGrid->z;
	double z_old = this21cmGrid->z_old;
	
	double h = thisCosmology->h;
	double omega_m = thisCosmology->omega_m;
	double omega_b = thisCosmology->omega_b;
	double Hubble = h*100.*km_cm/Mpc_cm*sqrt(omega_m);
	double n = 3.*SQR(h*100.*km_cm/Mpc_cm)/(8.*M_PI*G)/mp_g*omega_b;

	double factor = SQR((1.+z)/(1.+z_old));
	double sqrt_cub_z = pow(1.+z,-1.5);
	double sqrt_cub_z_old = pow(1.+z_old,-1.5);
	
	printf("factor = %e\t sqrt_cub_z = %e\t sqrt_cub_z_old = %e\n", factor, sqrt_cub_z, sqrt_cub_z_old);
	
	double gamma;
	double xray_heating;
	double D = 1.;
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				gamma = (thisXray_grid->xray_ionization_HI[i*nbins*nbins+j*nbins+k] + thisXray_grid->xray_ionization_HeI[i*nbins*nbins+j*nbins+k])*(1.-this21cmGrid->Xe[i*nbins*nbins+j*nbins+k]);
				
				xray_heating = (thisXray_grid->xray_heating_HI[i*nbins*nbins+j*nbins+k] + thisXray_grid->xray_heating_HeI[i*nbins*nbins+j*nbins+k])*(1.-this21cmGrid->Xe[i*nbins*nbins+j*nbins+k]);
				
// 				printf("gamma = %e\t heating = %e\n", gamma, xray_heating);
				
				index = i*nbins*nbins+j*nbins+k;
				
				this21cmGrid->Xe[i*nbins*nbins+j*nbins+k] = update_ion(this21cmGrid, thisRecombRates, index, z, z_old, gamma, n, Hubble);
				this21cmGrid->temp[i*nbins*nbins+j*nbins+k] = update_temp(this21cmGrid, index, sqrt_cub_z, sqrt_cub_z_old, factor, xray_heating, Hubble);
				
				this21cmGrid->dens[i*nbins*nbins+j*nbins+k] = update_dens(this21cmGrid, index, D);
			}
		}
	}
}