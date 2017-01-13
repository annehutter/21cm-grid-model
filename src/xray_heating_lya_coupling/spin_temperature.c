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

#include "../ion_and_temp_evolution/fractions_heat_ion.h"

#include "phys_const.h"
#include "chem_const.h"
#include "cosmology.h"
#include "collisional_coupling.h"
#include "xray_heating.h"
#include "wouthuysen_effect.h"
#include "21cm_grid.h"
#include "spin_temperature.h"

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

double tau_gp(double Hubble_z_inv, double nHI)
{
	return nHI*Hubble_z_inv*SQR(e)*falpha*lambda_lya/(3.*me_g*1.e-3*clight_cm*epsilon0);
}

double T_CMB(double z)
{
	return 2.73*(1.+z);
}

double coupling_coll(double z, double nH, k10_t* k10_table, double Tk)
{
	double Tbg = T_CMB(z);
	
	return (get_k10_HH(k10_table, Tk) + get_k10_eH(k10_table, Tk) + get_k10_pH(k10_table, Tk))*nH*T21cm/(A10*Tbg);
}

double coupling_alpha(double prefac_z, double modSalpha, double Jalpha)
{
	return prefac_z*modSalpha*Jalpha;
}

double coupling_alpha_prefac_z(double z, double prefac)
{
	return prefac/T_CMB(z);
}

double calc_modSalpha(double Hubble_z_inv, double nHI, double Tk_inv, double Ts_inv)
{
	return (1. - 0.0631789*Tk_inv + 0.115995*SQR(Tk_inv) - 0.401403*Ts_inv*Tk_inv + 0.336463*Ts_inv*SQR(Tk_inv))/(1. + 2.98394*calc_xi(Hubble_z_inv, nHI, Tk_inv) + 1.53583*SQR(calc_xi(Hubble_z_inv, nHI, Tk_inv)) + 3.85289*CUB(calc_xi(Hubble_z_inv, nHI, Tk_inv)));
}

double calc_xi(double Hubble_z_inv, double nHI, double Tk_inv)
{
	return pow(1.e-7*tau_gp(Hubble_z_inv, nHI)*SQR(Tk_inv), 1./3.);
}

double calc_Teff_inv(double Tk_inv, double Ts_inv)
{
	return Tk_inv + 0.405535*Tk_inv * (Ts_inv - Tk_inv);
}

double calc_Ts_inv(double xa_mod, double xc, double Teff_inv, double Tk_inv, double Tbg_inv)
{
	return (Tbg_inv + xa_mod*Teff_inv + xc*Tk_inv)/(1. + xa_mod + xc);
}

double compute_Ts_inv(double Hubble_z_inv, double nH, double XHI, double coupling_alpha_prefac_z, double Jalpha, double xc, double Tk_inv, double Tbg_inv)
{
	double nHI = nH*XHI;
	double Ts_inv = Tbg_inv;
	double Teff_inv;
	double xa_mod;
	
	int N = 5;
	
	for(int i=0; i<N; i++)
	{
		xa_mod = coupling_alpha(coupling_alpha_prefac_z, calc_modSalpha(Hubble_z_inv, nHI, Tk_inv, Ts_inv), Jalpha);
		Teff_inv = calc_Teff_inv(Tk_inv, Ts_inv);
		Ts_inv = calc_Ts_inv(xa_mod, xc, Teff_inv, Tk_inv, Tbg_inv);
// 		printf("%d: \txa_mod = %e\t Teff = %e\t xc = %e\t Ts = %e\t Tk = %e\n", i, xa_mod, 1./Teff_inv, xc,  1./Ts_inv, 1./Tk_inv);
	}
// 	printf("xa_mod = %e\n", xa_mod);
	
	return Ts_inv;
}

void compute_Ts_on_grid(lya_grid_t *thisLya_grid, k10_t *k10_table, grid_21cm_t *this21cmGrid, cosmology_t *thisCosmology)
{
	double Hubble_z_inv = 1./thisCosmology->Hubble_z;
	double z = thisCosmology->z;
	
	double coupling_alpha_prefac = (4.*M_PI*T21cm*SQR(e)*falpha)/(27.*epsilon0*A10*me_g*1.e-3*clight_cm);
	double xalpha_prefac_z = coupling_alpha_prefac_z(z, coupling_alpha_prefac);

	double Tbg_inv = 1./T_CMB(z);
	double Tk, Tk_inv;
	double Ts_inv;
	
	double dens, nH, XHI;
	double Jalpha;
	double xc;
	
	int nbins;
	int local_n0;
	
	nbins = this21cmGrid->nbins;
	local_n0 = this21cmGrid->nbins;
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				dens = creal(this21cmGrid->dens[i*nbins*nbins+j*nbins+k]);
				nH = thisCosmology->nH_z*dens;
				XHI = 1.-creal(this21cmGrid->Xe[i*nbins*nbins+j*nbins+k]);
				Jalpha = creal(thisLya_grid->lya[i*nbins*nbins+j*nbins+k]);
				Tk = creal(this21cmGrid->temp[i*nbins*nbins+j*nbins+k]);
				xc = coupling_coll(z, nH, k10_table, Tk);
				Tk_inv = 1./Tk;
				Ts_inv = compute_Ts_inv(Hubble_z_inv, nH, XHI, xalpha_prefac_z, Jalpha, xc, Tk_inv, Tbg_inv);
				this21cmGrid->Ts_inv[i*nbins*nbins+j*nbins+k] = Ts_inv + 0.*I;
// 				printf("dens = %e\t nH = %e\t XHI = %e\t Jalpha = %e\t xc = %e\t Tk = %e\n", dens, nH, XHI, Jalpha, xc, Tk);
// 				printf("Ts = %e\n", 1./Ts_inv);
			}
		}
	}
}

/*--------------------------------------------------------------------------------------------------*/
/* Wouthuysen coupling */

void lya_wouthuysen_coupling(lya_grid_t *thisLya_grid, lya_spectrum_t *thisSpectrum, xray_grid_t *thisXray_grid, grid_21cm_t *this21cmGrid, cosmology_t *thisCosmology)
{
    int nbins = thisLya_grid->nbins;
	int local_n0 = thisLya_grid->local_n0;
	
	/* compute comoving Lya flux */
	lya_bg_source_emission(thisLya_grid, thisSpectrum, thisCosmology);
	lya_bg_lya_excitation(thisLya_grid, thisXray_grid, this21cmGrid, thisCosmology);
}

void lya_bg_lya_excitation(lya_grid_t *thisLya_grid, xray_grid_t *thisXray_grid, grid_21cm_t *this21cmGrid, cosmology_t *thisCosmology)
{
	int nbins = thisLya_grid->nbins;
	int local_n0 = thisLya_grid->local_n0;
	
	double factor = SQR(lambda_lya)/(4.*M_PI*planck_cgs*clight_cm*thisCosmology->Hubble_z);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
                double const Xe = this21cmGrid->Xe[i*nbins*nbins+j*nbins+k];
                
				double const heating = creal(thisXray_grid->xray_heating_HI[i*nbins*nbins+j*nbins+k]) + creal(thisXray_grid->xray_heating_HeI[i*nbins*nbins+j*nbins+k]) + creal(thisXray_grid->xray_heating_HeII[i*nbins*nbins+j*nbins+k]);
				
				if(heating<0.) 
				{
					fprintf(stderr, "heating is negative!!! and has value = %e\nCheckyour x-ray heating!\nStopping execution.\n", heating);
					exit(EXIT_FAILURE);
				}
				thisLya_grid->lya[i*nbins*nbins+j*nbins+k] = creal(thisLya_grid->lya[i*nbins*nbins+j*nbins+k]) + heating*flya(Xe)*factor + 0.*I;
            }
		}
	}
}
