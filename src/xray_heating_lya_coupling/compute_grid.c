#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>    //included because M_PI is not defined in <math.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "../settings.h"

#include "phys_const.h"
#include "chem_const.h"
#include "cosmology.h"
#include "collisional_coupling.h"
#include "xray_heating.h"
#include "wouthuysen_effect.h"
#include "21cm_grid.h"
#include "3cm_grid.h"
#include "spin_temperature.h"
#include "spin_temperature_He.h"
#include "compute_grid.h"

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

/* compute 21cm brightness temperature */
void compute_Tb_grid(grid_21cm_t *this21cmGrid, cosmology_t *thisCosmology)
{
    int nbins = this21cmGrid->nbins;
    int local_n0 = this21cmGrid->local_n0;
    
    double factor;
    double Tbg;
    
    double Hubble_z = thisCosmology->Hubble_z;
    double nH = thisCosmology->nH_z;
    double z = thisCosmology->z;
    
    nbins = this21cmGrid->nbins;
    local_n0 = this21cmGrid->local_n0;
    
    factor = (3.*clight_cm*SQR(lambda_21cm)*planck_cgs*A10*nH)/(32.*M_PI*boltzman_cgs*(1.+z)*Hubble_z);
    Tbg = T_CMB(z);
        
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                this21cmGrid->Tb[i*nbins*nbins+j*nbins+k] = factor*creal(this21cmGrid->XHI[i*nbins*nbins+j*nbins+k])*creal(this21cmGrid->dens[i*nbins*nbins+j*nbins+k])*(1. - Tbg * creal(this21cmGrid->Ts_inv[i*nbins*nbins+j*nbins+k])) + 0.*I;
                debug_printf(DEBUG_21CM_TB_CALCULATION, "+DEBUG+ 3 c lambda21^2 h A10 nH / (32 pi kB (1+z) H(z)) = %e K\t TCMB = %e K\t Tb = %e K\t Ts = %e K\n", factor, Tbg, creal(this21cmGrid->Tb[i*nbins*nbins+j*nbins+k]), 1./creal(this21cmGrid->Ts_inv[i*nbins*nbins+j*nbins+k]));
            }
        }
    }
}

/* compute xray and lya fields, get ionization and temperature & compute spin temperature */
void do_step_21cm_emission(cosmology_t *thisCosmology, xray_grid_t *thisXray_grid, xray_spectrum_t *thisXray_spectrum, lya_grid_t * thisLya_grid, lya_spectrum_t *thisLya_spectrum, grid_21cm_t *this21cmGrid, k10_t *k10_table, double Xe, int myRank)
{
    /* compute xray heating and ionization */
    xray_heating_and_ionization(thisXray_grid, thisCosmology, Xe, thisXray_spectrum, myRank);
      if(myRank == 0) printf("  21cm emission (do_step_21cm_emission): xray heating and ionization done\n");

    /* compute wouthuysen coupling */
    lya_wouthuysen_coupling(thisLya_grid, thisLya_spectrum, thisXray_grid, this21cmGrid, thisCosmology);
    if(myRank == 0) printf("  21cm emission (do_step_21cm_emission): lya coupling done\n");

    /* compute the spin temperature */
    compute_Ts_on_grid(thisLya_grid, k10_table, this21cmGrid, thisCosmology);
    if(myRank == 0) printf("  21cm emission (do_step_21cm_emission): Ts calclated\n");

    /* compute 21cm emission/absorption */
    compute_Tb_grid(this21cmGrid, thisCosmology);
    if(myRank == 0) printf("  21cm emission (do_step_21cm_emission): Tb calculated\n");
}


/* compute 3cm brightness temperature */
void compute_Tb_He_grid(grid_3cm_t *this3cmGrid, grid_21cm_t *this21cmGrid, cosmology_t *thisCosmology)
{
    int nbins = this21cmGrid->nbins;
    int local_n0 = this21cmGrid->local_n0;
    
    double factor;
    double Tbg;
    
    double Hubble_z = thisCosmology->Hubble_z;
    double nHe = thisCosmology->nHe_z;
    double f3He = thisCosmology->f3He;
    double z = thisCosmology->z;
    
    nbins = this21cmGrid->nbins;
    local_n0 = this21cmGrid->local_n0;
    
    z = thisCosmology->z;
    factor = (clight_cm*SQR(lambda_3cm)*planck_cgs*A10_He*nHe*f3He)/(32.*M_PI*boltzman_cgs*(1.+z)*Hubble_z);
    Tbg = T_CMB(z);
        
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                this3cmGrid->Tb[i*nbins*nbins+j*nbins+k] = factor*creal(this21cmGrid->dens[i*nbins*nbins+j*nbins+k])*(1. - Tbg * creal(this3cmGrid->Ts_inv[i*nbins*nbins+j*nbins+k])) + 0.*I;
                debug_printf(DEBUG_3CM_TB_CALCULATION, "+DEBUG+ c lambdas^2 h A10He nHe f3He / (32 pi kB (1+z) H(z)) = %e K\t TCMB = %e K\t Tb = %e K\t Ts = %e K\n", factor, Tbg, creal(this3cmGrid->Tb[i*nbins*nbins+j*nbins+k]), 1./creal(this3cmGrid->Ts_inv[i*nbins*nbins+j*nbins+k]));
            }
        }
    }
}

void do_step_3cm_emission(cosmology_t *thisCosmology, xray_spectrum_t *thisXraySpectrum, grid_21cm_t *this21cmGrid, grid_3cm_t *this3cmGrid, int myRank)
{
    /* compute the spin temperature */
    compute_Ts_He_on_grid(thisXraySpectrum, this3cmGrid, this21cmGrid, thisCosmology);
    if(myRank == 0) printf("  3cm emission (do_step_3cm_emission): Ts calclated\n");

    /* compute 3cm emission/absorption */
    compute_Tb_He_grid(this3cmGrid, this21cmGrid, thisCosmology);
    if(myRank == 0) printf("  3cm emission (do_step_3cm_emission): Tb calculated\n");
}
