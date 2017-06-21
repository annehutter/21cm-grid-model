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
#include "spin_temperature.h"
#include "3cm_grid.h"
#include "spin_temperature_He.h"

double coupling_coll_He(double z, double ne, double Tk)
{
    double Tbg = T_CMB(z);
    
    double tmp = 14.3*ev_to_erg *bohrHe_cm*bohrHe_cm / sqrt(boltzman_cgs*Tk*M_PI*me_g/8.);
    
    double result = tmp * ne*T3cm/(A10_He*Tbg);

    return result;
}

double calc_Tk_He(xray_spectrum_t *thisXraySpectrum)
{
    double const alphaX = thisXraySpectrum->alphaX;
    double const mean_energy = planck_cgs*nu_HI*alphaX/(alphaX+1.);
    
    return 5.5e4*pow(mean_energy/(40.8*ev_to_erg), 0.25);
}

void compute_Ts_He_on_grid(xray_spectrum_t *thisXraySpectrum, grid_3cm_t *this3cmGrid, grid_21cm_t *this21cmGrid, cosmology_t *thisCosmology)
{
    double z = thisCosmology->z;
    
    double Tbg_inv = 1./T_CMB(z);
    double Tk, Tk_inv;
    double Ta_inv;
    double Ts_inv;
    
    double dens, nH, nHe, Xe, ne;
    double xa, xc;
    
    int nbins;
    int local_n0;
    
    nbins = this21cmGrid->nbins;
    local_n0 = this21cmGrid->local_n0;

    xa = 0.;
    Tk = calc_Tk_He(thisXraySpectrum);
    Tk_inv = 1./Tk;
    Ta_inv = Tk_inv;   
    
    debug_printf(DEBUG_TS_HE_CALCULATION, "+DEBUG+ nHe = %e cm^-3\t nH = %e cm^-3\t xc = %e\n",  calc_nHe_z(thisCosmology->h, thisCosmology->omega_b, z, thisCosmology->Y), calc_nH_z(thisCosmology->h, thisCosmology->omega_b, z, thisCosmology->Y), coupling_coll_He(z, calc_nHe_z(thisCosmology->h, thisCosmology->omega_b, z, thisCosmology->Y)+calc_nH_z(thisCosmology->h, thisCosmology->omega_b, z, thisCosmology->Y), 1.e4));
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                dens = creal(this21cmGrid->dens[i*nbins*nbins+j*nbins+k]);
                nHe = thisCosmology->nHe_z*dens;
                nH = thisCosmology->nH_z*dens;
                Xe = creal(this21cmGrid->Xe[i*nbins*nbins+j*nbins+k]);
                ne = (nHe+nH)*Xe;
                xc = coupling_coll_He(z, ne, Tk);
                Ts_inv = calc_Ts_inv(xa, xc, Ta_inv, Tk_inv, Tbg_inv);
                this3cmGrid->Ts_inv[i*nbins*nbins+j*nbins+k] = Ts_inv + 0.*I;
            }
        }
    }
}
