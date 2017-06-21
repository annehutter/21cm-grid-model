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

#include "../xray_heating_lya_coupling/phys_const.h"
#include "../xray_heating_lya_coupling/chem_const.h"
#include "../xray_heating_lya_coupling/cosmology.h"
#include "../xray_heating_lya_coupling/collisional_coupling.h"
#include "../xray_heating_lya_coupling/xray_heating.h"
#include "../xray_heating_lya_coupling/wouthuysen_effect.h"
#include "../xray_heating_lya_coupling/21cm_grid.h"
#include "../xray_heating_lya_coupling/spin_temperature.h"
#include "../xray_heating_lya_coupling/3cm_grid.h"
#include "../xray_heating_lya_coupling/spin_temperature_He.h"
#include "../xray_heating_lya_coupling/compute_grid.h"

#include "recomb_rates.h"
#include "fractions_heat_ion.h"
#include "solve_ion_temp.h"

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

void update_ion_temp(recomb_t *thisRecombRates, double dens_old, double dens, double *Xe_old, double *Xe, double *temp_old, double *temp, double z_old, double z, double gamma, double n, double xray_heating, double Hubble)
{
    double recomb = thisRecombRates->recHII[temp_recomb_index(thisRecombRates, *temp_old)];
    
    *Xe = update_ion_smallXe(n*dens, *Xe_old, recomb, z, z_old, gamma, Hubble);
    *temp = update_temp(dens_old, dens, *Xe_old, *Xe, *temp_old, z, z_old, xray_heating, Hubble);
    
//     printf("Xe_old = %e\t Xe = %e\t z_old = %e\t z = %e\n", *Xe_old, *Xe, z_old, 0.5*(z_old+z));
//     if((*Xe)-(*Xe_old) > 0.1)
//     {
//         update_ion_temp(thisRecombRates, dens_old, dens, Xe_old, Xe, temp_old, temp, z_old, 0.5*(z_old+z), gamma, n, xray_heating, Hubble);
//         update_ion_temp(thisRecombRates, dens_old, dens, Xe_old, Xe, temp_old, temp, 0.5*(z_old+z), z, gamma, n, xray_heating, Hubble);
//     }
    
    if((*Xe) > 0.01)
    {
        *Xe = update_ion_constdens(n*dens*CUB(1.+z_old), *Xe_old, recomb, z, z_old, gamma, Hubble);
        *temp = update_temp(dens_old, dens, *Xe_old, *Xe, *temp_old, z, z_old, xray_heating, Hubble);
    }
    
    Xe_old = Xe;
    temp_old = temp;
}

double update_ion_smallXe(double dens, double Xe_old, double recomb, double z, double z_old, double gamma, double Hubble)
{    
    double z_old_CUB = sqrt(CUB(1.+z_old));
    
    double exponent = sqrt(4.*dens*recomb*gamma/SQR(Hubble) + 2.25);
    
    double a_div_b = pow((1.+z)/(1.+z_old), exponent);
    
    
    double factor = exponent*2.*Hubble;
    
    double tmp = (1.+ a_div_b)*factor;
    double tmp2 = 1.-a_div_b;
    
    double numerator = tmp + tmp2 * (-3.*Hubble + 4.*gamma/(Xe_old*z_old_CUB));
    double denominator = tmp + tmp2 * (3.*Hubble + 4.*dens*recomb*Xe_old*z_old_CUB);
    
    double Xe_new = Xe_old/sqrt(CUB(1.+z))*z_old_CUB * numerator/denominator;
    
    return Xe_new;
}

double update_ion_constdens(double dens, double Xe_old, double recomb, double z, double z_old, double gamma, double Hubble)
{
    double z_old_CUB = sqrt(CUB(1.+z_old));
    double z_CUB = sqrt(CUB(1.+z));
    
    double tmp = gamma/(2.*dens*recomb);
    double tmp2 = sqrt(1.+2./tmp);
    
    double atanh_argument = Xe_old+tmp/(tmp*tmp2);
    double tanh_argument = gamma*tmp2*(1./z_old_CUB - 1./z_CUB)/(3.*Hubble) - atanh(atanh_argument);
    
    double Xe_new = - tmp * (1. + tmp2 * tanh(tanh_argument));
    
    return Xe_new;
}

double update_temp(double dens_old, double dens, double Xe_old, double Xe, double temp_old, double z, double z_old, double xray_heating, double Hubble)
{
    double factor = SQR((1.+z)/(1.+z_old));
    double sqrt_cub_z = pow(1.+z,-1.5);
    double sqrt_cub_z_old = pow(1.+z_old,-1.5);
    
    double temp_new = temp_old*factor*SQR(1.+dens)*(1.+Xe)/(SQR(1.+dens_old)*(1.+Xe_old)) + 4.*xray_heating/(21.*boltzman_cgs*Hubble)*(sqrt_cub_z-sqrt_cub_z_old*factor);
    
    return temp_new;
}


/* =====================================================================================================*/

double update_dens(grid_21cm_t *this21cmGrid, int index, double D)
{
    double dens_old = this21cmGrid->dens_old[index];
    
    return dens_old*D;
}


/* =====================================================================================================*/

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

//     double factor = SQR((1.+z)/(1.+z_old));
//     double sqrt_cub_z = pow(1.+z,-1.5);
//     double sqrt_cub_z_old = pow(1.+z_old,-1.5);
        
    double dens, dens_old;
    double Xe, Xe_old;
    double temp, temp_old;
    
    double gamma=0.;
//     double gamma_heat=0.;
    double xray_heating=0.;
    double D = 1.;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                index = i*nbins*nbins+j*nbins+k;

                dens_old = creal(this21cmGrid->dens_old[index]);
                dens = creal(this21cmGrid->dens[index]);
                
                Xe_old = creal(this21cmGrid->Xe_old[index]);
                Xe = creal(this21cmGrid->Xe[index]);
                
                temp_old = creal(this21cmGrid->temp_old[index]);
                temp = creal(this21cmGrid->temp[index]);
                
                if(isnan(Xe) || isinf(Xe) || Xe <= 0. || Xe > 1.)
                {
                    printf("Xe = %e\n", Xe);
                }
                
                gamma = (creal(thisXray_grid->xray_ionization_HI[index]) + creal(thisXray_grid->xray_ionization_HeI[index]))*(1.-Xe);
//                 gamma_heat = (creal(thisXray_grid->xray_heating_HI[index]) + creal(thisXray_grid->xray_heating_HeI[index]))*(1.-Xe)*(fion_HI(Xe)/(planck_cgs*nu_HI) + fion_HeI(Xe)/(planck_cgs*nu_HeI));

                xray_heating = (creal(thisXray_grid->xray_heating_HI[index]) + creal(thisXray_grid->xray_heating_HeI[index]))*(1.-Xe)*fheat(Xe);
                                
                if(isnan(gamma))
                {
                    printf("index %d: gamma = %e\t heating = %e\t xray_ion_HI = %e\txray_ion_HeI = %e\t Xe = %e\n", index, gamma, xray_heating, creal(thisXray_grid->xray_ionization_HI[index]), creal(thisXray_grid->xray_ionization_HeI[index]), Xe);
                }
                
                this21cmGrid->dens[index] = update_dens(this21cmGrid, index, D) + I*0.;
                
                update_ion_temp(thisRecombRates, dens_old, dens, &Xe_old, &Xe, &temp_old, &temp, z_old, z, gamma, n, xray_heating, Hubble);

                this21cmGrid->Xe[index] = Xe + I*0.;
                this21cmGrid->temp[index] = temp + I*0.;
            }
        }
    }
}
