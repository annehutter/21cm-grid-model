#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>	//included because M_PI is not defined in <math.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "xray_heating_lya_coupling/phys_const.h"
#include "xray_heating_lya_coupling/chem_const.h"
#include "xray_heating_lya_coupling/cosmology.h"
#include "xray_heating_lya_coupling/collisional_coupling.h"
#include "xray_heating_lya_coupling/xray_heating.h"
#include "xray_heating_lya_coupling/wouthuysen_effect.h"
#include "xray_heating_lya_coupling/21cm_grid.h"
#include "xray_heating_lya_coupling/spin_temperature.h"
#include "xray_heating_lya_coupling/compute_grid.h"

#include "ion_and_temp_evolution/cell.h"
#include "ion_and_temp_evolution/recomb_rates.h"
#include "ion_and_temp_evolution/evolution_loop.h"
#include "ion_and_temp_evolution/solve_ionfraction.h"
#include "ion_and_temp_evolution/solve_temperature.h"
#include "ion_and_temp_evolution/evolution_loop.h"

#include "evolution.h"

void update_21cmgrid(grid_21cm_t *this21cmGrid, cell_t *thisCell, int x, int y, int z)
{
	int nbins = this21cmGrid->nbins;
	
	this21cmGrid->temp[x*nbins*nbins+y*nbins+z] = thisCell->temp;
	this21cmGrid->XHI[x*nbins*nbins+y*nbins+z] = thisCell->XHI;
	this21cmGrid->XHeI[x*nbins*nbins+y*nbins+z] = thisCell->XHeI;
	this21cmGrid->XHeII[x*nbins*nbins+y*nbins+z] = thisCell->XHeII;
	
	this21cmGrid->Xe[x*nbins*nbins+y*nbins+z] = 1.-thisCell->XHI;
}

void evolve()
{
/*-------------------------------------------------------------------------------------*/
/* COSMOLOGY */
/*-------------------------------------------------------------------------------------*/
	double h = 0.7;
	double omega_m = 0.27;
	double omega_b = 0.045;
	double omega_l = 0.73;
	
	double Y = 0.24;
	
/*-------------------------------------------------------------------------------------*/
/* XRAY Spectrum */
/*-------------------------------------------------------------------------------------*/
	double alphaX = -1.5;
	double nuX_min = 0.1e3*ev_to_erg/planck_cgs;
	double lumX = 3.4e40/(1.e3*ev_to_erg);
	printf("alphaX = %e\t nuX_min = %e\t lumX = %e\n", alphaX, nuX_min, lumX);
/*-------------------------------------------------------------------------------------*/
/* LYA Spectrum */
/*-------------------------------------------------------------------------------------*/
	double alphaLya = -1.5;
	double nuLya_min = 10.3*ev_to_erg/planck_cgs;
	double lumLya = 3.4e40/(1.e3*ev_to_erg);

/*-------------------------------------------------------------------------------------*/
/* GRID */
/*-------------------------------------------------------------------------------------*/
	int nbins = 128;
	int local_n0 = nbins;
	double box_size = 80.;
	
	char densfile[128];
	char xray_sourcefile[128];
	char lya_sourcefile[128];
	
	for (int i=0; i<128; i++)
	{
		densfile[i] = '\0';
		xray_sourcefile[i] = '\0';
		lya_sourcefile[i] = '\0';
	}
	strcat(densfile, "/home/anne/PostdocSwinburne/grid_model_mpi/inputfiles/grid128/density11_ic.in");
	strcat(xray_sourcefile, "/home/anne/PostdocSwinburne/grid_model_mpi/inputfiles/grid128/density11_ic.in");
	strcat(lya_sourcefile, "/home/anne/PostdocSwinburne/grid_model_mpi/inputfiles/grid128/density11_ic.in");
	
/*-------------------------------------------------------------------------------------*/
/* Initialize for ionization and temperature evolution */
/*-------------------------------------------------------------------------------------*/
	double zstart = 25.;
	double zend = 10.;
	double dz = 0.5;
  
	cell_t *cell;
	cell = initCell_temp(14., 1.);
	
	recomb_t * recomb_rates;
	recomb_rates = calcRecRate();
	
	int dim_matrix = 2;
	
/*-------------------------------------------------------------------------------------*/
/* Initialize for 21cm computation */
/*-------------------------------------------------------------------------------------*/

	cosmology_t *thisCosmology;
	
	xray_grid_t *thisXray_grid;
	xray_spectrum_t *thisXray_spectrum;
	
	lya_grid_t *thisLya_grid;
	lya_spectrum_t *thisLya_spectrum;
	
	grid_21cm_t *this21cmGrid;
	
	k10_t *k10_table;
	
	int double_precision = 0;
	
	thisCosmology = allocate_cosmology(h, omega_m, omega_l, omega_b, zstart, Y);
	
	thisXray_grid = allocate_xray_grid(nbins, box_size);
	thisXray_spectrum = allocate_xray_spectrum(lumX, alphaX, nuX_min);
	
	thisLya_grid = allocate_lya_grid(nbins, box_size);
	thisLya_spectrum = allocate_lya_spectrum(lumLya, alphaLya, nuLya_min);
	
	this21cmGrid = allocate_21cmgrid(nbins, box_size);
	
	read_density_21cmgrid(this21cmGrid, densfile, double_precision);
	read_lum_xraygrid(thisXray_grid, xray_sourcefile, double_precision);
	read_lum_lyagrid(thisLya_grid, lya_sourcefile, double_precision);
	
	set_temperature_21cmgrid(this21cmGrid, 14.);
	
	k10_table = create_k10_table();
	
/*-------------------------------------------------------------------------------------*/
/* REDSHIFT EVOLUTION */
/*-------------------------------------------------------------------------------------*/

	const double tmp_HI = 1./(boltzman_cgs*nu_HI);
	const double tmp_HeI = 1./(boltzman_cgs*nu_HeI);
	const double tmp_HeII = 1./(boltzman_cgs*nu_HeII);
	
	double fion_HI = 1.;
	double fion_HeI = 1.;
	double fion_HeII = 1.;
	
	double fheat = 1.;
	
	double Xe = 0.;
	
	double dz_tmp = dz;
	
	for(double z = zstart; z>zend; z = z-dz_tmp)
	{
		printf("z = %e\t %e\n", z, fmod(z,dz));
		cosmology_update_z(thisCosmology, z);
		if(fmod(z,dz) == 0.)
		{
			Xe = get_mean_Xe_21cmgrid(this21cmGrid);
			printf("Xe = %e\n", Xe);
			do_step_21cm_emission(thisCosmology, thisXray_grid, thisXray_spectrum, thisLya_grid, thisLya_spectrum, this21cmGrid, k10_table, Xe);
		}
		for(int i=0; i<local_n0; i++)
		{
			for(int j=0; j<nbins; j++)
			{
				for(int k=0; k<nbins; k++)
				{
					if(fmod(z,dz) == 0.)
					{
						double XHII = 1.-creal(this21cmGrid->XHI[i*nbins*nbins+j*nbins+k]);
						double XHeII = creal(this21cmGrid->XHeII[i*nbins*nbins+j*nbins+k]);
						double XHeIII = 1.-XHeII-creal(this21cmGrid->XHeI[i*nbins*nbins+j*nbins+k]);
						update_X_cell(cell, XHII, XHeII, XHeIII);
						
						double dens = creal(this21cmGrid->dens[i*nbins*nbins+j*nbins+k]);
						update_dens_cell(cell, dens);
						
						double temp = creal(this21cmGrid->temp[i*nbins*nbins+j*nbins+k]);
						update_temp_cell(cell, temp);
						
						double heatHI = fheat*creal(thisXray_grid->xray_heating_HI[i*nbins*nbins+j*nbins+k]);
						double heatHeI = fheat*creal(thisXray_grid->xray_heating_HeI[i*nbins*nbins+j*nbins+k]);
						double heatHeII = fheat*creal(thisXray_grid->xray_heating_HeII[i*nbins*nbins+j*nbins+k]);
						update_heat_cell(cell, heatHI, heatHeI, heatHeII);
						
						double photIonHI = creal(thisXray_grid->xray_ionization_HI[i*nbins*nbins+j*nbins+k]) + heatHI*fion_HI*tmp_HI;
						double photIonHeI = creal(thisXray_grid->xray_ionization_HeI[i*nbins*nbins+j*nbins+k]) + heatHeI*fion_HeI*tmp_HeI;
						double photIonHeII = creal(thisXray_grid->xray_ionization_HeII[i*nbins*nbins+j*nbins+k]) + heatHeII*fion_HeII*tmp_HeII;
						update_photIon_cell(cell, photIonHI, photIonHeI, photIonHeII);
						
// 						printf("photHI = %e\n", cell->photIonHI);
					}

					calc_step(recomb_rates, cell, dim_matrix, z, dz_tmp);
					printf("%d %d %d: z = %e:\tT = %e\t T = T0/(1+z)^2 = %e\t XHII = %e\t XHeII = %e\t XHeIII = %e\n", i,j,k,z, cell->temp, T_CMB(zstart)/((1.+zstart)*(1.+zstart))*((1.+z-dz_tmp)*(1.+z-dz_tmp)), cell->XHII, cell->XHeII, cell->XHeIII);
					
					if(fmod(z,dz) == 0.)
					{
						update_21cmgrid(this21cmGrid, cell, i, j, k);
					}
				}
			}
		}
	}
	  
	  
/*-------------------------------------------------------------------------------------*/
/* Deallocate ionization and temperature evolution structs */
/*-------------------------------------------------------------------------------------*/
	deallocate_recomb(recomb_rates);
	deallocate_cell(cell);
	
/*-------------------------------------------------------------------------------------*/
/* Deallocate 21cm computation structs */
/*-------------------------------------------------------------------------------------*/
	deallocate_cosmology(thisCosmology);
	
	deallocate_xray_grid(thisXray_grid);
	deallocate_xray_spectrum(thisXray_spectrum);
	
	deallocate_lya_grid(thisLya_grid);
	deallocate_lya_spectrum(thisLya_spectrum);
	
	deallocate_21cmgrid(this21cmGrid);
	
	deallocate_k10_table(k10_table);
}

