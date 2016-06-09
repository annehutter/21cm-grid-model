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

#include "fftw_array_tools.h"
#include "phys_const.h"
#include "chem_const.h"
#include "cosmology.h"
#include "xray_heating.h"
#include "wouthuysen_effect.h"
#include "redshift_tools.h"
#include "cross_sections.h"

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

double lya_frecycle(double n)
{
	if(n == 2.){
		return 1.;
	}else if(n == 3.){
		return 0.;
	}else if(n == 4.){
		return 0.26;
	}else if(n == 5.){
		return 0.31;
	}else if(n == 6.){
		return 0.32;
	}else if(n == 7.){
		return 0.33;
	}else if(n < 10.){
		return 0.34;
	}else{
		return 0.35;
	}
}

double lya_spectrum(lya_spectrum_t *thisSpectrum, double nu)
{
	if(nu>thisSpectrum->nu_min)
	{
		return thisSpectrum->lum*pow(nu, thisSpectrum->alpha);
	}else{
		return 0.;
	}
}

void lya_filter(lya_grid_t *thisLya_grid, double n, fftw_complex *filter, double h, double omega_m, double z, lya_spectrum_t *thisSpectrum)
{
	const int nbins = thisLya_grid->nbins;
	const int local_n0 = thisLya_grid->local_n0;
	const int half_nbins = nbins/2;
	
	const double factor = (thisLya_grid->box_size/h)/nbins*Mpc_cm;	//in comoving cm
	const double sq_factor = factor*factor;

	for(int i=0; i<local_n0; i++)
	{
		const double i_expr = half_nbins-abs(i - half_nbins);		
		const double sq_i_expr = SQR(i_expr);
		for(int j=0; j<nbins; j++)
		{
			const double j_expr = half_nbins - abs(j - half_nbins);
			const double sq_j_expr = SQR(j_expr);
			for(int k=0; k<nbins; k++)
			{
				const double k_expr = half_nbins - abs(k - half_nbins);
				const double sq_k_expr = SQR(k_expr);
				
				const double expr = (sq_i_expr + sq_j_expr + sq_k_expr + epsilon)*sq_factor;
				
				double nu = nu_lyalimit*(1.-1./SQR(n))/SQR(1.-h*100.*km_cm/Mpc_cm*sqrt(omega_m*(1.+z))/(2.*clight_cm)*sqrt(expr));
// 				printf("nu = %e\t expr = %e\t%e\n", nu, sqrt(expr), h*100.*km_cm/Mpc_cm*sqrt(omega_m*(1.+z))/(2.*clight_cm)*sqrt(expr));
				
				filter[i*nbins*nbins+j*nbins+k] = lya_frecycle(n)*lya_spectrum(thisSpectrum, nu)*SQR(1.+z)/(4.*M_PI*expr) + 0.*I;
// 				if(i==0 && j==0) printf("%d: lya_filter = %e\n", k, creal(filter[i*nbins*nbins+j*nbins+k]));
			}
		}
	}
}

void lya_convolve_fft(lya_grid_t *thisLya_grid, fftw_complex *filter, fftw_complex *output, fftw_complex *input)
{
	int nbins;
	double factor;
#ifdef __MPI
	ptrdiff_t alloc_local, local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	
	fftw_complex *input_ft, *filter_ft;
	fftw_plan plan_input, plan_filter, plan_back;
	
	nbins = thisLya_grid->nbins;
	local_n0 = thisLya_grid->local_n0;
	
#ifdef __MPI
	local_0_start = thisLya_grid->local_0_start;
	
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	assert(local_0_start == thisLya_grid->local_0_start);
	assert(local_n0 == thisLya_grid->local_n0);
	
	input_ft = fftw_alloc_complex(alloc_local);
	filter_ft = fftw_alloc_complex(alloc_local);
	
	plan_input = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, input, input_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MPI_TRANSPOSED_OUT); //FFTW_MPI_TRANSPOSED_OUT
	plan_filter = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MPI_TRANSPOSED_OUT);
#else 
	input_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(input_ft == NULL)
	{
		fprintf(stderr, "input_ft in convolve_fft (filtering.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	filter_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(filter_ft == NULL)
	{
		fprintf(stderr, "filter_ft in convolve_fft (filtering.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	plan_input = fftw_plan_dft_3d(nbins, nbins, nbins, input, input_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_filter = fftw_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, FFTW_FORWARD, FFTW_ESTIMATE);
#endif
	
	fftw_execute(plan_input);
	fftw_execute(plan_filter);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				input_ft[i*nbins*nbins+j*nbins+k] = input_ft[i*nbins*nbins+j*nbins+k]*filter_ft[i*nbins*nbins+j*nbins+k];
			}
		}
	}
	
#ifdef __MPI
	plan_back = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, input_ft, output, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MPI_TRANSPOSED_IN);
#else
	plan_back = fftw_plan_dft_3d(nbins, nbins, nbins, input_ft, output, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
	fftw_execute(plan_back);
	
	factor = 1./(nbins*nbins*nbins);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				output[i*nbins*nbins+j*nbins+k] = factor*output[i*nbins*nbins+j*nbins+k];
			}
		}
	}	
	
	fftw_destroy_plan(plan_input);
	fftw_destroy_plan(plan_filter);
	fftw_destroy_plan(plan_back);
	
	fftw_free(input_ft);
	fftw_free(filter_ft);
}

void lya_add_modes(lya_grid_t *thisLya_grid, fftw_complex *thisLya_mode)
{
	int nbins = thisLya_grid->nbins;
	int local_n0 = thisLya_grid->local_n0;
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++) 
			{
				thisLya_grid->lya[i*nbins*nbins+j*nbins+k] = (creal(thisLya_grid->lya[i*nbins*nbins+j*nbins+k]) + thisLya_mode[i*nbins*nbins+j*nbins+k]) + 0.*I;
			}
		}
	}
}

void lya_bg_source_emission(lya_grid_t *thisLya_grid, lya_spectrum_t *thisSpectrum, cosmology_t *thisCosmology)
{
	int nbins;
	
  	/* construct grid for mode in lya bg */
#ifdef __MPI
	ptrdiff_t alloc_local, local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	fftw_complex *thisLya_mode;
	fftw_complex *filter;
	nbins = thisLya_grid->nbins;
	local_n0 = nbins;
	
#ifdef __MPI
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	
	thisLya_mode = fftw_alloc_complex(alloc_local);
	filter = fftw_alloc_complex(alloc_local);
#else 
	thisLya_mode = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(thisLya_mode == NULL)
	{
		fprintf(stderr, "thisLya_mode in lya_bg_source_emission (wouthuysen.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	filter = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(filter == NULL)
	{
		fprintf(stderr, "filter in lya_bg_source_emission (wouthuysen.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
#endif
	
	for(int n=2; n<24; n++)		//nmax from ionized sphere around sources (10pkpc at z=20 and 35pkpc at z=10)
	{
		/* construct filter */
		lya_filter(thisLya_grid, (double)n, filter, thisCosmology->h, thisCosmology->omega_m, thisCosmology->z, thisSpectrum);
		
		/* lya emission field = convolution between filter and source grid */
		lya_convolve_fft(thisLya_grid, filter, thisLya_mode, thisLya_grid->lya_lum);
		
		lya_add_modes(thisLya_grid, thisLya_mode);
	}
	
	fftw_free(thisLya_mode);
	fftw_free(filter);
}

void lya_bg_lya_excitation(lya_grid_t *thisLya_grid, xray_grid_t *thisXray_grid, cosmology_t *thisCosmology)
{
	int nbins = thisLya_grid->nbins;
	int local_n0 = thisLya_grid->local_n0;
	
	double factor = SQR(lambda_lya)/(4.*M_PI*planck_cgs*clight_cm*thisCosmology->Hubble_z);
	printf("factor = %e\n", factor);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				double const heating = creal(thisXray_grid->xray_heating_HI[i*nbins*nbins+j*nbins+k]) + creal(thisXray_grid->xray_heating_HeI[i*nbins*nbins+j*nbins+k]) + creal(thisXray_grid->xray_heating_HeII[i*nbins*nbins+j*nbins+k]);
				
				if(heating<0.) 
				{
					fprintf(stderr, "heating is negative!!! and has value = %e\nCheckyour x-ray heating!\nStopping execution.\n", heating);
					exit(EXIT_FAILURE);
				}
				thisLya_grid->lya[i*nbins*nbins+j*nbins+k] = creal(thisLya_grid->lya[i*nbins*nbins+j*nbins+k]) + heating*factor + 0.*I;
			}
		}
	}
}

void lya_wouthuysen_coupling(lya_grid_t *thisLya_grid, lya_spectrum_t *thisSpectrum, xray_grid_t *thisXray_grid, cosmology_t *thisCosmology)
{
	int nbins = thisLya_grid->nbins;
	int local_n0 = thisLya_grid->local_n0;
	
	/* compute comoving Lya flux */
	lya_bg_source_emission(thisLya_grid, thisSpectrum, thisCosmology);
	lya_bg_lya_excitation(thisLya_grid, thisXray_grid, thisCosmology);
	
	/* compute prefactor for Lya coupling */
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				thisLya_grid->lya[i*nbins*nbins+j*nbins+k] = creal(thisLya_grid->lya[i*nbins*nbins+j*nbins+k]) + 0.*I;
				
// 				printf("lya = %e\n", creal(thisLya_grid->lya[i*nbins*nbins+j*nbins+k]));
			}
		}
	}
}



/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR LYA_PARAMS_T */
/*-------------------------------------------------------------------------------------*/

lya_params_t *initLya_params()
{
	lya_params_t *newLya_params;
	
	newLya_params = malloc(sizeof(lya_params_t));
	if(newLya_params == NULL)
	{
		fprintf(stderr, "newLya_params in initLya_params (lya_heating.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	return newLya_params;
}

void deallocate_lya_params(lya_params_t * thisLya_params)
{
	free(thisLya_params);
}


/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR LYA_SPECTRUM_T */
/*-------------------------------------------------------------------------------------*/

lya_spectrum_t *initLya_spectrum()
{
	lya_spectrum_t *newLya_spectrum;
	
	newLya_spectrum = malloc(sizeof(lya_spectrum_t));
	if(newLya_spectrum == NULL)
	{
		fprintf(stderr, "newLya_spectrum in initLya_spectrum (lya_heating.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	return newLya_spectrum;
}

void deallocate_lya_spectrum(lya_spectrum_t * thisLya_spectrum)
{
	free(thisLya_spectrum);
}

lya_spectrum_t *allocate_lya_spectrum(double lum, double alpha, double nu_min)
{
	lya_spectrum_t *thisLya_spectrum;
	
	thisLya_spectrum = initLya_spectrum();
	
	thisLya_spectrum->lum = lum;
	thisLya_spectrum->alpha = alpha;
	thisLya_spectrum->nu_min = nu_min;
	
	return thisLya_spectrum;
}

/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR LYA_GRID_T */
/*-------------------------------------------------------------------------------------*/

lya_grid_t *initLya_grid()
{
	lya_grid_t *newLya_grid;
	
	newLya_grid = malloc(sizeof(lya_grid_t));
	if(newLya_grid == NULL)
	{
		fprintf(stderr, "newLya_grid in initLya_grid (lya_heating.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	newLya_grid->nbins = 0;
	newLya_grid->box_size = 0.;
	
	newLya_grid->lya = NULL;
	newLya_grid->mean_lya = 0.;
	
	newLya_grid->lya_lum = NULL;
	
	newLya_grid->local_n0 = 0;
	newLya_grid->local_0_start = 0;
	
	return newLya_grid;
}

void deallocate_lya_grid(lya_grid_t * thisLya_grid)
{
	if(thisLya_grid->lya != NULL) fftw_free(thisLya_grid->lya);
	if(thisLya_grid->lya_lum != NULL) fftw_free(thisLya_grid->lya_lum);
	free(thisLya_grid);
}

lya_grid_t *allocate_lya_grid(int nbins, float box_size)
{
#ifdef __MPI
	ptrdiff_t alloc_local, local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	lya_grid_t *thisLya_grid;

	thisLya_grid = initLya_grid();
	
	thisLya_grid->nbins = nbins;
	thisLya_grid->box_size = box_size;
	
	thisLya_grid->local_n0 = nbins;
#ifdef __MPI
	fftw_mpi_init();
	
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	
	thisLya_grid->local_n0 = local_n0;
	thisLya_grid->local_0_start = local_0_start;
	
	thisLya_grid->lya = fftw_alloc_complex(alloc_local);
	thisLya_grid->lya_lum = fftw_alloc_complex(alloc_local);
	
#else
	thisLya_grid->lya = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	thisLya_grid->lya_lum = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
#endif
	initialize_lya_grid(thisLya_grid);
	
	return thisLya_grid;
}

void initialize_lya_grid(lya_grid_t *thisLya_grid)
{
	int nbins = thisLya_grid->nbins;
	int local_n0 = thisLya_grid->local_n0;
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				thisLya_grid->lya[i*nbins*nbins+j*nbins+k] = 0. + 0.*I;
				thisLya_grid->lya_lum[i*nbins*nbins+j*nbins+k] = 0. + 0.*I;
			}
		}
	}
}

void read_lum_lyagrid(lya_grid_t *thisLya_grid, char *filename, int double_precision)
{
#ifdef __MPI
	ptrdiff_t local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	int nbins;
	
	nbins = thisLya_grid->nbins;
	local_n0 = thisLya_grid->local_n0;
	
	if(double_precision == 1)
	{
#ifdef __MPI
	local_0_start = thisLya_grid->local_0_start;
	read_grid_doubleprecision(thisLya_grid->lya_lum, nbins, local_n0, local_0_start, filename);
#else
	read_grid_doubleprecision(thisLya_grid->lya_lum, nbins, local_n0, filename);
#endif
	}else{
#ifdef __MPI
	local_0_start = thisLya_grid->local_0_start;
	read_grid(thisLya_grid->lya_lum, nbins, local_n0, local_0_start, filename);
#else
	read_grid(thisLya_grid->lya_lum, nbins, local_n0, filename);
#endif
	}
}