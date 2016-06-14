#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>	//included because M_PI is not defined in <math.h>
#include <assert.h>
#include <complex.h>
#include <gsl/gsl_integration.h>

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
#include "redshift_tools.h"
#include "cross_sections.h"
#include "xray_heating.h"
#include "convolution_fftw.h"

#define SQR(X) ((X) * (X))
#define CUB(X) ((X) * (X) * (X))

// /*-------------------------------------------------------------------------------------*/
// /*-------------------------------------------------------------------------------------*/
// 
// double freq_tau_equal_unity(double f_HI, double f1_HeI, double f2_HeI, double f_HeII)
// {
// 	const double tmp, tmp2;
// 	
// // 	tmp = (f_HI*CUB(nu_HI) + f2_HeI*CUB(nu_HeI));
// // 	tmp2 = pow(-27.*tmp + sqrt(729.*SQR(tmp)-108.*CUB(f1_HeI)*pow(nu_HeI,6)),1./3.);
// // 	
// // 	return -pow(2.,1.3.)*f1_HeI*SQR(nu_HeI)/tmp2 - tmp2/(3.*pow(2.,1./3.));
// 	
// 	
// 	tmp = (f_HI*CUB(nu_HI) + f2_HeI*CUB(nu_HeI) + f_HeII*CUB(nu_HeIII));
// 	tmp2 = sqrt(-108.*CUB(f1_HeI)*pow(nu_HeI,6) + 729.*SQR(tmp))/3. - 9.*tmp;
// 	
// 	return (2.*pow(3.,1./3.)*f1_HeI*SQR(nu_HeI) + pow(2.*SQR(tmp2),1./3.))/(pow(SQR(6.)*tmp2,1./3.));	//this solution is approximative as it assumes for cross_sec_HeI the exponents -2 and -3, instead of -2.05 and -3.05
// }
// 
// /* computes 1+z', whereas z' is redshift that has a comoving distance x from redshift z */
// double distance_to_redshift(double z, double x, double omega_m)
// {
// 	pow(1/sqrt(1.+z)-H0*sqrt(omega_m)/(2.*clight_cm)*x ,-2);
// }
// 
// /* function to compute the filters for xray heating and ionization */
// void xray_build_filter_functions(double *xray_filter_heating, double *xray_filter_ionization, double z)
// {
// 	for(int i=0; i<nbins; i++)
// 	{
// 		/* compute comoving transverse distance from source */
// 		double x = (thisXray_grid->box_size/simParam->h/(1.+z))*i/nbins;
// 		
// 		/* compute 1+z', z' is redshift at x */
// 		const double tmp = distance_to_redshift(z, x, omega_m);
// 		
// 		/* compute frequency at which tau == 1 */
// 		const double tmp_HI = pow(1.+z,1.5) * (pow(tmp/(1.+z), -1.5) - 1.);
// 		const double tmp_HeI_1 = pow(1.+z,1.5) * (pow(tmp/(1.+z), -0.55) - 1.);
// 		const double tmp_HeI_2 = pow(1.+z,1.5) * (pow(tmp/(1.+z), -1.55) - 1.);
// 
// 		const double f_HI = (2.*nHI*cross_sec_HI(nu_HI))/(H0*sqrt(omega_m)*3.)*tmp_HI;
// 		const double f_HeI_1 = (2.*nHeI*7.2e-18*1.66)/(H0*sqrt(omega_m)*1.05)*tmp_HeI_1;
// 		const double f_HeI_2 = (2.*nHeI*7.2e-18*0.66)/(H0*sqrt(omega_m)*3.05)*tmp_HeI_2;
// 		const double f_HeII = (2.*nHeII*cross_sec_HeII(nu_HeII))/(H0*sqrt(omega_m)*3.)*tmp_HI;
// 
// 		double freq_tau_unity = freq_tau_equal_unity(f_HI, f_HeI_1, f_HeI_2, f_HeII);
// 		double nu_min;
// 		
// 		/* compute the contributions of all species to heating & primary and secondary ionization*/
// 		if(freq_tau_unity > nu_HI) nu_min = freq_tau_unity;
// 		else nu_min = nu_HI;
// 		const double contribution_HI = cross_sec_HI(nu_HI)*pow(nu_min/nu_HI,-3) * planck_cgs*pow(nu_min,1.-alphaX)*emission*((alphaX-4.)*nu_min-(alphaX-5.)*nu_HI)/((alphaX-4)*(alphaX-5.));
// 		
// 		const double contribution_ion_HI = cross_sec_HI(nu_HI)*pow(nu_min/nu_HI,-3) * planck_cgs*pow(nu_min,1.-alphaX)*emission/(alphaX-2.);
// 
// 		if(freq_tau_unity > nu_HeI) nu_min = freq_tau_unity;
// 		else nu_min = nu_HeI;
// 		const double contribution_HeI_1 = 7.2e-18*1.66*pow(nu_min/nu_HeI,-2.05) * planck_cgs*pow(nu_min,1.-alphaX)*emission*((alphaX-3.05)*nu_min-(alphaX-4.05)*nu_HeI)/((alphaX-3.05)*(alphaX-4.05));
// 
// 		const double contribution_ion_HeI_1 = 7.2e-18*1.66*pow(nu_min/nu_HeI,-2.05) * planck_cgs*pow(nu_min,1.-alphaX)*emission*/(alphaX-3.05);
// 		
// 		const double contribution_HeI_2 = 7.2e-18*0.66*pow(nu_min/nu_HeI,-3.05) * planck_cgs*pow(nu_min,1.-alphaX)*emission*((alphaX-4.05)*nu_min-(alphaX-5.05)*nu_HeI)/((alphaX-4.05)*(alphaX-5.05));
// 		
// 		const double contribution_ion_HeI_2 = 7.2e-18*0.66*pow(nu_min/nu_HeI,-3.05) * planck_cgs*pow(nu_min,1.-alphaX)*emission*(alphaX-4.05);
// 		
// 		if(freq_tau_unity > nu_HeII) nu_min = freq_tau_unity;
// 		else nu_min = nu_HeII;
// 		const double contribution_HeII = cross_sec_HeII(nu_HeII)*pow(nu_min/nu_HeII,-3) * planck_cgs*pow(nu_min,1.-alphaX)*emission*((alphaX-4.)*nu_min-(alphaX-5.)*nu_HeII)/((alphaX-4)*(alphaX-5.));
// 		
// 		xray_filter_heating[i] = (contribution_HI + contribution_HeI_1 + contribution_HeI_2 + contribution_HeII)/(4.*M_PI*SQR(x));
// 		
// 		xray_filter_ionization[i] = ((contribution_HI + contribution_HeI_1 + contribution_HeI_2 + contribution_HeII) + contribution_ion_HI + contribution_ion_HeI_1 + contribution_ion_HeI_2 + contribution_ion_HeII)/(4.*M_PI*SQR(x));
// 	}
// }

/*-------------------------------------------------------------------------------------*/
/* compute x-ray mean free path -------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/

double xray_mean_free_path_function(double zd, void *p)
{
	xray_params_t * params = (xray_params_t *)p;
	
	double h = params->h;
	double omega_b = params->omega_b;
	double omega_m = params->omega_m;
	double omega_l = params->omega_l;
	
	double nu = params->nu;
	double zp = params->zp;
	double Y = params->Y;
	double Xe = params->Xe;
	
	double tmp = 3.*clight_cm*(h*100.*km_cm/Mpc_cm)*omega_b/((8.*M_PI*G)*mp_g*sqrt(omega_m*CUB(1.+zd)+omega_l));
	double nu_zd = nu*(1.+zd)/(1.+zp);
	
	return tmp*SQR(1.+zd)*((1.-Y)*(1.-Xe)*cross_sec_HI(nu_zd) + Y*( (1.-Xe)*cross_sec_HeI(nu_zd) +  Xe*cross_sec_HeII(nu_zd)));
}

double xray_mfp_calc_integral(xray_params_t params, double lowLim, double upLim)
{
	gsl_function F;
	F.function = &xray_mean_free_path_function;
	F.params = &params;
	double result, error;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	gsl_integration_qags(&F, lowLim, upLim, 1.e-9, 1.e-9, 10000, w, &result, &error);
	
	gsl_integration_workspace_free(w);
	return result;
}

// double xray_mean_free_path(double nu, double x, double dens, double Xe, double Y)
// {
// 	return x*dens*((1.-Y)*(1.-Xe)*cross_sec_HI(nu) + Y*( (1.-Xe)*cross_sec_HeI(nu) +  Xe*cross_sec_HeII(nu)));
// }

/*-------------------------------------------------------------------------------------*/
/* Integral over frequency for photoionization heating rate ---------------------------*/
/*-------------------------------------------------------------------------------------*/

double xray_heating_function_HI(double nu, void *p)
{
	xray_params_t *params = (xray_params_t *)p;
	
	params->nu = nu;	//nu is in reference frame of zp (at absorption site)
	double factor = (1.+params->z)/(1.+params->zp);
	double mfp = xray_mfp_calc_integral(*params, params->z, params->zp);
	
	if(mfp <= 1. && nu >= params->nu_min/factor){	// nu_min is in reference frame of emitter, i.e. z
		return params->emission*pow(nu*factor,params->alphaX)*planck_cgs*(nu-nu_HI)*cross_sec_HI(nu);
	}else{
		return 0.;
	}
}

double xray_heating_function_HeI(double nu, void *p)
{
	xray_params_t *params = (xray_params_t *)p;
	
	params->nu = nu;	//nu is in reference frame of zp (at absorption site)
	double factor = (1.+params->z)/(1.+params->zp);
	double mfp = xray_mfp_calc_integral(*params, params->z, params->zp);
	
	if(mfp <= 1. && nu >= params->nu_min/factor){	// nu_min is in reference frame of emitter, i.e. z
		return params->emission*pow(nu*factor,params->alphaX)*planck_cgs*(nu-nu_HeI)*cross_sec_HeI(nu);
	}else{
		return 0.;
	}
}

double xray_heating_function_HeII(double nu, void *p)
{
	xray_params_t *params = (xray_params_t *)p;
	
	params->nu = nu;	//nu is in reference frame of zp (at absorption site)
	double factor = (1.+params->z)/(1.+params->zp);
	double mfp = xray_mfp_calc_integral(*params, params->z, params->zp);
	
	if(mfp <= 1. && nu >= params->nu_min/factor){	// nu_min is in reference frame of emitter, i.e. z
		return params->emission*pow(nu*factor,params->alphaX)*planck_cgs*(nu-nu_HeII)*cross_sec_HeII(nu);
	}else{
		return 0.;
	}
}

double xray_heating_HI_calc_integral(xray_params_t params, double lowLim, double upLim)
{
	gsl_function F;
	F.function = &xray_heating_function_HI;
	F.params = &params;
	double result, error;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	gsl_integration_qag(&F, lowLim, upLim, 1.e-9, 1.e-9, 10000, 1, w, &result, &error);
	
	gsl_integration_workspace_free(w);
	return result;
}

double xray_heating_HeI_calc_integral(xray_params_t params, double lowLim, double upLim)
{
	gsl_function F;
	F.function = &xray_heating_function_HeI;
	F.params = &params;
	double result, error;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	gsl_integration_qag(&F, lowLim, upLim, 1.e-9, 1.e-9, 10000, 1, w, &result, &error);
	
	gsl_integration_workspace_free(w);
	return result;
}

double xray_heating_HeII_calc_integral(xray_params_t params, double lowLim, double upLim)
{
	gsl_function F;
	F.function = &xray_heating_function_HeII;
	F.params = &params;
	double result, error;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	gsl_integration_qag(&F, lowLim, upLim, 1.e-9, 1.e-9, 10000, 1, w, &result, &error);
	
	gsl_integration_workspace_free(w);
	return result;
}

/*-------------------------------------------------------------------------------------*/
/* Integral over frequency for photoionization and secondary ionization ---------------*/
/*-------------------------------------------------------------------------------------*/

double xray_ionization_function_HI(double nu, void *p)
{
	xray_params_t *params = (xray_params_t *)p;
	
	params->nu = nu;	//nu is in reference frame of zp (at absorption site)
	double factor = (1.+params->z)/(1.+params->zp);
	double mfp = xray_mfp_calc_integral(*params, params->z, params->zp);
	
	if(mfp <= 1. && nu >= params->nu_min/factor){	// nu_min is in reference frame of emitter, i.e. z
		return params->emission*pow(nu*factor,params->alphaX)*cross_sec_HI(nu);
	}else{
		return 0.;
	}
}

double xray_ionization_function_HeI(double nu, void *p)
{
	xray_params_t *params = (xray_params_t *)p;
	
	params->nu = nu;	//nu is in reference frame of zp (at absorption site)
	double factor = (1.+params->z)/(1.+params->zp);
	double mfp = xray_mfp_calc_integral(*params, params->z, params->zp);
	
	if(mfp <= 1. && nu >= params->nu_min/factor){	// nu_min is in reference frame of emitter, i.e. z
		return params->emission*pow(nu*factor,params->alphaX)*cross_sec_HeI(nu);
	}else{
		return 0.;
	}
}

double xray_ionization_function_HeII(double nu, void *p)
{
	xray_params_t *params = (xray_params_t *)p;
	
	params->nu = nu;	//nu is in reference frame of zp (at absorption site)
	double factor = (1.+params->z)/(1.+params->zp);
	double mfp = xray_mfp_calc_integral(*params, params->z, params->zp);
	
	if(mfp <= 1. && nu >= params->nu_min/factor){	// nu_min is in reference frame of emitter, i.e. z
		return params->emission*pow(nu*factor,params->alphaX)*cross_sec_HeII(nu);
	}else{
		return 0.;
	}
}

double xray_ionization_HI_calc_integral(xray_params_t params, double lowLim, double upLim)
{
	gsl_function F;
	F.function = &xray_ionization_function_HI;
	F.params = &params;
	double result, error;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	gsl_integration_qag(&F, lowLim, upLim, 1.e-9, 1.e-9, 10000, 1, w, &result, &error);
	
	gsl_integration_workspace_free(w);
	return result;
}

double xray_ionization_HeI_calc_integral(xray_params_t params, double lowLim, double upLim)
{
	gsl_function F;
	F.function = &xray_ionization_function_HeI;
	F.params = &params;
	double result, error;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	gsl_integration_qag(&F, lowLim, upLim, 1.e-9, 1.e-9, 10000, 1, w, &result, &error);
	
	gsl_integration_workspace_free(w);
	return result;
}

double xray_ionization_HeII_calc_integral(xray_params_t params, double lowLim, double upLim)
{
	gsl_function F;
	F.function = &xray_ionization_function_HeII;
	F.params = &params;
	double result, error;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	gsl_integration_qag(&F, lowLim, upLim, 1.e-9, 1.e-9, 10000, 1, w, &result, &error);
	
	gsl_integration_workspace_free(w);
	return result;
}

/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/

void xray_build_filter_functions(xray_grid_t *thisXray_grid, double *xray_filter_heating, double *xray_filter_ionization, xray_params_t *xray_params, int type)
{
	int nbins = thisXray_grid->nbins;
	
	double h = xray_params->h;
	double omega_m = xray_params->omega_m;
	double z = xray_params->zp;	//redshift of frame of absorption
	
	double eps = (thisXray_grid->box_size/h)*epsilon/nbins*Mpc_cm;
	
	for(int i=0; i<nbins; i++)
	{
		const double x = (thisXray_grid->box_size/h)*i/nbins*Mpc_cm;	//in comoving cm
		xray_params->x = x;
			
		/* compute 1+z', z' is redshift at x */
		const double tmp = z_distance_to_redshift(z, x, h, omega_m);
		xray_params->z = tmp - 1.;	//redshift of emitter frame
		
		if(type == 0)
		{
			xray_filter_heating[i] = xray_heating_HI_calc_integral(*xray_params, nu_HI, nu_HI*1.e4)*SQR(1.+z)/(4.*M_PI*SQR(x + eps));
			xray_filter_ionization[i] = xray_ionization_HI_calc_integral(*xray_params, nu_HI, nu_HI*1.e4)*SQR(1.+z)/(4.*M_PI*SQR(x + eps));
		}else if(type == 1)
		{
			xray_filter_heating[i] = xray_heating_HeI_calc_integral(*xray_params, nu_HeI, nu_HI*1.e4)*SQR(1.+z)/(4.*M_PI*SQR(x + eps));
			xray_filter_ionization[i] = xray_ionization_HeI_calc_integral(*xray_params, nu_HeI, nu_HI*1.e4)*SQR(1.+z)/(4.*M_PI*SQR(x + eps));
		}else if(type ==2)
		{
			xray_filter_heating[i] = xray_heating_HeII_calc_integral(*xray_params, nu_HeII, nu_HI*1.e4)*SQR(1.+z)/(4.*M_PI*SQR(x + eps));
			xray_filter_ionization[i] = xray_ionization_HeII_calc_integral(*xray_params, nu_HeII, nu_HI*1.e4)*SQR(1.+z)/(4.*M_PI*SQR(x + eps));
		}
// 		printf("x = %e\t z_em = %e\t heating[%d] = %e\t ionization[%d] = %e\n", x, tmp-1., i, xray_filter_heating[i], i, xray_filter_ionization[i]);
	}
}

/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/

void xray_filter(xray_grid_t *thisXray_grid, double *xray_filter_function, fftw_complex *filter)
{
	const int nbins = thisXray_grid->nbins;
	const int local_n0 = thisXray_grid->local_n0;
	const int half_nbins = nbins/2;

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
				
				const double expr = sqrt(sq_i_expr + sq_j_expr + sq_k_expr);
				const int expr_int = floor(expr);
				
				const double value = xray_filter_function[expr_int] + (xray_filter_function[expr_int+1]-xray_filter_function[expr_int])*(expr - floor(expr));
// 				if(i==0 && j==0) printf("expr_int = %d \tf(expr_int) = %e\t value = %e\n", expr_int, xray_filter_function[expr_int], value);
				
				filter[i*nbins*nbins+j*nbins+k] = value + 0.*I;
			}
		}
	}
}

/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/

void xray_heating_and_ionization(xray_grid_t *thisXray_grid, cosmology_t *thisCosmology, double Xe, xray_spectrum_t *thisSpectrum)
{
	fftw_complex *xray_filter_heating;
	fftw_complex *xray_filter_ionization;
	
	double *xray_filter_function_heating;
	double *xray_filter_function_ionization;
	
	int nbins = thisXray_grid->nbins;
	
	xray_params_t *xray_params;
	
	xray_params = initXray_params();
	
	xray_params->emission = thisSpectrum->lumX;
	xray_params->nu_min = thisSpectrum->nu_min;
	xray_params->alphaX = thisSpectrum->alphaX;
	xray_params->dens = thisCosmology->nH_z + 4.*thisCosmology->nHe_z;
	xray_params->Xe = Xe;
	xray_params->Y = thisCosmology->Y;

	xray_params->x = 0.;
	xray_params->nu = 0.;
	
	xray_params->z = thisCosmology->z;
	xray_params->zp = thisCosmology->z;
	
	xray_params->h = thisCosmology->h;
	xray_params->omega_b = thisCosmology->omega_b;
	xray_params->omega_m = thisCosmology->omega_m;
	xray_params->omega_l = thisCosmology->omega_l;

// 	printf("xray heating and ionization: initialisation done\n");
	
	/* allocating space for filters for heating and ionization */
	xray_filter_heating = (fftw_complex*)malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(xray_filter_heating == NULL)
	{
		fprintf(stderr, "xray_filter_heating in xray_heating (xray_heating.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	xray_filter_ionization = (fftw_complex*)malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(xray_filter_heating == NULL)
	{
		fprintf(stderr, "xray_filter_ionization in xray_heating (xray_heating.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	
	xray_filter_function_heating = (double*)malloc(sizeof(double)*nbins);
	if(xray_filter_function_heating == NULL)
	{
		fprintf(stderr, "xray_filter_function_heating in xray_heating (xray_heating.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	xray_filter_function_ionization = (double*)malloc(sizeof(double)*nbins);
	if(xray_filter_function_ionization == NULL)
	{
		fprintf(stderr, "xray_filter_function_ionization in xray_heating (xray_heating.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	for(int type=0; type<3; type++)
	{
// 		printf("building filter function for type %d\n", type);
		/* build the filter function */
		xray_build_filter_functions(thisXray_grid, xray_filter_function_heating, xray_filter_function_ionization, xray_params, type);
		
// 		printf("filter functions built\n");

		/* generate 3D grid of filter */
		xray_filter(thisXray_grid, xray_filter_function_heating, xray_filter_heating);
		xray_filter(thisXray_grid, xray_filter_function_ionization, xray_filter_ionization);
		
// 		printf("filter array built\n");
		
		if(type == 0){
			convolve_fft(thisXray_grid, xray_filter_heating, thisXray_grid->xray_heating_HI, thisXray_grid->xray_lum);
			convolve_fft(thisXray_grid, xray_filter_ionization, thisXray_grid->xray_ionization_HI, thisXray_grid->xray_lum);
// #ifdef __MPI
// 			write_grid_to_file_float(xray_filter_heating, nbins, local_n0, thisXray_grid->local_0_start, "xray_heating_filter.dat");
// 			write_grid_to_file_float(thisXray_grid->xray_lum, nbins, local_n0, thisXray_grid->local_0_start, "xray_lum_HI.dat");
// #else
// 			write_grid_to_file_float(xray_filter_heating, nbins, local_n0, "xray_heating_filter.dat");
// 			write_grid_to_file_float(thisXray_grid->xray_lum, nbins, local_n0, "xray_heating_lum.dat");
// #endif
		}else if(type == 1){
			convolve_fft(thisXray_grid, xray_filter_heating, thisXray_grid->xray_heating_HeI, thisXray_grid->xray_lum);
			convolve_fft(thisXray_grid, xray_filter_ionization, thisXray_grid->xray_ionization_HeI, thisXray_grid->xray_lum);
		}else if(type == 2){
			convolve_fft(thisXray_grid, xray_filter_heating, thisXray_grid->xray_heating_HeII, thisXray_grid->xray_lum);
			convolve_fft(thisXray_grid, xray_filter_ionization, thisXray_grid->xray_ionization_HeII, thisXray_grid->xray_lum);
		}
// 		printf("convolution done\n");
	}
	
// #ifdef __MPI
// 	write_grid_to_file_float(thisXray_grid->xray_heating_HI, nbins, local_n0, thisXray_grid->local_0_start, "xray_heating_HI.dat");
// #else
// 	write_grid_to_file_float(thisXray_grid->xray_heating_HI, nbins, local_n0, "xray_heating_HI.dat");
// #endif

	fftw_free(xray_filter_heating);
	fftw_free(xray_filter_ionization);
	
	free(xray_filter_function_heating);
	free(xray_filter_function_ionization);
	
	deallocate_xray_params(xray_params);
	
// 	for(int i=0; i<local_n0; i++)
// 	{
// 		for(int j=0; i<nbins; j++)
// 		{
// 			for(int k=0; k<nbins; k++)
// 			{
// 				printf("xray_heating_HI = %e\t xray_ionization_HI = %e\n", creal(thisXray_grid->xray_heating_HI[i*nbins*nbins+j*nbins+k]), creal(thisXray_grid->xray_ionization_HI[i*nbins*nbins+j*nbins+k]));
// 				
// 				printf("xray_heating_HeI = %e\t xray_ionization_HeI = %e\n", creal(thisXray_grid->xray_heating_HeI[i*nbins*nbins+j*nbins+k]), creal(thisXray_grid->xray_ionization_HeI[i*nbins*nbins+j*nbins+k]));
// 
// 				printf("xray_heating_HeII = %e\t xray_ionization_HeII = %e\n\n", creal(thisXray_grid->xray_heating_HeII[i*nbins*nbins+j*nbins+k]), creal(thisXray_grid->xray_ionization_HeII[i*nbins*nbins+j*nbins+k]));
// 
// 			}
// 		}
// 	}
}

void xray_heating_and_ionization_global(xray_grid_t *thisXray_grid, cosmology_t *thisCosmology, double Xe, xray_spectrum_t *thisSpectrum)
{
	int nbins = thisXray_grid->nbins;
	int local_n0 = thisXray_grid->local_n0;
	
	xray_heating_and_ionization(thisXray_grid, thisCosmology, Xe, thisSpectrum);

	double sum_heating_HI = 0.;
	double sum_heating_HeI = 0.;
	double sum_heating_HeII = 0.;

	double sum_ionization_HI = 0.;
	double sum_ionization_HeI = 0.;
	double sum_ionization_HeII = 0.;

	/* compute mean across all cells */
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; i<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				sum_heating_HI += creal(thisXray_grid->xray_heating_HI[i*nbins*nbins+j*nbins+k]);
				sum_heating_HeI += creal(thisXray_grid->xray_heating_HeI[i*nbins*nbins+j*nbins+k]);
				sum_heating_HeII += creal(thisXray_grid->xray_heating_HeII[i*nbins*nbins+j*nbins+k]);

				sum_ionization_HI += creal(thisXray_grid->xray_ionization_HI[i*nbins*nbins+j*nbins+k]);
				sum_ionization_HeI += creal(thisXray_grid->xray_ionization_HeI[i*nbins*nbins+j*nbins+k]);
				sum_ionization_HeII += creal(thisXray_grid->xray_ionization_HeII[i*nbins*nbins+j*nbins+k]);
			}
		}
	}
#ifdef __MPI
	MPI_Allreduce(&sum_heating_HI, &thisXray_grid->mean_xray_heating_HI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sum_heating_HeI, &thisXray_grid->mean_xray_heating_HeI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sum_heating_HeII, &thisXray_grid->mean_xray_heating_HeII, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	MPI_Allreduce(&sum_ionization_HI, &thisXray_grid->mean_xray_ionization_HI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sum_ionization_HeI, &thisXray_grid->mean_xray_ionization_HeI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sum_ionization_HeII, &thisXray_grid->mean_xray_ionization_HeII, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	thisXray_grid->mean_xray_heating_HI /= (nbins*nbins*nbins);
	thisXray_grid->mean_xray_heating_HeI /= (nbins*nbins*nbins);
	thisXray_grid->mean_xray_heating_HeII /= (nbins*nbins*nbins);

	thisXray_grid->mean_xray_ionization_HI /= (nbins*nbins*nbins);
	thisXray_grid->mean_xray_ionization_HeI /= (nbins*nbins*nbins);
	thisXray_grid->mean_xray_ionization_HeII /= (nbins*nbins*nbins);
}

/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR XRAY_PARAMS_T */
/*-------------------------------------------------------------------------------------*/

xray_params_t *initXray_params()
{
	xray_params_t *newXray_params;
	
	newXray_params = malloc(sizeof(xray_params_t));
	if(newXray_params == NULL)
	{
		fprintf(stderr, "newXray_params in initXray_params (xray_heating.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	return newXray_params;
}

void deallocate_xray_params(xray_params_t * thisXray_params)
{
	free(thisXray_params);
}


/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR XRAY_SPECTRUM_T */
/*-------------------------------------------------------------------------------------*/

xray_spectrum_t *initXray_spectrum()
{
	xray_spectrum_t *newXray_spectrum;
	
	newXray_spectrum = malloc(sizeof(xray_spectrum_t));
	if(newXray_spectrum == NULL)
	{
		fprintf(stderr, "newXray_spectrum in initXray_spectrum (xray_heating.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	return newXray_spectrum;
}

void deallocate_xray_spectrum(xray_spectrum_t * thisXray_spectrum)
{
	free(thisXray_spectrum);
}

xray_spectrum_t *allocate_xray_spectrum(double lumX, double alphaX, double nu_min)
{
	xray_spectrum_t * thisXray_spectrum;
	
	thisXray_spectrum = initXray_spectrum();
	
	thisXray_spectrum->lumX = lumX;
	thisXray_spectrum->alphaX = alphaX;
	thisXray_spectrum->nu_min = nu_min;
	
	return thisXray_spectrum;
}


/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR XRAY_GRID_T */
/*-------------------------------------------------------------------------------------*/

xray_grid_t *initXray_grid()
{
	xray_grid_t *newXray_grid;
	
	newXray_grid = malloc(sizeof(xray_grid_t));
	if(newXray_grid == NULL)
	{
		fprintf(stderr, "newXray_grid in initXray_grid (xray_heating.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	newXray_grid->nbins = 0;
	newXray_grid->box_size = 0.;
	
	newXray_grid->xray_heating_HI = NULL;
	newXray_grid->mean_xray_heating_HI = 0.;
	newXray_grid->xray_heating_HeI = NULL;
	newXray_grid->mean_xray_heating_HeI = 0.;
	newXray_grid->xray_heating_HeII = NULL;
	newXray_grid->mean_xray_heating_HeII = 0.;
	
	newXray_grid->xray_ionization_HI = NULL;
	newXray_grid->mean_xray_ionization_HI = 0.;
	newXray_grid->xray_ionization_HeI = NULL;
	newXray_grid->mean_xray_ionization_HeI = 0.;
	newXray_grid->xray_ionization_HeII = NULL;
	newXray_grid->mean_xray_ionization_HeII = 0.;
	
	newXray_grid->xray_lum = NULL;
	
	newXray_grid->local_n0 = 0;
	newXray_grid->local_0_start = 0;
	
	return newXray_grid;
}

void deallocate_xray_grid(xray_grid_t * thisXray_grid)
{
	if(thisXray_grid->xray_heating_HI != NULL) fftw_free(thisXray_grid->xray_heating_HI);
	if(thisXray_grid->xray_heating_HeI != NULL) fftw_free(thisXray_grid->xray_heating_HeI);
	if(thisXray_grid->xray_heating_HeII != NULL) fftw_free(thisXray_grid->xray_heating_HeII);

	if(thisXray_grid->xray_ionization_HI != NULL) fftw_free(thisXray_grid->xray_ionization_HI);
	if(thisXray_grid->xray_ionization_HeI != NULL) fftw_free(thisXray_grid->xray_ionization_HeI);
	if(thisXray_grid->xray_ionization_HeII != NULL) fftw_free(thisXray_grid->xray_ionization_HeII);

	if(thisXray_grid->xray_lum != NULL) fftw_free(thisXray_grid->xray_lum);
	free(thisXray_grid);
}

xray_grid_t *allocate_xray_grid(int nbins, float box_size)
{
#ifdef __MPI
	ptrdiff_t alloc_local, local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	xray_grid_t *thisXray_grid;

	thisXray_grid = initXray_grid();
	
	thisXray_grid->nbins = nbins;
	thisXray_grid->box_size = box_size;
	
	thisXray_grid->local_n0 = nbins;
	thisXray_grid->local_0_start = 0;
#ifdef __MPI
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	
	thisXray_grid->local_n0 = local_n0;
	thisXray_grid->local_0_start = local_0_start;
	
	thisXray_grid->xray_heating_HI = fftw_alloc_complex(alloc_local);
	thisXray_grid->xray_heating_HeI = fftw_alloc_complex(alloc_local);
	thisXray_grid->xray_heating_HeII = fftw_alloc_complex(alloc_local);

	thisXray_grid->xray_ionization_HI = fftw_alloc_complex(alloc_local);
	thisXray_grid->xray_ionization_HeI = fftw_alloc_complex(alloc_local);
	thisXray_grid->xray_ionization_HeII = fftw_alloc_complex(alloc_local);

	thisXray_grid->xray_lum = fftw_alloc_complex(alloc_local);
	
#else
	thisXray_grid->xray_heating_HI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	thisXray_grid->xray_heating_HeI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	thisXray_grid->xray_heating_HeII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);

	thisXray_grid->xray_ionization_HI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	thisXray_grid->xray_ionization_HeI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
	thisXray_grid->xray_ionization_HeII = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);

	thisXray_grid->xray_lum = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*local_n0*nbins*nbins);
#endif
	
	return thisXray_grid;
}


void read_lum_xraygrid(xray_grid_t *thisXray_grid, char *filename, int double_precision)
{
#ifdef __MPI
	ptrdiff_t local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	int nbins;
	
	nbins = thisXray_grid->nbins;
	local_n0 = thisXray_grid->local_n0;
	
	if(double_precision == 1)
	{
#ifdef __MPI
	local_0_start = thisXray_grid->local_0_start;
	read_grid_doubleprecision(thisXray_grid->xray_lum, nbins, local_n0, local_0_start, filename);
#else
	read_grid_doubleprecision(thisXray_grid->xray_lum, nbins, local_n0, filename);
#endif
	}else{
#ifdef __MPI
	local_0_start = thisXray_grid->local_0_start;
	read_grid(thisXray_grid->xray_lum, nbins, local_n0, local_0_start, filename);
#else
	read_grid(thisXray_grid->xray_lum, nbins, local_n0, filename);
#endif
	}
}