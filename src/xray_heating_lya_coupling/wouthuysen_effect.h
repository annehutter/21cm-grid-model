#ifndef WOUTHUYSEN_EFFECT_H
#define WOUTHUYSEN_EFFECT_H
#endif

typedef struct
{
	double emission;
	double alpha;
	double nu_min;
	
	double dens;
	double Xe;
	double Y;
	
	double x;
	
	double fheat;
} lya_params_t;

typedef struct
{
	double lum;
	double alpha;
	double nu_min;
} lya_spectrum_t;

typedef struct
{
	int nbins;
	float box_size;
	
	fftw_complex *lya;
	double mean_lya;
	
	fftw_complex *lya_lum;
	
	int local_n0;
	int local_0_start;
} lya_grid_t;


double lya_spectrum(lya_spectrum_t *thisSpectrum, double nu);
void lya_filter(lya_grid_t *thisLya_grid, double n, fftw_complex *filter, double h, double omega_m, double z, lya_spectrum_t *thisSpectrum);
void lya_convolve_fft(lya_grid_t *thisLya_grid, fftw_complex *filter, fftw_complex *output, fftw_complex *input);
void lya_add_modes(lya_grid_t *thisLya_grid, fftw_complex *thisLya_mode);
void lya_bg_source_emission(lya_grid_t *thisLya_grid, lya_spectrum_t*thisSpectrum, cosmology_t *thisCosmology);
void lya_bg_lya_excitation(lya_grid_t *thisLya_grid, xray_grid_t *thisXray_grid, cosmology_t *thisCosmology);
void lya_wouthuysen_coupling(lya_grid_t *thisLya_grid, lya_spectrum_t *thisSpectrum, xray_grid_t *thisXray_grid, cosmology_t *thisCosmology);

double TCMB(double z);

/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR LYA_PARAMS_T */
/*-------------------------------------------------------------------------------------*/

lya_params_t *initLya_params();
void deallocate_lya_params(lya_params_t * thisLya_params);
lya_spectrum_t *allocate_lya_spectrum(double lum, double alpha, double nu_min);

/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR LYA_SPECTRUM_T */
/*-------------------------------------------------------------------------------------*/

lya_spectrum_t *initLya_spectrum();
void deallocate_lya_spectrum(lya_spectrum_t * thisLya_spectrum);

/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR LYA_GRID_T */
/*-------------------------------------------------------------------------------------*/

lya_grid_t *initLya_grid();
void deallocate_lya_grid(lya_grid_t * thisLya_grid);
lya_grid_t *allocate_lya_grid(int nbins, float box_size);
void initialize_lya_grid(lya_grid_t *thisLya_grid);
void read_lum_lyagrid(lya_grid_t *thisLya_grid, char *filename, int double_precision);