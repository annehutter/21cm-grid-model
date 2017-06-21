#ifndef XRAY_HEATING_H
#define XRAY_HEATING_H
#endif

typedef struct
{
	double emission;
	double nu_min;
	double alphaX;
	
	double dens;
	double Xe;
	double Y;
	
	double x;
	double nu;
	
	double z;
	double zp;
	
	double h;
	double omega_b;
	double omega_m;
	double omega_l;
} xray_params_t;

typedef struct
{
	double lumX;
	double alphaX;
	double nu_min;
} xray_spectrum_t;


typedef struct
{
	int nbins;
	float box_size;
	
	fftw_complex *xray_heating_HI;
	double mean_xray_heating_HI;
	fftw_complex *xray_heating_HeI;
	double mean_xray_heating_HeI;
	fftw_complex *xray_heating_HeII;
	double mean_xray_heating_HeII;
	
	fftw_complex *xray_ionization_HI;
	double mean_xray_ionization_HI;
	fftw_complex *xray_ionization_HeI;
	double mean_xray_ionization_HeI;
	fftw_complex *xray_ionization_HeII;
	double mean_xray_ionization_HeII;
	
	fftw_complex *xray_lum;
	
	int local_n0;
	int local_0_start;
} xray_grid_t;

// double freq_tau_equal_unity(double f_HI, double f1_HeI, double f2_HeI, double f_HeII);
// double distance_to_redshift(double z, double x, double omega_m);
// void xray_build_filter_functions(double *xray_filter_heating, double *xray_filter_ionization, double z);

double xray_mean_free_path_function(double zd, void *p);
double xray_mfp_calc_integral(xray_params_t params, double lowLim, double upLim);
// double xray_mean_free_path(double nu, double x, double dens, double Xe, double Y);

double xray_heating_function_HI(double nu, void *p);
double xray_heating_function_HeI(double nu, void *p);
double xray_heating_function_HeII(double nu, void *p);

double xray_heating_HI_calc_integral(xray_params_t params, double lowLim, double upLim);
double xray_heating_HeI_calc_integral(xray_params_t params, double lowLim, double upLim);
double xray_heating_HeII_calc_integral(xray_params_t params, double lowLim, double upLim);

double xray_ionization_function_HI(double nu, void *p);
double xray_ionization_function_HeI(double nu, void *p);
double xray_ionization_function_HeII(double nu, void *p);

double xray_ionization_HI_calc_integral(xray_params_t params, double lowLim, double upLim);
double xray_ionization_HeI_calc_integral(xray_params_t params, double lowLim, double upLim);
double xray_ionization_HeII_calc_integral(xray_params_t params, double lowLim, double upLim);


void xray_build_filter_functions(xray_grid_t *thisXray_grid, double *xray_filter_heating, double *xray_filter_ionization, xray_params_t *xray_params, int type);
void xray_filter(xray_grid_t *thisXray_grid, double *xray_filter_function, fftw_complex *filter);

void xray_heating_and_ionization(xray_grid_t *thisXray_grid, cosmology_t *thisCosmology, double Xe, xray_spectrum_t *thisSpectrum, int myRank);

void xray_heating_and_ionization_global(xray_grid_t *thisXray_grid, cosmology_t *thisCosmology, double Xe, xray_spectrum_t *thisSpectrum, int myRank);


/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR XRAY_PARAMS_T */
/*-------------------------------------------------------------------------------------*/

xray_params_t *initXray_params();
void deallocate_xray_params(xray_params_t * thisXray_params);


/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR XRAY_SPECTRUM_T */
/*-------------------------------------------------------------------------------------*/

xray_spectrum_t *initXray_spectrum();
void deallocate_xray_spectrum(xray_spectrum_t * thisXray_spectrum);
xray_spectrum_t *allocate_xray_spectrum(double lumX, double alphaX, double nu_min);


/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR XRAY_GRID_T */
/*-------------------------------------------------------------------------------------*/

xray_grid_t *initXray_grid();
void deallocate_xray_grid(xray_grid_t * thisXray_grid);
xray_grid_t *allocate_xray_grid(int nbins, float box_size);
void initialize_xray_grid(fftw_complex *thisArray, int nbins, int local_n0, double value);

void read_lum_xraygrid(xray_grid_t *thisXray_grid, char *filename, int double_precision);
