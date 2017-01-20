#ifndef SPIN_TEMPERATURE_HE
#define SPIN_TEMPERATURE_HE
#endif

double coupling_coll_He(double z, double ne, double Tk);
double calc_Tk_He(xray_spectrum_t *thisXraySpectrum);
void compute_Ts_He_on_grid(xray_spectrum_t *thisXraySpectrum, grid_3cm_t *this3cmGrid, grid_21cm_t *this21cmGrid, cosmology_t *thisCosmology);

