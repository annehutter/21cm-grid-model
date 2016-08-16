#ifndef SPIN_TEMPERATURE_H
#define SPIN_TEMPERATURE_H
#endif

double tau_gp(double Hubble_z_inv, double nHI);

double T_CMB(double z);

double coupling_coll(double z, double nH, k10_t *k10_table, double Tk);

double coupling_alpha(double prefac_z, double modSalpha, double Jalpha);
double coupling_alpha_prefac_z(double z, double prefac);

double calc_modSalpha(double Hubble_z_inv, double nHI, double Tk_inv, double Ts_inv);
double calc_xi(double Hubble_z_inv, double nHI, double Tk_inv);
double calc_Teff_inv(double Tk_inv, double Ts_inv);

double calc_Ts_inv(double xa_mod, double xc, double Teff_inv, double Tk_inv, double Tbg_inv);
double compute_Ts_inv(double Hubble_z_inv, double nH, double XHI, double coupling_alpha_prefac_z, double Jalpha, double xc, double Tk_inv, double Tbg_inv);
void compute_Ts_on_grid(lya_grid_t *thisLya_grid, k10_t *k10_table, grid_21cm_t *this21cmGrid, cosmology_t *thisCosmology);


/* wouthuysen coupling */
void lya_wouthuysen_coupling(lya_grid_t *thisLya_grid, lya_spectrum_t *thisSpectrum, xray_grid_t *thisXray_grid, grid_21cm_t *this21cmGrid, cosmology_t *thisCosmology);
void lya_bg_lya_excitation(lya_grid_t *thisLya_grid, xray_grid_t *thisXray_grid, grid_21cm_t *this21cmGrid, cosmology_t *thisCosmology);
