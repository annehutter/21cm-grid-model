#ifndef SOLVE_ION_TEMP_H
#define SOLVE_ION_TEMP_H
#endif

double update_ion(grid_21cm_t *this21cmGrid, recomb_t *thisRecombRates, int index, double z, double z_old, double gamma, double n, double Hubble);
double update_temp(grid_21cm_t *this21cmGrid, int index, double sqrt_cub_z, double sqrt_cub_z_old, double factor, double xray_heating, double Hubble);
double update_dens(grid_21cm_t *this21cmGrid, int index, double D);
void overwrite_old_values_21cmgrid(grid_21cm_t *this21cmGrid);
void compute_temp_ion_grid(grid_21cm_t *this21cmGrid, recomb_t *thisRecombRates, cosmology_t *thisCosmology, xray_grid_t *thisXray_grid);