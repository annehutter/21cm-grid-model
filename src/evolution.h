#ifndef EVOLUTION_H
#define EVOLUTION_H
#endif


void read_update_density_grid(int snap, int double_precision, grid_21cm_t *this21cmGrid, confObj_t simParam);
void read_update_xray_grid(int snap, int double_precision, xray_grid_t *thisXray_grid, confObj_t simParam);
void update_xray_spectrum(xray_spectrum_t *thisXray_spectrum, double lumX, double alphaX, double nuX_min);
void read_update_lya_grid(int snap, int double_precision, lya_grid_t *thisLya_grid, confObj_t simParam);
void update_lya_spectrum(lya_spectrum_t *thisLya_spectrum, double lumLya, double alphaLya, double nuLya_min);
void update_21cmgrid(grid_21cm_t *this21cmGrid, cell_t *thisCell, int x, int y, int z);
void BubbleSort(double a[], int array_size);
double *create_redshift_table(inputlist_t *thisInputlist, double zstart, double zend, double dz);
void evolve(confObj_t simParam, int myRank);
