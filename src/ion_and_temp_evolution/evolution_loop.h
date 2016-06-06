#ifndef EVOLUTION_LOOP_H
#define EVOLUTION_LOOP_H
#endif

void do_ion_integration_step(int dim_matrix, double X_old[], double dt, double X[], recomb_t *thisRecombRates, cell_t *thisCell, int N);
void do_temp_integration_step(double *temp_old, double dt_sub, double dt, double n, double dn, cell_t *thisCell, int N);
void calc_step(recomb_t *thisRecombRates, cell_t *thisCell, int dim_matrix, double z, double dz);
double dt_from_dz(double z, double dz);
double t_at_z(double z);