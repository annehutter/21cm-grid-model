#ifndef REDSHIFT_TOOLS_H
#define REDSHIFT_TOOLS_H
#endif

 
double z_hubble(double h, double omega_m, double omega_l, double z);

double z_time_from_redshiftinterval_flatuniverse(double h, double omega_m, double omega_l, double zmin, double zmax);

double z_distance_to_redshift(double z, double x, double h, double omega_m);

double z_calc_numdensity_per_mp_at_redshift(double h, double omega_b, double z);
double z_calc_numdensity_H_at_redshift(double h, double omega_b, double z, double Y);
double z_calc_numdensity_He_at_redshift(double h, double omega_b, double z, double Y);