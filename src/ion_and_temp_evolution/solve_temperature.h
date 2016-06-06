#ifndef SOLVE_TEMPERATURE_H
#define SOLVE_TEMPERATURE_H
#endif

double create_prod_rate_temp();
double create_destr_rate_temp(double temperature, double n, double dn, double dt);
double Patankar_step(double oldvalue, double intervalue, double dt, double prod_rate, double destr_rRate);
double PatankarRKnoncon(double oldvalue, double dt, double dt_max, double n, double dn, cell_t *thisCell);

