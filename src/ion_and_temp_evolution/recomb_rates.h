#ifndef RECOMB_RATES_H
#define RECOMB_RATES_H
#endif

typedef struct
{
	double temp_max;
	double saturation;
	int lower_steps;
	int upper_steps;
	int N;
	
	double *recHII;
	double *recHeII;
	double *recHeIII;
} recomb_t;

double recHII(double temp);
double recHeII(double temp);
double recHeIII(double temp);

recomb_t *calcRecRate();
void deallocate_recomb(recomb_t *thisRecombRates);
void extend_recomb(recomb_t *thisRecombRates, double new_temp_max);
int temp_recomb_index(recomb_t *thisRecombRates, double temperature);