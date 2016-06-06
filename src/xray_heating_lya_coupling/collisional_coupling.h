#ifndef COLLISIONAL_COUPLING_H
#define COLLISIONAL_COUPLING_H
#endif

typedef struct
{
	int num_HH;
	double *T_HH;
	double *k10_HH;
	
	int num;
	double *T;
	double *k10_eH;
	double *k10_pH;
} k10_t;


k10_t *create_k10_table();

void deallocate_k10_table(k10_t *thisTable);

double get_k10_HH(k10_t *thisTable, double temp);
double get_k10_eH(k10_t *thisTable, double temp);
double get_k10_pH(k10_t *thisTable, double temp);