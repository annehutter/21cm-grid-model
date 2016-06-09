#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>	//included because M_PI is not defined in <math.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "collisional_coupling.h"

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

k10_t *create_k10_table()
{
	k10_t *thisTable;
	
	thisTable = malloc(sizeof(k10_t));
	
	thisTable->num_HH= 19;
	thisTable->T_HH = malloc(thisTable->num*sizeof(double));
	thisTable->k10_HH = malloc(thisTable->num*sizeof(double));
	
	thisTable->T_HH[0] = 1.;
	thisTable->T_HH[1] = 2.;
	thisTable->T_HH[2] = 4.;
	thisTable->T_HH[3] = 6.;
	thisTable->T_HH[4] = 8.;
	thisTable->T_HH[5] = 10.;
	thisTable->T_HH[6] = 15.;
	thisTable->T_HH[7] = 20.;
	thisTable->T_HH[8] = 25.;
	thisTable->T_HH[9] = 30.;
	thisTable->T_HH[10] = 40.;
	thisTable->T_HH[11] = 50.;
	thisTable->T_HH[12] = 60.;
	thisTable->T_HH[13] = 70.;
	thisTable->T_HH[14] = 80.;
	thisTable->T_HH[15] = 90.;
	thisTable->T_HH[16] = 100.;
	thisTable->T_HH[17] = 200.;
	thisTable->T_HH[18] = 300.;

	thisTable->k10_HH[0] = 1.38e-13;
	thisTable->k10_HH[1] = 1.43e-13;
	thisTable->k10_HH[2] = 2.71e-13;
	thisTable->k10_HH[3] = 6.60e-13;
	thisTable->k10_HH[4] = 1.47e-12;
	thisTable->k10_HH[5] = 2.88e-12;
	thisTable->k10_HH[6] = 9.10e-12;
	thisTable->k10_HH[7] = 1.78e-11;
	thisTable->k10_HH[8] = 2.73e-11;
	thisTable->k10_HH[9] = 3.67e-11;
	thisTable->k10_HH[10] = 5.38e-11;
	thisTable->k10_HH[11] = 6.86e-11;
	thisTable->k10_HH[12] = 8.14e-11;
	thisTable->k10_HH[13] = 9.25e-11;
	thisTable->k10_HH[14] = 1.02e-10;
	thisTable->k10_HH[15] = 1.11e-10;
	thisTable->k10_HH[16] = 1.19e-10;
	thisTable->k10_HH[17] = 1.75e-10;
	thisTable->k10_HH[18] = 2.09e-10;
	
// 	for(int i=0; i<19; i++)
// 	{
// 		thisTable->T_HH[i] = log10(thisTable->T_HH[i]);
// 		thisTable->k10_HH[i] = log10(thisTable->k10_HH[i]);
// 	}
	
	thisTable->num = 17;
	thisTable->T = malloc(thisTable->num*sizeof(double));
	thisTable->k10_eH = malloc(thisTable->num*sizeof(double));
	thisTable->k10_pH = malloc(thisTable->num*sizeof(double));
	
	thisTable->T[0] = 1.;
	thisTable->T[1] = 2.;
	thisTable->T[2] = 5.;
	thisTable->T[3] = 10.;
	thisTable->T[4] = 20.;
	thisTable->T[5] = 50.;
	thisTable->T[6] = 100.;
	thisTable->T[7] = 200.;
	thisTable->T[8] = 500.;
	thisTable->T[9] = 1000.;
	thisTable->T[10] = 2000.;
	thisTable->T[11] = 3000.;
	thisTable->T[12] = 5000.;
	thisTable->T[13] = 7000.;
	thisTable->T[14] = 10000.;
	thisTable->T[15] = 15000.;
	thisTable->T[16] = 20000.;

	thisTable->k10_eH[0] = 0.2389e-9;
	thisTable->k10_eH[1] = 0.3371e-9;
	thisTable->k10_eH[2] = 0.5303e-9;
	thisTable->k10_eH[3] = 0.7459e-9;
	thisTable->k10_eH[4] = 1.047e-9;
	thisTable->k10_eH[5] = 1.629e-9;
	thisTable->k10_eH[6] = 2.260e-9;
	thisTable->k10_eH[7] = 3.106e-9;
	thisTable->k10_eH[8] = 4.595e-9;
	thisTable->k10_eH[9] = 5.917e-9;
	thisTable->k10_eH[10] = 7.153e-9;
	thisTable->k10_eH[11] = 7.712e-9;
	thisTable->k10_eH[12] = 8.170e-9;
	thisTable->k10_eH[13] = 8.321e-9;
	thisTable->k10_eH[14] = 8.366e-9;
	thisTable->k10_eH[15] = 8.285e-9;
	thisTable->k10_eH[16] = 8.114e-9;
	
	thisTable->k10_pH[0] = 0.40e-9;
	thisTable->k10_pH[1] = 0.45e-9;
	thisTable->k10_pH[2] = 0.43e-9;
	thisTable->k10_pH[3] = 0.369e-9;
	thisTable->k10_pH[4] = 0.317e-9;
	thisTable->k10_pH[5] = 0.3047e-9;
	thisTable->k10_pH[6] = 0.3379e-9;
	thisTable->k10_pH[7] = 0.4043e-9;
	thisTable->k10_pH[8] = 0.5471e-9;
	thisTable->k10_pH[9] = 0.7051e-9;
	thisTable->k10_pH[10] = 0.9167e-9;
	thisTable->k10_pH[11] = 1.070e-9;
	thisTable->k10_pH[12] = 1.301e-9;
	thisTable->k10_pH[13] = 1.48e-9;
	thisTable->k10_pH[14] = 1.695e-9;
	thisTable->k10_pH[15] = 1.975e-9;
	thisTable->k10_pH[16] = 2.201e-9;
	
// 	for(int i=0; i<17; i++)
// 	{
// 		thisTable->T[i] = log10(thisTable->T[i]);
// 		thisTable->k10_eH[i] = log10(thisTable->k10_eH[i]);
// 		thisTable->k10_pH[i] = log10(thisTable->k10_pH[i]);
// 	}
	
	return thisTable;
}

void deallocate_k10_table(k10_t *thisTable)
{
	free(thisTable->T_HH);
	free(thisTable->k10_HH);
	free(thisTable->T);
	free(thisTable->k10_eH);
	free(thisTable->k10_pH);
	free(thisTable);
}

double get_k10_HH(k10_t *thisTable, double temp)	//temperature in log10
{
	int i = 0;
	while(temp > thisTable->T_HH[i] && i<18)
	{
		i++;
	}
	
	if(i >= 19){
	      return 0.44*temp - 10.76;
	}else{
	      return thisTable->k10_HH[i-1] + (temp - thisTable->T_HH[i-1])*(thisTable->k10_HH[i] - thisTable->k10_HH[i-1])/(thisTable->T_HH[i] - thisTable->T_HH[i-1]);
	}
}

double get_k10_eH(k10_t *thisTable, double temp)
{
	int i = 0;
	while(temp > thisTable->T[i] && i<16)
	{
		i++;
	}
	if(i >= 17){
		return thisTable->k10_eH[16];
	}else{
		return thisTable->k10_eH[i-1] + (temp - thisTable->T[i-1])*(thisTable->k10_eH[i] - thisTable->k10_eH[i-1])/(thisTable->T[i] - thisTable->T[i-1]);
	}
}

double get_k10_pH(k10_t *thisTable, double temp)	//temperature in log10
{
	int i = 0;
	while(temp > thisTable->T[i] && i<16)
	{
		i++;
	}
	if(i >= 17){
		return 0.38*temp - 8.01;
	}else{
		return thisTable->k10_pH[i-1] + (temp - thisTable->T[i-1])*(thisTable->k10_pH[i] - thisTable->k10_pH[i-1])/(thisTable->T[i] - thisTable->T[i-1]);
	}
}