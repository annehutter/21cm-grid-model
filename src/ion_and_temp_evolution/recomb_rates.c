#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "recomb_rates.h"

double recHII(double temp){
	return 6.28e-11/sqrt(temp)*pow(temp/1.e3,-0.2)/(1.+pow(temp/1.e6,0.7));	//case A
// 	return 2.6e-13*pow(temp/1.e4,-0.85);	//case B
}

double recHeII(double temp){
	return 1.5e-10*pow(temp,-0.6353);
}

double recHeIII(double temp){
	return 3.36e-10/sqrt(temp)*pow(temp/1.e3,-0.2)/(1.+pow(temp/4.e6,0.7));
}

recomb_t *calcRecRate(){
	double temp, temp_max;
	double saturation;
	int lower_steps, upper_steps;
	int N;
	
	recomb_t *thisRecombRates;
	
	thisRecombRates = malloc(sizeof(recomb_t));
	if(thisRecombRates == NULL)
	{
		fprintf(stderr, "thisRecombRates in calcRecRate (ion_recomb_rates.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	temp_max = 1.e5;
	thisRecombRates->temp_max = temp_max;
	
	saturation = 500.;
	lower_steps = 2000;
	upper_steps = 1990;
	thisRecombRates->saturation = saturation;
	thisRecombRates->lower_steps = lower_steps;
	thisRecombRates->upper_steps = upper_steps;
	
	N = thisRecombRates->lower_steps + thisRecombRates->upper_steps;
	thisRecombRates->N = N;
	
	thisRecombRates->recHII = malloc(sizeof(double)*(N+1));
	thisRecombRates->recHeII = malloc(sizeof(double)*(N+1));
	thisRecombRates->recHeIII = malloc(sizeof(double)*(N+1));

	for(int i=0; i<=lower_steps; i++){
		if(i==0){
			temp = saturation/(double)lower_steps;
		}
		else{
			temp = saturation/(double)lower_steps*i;
		}
		thisRecombRates->recHII[i] = recHII(temp);
		thisRecombRates->recHeII[i] = recHeII(temp);
		thisRecombRates->recHeIII[i] = recHeIII(temp);
//         printf("%d: \t%e\t%e\n", i, temp, thisRecombRates->recHII[i]);
	}
	for(int i=1; i<=upper_steps; i++){
		temp = saturation + (temp_max - saturation)/upper_steps*i;
		thisRecombRates->recHII[i + lower_steps] = recHII(temp);
		thisRecombRates->recHeII[i + lower_steps] = recHeII(temp);
		thisRecombRates->recHeIII[i + lower_steps] = recHeIII(temp);
//         printf("%d: \t%e\t%e\n", i+lower_steps, temp, thisRecombRates->recHII[i + lower_steps]);
	}
	
	return thisRecombRates;
}

void deallocate_recomb(recomb_t *thisRecombRates)
{
	free(thisRecombRates->recHII);
	free(thisRecombRates->recHeII);
	free(thisRecombRates->recHeIII);
	free(thisRecombRates);
}

void extend_recomb(recomb_t *thisRecombRates, double new_temp_max)
{
	double temp;
	double temp_max_old = thisRecombRates->temp_max;
	double saturation = thisRecombRates->saturation;
	int upper_steps_old = thisRecombRates->upper_steps;
	int upper_steps;
	int lower_steps = thisRecombRates->lower_steps;
	int N;
	
	thisRecombRates->upper_steps += (new_temp_max-temp_max_old)/((temp_max_old-saturation)/upper_steps_old);
	upper_steps = thisRecombRates->upper_steps;
	
	thisRecombRates->temp_max = new_temp_max;
	N = thisRecombRates->lower_steps + thisRecombRates->upper_steps;
	thisRecombRates->N = N;
	
	thisRecombRates->recHII = realloc(thisRecombRates->recHII,(N+1));
	thisRecombRates->recHeII = realloc(thisRecombRates->recHeII,(N+1));
	thisRecombRates->recHeIII = realloc(thisRecombRates->recHeIII,(N+1));
	
	for(int i=upper_steps; i<=thisRecombRates->upper_steps; i++)
	{
		temp = saturation + (temp_max_old - saturation)/upper_steps_old*i;
		thisRecombRates->recHII[i + lower_steps] = recHII(temp);
		thisRecombRates->recHeII[i + lower_steps] = recHeII(temp);
		thisRecombRates->recHeIII[i + lower_steps] = recHeIII(temp);
	}
}

int temp_recomb_index(recomb_t *thisRecombRates, double temperature)
{
	int index = 0;
	if(temperature <= thisRecombRates->saturation) {
		index = (int)(temperature*thisRecombRates->lower_steps)/thisRecombRates->saturation;
	}
	else if(temperature > thisRecombRates->saturation) {
		index = (int)(thisRecombRates->lower_steps + (temperature - thisRecombRates->saturation)*thisRecombRates->upper_steps/ (thisRecombRates->temp_max - thisRecombRates->saturation));
		if(temperature > thisRecombRates->temp_max) {
			extend_recomb(thisRecombRates, temperature);
			index = thisRecombRates->N;
		}
	}
	return index;
}
