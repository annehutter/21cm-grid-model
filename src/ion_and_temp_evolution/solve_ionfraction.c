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

#include "cell.h"
#include "recomb_rates.h"
#include "solve_ionfraction.h"

void initializeVector(int dim_matrix, double vector[], double x[]){
	for(int i=0; i<dim_matrix; i++){
		vector[i] = x[i];
	}
}

void checkVector(int dim_matrix, double vector[]){
	for(int i=0; i<dim_matrix;i++){
		if(vector[i]==0.) vector[i] = 1.e-30;
		else if(vector[i]<1.e-30) vector[i] = 1.e-30;
		else if(vector[i] >= 1.e-30) continue;
		else {
			for(int j=0; j<dim_matrix; j++) printf("vector[%d] = %e\t", j, vector[j]);
			printf("\n");
			fprintf(stderr, "HI/HII or HI/HII & HeI/HeII/HeIII have not allowed values\n");
			exit(EXIT_FAILURE);
		}
	}
}

void create_rate_functions(int dim_matrix, double d[][dim_matrix], double vector[], recomb_t *thisRecombRates, cell_t *thisCell){

	/* initialization */
	for (int i=0; i<dim_matrix; i++){
		for(int j=0; j<dim_matrix; j++){
			d[i][j]=0.;
		}
	}

	double ne = (vector[1] + vector[3] + 2.*vector[4])*thisCell->dens*1.9e-7;
	double temperature = thisCell->temp;
	int temp_rec_index = temp_recomb_index(thisRecombRates, temperature);
	
// 	printf("temperature = %e\t %d\n", temperature, temp_rec_index);
// 	printf("recombHII = %e\t recombHeII = %e\t recombHeIII = %e\n", thisRecombRates->recHII[temp_rec_index], thisRecombRates->recHeII[temp_rec_index], thisRecombRates->recHeIII[temp_rec_index]);
// 	printf("v[1] = %e\t v[3] = %e\t v[4] = %e\n", vector[1], vector[3], vector[4]);
// 	printf("dens = %e\n", thisCell->dens);
// 	printf("ne = %e\n",ne);
	if(dim_matrix == 5) {
		d[0][1] = thisCell->photIonHI*vector[0];
		d[1][0] = thisRecombRates->recHII[temp_rec_index]*vector[1]*ne;
		d[2][3] = thisCell->photIonHeI*vector[2];
		d[3][2] = thisRecombRates->recHeII[temp_rec_index]*vector[3]*ne;
		d[3][4] = thisCell->photIonHeII*vector[3];
		d[4][3] = thisRecombRates->recHeIII[temp_rec_index]*vector[4]*ne;
	}

	if(dim_matrix == 2) {
		d[0][1] = thisCell->photIonHI*vector[0];
		d[1][0] = thisRecombRates->recHII[temp_rec_index]*vector[1]*ne;
	}
	
// 	for (int i=0; i<dim_matrix; i++){
// 		for(int j=0; j<dim_matrix; j++){
// // 			printf("desctr[%d][%d] = %e\n", i, j, d[i][j]);
// 		}
// 	}
}

void add_rate_functions(int dim_matrix, double d[][dim_matrix], double d1[][dim_matrix]){
	for(int i=0; i < dim_matrix; i++){
		for(int j=0; j<dim_matrix; j++){
			d[i][j]=d1[i][j]+d[i][j];
		}
	}
}

void calculateMatrix(int dim_matrix, double delta_t, double matrix[][dim_matrix], double d[][dim_matrix], double vector[]){
	double value;

	for(int i=0; i<dim_matrix; i++){
        
        double tmp1 = delta_t/(vector[i]+1.e-30);
        
		for(int j=0; j<dim_matrix; j++){
			if(i==j) {
				value=0.;
				for(int k=0; k<dim_matrix; k++){
					value=value+d[i][k];
				}
				matrix[j][i] = 1. + (value-d[j][i])*tmp1;
			}
			else {
				matrix[j][i] = -d[i][j]*tmp1;
			}
// 			printf("matrix[%d][%d] = %e\n", i,j,matrix[i][j]);
		}
	}
}

void solveSetOfEquations(int dim_matrix, double matrix[][dim_matrix], double vector[], double x[]){
	if(dim_matrix ==2 || dim_matrix ==5){
		if(matrix[0][0]*matrix[1][1] != matrix[0][1]*matrix[1][0]){
			x[1] = (matrix[0][0]*vector[1]-matrix[1][0]*vector[0])/(matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]);
			if(matrix[0][0]!=0) {
				x[0] = (vector[0]-matrix[0][1]*x[1])/matrix[0][0];
			}else{
				x[0] = (-matrix[0][1]*vector[1]+matrix[1][1]*vector[0])/(matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]);
			}
		}else{
			fprintf(stderr, "ERROR in solving chemical equations!\n");
			exit(EXIT_FAILURE);
		}
	}
	if(dim_matrix ==5){
		double temp1 = (matrix[4][2]*(matrix[2][4]*matrix[3][3]-matrix[2][3]*matrix[3][4])-matrix[4][3]*(matrix[2][4]*matrix[3][2]-matrix[2][2]*matrix[3][4])+matrix[4][4]*(matrix[2][3]*matrix[3][2]-matrix[2][2]*matrix[3][3]));
		double temp2 = 1./temp1;
		if(temp1 != 0){
			x[2] = ((matrix[3][4]*matrix[4][3]-matrix[3][3]*matrix[4][4])*vector[2] - (matrix[2][4]*matrix[4][3]-matrix[2][3]*matrix[4][4])*vector[3] + (matrix[2][4]*matrix[3][3]-matrix[2][3]*matrix[3][4])*vector[4])*temp2;
			x[3] = -((matrix[3][4]*matrix[4][2]-matrix[3][2]*matrix[4][4])*vector[2] - (matrix[2][4]*matrix[4][2]-matrix[2][2]*matrix[4][4])*vector[3] + (matrix[2][4]*matrix[3][2]-matrix[2][2]*matrix[3][4])*vector[4])*temp2;
			x[4] = ((matrix[3][3]*matrix[4][1]-matrix[3][2]*matrix[4][3])*vector[2] - (matrix[2][3]*matrix[4][2]-matrix[2][2]*matrix[4][3])*vector[3] + (matrix[2][3]*matrix[3][2]-matrix[2][2]*matrix[3][3])*vector[4])*temp2;
		}else{
			fprintf(stderr, "ERROR in solving chemical equations!\n");
			exit(EXIT_FAILURE);
		}
	}
}

void ModifiedPatankarRK(int dim_matrix, double oldvalue[], double dt, double newvalue[], recomb_t *thisRecombRates, cell_t *thisCell){
	double intervalue[dim_matrix];
	double oldvalue_local[dim_matrix];
	double newvalue_local[dim_matrix];
	double destr_rate[dim_matrix][dim_matrix];
	double inter_destr_rate[dim_matrix][dim_matrix];
	double matrix[dim_matrix][dim_matrix];

	initializeVector(dim_matrix, oldvalue_local, oldvalue);
	initializeVector(dim_matrix, newvalue_local, oldvalue);
	
// 	printf("oldvalue_local[0] = %e\toldvalue_local[1] = %e\n", oldvalue_local[0], oldvalue_local[1]);

// 	First Runge Kutta step
	initializeVector(dim_matrix, oldvalue_local, newvalue_local);	//oldvalue = x

	create_rate_functions(dim_matrix, destr_rate, oldvalue_local, thisRecombRates, thisCell);	//calculate destruction rates

	calculateMatrix(dim_matrix, dt, matrix, destr_rate, oldvalue_local);	//calculate matrix for destr_rate and oldvalue

	solveSetOfEquations(dim_matrix, matrix, oldvalue_local, intervalue);	//Solve set of equations for intervalue

	checkVector(dim_matrix,intervalue);

// 	Second Runge Kutta step
	initializeVector(dim_matrix, oldvalue_local, oldvalue);		//oldvalue = newvalue

	create_rate_functions(dim_matrix, inter_destr_rate, intervalue, thisRecombRates, thisCell);	//calculate new destruction rates from intervalue

	add_rate_functions(dim_matrix, inter_destr_rate, destr_rate);	//add rate functions destr_rate and interdestr_rate, output on interdestr_rate

	calculateMatrix(dim_matrix, dt*0.5, matrix, inter_destr_rate, intervalue);	//calculate matrix for interdestr_rate and intervalue

	solveSetOfEquations(dim_matrix, matrix, oldvalue_local, newvalue_local);	//solve set of equations for newvalue

	checkVector(dim_matrix, newvalue_local);

	initializeVector(dim_matrix, newvalue, newvalue_local);
}