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

#include "phys_const.h"
#include "chem_const.h"

double cross_sec_HI(double nu)
{
	if(nu >= nu_HI){
		return 6.3e-18*pow(nu/nu_HI,-3);
	}else{
		return 0.;
	}
}

double cross_sec_HeI(double nu)
{
	const double tmp = nu/nu_HeI;
	
	if(nu >= nu_HeI){
		return 7.2e-18*(1.66*pow(tmp,-2.05) + 0.66*pow(tmp,-3.05));
	}else{
		return 0.;
	}
}

double cross_sec_HeII(double nu)
{
	if(nu >= nu_HeII){
		return 1.58e-18*pow(nu/nu_HeII,-3);
	}else{
		return 0.;
	}
}