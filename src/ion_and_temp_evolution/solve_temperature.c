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

#include "../xray_heating_lya_coupling/phys_const.h"
#include "cell.h"
#include "solve_temperature.h"

double create_prod_rate_temp(cell_t *thisCell)
{
	double tmp = thisCell->heatHI*thisCell->XHI + thisCell->heatHeI*thisCell->XHeI + thisCell->heatHeII*thisCell->XHeII;
	
	return (2.*tmp)/(3.*boltzman_cgs);
}

double create_destr_rate_temp(double temperature, double n, double dn, double dt)
{
	const double ctn = 2.*temperature*dn/(3.*n*dt);
	
	return ctn;
}

double Patankar_step(double oldvalue, double intervalue, double dt, double prod_rate, double destr_rate){
	return intervalue*(oldvalue + dt*prod_rate)/(intervalue + dt*destr_rate);
}

double PatankarRKnoncon(double oldvalue, double dt, double dt_max, double n, double dn, cell_t *thisCell){
	double prod_rate, destr_rate, help_rate;
	double new_prod_rate, new_destr_rate;
	double intervalue, newvalue;
	
	prod_rate = create_prod_rate_temp(thisCell);
	destr_rate = create_destr_rate_temp(oldvalue, n, dn, dt_max);
	
	if(destr_rate < 0.){
		if(prod_rate >= 0.){
			prod_rate = prod_rate - destr_rate;
			destr_rate = 0.;
		}else{
			help_rate = - prod_rate;
			prod_rate = - destr_rate;
			destr_rate = help_rate;
		}
	}else{
		if(prod_rate < 0.){
			destr_rate = destr_rate - prod_rate;
			prod_rate = 0.;
		}
	}
	intervalue = Patankar_step(oldvalue, oldvalue, dt, prod_rate, destr_rate);

	new_prod_rate = prod_rate;
	new_destr_rate = create_destr_rate_temp(intervalue, n, dn, dt_max);
	
	if(new_destr_rate < 0.){
		if(new_prod_rate >= 0.){
			new_prod_rate = new_prod_rate - new_destr_rate;
			new_destr_rate = 0.;
		}else{
			help_rate = - new_prod_rate;
			new_prod_rate = - new_destr_rate;
			new_destr_rate = help_rate;
		}
	}else{
		if(new_prod_rate < 0.){
			new_destr_rate = new_destr_rate - new_prod_rate;
			new_prod_rate = 0.;
		}
	}
	newvalue = Patankar_step(oldvalue, intervalue, dt*0.5, prod_rate + new_prod_rate, destr_rate + new_destr_rate);
	return newvalue;
}



