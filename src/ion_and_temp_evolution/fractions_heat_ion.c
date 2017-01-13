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

#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))

/*--------------------------------------------------------------------------------------------*/
/* fraction of electron's energy going into secondary ionizations (Shull & van Steenberg 1985)*/

double fheat(double Xe)
{
    return 0.9971*(1.-pow(1.-pow(Xe,0.2663), 1.3163));
}

double fion_HI(double Xe)
{
    return 0.3908*pow(1.-pow(Xe, 0.4092), 1.7592);
}

double fion_HeI(double Xe)
{
    return 0.0554*pow(1.-pow(Xe, 0.4614), 1.6660);
}

double flya(double Xe)
{
    double const tmp = (1.-fheat(Xe)-fion_HI(Xe)-fion_HeI(Xe))*0.78;
    if(tmp>0. && tmp<1.)
    {
        return tmp;
    }else{
        return 0.;
    }
}

/*--------------------------------------------------------------------------------------------*/
