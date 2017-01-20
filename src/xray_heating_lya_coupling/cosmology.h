#ifndef COSMOLOGY_H
#define COSMOLOGY_H
#endif

typedef struct
{
	double h;
	double omega_m;
	double omega_l;
	double omega_b;
	
	double z;
	double Hubble_z;
	
	double Y;
	double nH_z;
	double nHe_z;
    
    double f3He;
} cosmology_t;


cosmology_t *initCosmology();
void deallocate_cosmology(cosmology_t *thisCosmology);
cosmology_t *allocate_cosmology(double h, double omega_m, double omega_l, double omega_b, double z, double Y);
void cosmology_update_z(cosmology_t *thisCosmology, double z);


/*-------------------------------------------------------------------------------------*/
/* QUANTITIES AT REDSHIFT Z */
/*-------------------------------------------------------------------------------------*/

double calc_Hubble_z(double h, double omega_m, double omega_l, double z);
double calc_nH_z(double h, double omega_b, double z, double Y);		//number density per proton / H-atom
double calc_nHe_z(double h, double omega_b, double z, double Y);		//number density per He-atom

