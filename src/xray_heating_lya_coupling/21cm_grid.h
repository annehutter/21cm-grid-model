#ifndef L21CM_GRID_H
#define L21CM_GRID_H
#endif

typedef struct
{
	int nbins;
	double box_size;
	
	double z_old;
	fftw_complex *temp_old;
	fftw_complex *dens_old;
	fftw_complex *Xe_old;
	
	double z;
	fftw_complex *temp;
	fftw_complex *dens;
	fftw_complex *Xe;
	
	fftw_complex *XHI;
	fftw_complex *XHeI;
	fftw_complex *XHeII;
	fftw_complex *Ts_inv;
	fftw_complex *Tb;
	
	int local_n0;
	int local_0_start;
} grid_21cm_t;

/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR 21CM GRID */
/*-------------------------------------------------------------------------------------*/

grid_21cm_t *init21cmgrid();
void deallocate_21cmgrid(grid_21cm_t *this21cmGrid);
grid_21cm_t *allocate_21cmgrid(int nbins, float box_size);

void initialize_21cmgrid(fftw_complex *thisArray, int nbins, int local_n0, double value);

void set_temperature_21cmgrid(grid_21cm_t *this21cmGrid, double value);
void read_density_21cmgrid(grid_21cm_t *this21cmGrid, char *filename, int double_precision);

double get_mean_Xe_21cmgrid(grid_21cm_t *this21cmGrid);
double get_mean_temp_21cmgrid(grid_21cm_t *this21cmGrid);

void write_Tb_field_file(grid_21cm_t *this21cmGrid, char *filename);
void write_Ts_field_file(grid_21cm_t *this21cmGrid, char *filename);
void write_Xe_field_file(grid_21cm_t *this21cmGrid, char *filename);
void write_temp_field_file(grid_21cm_t *this21cmGrid, char *filename);
