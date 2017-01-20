#ifndef L3CM_GRID_H
#define L3CM_GRID_H
#endif

typedef struct
{
	int nbins;
	double box_size;
	
	double z_old;
	double temp_old;
	
	double z;
	double temp;
	
	fftw_complex *Ts_inv;
	fftw_complex *Tb;
	
	int local_n0;
	int local_0_start;
} grid_3cm_t;

/*-------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR 21CM GRID */
/*-------------------------------------------------------------------------------------*/

grid_3cm_t *init3cmgrid();
void deallocate_3cmgrid(grid_3cm_t *this3cmGrid);
grid_3cm_t *allocate_3cmgrid(int nbins, float box_size);

void initialize_3cmgrid(fftw_complex *thisArray, int nbins, int local_n0, double value);
void set_temperature_3cmgrid(grid_3cm_t *this3cmGrid, double value);

void write_Tb_3cm_field_file(grid_3cm_t *this3cmGrid, char *filename);
void write_Ts_3cm_field_file(grid_3cm_t *this3cmGrid, char *filename);
