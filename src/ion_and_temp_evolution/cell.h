#ifndef CELL_H
#define CELL_H
#endif

typedef struct
{
	double XHI;
	double XHII;
	double XHeI;
	double XHeII;
	double XHeIII;
	
	double temp;
	
	double heatHI;
	double heatHeI;
	double heatHeII;
	
	double photIonHI;
	double photIonHeI;
	double photIonHeII;
	
	double dens;
} cell_t;

cell_t *initCell();
void deallocate_cell(cell_t *thisCell);
cell_t *initCell_temp(double temp, double dens);
void update_X_cell(cell_t *thisCell, double XHII, double XHeII, double XHeIII);
void update_photIon_cell(cell_t *thisCell, double photIonHI, double photIonHeI, double photIonHeII);
void update_dens_cell(cell_t *thisCell, double dens);
void update_temp_cell(cell_t *thisCell, double temp);
void update_heat_cell(cell_t *thisCell, double heatHI, double heatHeI, double heatHeII);