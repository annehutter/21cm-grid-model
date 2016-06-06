#ifndef SOLVE_IONFRACTION_H
#define SOLVE_IONFRACTION_H
#endif

void initializeVector(int dim_matrix, double vector[], double x[]);
void checkVector(int dim_matrix, double vector[]);
void create_rate_functions(int dim_matrix, double d[][dim_matrix], double vector[], recomb_t *thisRecombRates, cell_t *thisCell);
void add_rate_functions(int dim_matrix, double d[][dim_matrix], double d1[][dim_matrix]);
void calculateMatrix(int dim_matrix, double delta_t, double matrix[][dim_matrix], double d[][dim_matrix], double vector[]);
void solveSetOfEquations(int dim_matrix, double matrix[][dim_matrix], double vector[], double x[]);
void ModifiedPatankarRK(int dim_matrix, double oldvalue[], double dt, double newvalue[], recomb_t *thisRecombRates, cell_t *thisCell);
