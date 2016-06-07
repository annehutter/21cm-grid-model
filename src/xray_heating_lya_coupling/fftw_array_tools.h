#ifndef FFTW_ARRAY_TOOLS_H
#define FFTW_ARRAY_TOOLS_H
#endif

/*-------------------------------------------------------------------------------------*/
/* Reading files to fftw arrays */
/*-------------------------------------------------------------------------------------*/


#ifdef __MPI
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename);
#else
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, char *filename);
#endif


#ifdef __MPI
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename);
#else
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, char *filename);
#endif

/*-------------------------------------------------------------------------------------*/
/* Writing fftw arrays to file */
/*-------------------------------------------------------------------------------------*/

#ifdef __MPI
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename);
#else
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, char *filename);
#endif


#ifdef __MPI
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename);
#else
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, char *filename);
#endif
