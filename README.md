# 21cm_grid
Code to compute HI 21cm & 3HeII 3.5cm brightness temperature

Requirements:
- Serial run: fftw3 library, gsl library (math, integrate)
- Parallel run with MPI: fftw3, fftw3-mpi and mpi libraries, gsl library (math, integrate)

Compilation:
make

Usage:
./21CMGRID iniFile.ini

iniFile.ini contains the necessary input files and values for a run. 


Parameters in iniFile.ini

General:
- inputFilesAreInDoublePrecision: 0 for single, 1 for double precision of data files to be read in

- inputList: list with redshift steps: z, dens_id, xray_id, lumX, alphaX, nuX_min, lya_id, lumLya, alphaLya, nuLya_min
- inputIgmDensityFile: path to density grids (without "_<dens_id>" at the end)
- inputXraySourcesFile: path to Xray sources grids (without "_<xray_id>" at the end)
- inputLyaSourcesFile: path to Lya sources grids (without "_<lya_id>" at the end)

- gridsize: size of the grid (should be a power of 2)
- boxsize: comoving boxsize in Mpc/h

- initialRedshift: redshift at which simulation starts
- finalRedshift: redshift at which simulation ends
- dz_max: maximum decrease in redshift when evolving

- h: H = 100*h km/s/Mpc
- omega_b: baryon density parameter
- omega_m: matter density parameter
- omega_l: lambda density parameter
- sigma8: sigma8
- Y: mass fraction of Helium in the primordial gas (assumed to consist of H and He)

- solve_3He: 1 for solving brightness temperature for 3HeII, 0 for only HI

- initialIGMtemperature: initial IGM temperature at "initialRedshift"
