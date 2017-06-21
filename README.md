Description
===========

Code to compute HI 21cm & 3HeII 3.5cm brightness temperature

When you should use this code
=============================

If you want to compute the 21cm signal in the early stages of reionization from a cosmological simulation. You will need

- cosmological box with DM/gas overdensities (grid)
- a box with ovderdensity of number of ionizing photons (grid)

Why should you use it
=====================

1. **MPI Parallel** The code can be run on multiple cores and distributed memory.
2. **redshift-variable spectral index** An input list needs to be provided where the spectral index of the Lyman-alpha and x-ray sources can be changed over time.
3. **3HeII 3cm emission** The code includes the option to compute the 3HeII 3cm emission line.

Installation
============

Pre-requisities
---------------

Serial run
``````````

1. fftw3 library: ``fftw3 >= 3.3.3``
2. gsl library (math, integrate): ``gsl >= 1.16``

Parallel run
````````````

1. MPI library
2. fftw3 & fftw3-mpi library: ``fftw3 >= 3.3.3``
3. gsl library (math, integrate): ``gsl >= 1.16``

FFTW3
'''''

Go to the `FFTW webpage <http://www.fftw.org/download.html>`__ to install fftw3. Ensure to compile the library with the ``enable-mpi`` flag for parallel runs
::
    
    $ ./configure --enable-mpi
    $ make
    $ make install
    
GSL
'''

Go to the `GSL webpage <https://www.gnu.org/software/gsl/>`__ to install gsl and follow the instructions there. 


Download & Build
----------------

::

    $ git clone https://github.com/annehutter/21cm_grid.git
    $ make

This will download the code and first test case from the github directory and compile the source code.

Execution
---------

The code is executed by
::

    $ ./21CMGRID iniFile.ini

``iniFile.ini`` contains all input parameters that are needed for any runs. For a different simulation the code does not need to be recompiled but just this parameter file iniFile.ini to be adapted.


Parameter file
''''''''''''''

**General**
...........

- ``inputFilesAreInDoublePrecision``: 0 for single, 1 for double precision of data files to be read in

- ``inputList``: list with redshift steps: z, dens_id, xray_id, lumX, alphaX, nuX_min, lya_id, lumLya, alphaLya, nuLya_min
- ``inputIgmDensityFile``: path to density grids (without "_<dens_id>" at the end)
- ``inputXraySourcesFile``: path to Xray sources grids (without "_<xray_id>" at the end)
- ``inputLyaSourcesFile``: path to Lya sources grids (without "_<lya_id>" at the end)

- ``gridsize``: size of the grid (should be a power of 2)
- ``boxsize``: comoving boxsize in Mpc/h

- ``initialRedshift``: redshift at which simulation starts
- ``finalRedshift``: redshift at which simulation ends
- ``dz_max``: maximum decrease in redshift when evolving

- ``h``: H = 100*h km/s/Mpc
- ``omega_b``: baryon density parameter
- ``omega_m``: matter density parameter
- ``omega_l``: lambda density parameter
- ``sigma8``: sigma8
- ``Y``: mass fraction of Helium in the primordial gas (assumed to consist of H and He)

- ``solve_3He``: 1 for solving brightness temperature for 3HeII, 0 for only HI

- ``initialIGMtemperature``: initial IGM temperature at "initialRedshift"
