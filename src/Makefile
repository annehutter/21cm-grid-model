SOURCES := 	./main.c \
		./evolution.c \
		./ion_and_temp_evolution/cell.c \
		./ion_and_temp_evolution/evolution_loop.c \
		./ion_and_temp_evolution/recomb_rates.c \
		./ion_and_temp_evolution/solve_ionfraction.c \
		./ion_and_temp_evolution/solve_temperature.c \
		./xray_heating_lya_coupling/21cm_grid.c \
		./xray_heating_lya_coupling/chem_const.c \
		./xray_heating_lya_coupling/collisional_coupling.c \
		./xray_heating_lya_coupling/compute_grid.c \
		./xray_heating_lya_coupling/convolution_fftw.c \
		./xray_heating_lya_coupling/cosmology.c \
		./xray_heating_lya_coupling/cross_sections.c \
		./xray_heating_lya_coupling/fftw_array_tools.c \
		./xray_heating_lya_coupling/phys_const.c \
		./xray_heating_lya_coupling/redshift_tools.c \
		./xray_heating_lya_coupling/spin_temperature.c \
		./xray_heating_lya_coupling/wouthuysen_effect.c \
		./xray_heating_lya_coupling/xray_heating.c 
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := 21CMGRID

USE-MPI = YES

OPTIMIZE = -O3 -ftree-vectorize
WARNING = -Wall -Wextra -Wshadow -g

ifdef USE-MPI
	CC := mpicc
	CFLAGS := -c -std=c99 -march=native -lm $(WARNING) $(OPTIMIZE) -D __MPI
	LDFLAGS := -lfftw3_mpi -lfftw3 -lm -lmpich -lgsl -lgslcblas

else
	CC := gcc
	CFLAGS := -c -std=c99 -march=native -lm $(WARNING) $(OPTIMIZE) 
	LDFLAGS := -lfftw3 -lm -lgsl -lgslcblas
endif

.PHONY: all clean

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
