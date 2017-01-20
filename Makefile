SOURCES := 	./src/main.c \
		./src/confObj.c \
		./src/evolution.c \
		./src/inputs.c \
		./src/outputs.c \
		./src/parse_ini.c \
		./src/utils.c \
		./src/xmem.c \
		./src/xstring.c \
		./src/ion_and_temp_evolution/cell.c \
		./src/ion_and_temp_evolution/evolution_loop.c \
		./src/ion_and_temp_evolution/recomb_rates.c \
		./src/ion_and_temp_evolution/solve_ionfraction.c \
		./src/ion_and_temp_evolution/solve_temperature.c \
		./src/ion_and_temp_evolution/solve_ion_temp.c \
		./src/ion_and_temp_evolution/fractions_heat_ion.c \
		./src/xray_heating_lya_coupling/21cm_grid.c \
		./src/xray_heating_lya_coupling/chem_const.c \
		./src/xray_heating_lya_coupling/collisional_coupling.c \
		./src/xray_heating_lya_coupling/compute_grid.c \
		./src/xray_heating_lya_coupling/convolution_fftw.c \
		./src/xray_heating_lya_coupling/cosmology.c \
		./src/xray_heating_lya_coupling/cross_sections.c \
		./src/xray_heating_lya_coupling/fftw_array_tools.c \
		./src/xray_heating_lya_coupling/phys_const.c \
		./src/xray_heating_lya_coupling/redshift_tools.c \
		./src/xray_heating_lya_coupling/spin_temperature.c \
		./src/xray_heating_lya_coupling/wouthuysen_effect.c \
		./src/xray_heating_lya_coupling/xray_heating.c \
		./src/xray_heating_lya_coupling/3cm_grid.c \
		./src/xray_heating_lya_coupling/spin_temperature_He.c
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := 21CMGRID

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
