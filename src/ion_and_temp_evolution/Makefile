SOURCES := 	./main.c \
		./cell.c \
		./evolution_loop.c \
		./recomb_rates.c \
		./solve_ionfraction.c \
		./solve_temperature.c 
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := evoltemp

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
