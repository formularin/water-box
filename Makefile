.PHONY: all clean setup

F90 = gfortran
F77 = gfortran -ffixed-form
RM = trash
LIBGMXFLAGS = -I/usr/local/include `pkg-config --libs libgmxfort`
F90FLAGS = -Wall -fbackslash
F77FLAGS = 
PROGRAMS = bin/diffusion bin/rdf bin/count_waters

all: setup $(PROGRAMS)

bin/%: build/utils.o build/spherenoncube.o scripts/%.f90
	$(F90) $^ -o $@ $(F90FLAGS) $(LIBGMXFLAGS)

build/spherenoncube.o: scripts/spherenoncube.f
	$(F77) -c $^ $(F77FLAGS) -o $@

build/utils.o: scripts/utils.f90
	$(F90) -c $^ $(F90FLAGS) -o $@

setup:
	mkdir -p bin build

clean:
	$(RM) -rf bin/* build/*