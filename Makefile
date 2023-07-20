.PHONY: all clean setup

F90 = gfortran
F77 = gfortran -ffixed-form
RM = trash
LIBGMXFLAGS = -I/usr/local/include `pkg-config --libs libgmxfort`
F90FLAGS = -Wall -fbackslash
F77FLAGS = 
PROGRAMS = bin/diffusion bin/linreg bin/rdf

all: setup $(PROGRAMS)

bin/%: scripts/%.f90 build/spherenoncube.o
	$(F90) $^ -o $@ $(F90FLAGS) $(LIBGMXFLAGS)

build/spherenoncube.o: scripts/spherenoncube.f
	$(F77) -c $^ $(F77FLAGS) -o $@

setup:
	mkdir -p bin build

clean:
	$(RM) -rf bin/* build/*