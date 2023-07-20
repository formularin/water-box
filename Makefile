.PHONY: all clean setup

F90 = gfortran
RM = trash
LIBGMXFLAGS = -I/usr/local/include `pkg-config --libs libgmxfort`
F90FLAGS = -Wall -fbackslash
PROGRAMS = bin/diffusion bin/linreg

all: setup $(PROGRAMS)

bin/%: scripts/%.f90
	$(F90) $^ -o $@ $(F90FLAGS) $(LIBGMXFLAGS)

setup:
	mkdir -p bin

clean:
	$(RM) -rf bin/*