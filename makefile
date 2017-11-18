



FC = gfortran

UTILS_DIR = $(HOME)/fortran-utils/src/


FCFLAGS = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -I$(UTILS_DIR) -ffpe-trap=zero,underflow,overflow
FLDFLAGS =  -llapack -L$(UTILS_DIR) -lfortran_utils -lblas


FSOURCES =    sudo.f90


all: $(FSOURCES)
	$(FC)  $(FCFLAGS)  $? -o sudo.x $(FLDFLAGS)
test1: tests/test_interpolation.f90
	$(FC)  $(FCFLAGS)  $? -o test1.x $(FLDFLAGS)

 %.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<
