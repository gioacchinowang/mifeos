.SUFFIXES:
.SUFFIXES:  .out .o .f90 .for .f
# Edit this to be correct on your system
F90 = gfortran
F90OPTS = -O3 -fopenmp

HEALPIX = /home/lab/healpix

HEALPIXMOD = -I$(HEALPIX)/include
HEALPIXLIB = -L$(HEALPIX)/lib -lhealpix -lcfitsio

LAPACKLIB = -llapack -lblas

#--------------------------------------------------------------

all : CND_REG2D

#--------------------------------------------------------------
# Dependencies for object files
CND_REG2D.o: CND_REG2D_mod.o CND_REG2D.f90

#--------------------------------------------------------------
# test driver and executables linkage

CND_REG2D: CND_REG2D.o CND_REG2D_mod.o
	$(F90) $(F90OPTS) $(LDOPTS) -o CND_REG2D CND_REG2D.o CND_REG2D_mod.o \
		$(HEALPIXLIB) $(LAPACKLIB)

#-----------------------------------------------------

OBJ = CND_REG2D.o CND_REG2D_mod.o

objects: $(OBJ)
#----------------------------------------------------

CND_REG2D_mod.o : CND_REG2D_mod.f90
	$(F90) $(F90OPTS) -c CND_REG2D_mod.f90 -o CND_REG2D_mod.o $(HEALPIXMOD)

CND_REG2D.o : CND_REG2D.f90
	$(F90) $(F90OPTS) -c CND_REG2D.f90 -o CND_REG2D.o $(HEALPIXMOD)

#----------------------------------------------------------

clean:
	/bin/rm -f *.o *.mod CND_REG2D
