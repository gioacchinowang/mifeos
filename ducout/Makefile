.SUFFIXES:
.SUFFIXES:  .out .o .f90 .for .f
# Edit this to be correct on your system
F90 = ifort
MODULEDIR= include
#F90OPTS = -check bounds -check format -warn declarations -warn argument_checking
#F90OPTS = -O4 -static-libcxa -module $(MODULEDIR) -I$(MODULEDIR) -qp
F90OPTS = -O3 -module $(MODULEDIR) -I$(MODULEDIR) -openmp
LDOPTS = 

# We rely on HEALPIX variable being set by standard HealPix config script. 
# If not, set it manually
# HEALPIX = /PATH_TO_INSTALLATION_DIR/Healpix_2.12a

HEALPIXMOD = -I$(HEALPIX)/include
HEALPIXLIB = -L$(HEALPIX)/lib -lhealpix -lcfitsio

# Now we use LAPACK for Intel mkl.  For my machine MKLROOT is
# set by startup scripts

MKLARCH = em64t
MKLPATH=$(MKLROOT)/lib/$(MKLARCH)
                          # Other choices:  32 64
			  #

LAPACK  = -I$(MKLROOT)/include -mkl 

LAPACKLIB = $(LAPACK) -lguide -lpthread

#--------------------------------------------------------------

all : CND_REG2D

#--------------------------------------------------------------
# Dependencies for object files
CND_REG2D.o:		CND_REG2D_mod.o CND_REG2D.f90

#--------------------------------------------------------------
# test driver and executables linkage

CND_REG2D: CND_REG2D.o CND_REG2D_mod.o
	$(F90) $(F90OPTS) $(LDOPTS) -o CND_REG2D CND_REG2D.o CND_REG2D_mod.o \
		$(HEALPIXLIB) $(LAPACKLIB)

#-----------------------------------------------------

OBJ = CND_REG2D.o CND_REG2D_mod.o

objects: $(OBJ)
	echo "Object files made"
#----------------------------------------------------

CND_REG2D_mod.o : CND_REG2D_mod.f90
	$(F90) $(F90OPTS) -c CND_REG2D_mod.f90 $(HEALPIXMOD)
CND_REG2D.o : CND_REG2D.f90
	$(F90) $(F90OPTS) -c CND_REG2D.f90 $(HEALPIXMOD)

#----------------------------------------------------------

clean:
	/bin/rm -f *.o CND_REG2D 

cleanall:
	/bin/rm -f *.o CND_REG2D $(MODULEDIR)/[!n]*.mod	
