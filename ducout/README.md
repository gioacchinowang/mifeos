Code : 3 Minkowski functionals + number of isolated connected clusters

This is a scientific code which has been developed to work in our environment.
Please let us know if you have any problems running it for your task and,
of course, if you find any bugs. 

If used in publications, please cite:

1) arXiv:1209.1223 (or updated journal reference as it comes out)
2) Gay, C., Pichon, C., Pogosyan, D. Physical Review D, vol. 85, Issue 2, id. 023011.

CONTACT

Anne Ducout (anne at ducout.com),  Dmitri Pogosyan (pogosyan at ualberta.ca)

INSTALLATION

Code is in fortran90/95 and has been checked with Intel ifort compiler.
Reguirements: HEALPix, CFITSIO and LAPACK available

supplied Makefile should compile the code if ${HEALPIX} environment
variable is set,  cfistio is a system library (i.e can be found in standard
library locations) and one uses LAPACK/BLAS that comes with 
Intel Composer XE suite.

As an example of a less standard installation, Makefile_macosx_sample is
provided.  This makefile will not work out-of-the-box since it contains
references to private directories, but may serve as a guidance for setup
under MacOSX.


STRUCTURE of the code is the following:

1) CND_REG2D.f90 contains main program + subroutine mink

main program : constructs/reads data (map+mask) calls subroutine mink
mink         : loop on thresholds, writes results to a text file,
               calls subroutine CND_REG2D, in module CND_REG2D_mod.f90

This part of the code is expected to be modified by a user to suit
his/her particular needs.


2) CND_REG2D_mod.f90 : the core algorithm for MFs, computes the 4 functionals
                       at a given threshold.

The package also contains  a sample nside=64 map "map_cmb_ns64.fits"
and the mask "mask_gal_fsky0.80_ns64.fits" as a test case
as well as the expected output of the code for this data in "mf.dat".


GIT REPOSITORY

For git users, the code (and perhaps bugfixes) is available
from the git repository, try

git clone git://quad-opteron.nic.ualberta.ca/genus.git


POSSIBLE ISSUES

for very high resolution data, some 4-byte INTEGER(I4B) variables may become
insufficient, and will need to be replaced by 8-byte INTEGER(I8B)
