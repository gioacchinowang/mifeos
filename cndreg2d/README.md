## INTRODUCTION

3 Minkowski functionals + number of isolated connected clusters

This is a scientific code which has been developed to work in our environment.
Please let us know if you have any problems running it for your task and,
of course, if you find any bugs. 

If used in publications, please cite:

1) [Ducout, A., Bouchet, F., Colombi, S., Pogosyan, D., Prunet, S., MNRAS, 429, 3, 2104-2126](https://academic.oup.com/mnras/article/429/3/2104/1001307) 
2) [Gay, C., Pichon, C., Pogosyan, D., PhyRevD, 85, 023011](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.85.023011)

## INSTALLATION

Code is in fortran90/95 and has been checked with Intel ifort compiler.
Reguirements: HEALPix, CFITSIO, LAPACK and BLAS available.

STRUCTURE of the code is the following:

1) CND_REG2D.f90 contains main program + subroutine mink

main program : 
    constructs/reads data (map+mask) calls subroutine mink
mink : 
    loop on thresholds, writes results to a text file,
    calls subroutine CND_REG2D, in module CND_REG2D_mod.f90

This part of the code is expected to be modified by a user to suit
his/her particular needs.

2) CND_REG2D_mod.f90 : 
    the core algorithm for MFs, computes the 3 functionals at a given threshold.
    
## POSSIBLE ISSUES

for very high resolution data, some 4-byte INTEGER(I4B) variables may become
insufficient, and will need to be replaced by 8-byte INTEGER(I8B)


## AUTHOR

- Dr. Anne Ducout
- Dr. Dmitri Pogosyan
