# **modern\_fishpack**
 
This Fortran project is a modernization of NCAR's FISHPACK90 library. The original work, written in both FORTRAN 77 and Fortran 90, was heavily refactored to incorporate features of modern Fortran. That is, all the procedures are now encapsulated in modules. Most importantly, the potential memory leak in the derived type **fish**, now rebaptized as **FishpackWorkspace**, was addressed by replacing pointers with allocatable arrays. 

This project is still a work in progress.

-----------------------------------------------------------------------------

# TODO
* Introduce interfaces to replace assumed-size arrays with assumed shape arrays. 
* Implement object-oriented features to hide workspace arrays.
* Replace **COMMON** blocks with module variables.

-----------------------------------------------------------------------------

## Requirements
* The GNU Make tool https://www.gnu.org/software/make/
* The GNU gfortran compiler https://gcc.gnu.org/wiki/GFortran

-----------------------------------------------------------------------------

## To build the project

Type the following command line arguments
```
git clone https://github.com/jlokimlin/modern_fishpack.git

cd modern_fishpack; make all
```
