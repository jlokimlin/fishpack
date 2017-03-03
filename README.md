# **fishpack - A Fortran library for solving of separable elliptic partial differential equations**
 
A modernization of NCAR's FISHPACK90.

* The original work, written in both FORTRAN 77 and Fortran 90, was heavily refactored to incorporate features of modern Fortran (2008+). 
* The library is now fully Fortran 2008 (ISO/IEC 1539-1:2010) compliant.
* Every common block, subroutine, and function is now encapsulated in a module. 
* Numerous memory leaks in the derived data type `fishworkspace`, now refactored as `FishpackWorkspace`, are resolved by replacing pointers with allocatable arrays. A valgrind run confirms the fix.
* This project is still a work in progress and every refactored solver passes their original unit test.

-----------------------------------------------------------------------------


## What is fishpack?

A collection of Fortran programs and subroutines that solve 2nd- and 4th-order finite difference approximations to separable elliptic Partial Differential Equations (PDEs). 

These include Helmholtz equations in cartesian, polar, cylindrical, and spherical coordinates, as well as more general separable elliptic equations. The solvers use the cyclic reduction algorithm. When the problem is singular, a least-squares solution is computed. Singularities induced by the coordinate system are handled, including at the origin *r=0* in cylindrical coordinates, and at the poles in spherical coordinates.

Test programs are provided for the 19 solvers in the `examples` folder. Each serves two purposes: as a template to guide you in writing your own codes utilizing the fishpack solvers, and as a demonstration that you can correctly produce the executables. 

-----------------------------------------------------------------------------

## Contributing

This project is still a work in progress and anyone is free to contribute under the proviso that they abstain from using the dreaded **go to**. 

For bug reports or feature requests please open an issue on github.

-----------------------------------------------------------------------------

## Requirements
* The GNU Make tool https://www.gnu.org/software/make/
* The GNU gfortran compiler https://gcc.gnu.org/wiki/GFortran

-----------------------------------------------------------------------------

## To build the project

Type the following command line arguments
```
git clone https://github.com/jlokimlin/fishpack.git; cd fishpack; make all
```

-----------------------------------------------------------------------------

## Result

```

     blktri *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 1.6478E-05
     The output from your computer is: 
     ierror =  0 discretization error =    1.647786E-05

     cblktri *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 1.6457E-05
     The output from your computer is: 
     ierror =  0 discretization error =    1.645720E-05

     cmgnbn *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 9.1620E-3
     The output from your computer is: 
     ierror =  0 discretization error =    9.161996E-03

     genbun *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 9.6406E-3
     The output from your computer is: 
     ierror =  0 discretization error =    9.640629E-03

     hstcrt *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 1.2600E-3
     The output from your computer is: 
     ierror =  0 discretization error =    1.260008E-03

     hstcsp *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 5.5843E-3
     The output from your computer is: 
     ierror =  0 discretization error =    5.584324E-03

     hstcyl *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  PERTRB =-4.4311E-4
     discretization error = 7.5280E-5 
     The output from your computer is: 
     ierror =  0 PERTRB =   -4.431139E-04
     discretization error =    7.527963E-05

     hstplr *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 1.1303E-3
     The output from your computer is: 
     ierror =  0 discretization error =    1.130379E-03

     hstssp *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  PERTRB = 6.35830E-4
     discretization error = 3.37523E-3
     The output from your computer is: 
     ierror =  0 PERTRB =    6.358300E-04
     discretization error =    3.375232E-03

     hw3crt *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 9.6480E-3
     The output from your computer is: 
     ierror =  0 discretization error =    9.648020E-03

     hwscrt *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 5.36508-4
     The output from your computer is: 
     ierror =  0 discretization error =    5.365082E-04

     hwscsp *** TEST RUN, EXAMPLE 1 *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 7.9984E-4 
     The output from your computer is: 
     ierror =  0     discretization error =   7.998416E-04

     hwscsp *** TEST RUN, EXAMPLE 2 *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0, discretization error = 5.8682E-5 
     The output from your computer is: 
     ierror =  0     discretization error =   5.868243E-05

     hwscyl *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  PERTRB = 2.2674E-4
     discretization error = 3.7367E-4 
     The output from your computer is: 
     ierror =  0 PERTRB =    2.267427E-04
     discretization error =    3.736722E-04

     hwsplr *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 6.19134E-4
     The output from your computer is: 
     ierror =  0 discretization error =    6.191342E-04

     hwsssp *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 3.38107E-3
     The output from your computer is: 
      ierror =  0 discretization error =    3.381073E-03

     pois3d *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 2.93277E-2
     The output from your computer is: 
     ierror =  0 discretization error =    2.932770E-02

     poistg *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 5.6417E-4
     The output from your computer is: 
     ierror =  0 discretization error =    5.641706E-04

     sepeli *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0
     Second Order discretization error = 9.7891E-5
     Fourth Order discretization error = 1.4735E-6
     The output from your computer is: 
     ierror =  0
     Second Order discretization error =   9.789104E-05
     Fourth Order discretization error =   1.473507E-06

     sepx4 *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0
     Second Order discretization error = 1.5985E-4
     Fourth Order discretization error = 1.8575E-6
     The output from your computer is: 
     ierror =  0
     Second Order discretization error =   1.598503E-04
     Fourth Order discretization error =   1.857488E-06


```
