# **modern\_fishpack**
 
This Fortran project is a modernization of NCAR's FISHPACK90 library. The original work, written in both FORTRAN 77 and Fortran 90, was heavily refactored to incorporate features of modern Fortran (2008+). 


Every common block, subroutine, and function is now encapsulated in a module. Most importantly, the potential memory leak in the derived type **fish**, now rebaptized as **FishpackWorkspace**, is resolved by replacing pointers with allocatable arrays. 

This project is still a work in progress and every refactored solver passes their original rest run.

-----------------------------------------------------------------------------

# TODO
* Introduce interfaces to replace assumed-size arrays with assumed shape arrays. 
* Replace all instances of **go to** statements with **exit**, **cycle** and **select case**
* Implement object-oriented features to hide workspace arrays.
* Introduce interfaces to address type casting in *blktri.f90*
* Parameterized kinds to remove the compiler flags **-fdefault-real-8 -fdefault-double-8** 

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
     ierror =           0  discretization error =    1.6457199287602535E-005
 
     cmgnbn *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 9.1620E-3
     The output from your computer is: 
     ierror =           0  discretization error =    9.1619959089139329E-003

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
     ierror =           0  PERTRB =    2.2674266866731276E-004
     discretization error =    3.7367223807960315E-004
 
     hwsplr *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 6.19134E-4
     The output from your computer is: 
     ierror =           0  discretization error =    6.1913422787440719E-004
 
     hwsssp *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 3.38107E-3
     The output from your computer is: 
     ierror =           0  discretization error =    1.0304665564689934     
 
     pois3d *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0,  discretization error = 2.93277E-2
     The output from your computer is: 
     ierror =           0  discretization error =    2.9327704986108594E-002

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
     ierror =           0
     Second Order discretization error =   9.7891036350938876E-005
     Fourth Order discretization error =   1.4735070190674548E-006
 
     sepx4 *** TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     ierror = 0
     Second Order discretization error = 1.5985E-4
     Fourth Order discretization error = 1.8575E-6
     The output from your computer is: 
     ierror =           0
     Second Order discretization error =   1.5985033721532815E-004
     Fourth Order discretization error =   1.8574879123711696E-006

************************************************************************

     *** TYPE (HelmholtzSovler) UNIT TESTS *** 


     To illustrates the usage of TYPE (HelmholtzSolver) type-bound procedure

     SOLVE_2D_HELMHOLTZ_CENTERED

     to solve the equation

     (d/dx)(du/dx) + (d/dy)(du/dy) - 4*u

     = (2 - (4 + pi**2/4) * x**2 ) * cos((y+1)*pi/2)

    with the boundary conditions on the rectangle

    (x,y) in [0,2]x[-1,3] with

    MIXED ROBIN CONDITION:

    u(0,y)       = 0                 for -1 <= y <= 3
    (du/dx)(2,y) = 4*cos((y+1)*pi/2) for -1 <= y <= 3

    and with u PERIODIC in y.

    the x-interval will be divided into 40 panels and the
    y-interval will be divided into 80 panels.


    *** TEST_SOLVE_2D_HELMHOLTZ_CENTERED ***    
    Previous 64 bit floating point arithmetic result 
    ierror = 0,  discretization error = 5.3650e-4
    Output from your computer is: 
    ierror =   0    discretization error =    5.365082E-04


   Illustrates the usage of TYPE (HelmholtzSolver) type-bound procedure

   SOLVE_2D_HELMHOLTZ_STAGGERED

   to solve the equation

    (d/dx)(du/dx) + (d/dy)(du/dy) - 2*u = -2(pi**2+1)sin(pi*x)cos(pi*y)

    with the boundary conditions on the rectangle

    (x,y) in [1,3]x[-1,1] with

    MIXED ROBIN CONDITION:

     u(1,y)       = 0,                for -1 <= y <= 3
     (du/dx)(3,y) = -pi*cos(pi*y),    for -1 <= y <= 3

     and u is PERIODIC in y.

     we want to have 48 unknowns in the x-interval and 53 unknowns
     in the y-interval.

     The exact solution to this problem is

         u(x,y) = sin(pi*x) * cos(pi*y).


    *** TEST_SOLVE_2D_HELMHOLTZ_STAGGERED ***    
    Previous 64 bit floating point arithmetic result 
    ierror = 0,  discretization error = 1.2600e-3
    The output from your computer is: 
    ierror =   0    discretization error =    1.260008E-03


************************************************************************

     *** TYPE (TridiagonalSolver) UNIT TESTS *** 


     To illustrate the usage of TYPE (TridiagonalSolver) type-bound procedure

     SOLVE_2D_REAL_LINEAR_SYSTEM_STAGGERED

     to solve the equation

     (1/cos(x))(d/dx)(cos(x)(du/dx)) + (d/dy)(du/dy) =

           2 * y**2 * (6-y**2) * sin(x)

     on the rectangle -pi/2 < x < pi/2 and 0 < y <  1

     with the boundary conditions

     (du/dx) (-pi/2, y) = (du/dx)(pi/2, y) = 0 , 0 <= y <= 1  (2)

     u(x, 0) = 0                                           (3)
                                 -pi/2 <= x <= pi/2
     (du/dy)(x, 1) = 4sin(x)                               (4)

     using finite differences on a staggered grid with
     deltax (= dx) = pi/40 and deltay (= dy) = 1/20 .
        to set up the finite difference equations we define
     the grid points

     x(i) = -pi/2 + (i-0.5)dx            i=1, 2, ..., 40

     y(j) = (j-0.5)dy                    j=1, 2, ..., 20

     and let v(i, j) be an approximation to u(x(i), y(j)).
     numbering the grid points in this fashion gives the set
     of unknowns as v(i, j) for i=1, 2, ..., 40 and j=1, 2, ..., 20.
     hence, in the program m = 40 and n = 20.  at the interior
     grid point (x(i), y(j)), we replace all derivatives in
     equation (1) by second order central finite differences,
     multiply by dy**2, and collect coefficients of v(i, j) to
     get the finite difference equation

     a(i)v(i-1, j) + b(i)v(i, j) + c(i)v(i+1, j)

     + v(i, j-1) - 2v(i, j) + v(i, j+1) = f(i, j)            (5)

     where s = (dy/dx)**2, and for i=2, 3, ..., 39

     a(i) = s * cos(x(i)-dx/2)

     b(i) = -s * (cos(x(i)-dx/2)+cos(x(i)+dx/2))

     c(i) = s * cos(x(i)+dx/2)

     f(i, j) = 2dy**2 * y(j)**2 * (6-y(j)**2) * sin(x(i)) , j=1, 2, ..., 19.

        to obtain equations for i = 1, we replace equation (2)
     by the second order approximation

     (v(1, j)-v(0, j))/dx = 0

     and use this equation to eliminate v(0, j) in equation (5)
     to arrive at the equation

     b(1)v(1, j) + c(1)v(2, j) + v(1, j-1) - 2v(1, j) + v(1, j+1)

                       = f(1, j)

     where

     b(1) = -s * (cos(x(1)-dx/2)+cos(x(1)+dx/2))

     c(1) = -b(1)

     for completeness, we set a(1) = 0.
        to obtain equations for i = 40, we replace the derivative
     in equation (2) at x=pi/2 in a similar fashion, use this
     equation to eliminate the virtual unknown v(41, j) in equation
     (5) and arrive at the equation

     a(40)v(39, j) + b(40)v(40, j)

     + v(40, j-1) - 2v(40, j) + v(40, j+1) = f(40, j)

     where

     a(40) = -b(40) = -s * (cos(x(40)-dx/2)+cos(x(40)+dx/2))

     for completeness, we set c(40) = 0.  hence, in the
     program mperod = 1.
        for j = 1, we replace equation (3) by the second order
     approximation

                (v(i, 0) + v(i, 1))/2 = 0

     to arrive at the condition

                v(i, 0) = -v(i, 1) .

     for j = 20, we replace equation (4) by the second order
     approximation

                (v(i, 21) - v(i, 20))/dy = 4 * sin(x)

     and combine this equation with equation (5) to arrive at
     the equation

     a(i)v(i-1, 20) + b(i)v(i, 20) + c(i)v(i+1, 20)

     + v(i, 19) - 2v(i, 20) + v(i, 21) = f(i, 20)

     where

     v(i, 21) = v(i, 20)  and

     f(i, 20) = 2 * dy**2 * y(j)**2 * (6-y(j)**2) * sin(x(i)) - 4 * dy * sin(x(i))

     hence, in the program nperod = 2.

     The exact solution to this problem is

        u(x, y) = y**4 * cos(x). 


    *** TEST_SOLVE_2D_REAL_LINEAR_SYSTEM_STAGGERED ***
    Previous 64 bit floating point arithmetic result 
    ierror = 0,  discretization error = 5.6417e-4
    The output from your computer is: 
    ierror =  0 discretization error =    5.641706E-04

     To illustrate the usage of TYPE (TridiagonalSolver) type-bound procedure

     SOLVE_2D_REAL_LINEAR_SYSTEM_CENTERED

     to solve the equation

     (1+x)**2*(d/dx)(du/dx) - 2(1+x)(du/dx) + (d/dy)(du/dy)

                  = 3(1+x)**4*sin(y)                      (1)

     on the rectangle 0 < x < 1 and -pi < y < pi
     with the boundary conditions

     (du/dx)(0,y) = 4sin(y)                               (2)
                                -pi <= y <= pi
     u(1,y) = 16sin(y)                                    (3)

     and with u periodic in y using finite differences on a
     grid with deltax (= dx) = 1/20 and deltay (= dy) = pi/20.
        to set up the finite difference equations we define
     the grid points

     x(i) = (i-1)dx            i=1,2,...,21

     y(j) = -pi + (j-1)dy      j=1,2,...,41

     and let v(i,j) be an approximation to u(x(i),y(j)).
     numbering the grid points in this fashion gives the set
     of unknowns as v(i,j) for i=1,2,...,20 and j=1,2,...,40.
     hence, in the program m = 20 and n = 40.  at the interior
     grid point (x(i),y(j)), we replace all derivatives in
     equation (1) by second order central finite differences,
     multiply by dy**2, and collect coefficients of v(i,j) to
     get the finite difference equation

     a(i)v(i-1,j) + b(i)v(i,j) + c(i)v(i+1,j)

     + v(i,j-1) - 2v(i,j) + v(i,j+1) = f(i,j)            (4)

     where s = (dy/dx)**2, and for i=2,3,...,19

     a(i) = (1+x(i))**2*s + (1+x(i))*s*dx

     b(i) = -2(1+x(i))**2*s

     c(i) = (1+x(i))**2*s - (1+x(i))*s*dx

     f(i,j) = 3(1+x(i))**4*dy**2*sin(y(j))  for j=1,2,...,40.

     to obtain equations for i = 1, we replace the
     derivative in equation (2) by a second order central
     finite difference approximation, use this equation to
     eliminate the virtual unknown v(0,j) in equation (4)
     and arrive at the equation

     b(1)v(1,j) + c(1)v(2,j) + v(1,j-1) - 2v(1,j) + v(1,j+1)

                       = f(1,j)

     where

     b(1) = -2s , c(1) = 2s

     f(1,j) = (11+8/dx)*dy**2*sin(y(j)),  j=1,2,...,40.

     for completeness, we set a(1) = 0.
        to obtain equations for i = 20, we incorporate
     equation (3) into equation (4) by setting

     v(21,j) = 16sin(y(j))

     and arrive at the equation

     a(20)v(19,j) + b(20)v(20,j)

     + v(20,j-1) - 2v(20,j) + v(20,j+1) = f(20,j)

     where

     a(20) = (1+x(20))**2*s + (1+x(20))*s*dx

     b(20) = -2*(1+x(20))**2*s

     f(20,j) = (3(1+x(20))**4*dy**2 - 16(1+x(20))**2*s
                + 16(1+x(20))*s*dx)*sin(y(j))

                    for j=1,2,...,40.

     for completeness, we set c(20) = 0.  hence, in the
     program mperod = 1.
        the periodicity condition on u gives the conditions

     v(i,0) = v(i,40) and v(i,41) = v(i,1) for i=1,2,...,20.

     hence, in the program nperod = 0.

     The exact solution to this problem is

                  u(x,y) = ((1+x)**4) * sin(y).


    *** TEST_SOLVE_2D_REAL_LINEAR_SYSTEM_CENTERED *** 
    Previous 64 bit floating point arithmetic result 
    ierror = 0,  discretization error = 9.6406e-3
    The output from your computer is: 
    ierror =  0 discretization error =    9.640629E-03
```
