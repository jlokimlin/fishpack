!
!     file thwsplr.f
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  Version 1.1                    *
!     *                                                               *
!     *                 A Package of Fortran 77 and 90                *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *               for Modeling Geophysical Processes              *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        John Adams, Paul Swarztrauber and Roland Sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!     Purpose:
!
!     Program to illustrate the use of subroutine hwsplr to solve
!     the equation
!
!     (1/r)(d/dr)(r*(du/dr)) + (1/r**2)(d/dtheta)(du/dtheta) = 16*r**2
!
!     on the quarter-disk 0 .lt. r .lt. 1, 0 .lt. theta .lt. pi/2 with
!     with the boundary conditions
!
!     u(1, theta) = 1 - cos(4*theta), 0 .le. theta .le. 1
!
!     and
!
!     (du/dtheta)(r, 0) = (du/dtheta)(r, pi/2) = 0,  0 .le. r .le. 1.
!
!     (note that the solution u is unspecified at r = 0.)
!          the r-interval will be divided into 50 panels and the
!     theta-interval will be divided into 48 panels.
!
!
!     from dimension statement we get value of idimf.
!
program thwsplr

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        hwsplr

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer :: idimf, m, mbdcnd, n, nbdcnd, mp1, np1, i, j, ierror
    real , dimension(100, 50) :: f
    real , dimension(51) :: bdc, bdd, r, bda, bdb
    real , dimension(49) :: theta
    real :: a, b, c, pi, dum, d, elmbda, pertrb, discretization_error, z
    !-----------------------------------------------

    idimf = 100
    a = 0.
    b = 1.
    m = 50
    mbdcnd = 5
    c = 0.
    pi = acos( -1.0 )
    d = pi/2.
    n = 48
    nbdcnd = 3
    elmbda = 0.
    !
    !     auxiliary quantities.
    !
    mp1 = m + 1
    np1 = n + 1
    !
    !     generate and store grid points for the purpose of computing
    !     boundary data and the right side of the poisson equation.
    !
    do i = 1, mp1
        r(i) = real(i - 1)/50.
    end do
    do j = 1, np1
        theta(j) = real(j - 1)*pi/96.
    end do
    !
    !     generate boundary data.
    !
    bdc(:mp1) = 0.
    bdd(:mp1) = 0.
    !
    !     bda and bdb are dummy variables.
    !
    do j = 1, np1
        f(mp1, j) = 1. - cos(4.*theta(j))
    end do
    !
    !     generate right side of equation.
    !
    do i = 1, m
        f(i, :np1) = 16.*r(i)**2
    end do

    ! Solve system
    call hwsplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error.  the exact solution is
    !                u(r, theta) = r**4*(1 - cos(4*theta))
    !
    discretization_error = 0.
    do i = 1, mp1
        do j = 1, np1
            z = abs(f(i, j)-r(i)**4*(1.-cos(4.*theta(j))))
            discretization_error = max(z, discretization_error)
        end do
    end do

    !     Print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     hwsplr *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  discretization error = 6.19134e-4'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') &
        '     ierror =', ierror, ' discretization error = ', discretization_error

end program thwsplr
