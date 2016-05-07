!
!     file thstssp.f
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
program thstssp

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        FishpackSolver

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    type (FishpackSolver)    :: solver
    integer :: m, mbdcnd, n, nbdcnd, idimf, i, j, ierror
    real (wp), dimension(18, 72) :: f
    real (wp), dimension(72) :: bda, bdb, bdc, bdd
    real (wp), dimension(18) :: sint
    real (wp), dimension(72) :: sinp
    real::pi, a, b, c, d, elmbda, dtheta, dphi, pertrb, err, z
    !-----------------------------------------------
    !
    !     the value of idimf is the first dimension of f.
    !
    pi = acos( -1.0 )
    a = 0.
    b = pi/2.
    m = 18
    mbdcnd = 6
    c = 0.
    d = 2.*pi
    n = 72
    nbdcnd = 0
    elmbda = 0.
    idimf = 18
    !
    !     generate sines for use in subsequent computations
    !
    dtheta = b/m
    do i = 1, m
        sint(i) = sin((real(i) - 0.5)*dtheta)
    end do
    dphi = d/n
    do j = 1, n
        sinp(j) = sin((real(j) - 0.5)*dphi)
    end do
    !
    !     compute right side of equation and store in f
    !
    do j = 1, n
        f(:m, j) = 2. - 6.*(sint(:m)*sinp(j))**2
    end do
    !
    !     store derivative data at the equator
    !
    bdb(:n) = 0.
    !
    !     bda, bdc, and bdd are dummy variables.
    !
    call solver%hstssp(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error. since problem is singular, the
    !     solution must be normalized.
    !
    err = 0.
    do j = 1, n
        do i = 1, m
            z = abs(f(i, j)-(sint(i)*sinp(j))**2-f(1, 1))
            err = max(z, err)
        end do
    end do

    !     Print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     hstssp *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  pertrb = 6.35830e-4'
    write( stdout, '(A)') '     discretization error = 3.37523e-3'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') '     ierror =', ierror, ' pertrb = ', pertrb
    write( stdout, '(A,1pe15.6)') '     discretization error = ', err

end program thstssp
