!
!     file thstssp.f90
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
!     *                      A Package of Fortran                     *
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
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        FishpackSolver

    ! Explicit typing only
    implicit None

    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    type (FishpackSolver)    :: solver
    integer (ip) :: m, mbdcnd, n, nbdcnd, idimf, i, j, ierror
    real (wp), dimension(18, 72) :: f
    real (wp), dimension(72) :: bda, bdb, bdc, bdd
    real (wp), dimension(18) :: sint
    real (wp), dimension(72) :: sinp
    real (wp) ::pi, a, b, c, d, elmbda, dtheta, dphi, pertrb, discretization_error, z
    !-----------------------------------------------
    !
    !     the value of idimf is the first dimension of f.
    !
    pi = acos(-1.0)
    a = 0.
    b = acos(0.0_wp)
    m = 18
    mbdcnd = 6
    c = 0.
    d = 2.0_wp*pi
    n = 72
    nbdcnd = 0
    elmbda = 0.0_wp
    idimf = 18
    !
    !     generate sines for use in subsequent computations
    !
    dtheta = b/m
    do i = 1, m
        sint(i) = sin((real(i, kind=wp) - 0.5_wp)*dtheta)
    end do
    dphi = d/n
    do j = 1, n
        sinp(j) = sin((real(j, kind=wp) - 0.5_wp)*dphi)
    end do
    !
    !     compute right side of equation and store in f
    !
    do j = 1, n
        f(:m, j) = 2.0_wp - 6.0_wp*(sint(:m)*sinp(j))**2
    end do
    !
    !     store derivative data at the equator
    !
    bdb(:n) = 0.0_wp
    !
    !     bda, bdc, and bdd are dummy variables.
    !
    call solver%hstssp(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error. since problem is singular, the
    !     solution must be normalized.
    !
    discretization_error = 0.0_wp
    do j = 1, n
        do i = 1, m
            z = abs(f(i, j)-(sint(i)*sinp(j))**2-f(1, 1))
            discretization_error = max(z, discretization_error)
        end do
    end do

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     hstssp *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  pertrb = 6.35830e-4'
    write( stdout, '(a)') '     discretization error = 3.37523e-3'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6)') '     ierror =', ierror, ' pertrb = ', pertrb
    write( stdout, '(a,1pe15.6/)') '     discretization error = ', discretization_error


end program thstssp
