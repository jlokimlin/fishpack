!     file thstplr.f
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
program thstplr

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        hstplr

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer (ip) :: idimf, m, mbdcnd, n, nbdcnd, i, j, ierror
    real (wp), dimension(51, 50) :: f
    real (wp), dimension(48) :: bda, bdb
    real (wp), dimension(50) :: bdc, bdd, r
    real (wp), dimension(48) :: theta
    real (wp) :: a, b, c, pi, d, elmbda, pertrb, discretization_error, z
    !-----------------------------------------------
    !
    !     from dimension statement we get value of idimf.
    !
    idimf = 51
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
    !     generate and store grid points for the purpose of computing
    !     boundary data and the right side of the poisson equation.
    !
    do i = 1, m
        r(i) = (real(i) - 0.5)/50.
    end do
    do j = 1, n
        theta(j) = (real(j) - 0.5)*pi/96.
    end do
    !
    !     generate boundary data.
    !
    do j = 1, n
        bdb(j) = 1. - cos(4.*theta(j))
    end do
    !
    !     generate boundary data.
    !
    bdc(:m) = 0.
    bdd(:m) = 0.
    !
    !     bda is a dummy variable.
    !
    !
    !     generate right side of equation.
    !
    do i = 1, m
        f(i, :n) = 16.*r(i)**2
    end do

    ! Solve system
    call hstplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error.  the exact solution is
    !
    !                u(r, theta) = r**4*(1 - cos(4*theta))
    !
    discretization_error = 0.0_wp
    do i = 1, m
        do j = 1, n
            z = abs(f(i, j)-(r(i)**4) * (1.0_wp - cos(4.0_wp * theta(j))))
            discretization_error = max(z, discretization_error)
        end do
    end do

    !     Print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     hstplr *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  discretization error = 1.1303e-3'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') &
        '     ierror =', ierror, ' discretization error = ', discretization_error

end program thstplr
