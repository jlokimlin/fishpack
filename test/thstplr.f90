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
program thstplr

    use, intrinsic :: iso_fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        FishpackSolver

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    type (FishpackSolver)    :: solver
    integer (ip)             :: idimf, m, mbdcnd, n, nbdcnd, i, j, ierror
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
    pi = acos(-1.0_wp)
    d = acos(0.0_wp)
    n = 48
    nbdcnd = 3
    elmbda = 0.0_wp
    !
    !     generate and store grid points for the purpose of computing
    !     boundary data and the right side of the poisson equation.
    !
    do i = 1, m
        r(i) = (real(i, kind=wp) - 0.5_wp)/50
    end do

    do j = 1, n
        theta(j) = (real(j, kind=wp) - 0.5_wp)*pi/96
    end do
    !
    !     generate boundary data.
    !
    do j = 1, n
        bdb(j) = 1.0_wp - cos(4.0_wp*theta(j))
    end do
    !
    !     generate boundary data.
    !
    bdc(:m) = 0.0_wp
    bdd(:m) = 0.0_wp
    !
    !     bda is a dummy variable.
    !
    !
    !     generate right side of equation.
    !
    do i = 1, m
        f(i, :n) = 16.0_wp*r(i)**2
    end do

    ! Solve system
    call solver%hstplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
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

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     hstplr *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 1.1303e-3'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        ' discretization error = ', discretization_error

end program thstplr
