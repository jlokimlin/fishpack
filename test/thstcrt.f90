!     file thstcrt.f90
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
program thstcrt

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        FishpackSolver

    ! Explicit typing only
    implicit None

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    type (FishpackSolver)   :: solver
    integer (ip)            :: idimf, m, mbdcnd, n, nbdcnd, i, j, ierror
    real (wp), dimension(50, 53) :: f
    real (wp), dimension(53) :: bda, bdb, bdc, bdd
    real (wp), dimension(48) :: x
    real (wp), dimension(53) :: y
    real (wp) :: a, b, dx, c, d, dy, elmbda, pi, pi2, t, pertrb, discretization_error
    !-----------------------------------------------
    !
    !     from the dimension statement we get idimf = 50.
    !
    idimf = 50
    a = 1.
    b = 3.
    m = 48
    dx = (b - a)/m
    mbdcnd = 2
    c = -1.
    d = 1.
    n = 53
    dy = (d - c)/n
    nbdcnd = 0
    elmbda = -2.
    !
    !     auxiliary quantities
    !
    pi = acos(-1.0)
    pi2 = pi**2
    !
    !     generate and store grid points for computation of boundary data
    !     and the right side of the helmholtz equation.
    !
    do i = 1, m
        x(i) = a + (real(i, kind=wp) - 0.5_wp)*dx
    end do

    do j = 1, n
        y(j) = c + (real(j, kind=wp) - 0.5_wp)*dy
    end do
    !
    !     generate boundary data.
    !
    do j = 1, n
        bda(j) = 0.0_wp
        bdb(j) = -pi*cos(pi*y(j))
    end do
    !
    !     bdc and bdd are dummy arguments in this example.
    !
    !     generate right side of equation.
    !
    t = -2.0_wp*(pi2 + 1.0_wp)
    do i = 1, m
        do j = 1, n
            f(i, j) = t*sin(pi*x(i))*cos(pi*y(j))
        end do
    end do

    call solver%hstcrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, &
        bdc, bdd, elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error.  the exact solution is
    !
    !               u(x, y) = sin(pi*x)*cos(pi*y) .
    !
    discretization_error = 0.0_wp
    do i = 1, m
        do j = 1, n
            t = abs(f(i, j)-sin(pi*x(i))*cos(pi*y(j)))
            discretization_error = max(t, discretization_error)
        end do
    end do

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     hstcrt *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 1.2600e-3'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        ' discretization error = ', discretization_error

end program thstcrt
