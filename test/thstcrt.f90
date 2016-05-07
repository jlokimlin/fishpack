!     file thstcrt.f
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
program thstcrt

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        hstcrt

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer (ip) :: idimf, m, mbdcnd, n, nbdcnd, i, j, ierror
    real (wp), dimension(50, 53) :: f
    real (wp), dimension(53) :: bda, bdb, bdc, bdd
    real (wp), dimension(48) :: x
    real (wp), dimension(53) :: y
    real (wp) :: a, b, dx, c, d, dy, elmbda, pi, pisq, t, pertrb, discretization_error
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
    pi = acos( -1.0 )
    pisq = pi*pi
    !
    !     generate and store grid points for computation of boundary data
    !     and the right side of the helmholtz equation.
    !
    do i = 1, m
        x(i) = a + (real(i) - 0.5)*dx
    end do
    do j = 1, n
        y(j) = c + (real(j) - 0.5)*dy
    end do
    !
    !     generate boundary data.
    !
    do j = 1, n
        bda(j) = 0.
        bdb(j) = -pi*cos(pi*y(j))
    end do
    !
    !     bdc and bdd are dummy arguments in this example.
    !
    !     generate right side of equation.
    !
    t = -2.*(pisq + 1.)
    do i = 1, m
        do j = 1, n
            f(i, j) = t*sin(pi*x(i))*cos(pi*y(j))
        end do
    end do

    call hstcrt (a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd &
        , elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error.  the exact solution is
    !
    !               u(x, y) = sin(pi*x)*cos(pi*y) .
    !
    discretization_error = 0.
    do i = 1, m
        do j = 1, n
            t = abs(f(i, j)-sin(pi*x(i))*cos(pi*y(j)))
            discretization_error = max(t, discretization_error)
        end do
    end do
    !     Print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     hstcrt *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  discretization error = 1.2600e-3'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') &
        '     ierror =', ierror, ' discretization error = ', discretization_error

end program thstcrt
