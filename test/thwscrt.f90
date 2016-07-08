!
!     file thwscrt.f90
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
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
program thwscrt

    use, intrinsic :: iso_fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        PI, &
        FishpackSolver

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    type (FishpackSolver)    :: solver
    integer (ip) :: idimf, m, mbdcnd, n, nbdcnd, mp1, np1, i, j, ierror
    real (wp), dimension(45, 82) :: f
    real (wp), dimension(81) :: bdb, bda, bdc, bdd, y
    real (wp), dimension(41) :: x
    real (wp) :: a, b, c, d, elmbda, HALF_PI, PI2, pertrb, discretization_error, z
    !-----------------------------------------------
    !
    !     from dimension statement we get value of idimf.
    !
    idimf = 45
    a = 0.0_wp
    b = 2.0_wp
    m = 40
    mbdcnd = 2
    c = -1.0_wp
    d = 3.0_wp
    n = 80
    nbdcnd = 0
    elmbda = -4.0_wp
    !
    !     auxiliary quantities.
    !
    HALF_PI = PI/2
    PI2 = PI**2
    mp1 = m + 1
    np1 = n + 1
    !
    !     generate and store grid points for the purpose of computing
    !     boundary data and the right side of the helmholtz equation.
    !
    do i = 1, mp1
        x(i) = real(i - 1, kind=wp)/20
    end do

    do j = 1, np1
        y(j) = -1.0_wp + real(j - 1, kind=wp)/20
    end do
    !
    !     generate boundary data.
    !
    do j = 1, np1
        bdb(j) = 4.0_wp * cos((y(j) + 1.0_wp)*HALF_PI)
    end do
    !
    !     bda, bdc, and bdd are dummy variables.
    !
    f(1, :np1) = 0.0_wp
    !
    !     generate right side of equation.
    !
    do i = 2, mp1
        do j = 1, np1
            f(i, j) = (2.0_wp - (4.0_wp + PI2/4)*(x(i)**2)) * cos((y(j) + 1.0_wp)*HALF_PI)
        end do
    end do

    call solver%hwscrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error.  the exact solution is
    !                u(x, y) = x**2*cos((y+1)*(pi/2))
    !
    discretization_error = 0.0_wp
    do i = 1, mp1
        do j = 1, np1
            z = abs(f(i, j)-(x(i)**2) * cos((y(j)+1.0_wp) * HALF_PI))
            discretization_error = max(z, discretization_error)
        end do
    end do

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     hwscrt *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 5.36508e-4'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        ' discretization error = ', discretization_error

end program thwscrt
