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

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        ip, wp, HALF_PI, hstplr

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer(ip), parameter :: M = 50, N = 48
    integer(ip), parameter :: IDIMF = M + 1, NP1 = N + 1
    integer(ip)            :: mbdcnd, nbdcnd, i, j, ierror
    real(wp)               :: f(IDIMF, NP1)
    real(wp), dimension(N) :: bdb, theta
    real(wp), dimension(M) :: bdc, bdd, r
    real(wp)               :: a, b, c, d, bda(1), dr, dtheta
    real(wp)               :: elmbda, pertrb
    real(wp), parameter    :: ZERO = 0.0_wp, HALF = 0.5_wp
    real(wp), parameter    :: ONE = 1.0_wp, FOUR = 4.0_wp
    !-----------------------------------------------

    ! Set domain
    a = ZERO
    b = ONE
    c = ZERO
    d = HALF_PI

    ! Set boundary conditions
    mbdcnd = 5
    nbdcnd = 3

    ! Set helmholtz constant
    elmbda = ZERO

    ! Set mesh sizes
    dr = (b - a)/M
    dtheta = (d - c)/N

    ! Generate and store grid points for the purpose of computing
    ! boundary data and the right side of the poisson equation.
    do i = 1, M
        r(i) = (real(i, kind=wp) - HALF) * dr
    end do

    do j = 1, N
        theta(j) = (real(j, kind=wp) - HALF) * dtheta
    end do

    ! Generate boundary data. In our example, bda is a 1-dimensional dummy variable.
    bdb = ONE - cos(FOUR*theta)
    bdc = ZERO
    bdd = ZERO

    ! Generate right side of equation.
    block
        real(wp), parameter :: SIXTEEN = 16.0_wp

        do i = 1, M
            f(i, :N) = SIXTEEN * r(i)**2
        end do
    end block

    ! Solve 2D Helmholtz in polar coordinates on staggered grid
    call hstplr(a, b, M, mbdcnd, bda, bdb, c, d, N, nbdcnd, bdc, bdd, &
        elmbda, f, IDIMF, pertrb, ierror)

    ! Compute discretization error.  the exact solution is
    !
    ! u(r, theta) = r**4*(1 - cos(4*theta))
    !
    block
        real(wp) :: discretization_error
        real(wp) :: exact_solution(M,N)

        do j = 1, N
            do i = 1, M
                exact_solution(i,j) = (r(i)**4) * (ONE - cos(FOUR * theta(j)))
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:M,:N)))

        ! Print earlier output from platforms with 64-bit floating point
        ! arithmetic followed by the output from this computer
        write( stdout, '(/a)') '     hstplr *** TEST RUN *** '
        write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(a)') '     ierror = 0,  discretization error = 1.1303e-3'
        write( stdout, '(a)') '     The output from your computer is: '
        write( stdout, '(a,i3,a,1pe15.6/)') &
            '     ierror =', ierror, &
            ' discretization error = ', discretization_error
    end block

end program thstplr
