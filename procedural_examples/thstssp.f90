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

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        ip, wp, HALF_PI, TWO_PI, hstssp

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    integer(ip), parameter :: M = 18
    integer(ip), parameter :: N = 72
    integer(ip), parameter :: IDIMF = M
    integer(ip)            :: mbdcnd, nbdcnd, i, j, ierror
    real(wp)               :: f(IDIMF,N), sint(M)
    real(wp), dimension(1) :: bda, bdc, bdd
    real(wp), dimension(N) :: sinp, bdb
    real(wp)               :: a, b, c, d, elmbda, dtheta, dphi
    real(wp)               :: pertrb
    real(wp), parameter    :: ZERO = 0.0_wp, HALF = 0.5_wp
    real(wp), parameter    :: TWO = 2.0_wp
    !-----------------------------------------------

    ! Set domain
    a = ZERO
    b = HALF_PI
    c = ZERO
    d = TWO_PI

    ! Set boundary conditions
    nbdcnd = 0
    mbdcnd = 6

    ! Set helmholtz constant
    elmbda = ZERO

    ! Set mesh sizes
    dtheta = (b - a)/M
    dphi = (d - c)/N

    ! Generate sines for use in subsequent computations
    do i = 1, M
        sint(i) = sin((real(i, kind=wp) - HALF) * dtheta)
    end do

    do j = 1, N
        sinp(j) = sin((real(j, kind=wp) - HALF) * dphi)
    end do

    ! Compute right hand side of equation and store in f
    block
        real(wp), parameter :: SIX = 6.0_wp
        do j = 1, N
            f(:, j) = TWO - SIX * (sint * sinp(j))**2
        end do
    end block
    !
    ! Set boundary data; store derivative data at the equator
    ! In our example, bda, bdc, and bdd are 1-dimensional dummy variables.
    bdb = ZERO

    ! Solve 2D Helmholtz in spherical coordinates on staggered grid
    call hstssp(a, b, M, mbdcnd, bda, bdb, c, d, N, nbdcnd, bdc, bdd, &
        elmbda, f, IDIMF, pertrb, ierror)

    ! Compute discretization error. since problem is singular, the
    ! solution must be normalized.
    block
        real(wp) :: discretization_error
        real(wp) :: exact_solution(M,N)

        do j = 1, N
            do i = 1, M
                exact_solution(i,j) = (sint(i) * sinp(j))**2 + f(1, 1)
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:M,:N)))

        ! Print earlier output from platforms with 64-bit floating point
        ! arithmetic followed by the output from this computer
        write( stdout, '(/a)') '     hstssp *** TEST RUN *** '
        write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(a)') '     ierror = 0,  pertrb = 6.35830e-4'
        write( stdout, '(a)') '     discretization error = 3.37523e-3'
        write( stdout, '(a)') '     The output from your computer is: '
        write( stdout, '(a,i3,a,1pe15.6)') '     ierror =', ierror, ' pertrb = ', pertrb
        write( stdout, '(a,1pe15.6/)') '     discretization error = ', discretization_error
    end block

end program thstssp
