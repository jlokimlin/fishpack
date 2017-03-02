!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                         Fishpack                              *
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
program test_hstcsp

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    type(FishpackWorkspace) :: workspace
    integer(ip), parameter  :: M = 45
    integer(ip), parameter  :: N = 15
    integer(ip), parameter  :: IDIMF = M + 2, NP1 = N + 1
    integer(ip)             :: mbdcnd, nbdcnd, i, j, intl, ierror
    real(wp)                :: f(IDIMF, NP1), r(N)
    real(wp), dimension(M)  :: bdd, theta, cost
    real(wp), dimension(1)  :: bda, bdb, bdc
    real(wp)                :: a, b, c, d
    real(wp)                :: dt, dr, elmbda, pertrb
    real(wp), parameter     :: ZERO = 0.0_wp, HALF = 0.5_wp, ONE = 1.0_wp
    !-----------------------------------------------

    ! Set domaiin
    a = ZERO
    b = PI
    c = ZERO
    d = ONE

    ! Set boundary conditions
    mbdcnd = 9
    nbdcnd = 5

    ! Set helmholtz constant
    elmbda = ZERO

    ! Define mesh sizes
    dt = (b - a)/M
    dr = (d - c)/N

    ! Define grid points theta(i) and cos(theta(i))
    do i = 1, M
        theta(i) = a + (real(i, kind=wp) - HALF) * dt
        cost(i) = cos(theta(i))
    end do

    ! Define grid points r(j)
    do j = 1, N
        r(j) = c + (real(j) - HALF) * dr
    end do

    ! Define boundary array bdd.
    ! bda, bdb, and bdc are one-dimensional dummy variables in this example
    bdd = cost**4

    ! Define right hand side of equation
    block
        real(wp), parameter :: TWELVE = 12.0_wp

        do i = 1, M
            f(i,:N) = TWELVE * (r(:N)*cost(i))**2
        end do
    end block

    ! Initialize call flag
    intl = 0

    ! Solve 2D axisymmetric Helmholtz equation on staggered grid
    call hstcsp(intl, a, b, M, mbdcnd, bda, bdb, c, d, N, nbdcnd, bdc, &
        bdd, elmbda, f, IDIMF, pertrb, ierror, workspace)

    ! Compute discretization error. The exact solution is
    !
    ! u(theta, r) = (r*cos(theta))**4
    !
    block
        real(wp), parameter :: KNOWN_ERROR = 0.558432375660800e-2_wp
        real(wp) :: discretization_error
        real(wp) :: exact_solution(M,N)

        do j = 1, N
            do i = 1, M
                exact_solution(i,j) = (r(j)*cost(i))**4
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:M,:N)))

        call check_output('hstcsp', ierror, KNOWN_ERROR, discretization_error)
    end block

    ! Release memory
    call workspace%destroy()

end program test_hstcsp
