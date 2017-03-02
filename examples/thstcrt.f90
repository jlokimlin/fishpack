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
program test_hstcrt

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack

    ! Explicit typing only
    implicit none

    ! Dictionary
    integer(ip), parameter :: M = 48
    integer(ip), parameter :: N = 53
    integer(ip), parameter :: IDIMF = M + 2
    integer(ip)            :: mbdcnd, nbdcnd, i, j, ierror
    real(wp)               :: f(IDIMF, N)
    real(wp), dimension(N) :: bda, bdb, y
    real(wp)               :: a, b, c, d
    real(wp)               :: dx, dy, x(M), bdc(1), bdd(1)
    real(wp)               :: elmbda, pertrb
    real(wp), parameter    :: PI2 = PI**2
    real(wp), parameter    :: ZERO = 0.0_wp, HALF = 0.5_wp, ONE = 1.0_wp
    real(wp), parameter    :: TWO = 2.0_wp, THREE = 3.0_wp

    ! Set domain
    a = ONE
    b = THREE
    c = -ONE
    d = ONE

    ! Set boundary conditions
    mbdcnd = 2
    nbdcnd = 0

    ! Set helmholtz constant
    elmbda = -TWO

    ! Set mesh sizes
    dx = (b - a)/M
    dy = (d - c)/N

    ! Generate and store grid points for computation of boundary data
    ! and the right side of the helmholtz equation.
    do i = 1, M
        x(i) = a + (real(i, kind=wp) - HALF) * dx
    end do

    do j = 1, N
        y(j) = c + (real(j, kind=wp) - HALF) * dy
    end do

    ! Generate boundary data. bdc and bdd are dummy arguments in this example.
    bda = ZERO
    bdb = -PI * cos(PI*y)

    ! Generate right side of equation
    block
        real(wp), parameter :: CONST = -TWO*(PI2 + ONE)

        do i = 1, M
            f(i,:) = CONST*sin(PI*x(i))*cos(PI*y(:))
        end do
    end block

    ! Solve 2D Helmholtz in cartesian coordinates on staggered grid
    call hstcrt(a, b, M, mbdcnd, bda, bdb, c, d, N, nbdcnd, &
        bdc, bdd, elmbda, f, IDIMF, pertrb, ierror)

    ! Compute discretization error. The exact solution is
    !
    ! u(x, y) = sin(pi*x) * cos(pi*y).
    block
        real(wp), parameter :: KNOWN_ERROR = 0.126000790151959e-2_wp
        real(wp) :: discretization_error
        real(wp) :: exact_solution(M,N)

        do j = 1, N
            do i = 1, M
                exact_solution(i,j) = sin(PI*x(i))*cos(PI*y(j))
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:M,:N)))

        call check_output('hstcrt', ierror, KNOWN_ERROR, discretization_error)
    end block

end program test_hstcrt
