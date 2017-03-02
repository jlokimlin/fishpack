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
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
program test_hwscrt

    use fishpack

    ! Explicit typing only
    implicit none

    ! Dictionary
    integer(ip), parameter   :: M = 40, N = 80
    integer(ip), parameter   :: MP1 = M + 1, NP1 = N + 1
    integer(ip), parameter   :: IDIMF = M + 5
    integer(ip)              :: mbdcnd, nbdcnd, i, j, ierror
    real(wp)                 :: f(IDIMF, NP1), x(MP1)
    real(wp), dimension(NP1) :: bdb, y
    real(wp), dimension(1)   :: bda, bdc, bdd
    real(wp)                 :: a, b, c, d, elmbda, pertrb
    real(wp), parameter      :: PI2 = PI**2
    real(wp), parameter      :: ZERO = 0.0_wp, ONE = 1.0_wp
    real(wp), parameter      :: TWO = 2.0_wp, THREE = 3.0_wp, FOUR = 4.0_wp

    ! Set domain
    a = ZERO
    b = TWO
    c = -ONE
    d = THREE

    ! Set boundary conditions
    mbdcnd = 2
    nbdcnd = 0

    ! Set helmholtz constant
    elmbda = -FOUR

    ! Generate and store grid points for the purpose of computing
    ! boundary data and the right side of the helmholtz equation.
    do i = 1, MP1
        x(i) = TWO * real(i - 1, kind=wp)/M
    end do

    do j = 1, NP1
        y(j) = -ONE + FOUR * real(j - 1, kind=wp)/N
    end do

    ! Generate boundary data.
    ! In our example, bda, bdc, and bdd are 1-dimensional dummy variables.
    bdb = FOUR * cos((y + ONE) * HALF_PI)

    ! Generate right side of equation.
    f = ZERO
    do j = 1, NP1
        do i = 2, MP1
            f(i, j) = (TWO - (FOUR + PI2/4)*(x(i)**2)) * cos((y(j) + ONE)*HALF_PI)
        end do
    end do

    ! Solve 2D Helmholtz in cartesian coordinates on a centered grid
    call hwscrt(a, b, M, mbdcnd, bda, bdb, c, d, N, nbdcnd, bdc, bdd, &
        elmbda, f, IDIMF, pertrb, ierror)

    ! Compute discretization error. The exact solution is
    !
    ! u(x, y) = (x**2) * cos((y+1)*(pi/2))
    block
        real(wp), parameter :: KNOWN_ERROR = 0.536508246868017e-3_wp
        real(wp) :: discretization_error
        real(wp) :: exact_solution(MP1, NP1)

        do j = 1, NP1
            do i = 1, MP1
                exact_solution(i,j) = (x(i)**2) * cos((y(j)+ONE) * HALF_PI)
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:MP1,:NP1)))

        call check_output('hwscrt', ierror, KNOWN_ERROR, discretization_error)
    end block

end program test_hwscrt
