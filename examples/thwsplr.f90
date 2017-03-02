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
!
!     Purpose:
!
!     Program to illustrate the use of subroutine hwsplr to solve
!     the equation
!
!     (1/r)(d/dr)(r*(du/dr)) + (1/r**2)(d/dtheta)(du/dtheta) = 16*r**2
!
!     on the quarter-disk 0 < r < 1, 0 < theta < pi/2 with
!     with the boundary conditions
!
!     u(1, theta) = 1 - cos(4*theta), 0 <= theta <= 1
!
!     and
!
!     (du/dtheta)(r, 0) = (du/dtheta)(r, pi/2) = 0,  0 <= r <= 1.
!
!     (note that the solution u is unspecified at r = 0.)
!          the r-interval will be divided into 50 panels and the
!     theta-interval will be divided into 48 panels.
!
!
!     from dimension statement we get value of idimf.
!
program test_hwsplr

    use fishpack

    ! Explicit typing only
    implicit none

    ! Dictionary
    integer(ip), parameter    :: m = 50, n = 48
    integer(ip), parameter    :: idimf = 100
    integer(ip), parameter    :: mp1 = m + 1, np1 = n + 1
    integer(ip)               :: mbdcnd, nbdcnd, i, j, ierror
    real(wp)                  :: f(idimf, m), theta(np1)
    real(wp), dimension (mp1) :: bdc, bdd, r
    real(wp), dimension(1)    :: bda, bdb
    real(wp)                  :: a, b, c, d, elmbda, pertrb
    real(wp), parameter       :: ZERO = 0.0_wp, ONE = 1.0_wp, FOUR = 4.0_wp

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

    ! Generate and store grid points for the purpose of computing
    ! boundary data and the right side of the poisson equation.
    do i = 1, mp1
        r(i) = real(i - 1, kind=wp)/m
    end do

    do j = 1, np1
        theta(j) = real(j - 1, kind=wp) * (HALF_PI/n)
    end do

    ! Generate boundary data.
    ! In our example, bda and bdb are 1-dimensional dummy variables.
    bdc = ZERO
    bdd = ZERO

    ! Generate right hand side of equation
    block
        real(wp), parameter :: SIXTEEN = 16.0_wp

        do j = 1, np1
            f(mp1, j) = ONE - cos(FOUR * theta(j))
        end do

        do i = 1, m
            f(i, :np1) = SIXTEEN * r(i)**2
        end do
    end block

    ! Solve 2D Helmholtz in polar coordinates on centered grid
    call hwsplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)

    ! Compute discretization error. The exact solution is
    !
    ! u(r,theta) = (r**4)*(1 - cos(4*theta))
    block
        real(wp), parameter :: KNOWN_ERROR = 0.619134227874629e-003
        real(wp) :: discretization_error, exact_solution(mp1, np1)

        do j = 1, np1
            do i = 1, mp1
                exact_solution(i,j) = (r(i)**4)*(ONE-cos(FOUR * theta(j)))
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:mp1,:np1)))

        call check_output('hwsplr', ierror, KNOWN_ERROR, discretization_error)
    end block

end program test_hwsplr
