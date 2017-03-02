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
!     Purpose:
!
!     To illustrate the usage of poistg to solve the equation
!
!     (1/cos(x))(d/dx)(cos(x)(du/dx)) + (d/dy)(du/dy) =
!
!           2 * y**2 * (6-y**2) * sin(x)
!
!     on the rectangle -pi/2 < x < pi/2 and 0 < y <  1
!
!     with the boundary conditions
!
!     (du/dx) (-pi/2, y) = (du/dx)(pi/2, y) = 0 , 0 <= y <= 1  (2)
!
!     u(x, 0) = 0                                           (3)
!                                 -pi/2 <= x <= pi/2
!     (du/dy)(x, 1) = 4sin(x)                               (4)
!
!     using finite differences on a staggered grid with
!     deltax (= dx) = pi/40 and deltay (= dy) = 1/20 .
!        to set up the finite difference equations we define
!     the grid points
!
!     x(i) = -pi/2 + (i-0.5)dx            i=1, 2, ..., 40
!
!     y(j) = (j-0.5)dy                    j=1, 2, ..., 20
!
!     and let v(i, j) be an approximation to u(x(i), y(j)).
!     numbering the grid points in this fashion gives the set
!     of unknowns as v(i, j) for i=1, 2, ..., 40 and j=1, 2, ..., 20.
!     hence, in the program m = 40 and n = 20.  at the interior
!     grid point (x(i), y(j)), we replace all derivatives in
!     equation (1) by second order central finite differences,
!     multiply by dy**2, and collect coefficients of v(i, j) to
!     get the finite difference equation
!
!     a(i)v(i-1, j) + b(i)v(i, j) + c(i)v(i+1, j)
!
!     + v(i, j-1) - 2v(i, j) + v(i, j+1) = f(i, j)            (5)
!
!     where s = (dy/dx)**2, and for i=2, 3, ..., 39
!
!     a(i) = s * cos(x(i)-dx/2)
!
!     b(i) = -s * (cos(x(i)-dx/2)+cos(x(i)+dx/2))
!
!     c(i) = s * cos(x(i)+dx/2)
!
!     f(i, j) = 2dy**2 * y(j)**2 * (6-y(j)**2) * sin(x(i)) , j=1, 2, ..., 19.
!
!        to obtain equations for i = 1, we replace equation (2)
!     by the second order approximation
!
!     (v(1, j)-v(0, j))/dx = 0
!
!     and use this equation to eliminate v(0, j) in equation (5)
!     to arrive at the equation
!
!     b(1)v(1, j) + c(1)v(2, j) + v(1, j-1) - 2v(1, j) + v(1, j+1)
!
!                       = f(1, j)
!
!     where
!
!     b(1) = -s * (cos(x(1)-dx/2)+cos(x(1)+dx/2))
!
!     c(1) = -b(1)
!
!     for completeness, we set a(1) = 0.
!        to obtain equations for i = 40, we replace the derivative
!     in equation (2) at x=pi/2 in a similar fashion, use this
!     equation to eliminate the virtual unknown v(41, j) in equation
!     (5) and arrive at the equation
!
!     a(40)v(39, j) + b(40)v(40, j)
!
!     + v(40, j-1) - 2v(40, j) + v(40, j+1) = f(40, j)
!
!     where
!
!     a(40) = -b(40) = -s * (cos(x(40)-dx/2)+cos(x(40)+dx/2))
!
!     for completeness, we set c(40) = 0.  hence, in the
!     program nperod = 1.
!
!     For j = 1, we replace equation (3) by the second order
!     approximation
!
!                (v(i, 0) + v(i, 1))/2 = 0
!
!     to arrive at the condition
!
!                v(i, 0) = -v(i, 1) .
!
!     For j = 20, we replace equation (4) by the second order
!     approximation
!
!                (v(i, 21) - v(i, 20))/dy = 4 * sin(x)
!
!     and combine this equation with equation (5) to arrive at
!     the equation
!
!     a(i)v(i-1, 20) + b(i)v(i, 20) + c(i)v(i+1, 20)
!
!     + v(i, 19) - 2v(i, 20) + v(i, 21) = f(i, 20)
!
!     where
!
!     v(i, 21) = v(i, 20)  and
!
!     f(i, 20) = 2 * dy**2 * y(j)**2 * (6-y(j)**2) * sin(x(i)) - 4 * dy * sin(x(i))
!
!     hence, in the program nperod = 2.
!
!     The exact solution to this problem is
!
!     u(x, y) = y**4 * cos(x) .
!
!     From dimension statement we get that size(f, dim=1) = 42
!
program test_poistg

    use fishpack

    ! Explicit typing only
    implicit none

    ! Dictionary
    integer(ip), parameter        :: M = 40, N = 20
    integer(ip), parameter        :: IDIMF = M + 2
    integer(ip)                   :: mperod, nperod, i, j, ierror
    real(wp), dimension(IDIMF, N) :: f
    real(wp), dimension(M)        :: a, b, c, x
    real(wp)                      :: y(N), dx, dy
    real(wp), parameter           :: ZERO = 0.0_wp, HALF = 0.5_wp
    real(wp), parameter           :: ONE = 1.0_wp, TWO = 2.0_wp

    ! Set boundary conditions
    mperod = 1
    nperod = 2

    ! Set mesh sizes
    dx = PI/M
    dy = ONE/N

    ! Generate and store grid points for computation
    do i = 1, M
        x(i) = -HALF_PI + (real(i, kind=wp) - HALF)*dx
    end do

    do j = 1, N
        y(j) = (real(j, kind=wp) - HALF)*dy
    end do

    ! Generate coefficients
    block
        real(wp) :: s, half_dx

        s = (dy/dx)**2
        half_dx = dx/TWO
        a(1) = ZERO
        b(1) = -s * cos(-HALF_PI + dx)/cos(x(1))
        c(1) = -b(1)

        do i = 2, M
            a(i) = s * cos(x(i)-half_dx)/cos(x(i))
            c(i) = s * cos(x(i)+half_dx)/cos(x(i))
            b(i) = -(a(i) + c(i))
        end do

        a(M) = -b(M)
        c(M) = ZERO
    end block

    ! Generate right hand side of equation
    block
        real(wp), parameter :: FOUR = 4.0_wp
        real(wp), parameter :: SIX = 6.0_wp

        do j = 1, N
            do i = 1, M
                f(i, j) = TWO * (dy**2) * (y(j)**2) * (SIX - y(j)**2) * sin(x(i))
            end do
        end do
        f(1:M,N) = f(1:M,N) - FOUR * dy * sin(x)
    end block

    ! Solve 2D real linear system on staggered grid
    call poistg(nperod, N, mperod, M, a, b, c, IDIMF, f, ierror)

    ! Compute discretization error. The exact solution is
    !
    !  u(x,y) = (y**4) * sin(x)
    !
    block
        real(wp), parameter :: KNOWN_ERROR = 0.564170618941665e-003
        real(wp) :: discretization_error
        real(wp) :: exact_solution(M,N)

        do j = 1, N
            do i = 1, M
                exact_solution(i,j) = (y(j)**4) * sin(x(i))
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:M,:N)))

        call check_output('poistg', ierror, KNOWN_ERROR, discretization_error)
    end block

end program test_poistg
