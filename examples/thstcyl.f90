!     file thstcyl.f90
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
!     program to illustrate the use of hstcyl to solve the equation
!
!    (1/r)(d/dr)(r*du/dr) + (d/dz)(du/dz) = (2*r*z)**2*(4*z**2 + 3*r**2)
!
!     on the rectangle 0 < r < 1 , 0 < z < 1 with the
!     boundary conditions
!
!     (du/dr)(1, z) = 4*z**2  for  0 <= z <= 1
!
!     and
!
!     (du/dz)(r, 0) = 0 and (du/dz)(r, 1) = 4*r**2  for  0 <= r <= 1 .
!
!     the solution to this problem is not unique.  it is a
!     one-parameter family of solutions given by
!
!            u(r, z) = (r*z)**4 + arbitrary constant .
!
!     the r-interval will contain 50 unknowns and the z-interval will
!     contain 52 unknowns.
!
program test_hstcyl

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer(ip), parameter :: M = 50
    integer(ip), parameter :: N = 52
    integer(ip), parameter :: IDIMF = M + 1
    integer(ip)            :: mbdcnd, nbdcnd, i, j, ierror
    real(wp)               :: f(IDIMF,N), pertrb
    real(wp), dimension(N) :: bdb, z
    real(wp), dimension(M) :: bdc, bdd, r
    real(wp)               :: a, b, c, d, elmbda, bda(1)
    real(wp), parameter    :: ZERO = 0.0_wp, HALF = 0.5_wp
    real(wp), parameter    :: ONE = 1.0_wp, FOUR = 4.0_wp
    !-----------------------------------------------

    ! Set domain
    a = ZERO
    b = ONE
    c = ZERO
    d = ONE

    ! Set boundary conditions
    mbdcnd = 6
    nbdcnd = 3

    ! Set helmholtz constant
    elmbda = ZERO

    ! Generate and store grid points for the purpose of computing
    ! boundary data and the right side of the poisson equation.
    do i = 1, M
        r(i) = (real(i, kind=wp) - HALF)/M
    end do

    do j = 1, N
        z(j) = (real(j, kind=wp) - HALF)/N
    end do

    ! Generate boundary data. bda is a 1-dimensional dummy variable
    bdb = FOUR * (z**4)
    bdc = ZERO
    bdd = FOUR * (r**4)

    ! Generate right side of equation.
    block
        real(wp), parameter :: THREE = 3.0_wp

        do i = 1, M
            f(i,:N) = FOUR * (r(i)**2) * (z**2) * (FOUR * (z**2) + THREE * (r(i)**2) )
        end do
    end block

    ! Solve 2D Helmholtz in cylindrical coordinates on staggered grid
    call hstcyl(a, b, M, mbdcnd, bda, bdb, c, d, N, nbdcnd, bdc, bdd, &
        elmbda, f, IDIMF, pertrb, ierror)

    ! Compute discretization error by minimizing over all a the function
    ! norm(f(i, j) - a*1 - u(r(i), z(j))). The exact solution is
    !
    ! u(r, z) = (r*z)**4 + arbitrary constant.
    block
        real(wp), parameter :: KNOWN_PERTRB = -0.443113920336705e-3_wp
        real(wp), parameter :: KNOWN_ERROR = 0.752796331450795e-4_wp
        real(wp) :: x, discretization_error
        real(wp) :: exact_solution(M,N)

        x = ZERO
        do i = 1, M
            x = x + sum(f(i,:)-(r(i)*z)**4)
        end do
        x = x/(M*N)
        f(:M,:) = f(:M,:) - x

        do j = 1, N
            do i = 1, M
                exact_solution(i,j) = (r(i)*z(j))**4
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:M,:N)))

        call check_output('hstcyl', &
            ierror, KNOWN_ERROR, discretization_error, KNOWN_PERTRB, pertrb)
    end block

end program test_hstcyl
