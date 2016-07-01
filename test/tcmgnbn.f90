!
!     file tcmgnbn.f90
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
!
! Purpose
!
!     program to illustrate the use of subroutine cmgnbn to solve
!     the equation
!
!     (1+x)**2*(d/dx)(du/dx) - 2(1+x)(du/dx) + (d/dy)(du/dy)
!
!             - sqrt(-1)*u = (3 - sqrt(-1))*(1+x)**4*sin(y)         (1)
!
!     on the rectangle 0 < x < 1 and -pi < y < pi
!     with the boundary conditions
!
!     (du/dx)(0, y) = 4sin(y)                               (2)
!                                -pi <= y <= pi
!     u(1, y) = 16sin(y)                                    (3)
!
!     and with u periodic in y using finite differences on a
!     grid with deltax (= dx) = 1/20 and deltay (= dy) = pi/20.
!        to set up the finite difference equations we define
!     the grid points
!
!     x(i) = (i-1)dx            i=1, 2, ..., 21
!
!     y(j) = -pi + (j-1)dy      j=1, 2, ..., 41
!
!     and let v(i, j) be an approximation to u(x(i), y(j)).
!     numbering the grid points in this fashion gives the set
!     of unknowns as v(i, j) for i=1, 2, ..., 20 and j=1, 2, ..., 40.
!     hence, in the program m = 20 and n = 40.  at the interior
!     grid point (x(i), y(j)), we replace all derivatives in
!     equation (1) by second order central finite differences,
!     multiply by dy**2, and collect coefficients of v(i, j) to
!     get the finite difference equation
!
!     a(i)v(i-1, j) + b(i)v(i, j) + c(i)v(i+1, j)
!
!     + v(i, j-1) - 2v(i, j) + v(i, j+1) = f(i, j)            (4)
!
!     where s = (dy/dx)**2, and for i=2, 3, ..., 19
!
!     a(i) = (1+x(i))**2*s + (1+x(i))*s*dx
!
!     b(i) = -2(1+x(i))**2*s - sqrt(-1)*(dy**2)
!
!     c(i) = (1+x(i))**2*s - (1+x(i))*s*dx
!
!     f(i, j) = (3 - sqrt(-1))*(1+x(i))**4*(dy**2)*sin(y(j))
!              for j=1, 2, ..., 40.
!
!        to obtain equations for i = 1, we replace the
!     derivative in equation (2) by a second order central
!     finite difference approximation, use this equation to
!     eliminate the virtual unknown v(0, j) in equation (4)
!     and arrive at the equation
!
!     b(1)v(1, j) + c(1)v(2, j) + v(1, j-1) - 2v(1, j) + v(1, j+1)
!
!                       = f(1, j)
!
!     where
!
!     b(1) = -2s - sqrt(-1)*(dy**2) , c(1) = 2s
!
!     f(1, j) = (11-sqrt(-1)+8/dx)*(dy**2)*sin(y(j)),  j=1, 2, ..., 40.
!
!     for completeness, we set a(1) = 0.
!        to obtain equations for i = 20, we incorporate
!     equation (3) into equation (4) by setting
!
!     v(21, j) = 16sin(y(j))
!
!     and arrive at the equation
!
!     a(20)v(19, j) + b(20)v(20, j)
!
!     + v(20, j-1) - 2v(20, j) + v(20, j+1) = f(20, j)
!
!     where
!
!     a(20) = (1+x(m))**2*s + (1+x(m))*s*dx
!
!     b(20) = -2*(1+x(m))**2*s - sqrt(-1)*(dy**2)
!
!     f(20, j) = ((3-sqrt(-1))*(1+x(m))**4*(dy**2) - 16(1+x(m))**2*s
!                + 16(1+x(m))*s*dx)*sin(y(j))
!
!                    for j=1, 2, ..., 40.
!
!     for completeness, we set c(20) = 0.  hence, in the
!     program mperod = 1.
!        the periodicity condition on u gives the conditions
!
!     v(i, 0) = v(i, 40) and v(i, 41) = v(i, 1) for i=1, 2, ..., 20.
!
!     hence, in the program nperod = 0.
!
!          the exact solution to this problem is
!
!                  u(x, y) = (1+x)**4*sin(y) .
!
!
!     from the dimension statement we get that idimf = 22
!
program tcmgnbn

    use, intrinsic :: iso_fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        PI, &
        FishpackSolver, &
        FishpackGrid

    ! Explicit typing only
    implicit None

    !------------------------------------------------------------------
    ! Dictionary
    !------------------------------------------------------------------
    type (FishpackSolver)     :: solver
    type (FishpackGrid)       :: grid
    integer (ip), parameter   :: M = 20
    integer (ip), parameter   :: N = 40
    integer (ip), parameter   :: IDIMF = M + 2
    integer (ip)              :: mp1, mperod, nperod, k, j, ierror
    real (wp), allocatable    :: x(:), y(:)
    real (wp)                 :: discretization_error
    complex (wp)              :: f(IDIMF, N)
    complex (wp), allocatable :: a(:), b(:), c(:)
    !------------------------------------------------------------------

    mp1 = M + 1
    mperod = 1
    nperod = 0
    !
    !==> Generate grid points for later use.
    !
    x = grid%get_centered_grid(start=0.0_wp, stop=1.0_wp, num=M)
    y = grid%get_centered_grid(start=-PI, stop=PI, num=N)

    !
    !==> Generate coefficients.
    !
    call get_coefficients(x, y, a, b, c)

    !
    !==> Generate right side.
    !
    call get_right_hand_side(x, y, f)

    !
    !==> Solve system
    !
    call solver%cmgnbn(nperod, N, mperod, M, a, b, c, IDIMF, f, ierror)

    discretization_error = 0.0_wp
    do k = 1, M
        do j = 1, N
            associate( local_error =>  abs(f(k, j)-(1.0_wp+x(k))**4*sin(y(j))) )
                !
                !==> Compute discretization error. The exact solution is
                !
                !            u(x, y) = (1+x)**4*sin(y).
                !
                discretization_error = max(local_error, discretization_error)

            end associate
        end do
    end do

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     cmgnbn *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 9.1620e-3'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        ' discretization error = ', discretization_error

    !
    !==> Release memory
    !
    deallocate( x, y, a, b, c )


contains


    pure subroutine get_coefficients(x, y, a, b, c)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        real (wp),                 intent (in)  :: x(:)
        real (wp),                 intent (in)  :: y(:)
        complex (wp), allocatable, intent (out) :: a(:)
        complex (wp), allocatable, intent (out) :: b(:)
        complex (wp), allocatable, intent (out) :: c(:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip)            :: m, k
        real (wp)               :: dx, dy, s, t, t2
        complex (wp), parameter :: I = (0.0_wp, 1.0_wp)
        !-----------------------------------------------

        m = size(x) - 1
        !
        !==> Allocate memory
        !
        allocate( a(m), b(m), c(m) )

        dx = x(2)-x(1)
        dy = y(2)-y(1)
        s = (dy/dx)**2
        do k = 2, m - 1
            t = 1.0_wp + x(k)
            t2 = t**2
            a(k) = cmplx((t2 + t*dx)*s, 0.0_wp, kind=wp)
            b(k) = -2.0_wp * t2*s - I*(dy**2)
            c(k) = cmplx((t2 - t*dx)*s, 0.0_wp, kind=wp)
        end do
        a(1) = 0.0_wp
        b(1) = -2.0_wp * s - I*(dy**2)
        c(1) = cmplx(2.0_wp * s, 0.0_wp, kind=wp)
        b(m) = (-2.0_wp * s*(1.0_wp + x(m))**2) - I*(dy**2)
        a(m) = cmplx(s*(1.0_wp + x(m))**2 + (1.0_wp + x(m))*dx*s, 0.0_wp, kind=wp)
        c(m) = 0.0_wp

    end subroutine get_coefficients



    pure subroutine get_right_hand_side(x, y, f)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        real (wp),    intent (in)  :: x(:)
        real (wp),    intent (in)  :: y(:)
        complex (wp), intent (out) :: f(:,:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: n, m, k, j
        real (wp)    :: dx, dy, s, t, t2, t4
        !-----------------------------------------------

        m = size(x) - 1
        n = size(y) - 1
        dx = x(2)-x(1)
        dy = y(2)-y(1)
        s = (dy/dx)**2

        do k = 2, m - 1
            do j = 1, n
                f(k, j) = (3.0_wp, -1.0_wp)*(1.0_wp + x(k))**4*(dy**2)*sin(y(j))
            end do
        end do

        s = (dy/dx)**2
        t = 1.0_wp + x(m)
        t2 = t**2
        t4 = t**4
        do j = 1, n
            f(1, j) = ((11.0_wp, -1.0_wp) + 8.0_wp/dx)*(dy**2)*sin(y(j))
            f(m, j)=((3.0_wp, -1.0_wp)*t4*(dy**2) &
                -16.0_wp*t2*s+16.0_wp*t*s*dx)*sin(y(j))
        end do

    end subroutine get_right_hand_side

end program tcmgnbn
