!
!     file tcmgnbn.f
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
!     b(i) = -2(1+x(i))**2*s - sqrt(-1)*dy**2
!
!     c(i) = (1+x(i))**2*s - (1+x(i))*s*dx
!
!     f(i, j) = (3 - sqrt(-1))*(1+x(i))**4*dy**2*sin(y(j))
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
!     b(1) = -2s - sqrt(-1)*dy**2 , c(1) = 2s
!
!     f(1, j) = (11-sqrt(-1)+8/dx)*dy**2*sin(y(j)),  j=1, 2, ..., 40.
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
!     a(20) = (1+x(20))**2*s + (1+x(20))*s*dx
!
!     b(20) = -2*(1+x(20))**2*s - sqrt(-1)*dy**2
!
!     f(20, j) = ((3-sqrt(-1))*(1+x(20))**4*dy**2 - 16(1+x(20))**2*s
!                + 16(1+x(20))*s*dx)*sin(y(j))
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
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        FishpackSolver

    ! Explicit typing only
    implicit None

    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    type (FishpackSolver)    :: solver
    integer (ip) :: idimf, m, mp1, mperod, n, nperod, i, j, ierror
    real (wp), dimension(21) :: x
    real (wp), dimension(41) :: y
    real (wp) :: dx, pi, dy, s, t, tsq, t4, discretization_error
    complex (wp), dimension(22, 40) :: f
    complex, dimension(20) :: a, b, c
    !-----------------------------------------------

    idimf = 22
    m = 20
    mp1 = m + 1
    mperod = 1
    dx = 0.05
    n = 40
    nperod = 0
    pi = acos( -1.0 )
    dy = pi/20.
    !
    !     generate grid points for later use.
    !
    do i = 1, mp1
        x(i) = real(i - 1)*dx
    end do
    do j = 1, n
        y(j) = (-pi) + real(j - 1)*dy
    end do
    !
    !     generate coefficients.
    !
    s = (dy/dx)**2
    do i = 2, 19
        t = 1. + x(i)
        tsq = t**2
        a(i) = cmplx((tsq + t*dx)*s, 0.)
        b(i) = (-2.0_wp * tsq*s) - (0., 1.)*dy**2
        c(i) = cmplx((tsq - t*dx)*s, 0.)
    end do
    a(1) = (0., 0.)
    b(1) = (-2.0_wp * s) - (0., 1.)*dy**2
    c(1) = cmplx(2.0_wp * s, 0.)
    b(20) = (-2.0_wp * s*(1. + x(20))**2) - (0., 1.)*dy**2
    a(20) = cmplx(s*(1. + x(20))**2+(1.+x(20))*dx*s, 0.)
    c(20) = (0., 0.)
    !
    !     generate right side.
    !
    do i = 2, 19
        do j = 1, n
            f(i, j) = (3., -1.)*(1. + x(i))**4*dy**2*sin(y(j))
        end do
    end do
    t = 1. + x(20)
    tsq = t**2
    t4 = tsq**2
    do j = 1, n
        f(1, j) = ((11., -1.) + 8./dx)*dy**2*sin(y(j))
        f(20, j)=((3., -1.)*t4*dy**2-16.*tsq*s+16.*t*s*dx)*sin(y(j))
    end do

    ! Solve system
    call solver%cmgnbn(nperod, n, mperod, m, a, b, c, idimf, f, ierror)

    !     compute discretization error.  the exact solution is
    !
    !            u(x, y) = (1+x)**4*sin(y) .
    !
    discretization_error = 0.
    do i = 1, m
        do j = 1, n
            t = abs(f(i, j)-(1.+x(i))**4*sin(y(j)))
            discretization_error = max(t, discretization_error)
        end do
    end do

    !     Print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     cmgnbn *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  discretization error = 9.1620e-3'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') &
        '     ierror =', ierror, ' discretization error = ', &
        discretization_error

end program tcmgnbn
