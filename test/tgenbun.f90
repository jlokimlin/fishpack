!     file tgenbun.f90
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
!  Purpose:
!
!     To illustrate the usage of TYPE (TridiagonalSolver) type-bound procedure
!
!     SOLVE_2D_REAL_LINEAR_SYSTEM_CENTERED
!
!     to solve the equation
!
!     (1+x)**2*(d/dx)(du/dx) - 2(1+x)(du/dx) + (d/dy)(du/dy)
!
!                  = 3(1+x)**4*sin(y)                      (1)
!
!     on the rectangle 0 < x < 1 and -pi < y < pi
!     with the boundary conditions
!
!     (du/dx)(0,y) = 4sin(y)                               (2)
!                                -pi <= y <= pi
!     u(1,y) = 16sin(y)                                    (3)
!
!     and with u periodic in y using finite differences on a
!     grid with deltax (= dx) = 1/20 and deltay (= dy) = pi/20.
!        to set up the finite difference equations we define
!     the grid points
!
!     x(i) = (i-1)dx            i=1,2,...,21
!
!     y(j) = -pi + (j-1)dy      j=1,2,...,41
!
!     and let v(i,j) be an approximation to u(x(i),y(j)).
!     numbering the grid points in this fashion gives the set
!     of unknowns as v(i,j) for i=1,2,...,20 and j=1,2,...,40.
!     hence, in the program m = 20 and n = 40.  at the interior
!     grid point (x(i),y(j)), we replace all derivatives in
!     equation (1) by second order central finite differences,
!     multiply by dy**2, and collect coefficients of v(i,j) to
!     get the finite difference equation
!
!     a(i)v(i-1,j) + b(i)v(i,j) + c(i)v(i+1,j)
!
!     + v(i,j-1) - 2v(i,j) + v(i,j+1) = f(i,j)            (4)
!
!     where s = (dy/dx)**2, and for i=2,3,...,19
!
!     a(i) = (1+x(i))**2*s + (1+x(i))*s*dx
!
!     b(i) = -2(1+x(i))**2*s
!
!     c(i) = (1+x(i))**2*s - (1+x(i))*s*dx
!
!     f(i,j) = 3(1+x(i))**4*dy**2*sin(y(j))  for j=1,2,...,40.
!
!        to obtain equations for i = 1, we replace the
!     derivative in equation (2) by a second order central
!     finite difference approximation, use this equation to
!     eliminate the virtual unknown v(0,j) in equation (4)
!     and arrive at the equation
!
!     b(1)v(1,j) + c(1)v(2,j) + v(1,j-1) - 2v(1,j) + v(1,j+1)
!
!                       = f(1,j)
!
!     where
!
!     b(1) = -2s , c(1) = 2s
!
!     f(1,j) = (11+8/dx)*dy**2*sin(y(j)),  j=1,2,...,40.
!
!     for completeness, we set a(1) = 0.
!        to obtain equations for i = 20, we incorporate
!     equation (3) into equation (4) by setting
!
!     v(21,j) = 16sin(y(j))
!
!     and arrive at the equation
!
!     a(20)v(19,j) + b(20)v(20,j)
!
!     + v(20,j-1) - 2v(20,j) + v(20,j+1) = f(20,j)
!
!     where
!
!     a(20) = (1+x(20))**2*s + (1+x(20))*s*dx
!
!     b(20) = -2*(1+x(20))**2*s
!
!     f(20,j) = (3(1+x(20))**4*dy**2 - 16(1+x(20))**2*s
!                + 16(1+x(20))*s*dx)*sin(y(j))
!
!                    for j=1,2,...,40.
!
!     for completeness, we set c(20) = 0.  hence, in the
!
!     program Y_PERIODICITY = 1.
!
!        the periodicity condition on u gives the conditions
!
!     v(i,0) = v(i,40) and v(i,41) = v(i,1) for i=1,2,...,20.
!
!     hence, in the program X_PERIODICITY = 0.
!
!     The exact solution to this problem is
!
!                  u(x,y) = ((1+x)**4) * sin(y).
!
program tgenbun

    use, intrinsic :: iso_fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        TridiagonalData, &
        TridiagonalSolver, &
        Grid, &
        CenteredGrid

    ! Explicit typing only
    implicit none

    !--------------------------------------------------------------------------------
    ! Dictionary
    !--------------------------------------------------------------------------------
    type (TridiagonalSolver) :: solver
    type (CenteredGrid)      :: centered_grid
    integer (ip), parameter  :: NX = 20 !! number of horizontally staggered grid points
    integer (ip), parameter  :: NY = 40 !! number of vertically staggered grid points
    integer (ip)             :: i, j !! counter
    integer (ip)             :: error_flag
    real (wp),    parameter  :: PI = acos(-1.0_wp)
    real (wp)                :: approximate_solution(NX, NY)
    real (wp)                :: source(NX + 2, NY)
    real (wp)                :: discretization_error
    real (wp)                :: exact_solution
    !------------------------------------------------------------------------------


    associate( &
        x_interval => [ 0.0_wp, 1.0_wp ], &
        y_interval => [ -PI, PI ] &
        )
        !
        !==> Allocate memory
        !
        centered_grid = CenteredGrid(x_interval, y_interval, NX, NY)

    end associate

    !
    !==> Associate various quantities
    !
    associate( &
        dx => centered_grid%DX, &
        dy => centered_grid%DY, &
        x => centered_grid%x, &
        y => centered_grid%y, &
        f => source, &
        u => approximate_solution &
        )

        !
        !==> Initialize tridiagonal solver
        !
        associate( &
            X_PERIODICITY => 1, &
            Y_PERIODICITY => 0 &
            )
            call solver%create(NX, X_PERIODICITY, Y_PERIODICITY, coeff_procedure)
        end associate

        !
        !==> Assign coefficients
        !
        call solver%assign_coefficients(centered_grid)

        !
        !==> Generate the source, i.e., right hand side of the equation
        !
        associate( s => (dy/dx)**2 )

            do j = 1, NY
                do i = 2, NX - 1
                    f(i, j) = 3.0_wp * (1.0_wp + x(i))**4 * (dy**2) * sin( y(j) )
                end do
            end do

            associate( t => 1.0_wp + x(NX) )
                do j = 1, NY
                    f(1, j) = &
                        (11.0_wp + 8.0_wp/dx) * (dy**2) * sin( y(j) )

                    f(NX, j) = &
                        (3.0_wp * (t**4) * (dy**2) &
                        - 16.0_wp * (t**2) * s &
                        + 16.0_wp * t * s * dx) * sin(y(j))
                end do
            end associate
        end associate

        !
        !==> Solve system
        !
        call solver%solve_2d_real_linear_system_centered(f, u, error_flag)

        !
        !==> Compute the discretization error
        !
        discretization_error = 0.0_wp
        do j = 1, NY
            do i = 1, NX
                ! Set exact solution
                exact_solution = ( (1.0_wp + x(i))**4 ) * sin(y(j))
                ! Set local error
                associate( local_error => abs(u(i, j) - exact_solution))
                    ! Set discretization error
                    discretization_error = max(discretization_error, local_error)
                end associate
            end do
        end do
    end associate

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     genbun *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 9.6406e-3'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     error_flag =', error_flag, &
        ' discretization error = ', discretization_error

    !
    !==> Release memory
    !
    call solver%destroy()
    call centered_grid%destroy()

contains


    subroutine coeff_procedure(solver, centered_grid)
        !
        ! Purpose:
        !
        ! User-supplied subroutine to assign coefficients
        !
        !--------------------------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------------------------
        class (TridiagonalData), intent (in out)  :: solver
        class (Grid),            intent (in out)  :: centered_grid
        !--------------------------------------------------------------------------------
        integer (ip) :: i !! Counter
        !--------------------------------------------------------------------------------

        select type (solver)
            class is (TridiagonalSolver)
            select type (centered_grid)
                class is (CenteredGrid)
                associate( &
                    nx => centered_grid%NX, &
                    dx => centered_grid%DX, &
                    dy => centered_grid%DY, &
                    x => centered_grid%x &
                    )
                    associate( &
                        s => (dy/dx)**2, &
                        a => solver%subdiagonal, &
                        b => solver%diagonal, &
                        c => solver%superdiagonal &
                        )
                        do i = 2, nx - 1
                            associate( t => 1.0_wp + x(i) )
                                a(i) = ( (t**2) + t * dx) * s
                                b(i) = -2.0_wp * ( t**2 ) * s
                                c(i) = ( (t**2) - t * dx) * s
                            end associate
                        end do
                        a(1) = 0.0_wp
                        b(1) = -2.0_wp * s
                        c(1) = -b(1)
                        b(nx) = -2.0_wp * s * (1.0_wp + x(nx))**2
                        a(nx) = (-b(nx)/2.0_wp) + (1.0_wp + x(nx)) * dx * s
                        c(nx) = 0.0_wp
                    end associate
                end associate
            end select
        end select

    end subroutine coeff_procedure


end program tgenbun
