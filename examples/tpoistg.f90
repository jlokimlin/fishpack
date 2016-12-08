!
!     file tpoistg.f
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
!     Purpose:
!
!     To illustrate the usage of TYPE(TridiagonalSolver) type-bound procedure
!
!     SOLVE_2D_REAL_LINEAR_SYSTEM_STAGGERED
!
!     to solve the equation
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
!     program Y_PERIODICITY = 1.
!
!        for j = 1, we replace equation (3) by the second order
!     approximation
!
!                (v(i, 0) + v(i, 1))/2 = 0
!
!     to arrive at the condition
!
!                v(i, 0) = -v(i, 1) .
!
!     for j = 20, we replace equation (4) by the second order
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
!     hence, in the program X_PERIODICITY = 2 .
!
!     The exact solution to this problem is
!
!        u(x, y) = y**4 * cos(x) .
!
!
!     from dimension statement we get that size(f, dim=1) = 42
!
program tpoistg

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        TridiagonalData, &
        TridiagonalSolver, &
        Grid, &
        StaggeredGrid

    ! Explicit typing only
    implicit none

    !--------------------------------------------------------------
    ! Dictionary
    !--------------------------------------------------------------
    type(TridiagonalSolver) :: solver
    type(StaggeredGrid)     :: staggered_grid
    integer(ip), parameter  :: NX = 40 !! number of horizontally staggered grid points
    integer(ip), parameter  :: NY = 20 !! number of vertically staggered grid points
    integer(ip)             :: i, j  !! counters
    integer(ip)             :: error_flag
    real(wp),    parameter  :: PI = acos(-1.0_wp)
    real(wp)                :: approximate_solution(NX, NY)
    real(wp)                :: source(NX, NY)
    real(wp)                :: discretization_error
    real(wp)                :: exact_solution
    !------------------------------------------------------------------------------

    associate( &
        x_interval => [ -PI/2, PI/2 ], &
        y_interval => [ 0.0_wp, 1.0_wp ] &
        )
        !
        !==> Allocate memory
        !
        staggered_grid = StaggeredGrid(x_interval, y_interval, NX, NY)
    end associate

    !
    !==> Associate various quantities
    !
    associate( &
        dx => staggered_grid%DX, &
        dy => staggered_grid%DY, &
        x => staggered_grid%x, &
        y => staggered_grid%y, &
        f => source, &
        u => approximate_solution, &
        HALF_PI => PI/2 &
        )

        !
        !==> Initialize the tridiagonal solver
        !
        associate( &
            X_PERIODICITY => 1, &
            Y_PERIODICITY => 2 &
            )
            call solver%create(NX, X_PERIODICITY, Y_PERIODICITY, coeff_procedure)
        end associate

        !
        !==> Set system coefficients
        !
        call solver%assign_coefficients(staggered_grid)

        !
        !==> Generate the source, i.e., right side of equation
        !
        do j = 1, NY
            do i = 1, NX
                f(i, j) = 2.0_wp * ( dy**2 ) * y(j)**2 * ( 6.0_wp - y(j)**2 ) * sin(x(i))
            end do
        end do

        f(:, NY) = f(:, NY) - 4.0_wp * dy * sin(x(1:NX))

        !
        !==> Solve system
        !
        call solver%solve_2d_real_linear_system_staggered(f, u, error_flag)

        !
        !==> Compute discretization error
        !
        discretization_error = 0.0_wp
        do j = 1, NY
            do i = 1, NX
                ! Set exact solution
                exact_solution = y(j)**4 * sin(x(i))
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
    write( stdout, '(/a)') '     poistg *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 5.6417e-4'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)')&
        '     error_flag = ', error_flag, ' discretization error = ', discretization_error

    !
    !==> Release memory
    !
    call solver%destroy()
    call staggered_grid%destroy()

contains


    subroutine coeff_procedure(solver, staggered_grid)
        !
        ! Purpose:
        !
        ! User-supplied subroutine to assign coefficients
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(TridiagonalData), intent(inout)  :: solver
        class(Grid),            intent(inout)  :: staggered_grid
        !--------------------------------------------------------------
        integer(ip) :: i ! counter
        !--------------------------------------------------------------

        select type(solver)
            class is (TridiagonalSolver)
            select type(staggered_grid)
                class is (StaggeredGrid)
                associate( &
                    nx => staggered_grid%NX, &
                    dx => staggered_grid%DX, &
                    dy => staggered_grid%DY &
                    )
                    associate( &
                        s => (dy/dx)**2, &
                        a => solver%subdiagonal, &
                        b => solver%diagonal, &
                        c => solver%superdiagonal, &
                        x => staggered_grid%x, &
                        pi => acos(-1.0_wp) &
                        )
                        a(1) = 0.0_wp
                        b(1) = -s * cos((-pi/2.0_wp) + dx)/cos(x(1))
                        c(1) = -b(1)
                        do i = 2, nx
                            a(i) = s * cos( x(i) - dx/2.0_wp ) / cos(x(i))
                            c(i) = s * cos( x(i) + dx/2.0_wp ) / cos(x(i))
                            b(i) = -( a(i)+c(i) )
                        end do
                        a(nx) = -b(nx)
                        c(nx) = 0.0_wp
                    end associate
                end associate
            end select
        end select

    end subroutine coeff_procedure

end program tpoistg
