module type_TridiagonalSolver

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT, &
        stdout => OUTPUT_UNIT

    use type_TridiagonalData, only: &
        TriDiagonalData

    use type_Grid, only: &
        Grid

    use type_CenteredGrid, only: &
        CenteredGrid

    use type_StaggeredGrid, only: &
        StaggeredGrid

    use module_genbun, only: &
        genbun

    use module_poistg, only: &
        poistg

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: TridiagonalSolver

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    character (len=250)     :: error_message            !! Probably long enough
    integer (ip)            :: allocate_status          !! To check allocation status
    integer (ip)            :: deallocate_status        !! To check deallocation status
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type, extends( TridiagonalData ), public :: TridiagonalSolver
       !---------------------------------------------------------------------------------
       ! Class variables
       !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Class methods
        !---------------------------------------------------------------------------------
        procedure,         public :: solve_2d_real_linear_system_staggered
        procedure,         public :: solve_2d_real_linear_system_centered
        procedure, nopass, public :: unit_test => unit_test_tridiagonal_solver
        final                     :: finalize_tridiagonal_solver
        !---------------------------------------------------------------------------------
    end type TridiagonalSolver


contains


    subroutine finalize_tridiagonal_solver( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (TridiagonalSolver), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_tridiagonal_solver


    subroutine solve_2d_real_linear_system_staggered( & ! POISTG
        this, source, solution, error_flag )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (TridiagonalSolver), intent (in out) :: this
        real (wp), contiguous,     intent (in out) :: source(:,:)
        real (wp), contiguous,     intent (out)    :: solution(:,:)
        integer (ip), optional,    intent (out)    :: error_flag
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)   :: local_error_flag
        !--------------------------------------------------------------------------------

        ! Invoke procedural solver
        associate( &
            nperod => this%Y_BOUNDARY_CONDITION, &
            n => size( source, dim = 2 ), &
            mperod => this%X_BOUNDARY_CONDITION, &
            m => size( this%subdiagonal ), &
            a => this%subdiagonal, &
            b => this%diagonal, &
            c => this%superdiagonal, &
            idimy  => size( source, dim = 1 ), &
            y => source, &
            ierror => local_error_flag &
            )
            call poistg( nperod, n, mperod, m, a, b, c, idimy, y, ierror )
        end associate

        ! Address the error flag
        if ( local_error_flag /= 0 ) then
            write( stderr, '(A)') 'ERROR: 2D_REAL_LINEAR_SYSTEM_STAGGERED (POISTG)'
            select case (local_error_flag)
                case(1)
                    write( stderr, '(A)') 'Insufficient number of horizontally staggered grid points'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'NX >= 2'
                case(2)
                    write( stderr, '(A)') 'Insufficient number of vertically staggered grid points'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'NY >= 2'
                case(3)
                    write( stderr, '(A)') 'Invalid rank for SOURCE_TERM'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'size( SOURCE_TERM, dim = 1) >= NX'
                case(4)
                    write( stderr, '(A)') 'Invalid boundary condition type in y'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') '1 <= Y_BOUNDARY_CONDITION_TYPE <= 4'
                case(5)
                    write( stderr, '(A)') 'Invalid boundary condition type in x'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') '0 <= X_BOUNDARY_CONDITION_TYPE <= 1'
                case(6)
                    write( stderr, '(A)') 'Invalid tridiagonal system'
                    write( stderr, '(A)') 'X_BOUNDARY_TYPE = 0 and A(1) /= C(1)'
                    write( stderr, '(A)') 'or B(i) /= B(1) or C(i) /= C(1)'
                    write( stderr, '(A)') 'for some i = 1, 2, ..,NX'
                case(7)
                    write( stderr, '(A)') 'Invalid tridiagonal system'
                    write( stderr, '(A)') 'Y_BOUNDARY_TYPE = 1 and '
                    write( stderr, '(A)') 'A(1) /= 0 or C(NX) \= 0 '
                case(20)
                    write( stderr, '(A)') 'The dynamic allocation of real and or complex WORKSPACE failed'
                    write( stderr, '(A)') 'NX or NY may be too large for your computer'
                case default
                    write( stderr, '(A)') 'Undetermined error flag'
            end select
        end if

        ! Set solution
        associate( &
            nx => size( solution, dim = 1), &
            ny => size( solution, dim = 2) &
            )
            solution = source( 1:nx, 1:ny )
        end associate


        ! Address optional arguments
        if ( present (error_flag) ) then
            error_flag = local_error_flag
        end if

    end subroutine solve_2d_real_linear_system_staggered


    subroutine solve_2d_real_linear_system_centered( & ! genbun
        this, source, solution, error_flag )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (TridiagonalSolver), intent (in out) :: this
        real (wp), contiguous,     intent (in out) :: source(:,:)
        real (wp), contiguous,     intent (out)    :: solution(:,:)
        integer (ip), optional,    intent (out)    :: error_flag
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)   :: local_error_flag
        !--------------------------------------------------------------------------------

        ! Invoke procedural solver
        associate( &
            nperod => this%Y_BOUNDARY_CONDITION, &
            n => size( source, dim = 2), &
            mperod => this%X_BOUNDARY_CONDITION, &
            m => size( this%subdiagonal ), &
            a => this%subdiagonal, &
            b => this%diagonal, &
            c => this%superdiagonal, &
            idimy => size( source, dim = 1 ), &
            y => source, &
            ierror => local_error_flag &
            )
            call genbun( nperod, n, mperod, m, a, b, c, idimy, y, ierror )
        end associate

        ! Address the error flag
        if ( local_error_flag /= 0 ) then
            write( stderr, '(A)') 'ERROR: 2D_REAL_LINEAR_SYSTEM_CENTERED (GENBUN)'
            select case (local_error_flag)
                case(1)
                    write( stderr, '(A)') 'Insufficient number of horizontally staggered grid points'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'NX >= 2'
                case(2)
                    write( stderr, '(A)') 'Insufficient number of vertically staggered grid points'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'NY >= 2'
                case(3)
                    write( stderr, '(A)') 'Invalid rank for SOURCE_TERM'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'size( SOURCE_TERM, dim = 1) >= NX'
                case(4)
                    write( stderr, '(A)') 'Invalid boundary condition type in y'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') '0 <= Y_BOUNDARY_CONDITION_TYPE <= 4'
                case(5)
                    write( stderr, '(A)') 'Invalid boundary condition type in x'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') '0 <= X_BOUNDARY_CONDITION_TYPE <= 1'
                case(6)
                    write( stderr, '(A)') 'Invalid tridiagonal system'
                    write( stderr, '(A)') 'A(1) /= C(1) or C(i) /= C(1)'
                    write( stderr, '(A)') 'or B(i) /= B(1) '
                    write( stderr, '(A)') 'for some i = 1, 2, ..,NX'
                case(7)
                    write( stderr, '(A)') 'Invalid tridiagonal system'
                    write( stderr, '(A)') 'X_BOUNDARY_TYPE = 1 and '
                    write( stderr, '(A)') 'A(1) /= 0 or C(NX) \= 0 '
                case(20)
                    write( stderr, '(A)') 'The dynamic allocation of real and or complex WORKSPACE failed'
                    write( stderr, '(A)') 'NX or NY may be too large for your computer'
                case default
                    write( stderr, '(A)') 'Undetermined error flag'
            end select
        end if

        ! set solution
        associate( &
            nx => size( solution, dim = 1 ), &
            ny => size( solution, dim = 2 ) &
            )
            solution = source( 1:nx, 1:ny )
        end associate

        ! Address optional arguments
        if ( present (error_flag) ) then
            error_flag = local_error_flag
        end if

    end subroutine solve_2d_real_linear_system_centered


    subroutine test_solve_2d_real_linear_system_staggered() ! tpoistg
        !
        ! Purpose:
        !
        !     To illustrate the usage of TYPE (TridiagonalSolver) type-bound procedure
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
        !     program mperod = 1.
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
        !     hence, in the program nperod = 2 .
        !
        !     The exact solution to this problem is
        !
        !        u(x, y) = y**4 * cos(x) .
        !
        !
        !     from dimension statement we get that size(f, dim = 1) = 42
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (TridiagonalSolver) :: solver
        type (StaggeredGrid)     :: staggered_grid
        integer (ip), parameter  :: NX = 40 !! number of horizontally staggered grid points
        integer (ip), parameter  :: NY = 20 !! number of vertically staggered grid points
        integer (ip)             :: i, j    !! counter
        integer (ip)             :: error_flag
        real (wp),    parameter  :: PI = acos( -1.0_wp )
        real (wp)                :: approximate_solution( NX, NY )
        real (wp)                :: source( NX, NY )
        real (wp)                :: discretization_error
        real (wp)                :: exact_solution
        !------------------------------------------------------------------------------

        ! Print description to console
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     To illustrate the usage of TYPE (TridiagonalSolver) type-bound procedure'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     SOLVE_2D_REAL_LINEAR_SYSTEM_STAGGERED'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     to solve the equation'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     (1/cos(x))(d/dx)(cos(x)(du/dx)) + (d/dy)(du/dy) ='
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '           2 * y**2 * (6-y**2) * sin(x)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     on the rectangle -pi/2 < x < pi/2 and 0 < y <  1'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     with the boundary conditions'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     (du/dx) (-pi/2, y) = (du/dx)(pi/2, y) = 0 , 0 <= y <= 1  (2)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     u(x, 0) = 0                                           (3)'
        write( stdout, '(A)' ) '                                 -pi/2 <= x <= pi/2'
        write( stdout, '(A)' ) '     (du/dy)(x, 1) = 4sin(x)                               (4)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     using finite differences on a staggered grid with'
        write( stdout, '(A)' ) '     deltax (= dx) = pi/40 and deltay (= dy) = 1/20 .'
        write( stdout, '(A)' ) '        to set up the finite difference equations we define'
        write( stdout, '(A)' ) '     the grid points'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     x(i) = -pi/2 + (i-0.5)dx            i=1, 2, ..., 40'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     y(j) = (j-0.5)dy                    j=1, 2, ..., 20'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     and let v(i, j) be an approximation to u(x(i), y(j)).'
        write( stdout, '(A)' ) '     numbering the grid points in this fashion gives the set'
        write( stdout, '(A)' ) '     of unknowns as v(i, j) for i=1, 2, ..., 40 and j=1, 2, ..., 20.'
        write( stdout, '(A)' ) '     hence, in the program m = 40 and n = 20.  at the interior'
        write( stdout, '(A)' ) '     grid point (x(i), y(j)), we replace all derivatives in'
        write( stdout, '(A)' ) '     equation (1) by second order central finite differences,'
        write( stdout, '(A)' ) '     multiply by dy**2, and collect coefficients of v(i, j) to'
        write( stdout, '(A)' ) '     get the finite difference equation'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     a(i)v(i-1, j) + b(i)v(i, j) + c(i)v(i+1, j)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     + v(i, j-1) - 2v(i, j) + v(i, j+1) = f(i, j)            (5)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     where s = (dy/dx)**2, and for i=2, 3, ..., 39'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     a(i) = s * cos(x(i)-dx/2)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     b(i) = -s * (cos(x(i)-dx/2)+cos(x(i)+dx/2))'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     c(i) = s * cos(x(i)+dx/2)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     f(i, j) = 2dy**2 * y(j)**2 * (6-y(j)**2) * sin(x(i)) , j=1, 2, ..., 19.'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '        to obtain equations for i = 1, we replace equation (2)'
        write( stdout, '(A)' ) '     by the second order approximation'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     (v(1, j)-v(0, j))/dx = 0'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     and use this equation to eliminate v(0, j) in equation (5)'
        write( stdout, '(A)' ) '     to arrive at the equation'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     b(1)v(1, j) + c(1)v(2, j) + v(1, j-1) - 2v(1, j) + v(1, j+1)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '                       = f(1, j)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     where'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     b(1) = -s * (cos(x(1)-dx/2)+cos(x(1)+dx/2))'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     c(1) = -b(1)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     for completeness, we set a(1) = 0.'
        write( stdout, '(A)' ) '        to obtain equations for i = 40, we replace the derivative'
        write( stdout, '(A)' ) '     in equation (2) at x=pi/2 in a similar fashion, use this'
        write( stdout, '(A)' ) '     equation to eliminate the virtual unknown v(41, j) in equation'
        write( stdout, '(A)' ) '     (5) and arrive at the equation'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     a(40)v(39, j) + b(40)v(40, j)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     + v(40, j-1) - 2v(40, j) + v(40, j+1) = f(40, j)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     where'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     a(40) = -b(40) = -s * (cos(x(40)-dx/2)+cos(x(40)+dx/2))'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     for completeness, we set c(40) = 0.  hence, in the'
        write( stdout, '(A)' ) '     program mperod = 1.'
        write( stdout, '(A)' ) '        for j = 1, we replace equation (3) by the second order'
        write( stdout, '(A)' ) '     approximation'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '                (v(i, 0) + v(i, 1))/2 = 0'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     to arrive at the condition'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '                v(i, 0) = -v(i, 1) .'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     for j = 20, we replace equation (4) by the second order'
        write( stdout, '(A)' ) '     approximation'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '                (v(i, 21) - v(i, 20))/dy = 4 * sin(x)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     and combine this equation with equation (5) to arrive at'
        write( stdout, '(A)' ) '     the equation'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     a(i)v(i-1, 20) + b(i)v(i, 20) + c(i)v(i+1, 20)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     + v(i, 19) - 2v(i, 20) + v(i, 21) = f(i, 20)'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     where'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     v(i, 21) = v(i, 20)  and'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     f(i, 20) = 2 * dy**2 * y(j)**2 * (6-y(j)**2) * sin(x(i)) - 4 * dy * sin(x(i))'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     hence, in the program nperod = 2.'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '     The exact solution to this problem is'
        write( stderr, '(A)') ''
        write( stdout, '(A)' ) '        u(x, y) = y**4 * cos(x). '
        write( stderr, '(A)') ''

        ! Create staggered grid
        associate( &
            x_interval => [ -PI/2, PI/2 ], &
            y_interval => [ 0.0_wp, 1.0_wp ] &
            )
            call staggered_grid%create( x_interval, y_interval, NX, NY )
        end associate

        ! Associate various quantities
        associate( &
            dx => staggered_grid%DX, &
            dy => staggered_grid%DY, &
            x => staggered_grid%x, &
            y => staggered_grid%y, &
            f => source, &
            u => approximate_solution, &
            HALF_PI => PI/2 &
            )

            ! create the tridiagonal type
            associate( &
                MIXED_ROBIN_IN_X => 1, &
                PERIODIC_IN_Y => 2 &
                )
                call solver%create( NX, MIXED_ROBIN_IN_X, PERIODIC_IN_Y, proc )
            end associate

            ! Set system coefficients
            call solver%assign_coefficients( staggered_grid )

            ! Generate the source, i.e., right side of equation
            do j = 1, NY
                do i = 1, NX
                    f(i, j) = 2.0_wp * ( dy**2 ) * y(j)**2 * ( 6.0_wp - y(j)**2 ) * sin(x(i))
                end do
            end do

            f(:, NY) = f(:, NY) - 4.0_wp * dy * sin(x(1:NX))

            ! Solve system
            call solver%solve_2d_real_linear_system_staggered( f, u, error_flag )

            ! Compute discretization error
            discretization_error = 0.0_wp
            do j = 1, NY
                do i = 1, NX
                    exact_solution = y(j)**4 * sin(x(i))
                    associate( local_error => abs( u(i, j) - exact_solution ))
                        discretization_error = max( discretization_error, local_error )
                    end associate
                end do
            end do
        end associate

        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    *** TEST_SOLVE_2D_REAL_LINEAR_SYSTEM_STAGGERED ***'
        write( stdout, '(A)' ) '    Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)' ) '    ierror = 0,  discretization error = 5.6417e-4'
        write( stdout, '(A)' ) '    The output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)' ) &
            '    ierror =', error_flag, ' discretization error = ', discretization_error


    contains


        subroutine proc( solver, staggered_grid )
            !
            ! Purpose:
            !
            ! User-supplied subroutine to assign coefficients
            !
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            class (TridiagonalData), intent (in out)  :: solver
            class (Grid),            intent (in out)  :: staggered_grid
            !--------------------------------------------------------------------------------
            integer (ip) :: i ! counter
            !--------------------------------------------------------------------------------

            select type (solver)
                class is (TridiagonalSolver)
                select type (staggered_grid)
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
                            pi => acos( -1.0_wp ) &
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

        end subroutine proc

    end subroutine test_solve_2d_real_linear_system_staggered


    subroutine test_solve_2d_real_linear_system_centered() ! tgenbun
        !
        !< Purpose:
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
        !     program mperod = 1.
        !        the periodicity condition on u gives the conditions
        !
        !     v(i,0) = v(i,40) and v(i,41) = v(i,1) for i=1,2,...,20.
        !
        !     hence, in the program nperod = 0.
        !
        !     The exact solution to this problem is
        !
        !                  u(x,y) = ((1+x)**4) * sin(y).
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (TridiagonalSolver) :: solver
        type (CenteredGrid)      :: centered_grid
        integer (ip), parameter  :: NX = 20 !! number of horizontally staggered grid points
        integer (ip), parameter  :: NY = 40 !! number of vertically staggered grid points
        integer (ip)             :: i, j !! counter
        integer (ip)             :: error_flag
        real (wp),    parameter  :: PI = acos( -1.0_wp )
        real (wp)                :: approximate_solution( NX, NY )
        real (wp)                :: source( NX + 2, NY )            !! source term
        real (wp)                :: discretization_error
        real (wp)                :: exact_solution
        !------------------------------------------------------------------------------

        ! Print description to console
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     To illustrate the usage of TYPE (TridiagonalSolver) type-bound procedure'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     SOLVE_2D_REAL_LINEAR_SYSTEM_CENTERED'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     to solve the equation'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     (1+x)**2*(d/dx)(du/dx) - 2(1+x)(du/dx) + (d/dy)(du/dy)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '                  = 3(1+x)**4*sin(y)                      (1)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     on the rectangle 0 < x < 1 and -pi < y < pi'
        write( stdout, '(A)' ) '     with the boundary conditions'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     (du/dx)(0,y) = 4sin(y)                               (2)'
        write( stdout, '(A)' ) '                                -pi <= y <= pi'
        write( stdout, '(A)' ) '     u(1,y) = 16sin(y)                                    (3)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     and with u periodic in y using finite differences on a'
        write( stdout, '(A)' ) '     grid with deltax (= dx) = 1/20 and deltay (= dy) = pi/20.'
        write( stdout, '(A)' ) '        to set up the finite difference equations we define'
        write( stdout, '(A)' ) '     the grid points'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     x(i) = (i-1)dx            i=1,2,...,21'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     y(j) = -pi + (j-1)dy      j=1,2,...,41'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     and let v(i,j) be an approximation to u(x(i),y(j)).'
        write( stdout, '(A)' ) '     numbering the grid points in this fashion gives the set'
        write( stdout, '(A)' ) '     of unknowns as v(i,j) for i=1,2,...,20 and j=1,2,...,40.'
        write( stdout, '(A)' ) '     hence, in the program m = 20 and n = 40.  at the interior'
        write( stdout, '(A)' ) '     grid point (x(i),y(j)), we replace all derivatives in'
        write( stdout, '(A)' ) '     equation (1) by second order central finite differences,'
        write( stdout, '(A)' ) '     multiply by dy**2, and collect coefficients of v(i,j) to'
        write( stdout, '(A)' ) '     get the finite difference equation'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     a(i)v(i-1,j) + b(i)v(i,j) + c(i)v(i+1,j)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     + v(i,j-1) - 2v(i,j) + v(i,j+1) = f(i,j)            (4)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     where s = (dy/dx)**2, and for i=2,3,...,19'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     a(i) = (1+x(i))**2*s + (1+x(i))*s*dx'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     b(i) = -2(1+x(i))**2*s'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     c(i) = (1+x(i))**2*s - (1+x(i))*s*dx'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     f(i,j) = 3(1+x(i))**4*dy**2*sin(y(j))  for j=1,2,...,40.'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     to obtain equations for i = 1, we replace the'
        write( stdout, '(A)' ) '     derivative in equation (2) by a second order central'
        write( stdout, '(A)' ) '     finite difference approximation, use this equation to'
        write( stdout, '(A)' ) '     eliminate the virtual unknown v(0,j) in equation (4)'
        write( stdout, '(A)' ) '     and arrive at the equation'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     b(1)v(1,j) + c(1)v(2,j) + v(1,j-1) - 2v(1,j) + v(1,j+1)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '                       = f(1,j)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     where'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     b(1) = -2s , c(1) = 2s'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     f(1,j) = (11+8/dx)*dy**2*sin(y(j)),  j=1,2,...,40.'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     for completeness, we set a(1) = 0.'
        write( stdout, '(A)' ) '        to obtain equations for i = 20, we incorporate'
        write( stdout, '(A)' ) '     equation (3) into equation (4) by setting'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     v(21,j) = 16sin(y(j))'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     and arrive at the equation'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     a(20)v(19,j) + b(20)v(20,j)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     + v(20,j-1) - 2v(20,j) + v(20,j+1) = f(20,j)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     where'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     a(20) = (1+x(20))**2*s + (1+x(20))*s*dx'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     b(20) = -2*(1+x(20))**2*s'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     f(20,j) = (3(1+x(20))**4*dy**2 - 16(1+x(20))**2*s'
        write( stdout, '(A)' ) '                + 16(1+x(20))*s*dx)*sin(y(j))'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '                    for j=1,2,...,40.'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     for completeness, we set c(20) = 0.  hence, in the'
        write( stdout, '(A)' ) '     program mperod = 1.'
        write( stdout, '(A)' ) '        the periodicity condition on u gives the conditions'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     v(i,0) = v(i,40) and v(i,41) = v(i,1) for i=1,2,...,20.'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     hence, in the program nperod = 0.'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     The exact solution to this problem is'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '                  u(x,y) = ((1+x)**4) * sin(y).'
        write( stdout, '(A)' ) ''


        ! Create centered grid
        associate( &
            x_interval => [ 0.0_wp, 1.0_wp ], &
            y_interval => [ -PI, PI ] &
            )
            call centered_grid%create( x_interval, y_interval, NX, NY )
        end associate

        ! Associate various quantities
        associate( &
            dx => centered_grid%DX, &
            dy => centered_grid%DY, &
            x => centered_grid%x, &
            y => centered_grid%y, &
            f => source, &
            u => approximate_solution &
            )

            ! create the real tridiagonal system
            call solver%create( NX, 1, 0, proc )

            ! Assign coefficients
            call solver%assign_coefficients( centered_grid )

            ! Associate mesh ratio
            associate( s => (dy/dx)**2 )

                ! Generate the source, i.e., right hand side of the equation
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

            ! Solve the real tridiagonal system
            call solver%solve_2d_real_linear_system_centered( f, u, error_flag )

            ! Compute the discretization error
            discretization_error = 0.0_wp
            do j = 1, NY
                do i = 1, NX
                    exact_solution = ( (1.0_wp + x(i))**4 ) * sin( y(j) )
                    associate( local_error => abs( u(i, j) - exact_solution ))
                        discretization_error = max( discretization_error, local_error )
                    end associate
                end do
            end do
        end associate

        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    *** TEST_SOLVE_2D_REAL_LINEAR_SYSTEM_CENTERED *** '
        write( stdout, '(A)' ) '    Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)' ) '    ierror = 0,  discretization error = 9.6406e-3'
        write( stdout, '(A)' ) '    The output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)' ) &
            '    ierror =', error_flag, ' discretization error = ', discretization_error


    contains


        subroutine proc( solver, centered_grid )
            !
            ! Purpose:
            !
            ! User-supplied subroutine to assign coefficients
            !
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
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

        end subroutine proc


    end subroutine test_solve_2d_real_linear_system_centered


    subroutine unit_test_tridiagonal_solver()
        !
        ! Purpose:
        !
        ! *** UNIT TESTS ***
        !
        !--------------------------------------------------------------------------------

        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '************************************************************************'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     *** TYPE (TridiagonalSolver) UNIT TESTS *** '
        write( stdout, '(A)' ) ''

        call test_solve_2d_real_linear_system_staggered() ! poistg
        call test_solve_2d_real_linear_system_centered() ! genbun

    end subroutine unit_test_tridiagonal_solver


end module type_TridiagonalSolver
