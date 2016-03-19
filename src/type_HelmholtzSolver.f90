module type_HelmholtzSolver

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT, &
        stdout => OUTPUT_UNIT

    use module_hwscrt, only: &
        hwscrt

    use module_hstcrt, only: &
        hstcrt

    use type_HelmholtzData, only: &
        HelmholtzData

    use type_Grid, only: &
        Grid

    use type_CenteredGrid, only: &
        CenteredGrid

    use type_StaggeredGrid, only: &
        StaggeredGrid

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: HelmholtzSolver

    ! Declare derived data type
    type, extends( HelmholtzData ), public ::  HelmholtzSolver
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Class methods
        !---------------------------------------------------------------------------------
        procedure,         public :: solve_2d_helmholtz_centered
        procedure,         public :: solve_2d_helmholtz_staggered
        procedure, nopass, public :: unit_test => unit_test_helmholtz_solver
        final                     :: finalize_helmholtz_solver
        !---------------------------------------------------------------------------------
    end type HelmholtzSolver


contains


    subroutine solve_2d_helmholtz_centered( &  ! hwscrt
        this, helmholtz_constant, source_term, solution, perturbation, error_flag )
        !
        ! Purpose:
        !
        ! Solves the standard five-point finite difference approximation to the helmholtz
        ! equation in cartesian coordinates. This equation is
        !
        ! (d/dx)(du/dx) + (d/dy)(du/dy)+ lambda * u = f(x, y).
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (HelmholtzSolver), intent (in out) :: this
        real (wp),               intent (in)     :: helmholtz_constant
        real (wp), contiguous,   intent (in out) :: source_term(:,:)
        real (wp), contiguous,   intent (out)    :: solution(:,:)
        real (wp), optional,     intent (out)    :: perturbation
        integer (ip), optional,  intent (out)    :: error_flag
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)  :: local_error_flag
        real (wp)     :: local_perturbation
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if ( .not.this%initialized ) then
            error stop 'uninitialized object in SOLVE_2D_HELMHOLTZ_CENTERED!'
        end if

        ! Call procedural solver
        associate ( &
            a => this%domain%A, &
            b => this%domain%B, &
            m => size( solution, dim = 1) - 1, &
            mbdcnd => this%Y_BOUNDARY_CONDITION_TYPE, &
            bda => this%west, &
            bdb => this%east, &
            c => this%domain%C, &
            d => this%domain%D, &
            n => size( solution, dim = 2) - 1, &
            nbdcnd => this%X_BOUNDARY_CONDITION_TYPE, &
            bdc => this%south, &
            bdd => this%north, &
            elmbda => helmholtz_constant, &
            f => source_term, &
            idimf => size( source_term, dim = 1 ), &
            pertrb => local_perturbation, &
            ierror => local_error_flag &
            )
            call hwscrt( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror )
        end associate

        ! Address the error flag
        if ( local_error_flag /= 0 ) then
            write( stderr, '(A)') 'ERROR: SOLVE_2D_HELMHOLTZ_CENTERED (hwscrt)'
            select case (local_error_flag)
                case(1)
                    write( stderr, '(A)') 'Invalid X_INTERVAL'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'B >= A'
                case(2)
                    write( stderr, '(A)') 'Invalid boundary condition type in x'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') '0 <= X_BOUNDARY_CONDITION_TYPE <= 4'
                case(3)
                    write( stderr, '(A)') 'Invalid Y_INTERVAL'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'D >= C'
                case(4)
                    write( stderr, '(A)') 'Insufficient number of vertically staggered grid points'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'NY >= 3'
                case(5)
                    write( stderr, '(A)') 'Invalid boundary condition type in y'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') '0 <= Y_BOUNDARY_CONDITION_TYPE <= 4'
                case(6)
                    write( stderr, '(A)') 'Invalid helmholtz constant'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'HELMHOLTZ_CONSTANT <= 0'
                case(7)
                    write( stderr, '(A)') 'Invalid rank for SOURCE_TERM'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'size( SOURCE_TERM, dim = 1) >= NX + 1'
                case(8)
                    write( stderr, '(A)') 'Insufficient number of horizontally staggered grid points'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'NX >= 3'
                case(20)
                    write( stderr, '(A)') 'The dynamic allocation of real and or complex WORKSPACE failed'
                    write( stderr, '(A)') 'NX or NY may be too large for your computer'
                case default
                    write( stderr, '(A)') 'Undetermined error flag'
            end select
        end if

        ! set solution
        associate ( &
            nx => size( solution, dim = 1) - 1, &
            ny => size( solution, dim = 2) - 1 &
            )
            solution = source_term( 1:nx + 1, 1:ny + 1 )
        end associate

        ! Check optional argument
        if ( present (perturbation) ) then
            perturbation = local_perturbation
        end if

    end subroutine solve_2d_helmholtz_centered


    subroutine solve_2d_helmholtz_staggered( & ! hstcrt
        this, helmholtz_constant, source_term, solution, perturbation, error_flag )
        !
        ! Purpose:
        !
        ! Solves the standard five-point finite difference
        ! approximation to the helmholtz equation
        !
        ! (d/dx)(du/dx) + (d/dy)(du/dy) + lambda * u = f(x, y)
        !
        !  on a staggered grid in cartesian coordinates.
        !
        !--------------------------------------------------------------------------------
        class (HelmholtzSolver), intent (in out) :: this
        real (wp),               intent (in)     :: helmholtz_constant
        real (wp), contiguous,   intent (in out) :: source_term(:,:)
        real (wp), contiguous,   intent (out)    :: solution(:,:)
        real (wp), optional,     intent (out)    :: perturbation
        integer (ip), optional,  intent (out)    :: error_flag
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: local_error_flag
        real (wp)    :: local_perturbation
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if ( .not.this%initialized ) then
            error stop 'uninitialized object in SOLVE_2D_HELMHOLTZ_STAGGERED!'
        end if

        ! Call procedural solver
        associate ( &
            a => this%domain%A, &
            b => this%domain%B, &
            m => size( this%south ), &
            mbdcnd => this%X_BOUNDARY_CONDITION_TYPE, &
            bda => this%west, &
            bdb => this%east, &
            c => this%domain%C, &
            d => this%domain%D, &
            n => size( this%west ), &
            nbdcnd => this%Y_BOUNDARY_CONDITION_TYPE, &
            bdc => this%south, &
            bdd => this%north, &
            elmbda => helmholtz_constant, &
            f => source_term, &
            idimf => size( source_term, dim = 1 ), &
            pertrb => local_perturbation, &
            ierror => local_error_flag &
            )
            call hstcrt( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror )
        end associate

        ! Address the error flag
        if ( local_error_flag /= 0 ) then
            write( stderr, '(A)') 'ERROR: SOLVE_2D_HELMHOLTZ_STAGGERED (hstcrt)'
            select case (local_error_flag)
                case(1)
                    write( stderr, '(A)') 'Invalid X_INTERVAL'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'B >= A'
                case(2)
                    write( stderr, '(A)') 'Invalid boundary condition type in x'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') '0 <= X_BOUNDARY_CONDITION_TYPE <= 4'
                case(3)
                    write( stderr, '(A)') 'Invalid Y_INTERVAL'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'D >= C'
                case(4)
                    write( stderr, '(A)') 'Insufficient number of vertically staggered grid points'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'NY >= 2'
                case(5)
                    write( stderr, '(A)') 'Invalid boundary condition type in y'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') '0 <= Y_BOUNDARY_CONDITION_TYPE <= 4'
                case(6)
                    write( stderr, '(A)') 'Invalid helmholtz constant'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'HELMHOLTZ_CONSTANT <= 0'
                case(7)
                    write( stderr, '(A)') 'Invalid rank for SOURCE_TERM'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'size( SOURCE_TERM, dim = 1) >= NX + 1'
                case(8)
                    write( stderr, '(A)') 'Insufficient number of horizontally staggered grid points'
                    write( stderr, '(A)') 'Fails to satisfy:'
                    write( stderr, '(A)') 'NX >= 2'
                case(20)
                    write( stderr, '(A)') 'The dynamic allocation of real and or complex WORKSPACE failed'
                    write( stderr, '(A)') 'NX or NY may be too large for your computer'
                case default
                    write( stderr, '(A)') 'Undetermined error flag'
            end select
        end if

        ! set solution
        associate ( &
            nx => size( solution, dim = 1), &
            ny => size( solution, dim = 2) &
            )
            solution = source_term( 1:nx, 1:ny )
        end associate

        ! Check optional arguments
        if ( present (perturbation) ) then
            perturbation = local_perturbation
        end if

        if ( present (error_flag) ) then
            error_flag = local_error_flag
        end if

    end subroutine solve_2d_helmholtz_staggered


    subroutine finalize_helmholtz_solver( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (HelmholtzSolver), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_helmholtz_solver


    subroutine test_solve_2d_helmholtz_centered() ! Thwscrt
        !
        !< Purpose:
        !
        !     To illustrates the usage of TYPE (HelmholtzSolver) type-bound procedure'
        !
        !     SOLVE_2D_HELMHOLTZ_CENTERED
        !
        !     to solve the equation
        !
        !     (d/dx)(du/dx) + (d/dy)(du/dy) - 4*u
        !
        !     = (2 - (4 + pi**2/4) * x**2 ) * cos((y+1)*pi/2)
        !
        !    with the boundary conditions on the rectangle
        !
        !    (x,y) in [0,2]x[-1,3] with
        !
        !    MIXED ROBIN CONDITION:
        !
        !    u(0,y)       = 0                 for -1 <= y <= 3
        !
        !    (du/dx)(2,y) = 4*cos((y+1)*pi/2) for  -1 <= y <= 3
        !
        !    and with u PERIODIC in y.
        !
        !    the x-interval will be divided into 40 panels and the
        !    y-interval will be divided into 80 panels.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        type (HelmholtzSolver)  :: helmholtz_solver
        type (CenteredGrid)     :: centered_grid
        integer (ip)            :: error_flag
        integer (ip)            :: i, j    !! Counters
        integer (ip), parameter :: NX = 40 !! Number of horizontally staggered grid poins
        integer (ip), parameter :: NY = 80 !! Number of vertically staggered grid points
        real (wp),    parameter :: PI = acos( -1.0_wp )
        real (wp)               :: approximate_solution( NX + 1, NY + 1 )
        real (wp)               :: source( NX + 5, NY + 2 )
        real (wp)               :: discretization_error
        real (wp)               :: exact_solution
        !--------------------------------------------------------------------------------

        ! Print description to console
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     To illustrates the usage of TYPE (HelmholtzSolver) type-bound procedure'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     SOLVE_2D_HELMHOLTZ_CENTERED'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     to solve the equation'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     (d/dx)(du/dx) + (d/dy)(du/dy) - 4*u'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     = (2 - (4 + pi**2/4) * x**2 ) * cos((y+1)*pi/2)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    with the boundary conditions on the rectangle'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    (x,y) in [0,2]x[-1,3] with'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    MIXED ROBIN CONDITION:'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    u(0,y)       = 0                 for -1 <= y <= 3'
        write( stdout, '(A)' ) '    (du/dx)(2,y) = 4*cos((y+1)*pi/2) for -1 <= y <= 3'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    and with u PERIODIC in y.'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    the x-interval will be divided into 40 panels and the'
        write( stdout, '(A)' ) '    y-interval will be divided into 80 panels.'
        write( stdout, '(A)' ) ''

        ! Create centered grid
        associate( &
            x_interval => [  0.0_wp, 2.0_wp ], &
            y_interval => [ -1.0_wp, 3.0_wp ] &
            )
            call centered_grid%create( x_interval, y_interval, NX, NY )
        end associate


        ! Associate various quantities
        associate( &
            x => centered_grid%x, &
            y => centered_grid%y, &
            domain => centered_grid%domain, &
            ROBIN_CONDITIONS_IN_X => 2, &
            PERIODIC_IN_Y => 0, &
            HELMHOLTZ_CONSTANT => -4.0_wp, &
            f => source, &
            u => approximate_solution, &
            HALF_PI => PI/2, &
            PI_SQUARED => PI**2 &
            )


            ! Set the boundary conditions
            call helmholtz_solver%create( NX + 1, NY + 1, PERIODIC_IN_Y, ROBIN_CONDITIONS_IN_X, domain, func )


            ! generate boundary data
            call helmholtz_solver%assign_boundary_data( centered_grid )


            ! Set source, i.e., right side of helmholtz equation
            f(1, :NY + 1) = 0.0_wp
                        do j = 1, NY + 1
                do i = 2, NX + 1
                    f(i, j) = ( 2.0_wp  - ( 4.0_wp + PI_SQUARED/4.0_wp ) * ( x(i)**2 )) * cos( ( y(j) + 1.0_wp) * HALF_PI )
                end do
            end do


            ! Solve helmholtz equation
            call helmholtz_solver%solve_2d_helmholtz_centered( HELMHOLTZ_CONSTANT, f, u, error_flag = error_flag )


            ! Compute discretization error
            discretization_error = 0.0_wp
            do j = 1, NY + 1
                do i = 1, NX + 1
                    exact_solution = ( x(i)**2 ) * cos( ( y(j)+1.0_wp ) * HALF_PI )
                    associate( local_error => abs( u(i, j) - exact_solution ))
                        discretization_error = max( discretization_error, local_error)
                    end associate
                end do
            end do
        end associate

        !     Print earlier output from platforms with 64 bit floating point
        !     arithmetic followed by the output from this computer
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    *** TEST_SOLVE_2D_HELMHOLTZ_CENTERED ***    '
        write( stdout, '(A)' ) '    Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)' ) '    ierror = 0,  discretization error = 5.3650e-4'
        write( stdout, '(A)' ) '    Output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)') &
            '    ierror = ', error_flag, '    discretization error = ', discretization_error
        write( stdout, '(A)' ) ''


    contains


        subroutine func( solver, centered_grid )
            !
            ! Purpose:
            !
            ! User-supplied subroutine to assign boundary data
            !
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            class (HelmholtzData), intent (in out)  :: solver
            class (Grid),          intent (in out)  :: centered_grid
            !--------------------------------------------------------------------------------

            select type (solver)
                class is (HelmholtzSolver)
                select type (centered_grid)
                    class is (CenteredGrid)
                    associate( bdb => solver%east )
                        associate( &
                            y => centered_grid%y(1:size(bdb)), &
                            HALF_PI => acos(-1.0_wp)/2 &
                            )
                            bdb = 4.0_wp * cos( ( y + 1.0_wp ) * HALF_PI)
                            ! Remark: The west (bda), south(bdc), and north (bdd) components are dummy variables.
                        end associate
                    end associate
                end select
            end select

        end subroutine func


    end subroutine test_solve_2d_helmholtz_centered


    subroutine test_solve_2d_helmholtz_staggered() ! Thstcrt
        !
        !< Purpose:
        !
        !     To illustrate the usage of TYPE(HelmholtzSolver) type-bound procedure
        !
        !     SOLVE_2D_HELMHOLTZ_STAGGERED
        !
        !     to solve the equation
        !
        !    (d/dx)(du/dx) + (d/dy)(du/dy) - 2*u = -2(pi**2+1)sin(pi*x)cos(pi*y)
        !
        !    with the boundary conditions on the rectangle
        !
        !    (x,y) in [1,3]x[-1,1] with
        !
        !    MIXED ROBIN CONDITION:
        !
        !     u(1,y)       = 0,                for -1 <= y <= 3
        !     (du/dx)(3,y) = -pi*cos(pi*y),    for -1 <= y <= 3
        !
        !     and u is PERIODIC in y.
        !
        !     we want to have 48 unknowns in the x-interval and 53 unknowns
        !     in the y-interval.
        !
        !     The exact solution to this problem is
        !
        !                  u(x,y) = sin(pi*x) * cos(pi*y).
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (HelmholtzSolver)  :: helmholtz_solver
        type (StaggeredGrid)    :: staggered_grid
        integer (ip)            :: error_flag
        integer (ip)            :: i, j     !! Counters
        integer (ip), parameter :: NX = 48  !! Number of horizontally staggered grid poins
        integer (ip), parameter :: NY = 53  !! Number of vertically staggered grid points
        real (wp),    parameter :: PI = acos( -1.0_wp )
        real (wp)               :: approximate_solution( NX, NY )
        real (wp)               :: source( NX + 2, NY )
        real (wp)               :: discretization_error
        real (wp)               :: exact_solution
        !--------------------------------------------------------------------------------


        ! Print description to console
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '   Illustrates the usage of TYPE (HelmholtzSolver) type-bound procedure'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '   SOLVE_2D_HELMHOLTZ_STAGGERED'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '   to solve the equation'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    (d/dx)(du/dx) + (d/dy)(du/dy) - 2*u = -2(pi**2+1)sin(pi*x)cos(pi*y)'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    with the boundary conditions on the rectangle'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    (x,y) in [1,3]x[-1,1] with'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    MIXED ROBIN CONDITION:'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     u(1,y)       = 0,                for -1 <= y <= 3'
        write( stdout, '(A)' ) '     (du/dx)(3,y) = -pi*cos(pi*y),    for -1 <= y <= 3'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     and u is PERIODIC in y.'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     we want to have 48 unknowns in the x-interval and 53 unknowns'
        write( stdout, '(A)' ) '     in the y-interval.'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     The exact solution to this problem is'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '         u(x,y) = sin(pi*x) * cos(pi*y).'
        write( stdout, '(A)' ) ''


        ! Create staggered grid
        associate( &
            x_interval => [  1.0_wp, 3.0_wp ], &
            y_interval => [ -1.0_wp, 1.0_wp ] &
            )
            call staggered_grid%create( x_interval, y_interval,  NX, NY )
        end associate


        associate( &
            x => staggered_grid%x, &
            y => staggered_grid%y, &
            domain => staggered_grid%domain, &
            f => source, &
            u => approximate_solution, &
            MIXED_ROBIN_IN_X => 2, &
            PERIODIC_IN_Y => 0, &
            HELMHOLTZ_CONSTANT => -2.0_wp &
            )

            ! Set boundary conditions
            call helmholtz_solver%create( NX, NY, MIXED_ROBIN_IN_X, PERIODIC_IN_Y, domain, func )

            ! generate boundary data
            call helmholtz_solver%assign_boundary_data( staggered_grid )


            ! Generate the source, i.e., right side of equation
            associate( C => -2.0_wp * ( (PI **2 ) + 1.0_wp) )
                do i = 1, NX
                    do j = 1, NY
                        f(i, j) = C * sin( PI * x(i) ) * cos( PI * y(j) )
                    end do
                end do
            end associate


            ! Solve helmholtz equation
            call helmholtz_solver%solve_2d_helmholtz_staggered( HELMHOLTZ_CONSTANT, f, u, error_flag = error_flag )


            ! Compute discretization error
            discretization_error = 0.0_wp
            do i = 1, NX
                do j = 1, NY
                    exact_solution = sin( PI * x(i) ) * cos( PI * y(j) )
                    associate( local_error => abs( u(i, j) - exact_solution ))
                        discretization_error = max( discretization_error, local_error)
                    end associate
                end do
            end do

        end associate

        !     Print earlier output from platforms with 64 bit floating point
        !     arithmetic followed by the output from this computer
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '    *** TEST_SOLVE_2D_HELMHOLTZ_STAGGERED ***    '
        write( stdout, '(A)' ) '    Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)' ) '    ierror = 0,  discretization error = 1.2600e-3'
        write( stdout, '(A)' ) '    The output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)' ) &
            '    ierror = ', error_flag, '    discretization error = ', discretization_error
        write( stdout, '(A)' ) ''


    contains


        subroutine func( solver, staggered_grid )
            !
            ! Purpose:
            !
            ! User-supplied subroutine to assign boundary data
            !
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            class (HelmholtzData), intent (in out)  :: solver
            class (Grid),          intent (in out)  :: staggered_grid
            !--------------------------------------------------------------------------------

            select type (solver)
                class is (HelmholtzSolver)
                select type (staggered_grid)
                    class is (StaggeredGrid)
                    associate( &
                        bda => solver%west, &
                        bdb => solver%east &
                        )
                        associate( &
                            y => staggered_grid%y( 1:size(bdb) ), &
                            pi => acos(-1.0_wp) &
                            )
                            bda = 0.0_wp
                            bdb = -pi * cos( pi * y )
                            ! Remark: the south (bdc) and north (bdd) components are dummy variables.
                        end associate
                    end associate
                end select
            end select

        end subroutine func


    end subroutine test_solve_2d_helmholtz_staggered


    subroutine unit_test_helmholtz_solver()
        !
        ! Purpose:
        !
        ! *** UNIT TESTS ***
        !
        !--------------------------------------------------------------------------------

        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '************************************************************************'
        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) '     *** TYPE (HelmholtzSovler) UNIT TESTS *** '
        write( stdout, '(A)' ) ''

        call test_solve_2d_helmholtz_centered() ! compare with test_hwscrt
        call test_solve_2d_helmholtz_staggered() ! compare with test_hstcrt

    end subroutine unit_test_helmholtz_solver


end module type_HelmholtzSolver
