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
    type, extends (HelmholtzData), public ::  HelmholtzSolver
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Class methods
        !---------------------------------------------------------------------------------
        procedure, public :: solve_2d_helmholtz_centered
        procedure, public :: solve_2d_helmholtz_staggered
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
        real (wp),               intent (in out) :: source_term(:,:)
        real (wp),               intent (out)    :: solution(:,:)
        real (wp),    optional,  intent (out)    :: perturbation
        integer (ip), optional,  intent (out)    :: error_flag
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)  :: error_flag_op
        real (wp)     :: perturbation_op
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (this%initialized .eqv. .false.) then
            error stop 'Uninitialized object of class (HelmholtzSolver) '&
                //'in solve_2d_helmholtz_centered'
        end if

        !
        !==> Call procedural solver
        !
        associate ( &
            a => this%domain%A, &
            b => this%domain%B, &
            m => size(solution, dim=1) - 1, &
            mbdcnd => this%Y_BOUNDARY_CONDITION_TYPE, &
            bda => this%west, &
            bdb => this%east, &
            c => this%domain%C, &
            d => this%domain%D, &
            n => size(solution, dim=2) - 1, &
            nbdcnd => this%X_BOUNDARY_CONDITION_TYPE, &
            bdc => this%south, &
            bdd => this%north, &
            elmbda => helmholtz_constant, &
            f => source_term, &
            idimf => size(source_term, dim=1), &
            pertrb => perturbation_op, &
            ierror => error_flag_op &
            )
            call hwscrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror)
        end associate

        ! Address the error flag
        if ( error_flag_op /= 0 ) then
            write( stderr, '(A)') 'ERROR: SOLVE_2D_HELMHOLTZ_CENTERED (hwscrt)'
            select case (error_flag_op)
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
                    write( stderr, '(A)') 'size(SOURCE_TERM, dim=1) >= NX + 1'
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
            nx => size(solution, dim=1) - 1, &
            ny => size(solution, dim=2) - 1 &
            )
            solution = source_term(1:nx + 1, 1:ny + 1)
        end associate

        !
        !==> Check optional arguments
        !
        if (present(perturbation)) then
            perturbation = perturbation_op
        end if

        if (present(error_flag)) then
            error_flag = error_flag_op
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
        real (wp),               intent (in out) :: source_term(:,:)
        real (wp),               intent (out)    :: solution(:,:)
        real (wp),    optional,  intent (out)    :: perturbation
        integer (ip), optional,  intent (out)    :: error_flag
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: error_flag_op
        real (wp)    :: perturbation_op
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (this%initialized .eqv. .false.) then
            error stop 'uninitialized object of class (HelmholtzSolver) '&
                //'in solve_2d_helmholtz_staggered'
        end if

        !
        !==> Call procedural solver
        !
        associate ( &
            a => this%domain%A, &
            b => this%domain%B, &
            m => size(this%south ), &
            mbdcnd => this%X_BOUNDARY_CONDITION_TYPE, &
            bda => this%west, &
            bdb => this%east, &
            c => this%domain%C, &
            d => this%domain%D, &
            n => size(this%west), &
            nbdcnd => this%Y_BOUNDARY_CONDITION_TYPE, &
            bdc => this%south, &
            bdd => this%north, &
            elmbda => helmholtz_constant, &
            f => source_term, &
            idimf => size(source_term, dim=1), &
            pertrb => perturbation_op, &
            ierror => error_flag_op &
            )
            call hstcrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror)
        end associate

        ! Address the error flag
        if (error_flag_op /= 0) then
            write( stderr, '(A)') 'TYPE (HelmholtzSolver): '&
                //' in SOLVE_2D_HELMHOLTZ_STAGGERED (hstcrt)'
            select case (error_flag_op)
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
                    write( stderr, '(A)') 'size(SOURCE_TERM, dim=1) >= NX + 1'
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
            nx => size(solution, dim=1), &
            ny => size(solution, dim=2) &
            )
            solution = source_term(1:nx, 1:ny)
        end associate

        !
        !==> Check optional arguments
        !
        if (present(perturbation)) then
            perturbation = perturbation_op
        end if

        if (present(error_flag)) then
            error_flag = error_flag_op
        end if

    end subroutine solve_2d_helmholtz_staggered



    subroutine finalize_helmholtz_solver(this)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (HelmholtzSolver), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_helmholtz_solver



end module type_HelmholtzSolver
