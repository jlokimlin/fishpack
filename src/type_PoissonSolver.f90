module type_PoissonSolver

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_HelmholtzSolver, only: &
        HelmholtzSolver

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: PoissonSolver

    ! Declare derived data type
    type, extends (HelmholtzSolver), public :: PoissonSolver
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Class methods
        !---------------------------------------------------------------------------------
        procedure, non_overridable, public :: solve_2d_poisson_centered
        procedure, non_overridable, public :: solve_2d_poisson_staggered
        final                              :: finalize_poisson_solver
        !---------------------------------------------------------------------------------
    end type PoissonSolver


contains


    subroutine solve_2d_poisson_centered(this, source_term, solution)
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (PoissonSolver),  intent (in out) :: this
        real (wp),              intent (in out) :: source_term(:,:)
        real (wp),              intent (out)    :: solution(:,:)
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (this%initialized .eqv. .false.) then
            error stop 'TYPE (PoissonSolver): '&
                //'uninitialized object in SOLVE_2D_POISSON_CENTERED!'
        end if

        ! Solve poisson's equation
        associate( helmholtz_constant => 0.0_wp )
            call this%solve_2d_helmholtz_centered(helmholtz_constant, source_term, solution)
        end associate

    end subroutine solve_2d_poisson_centered



    subroutine solve_2d_poisson_staggered(this, source_term, solution, perturbation, error_flag)
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (PoissonSolver),  intent (in out) :: this
        real (wp),              intent (in out) :: source_term(:,:)
        real (wp),              intent (out)    :: solution(:,:)
        real (wp),    optional, intent (out)    :: perturbation
        integer (ip), optional, intent (out)    :: error_flag
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: error_flag_op
        real (wp)    :: perturbation_op
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (this%initialized .eqv. .false.) then
            error stop 'TYPE (PoissonSolver): '&
                //'uninitialized object in SOLVE_2D_POISSON_STAGGERED!'
        end if

        !
        !==> Solve poisson's equation
        !
        associate( helmholtz_constant => 0.0_wp )
            call this%solve_2d_helmholtz_staggered( &
                helmholtz_constant, source_term, solution, perturbation_op, error_flag_op)
        end associate

        !
        !==> Address optional arguments
        !
        if (present(perturbation)) then
            perturbation = perturbation_op
        end if

        if (present(error_flag)) then
            error_flag = error_flag_op
        end if

    end subroutine solve_2d_poisson_staggered



    subroutine finalize_poisson_solver( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (PoissonSolver), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_poisson_solver



end module type_PoissonSolver
