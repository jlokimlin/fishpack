module type_PoissonSolver

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

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
        ! Type components
        !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Type-bound procedures
        !---------------------------------------------------------------------------------
        procedure, public :: solve_2d_poisson_centered
        procedure, public :: solve_2d_poisson_staggered
        !final             :: finalize_poisson_solver
        !---------------------------------------------------------------------------------
    end type PoissonSolver



contains



    subroutine solve_2d_poisson_centered(this, source_term, solution)
        !
        !--------------------------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------------------------
        class (PoissonSolver), intent (in out) :: this
        real (wp),             intent (in out) :: source_term(:,:)
        real (wp),             intent (out)    :: solution(:,:)
        !--------------------------------------------------------------------------------
        real (wp), parameter :: HELMHOLTZ_CONSTANT = 0.0_wp
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (this%initialized .eqv. .false.) then
            error stop 'Uninitialized object of class (PoissonSolver) '&
                //'in solve_2d_poisson_centered'
        end if

        !
        !==> Solve poisson's equation on centered grid
        !
        associate( &
            elmbda => HELMHOLTZ_CONSTANT, &
            rhs => source_term, &
            f => solution &
            )
            call this%solve_2d_helmholtz_centered(elmbda, rhs, f)
        end associate

    end subroutine solve_2d_poisson_centered



    subroutine solve_2d_poisson_staggered(this, source_term, solution)
        !
        !--------------------------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------------------------
        class (PoissonSolver), intent (in out) :: this
        real (wp),             intent (in out) :: source_term(:,:)
        real (wp),             intent (out)    :: solution(:,:)
        !--------------------------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------------------------
        real (wp), parameter :: HELMHOLTZ_CONSTANT = 0.0_wp
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (this%initialized .eqv. .false.) then
            error stop 'Uninitialized object of class (PoissonSolver) '&
                //' in solve_2d_poisson_staggered'
        end if

        !
        !==> Solve poisson's equation on staggered grid
        !
        associate( &
            elmbda => HELMHOLTZ_CONSTANT, &
            rhs => source_term, &
            f => solution &
            )
            call this%solve_2d_helmholtz_staggered(elmbda, rhs, f)
        end associate

    end subroutine solve_2d_poisson_staggered



    subroutine finalize_poisson_solver(this)
        !--------------------------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------------------------
        type (PoissonSolver), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_poisson_solver



end module type_PoissonSolver
