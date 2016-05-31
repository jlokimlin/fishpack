module type_RectangularDomain

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT

    ! Explicit typing only
    implicit None

    ! Everything is private unless stated otherwise
    private
    public :: RectangularDomain

    ! Declare derived data type
    type, public :: RectangularDomain
        !-------------------------------------------------------------------------------
        ! Class variables
        !-------------------------------------------------------------------------------
        logical,   public  :: initialized = .false.
        real (wp), public  :: A = 0.0_wp  !! A <= x <= B
        real (wp), public  :: B = 0.0_wp  !! A <= x <= B
        real (wp), public  :: C = 0.0_wp  !! C <= y <= D
        real (wp), public  :: D = 0.0_wp  !! C <= y <= D
        !-------------------------------------------------------------------------------
    contains
        !-------------------------------------------------------------------------------
        ! Class methods
        !-------------------------------------------------------------------------------
        procedure, non_overridable, public  :: create => create_rectangular_domain
        procedure, non_overridable, public  :: destroy => destroy_rectangular_domain
        generic,                    public  :: assignment (=) => copy_rectangular_domain
        procedure,                  private :: copy_rectangular_domain
        final                               :: finalize_rectangular_domain
        !-------------------------------------------------------------------------------
    end type RectangularDomain


contains


    subroutine create_rectangular_domain(this, x_interval, y_interval )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (RectangularDomain), intent (in out)  :: this
        real (wp),                 intent (in)      :: x_interval(:) !! A <= x <= B
        real (wp),                 intent (in)      :: y_interval(:) !! C <= y <= D
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Set x interval: A <= x <= B
        this%A = x_interval(1)
        this%B = x_interval(2)

        ! Set y interval: A <= y <= B
        this%C = y_interval(1)
        this%D = y_interval(2)

        ! Set initialization flag
        this%initialized = .true.

    end subroutine create_rectangular_domain


    subroutine destroy_rectangular_domain( this )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (RectangularDomain), intent (in out)   :: this
        !--------------------------------------------------------------------------------
        
        ! Check if object is already usable
        if (.not.this%initialized) return
        
        ! Reset constants
        this%A = 0.0_wp
        this%B = 0.0_wp
        this%C = 0.0_wp
        this%D = 0.0_wp

        ! Reset initialization flag
        this%initialized = .false.

    end subroutine destroy_rectangular_domain


    subroutine copy_rectangular_domain(this, rectangle_to_be_copied )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (RectangularDomain), intent (out) :: this
        class (RectangularDomain), intent (in)  :: rectangle_to_be_copied
        !--------------------------------------------------------------------------------

        this%A = rectangle_to_be_copied%A
        this%B = rectangle_to_be_copied%B
        this%C = rectangle_to_be_copied%C
        this%D = rectangle_to_be_copied%D

        ! Set initialization flag
        this%initialized = .true.

    end subroutine copy_rectangular_domain


    subroutine finalize_rectangular_domain( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (RectangularDomain), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_rectangular_domain
    

end module type_RectangularDomain
