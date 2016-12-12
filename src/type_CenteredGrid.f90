module type_CenteredGrid

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    use, intrinsic :: ISO_Fortran_env, only: &
        stderr => ERROR_UNIT

    use type_Grid, only: &
        Grid

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: CenteredGrid

    !---------------------------------------------------------------
    ! global variables confined to the module
    !---------------------------------------------------------------
    integer(ip) :: deallocate_status  !! To check deallocation status
    !---------------------------------------------------------------

    
    type, extends (Grid), public :: CenteredGrid
        !---------------------------------------------------------------
        ! Type components
        !---------------------------------------------------------------
        real(wp), allocatable, public :: x(:)
        real(wp), allocatable, public :: y(:)
        !---------------------------------------------------------------
    contains
        !---------------------------------------------------------------
        ! Type-bound procedures
        !---------------------------------------------------------------
        procedure, public :: create => create_centered_grid
        procedure, public :: destroy => destroy_centered_grid
        !---------------------------------------------------------------
    end type CenteredGrid


    ! Declare constructor
    interface CenteredGrid
        module procedure centered_grid_constructor
    end interface



contains



    function centered_grid_constructor(x_interval, y_interval, nx, ny) result (return_value)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real(wp), contiguous, intent(in) :: x_interval(:) !! Interval: A <= x <= B
        real(wp), contiguous, intent(in) :: y_interval(:) !! Interval: C <= y <= D
        integer(ip),          intent(in) :: nx  !! Number of horizontally staggered grid points in x
        integer(ip),          intent(in) :: ny  !! Number of vertically staggered grid points in y
        type(CenteredGrid)                :: return_value
        !--------------------------------------------------------------

        call return_value%create(x_interval, y_interval, nx, ny)

    end function centered_grid_constructor


    subroutine create_centered_grid(this, x_interval, y_interval, nx, ny )
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(CenteredGrid),  intent(inout) :: this
        real(wp), contiguous, intent(in)     :: x_interval(:) !! Interval: A <= x <= B
        real(wp), contiguous, intent(in)     :: y_interval(:) !! Interval: C <= y <= D
        integer(ip),          intent(in)     :: nx !! Number of horizontally staggered grid points in x
        integer(ip),          intent(in)     :: ny !! Number of vertically staggered grid points in y
        !--------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Set constants
        this%NX = nx
        this%NY = ny

        ! Allocate parent type
        call this%create_grid( x_interval, y_interval, nx, ny )

        associate( &
            A => x_interval(1), &
            C => y_interval(1) &
            )
            call this%get_centered_grids( A, C, nx, ny, this%x, this%y )
        end associate

        ! Set status
        this%initialized = .true.

    end subroutine create_centered_grid

    subroutine destroy_centered_grid(this)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(CenteredGrid), intent(inout) :: this
        !--------------------------------------------------------------

        ! Check initialization flag
        if (this%initialized .eqv. .false. ) return

        ! Deallocate horizontally staggered grid in x
        if ( allocated( this%x ) ) then
            ! Deallocate grid
            deallocate ( this%x, stat=deallocate_status )
            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                error stop 'Object of class(CenteredGrid): '&
                    //'Deallocating X failed in DESTROY_CENTERED_GRID'
            end if
        end if

        ! Deallocate vertically staggered grid in y
        if ( allocated(this%y) ) then
            ! Deallocate grid
            deallocate ( this%y, stat=deallocate_status )
            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                error stop 'Object of class(CenteredGrid): '&
                    //'Deallocating Y failed in DESTROY_CENTERED_GRID'
            end if
        end if

        ! destroy component type
        call this%domain%destroy()

        ! Destroy parent type
        call this%destroy_grid()

        ! Reset initialization flag
        this%initialized = .false.

    end subroutine destroy_centered_grid

end module type_CenteredGrid
