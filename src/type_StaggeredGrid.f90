module type_StaggeredGrid

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT

    use type_Grid, only: &
        Grid

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: StaggeredGrid

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    character (len=250) :: error_message !! Probably long enough
    integer (ip)        :: deallocate_status !! To check deallocation status
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type, extends( Grid ), public :: StaggeredGrid
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
        real (wp), allocatable, public :: x(:)
        real (wp), allocatable, public :: y(:)
        !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Class methods
        !---------------------------------------------------------------------------------
        procedure, public :: create => create_staggered_grid
        procedure, public :: destroy => destroy_staggered_grid
        final             :: finalize_staggered_grid
        !---------------------------------------------------------------------------------
    end type StaggeredGrid


    ! Declare constructor
    interface StaggeredGrid
        module procedure staggered_grid_constructor
    end interface



contains



    function staggered_grid_constructor(x_interval, y_interval, nx, ny) result (return_value)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), contiguous, intent (in) :: x_interval(:) !! Interval: A <= x <= B
        real (wp), contiguous, intent (in) :: y_interval(:) !! Interval: C <= y <= D
        integer (ip),          intent (in) :: nx  !! Number of horizontally staggered grid points in x
        integer (ip),          intent (in) :: ny  !! Number of vertically staggered grid points in y
        type (StaggeredGrid)               :: return_value
        !--------------------------------------------------------------------------------

        call return_value%create(x_interval, y_interval, nx, ny)

    end function staggered_grid_constructor


    subroutine create_staggered_grid(this, x_interval, y_interval, nx, ny)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (StaggeredGrid), intent (in out) :: this
        real (wp), contiguous, intent (in)     :: x_interval(:) !! Interval: A <= x <= B
        real (wp), contiguous, intent (in)     :: y_interval(:) !! Interval: C <= y <= D
        integer (ip),          intent (in)     :: nx !! Number of horizontally staggered grid points in x
        integer (ip),          intent (in)     :: ny !! Number of vertically staggered grid points in y
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Set constants
        this%NX = nx
        this%NY = ny

        ! Create parent type
        call this%create_grid(x_interval, y_interval, nx, ny)

        associate( &
            A => x_interval(1), &
            C => y_interval(1) &
            )
            call this%get_staggered_grids( A, C, nx, ny, this%x, this%y )
        end associate

        ! Set status
        this%initialized = .true.

    end subroutine create_staggered_grid


    subroutine destroy_staggered_grid( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (StaggeredGrid), intent (in out) :: this
        !--------------------------------------------------------------------------------

        ! Deallocate horizontally staggered grid in x
        if ( allocated( this%x ) ) then

            ! Deallocate grid
            deallocate ( &
                this%x, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (StaggeredGrid)'
                write( stderr, '(A)' ) 'Deallocating X failed in DESTROY_STAGGERED_GRID'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Deallocate vertically staggered grid in y
        if ( allocated(this%y) ) then

            ! Deallocate grid
            deallocate ( &
                this%y, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (StaggeredGrid)'
                write( stderr, '(A)' ) 'Deallocating Y failed in DESTROY_STAGGERED_GRID'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! destroy component type
        call this%domain%destroy()

        ! Destroy parent type
        call this%destroy_grid()

        ! Reset initialization flag
        this%initialized = .false.

    end subroutine destroy_staggered_grid


    subroutine finalize_staggered_grid( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (StaggeredGrid), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_staggered_grid


end module type_StaggeredGrid
