module type_Grid

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    use, intrinsic :: ISO_Fortran_env, only: &
        stderr => ERROR_UNIT

    use type_RectangularDomain, only: &
        RectangularDomain

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: Grid

    !---------------------------------------------------------------
    ! global variables confined to the module
    !---------------------------------------------------------------
    character(len=250) :: error_message !! Probably long enough
    integer(ip)        :: allocate_status !! To check allocation status
    integer(ip)        :: deallocate_status !! To check deallocation status
    !---------------------------------------------------------------

    ! Declare derived data type
    type, public ::  Grid
        !---------------------------------------------------------------
        ! Type components
        !---------------------------------------------------------------
        logical,                 public :: initialized = .false.
        integer(ip),             public :: NX = 0 !! Number of horizontally staggered grid points
        integer(ip),             public :: NY = 0 !! Number of vertically staggered grid points
        real(wp),                public :: DX = 0.0_wp !! Horizontal mesh
        real(wp),                public :: DY = 0.0_wp !! Vertical mesh
        real(wp),                public :: DX_SQUARED = 0.0_wp
        real(wp),                public :: DY_SQUARED = 0.0_wp
        type(RectangularDomain), public :: domain  !! A <= x <= B, C <= y <= D
        !---------------------------------------------------------------
    contains
        !---------------------------------------------------------------
        ! Type-bound procedures
        !---------------------------------------------------------------
        procedure,         public :: create => create_grid
        procedure,         public :: destroy => destroy_grid
        procedure,         public :: create_grid
        procedure,         public :: destroy_grid
        procedure, nopass, public :: compute_one_dimensional_grid
        procedure,         public :: get_discretization_mesh
        procedure,         public :: get_centered_grids
        procedure,         public :: get_staggered_grids
        !final                     :: finalize_grid
        !---------------------------------------------------------------
    end type Grid

contains

    subroutine create_grid(this, x_interval, y_interval, nx, ny)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(Grid),          intent(inout) :: this
        real(wp), contiguous, intent(in)    :: x_interval(:) !! Interval: A <= x <= B
        real(wp), contiguous, intent(in)    :: y_interval(:) !! Interval: C <= y <= D
        integer(ip),          intent(in)    :: nx  !! Number of horizontally staggered grid points in x
        integer(ip),          intent(in)    :: ny  !! Number of vertically staggered grid points in y
        !--------------------------------------------------------------

        ! Ensure object is usable
        call this%destroy()

        ! Set constants
        this%NX = nx
        this%NY = ny

        associate( &
            A => x_interval(1), &
            B => x_interval(2), &
            C => y_interval(1), &
            D => y_interval(2) &
            )

            call this%domain%create(x_interval, y_interval)

            call this%get_discretization_mesh(A, B, C, D, nx, ny)

        end associate

        ! Set status
        this%initialized = .true.

    end subroutine create_grid


    subroutine destroy_grid(this)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(Grid), intent(inout) :: this
        !--------------------------------------------------------------

        ! Check if object is already usable
        if (.not.this%initialized) return
        
        ! Reset constants
        this%NX = 0
        this%DX = 0.0_wp
        this%DX_SQUARED = 0.0_wp
        this%NY = 0
        this%DY = 0.0_wp
        this%DY_SQUARED = 0.0_wp

        ! Reset initialization flag
        this%initialized = .false.
        
        ! destroy parent type
        call this%domain%destroy()
        
    end subroutine destroy_grid


    subroutine get_discretization_mesh(this, A, B, C, D, nx, ny )
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(Grid), intent(inout) :: this
        real(wp),    intent(in)     :: A  ! Interval: A <= x <= B
        real(wp),    intent(in)     :: B  ! Interval: A <= x <= B
        real(wp),    intent(in)     :: C  ! Interval: C <= y <= D
        real(wp),    intent(in)     :: D  ! Interval: C <= y <= D
        integer(ip), intent(in)     :: nx ! Number of horizontally staggered grid points
        integer(ip), intent(in)     :: ny ! Number of vertically staggered grid points
        !--------------------------------------------------------------

        associate( &
            dx => (B - A) / nx, &
            dy => (D - C) / ny &
            )

            ! Set horizontal mesh
            this%DX = dx
            this%DX_SQUARED = (dx ** 2)

            ! Set vertical mesh
            this%DY = dy
            this%DY_SQUARED = (dy ** 2)

        end associate

    end subroutine get_discretization_mesh


    subroutine compute_one_dimensional_grid( &
        lower_bound, upper_bound, left_interval_endpoint, uniform_mesh, grid, staggered )
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip),           intent(in)  :: lower_bound
        integer(ip),           intent(in)  :: upper_bound
        real(wp),              intent(in)  :: left_interval_endpoint
        real(wp),              intent(in)  :: uniform_mesh
        real(wp), allocatable, intent(out) :: grid(:)
        logical, optional,      intent(in)  :: staggered
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip):: n !! Counter
        !--------------------------------------------------------------

        ! Handle the case where array is already allocated
        if ( allocated( grid ) ) then

            ! Deallocate grid
            deallocate ( &
                grid, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(a)' ) 'type(Grid)'
                write( stderr, '(a)' ) 'Deallocation failed in COMPUTE_ONE_DIMENSIONAL_GRID'
                write( stderr, '(a)' ) trim( error_message )
            end if
        end if

        ! Allocate grid
        allocate ( &
            grid( lower_bound:upper_bound ), &
            stat=allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            write( stderr, '(a)' ) 'type(Grid)'
            write( stderr, '(a)' ) 'Allocation failed in COMPUTE_ONE_DIMENSIONAL_GRID'
            write( stderr, '(a)' ) trim( error_message )
        end if

        ! Set constants
        associate( &
            t0 => left_interval_endpoint, &
            dt => uniform_mesh &
            )
            ! Check if grid is staggered or centered
            if ( present(staggered) ) then
                do n = lower_bound, upper_bound
                    ! Compute staggered grid
                    grid(n) = t0 +  (real(n, kind=wp) - 0.5_wp) * dt
                end do
            else
                ! Compute centered grid
                do n = lower_bound, upper_bound
                    grid(n) = t0 + real(n - 1, kind=wp) * dt
                end do
            end if
        end associate

    end subroutine compute_one_dimensional_grid


    subroutine get_centered_grids(this, A, C, nx, ny, x, y )
        ! Remark:
        ! The subroutine "get_discretization_mesh" must be called first to initialize dx and dy
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(Grid),           intent(inout) :: this
        real(wp),              intent(in)     :: A  ! A <= x
        real(wp),              intent(in)     :: C  ! C <= y
        integer(ip),           intent(in)     :: nx ! Number of horizontal grid points in x
        integer(ip),           intent(in)     :: ny ! Number of vertical grid points in y
        real(wp), allocatable, intent(out)    :: x(:)
        real(wp), allocatable, intent(out)    :: y(:)
        !--------------------------------------------------------------

        associate ( &
            dx => this%DX, &
            dy => this%DY &
            )

             ! Set horizontally centered ghost grid
            call compute_one_dimensional_grid( -1, nx + 3, A, dx, x )

            ! Set vertically centered ghost grid
            call compute_one_dimensional_grid(  -1, ny + 3, C, dy, y )

        end associate

    end subroutine get_centered_grids


    subroutine get_staggered_grids(this, A, C, nx, ny, x, y )
        ! Remark:
        ! The subroutine "get_discretization_mesh" must be called first to initialize dx and dy
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(Grid),           intent(inout) :: this
        real(wp),              intent(in)     :: A  ! A <= x
        real(wp),              intent(in)     :: C  ! C <= y
        integer(ip),           intent(in)     :: nx ! Number of horizontal grid points in x
        integer(ip),           intent(in)     :: ny ! Number of vertical grid points in y
        real(wp), allocatable, intent(out)    :: x(:)
        real(wp), allocatable, intent(out)    :: y(:)
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        logical, parameter :: staggered = .true.
        !--------------------------------------------------------------

        associate ( &
            dx => this%DX, &
            dy => this%DY &
            )
             ! Set horizontally staggered ghost grid
            call compute_one_dimensional_grid( -1, nx + 2, A, dx, x, staggered )

            ! Set vertically staggered ghost grid
            call compute_one_dimensional_grid( -1, ny + 2, C, dy, y, staggered )
        end associate

    end subroutine get_staggered_grids

    subroutine finalize_grid(this)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        type(Grid), intent(inout) :: this
        !--------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_grid

end module type_Grid
