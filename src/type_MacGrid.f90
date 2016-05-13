!
! Purpose:
!
! Defines the (Marker-and-Cell (MAC) grid object to discretize functions for
! the projection method
!
! Reference:
!
! Harlow, Francis H., and J. Eddie Welch.
! "Numerical calculation of time-dependent viscous incompressible 
! flow of fluid with free surface." 
! Physics of fluids 8.12 (1965): 2182.
!
!*****************************************************************************************
!
! Components of the MacGrid object are:
!
! DOMAIN (RectangularDomain) ::
!
! Derived data type with components
!
! A (real) :: Lower horizontal endpoint
! B (real) :: Upper horizontal endpoint
!
! The range of the horizontal grid in x,
!
! i.e.,  A <= x <= B. Where A must be less than B.
!
! C (real) :: Lower vertical endpoint
! D (real) :: Upper vertical endpoint
!
! The range of the vertical grid in y,
!
! i.e.,  C <= x <= D. Where C must be less than D.
!
!--------------------------------------------------------------------------------
!
! NX (integer) ::
!
! Number of horizontally staggered grid points in x
!
!--------------------------------------------------------------------------------
!
! NY (integer) ::
!
! Number of vertically staggered grid points in y
!
!--------------------------------------------------------------------------------
!
! X_CENTERED (real 1-dimensional array) ::
!
! Contains the points into which the interval (A, B) is subdivided.
! Hence, there will be (nx + 1)-many  grid points in the x-direction
! given by
!
! x(i) = A + (i - 1)dx for i = 1, 2, ..., nx + 1,
!
! where dx = (B - A)/nx is the panel width.
!
! nx must be greater than 3.
!
!--------------------------------------------------------------------------------
!
! X_STAGGERED (real 1-dimensional array) ::
!
! Contains the points into which the interval (A, B) is subdivided.
! Hence, there will be nx-many  grid points in the x-direction
! given by
!
! x(i) = A + (i - 1/2)dx for i = 1, 2, ..., nx,
!
! where dx = (B - A)/nx is the panel width.
!
! nx must be greater than 2.
!
!--------------------------------------------------------------------------------
!
! Y_CENTERED (real 1-dimensional array) ::
!
! Contains the points into which the interval (C, D) is subdivided.
! Hence, there will be (ny + 1)-many  grid points in the y-direction
! given by
!
! y(j) = C + (j - 1)dy for i = 1, 2, ..., ny + 1,
!
! where dy = (D-C)/ny is the panel width.
!
! ny must be greater than 3.
!
!--------------------------------------------------------------------------------
!
! Y_STAGGERED (real 1-dimensional array) ::
!
! Contains the points into which the interval (C, D) is subdivided.
! Hence, there will be ny-many  grid points in the y-direction
! given by
!
! y(j) = C + (j - 1/2)dy for i = 1, 2, ..., ny,
!
! where dy = (D-C)/ny is the panel width.
!
! ny must be greater than 2.
!
!
module type_MacGrid

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
    public :: MacGrid

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    character (len=250) :: error_message      !! Probably long enough
    integer (ip)        :: deallocate_status  !! To check deallocation status
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type, extends( Grid ), public :: MacGrid
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
        real (wp), allocatable, public :: x_centered(:)
        real (wp), allocatable, public :: y_centered(:)
        real (wp), allocatable, public :: x_staggered(:)
        real (wp), allocatable, public :: y_staggered(:)
        !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Class methods
        !---------------------------------------------------------------------------------
        procedure, public :: create => create_mac_grid
        procedure, public :: destroy => destroy_mac_grid
        procedure, public :: unformatted_print => print_to_unformatted_binary_files
        final             :: finalize_mac_grid
        !---------------------------------------------------------------------------------
    end type MacGrid


    ! Declare constructor
    interface MacGrid
        module procedure mac_grid_constructor
    end interface



contains



    function mac_grid_constructor(x_interval, y_interval, nx, ny) result (return_value)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), contiguous, intent (in) :: x_interval(:) !! Interval: A <= x <= B
        real (wp), contiguous, intent (in) :: y_interval(:) !! Interval: C <= y <= D
        integer (ip),          intent (in) :: nx  !! Number of horizontally staggered grid points in x
        integer (ip),          intent (in) :: ny  !! Number of vertically staggered grid points in y
        type (MacGrid)                     :: return_value
        !--------------------------------------------------------------------------------

        call return_value%create(x_interval, y_interval, nx, ny)

    end function



    subroutine create_mac_grid(this, x_interval, y_interval, nx, ny )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (MacGrid),       intent (in out) :: this
        real (wp), contiguous, intent (in)     :: x_interval(:) !! Interval: A <= x <= B
        real (wp), contiguous, intent (in)     :: y_interval(:) !! Interval: C <= y <= D
        integer (ip),          intent (in)     :: nx  !! Number of horizontally staggered grid points in x
        integer (ip),          intent (in)     :: ny  !! Number of vertically staggered grid points in y
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Create parent type
        call this%create_grid( x_interval, y_interval, nx, ny )

        associate( &
            A => x_interval(1), &
            B => x_interval(2), &
            C => y_interval(1), &
            D => y_interval(2) &
            )

            call this%get_centered_grids( A, C, nx, ny, this%x_centered, this%y_centered )

            call this%get_staggered_grids( A, C, nx, ny, this%x_staggered, this%y_staggered )

        end associate

        ! Set status
        this%initialized = .true.

    end subroutine create_mac_grid


    subroutine destroy_mac_grid( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (MacGrid), intent (in out) :: this
        !--------------------------------------------------------------------------------

        ! Deallocate horizontally centered grid in x
        if ( allocated( this%x_centered ) ) then

            ! Deallocate grid
            deallocate ( &
                this%x_centered, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (MacGrid)'
                write( stderr, '(A)' ) 'Deallocating X_CENTERED failed in DESTROY_MAC_GRID'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Deallocate horizontally staggered grid in x
        if ( allocated( this%x_staggered ) ) then

            ! Deallocate grid
            deallocate ( &
                this%x_staggered, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (MacGrid)'
                write( stderr, '(A)' ) 'Deallocating X_STAGGERED failed in DESTROY_MAC_GRID'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Deallocate vertically centered grid in y
        if ( allocated( this%y_centered ) ) then

            ! Deallocate grid
            deallocate ( &
                this%y_centered, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (MacGrid)'
                write( stderr, '(A)' ) 'Deallocating Y_CENTERED failed in DESTROY_MAC_GRID'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Deallocate vertically staggered grid in y
        if ( allocated( this%y_staggered ) ) then

            ! Deallocate grid
            deallocate ( &
                this%y_staggered, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (MacGrid)'
                write( stderr, '(A)' ) 'Deallocating Y_STAGGERED failed in DESTROY_MAC_GRID'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Destroy parent type
        call this%destroy_grid()

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_mac_grid


    subroutine print_to_unformatted_binary_files(this, header )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (MacGrid),   intent (in out) :: this
        character (len=*), intent (in)     :: header
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)  :: file_unit
        integer (ip)  :: record_length
        !--------------------------------------------------------------------------------

        associate( &
            nx => this%NX, &
            ny => this%NY &
            )
            associate( &
                x_staggered => this%x_staggered( 1:nx ), &
                y_staggered => this%y_staggered( 1:ny ), &
                x_centered => this%x_centered( 1:nx + 1), &
                y_centered => this%y_centered( 1:ny + 1) &
                )

                ! Write x_staggered
                inquire( iolength = record_length ) x_staggered
                open( newunit=file_unit, file=header//'x_staggered.dat', status='replace', &
                    form='unformatted', access='direct', recl=record_length )
                write( file_unit, rec=1 ) x_staggered
                close( file_unit )

                ! Write y_staggered
                inquire( iolength = record_length ) y_staggered
                open( newunit=file_unit, file=header//'y_staggered.dat', status='replace', &
                    form='unformatted', access='direct', recl=record_length )
                write( file_unit, rec=1 ) y_staggered
                close( file_unit )

                ! Write x_centered
                inquire( iolength = record_length ) x_centered
                open( newunit=file_unit, file=header//'x_centered.dat', status='replace', &
                    form='unformatted', access='direct', recl=record_length )
                write( file_unit, rec=1 ) x_centered
                close( file_unit )

                ! Write y_centered
                inquire( iolength = record_length ) y_centered
                open( newunit=file_unit, file=header//'y_centered.dat', status='replace', &
                    form='unformatted', access='direct', recl=record_length )
                write( file_unit, rec=1 ) y_centered
                close( file_unit )

            end associate
        end associate

    end subroutine print_to_unformatted_binary_files


    subroutine finalize_mac_grid( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (MacGrid), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_mac_grid


end module type_MacGrid
