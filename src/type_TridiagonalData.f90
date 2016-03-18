module type_TridiagonalData

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: TridiagonalData

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    character (len=250) :: error_message !! Probably long enough
    integer (ip)        :: allocate_status  !! To check allocation status
    integer (ip)        :: deallocate_status !! To check deallocation status
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type, public :: TridiagonalData
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
        logical,                public :: initialized = .false.
        integer (ip),           public :: X_BOUNDARY_CONDITION = -1
        integer (ip),           public :: Y_BOUNDARY_CONDITION = -1
        real (wp), allocatable, public :: subdiagonal(:)
        real (wp), allocatable, public :: diagonal(:)
        real (wp), allocatable, public :: superdiagonal(:)
        !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Class methods
        !---------------------------------------------------------------------------------
        procedure,                  public  :: create => create_tridiagonal_data
        procedure,                  public  :: destroy => destroy_tridiagonal_data
        procedure, non_overridable, public  :: create_tridiagonal_data
        procedure, non_overridable, public  :: destroy_tridiagonal_data
        procedure, nopass,          private :: get_x_boundary_condition_type
        procedure, nopass,          private :: get_y_boundary_condition_type
        final                               :: finalize_tridiagonal_data
        !---------------------------------------------------------------------------------
    end type TridiagonalData


contains


    subroutine create_tridiagonal_data( this, nx, x_type, y_type )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (TridiagonalData), intent (in out)  :: this
        integer (ip),            intent (in)      :: nx
        integer (ip),            intent (in)      :: x_type
        integer (ip),            intent (in)      :: y_type
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Allocate arrays
        allocate ( &
            this%subdiagonal( nx ), &
            this%diagonal( nx ), &
            this%superdiagonal( nx ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            write( stderr, '(A)' ) 'TYPE (TridiagonalData)'
            write( stderr, '(A)' ) 'Allocation failed in CREATE_TRIDIAGONAL_DATA'
            write( stderr, '(A)' ) trim( error_message )
        end if

        ! Set boundary condition types
        call get_x_boundary_condition_type( x_type, this%X_BOUNDARY_CONDITION )
        call get_y_boundary_condition_type( y_type, this%Y_BOUNDARY_CONDITION )

        ! Set initialization status
        this%initialized = .true.

    end subroutine create_tridiagonal_data


    subroutine destroy_tridiagonal_data( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (TridiagonalData), intent (in out)   :: this
        !--------------------------------------------------------------------------------

        if (.not.this%initialized) return

        ! Check if array is allocated
        if ( allocated( this%subdiagonal ) ) then

            ! Deallocate array
            deallocate( &
                this%subdiagonal, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TridiagonalData)'
                write( stderr, '(A)' ) 'Deallocating SUBDIAGONAL failed in DESTROY_TRIDIAGONAL_DATA'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%diagonal ) ) then

            ! Deallocate array
            deallocate( &
                this%diagonal, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TridiagonalData)'
                write( stderr, '(A)' ) 'Deallocating DIAGONAL failed in DESTROY_TRIDIAGONAL_DATA'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%superdiagonal ) ) then

            ! Deallocate array
            deallocate( &
                this%superdiagonal, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TridiagonalData)'
                write( stderr, '(A)' ) 'Deallocating SUPERDIAGONAL failed in DESTROY_TRIDIAGONAL_DATA'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Reset constants
        this%X_BOUNDARY_CONDITION = -1
        this%Y_BOUNDARY_CONDITION = -1

        ! Reset initialization flag
        this%initialized = .false.

    end subroutine destroy_tridiagonal_data


    subroutine get_y_boundary_condition_type( type, return_value )
        !
        ! Purpose:
        !
        ! Indicates values which x(i,0) and x(i,ny + 1)
        ! are assumed to have.
        !
        ! That is,
        !
        !  = 1 if x(i,0) = -x(i,1) and x(i,ny+1) = -x(i,ny)
        !  = 2 if x(i,0) = -x(i,1) and x(i,ny+1) =  x(i,ny)
        !  = 3 if x(i,0) =  x(i,1) and x(i,ny+1) =  x(i,ny)
        !  = 4 if x(i,0) =  x(i,1) and x(i,ny+1) = -x(i,ny)
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)   :: type
        integer (ip), intent (out)  :: return_value
        !--------------------------------------------------------------------------------

        ! Initialize return value
        return_value = -1

        select case (type)
            case(4)
                return_value = type
            case(3)
                return_value = type
            case(2)
                return_value = type
            case(1)
                return_value = type
            case(0)
                return_value = type
            case default
                ! Handle invalid boundary type
                write( stderr, '(A)' ) 'TYPE (TridiagonalData)'
                write( stderr, '(A)' ) 'Invalid boundary type'
                write( stderr, '(A,I11)' ) 'The y boundary condition argument type = ', type
                write( stderr, '(A)' ) 'must be either 0, 1, ..., 4'
        end select

    end subroutine get_y_boundary_condition_type


    subroutine get_x_boundary_condition_type( type, return_value )
        !
        ! Purpose:
        !
        ! = 0 if a(1) and c(nx) are not zero
        ! = 1 if a(1) = c(nx) = 0
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)   :: type
        integer (ip), intent (out)  :: return_value
        !--------------------------------------------------------------------------------

        ! Initialize return value
        return_value = -1

        select case (type)
            case(1)
                return_value = type
            case(0)
                return_value = type
            case default
                ! Handle invalid boundary type
                write( stderr, '(A)' ) 'TYPE (TridiagonalData)'
                write( stderr, '(A)' ) 'Invalid boundary type'
                write( stderr, '(A,I11)' ) 'The x boundary condition argument type = ', type
                write( stderr, '(A)' ) 'must be either 0, or 1'
        end select

    end subroutine get_x_boundary_condition_type


    subroutine finalize_tridiagonal_data( this )
        !
        ! Purpose:
        !< Finalize object
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (TridiagonalData), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_tridiagonal_data


end module type_TridiagonalData
