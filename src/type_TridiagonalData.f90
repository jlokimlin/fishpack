module type_TridiagonalData

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT

    use type_Grid, only: &
        Grid

    ! Explicit typing only
    implicit None

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
        logical,                             public :: initialized = .false.
        integer (ip),                        public :: X_BOUNDARY_CONDITION = -1
        integer (ip),                        public :: Y_BOUNDARY_CONDITION = -1
        real (wp), allocatable,              public :: subdiagonal(:)
        real (wp), allocatable,              public :: diagonal(:)
        real (wp), allocatable,              public :: superdiagonal(:)
        procedure (proc_interface), pointer, public :: assign_coefficients => null()
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


    abstract interface
        subroutine proc_interface(this, grid_type )
            import :: TridiagonalData, Grid, wp
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            class (TridiagonalData), intent (in out) :: this
            class (Grid),            intent (in out) :: grid_type
            !--------------------------------------------------------------------------------
        end subroutine proc_interface
    end interface


contains


    subroutine create_tridiagonal_data(this, nx, x_type, y_type, proc )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (TridiagonalData), intent (in out)  :: this
        integer (ip),            intent (in)      :: nx
        integer (ip),            intent (in)      :: x_type
        integer (ip),            intent (in)      :: y_type
        procedure (proc_interface), optional      :: proc
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Allocate arrays
        allocate ( &
            this%subdiagonal( nx ), &
            this%diagonal( nx ), &
            this%superdiagonal( nx ), &
            stat=allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            write( stderr, '(a)' ) 'TYPE (TridiagonalData)'
            write( stderr, '(a)' ) 'Allocation failed in CREATE_TRIDIAGONAL_DATA'
            write( stderr, '(a)' ) trim( error_message )
        end if

        ! Set boundary condition types
        associate( &
            x_bc => this%X_BOUNDARY_CONDITION, &
            y_bc => this%Y_BOUNDARY_CONDITION &
            )
            call this%get_x_boundary_condition_type( x_type, x_bc )
            call this%get_y_boundary_condition_type( y_type, y_bc )
        end associate

        ! Associate pointer
        if (present(proc)) then
            this%assign_coefficients => proc
        end if

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
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(a)' ) 'TYPE (TridiagonalData)'
                write( stderr, '(a)' ) 'Deallocating SUBDIAGONAL failed in DESTROY_TRIDIAGONAL_DATA'
                write( stderr, '(a)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%diagonal ) ) then

            ! Deallocate array
            deallocate( &
                this%diagonal, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(a)' ) 'TYPE (TridiagonalData)'
                write( stderr, '(a)' ) 'Deallocating DIAGONAL failed in DESTROY_TRIDIAGONAL_DATA'
                write( stderr, '(a)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%superdiagonal ) ) then

            ! Deallocate array
            deallocate( &
                this%superdiagonal, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(a)' ) 'TYPE (TridiagonalData)'
                write( stderr, '(a)' ) 'Deallocating SUPERDIAGONAL failed in DESTROY_TRIDIAGONAL_DATA'
                write( stderr, '(a)' ) trim( error_message )
            end if
        end if

        ! Reset constants
        this%X_BOUNDARY_CONDITION = -1
        this%Y_BOUNDARY_CONDITION = -1

        ! Nullify pointer
        if (associated(this%assign_coefficients)) nullify( this%assign_coefficients )

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
                write( stderr, '(a)' ) 'TYPE (TridiagonalData)'
                write( stderr, '(a)' ) 'Invalid boundary type'
                write( stderr, '(A,I11)' ) 'The y boundary condition argument type = ', type
                write( stderr, '(a)' ) 'must be either 0, 1, ..., 4'
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
                write( stderr, '(a)' ) 'TYPE (TridiagonalData)'
                write( stderr, '(a)' ) 'Invalid boundary type'
                write( stderr, '(A,I11)' ) 'The x boundary condition argument type = ', type
                write( stderr, '(a)' ) 'must be either 0, or 1'
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
