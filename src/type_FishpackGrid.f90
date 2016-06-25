module type_FishpackGrid

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: FishpackGrid

    ! Declare derived data type
    type, public ::  FishpackGrid
        !--------------------------------------------------------------
        ! Class variables
        !--------------------------------------------------------------
    contains
        !--------------------------------------------------------------
        ! Class methods
        !--------------------------------------------------------------
        procedure, nopass :: linspace
        procedure, nopass :: get_staggered_grid
        procedure, nopass :: get_centered_grid
        !--------------------------------------------------------------
    end type FishpackGrid


contains


    function linspace(start, stop, num, endpoint) &
        result (return_value)
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        real (wp),              intent (in) :: start
        real (wp),              intent (in) :: stop
        integer (ip), optional, intent (in) :: num
        logical,      optional, intent (in) :: endpoint
        real (wp), allocatable              :: return_value(:)
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip)     :: num_op, div, i
        real (wp)        :: delta, step
        !--------------------------------------------------------------
        
        ! Number of samples to generate. Default is 50. Must be non-negative
        if (present(num)) then
            if (num < 0) error stop "linspace: num < 0"
            num_op = num
        else
            num_op = 50
        end if

        ! Allocate memory
        allocate( return_value(num_op) )

        ! If True, stop is the last sample. Otherwise, it is not included. Default is True.
        div = num_op
        if (present(endpoint)) then
            if (endpoint) div = (num_op-1)
        end if

        delta = stop - start
        step = delta / div

        do concurrent (i=1:num_op)
            return_value(i) = start + real(i,kind=wp)*step
        end do

    end function linspace



    function get_centered_grid(start, stop, num) &
        result (return_value)
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        real (wp),    intent (in) :: start
        real (wp),    intent (in) :: stop
        integer (ip), intent (in) :: num
        real (wp), allocatable    :: return_value(:)
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip) :: i
        real (wp)    :: delta, step
        !--------------------------------------------------------------

        if (num < 0) error stop "get_centered_grid: num < 0"

        ! Allocate memory
        allocate( return_value(num + 1) )

        delta = stop - start
        step = delta / num

        do concurrent (i=1:size(return_value))
            return_value(i) = start + real(i - 1,kind=wp)*step
        end do

    end function get_centered_grid



    function get_staggered_grid(start, stop, num) &
        result (return_value)
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        real (wp),    intent (in) :: start
        real (wp),    intent (in) :: stop
        integer (ip), intent (in) :: num
        real (wp), allocatable    :: return_value(:)
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip) :: i
        real (wp)    :: step
        !--------------------------------------------------------------

        if (num < 0) error stop "get_staggered_grid: num < 0"

        ! Allocate memory
        allocate( return_value(num + 1) )

        step = (stop - start) / num

        do concurrent (i=1:size(return_value))
            return_value(i) = start + (real(i,kind=wp) - 0.5_wp) * step
        end do

    end function get_staggered_grid

end module type_FishpackGrid
