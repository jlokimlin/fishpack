module fishpack_precision

    use, intrinsic :: iso_fortran_env, only: &
        INT32

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: wp, ip
    public :: PI, TWO_PI, EPS
    public :: pimach, epmach

    !-----------------------------------------------
    ! Dictionary: precision constants
    !-----------------------------------------------
    integer,   parameter :: ip = INT32
    integer,   parameter :: sp = selected_real_kind(p=6, r=37)
    integer,   parameter :: dp = selected_real_kind(p=15, r=307)
    integer,   parameter :: qp = selected_real_kind(p=33, r=4931)
    integer,   parameter :: wp = dp
    real (wp), parameter :: PI = acos(-1.0_wp)
    real (wp), parameter :: TWO_PI = 2.0_wp * PI
    real (wp), parameter :: EPS = epsilon(1.0_wp)
    !-----------------------------------------------

contains


    pure function pimach(dum) result (return_value)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        real (wp), optional, intent (in) :: dum
        real (wp)                        :: return_value
        !-----------------------------------------------

        return_value = 3.141592653589793238462643383279502884197169399375105820974_wp

    end function pimach


    pure function epmach(dum) result (return_value)
        !
        ! Purpose:
        !
        ! Computes an approximate machine epsilon (accuracy), i.e.,
        !
        ! the smallest number EPS of the kind wp such that 1 + EPS > 1.
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        real (wp), optional, intent (in) :: dum
        real (wp)                        :: return_value
        !-----------------------------------------------

        return_value = epsilon(1.0_wp)

    end function epmach


end module fishpack_precision
