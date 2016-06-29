module fishpack_precision

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: wp, ip, PI

    integer, parameter :: sp = selected_real_kind(p=6, r=37)
    integer, parameter :: dp = selected_real_kind(p=15, r=307)
    integer, parameter :: qp = selected_real_kind(p=33, r=4931)
    integer, parameter :: wp = dp

    real (wp), parameter :: PI = acos(-1.0_wp)


end module fishpack_precision
