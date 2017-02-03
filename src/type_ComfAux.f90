!
!     file type_ComfAux.f90
!
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                         Fishpack                              *
!     *                                                               *
!     *                      A Package of Fortran                     *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *               for Modeling Geophysical Processes              *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        John Adams, Paul Swarztrauber and Roland Sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! The entries in this package are low-level
!  entries, supporting fishpack entries blktri
!  and cblktri. that is, these routines are
!  not called directly by users, but rather
!  by entries within blktri and cblktri.
!
module type_ComfAux

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: ppsgf
    public :: ppspf
    public :: psgf
    public :: ComfAux
    public :: comf_interface

    ! Parameters confined to the module
    real(wp), parameter :: ONE = 1.0_wp
    
    type, public :: ComfAux
    contains
        ! Type-bound procedures
        procedure, nopass, public :: ppsgf
        procedure, nopass, public :: ppspf
        procedure, nopass, public :: psgf
    end type ComfAux

    ! Declare interface
    interface
        function comf_interface(x, iz, c, a, bh) &
            result (return_value)
            import :: ip, wp

            ! Dummy arguments
            integer(ip), intent(in) :: iz
            real(wp),    intent(in) :: x
            real(wp),    intent(in) :: c(*)
            real(wp),    intent(in) :: a(*)
            real(wp),    intent(in) :: bh(*)
            real(wp)                :: return_value
        end function comf_interface
    end interface

contains

    pure function ppsgf(x, iz, c, a, bh) &
        result (return_value)

        ! Dummy arguments
        integer(ip),    intent(in)    :: iz
        real(wp),       intent(in)    :: x
        real(wp),       intent(in)    :: c(*)
        real(wp),       intent(in)    :: a(*)
        real(wp),       intent(in)    :: bh(*)
        real(wp)                      :: return_value

        return_value = sum(ONE/(x - bh(1:iz))**2)

    end function ppsgf

    pure function ppspf(x, iz, c, a, bh) &
        result (return_value)

        ! Dummy arguments
        integer(ip),    intent(in)    :: iz
        real(wp),       intent(in)    :: x
        real(wp),       intent(in)    :: c(*)
        real(wp),       intent(in)    :: a(*)
        real(wp),       intent(in)    :: bh(*)
        real(wp)                      :: return_value

        return_value = sum(ONE/(x - bh(1:iz)))

    end function ppspf

    pure function psgf(x, iz, c, a, bh) &
        result (return_value)

        ! Dummy arguments
        integer(ip),    intent(in)    :: iz
        real(wp),       intent(in)    :: x
        real(wp),       intent(in)    :: c(*)
        real(wp),       intent(in)    :: a(*)
        real(wp),       intent(in)    :: bh(*)
        real(wp)                      :: return_value

        ! Local variables
        integer(ip) :: j
        real(wp)    :: fsg, hsg, dd

        fsg = ONE
        hsg = ONE

        do j = 1, iz
            dd = ONE/(x - bh(j))
            fsg = fsg * a(j) * dd
            hsg = hsg * c(j) * dd
        end do

        select case (mod(iz,2))
            case (0)
                return_value = ONE - fsg - hsg
            case default
                return_value = ONE + fsg + hsg
        end select

    end function psgf

end module type_ComfAux
!
! REVISION HISTORY
!
! September 1973    Version 1
! April     1976    Version 2
! January   1978    Version 3
! December  1979    Version 3.1
! February  1985    Documentation upgrade
! November  1988    Version 3.2, FORTRAN 77 changes
! June      2004    Version 5.0, Fortran 90 changes
! April     2016    Modern Fortran (2008+) changes
!
