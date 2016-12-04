!
!     file comf.f90
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
!     *                    FISHPACK90  Version 1.1                    *
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
! PACKAGE COMF           The entries in this package are low-level
!                        entries, supporting fishpack entries blktri
!                        and cblktri. that is, these routines are
!                        not called directly by users, but rather
!                        by entries within blktri and cblktri.
!
!
! LATEST REVISION        April 2016
!
! SPECIAL CONDITIONS     None
!
! I/O                    None
!
! PRECISION              Set in the module fishpack_precision.f90
!
! REQUIRED LIBRARY       None
! FILES
!
! STANDARD               Fortran 2008
!
module module_comf

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

    ! Declare derived data type
    type, public :: ComfAux
        !--------------------------------------------------
        ! Type components
        !--------------------------------------------------
    contains
        !--------------------------------------------------
        ! Type-bound procedures
        !--------------------------------------------------
        procedure, nopass, public :: ppsgf
        procedure, nopass, public :: ppspf
        procedure, nopass, public :: psgf
        !--------------------------------------------------
    end type ComfAux

    ! Declare interface
    interface
        function comf_interface(x, iz, c, a, bh) result (return_value)
            import :: ip, wp
            !-----------------------------------------------
            ! Dummy arguments
            !-----------------------------------------------
            integer(ip), intent(in) :: iz
            real(wp),    intent(in) :: x
            real(wp),    intent(in) :: c(*)
            real(wp),    intent(in) :: a(*)
            real(wp),    intent(in) :: bh(*)
            real(wp)                :: return_value
            !-----------------------------------------------
        end function comf_interface
    end interface


contains


    pure function ppsgf(x, iz, c, a, bh) result (return_value)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in) :: iz
        real(wp),    intent(in) :: x
        real(wp),    intent(in) :: c(*)
        real(wp),    intent(in) :: a(*)
        real(wp),    intent(in) :: bh(*)
        real(wp)                :: return_value
        !-----------------------------------------------

        return_value = sum(1.0_wp/(x - bh(1:iz))**2)

    end function ppsgf


    pure function ppspf(x, iz, c, a, bh) result (return_value)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in) :: iz
        real(wp),    intent(in) :: x
        real(wp),    intent(in) :: c(*)
        real(wp),    intent(in) :: a(*)
        real(wp),    intent(in) :: bh(*)
        real(wp)                 :: return_value
        !-----------------------------------------------

        return_value = sum(1.0_wp/(x - bh(1:iz)))

    end function ppspf


    pure function psgf(x, iz, c, a, bh) result (return_value)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in) :: iz
        real(wp),    intent(in) :: x
        real(wp),    intent(in) :: c(*)
        real(wp),    intent(in) :: a(*)
        real(wp),    intent(in) :: bh(*)
        real(wp)                 :: return_value
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: j
        real(wp)    :: fsg, hsg, dd
        !-----------------------------------------------

        fsg = 1.0_wp
        hsg = 1.0_wp

        do j = 1, iz
            dd = 1.0_wp/(x - bh(j))
            fsg = fsg * a(j) * dd
            hsg = hsg * c(j) * dd
        end do

        select case (mod(iz,2))
            case (0)
                return_value = 1.0_wp - fsg - hsg
            case default
                return_value = 1.0_wp + fsg + hsg
        end select

    end function psgf


end module module_comf
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
! April     2016    Replaced epmach with intrinsic function epsilon
!                   and pimach with acos(-1.0_wp) where wp = REAL64 from the
!                   intrinsic module ISO_Fortran_env
!
