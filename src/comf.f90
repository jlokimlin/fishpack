!
!     file comf.f
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
!     *                    FISHPACK90  version 1.1                    *
!     *                                                               *
!     *                 A Package of Fortran 77 and 90                *
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
! PACKAGE COMF           THE ENTRIES IN THIS PACKAGE ARE LOWLEVEL
!                        ENTRIES, SUPPORTING FISHPACK ENTRIES blktri
!                        AND cblktri. THAT IS, THESE ROUTINES ARE
!                        NOT CALLED DIRECTLY BY USERS, BUT RATHER
!                        BY ENTRIES WITHIN blktri AND cblktri.
!
!
! LATEST REVISION        April 2016
!
! SPECIAL CONDITIONS     NONE
!
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED LIBRARY       NONE
! FILES
!
! LANGUAGE               Fortran 2008+
!
module module_comf

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: ppsgf
    public :: ppspf
    public :: psgf
    public :: comf_interface

    interface
        pure function comf_interface(x, iz, c, a, bh) result (return_value)
            import :: ip, wp
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in) :: iz
            real (wp),    intent (in) :: x
            real (wp),    intent (in) :: c(*)
            real (wp),    intent (in) :: a(*)
            real (wp),    intent (in) :: bh(*)
            real (wp)                 :: return_value
            !-----------------------------------------------
        end function
    end interface


contains


    pure function ppsgf(x, iz, c, a, bh) result (return_value)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in) :: iz
        real (wp),    intent (in) :: x
        real (wp),    intent (in) :: c(*)
        real (wp),    intent (in) :: a(*)
        real (wp),    intent (in) :: bh(*)
        real (wp)                 :: return_value
        !-----------------------------------------------

        return_value = sum(1.0_wp/(x - bh(1:iz))**2)

    end function ppsgf


    pure function ppspf(x, iz, c, a, bh) result (return_value)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in) :: iz
        real (wp),    intent (in) :: x
        real (wp),    intent (in) :: c(*)
        real (wp),    intent (in) :: a(*)
        real (wp),    intent (in) :: bh(*)
        real (wp)                 :: return_value
        !-----------------------------------------------

        return_value = sum(1.0_wp/(x - bh(1:iz)))

    end function ppspf


    pure function psgf(x, iz, c, a, bh) result (return_value)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in) :: iz
        real (wp),    intent (in) :: x
        real (wp),    intent (in) :: c(*)
        real (wp),    intent (in) :: a(*)
        real (wp),    intent (in) :: bh(*)
        real (wp)                 :: return_value
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: j
        real (wp)    :: fsg, hsg
        !-----------------------------------------------

        fsg = 1.0_wp
        hsg = 1.0_wp

        do j = 1, iz
            associate( dd => 1.0_wp/(x - bh(j)) )
                fsg = fsg*a(j)*dd
                hsg = hsg*c(j)*dd
            end associate
        end do

        if (mod(iz, 2) == 0) then
            return_value = 1.0_wp - fsg - hsg
        else
            return_value = 1.0_wp + fsg + hsg
        end if

    end function psgf


end module module_comf
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 changes
! April     2016    Replaced epmach with intrinsic function epsilon
!                   and pimach with acos(-1.0_wp) where wp = REAL64 from the
!                   intrinsic module iso_fortran_env
!
