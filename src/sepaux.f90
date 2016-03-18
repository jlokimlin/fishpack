!
!     file sepaux.f
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
! PACKAGE SEPAUX         CONTAINS NO USER ENTRY POINTS.
!
! LATEST REVISION        June 2004
!
! PURPOSE                THIS PACKAGE CONTAINS AUXILIARY ROUTINES FOR
!                        THE FISHPACK SOLVERS SEPELI AND sepx4.
!
! USAGE                  SINCE THIS PACKAGE CONTAINS NO USER ENTRIES,
!                        NO USAGE INSTRUCTIONS OR ARGUMENT DESCRIPTIONS
!                        ARE GIVEN HERE.
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
! LANGUAGE               FORTRAN 90
!
! HISTORY                DEVELOPED IN THE LATE 1970'S BY JOHN C. ADAMS
!                        OF NCAR'S SCIENTTIFIC COMPUTING DIVISION.
!                        Revised in June 2004 incorporating fortran 90
!                        features
!
! PORTABILITY            FORTRAN 90
!
module module_sepaux

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: seport
    public :: sepmin
    public :: septri
    public :: sepdx
    public :: sepdy
    public :: get_coefficients

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables shared with sepeli.f90 and sepx4.f90
    !---------------------------------------------------------------------------------
    integer (ip), save, public :: kswx, kswy, k, l, mit, nit, is, ms, js, ns
    real (wp),    save, public :: ait, bit, cit, dit
    real (wp),    save, public :: dlx, dly, tdlx3, tdly3, dlx4, dly4
    !---------------------------------------------------------------------------------

    interface

        pure subroutine get_coefficients( grid, a, b, c)
            import :: wp
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            real (wp), intent (in)  :: grid
            real (wp), intent (out) :: a
            real (wp), intent (out) :: b
            real (wp), intent (out) :: c
            !-----------------------------------------------
        end subroutine get_coefficients

    end interface


contains


    subroutine seport(usol, idmn, zn, zm, pertrb)
        !
        ! Purpose:
        !
        !    this subroutine orthoganalizes the array usol with respect to
        !     the constant array in a weighted least squares norm
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: idmn
        real (wp),    intent (out)    :: pertrb
        real (wp),    intent (in out) :: usol(idmn, 1)
        real (wp),    intent (in)     :: zn(*)
        real (wp),    intent (in)     :: zm(*)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: istr, ifnl, jstr, jfnl, i, ii, j, jj
        real (wp)    :: ute, ete
        !--------------------------------------------------------------------------------

        istr = is
        ifnl = ms
        jstr = js
        jfnl = ns
        !
        !     compute weighted inner products
        !
        ute = 0.0
        ete = 0.0
        do i = is, ms
            ii = i - is + 1
            ete = ete + sum(zm(ii)*zn(:ns-js+1))
            ute = ute + sum(usol(i, js:ns)*zm(ii)*zn(:ns-js+1))
        end do
        !
        !     set perturbation parameter
        !
        pertrb = ute/ete
        !
        !     subtract off constant pertrb
        !
        usol(istr:ifnl, jstr:jfnl) = usol(istr:ifnl, jstr:jfnl) - pertrb

    end subroutine seport


    subroutine sepmin(usol, idmn, zn, zm, pertb)
        !
        ! Purpose:
        !
        !     this subroutine orhtogonalizes the array usol with respect to
        !     the constant array in a weighted least squares norm
        !
        !
        !     entry at sepmin occurrs when the final solution is
        !     to be minimized with respect to the weighted
        !     least squares norm
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: idmn
        real (wp),    intent (out)    :: pertb
        real (wp),    intent (in out) :: usol(idmn, 1)
        real (wp),    intent (in)     :: zn(*)
        real (wp),    intent (in)     :: zm(*)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: istr, ifnl, jstr, jfnl, i, ii, j, jj
        real (wp)    :: ute, ete, pertrb
        !-----------------------------------------------

        istr = 1
        ifnl = k
        jstr = 1
        jfnl = l

        ! compute weighted inner products
        ute = 0.0_wp
        ete = 0.0_wp
        do i = is, ms
            ii = i - is + 1
            ete = ete + sum(zm(ii) * zn(:ns-js+1))
            ute = ute + sum(usol(i, js:ns) * zm(ii) * zn(:ns-js+1))
        end do

        ! set perturbation parameter
        pertrb = ute/ete

        ! subtract off constant pertrb
        usol(istr:ifnl, jstr:jfnl) = usol(istr:ifnl, jstr:jfnl) - pertrb

    end subroutine sepmin


    subroutine septri(n, a, b, c, d, u, z)
            !
        !     this subroutine solves for a non-zero eigenvector corresponding
        !     to the zero eigenvalue of the transpose of the rank
        !     deficient one matrix with subdiagonal a, diagonal b, and
        !     superdiagonal c , with a(1) in the (1, n) position, with
        !     c(n) in the (n, 1) position, and all other elements zero.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: n
        real (wp),    intent (in)     :: a(n)
        real (wp),    intent (in)     :: b(n)
        real (wp),    intent (in)     :: c(n)
        real (wp),    intent (in out) :: d(n)
        real (wp),    intent (in out) :: u(n)
        real (wp),    intent (in out) :: z(n)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: nm2, j, nm1, k
        real (wp)    :: bn, v, den, an
        !--------------------------------------------------------------------------------

        bn = b(n)
        d(1) = a(2)/b(1)
        v = a(1)
        u(1) = c(n)/b(1)
        nm2 = n - 2
        do j = 2, nm2
            den = b(j) - c(j-1)*d(j-1)
            d(j) = a(j+1)/den
            u(j) = -c(j-1)*u(j-1)/den
            bn = bn - v*u(j-1)
            v = -v*d(j-1)
        end do
        den = b(n-1) - c(n-2)*d(n-2)
        d(n-1) = (a(n)-c(n-2)*u(n-2))/den
        an = c(n-1) - v*d(n-2)
        bn = bn - v*u(n-2)
        den = bn - an*d(n-1)
        !
        !     set last component equal to one
        !
        z(n) = 1.0_wp
        z(n-1) = -d(n-1)
        nm1 = n - 1
        do j = 2, nm1
            k = n - j
            z(k) = (-d(k)*z(k+1)) - u(k)*z(n)
        end do

    end subroutine septri


    pure subroutine sepdx(u, idmn, i, j, uxxx, uxxxx)

        !
        !     this program computes second order finite difference
        !     approximations to the third and fourth x
        !     partial derivatives of u at the (i, j) mesh point
        !
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip), intent (in) :: idmn
        integer (ip), intent (in) :: i
        integer (ip), intent (in) :: j
        real (wp),    intent (out) :: uxxx
        real (wp),    intent (out) :: uxxxx
        real (wp),    intent (in) :: u(idmn, 1)
        !--------------------------------------------------------------------------------

        !     compute partial derivative approximations at x=a
        !
        if (i == 1) then
            if (kswx /= 1) then
                uxxx = ((-5.0_wp * u(1, j))+18.0_wp * u(2, j)-24.0_wp * u(3, j)+14.0_wp * u(4, j)- &
                    3.0_wp * u(5, j))/tdlx3
                uxxxx = (3.0_wp * u(1, j)-14.0_wp * u(2, j)+26.0_wp * u(3, j)-24.0_wp * u(4, j)+11.0 &
                    *u(5, j)-2.0_wp * u(6, j))/dlx4
                return
            else
                !
                !     periodic at x=a
                !
                uxxx = ((-u(k-2, j))+2.0_wp * u(k-1, j)-2.0_wp * u(2, j)+u(3, j))/tdlx3
                uxxxx = (u(k-2, j)-4.0_wp * u(k-1, j)+6.0_wp * u(1, j)-4.0_wp * u(2, j)+u(3, j)) &
                    /dlx4
                return
            end if
        !
        !     compute partial derivative approximations at x=a+dlx
        !
        else if (i == 2) then
            if (kswx /= 1) then
                uxxx = ((-3.0_wp * u(1, j))+10.0_wp * u(2, j)-12.0_wp * u(3, j)+6.0_wp * u(4, j)-u(5 &
                    , j))/tdlx3
                uxxxx = (2.0_wp * u(1, j)-9.0_wp * u(2, j)+16.0_wp * u(3, j)-14.0_wp * u(4, j)+6.0_wp * u &
                    (5, j)-u(6, j))/dlx4
                return
            else
                !
                !     periodic at x=a+dlx
                !
                uxxx = ((-u(k-1, j))+2.0_wp * u(1, j)-2.0_wp * u(3, j)+u(4, j))/tdlx3
                uxxxx = (u(k-1, j)-4.0_wp * u(1, j)+6.0_wp * u(2, j)-4.0_wp * u(3, j)+u(4, j))/ &
                    dlx4
                return
            end if
        else if (i>2 .and. i<k-1) then
            !
            !     compute partial derivative approximations on the interior
            !
            uxxx = ((-u(i-2, j))+2.0_wp * u(i-1, j)-2.0_wp * u(i+1, j)+u(i+2, j))/tdlx3
            uxxxx = (u(i-2, j)-4.0_wp * u(i-1, j)+6.0_wp * u(i, j)-4.0_wp * u(i+1, j)+u(i+2, j) &
                )/dlx4
            return
        else if (i == k - 1) then
            !
            !     compute partial derivative approximations at x=b-dlx
            !
            if (kswx /= 1) then
                uxxx = (u(k-4, j)-6.0_wp * u(k-3, j)+12.0_wp * u(k-2, j)-10.0_wp * u(k-1, j)+ &
                    3.0_wp * u(k, j))/tdlx3
                uxxxx = ((-u(k-5, j))+6.0_wp * u(k-4, j)-14.0_wp * u(k-3, j)+16.0_wp * u(k-2, j &
                    )-9.0_wp * u(k-1, j)+2.0_wp * u(k, j))/dlx4
                return
            else
                !
                !     periodic at x=b-dlx
                !
                uxxx = ((-u(k-3, j))+2.0_wp * u(k-2, j)-2.0_wp * u(1, j)+u(2, j))/tdlx3
                uxxxx = (u(k-3, j)-4.0_wp * u(k-2, j)+6.0_wp * u(k-1, j)-4.0_wp * u(1, j)+u(2, j &
                    ))/dlx4
                return
            end if
        else if (i == k) then
            !
            !     compute partial derivative approximations at x=b
            !
            uxxx = -(3.0_wp * u(k-4, j)-14.0_wp * u(k-3, j)+24.0_wp * u(k-2, j)&
                -18.0_wp * u(k-1, j) + 5.0_wp * u(k, j))/tdlx3

            uxxxx = ((-2.0_wp * u(k-5, j))+11.0_wp * u(k-4, j)-24.0_wp * u(k-3, j)&
                +26.0_wp * u(k-2, j)-14.0_wp * u(k-1, j)+3.0_wp * u(k, j))/dlx4
            return
        end if

    end subroutine sepdx


    pure subroutine sepdy(u, idmn, i, j, uyyy, uyyyy)
        !
        ! Purpose:
        !
        !     this program computes second order finite difference
        !     approximations to the third and fourth y
        !     partial derivatives of u at the (i, j) mesh point
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)  :: idmn
        integer (ip), intent (in)  :: i
        integer (ip), intent (in)  :: j
        real (wp),    intent (out) :: uyyy
        real (wp),    intent (out) :: uyyyy
        real (wp),    intent (in)  :: u(idmn, 6)
        !--------------------------------------------------------------------------------


        !     compute partial derivative approximations at y=c
        !
        if (j == 1) then
            if (kswy /= 1) then
                uyyy = ((-5.0_wp * u(i, 1))+18.0_wp * u(i, 2)-24.0_wp * u(i, 3)+14.0_wp * u(i, 4)- &
                    3.0_wp * u(i, 5))/tdly3
                uyyyy = (3.0_wp * u(i, 1)-14.0_wp * u(i, 2)+26.0_wp * u(i, 3)-24.0_wp * u(i, 4)+11.0 &
                    *u(i, 5)-2.0_wp * u(i, 6))/dly4
                return
            else
                !
                !     periodic at x=a
                !
                uyyy = ((-u(i, l-2))+2.0_wp * u(i, l-1)-2.0_wp * u(i, 2)+u(i, 3))/tdly3
                uyyyy = (u(i, l-2)-4.0_wp * u(i, l-1)+6.0_wp * u(i, 1)-4.0_wp * u(i, 2)+u(i, 3)) &
                    /dly4
                return
            end if
        !
        !     compute partial derivative approximations at y=c+dly
        !
        else if (j == 2) then
            if (kswy /= 1) then
                uyyy = ((-3.0_wp * u(i, 1))+10.0_wp * u(i, 2)-12.0_wp * u(i, 3)+6.0_wp * u(i, 4)-u(i &
                    , 5))/tdly3
                uyyyy = (2.0_wp * u(i, 1)-9.0_wp * u(i, 2)+16.0_wp * u(i, 3)-14.0_wp * u(i, 4)+6.0_wp * u &
                    (i, 5)-u(i, 6))/dly4
                return
            else
                !
                !     periodic at y=c+dly
                !
                uyyy = ((-u(i, l-1))+2.0_wp * u(i, 1)-2.0_wp * u(i, 3)+u(i, 4))/tdly3
                uyyyy = (u(i, l-1)-4.0_wp * u(i, 1)+6.0_wp * u(i, 2)-4.0_wp * u(i, 3)+u(i, 4))/ &
                    dly4
                return
            end if
        !
        !     compute partial derivative approximations on the interior
        !
        else if (j>2 .and. j<l-1) then
            uyyy = ((-u(i, j-2))+2.0_wp * u(i, j-1)-2.0_wp * u(i, j+1)+u(i, j+2))/tdly3
            uyyyy = (u(i, j-2)-4.0_wp * u(i, j-1)+6.0_wp * u(i, j)-4.0_wp * u(i, j+1)+u(i, j+2) &
                )/dly4
            return
        else if (j == l - 1) then
            !
            !     compute partial derivative approximations at y=d-dly
            !
            if (kswy /= 1) then
                uyyy = (u(i, l-4)-6.0_wp * u(i, l-3)+12.0_wp * u(i, l-2)-10.0_wp * u(i, l-1)+ &
                    3.0_wp * u(i, l))/tdly3
                uyyyy = ((-u(i, l-5))+6.0_wp * u(i, l-4)-14.0_wp * u(i, l-3)+16.0_wp * u(i, l-2 &
                    )-9.0_wp * u(i, l-1)+2.0_wp * u(i, l))/dly4
                return
            else
                !
                !     periodic at y=d-dly
                !
                uyyy = ((-u(i, l-3))+2.0_wp * u(i, l-2)-2.0_wp * u(i, 1)+u(i, 2))/tdly3
                uyyyy = (u(i, l-3)-4.0_wp * u(i, l-2)+6.0_wp * u(i, l-1)-4.0_wp * u(i, 1)+u(i, 2 &
                    ))/dly4
                return
            end if
        else if (j == l) then
            !
            !     compute partial derivative approximations at y=d
            !
            uyyy = -(3.0_wp * u(i, l-4)-14.0_wp * u(i, l-3)+24.0_wp * u(i, l-2)-18.0_wp * u(i, l-1) &
                +5.0_wp * u(i, l))/tdly3
            uyyyy = ((-2.0_wp * u(i, l-5))+11.0_wp * u(i, l-4)-24.0_wp * u(i, l-3)+26.0_wp * u(i, l &
                -2)-14.0_wp * u(i, l-1)+3.0_wp * u(i, l))/dly4
            return
        end if

    end subroutine sepdy


end module module_sepaux
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    version 5.0, fortran 90 changes
!-----------------------------------------------------------------------
