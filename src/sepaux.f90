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
!                        THE FISHPACK SOLVERS SEPELI AND SEPX4.
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
! **********************************************************************
module module_sepaux

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: seport
    public :: sepmin
    public :: septri
    public :: sepdx
    public :: sepdy

contains

    subroutine SEPORT(usol, idmn, zn, zm, pertrb)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: idmn
        real , intent (out) :: pertrb
        real , intent (in out) :: usol(idmn, 1)
        real , intent (in) :: zn(*)
        real , intent (in) :: zm(*)
        !-----------------------------------------------
        !   C o m m o n   B l o c k s
        !-----------------------------------------------
        !...  /SPLP/
        common /SPLP/ kswx, kswy, k, l, ait, bit, cit, dit, mit, nit, is, &
            ms, js, ns, dlx, dly, tdlx3, tdly3, dlx4, dly4
        integer   kswx, kswy, k, l, mit, nit, is, ms, js, ns
        real   ait, bit, cit, dit, dlx, dly, tdlx3, tdly3, dlx4, dly4
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: istr, ifnl, jstr, jfnl, i, ii, j, jj
        real :: ute, ete
        !-----------------------------------------------
        !
        !     THIS SUBROUTINE ORTHOGANALIZES THE ARRAY USOL WITH RESPECT TO
        !     THE CONSTANT ARRAY IN A WEIGHTED LEAST SQUARES NORM
        !
        istr = is
        ifnl = ms
        jstr = js
        jfnl = ns
        !
        !     COMPUTE WEIGHTED INNER PRODUCTS
        !
        ute = 0.0
        ete = 0.0
        do i = is, ms
            ii = i - is + 1
            ete = ete + SUM(ZM(ii)*ZN(:ns-js+1))
            ute = ute + SUM(USOL(i, js:ns)*ZM(ii)*ZN(:ns-js+1))
        end do
        !
        !     SET PERTURBATION PARAMETER
        !
        pertrb = ute/ete
        !
        !     SUBTRACT OFF CONSTANT PERTRB
        !
        usol(istr:ifnl, jstr:jfnl) = USOL(istr:ifnl, jstr:jfnl) - pertrb

    end subroutine SEPORT

    subroutine SEPMIN(usol, idmn, zn, zm, pertb)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: idmn
        real  :: pertb
        real , intent (in out) :: usol(idmn, 1)
        real , intent (in) :: zn(*)
        real , intent (in) :: zm(*)
        !-----------------------------------------------
        !   C o m m o n   B l o c k s
        !-----------------------------------------------
        !...  /SPLP/
        common /SPLP/ kswx, kswy, k, l, ait, bit, cit, dit, mit, nit, is, &
            ms, js, ns, dlx, dly, tdlx3, tdly3, dlx4, dly4
        integer   kswx, kswy, k, l, mit, nit, is, ms, js, ns
        real   ait, bit, cit, dit, dlx, dly, tdlx3, tdly3, dlx4, dly4
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: istr, ifnl, jstr, jfnl, i, ii, j, jj
        real :: ute, ete, pertrb
        !-----------------------------------------------
        !
        !     THIS SUBROUTINE ORHTOGONALIZES THE ARRAY USOL WITH RESPECT TO
        !     THE CONSTANT ARRAY IN A WEIGHTED LEAST SQUARES NORM
        !
        !
        !     ENTRY AT SEPMIN OCCURRS WHEN THE FINAL SOLUTION IS
        !     TO BE MINIMIZED WITH RESPECT TO THE WEIGHTED
        !     LEAST SQUARES NORM
        !
        istr = 1
        ifnl = k
        jstr = 1
        jfnl = l
        !
        !     COMPUTE WEIGHTED INNER PRODUCTS
        !
        ute = 0.0
        ete = 0.0
        do i = is, ms
            ii = i - is + 1
            ete = ete + SUM(ZM(ii)*ZN(:ns-js+1))
            ute = ute + SUM(USOL(i, js:ns)*ZM(ii)*ZN(:ns-js+1))
        end do
        !
        !     SET PERTURBATION PARAMETER
        !
        pertrb = ute/ete
        !
        !     SUBTRACT OFF CONSTANT PERTRB
        !
        usol(istr:ifnl, jstr:jfnl) = USOL(istr:ifnl, jstr:jfnl) - pertrb

    end subroutine SEPMIN

    subroutine SEPTRI(n, a, b, c, d, u, z)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: n
        real , intent (in) :: a(n)
        real , intent (in) :: b(n)
        real , intent (in) :: c(n)
        real , intent (in out) :: d(n)
        real , intent (in out) :: u(n)
        real , intent (in out) :: z(n)
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: nm2, j, nm1, k
        real :: bn, v, den, an
        !-----------------------------------------------
        !
        !     THIS SUBROUTINE SOLVES FOR A NON-ZERO EIGENVECTOR CORRESPONDING
        !     TO THE ZERO EIGENVALUE OF THE TRANSPOSE OF THE RANK
        !     DEFICIENT ONE MATRIX WITH SUBDIAGONAL A, DIAGONAL B, AND
        !     SUPERDIAGONAL C , WITH A(1) IN THE (1, N) POSITION, WITH
        !     C(N) IN THE (N, 1) POSITION, AND ALL OTHER ELEMENTS ZERO.
        !
        bn = B(n)
        d(1) = A(2)/B(1)
        v = A(1)
        u(1) = C(n)/B(1)
        nm2 = n - 2
        do j = 2, nm2
            den = B(j) - C(j-1)*D(j-1)
            d(j) = A(j+1)/den
            u(j) = -C(j-1)*U(j-1)/den
            bn = bn - v*U(j-1)
            v = -v*D(j-1)
        end do
        den = B(n-1) - C(n-2)*D(n-2)
        d(n-1) = (A(n)-C(n-2)*U(n-2))/den
        an = C(n-1) - v*D(n-2)
        bn = bn - v*U(n-2)
        den = bn - an*D(n-1)
        !
        !     SET LAST COMPONENT EQUAL TO ONE
        !
        z(n) = 1.0
        z(n-1) = -D(n-1)
        nm1 = n - 1
        do j = 2, nm1
            k = n - j
            z(k) = (-D(k)*Z(k+1)) - U(k)*Z(n)
        end do

    end subroutine SEPTRI

    subroutine SEPDX(u, idmn, i, j, uxxx, uxxxx)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: idmn
        integer , intent (in) :: i
        integer , intent (in) :: j
        real , intent (out) :: uxxx
        real , intent (out) :: uxxxx
        real , intent (in) :: u(idmn, 1)
        !-----------------------------------------------
        !   C o m m o n   B l o c k s
        !-----------------------------------------------
        !...  /SPLP/
        common /SPLP/ kswx, kswy, k, l, ait, bit, cit, dit, mit, nit, is, &
            ms, js, ns, dlx, dly, tdlx3, tdly3, dlx4, dly4
        integer   kswx, kswy, k, l, mit, nit, is, ms, js, ns
        real   ait, bit, cit, dit, dlx, dly, tdlx3, tdly3, dlx4, dly4
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        !-----------------------------------------------
        !
        !     THIS PROGRAM COMPUTES SECOND ORDER FINITE DIFFERENCE
        !     APPROXIMATIONS TO THE THIRD AND FOURTH X
        !     PARTIAL DERIVATIVES OF U AT THE (I, J) MESH POINT
        !
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A
        !
        if (i == 1) then
            if (kswx /= 1) then
                uxxx = ((-5.0*U(1, j))+18.0*U(2, j)-24.0*U(3, j)+14.0*U(4, j)- &
                    3.0*U(5, j))/tdlx3
                uxxxx = (3.0*U(1, j)-14.0*U(2, j)+26.0*U(3, j)-24.0*U(4, j)+11.0 &
                    *U(5, j)-2.0*U(6, j))/dlx4
                return
            else
                !
                !     PERIODIC AT X=A
                !
                uxxx = ((-U(k-2, j))+2.0*U(k-1, j)-2.0*U(2, j)+U(3, j))/tdlx3
                uxxxx = (U(k-2, j)-4.0*U(k-1, j)+6.0*U(1, j)-4.0*U(2, j)+U(3, j)) &
                    /dlx4
                return
            end if
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A+DLX
        !
        else if (i == 2) then
            if (kswx /= 1) then
                uxxx = ((-3.0*U(1, j))+10.0*U(2, j)-12.0*U(3, j)+6.0*U(4, j)-U(5 &
                    , j))/tdlx3
                uxxxx = (2.0*U(1, j)-9.0*U(2, j)+16.0*U(3, j)-14.0*U(4, j)+6.0*U &
                    (5, j)-U(6, j))/dlx4
                return
            else
                !
                !     PERIODIC AT X=A+DLX
                !
                uxxx = ((-U(k-1, j))+2.0*U(1, j)-2.0*U(3, j)+U(4, j))/tdlx3
                uxxxx = (U(k-1, j)-4.0*U(1, j)+6.0*U(2, j)-4.0*U(3, j)+U(4, j))/ &
                    dlx4
                return
            end if
        else if (i>2 .and. i<k-1) then
            !
            !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
            !
            uxxx = ((-U(i-2, j))+2.0*U(i-1, j)-2.0*U(i+1, j)+U(i+2, j))/tdlx3
            uxxxx = (U(i-2, j)-4.0*U(i-1, j)+6.0*U(i, j)-4.0*U(i+1, j)+U(i+2, j) &
                )/dlx4
            return
        else if (i == k - 1) then
            !
            !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B-DLX
            !
            if (kswx /= 1) then
                uxxx = (U(k-4, j)-6.0*U(k-3, j)+12.0*U(k-2, j)-10.0*U(k-1, j)+ &
                    3.0*U(k, j))/tdlx3
                uxxxx = ((-U(k-5, j))+6.0*U(k-4, j)-14.0*U(k-3, j)+16.0*U(k-2, j &
                    )-9.0*U(k-1, j)+2.0*U(k, j))/dlx4
                return
            else
                !
                !     PERIODIC AT X=B-DLX
                !
                uxxx = ((-U(k-3, j))+2.0*U(k-2, j)-2.0*U(1, j)+U(2, j))/tdlx3
                uxxxx = (U(k-3, j)-4.0*U(k-2, j)+6.0*U(k-1, j)-4.0*U(1, j)+U(2, j &
                    ))/dlx4
                return
            end if
        else if (i == k) then
            !
            !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B
            !
            uxxx = -(3.0*U(k-4, j)-14.0*U(k-3, j)+24.0*U(k-2, j)-18.0*U(k-1, j) &
                +5.0*U(k, j))/tdlx3
            uxxxx = ((-2.0*U(k-5, j))+11.0*U(k-4, j)-24.0*U(k-3, j)+26.0*U(k-2 &
                , j)-14.0*U(k-1, j)+3.0*U(k, j))/dlx4
            return
        end if

    end subroutine SEPDX

    subroutine SEPDY(u, idmn, i, j, uyyy, uyyyy)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: idmn
        integer , intent (in) :: i
        integer , intent (in) :: j
        real , intent (out) :: uyyy
        real , intent (out) :: uyyyy
        real , intent (in) :: u(idmn, 6)
        !-----------------------------------------------
        !   C o m m o n   B l o c k s
        !-----------------------------------------------
        !...  /SPLP/
        common /SPLP/ kswx, kswy, k, l, ait, bit, cit, dit, mit, nit, is, &
            ms, js, ns, dlx, dly, tdlx3, tdly3, dlx4, dly4
        integer   kswx, kswy, k, l, mit, nit, is, ms, js, ns
        real   ait, bit, cit, dit, dlx, dly, tdlx3, tdly3, dlx4, dly4
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        !-----------------------------------------------
        !
        !     THIS PROGRAM COMPUTES SECOND ORDER FINITE DIFFERENCE
        !     APPROXIMATIONS TO THE THIRD AND FOURTH Y
        !     PARTIAL DERIVATIVES OF U AT THE (I, J) MESH POINT
        !
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C
        !
        if (j == 1) then
            if (kswy /= 1) then
                uyyy = ((-5.0*U(i, 1))+18.0*U(i, 2)-24.0*U(i, 3)+14.0*U(i, 4)- &
                    3.0*U(i, 5))/tdly3
                uyyyy = (3.0*U(i, 1)-14.0*U(i, 2)+26.0*U(i, 3)-24.0*U(i, 4)+11.0 &
                    *U(i, 5)-2.0*U(i, 6))/dly4
                return
            else
                !
                !     PERIODIC AT X=A
                !
                uyyy = ((-U(i, l-2))+2.0*U(i, l-1)-2.0*U(i, 2)+U(i, 3))/tdly3
                uyyyy = (U(i, l-2)-4.0*U(i, l-1)+6.0*U(i, 1)-4.0*U(i, 2)+U(i, 3)) &
                    /dly4
                return
            end if
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY
        !
        else if (j == 2) then
            if (kswy /= 1) then
                uyyy = ((-3.0*U(i, 1))+10.0*U(i, 2)-12.0*U(i, 3)+6.0*U(i, 4)-U(i &
                    , 5))/tdly3
                uyyyy = (2.0*U(i, 1)-9.0*U(i, 2)+16.0*U(i, 3)-14.0*U(i, 4)+6.0*U &
                    (i, 5)-U(i, 6))/dly4
                return
            else
                !
                !     PERIODIC AT Y=C+DLY
                !
                uyyy = ((-U(i, l-1))+2.0*U(i, 1)-2.0*U(i, 3)+U(i, 4))/tdly3
                uyyyy = (U(i, l-1)-4.0*U(i, 1)+6.0*U(i, 2)-4.0*U(i, 3)+U(i, 4))/ &
                    dly4
                return
            end if
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
        !
        else if (j>2 .and. j<l-1) then
            uyyy = ((-U(i, j-2))+2.0*U(i, j-1)-2.0*U(i, j+1)+U(i, j+2))/tdly3
            uyyyy = (U(i, j-2)-4.0*U(i, j-1)+6.0*U(i, j)-4.0*U(i, j+1)+U(i, j+2) &
                )/dly4
            return
        else if (j == l - 1) then
            !
            !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY
            !
            if (kswy /= 1) then
                uyyy = (U(i, l-4)-6.0*U(i, l-3)+12.0*U(i, l-2)-10.0*U(i, l-1)+ &
                    3.0*U(i, l))/tdly3
                uyyyy = ((-U(i, l-5))+6.0*U(i, l-4)-14.0*U(i, l-3)+16.0*U(i, l-2 &
                    )-9.0*U(i, l-1)+2.0*U(i, l))/dly4
                return
            else
                !
                !     PERIODIC AT Y=D-DLY
                !
                uyyy = ((-U(i, l-3))+2.0*U(i, l-2)-2.0*U(i, 1)+U(i, 2))/tdly3
                uyyyy = (U(i, l-3)-4.0*U(i, l-2)+6.0*U(i, l-1)-4.0*U(i, 1)+U(i, 2 &
                    ))/dly4
                return
            end if
        else if (j == l) then
            !
            !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D
            !
            uyyy = -(3.0*U(i, l-4)-14.0*U(i, l-3)+24.0*U(i, l-2)-18.0*U(i, l-1) &
                +5.0*U(i, l))/tdly3
            uyyyy = ((-2.0*U(i, l-5))+11.0*U(i, l-4)-24.0*U(i, l-3)+26.0*U(i, l &
                -2)-14.0*U(i, l-1)+3.0*U(i, l))/dly4
            return
        end if

    end subroutine SEPDY

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
