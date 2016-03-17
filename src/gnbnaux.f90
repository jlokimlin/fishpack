!
!     file gnbnaux.f
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
! PACKAGE GNBNAUX
!
! LATEST REVISION        June 2004
!
! PURPOSE                TO PROVIDE AUXILIARY ROUTINES FOR FISHPACK
!                        ENTRIES GENBUN AND POISTG.
!
! USAGE                  THERE ARE NO USER ENTRIES IN THIS PACKAGE.
!                        THE ROUTINES IN THIS PACKAGE ARE NOT INTENDED
!                        TO BE CALLED BY USERS, BUT RATHER BY ROUTINES
!                        IN PACKAGES GENBUN AND POISTG.
!
! SPECIAL CONDITIONS     NONE
!
! I/O                    NONE
!
! PRECISION              SINGLE
!
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN IN 1979 BY ROLAND SWEET OF NCAR'S
!                        SCIENTIFIC COMPUTING DIVISION.  MADE AVAILABLE
!                        ON NCAR'S PUBLIC LIBRARIES IN JANUARY, 1980.
!                        Revised by John Adams in June 2004 incorporating
!                        Fortran 90 features
!
! PORTABILITY            FORTRAN 90
! ********************************************************************
module module_gnbnaux

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: cosgen
    public :: merge_rename
    public :: trix
    public :: tri3

contains


    pure subroutine cosgen(n, ijump, fnum, fden, a)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip), intent (in) :: ijump
        real (wp),    intent (in) :: fnum
        real (wp),    intent (in) :: fden
        real (wp),    intent (out) :: a(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip)         :: k3, k4, k, k1, k5, i, k2, np1
        real (wp), parameter :: PI = acos( -1.0_wp)
        real (wp)            :: dum, pibyn, x, y
        !-----------------------------------------------
        !
        !
        !     THIS SUBROUTINE COMPUTES REQUIRED COSINE VALUES IN ASCENDING
        !     ORDER.  WHEN IJUMP .GT. 1 THE ROUTINE COMPUTES VALUES
        !
        !        2*COS(J*PI/L) , J=1, 2, ..., L AND J .NE. 0(MOD N/IJUMP+1)
        !
        !     WHERE L = IJUMP*(N/IJUMP+1).
        !
        !
        !     WHEN IJUMP = 1 IT COMPUTES
        !
        !            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... , N
        !
        !     WHERE
        !        FNUM = 0.5, FDEN = 0.0,  FOR REGULAR REDUCTION VALUES
        !        FNUM = 0.0, FDEN = 1.0, FOR B-R AND C-R WHEN ISTAG = 1
        !        FNUM = 0.0, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
        !        FNUM = 0.5, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
        !                                IN POISN2 ONLY.
        !
        !

        if (n /= 0) then
            if (ijump /= 1) then
                k3 = n/ijump + 1
                k4 = k3 - 1
                pibyn = PI/real(n + ijump)
                do k = 1, ijump
                    k1 = (k - 1)*k3
                    k5 = (k - 1)*k4
                    do i = 1, k4
                        x = k1 + i
                        k2 = k5 + i
                        a(k2) = -2.0_wp * cos(x*pibyn)
                    end do
                end do
            else
                np1 = n + 1
                y = PI/(real(n, kind = wp) + fden)

                do i = 1, n
                    x = real(np1 - i, kind = wp) - fnum
                    a(i) = 2.0_wp * cos(x*y)
                end do
            end if
        end if

    end subroutine cosgen


    subroutine trix(idegbr, idegcr, m, a, b, c, y, tcos, d, w)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in) :: idegbr
        integer (ip), intent (in) :: idegcr
        integer (ip), intent (in) :: m
        real (wp), intent (in) :: a(*)
        real (wp), intent (in) :: b(*)
        real (wp), intent (in) :: c(*)
        real (wp), intent (in out) :: y(*)
        real (wp), intent (in) :: tcos(*)
        real (wp), intent (in out) :: d(*)
        real (wp), intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: mm1, ifb, ifc, l, lint, k, i, ip
        real (wp)    :: x, xx, z
        !-----------------------------------------------
        !
        !     SUBROUTINE TO SOLVE A SYSTEM OF LINEAR EQUATIONS WHERE THE
        !     COEFFICIENT MAtrix IS A RATIONAL FUNCTION IN THE MAtrix GIVEN BY
        !     TRIDIAGONAL  ( . . . , A(I), B(I), C(I), . . . ).
        !
        mm1 = m - 1
        ifb = idegbr + 1
        ifc = idegcr + 1
        l = ifb/ifc
        lint = 1
        do k = 1, idegbr
            x = TCOS(k)
            if (k == l) then
                i = idegbr + lint
                xx = x - TCOS(i)
                w(:m) = Y(:m)
                y(:m) = xx*Y(:m)
            end if
            z = 1./(B(1)-x)
            d(1) = C(1)*z
            y(1) = Y(1)*z
            do i = 2, mm1
                z = 1./(B(i)-x-A(i)*D(i-1))
                d(i) = C(i)*z
                y(i) = (Y(i)-A(i)*Y(i-1))*z
            end do
            z = B(m) - x - A(m)*D(mm1)
            if (z == 0.) then
                y(m) = 0.
            else
                y(m) = (Y(m)-A(m)*Y(mm1))/z
            end if
            do ip = 1, mm1
                y(m-ip) = Y(m-ip) - D(m-ip)*Y(m+1-ip)
            end do
            if (k /= l) cycle
            y(:m) = Y(:m) + W(:m)
            lint = lint + 1
            l = (lint*ifb)/ifc
        end do

    end subroutine trix


    subroutine tri3(m, a, b, c, k, y1, y2, y3, tcos, d, w1, w2, w3)
        !
        ! Purpose:
        !
        ! subroutine to solve three linear systems whose common coefficient
        ! matrix is a rational function in the matrix given by
        !
        !  tridiagonal (..., a(i), b(i), c(i), ...)
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in) :: m
        integer (ip), intent (in) :: k(4)
        real (wp), intent (in) :: a(*)
        real (wp), intent (in) :: b(*)
        real (wp), intent (in) :: c(*)
        real (wp), intent (in out) :: y1(*)
        real (wp), intent (in out) :: y2(*)
        real (wp), intent (in out) :: y3(*)
        real (wp), intent (in) :: tcos(*)
        real (wp), intent (in out) :: d(*)
        real (wp), intent (in out) :: w1(*)
        real (wp), intent (in out) :: w2(*)
        real (wp), intent (in out) :: w3(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: mm1, k1, k2, k3, k4, if1, if2, if3, if4, k2k3k4, l1, l2
        integer (ip) :: l3, lint1, lint2, lint3, kint1, kint2, kint3, n, i, ipp
        real (wp)    :: x, z, xx
        !-----------------------------------------------
        !

        mm1 = m - 1
        k1 = K(1)
        k2 = K(2)
        k3 = K(3)
        k4 = K(4)
        if1 = k1 + 1
        if2 = k2 + 1
        if3 = k3 + 1
        if4 = k4 + 1
        k2k3k4 = k2 + k3 + k4
        if (k2k3k4 /= 0) then
            l1 = if1/if2
            l2 = if1/if3
            l3 = if1/if4
            lint1 = 1
            lint2 = 1
            lint3 = 1
            kint1 = k1
            kint2 = kint1 + k2
            kint3 = kint2 + k3
        end if
        do n = 1, k1
            x = TCOS(n)
            if (k2k3k4 /= 0) then
                if (n == l1) then
                    w1(:m) = Y1(:m)
                end if
                if (n == l2) then
                    w2(:m) = Y2(:m)
                end if
                if (n == l3) then
                    w3(:m) = Y3(:m)
                end if
            end if
            z = 1./(B(1)-x)
            d(1) = C(1)*z
            y1(1) = Y1(1)*z
            y2(1) = Y2(1)*z
            y3(1) = Y3(1)*z
            do i = 2, m
                z = 1./(B(i)-x-A(i)*D(i-1))
                d(i) = C(i)*z
                y1(i) = (Y1(i)-A(i)*Y1(i-1))*z
                y2(i) = (Y2(i)-A(i)*Y2(i-1))*z
                y3(i) = (Y3(i)-A(i)*Y3(i-1))*z
            end do
            do ipp = 1, mm1
                y1(m-ipp) = Y1(m-ipp) - D(m-ipp)*Y1(m+1-ipp)
                y2(m-ipp) = Y2(m-ipp) - D(m-ipp)*Y2(m+1-ipp)
                y3(m-ipp) = Y3(m-ipp) - D(m-ipp)*Y3(m+1-ipp)
            end do
            if (k2k3k4 == 0) cycle
            if (n == l1) then
                i = lint1 + kint1
                xx = x - TCOS(i)
                y1(:m) = xx*Y1(:m) + W1(:m)
                lint1 = lint1 + 1
                l1 = (lint1*if1)/if2
            end if
            if (n == l2) then
                i = lint2 + kint2
                xx = x - TCOS(i)
                y2(:m) = xx*Y2(:m) + W2(:m)
                lint2 = lint2 + 1
                l2 = (lint2*if1)/if3
            end if
            if (n /= l3) cycle
            i = lint3 + kint3
            xx = x - TCOS(i)
            y3(:m) = xx*Y3(:m) + W3(:m)
            lint3 = lint3 + 1
            l3 = (lint3*if1)/if4
        end do

    end subroutine tri3
    !
    !*****************************************************************************************
    !
    subroutine merge_rename(tcos, i1, m1, i2, m2, i3)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in) :: i1
        integer (ip), intent (in) :: m1
        integer (ip), intent (in) :: i2
        integer (ip), intent (in) :: m2
        integer (ip), intent (in) :: i3
        real (wp),    intent (in out) :: tcos(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: j11, j3, j1, j2, j, l, k, m
        real (wp)    :: x, y
        !-----------------------------------------------
        !
        !     THIS SUBROUTINE MERGES TWO ASCENDING STRINGS OF NUMBERS IN THE
        !     ARRAY TCOS.  THE FIRST STRING IS OF LENGTH M1 AND STARTS AT
        !     TCOS(I1+1).  THE SECOND STRING IS OF LENGTH M2 AND STARTS AT
        !     TCOS(I2+1).  THE MERGED STRING GOES INTO TCOS(I3+1).
        !
        !
        j1 = 1
        j2 = 1
        j = i3
        if (m1 == 0) go to 107
        if (m2 == 0) go to 104
101 continue
    j11 = j1
    j3 = max(m1, j11)
    do j1 = j11, j3
        j = j + 1
        l = j1 + i1
        x = TCOS(l)
        l = j2 + i2
        y = TCOS(l)
        if (x - y > 0.) go to 103
        tcos(j) = x
    end do
    go to 106
103 continue
    tcos(j) = y
    j2 = j2 + 1
    if (j2 <= m2) go to 101
    if (j1 > m1) go to 109
104 continue
    k = j - j1 + 1
    do j = j1, m1
        m = k + j
        l = j + i1
        tcos(m) = TCOS(l)
    end do
    go to 109
106 continue
    if (j2 > m2) go to 109
107 continue
    k = j - j2 + 1
    do j = j2, m2
        m = k + j
        l = j + i2
        tcos(m) = TCOS(l)
    end do
109 continue

end subroutine merge_rename


end module module_gnbnaux
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! OCTOBER   1980    CHANGED SEVERAL DIVIDES OF FLOATING INTEGERS
!                   TO INTEGER DIVIDES TO ACCOMODATE CRAY-1 ARITHMETIC.
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
!-----------------------------------------------------------------------
