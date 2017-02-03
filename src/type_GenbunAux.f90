!
!     file type_GenbunAux.f90
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
! Package gnbnaux
!
! Latest revision        April 2016
!
! Purpose                to provide auxiliary routines for fishpack
!                        entries genbun and poistg.
!
! Usage                  There are no user entries in this package.
!                        the routines in this package are not intended
!                        to be called by users, but rather by routines
!                        in packages genbun and poistg.
!
! Special conditions     None
!
! I/O                    None
!
! Precision              64-bit double precision
!
!
! Language               Fortran 2008
!
! History                * Written in 1979 by Roland Sweet of NCAR'S
!                        scientific computing division. Made available
!                        on NCAR's public libraries in January, 1980.
!                        * Revised by John Adams in June 2004 incorporating
!                        Fortran 90 features
!                        * Revised by Jon Lo Kim Lin in 2016 incorporating
!                        Fortran 2008 features
!
! Portability            Fortran 2008
!
module type_GenbunAux

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        PI

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: GenbunAux

    
    type, public :: GenbunAux
    contains
        !--------------------------------------------------
        ! Type-bound procedures
        !--------------------------------------------------
        procedure, nopass, public :: cosgen
        procedure, nopass, public :: merger
        procedure, nopass, public :: trix
        procedure, nopass, public :: tri3
        procedure, nopass, private :: swap_real
        procedure, nopass, private :: swap_complex
        procedure, nopass, private :: swap_integer
        procedure, nopass, private :: swap_character
        !--------------------------------------------------
        ! Generic type-bound procedures
        !--------------------------------------------------
        generic, public :: swap => swap_real, &
            swap_complex, &
            swap_integer, &
            swap_character
        !--------------------------------------------------
    end type GenbunAux

    !---------------------------------------------------------------
    ! Parameters confined to the module
    !---------------------------------------------------------------
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    !---------------------------------------------------------------

contains

    pure subroutine cosgen(n, ijump, fnum, fden, a)
        !
        ! Purpose:
        !
        !     this subroutine computes required cosine values in ascending
        !     order.  when ijump >= 1 the routine computes values
        !
        !        2*cos(j*pi/l) , j=1, 2, ..., l and j /= 0(mod n/ijump+1)
        !
        !     where l = ijump*(n/ijump+1).
        !
        !
        !     when ijump = 1 it computes
        !
        !            2*cos((j-fnum)*pi/(n+fden)) ,  j=1, 2, ... , n
        !
        !     where
        !        fnum = 0.5, fden = 0.0, for regular reduction values
        !        fnum = 0.0, fden = 1.0, for b-r and c-r when istag = 1
        !        fnum = 0.0, fden = 0.5, for b-r and c-r when istag = 2
        !        fnum = 0.5, fden = 0.5, for b-r and c-r when istag = 2
        !                                in poisn2 only.
        !
        !
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: ijump
        real(wp),    intent(in)  :: fnum
        real(wp),    intent(in)  :: fden
        real(wp),    intent(out) :: a(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: k3, k4, k, k1, k5, i, k2, np1
        real(wp)    :: pibyn, x, y
        !-----------------------------------------------

        if (n == 0) return

        select case (ijump)
            case (1)
                !
                ! ijump == 1
                !
                np1 = n + 1
                y = PI/(real(n, kind=wp) + fden)

                do i = 1, n
                    x = real(np1 - i, kind=wp) - fnum
                    a(i) = TWO * cos(x*y)
                end do

            case default
                !
                ! ijump /= 1
                !
                k3 = n/ijump + 1
                k4 = k3 - 1
                pibyn = PI/(n + ijump)
                do k = 1, ijump
                    k1 = (k - 1)*k3
                    k5 = (k - 1)*k4
                    do i = 1, k4
                        x = k1 + i
                        k2 = k5 + i
                        a(k2) = -TWO * cos(x*pibyn)
                    end do
                end do
        end select

    end subroutine cosgen

    subroutine trix(idegbr, idegcr, m, a, b, c, y, tcos, d, w)
        !
        ! Purpose:
        !
        !     subroutine to solve a system of linear equations where the
        !     coefficient matrix is a rational function in the matrix given by
        !     tridiagonal  ( . . . , a(i), b(i), c(i), . . . ).
        !
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in)     :: idegbr
        integer(ip), intent(in)     :: idegcr
        integer(ip), intent(in)     :: m
        real(wp),    intent(in)     :: a(:)
        real(wp),    intent(in)     :: b(:)
        real(wp),    intent(in)     :: c(:)
        real(wp),    intent(inout)  :: y(:)
        real(wp),    intent(in)     :: tcos(:)
        real(wp),    intent(inout)  :: d(:)
        real(wp),    intent(inout)  :: w(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: mm1, ifb, ifc, l, lint, k, i, ip
        real(wp)    :: x, xx, z
        !-----------------------------------------------

        mm1 = m - 1
        ifb = idegbr + 1
        ifc = idegcr + 1
        l = ifb/ifc
        lint = 1

        outer_loop: do k = 1, idegbr

            x = tcos(k)

            if (k == l) then
                i = idegbr + lint
                xx = x - tcos(i)
                w = y
                y = xx*y
            end if

            z = ONE/(b(1)-x)
            d(1) = c(1)*z
            y(1) = y(1)*z

            do i = 2, mm1
                z = ONE/(b(i)-x-a(i)*d(i-1))
                d(i) = c(i)*z
                y(i) = (y(i)-a(i)*y(i-1))*z
            end do

            z = b(m) - x - a(m)*d(mm1)

            if (z == ZERO) then
                y(m) = ZERO
            else
                y(m) = (y(m)-a(m)*y(mm1))/z
            end if

            do ip = 1, mm1
                y(m-ip) = y(m-ip) - d(m-ip)*y(m+1-ip)
            end do

            if (k /= l) cycle outer_loop

            y = y + w
            lint = lint + 1
            l = (lint*ifb)/ifc

        end do outer_loop

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
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: k(4)
        real(wp),    intent(in)     :: a(:)
        real(wp),    intent(in)     :: b(:)
        real(wp),    intent(in)     :: c(:)
        real(wp),    intent(inout) :: y1(:)
        real(wp),    intent(inout) :: y2(:)
        real(wp),    intent(inout) :: y3(:)
        real(wp),    intent(in)    :: tcos(:)
        real(wp),    intent(inout) :: d(:)
        real(wp),    intent(inout) :: w1(:)
        real(wp),    intent(inout) :: w2(:)
        real(wp),    intent(inout) :: w3(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: mm1, k1, k2, k3, k4, if1, if2, if3, if4, k2k3k4, l1, l2
        integer(ip) :: l3, lint1, lint2, lint3, kint1, kint2, kint3, n, i, ipp
        real(wp)    :: x, z, xx
        !-----------------------------------------------

        mm1 = m - 1
        k1 = k(1)
        k2 = k(2)
        k3 = k(3)
        k4 = k(4)
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

        outer_loop: do n = 1, k1

            x = tcos(n)

            if (k2k3k4 /= 0) then
                if (n == l1) w1 = y1
                if (n == l2) w2 = y2
                if (n == l3) w3 = y3
            end if

            z = ONE/(b(1)-x)
            d(1) = c(1)*z
            y1(1) = y1(1)*z
            y2(1) = y2(1)*z
            y3(1) = y3(1)*z

            do i = 2, m
                z = ONE/(b(i)-x-a(i)*d(i-1))
                d(i) = c(i)*z
                y1(i) = (y1(i)-a(i)*y1(i-1))*z
                y2(i) = (y2(i)-a(i)*y2(i-1))*z
                y3(i) = (y3(i)-a(i)*y3(i-1))*z
            end do

            do ipp = 1, mm1
                y1(m-ipp) = y1(m-ipp) - d(m-ipp)*y1(m+1-ipp)
                y2(m-ipp) = y2(m-ipp) - d(m-ipp)*y2(m+1-ipp)
                y3(m-ipp) = y3(m-ipp) - d(m-ipp)*y3(m+1-ipp)
            end do

            if (k2k3k4 == 0) cycle outer_loop

            if (n == l1) then
                i = lint1 + kint1
                xx = x - tcos(i)
                y1 = xx*y1 + w1
                lint1 = lint1 + 1
                l1 = (lint1*if1)/if2
            end if

            if (n == l2) then
                i = lint2 + kint2
                xx = x - tcos(i)
                y2 = xx*y2 + w2
                lint2 = lint2 + 1
                l2 = (lint2*if1)/if3
            end if

            if (n /= l3) cycle outer_loop

            i = lint3 + kint3
            xx = x - tcos(i)
            y3 = xx*y3 + w3
            lint3 = lint3 + 1
            l3 = (lint3*if1)/if4

        end do outer_loop

    end subroutine tri3

    subroutine merger(tcos, i1, m1, i2, m2, i3)
        !
        ! Purpose:
        !
        !     this subroutine merges two ascending strings of numbers in the
        !     array tcos.  the first string is of length m1 and starts at
        !     tcos(i1+1).  the second string is of length m2 and starts at
        !     tcos(i2+1).  the merged string goes into tcos(i3+1).
        !
        !
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in)     :: i1
        integer(ip), intent(in)     :: m1
        integer(ip), intent(in)     :: i2
        integer(ip), intent(in)     :: m2
        integer(ip), intent(in)     :: i3
        real(wp),    intent(inout) :: tcos(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: j11, j3, j1, j2, j, l, k, m
        real(wp)    :: x, y
        !-----------------------------------------------

        j1 = 1
        j2 = 1
        j = i3

        if_construct: if (m1 /= 0) then
            if (m2 /= 0) then

                outer_loop: do

                    j11 = j1
                    j3 = max(m1, j11)

                    block_construct: block

                        do j1 = j11, j3
                            j = j + 1
                            l = j1 + i1
                            x = tcos(l)
                            l = j2 + i2
                            y = tcos(l)
                            if (x - y > ZERO) exit block_construct
                            tcos(j) = x
                        end do

                        if (j2 > m2) return

                        exit if_construct

                    end block block_construct

                    tcos(j) = y
                    j2 = j2 + 1

                    if (j2 > m2) exit outer_loop

                end do outer_loop

                if (j1 > m1) return

            end if

            k = j - j1 + 1

            do j = j1, m1
                m = k + j
                l = j + i1
                tcos(m) = tcos(l)
            end do

            return

        end if if_construct

        k = j - j2 + 1

        do j = j2, m2
            m = k + j
            l = j + i2
            tcos(m) = tcos(l)
        end do

    end subroutine merger

    pure subroutine swap_real(a, b)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        real(wp), intent(inout) :: a, b
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        real(wp) temp
        !-----------------------------------------------

        temp = a ; a = b ; b = temp

    end subroutine swap_real

    pure subroutine swap_complex(a, b)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        complex(wp), intent(inout) :: a, b
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        complex(wp) temp
        !-----------------------------------------------

        temp = a ; a = b ; b = temp

    end subroutine swap_complex

    pure subroutine swap_integer(a, b)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(inout) :: a, b
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) temp
        !-----------------------------------------------

        temp = a ; a = b ; b = temp

    end subroutine swap_integer

    pure subroutine swap_character(a, b)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        character, intent(inout) :: a, b
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        character temp
        !-----------------------------------------------

        temp = a
        a = b
        b = temp

    end subroutine swap_character

end module type_GenbunAux
!
! REVISION HISTORY
!
! September 1973    Version 1
! April     1976    Version 2
! January   1978    Version 3
! December  1979    Version 3.1
! October   1980    Changed several divides of floating integers
!                   to integer divides to accomodate cray-1 arithmetic.
! February  1985    Documentation upgrade
! November  1988    Version 3.2, FORTRAN 77 changes
! May       2016    Fortran 2008 changes
