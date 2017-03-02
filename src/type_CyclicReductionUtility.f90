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
! Purpose                to provide auxiliary routines for fishpack
!                        entries genbun and poistg.
!
! Usage                  There are no user entries in this package.
!                        the routines in this package are not intended
!                        to be called by users, but rather by routines
!                        in genbun and poistg.
!
! History                * Written in 1979 by Roland Sweet of NCAR'S
!                        scientific computing division. Made available
!                        on NCAR's public libraries in January, 1980.
!                        * Revised by John Adams in June 2004 incorporating
!                        Fortran 90 features
!
module type_CyclicReductionUtility

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        PI

        ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: CyclicReductionUtility
    
    type, public :: CyclicReductionUtility
    contains
        ! Public type-bound procedures
        procedure, nopass, public :: generate_cosines
        procedure, nopass, public :: merge_cosines
        procedure, nopass, public :: solve_tridiag
        procedure, nopass, public :: solve_tridiag3
        ! Private type-bound procedures
        procedure, nopass, private :: swap_real
        procedure, nopass, private :: swap_complex
        procedure, nopass, private :: swap_integer
        procedure, nopass, private :: swap_character
        ! Generic type-bound procedures
        generic, public :: swap => swap_real, &
            swap_complex, &
            swap_integer, &
            swap_character
    end type CyclicReductionUtility

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: TWO = 2.0_wp

contains

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
    pure subroutine generate_cosines(n, ijump, fnum, fden, a)

        ! Dummy arguments
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: ijump
        real(wp),    intent(in)  :: fnum
        real(wp),    intent(in)  :: fden
        real(wp),    intent(out) :: a(:)

        ! Local variables
        integer(ip) :: k3, k4, k, k1, k5, i, k2, np1
        real(wp)    :: pibyn, x, y

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

    end subroutine generate_cosines

    !
    ! Purpose:
    !
    !     subroutine to solve a system of linear equations where the
    !     coefficient matrix is a rational function in the matrix given by
    !     tridiagonal  ( . . . , a(i), b(i), c(i), . . . ).
    !
    subroutine solve_tridiag(idegbr, idegcr, m, a, b, c, y, tcos, d, w)

        ! Dummy arguments
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

        ! Local variables
        integer(ip) :: ifb, ifc, j, lint, k, i
        real(wp)    :: x, xx, temp

        ifb = idegbr + 1
        ifc = idegcr + 1
        j = ifb/ifc
        lint = 1

        outer_loop: do k = 1, idegbr

            x = tcos(k)

            if (k == j) then
                i = idegbr + lint
                xx = x - tcos(i)
                w = y
                y = xx*y
            end if

            temp = b(1) - x
            d(1) = c(1)/temp
            y(1) = y(1)/temp

            do i = 2, m - 1
                temp = b(i) - x - a(i) * d(i-1)
                d(i) = c(i)/temp
                y(i) = (y(i) - a(i) * y(i-1))/temp
            end do

            temp = b(m) - x - a(m) * d(m - 1)

            if (temp == ZERO) then
                y(m) = ZERO
            else
                y(m) = (y(m) - a(m) * y(m - 1))/temp
            end if

            do i = m - 1, 1, -1
                y(i) = y(i) - d(i) * y(i + 1)
            end do

            if (k /= j) cycle outer_loop

            y = y + w
            lint = lint + 1
            j = (lint*ifb)/ifc

        end do outer_loop

    end subroutine solve_tridiag

    !
    ! Purpose:
    !
    ! subroutine to solve three linear systems whose common coefficient
    ! matrix is a rational function in the matrix given by
    !
    !  tridiagonal (..., a(i), b(i), c(i), ...)
    !
    subroutine solve_tridiag3(m, a, b, c, k, y1, y2, y3, tcos, d, w1, w2, w3)

        ! Dummy arguments
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

        ! Local variables
        integer(ip) :: mm1, k1, k2, k3, k4, if1, if2, if3, if4, k2k3k4, l1, l2
        integer(ip) :: l3, lint1, lint2, lint3, kint1, kint2, kint3, n, i
        real(wp)    :: x, temp, xx

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

            temp = b(1) - x
            d(1) = c(1)/temp
            y1(1) = y1(1)/temp
            y2(1) = y2(1)/temp
            y3(1) = y3(1)/temp

            do i = 2, m
                temp = b(i)-x-a(i)*d(i-1)
                d(i) = c(i)/temp
                y1(i) = (y1(i)-a(i)*y1(i-1))/temp
                y2(i) = (y2(i)-a(i)*y2(i-1))/temp
                y3(i) = (y3(i)-a(i)*y3(i-1))/temp
            end do

            do i = m - 1, 1, -1
                y1(i) = y1(i) - d(i) * y1(i + 1)
                y2(i) = y2(i) - d(i) * y2(i + 1)
                y3(i) = y3(i) - d(i) * y3(i + 1)
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

    end subroutine solve_tridiag3

    ! Purpose:
    !
    !     this subroutine merges two ascending strings of numbers in the
    !     array tcos.  the first string is of length m1 and starts at
    !     tcos(i1+1).  the second string is of length m2 and starts at
    !     tcos(i2+1).  the merged string goes into tcos(i3+1).
    !
    subroutine merge_cosines(tcos, i1, m1, i2, m2, i3)

        ! Dummy arguments
        integer(ip), intent(in)    :: i1
        integer(ip), intent(in)    :: m1
        integer(ip), intent(in)    :: i2
        integer(ip), intent(in)    :: m2
        integer(ip), intent(in)    :: i3
        real(wp),    intent(inout) :: tcos(:)

        ! Local variables
        integer(ip) :: j1, j2, j3, istop

        if (m1 == 0 .and. m2 == 0) then
            return
        else if (m1 == 0 .and. m2 /= 0) then
            istop = m2
            tcos(i3+1:istop) = tcos(i2+1:istop)
        else if (m1 /= 0 .and. m2 == 0) then
            istop = m1
            tcos(i3+1:istop) = tcos(i1+1:istop)
        else
            j1 = 1
            j2 = 1
            j3 = 1
            do
                if (tcos(i1+j1)  <=  tcos(i2+j2)) then
                    tcos(i3+j3) = tcos(i1+j1)
                    j1 = j1+1
                    if (j1  >  m1) then
                        istop = m2-j2+1
                        tcos(i3+j3+1:istop) = tcos(i2+j2:istop)
                        return
                    end if
                else
                    tcos(i3+j3) = tcos(i2+j2)
                    j2 = j2+1
                    if (j2  >  m2) then
                        istop = m1-j1+1
                        tcos(i3+j3+1:istop) = tcos(i1+j1:istop)
                        return
                    end if
                end if
                j3 = j3+1
            end do
        end if

    end subroutine merge_cosines

    pure subroutine swap_real(a, b)

        ! Dummy arguments
        real(wp), intent(inout) :: a, b

        ! Local variables
        real(wp) temp

        temp = a ; a = b ; b = temp

    end subroutine swap_real

    pure subroutine swap_complex(a, b)

        ! Dummy arguments
        complex(wp), intent(inout) :: a, b

        ! Local variables
        complex(wp) temp

        temp = a ; a = b ; b = temp

    end subroutine swap_complex

    pure subroutine swap_integer(a, b)

        ! Dummy arguments
        integer(ip), intent(inout) :: a, b

        ! Local variables
        integer(ip) temp

        temp = a ; a = b ; b = temp

    end subroutine swap_integer

    pure subroutine swap_character(a, b)

        ! Dummy arguments
        character, intent(inout) :: a, b

        ! Local variables
        character temp

        temp = a
        a = b
        b = temp

    end subroutine swap_character

end module type_CyclicReductionUtility
