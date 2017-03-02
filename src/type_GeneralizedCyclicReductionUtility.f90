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
module type_GeneralizedCyclicReductionUtility

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        MACHINE_EPSILON ! Machine epsilon

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: GeneralizedCyclicReductionUtility
    public :: ppsgf
    public :: ppspf
    public :: psgf
    public :: comf_interface

    ! Parameters confined to the module
    real(wp),    parameter :: ZERO = 0.0_wp
    real(wp),    parameter :: ONE = 1.0_wp
    real(wp),    parameter :: TWO = 2.0_wp
    integer(ip), parameter :: SIZE_OF_WORKSPACE_INDICES = 6

    type, public :: GeneralizedCyclicReductionUtility
        ! Type components
        integer(ip) :: indices(SIZE_OF_WORKSPACE_INDICES)
        integer(ip) :: nl, npp, k, nm, ncmplx, ik
        real(wp)    :: cnv
    contains
        ! Public type-bound procedures
        procedure,         public :: blktrii
        procedure, nopass, public :: check_input_arguments
        procedure, nopass, public :: ppsgf
        procedure, nopass, public :: ppspf
        procedure, nopass, public :: psgf
        ! Private type-bound procedures
        procedure, private :: blktri_lower_routine
        procedure, private :: bsrh
        procedure, private :: compute_roots_of_b_polynomials
        procedure, private :: compute_eigen_values
        procedure, private :: tevls
        procedure, private :: compute_index_a_coeff
        procedure, private :: compute_index_b_coeff
        procedure, private :: compute_index_c_coeff
    end type GeneralizedCyclicReductionUtility

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

    pure subroutine check_input_arguments(n, m, idimy, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: n, m, idimy
        integer(ip), intent(out) :: ierror

        if (m < 5) then
            ierror = 1
        else if (n < 3) then
            ierror = 2
        else if (idimy < m) then
            ierror = 3
        else
            ierror = 0
        end if

    end subroutine check_input_arguments

    subroutine blktrii(self, iflg, np, n, an, bn, cn, mp, m, am, bm, cm, &
        idimy, y, ierror, w, wc)

        ! Dummy arguments
        class(GeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(in)    :: iflg
        integer(ip),      intent(in)    :: np
        integer(ip),      intent(in)    :: n
        integer(ip),      intent(in)    :: mp
        integer(ip),      intent(in)    :: m
        integer(ip),      intent(in)    :: idimy
        integer(ip),      intent(out)   :: ierror
        real(wp),         intent(in)    :: an(:)
        real(wp),         intent(in)    :: bn(:)
        real(wp),         intent(in)    :: cn(:)
        real(wp),         intent(in)    :: am(:)
        real(wp),         intent(in)    :: bm(:)
        real(wp),         intent(in)    :: cm(:)
        real(wp),         intent(inout) :: y(:,:)
        real(wp),         intent(inout) :: w(:)
        complex(wp),      intent(inout) :: wc(:)

        ! Local variables
        integer(ip) :: nh, iwah, iwbh

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ik => self%ik &
            )

            ! test m and n for the proper form
            nm = n

            ! Check again for solvers which call blktrii directly
            call check_input_arguments(nm, m, idimy, ierror)

            ! Check error flag
            if (ierror /= 0) return

            associate( &
                iw1 => self%indices(1), &
                iw2 => self%indices(2), &
                iw3 => self%indices(3), &
                iww => self%indices(4), &
                iwu => self%indices(5), &
                iwd => self%indices(6), &
                nl => self%nl &
                )

                select case (iflg)
                    case (0)
                        nh = n
                        npp = np
                        if (npp /= 0)  nh = nh + 1
                        ik = 4
                        k = 2

                        do
                            if (nh <= ik) exit
                            ik = 2*ik
                            k = k + 1
                        end do

                        nl = ik
                        ik = 2*ik
                        nl = nl - 1
                        iwah = (k - 2)*ik + k + 5

                        if (npp == 0) then
                            iwbh = iwah + 2*nm
                            iw1 = iwbh
                            nm = n - 1
                        else
                            iw1 = iwah
                            iwbh = iw1 + nm
                        end if
                        !
                        ! Set workspace indices
                        !
                        iw2 = iw1 + m
                        iw3 = iw2 + m
                        iwd = iw3 + m
                        iww = iwd + m
                        iwu = iww + m

                        call self%compute_roots_of_b_polynomials(nl, ierror, an, bn, cn, w, wc, w(iwah:), w(iwbh:))

                    case default

                        ! *** Important to reset nm for np = 0
                        if (npp == 0) nm = n - 1

                        select case (mp)
                            case (0)
                                call self%blktri_lower_routine(an, cn, m, am, bm, cm, y, w, wc, &
                                    w(iw1:), w(iw2:), w(iw3:), w(iwd:), w(iww:), w(iwu:), wc(iw1:), &
                                    wc(iw2:), wc(iw3:), prodp, cprodp)
                            case default
                                call self%blktri_lower_routine(an, cn, m, am, bm, cm, y, w, wc, &
                                    w(iw1:), w(iw2:), w(iw3:), w(iwd:), w(iww:), w(iwu:), wc(iw1:), &
                                    wc(iw2:), wc(iw3:), prod, cprod)
                        end select
                end select

            end associate

        end associate common_variables

    end subroutine blktrii

    ! Purpose:
    !
    ! Solves the linear system
    !
    ! b  contains the roots of all the b polynomials
    ! w1, w2, w3, wd, ww, wu  are all working arrays
    ! prdct  is either prodp or prod depending on whether the boundary
    ! conditions in the m direction are periodic or not
    ! cprdct is either cprodp or cprod which are the complex versions
    ! of prodp and prod. these are called in the event that some
    ! of the roots of the b sub p polynomial are complex
    !
    subroutine blktri_lower_routine(self, an, cn, m, am, bm, cm, y, b, bc, &
        w1, w2, w3, wd, ww, wu, cw1, cw2, cw3, prdct, cprdct)

        ! Dummy arguments
        class(GeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip), intent(in)     :: m
        real(wp),    intent(in)     :: an(:)
        real(wp),    intent(in)     :: cn(:)
        real(wp),    intent(in)     :: am(:)
        real(wp),    intent(in)     :: bm(:)
        real(wp),    intent(in)     :: cm(:)
        real(wp),    intent(inout)  :: y(:,:)
        real(wp),    intent(in)     :: b(*)
        real(wp),    intent(inout)  :: w1(*)
        real(wp),    intent(inout)  :: w2(*)
        real(wp),    intent(inout)  :: w3(*)
        real(wp),    intent(in)     :: wd(*)
        real(wp),    intent(in)     :: ww(*)
        real(wp),    intent(in)     :: wu(*)
        complex(wp), intent(in)     :: bc(*)
        complex(wp), intent(in)     :: cw1(*)
        complex(wp), intent(in)     :: cw2(*)
        complex(wp), intent(in)     :: cw3(*)
        external :: prdct, cprdct

        ! Local variables
        integer(ip) :: kdo, l, ir, i2, i1, i3, i4, irm1, im2, nm2, im3, nm3
        integer(ip) :: im1, nm1, i0, iif, i, ipi1, ipi2, ipi3, idxc, nc, idxa, na, ip2
        integer(ip) :: np2, ip1, np1, ip3, np3, iz, nz, izr, ll, ifd, iip, np
        integer(ip) :: imi1, imi2
        real(wp)    :: dummy_variable

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv &
            )
            ! begin reduction phase
            kdo = k - 1
            do l = 1, kdo
                ir = l - 1
                i2 = 2**ir
                i1 = i2/2
                i3 = i2 + i1
                i4 = i2 + i2
                irm1 = ir - 1
                call self%compute_index_b_coeff(i2, ir, im2, nm2)
                call self%compute_index_b_coeff(i1, irm1, im3, nm3)
                call self%compute_index_b_coeff(i3, irm1, im1, nm1)

                i0 = 0

                call prdct(nm2, b(im2), nm3, b(im3), nm1, b(im1), i0, dummy_variable, &
                    y(1, i2), w3, m, am, bm, cm, wd, ww, wu)

                iif = 2**k

                do i = i4, iif, i4

                    if (i > nm) cycle

                    ipi1 = i + i1
                    ipi2 = i + i2
                    ipi3 = i + i3

                    call self%compute_index_c_coeff(i, ir, idxc, nc)

                    if (i >= iif) cycle

                    call self%compute_index_a_coeff(i, ir, idxa, na)
                    call self%compute_index_b_coeff(i - i1, irm1, im1, nm1)
                    call self%compute_index_b_coeff(ipi2, ir, ip2, np2)
                    call self%compute_index_b_coeff(ipi1, irm1, ip1, np1)
                    call self%compute_index_b_coeff(ipi3, irm1, ip3, np3)
                    call prdct(nm1, b(im1), 0, dummy_variable, 0, dummy_variable, na, an(idxa), w3, &
                        w1, m, am, bm, cm, wd, ww, wu)

                    if (ipi2 > nm) then
                        w3(:m) = ZERO
                        w2(:m) = ZERO
                    else
                        call prdct(np2, b(ip2), np1, b(ip1), np3, b(ip3), 0, dummy_variable, &
                            y(1, ipi2), w3, m, am, bm, cm, wd, ww, wu)
                        call prdct(np1, b(ip1), 0, dummy_variable, 0, dummy_variable, nc, cn(idxc), w3, &
                            w2, m, am, bm, cm, wd, ww, wu)
                    end if
                    y(:m, i) = w1(:m) + w2(:m) + y(:m, i)
                end do
            end do

            if (npp == 0) then
                iif = 2**k
                i = iif/2
                i1 = i/2

                call self%compute_index_b_coeff(i - i1, k - 2, im1, nm1)
                call self%compute_index_b_coeff(i + i1, k - 2, ip1, np1)
                call self%compute_index_b_coeff(i, k - 1, iz, nz)
                call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dummy_variable, &
                    y(1, i), w1, m, am, bm, cm, wd, ww, wu)

                izr = i
                w2(:m) = w1(:m)

                do ll = 2, k
                    l = k - ll + 1
                    ir = l - 1
                    i2 = 2**ir
                    i1 = i2/2
                    i = i2

                    call self%compute_index_c_coeff(i, ir, idxc, nc)
                    call self%compute_index_b_coeff(i, ir, iz, nz)
                    call self%compute_index_b_coeff(i - i1, ir - 1, im1, nm1)
                    call self%compute_index_b_coeff(i + i1, ir - 1, ip1, np1)
                    call prdct(np1, b(ip1), 0, dummy_variable, 0, dummy_variable, nc, cn(idxc), w1, &
                        w1, m, am, bm, cm, wd, ww, wu)

                    w1(:m) = y(:m, i) + w1(:m)

                    call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, &
                        dummy_variable, w1, w1, m, am, bm, cm, wd, ww, wu)
                end do

                outer_loop: do ll = 2, k
                    l = k - ll + 1
                    ir = l - 1
                    i2 = 2**ir
                    i1 = i2/2
                    i4 = i2 + i2
                    ifd = iif - i2
                    inner_loop: do i = i2, ifd, i4

                        if (i - i2 - izr /= 0) cycle inner_loop

                        if (i > nm) cycle outer_loop

                        call self%compute_index_a_coeff(i, ir, idxa, na)
                        call self%compute_index_b_coeff(i, ir, iz, nz)
                        call self%compute_index_b_coeff(i - i1, ir - 1, im1, nm1)
                        call self%compute_index_b_coeff(i + i1, ir - 1, ip1, np1)
                        call prdct(nm1, b(im1), 0, dummy_variable, 0, dummy_variable, na, an(idxa), &
                            w2, w2, m, am, bm, cm, wd, ww, wu)

                        w2(:m) = y(:m, i) + w2(:m)

                        call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dummy_variable, &
                            w2, w2, m, am, bm, cm, wd, ww, wu)

                        izr = i

                        if (i == nm) exit outer_loop

                    end do inner_loop
                end do outer_loop

                y(:m, nm+1) = y(:m, nm+1) - cn(nm+1)*w1(:m) - an(nm+1)*w2(:m)

                call self%compute_index_b_coeff(iif/2, k - 1, im1, nm1)
                call self%compute_index_b_coeff(iif, k - 1, iip, np)

                if (ncmplx /= 0) then
                    call cprdct(nm + 1, bc(iip), nm1, b(im1), 0, dummy_variable, 0, dummy_variable, &
                        y(1, nm+1), y(1, nm+1), m, am, bm, cm, cw1, cw2, cw3)
                else
                    call prdct(nm + 1, b(iip), nm1, b(im1), 0, dummy_variable, 0, dummy_variable, &
                        y(1, nm+1), y(1, nm+1), m, am, bm, cm, wd, ww, wu)
                end if

                w1(:m) = an(1)*y(:m, nm+1)
                w2(:m) = cn(nm)*y(:m, nm+1)
                y(:m, 1) = y(:m, 1) - w1(:m)
                y(:m, nm) = y(:m, nm) - w2(:m)

                do l = 1, kdo
                    ir = l - 1
                    i2 = 2**ir
                    i4 = i2 + i2
                    i1 = i2/2
                    i = i4

                    call self%compute_index_a_coeff(i, ir, idxa, na)
                    call self%compute_index_b_coeff(i - i2, ir, im2, nm2)
                    call self%compute_index_b_coeff(i - i2 - i1, ir - 1, im3, nm3)
                    call self%compute_index_b_coeff(i - i1, ir - 1, im1, nm1)
                    call prdct(nm2, b(im2), nm3, b(im3), nm1, b(im1), 0, dummy_variable, &
                        w1, w1, m, am, bm, cm, wd, ww, wu)
                    call prdct(nm1, b(im1), 0, dummy_variable, 0, dummy_variable, na, an(idxa), w1, &
                        w1, m, am, bm, cm, wd, ww, wu)

                    y(:m, i) = y(:m, i) - w1(:m)
                end do

                izr = nm

                loop_131: do l = 1, kdo
                    ir = l - 1
                    i2 = 2**ir
                    i1 = i2/2
                    i3 = i2 + i1
                    i4 = i2 + i2
                    irm1 = ir - 1
                    loop_132: do i = i4, iif, i4
                        ipi1 = i + i1
                        ipi2 = i + i2
                        ipi3 = i + i3

                        if (ipi2 /= izr) then

                            if (i /= izr) cycle loop_132

                            cycle loop_131

                        end if

                        call self%compute_index_c_coeff(i, ir, idxc, nc)
                        call self%compute_index_b_coeff(ipi2, ir, ip2, np2)
                        call self%compute_index_b_coeff(ipi1, irm1, ip1, np1)
                        call self%compute_index_b_coeff(ipi3, irm1, ip3, np3)

                        call prdct(np2, b(ip2), np1, b(ip1), np3, b(ip3), 0, &
                            dummy_variable, w2, w2, m, am, bm, cm, wd, ww, wu)

                        call prdct(np1, b(ip1), 0, dummy_variable, 0, dummy_variable, nc, cn(idxc), &
                            w2, w2, m, am, bm, cm, wd, ww, wu)

                        y(:m, i) = y(:m, i) - w2(:m)
                        izr = i

                        cycle loop_131
                    end do loop_132
                end do loop_131
            end if
            !
            ! begin back substitution phase
            !
            do ll = 1, k
                l = k - ll + 1
                ir = l - 1
                irm1 = ir - 1
                i2 = 2**ir
                i1 = i2/2
                i4 = i2 + i2
                ifd = iif - i2
                inner_back_sub: do i = i2, ifd, i4

                    if (i > nm) cycle inner_back_sub

                    imi1 = i - i1
                    imi2 = i - i2
                    ipi1 = i + i1
                    ipi2 = i + i2

                    call self%compute_index_a_coeff(i, ir, idxa, na)
                    call self%compute_index_c_coeff(i, ir, idxc, nc)
                    call self%compute_index_b_coeff(i, ir, iz, nz)
                    call self%compute_index_b_coeff(imi1, irm1, im1, nm1)
                    call self%compute_index_b_coeff(ipi1, irm1, ip1, np1)

                    if (i <= i2) then
                        w1(:m) = ZERO
                    else
                        call prdct(nm1, b(im1), 0, dummy_variable, 0, dummy_variable, na, an(idxa), &
                            y(1,imi2), w1, m, am, bm, cm, wd, ww, wu)
                    end if

                    if (ipi2 > nm) then
                        w2(:m) = ZERO
                    else
                        call prdct(np1, b(ip1), 0, dummy_variable, 0, dummy_variable, nc, cn(idxc), &
                            y(1,ipi2), w2, m, am, bm, cm, wd, ww, wu)
                    end if

                    w1(:m) = y(:m, i) + w1(:m) + w2(:m)

                    call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dummy_variable, &
                        w1, y(1, i), m, am, bm, cm, wd, ww, wu)
                end do inner_back_sub
            end do
        end associate common_variables

    end subroutine blktri_lower_routine

    function bsrh(self, xll, xrr, iz, c, a, bh, sgn, f) &
        result (return_value)

        ! Dummy arguments
        class(GeneralizedCyclicReductionUtility), intent(inout) :: self
        real(wp),         intent(in)    :: xll
        real(wp),         intent(in)    :: xrr
        integer(ip),      intent(in)    :: iz
        real(wp),         intent(in)    :: c(:)
        real(wp),         intent(in)    :: a(:)
        real(wp),         intent(in)    :: bh(:)
        real(wp),         intent(in)    :: sgn
        procedure(comf_interface)       :: f
        real(wp)                        :: return_value

        ! Local variables
        real(wp) :: r1, xl, xr, dx, x

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv &
            )

            xl = xll
            xr = xrr
            dx = abs(xr - xl)/2
            x = (xl + xr)/2
            r1 = sgn * f(x, iz, c, a, bh)

            if (r1 >= ZERO) then
                if (r1 == ZERO) then
                    return_value = (xl + xr)/2
                    return
                end if
                xr = x
            else
                xl = x
            end if

            dx = dx/2

            do
                if (dx - cnv > ZERO) exit
                x = (xl + xr)/2
                r1 = sgn * f(x, iz, c, a, bh)
                if (r1 >= ZERO) then
                    if (r1 == ZERO) then
                        return_value = (xl + xr)/2
                        return
                    end if
                    xr = x
                else
                    xl = x
                end if
                dx = dx/2
            end do

            return_value = (xl + xr)/2
        end associate common_variables

    end function bsrh

    ! Purpose:
    !
    ! Computes the roots of the b polynomials using subroutine
    ! tevls which is a modification the eispack program tqlrat.
    ! ierror is set to 4 if either tevls fails or if a(j+1)*c(j) is
    ! less than zero for some j.  ah, bh are temporary work arrays.
    !
    subroutine compute_roots_of_b_polynomials(self, n, ierror, an, bn, cn, b, bc, ah, bh)

        ! Dummy arguments
        class(GeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(in)    :: n
        integer(ip),      intent(out)   :: ierror
        real(wp),         intent(in)    :: an(:)
        real(wp),         intent(in)    :: bn(:)
        real(wp),         intent(in)    :: cn(:)
        real(wp),         intent(inout) :: b(:)
        real(wp),         intent(inout) :: ah(:)
        real(wp),         intent(inout) :: bh(:)
        complex(wp),      intent(inout) :: bc(:)

        ! Local variables

        integer(ip)  :: j, iif, kdo, l, ir, i2, i4, ipl, ifd, i, ib, nb, js, jf
        integer(ip)  :: ls, lh, nmp, l1, l2, j2, j1, n2m2
        real(wp)     :: bnorm, arg, d1, d2, d3


        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv &
            )

            bnorm = abs(bn(1))

            do j = 2, nm
                bnorm = max(bnorm, abs(bn(j)))
                arg = an(j)*cn(j-1)

                if (arg < ZERO) then
                    ierror = 5
                    return
                end if

                b(j) = sign(sqrt(arg), an(j))
            end do

            cnv = MACHINE_EPSILON*bnorm
            iif = 2**k
            kdo = k - 1

            outer_loop: do l = 1, kdo
                ir = l - 1
                i2 = 2**ir
                i4 = i2 + i2
                ipl = i4 - 1
                ifd = iif - i4
                do i = i4, ifd, i4
                    call self%compute_index_b_coeff(i, l, ib, nb)

                    if (nb <= 0) cycle outer_loop

                    js = i - ipl
                    jf = js + nb - 1
                    ls = 0
                    bh(:jf-js+1) = bn(js:jf)
                    ah(:jf-js+1) = b(js:jf)

                    associate( order => nb )
                        call self%tevls(bh(1:order), ah(1:order), ierror)
                    end associate

                    if (ierror /= 0) then
                        ierror = 4
                        return
                    end if

                    lh = ib - 1
                    if (nb > 0) then
                        b(lh+1:nb+lh) = -bh(:nb)
                        lh = nb + lh
                    end if
                end do
            end do outer_loop

            b(:nm) = -bn(:nm)

            if (npp == 0) then
                nmp = nm + 1
                nb = nm + nmp
                do j = 1, nb
                    l1 = mod(j - 1, nmp) + 1
                    l2 = mod(j + nm - 1, nmp) + 1
                    arg = an(l1)*cn(l2)

                    if (arg < ZERO) then
                        ierror = 5
                        return
                    end if

                    bh(j) = sign(sqrt(arg), (-an(l1)))
                    ah(j) = -bn(l1)
                end do

                associate( order => nb )
                    call self%tevls(ah(1:order), bh(1:order), ierror)
                end associate

                if (ierror /= 0) then
                    ierror = 4
                    return
                end if

                call self%compute_index_b_coeff(iif, k - 1, j2, lh)
                call self%compute_index_b_coeff(iif/2, k - 1, j1, lh)

                j2 = j2 + 1
                lh = j2
                n2m2 = j2 + 2*nm - 2

                do
                    if (j2 <= n2m2) exit

                    d1 = abs(b(j1)-b(j2-1))
                    d2 = abs(b(j1)-b(j2))
                    d3 = abs(b(j1)-b(j2+1))

                    if (d1 <= d2 .or. d3 <= d2) then
                        b(lh) = b(j2)
                        j2 = j2 + 1
                        lh = lh + 1
                    else
                        j2 = j2 + 1
                        j1 = j1 + 1
                    end if
                end do

                b(lh) = b(n2m2+1)

                call self%compute_index_b_coeff(iif, k - 1, j1, j2)

                j2 = j1 + nmp + nmp

                associate( order => nm + 1 )
                    call self%compute_eigen_values(ierror, an(1:order), cn(1:order), bc(j1:order), b(j1:order), b(j2:order))
                end associate

            end if

        end associate common_variables

    end subroutine compute_roots_of_b_polynomials

    pure subroutine cprod(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, yy, &
        m, a, b, c, d, w, y)
        !
        ! Purpose:
        !
        ! cprod applies a sequence of matrix operations to the vector x and
        ! stores the result in yy (complex case)
        !
        ! aa   array containing scalar multipliers of the vector x
        !
        ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
        !
        ! bd, bm1, bm2 are arrays containing roots of certian b polynomials
        !
        ! na is the length of the array aa
        !
        ! x, yy the matrix operations are applied to x and the result is yy
        !
        ! a, b, c  are arrays which contain the tridiagonal matrix
        !
        ! m  is the order of the matrix
        !
        ! d, w, y are working arrays
        !
        ! isgn  determines whether or not a change in sign is made
        !

        ! Dummy arguments

        integer(ip), intent(in)  :: nd
        integer(ip), intent(in)  :: nm1
        integer(ip), intent(in)  :: nm2
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: bm1(nm1)
        real(wp),    intent(in)  :: bm2(nm2)
        real(wp),    intent(in)  :: aa(na)
        real(wp),    intent(in)  :: x(m)
        real(wp),    intent(out) :: yy(m)
        real(wp),    intent(in)  :: a(m)
        real(wp),    intent(in)  :: b(m)
        real(wp),    intent(in)  :: c(m)
        complex(wp), intent(in)  :: bd(nd)
        complex(wp), intent(out) :: d(m)
        complex(wp), intent(out) :: w(m)
        complex(wp), intent(out) :: y(m)

        ! Local variables

        integer(ip) :: j, mm, id, m1, m2, ia, iflg, k
        real(wp)    :: rt
        complex(wp) :: crt, den, y1, y2


        y = cmplx(x, ZERO, kind=wp)

        mm = m - 1
        id = nd
        m1 = nm1
        m2 = nm2
        ia = na

        main_loop: do

            iflg = 0

            if (id > 0) then
                crt = bd(id)
                id = id - 1
                !
                ! begin solution to system
                !
                d(m) = a(m)/(b(m)-crt)
                w(m) = y(m)/(b(m)-crt)

                do j = 2, mm
                    k = m - j
                    den = b(k+1) - crt - c(k+1)*d(k+2)
                    d(k+1) = a(k+1)/den
                    w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
                end do

                den = b(1) - crt - c(1)*d(2)

                if (abs(den) /= ZERO) then
                    y(1) = (y(1)-c(1)*w(2))/den
                else
                    y(1) = cmplx(ONE, ZERO, kind=wp)
                end if

                y(2:m) = w(2:m) - d(2:m)*y(1:m-1)

            end if

            if (m1 > 0 .and. m2 > 0) then
                if (m1 <= 0) then
                    rt = bm2(m2)
                    m2 = m2 - 1
                else
                    if (m2 <= 0) then
                        rt = bm1(m1)
                        m1 = m1 - 1
                    else
                        if (abs(bm1(m1)) - abs(bm2(m2)) > ZERO) then
                            rt = bm1(m1)
                            m1 = m1 - 1
                        else
                            rt = bm2(m2)
                            m2 = m2 - 1
                        end if
                    end if
                end if
                y1 = (b(1)-rt)*y(1) + c(1)*y(2)
                if (mm >= 2) then
                    do j = 2, mm
                        y2 = a(j)*y(j-1) + (b(j)-rt)*y(j) + c(j)*y(j+1)
                        y(j-1) = y1
                        y1 = y2
                    end do
                end if
                y(m) = a(m)*y(m-1) + (b(m)-rt)*y(m)
                y(m-1) = y1
                iflg = 1
                cycle main_loop
            end if

            if (ia > 0) then
                rt = aa(ia)
                ia = ia - 1
                iflg = 1
                !
                ! scalar multiplication
                !
                y = rt*y
            end if

            if (iflg <= 0) exit main_loop

        end do main_loop

        yy = real(y, kind=wp)

    end subroutine cprod

    pure subroutine cprodp(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, yy, m, a, &
        b, c, d, u, y)
        !
        ! Purpose:
        !
        ! cprodp applies a sequence of matrix operations to the vector x and
        ! stores the result in yy       periodic boundary conditions
        ! and  complex  case
        !
        ! bd, bm1, bm2 are arrays containing roots of certian b polynomials
        ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
        ! aa   array containing scalar multipliers of the vector x
        ! na is the length of the array aa
        ! x, yy the matrix operations are applied to x and the result is yy
        ! a, b, c  are arrays which contain the tridiagonal matrix
        ! m  is the order of the matrix
        ! d, u, y are working arrays
        ! isgn  determines whether or not a change in sign is made
        !

        ! Dummy arguments

        integer(ip), intent(in)  :: nd
        integer(ip), intent(in)  :: nm1
        integer(ip), intent(in)  :: nm2
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: bm1(nm1)
        real(wp),    intent(in)  :: bm2(nm2)
        real(wp),    intent(in)  :: aa(na)
        real(wp),    intent(in)  :: x(m)
        real(wp),    intent(out) :: yy(m)
        real(wp),    intent(in)  :: a(m)
        real(wp),    intent(in)  :: b(m)
        real(wp),    intent(in)  :: c(m)
        complex(wp), intent(in)  :: bd(nd)
        complex(wp), intent(out) :: d(m)
        complex(wp), intent(out) :: u(m)
        complex(wp), intent(out) :: y(m)

        ! Local variables

        integer(ip) :: j, mm, mm2, id, m1, m2, ia, iflg, k
        real(wp)    :: rt
        complex(wp) :: v, den, bh, ym, am, y1, y2, yh, crt


        y = cmplx(x, ZERO, kind=wp)

        mm = m - 1
        mm2 = m - 2
        id = nd
        m1 = nm1
        m2 = nm2
        ia = na

        main_loop: do

            iflg = 0

            if (id > 0) then
                crt = bd(id)
                id = id - 1
                iflg = 1
                !
                ! begin solution to system
                !
                bh = b(m) - crt
                ym = y(m)
                den = b(1) - crt
                d(1) = c(1)/den
                u(1) = a(1)/den
                y(1) = y(1)/den
                v = cmplx(c(m), ZERO, kind=wp)

                if (mm2 >= 2) then
                    do j = 2, mm2
                        den = b(j) - crt - a(j)*d(j-1)
                        d(j) = c(j)/den
                        u(j) = -a(j)*u(j-1)/den
                        y(j) = (y(j)-a(j)*y(j-1))/den
                        bh = bh - v*u(j-1)
                        ym = ym - v*y(j-1)
                        v = -v*d(j-1)
                    end do
                end if

                den = b(m-1) - crt - a(m-1)*d(m-2)
                d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
                y(m-1) = (y(m-1)-a(m-1)*y(m-2))/den
                am = a(m) - v*d(m-2)
                bh = bh - v*u(m-2)
                ym = ym - v*y(m-2)
                den = bh - am*d(m-1)

                if (abs(den) /= ZERO) then
                    y(m) = (ym - am*y(m-1))/den
                else
                    y(m) = cmplx(ONE, ZERO, kind=wp)
                end if

                y(m-1) = y(m-1) - d(m-1)*y(m)

                do j = 2, mm
                    k = m - j
                    y(k) = y(k) - d(k)*y(k+1) - u(k)*y(m)
                end do
            end if

            if (m1 > 0 .and. m2 > 0) then
                if (m1 <= 0) then
                    rt = bm2(m2)
                    m2 = m2 - 1
                else
                    if (m2 <= 0) then
                        rt = bm1(m1)
                        m1 = m1 - 1
                    else
                        if (abs(bm1(m1)) - abs(bm2(m2)) > ZERO) then
                            rt = bm1(m1)
                            m1 = m1 - 1
                        else
                            rt = bm2(m2)
                            m2 = m2 - 1
                        !
                        ! matrix multiplication
                        !
                        end if
                    end if
                end if
                yh = y(1)
                y1 = (b(1)-rt)*y(1) + c(1)*y(2) + a(1)*y(m)
                if (mm >= 2) then
                    do j = 2, mm
                        y2 = a(j)*y(j-1) + (b(j)-rt)*y(j) + c(j)*y(j+1)
                        y(j-1) = y1
                        y1 = y2
                    end do
                end if
                y(m) = a(m)*y(m-1) + (b(m)-rt)*y(m) + c(m)*yh
                y(m-1) = y1
                iflg = 1
                cycle main_loop
            end if

            if (ia > 0) then
                rt = aa(ia)
                ia = ia - 1
                iflg = 1
                !
                ! scalar multiplication
                !
                y = rt*y
            end if

            if (iflg <= 0) exit main_loop

        end do main_loop

        yy = real(y, kind=wp)


    end subroutine cprodp

    pure subroutine compute_index_a_coeff(self, i, ir, idxa, na)

        ! Dummy arguments

        class(GeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(in)    :: i
        integer(ip),      intent(in)    :: ir
        integer(ip),      intent(out)   :: idxa
        integer(ip),      intent(out)   :: na


        associate( nm => self%nm )
            na = 2**ir
            idxa = i - na + 1
            if (i > nm) na = 0
        end associate

    end subroutine compute_index_a_coeff

    pure subroutine compute_index_b_coeff(self, i, ir, idx, idp)

        ! Dummy arguments

        class(GeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(in)    :: i
        integer(ip),      intent(in)    :: ir
        integer(ip),      intent(out)   :: idx
        integer(ip),      intent(out)   :: idp

        ! Local variables

        integer(ip) :: izh, id, ipl


        associate( &
            nm => self%nm, &
            ik => self%ik &
            )
            !
            ! b(idx) is the location of the first root of the b(i, ir) polynomial
            !
            idp = 0

            select case (ir)
                case (0)
                    if (i > nm) then
                        return
                    else
                        idx = i
                        idp = 1
                    end if
                case default
                    if (ir > 0) then
                        izh = 2**ir
                        id = i - 2*izh
                        idx = 2*id + (ir - 1)*ik + ir + (ik - i)/izh + 4
                        ipl = izh - 1
                        idp = 2*izh - 1
                        if (i - ipl - nm > 0) then
                            idp = 0
                            return
                        end if
                        if (i + ipl - nm > 0) idp = nm + ipl - i + 1
                    end if
            end select
        end associate

    end subroutine compute_index_b_coeff

    pure subroutine compute_index_c_coeff(self, i, ir, idxc, nc)

        ! Dummy arguments

        class(GeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(in)    :: i
        integer(ip),      intent(in)    :: ir
        integer(ip),      intent(out)   :: idxc
        integer(ip),      intent(out)   :: nc


        associate( nm => self%nm )
            nc = 2**ir
            idxc = i
            if (idxc + nc - 1 - nm > 0) nc = 0
        end associate

    end subroutine compute_index_c_coeff

    ! Purpose
    !
    ! Computes the eigenvalues of the periodic tridiagonal matrix
    ! with coefficients an, bn, cn
    !
    ! n is the order of the bh and bp polynomials
    ! on output bp contains the eigenvalues
    ! cbp is the same as bp except type complex
    ! bh is used to temporarily store the roots of the b hat polynomial
    ! which enters through bp
    !
    subroutine compute_eigen_values(self, ierror, a, c, cbp, bp, bh)

        ! Dummy arguments
        class(GeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(out)   :: ierror
        real(wp),         intent(in)    :: a(:)
        real(wp),         intent(in)    :: c(:)
        real(wp),         intent(inout) :: bp(:)
        real(wp),         intent(inout) :: bh(:)
        complex(wp),      intent(inout) :: cbp(:)

        ! Local variables
        integer(ip)   :: iz, izm, izm2, j, nt, modiz, is
        integer(ip)   :: iif, ig, it, icv, i3, i2, nhalf
        real(wp)      :: r4, r5, r6, scnv, xl, db, sgn, xr, xm, psg
        real(wp)      :: temp
        complex(wp)   :: cx, fsg, hsg, dd, f, fp, fpp, cdis, r1, r2, r3

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv &
            )

            associate( n => size(a) )

                scnv = sqrt(cnv)
                iz = n
                izm = iz - 1
                izm2 = iz - 2

                if (bp(n) <= bp(1)) then
                    if (bp(n) == bp(1)) then
                        ierror = 4
                        return
                    else
                        bh(:n) = bp(n:1:(-1))
                    end if
                else
                    bh(:n) = bp(:n)
                end if

                ncmplx = 0
                modiz = mod(iz, 2)
                is = 1

                if (modiz /= 0) then
                    if (a(1) >= ZERO) then
                        if (a(1) == ZERO) then
                            ierror = 4
                            return
                        end if
                    end if

                    xl = bh(1)
                    db = bh(3) - bh(1)
                    xl = xl - db
                    r4 = self%psgf(xl, iz, c, a, bh)

                    do
                        if (ZERO < r4) exit
                        xl = xl - db
                        r4 = self%psgf(xl, iz, c, a, bh)
                    end do

                    sgn = -ONE

                    temp = self%bsrh(xl, bh(1), iz, c, a, bh, sgn, psgf)

                    cbp(1) = cmplx(temp, ZERO, kind=wp)

                    bp(1) = real(cbp(1), kind=wp)

                    is = 2

                end if

                iif = iz - 1

                if (modiz /= 0) then

                    if (a(1) == ZERO) then
                        ierror = 4
                        return
                    end if

                    xr = bh(iz)
                    db = bh(iz) - bh(iz-2)
                    xr = xr + db
                    r5 = self%psgf(xr, iz, c, a, bh)

                    do
                        if (ZERO <= r5) exit
                        xr = xr + db
                        r5 = self%psgf(xr, iz, c, a, bh)
                    end do

                    sgn = ONE
                    temp = self%bsrh(bh(iz), xr, iz, c, a, bh, sgn, psgf)

                    cbp(iz) = cmplx(temp, ZERO, kind=wp)
                    iif = iz - 2

                end if

                main_loop: do ig = is, iif, 2

                    xl = bh(ig)
                    xr = bh(ig+1)
                    sgn = -ONE
                    xm = self%bsrh(xl, xr, iz, c, a, bh, sgn, ppspf)
                    psg = self%psgf(xm, iz, c, a, bh)

                    block_construct: block

                        if (abs(psg) > MACHINE_EPSILON) then

                            r6 = psg*self%psgf(xm, iz, c, a, bh)

                            if (r6 > ZERO) exit block_construct

                            if (r6 /= ZERO) then

                                sgn = ONE
                                temp = self%bsrh(bh(ig), xm, iz, c, a, bh, sgn, psgf)
                                cbp(ig) = cmplx(temp, ZERO, kind=wp)
                                sgn = -ONE
                                temp = self%bsrh(xm, bh(ig+1), iz, c, a, bh, sgn, psgf)
                                cbp(ig+1) = cmplx(temp, ZERO, kind=wp)

                                cycle main_loop
                            !
                            ! Case of a multiple zero
                            !
                            end if
                        end if

                        cbp(ig) = cmplx(xm, ZERO, kind=wp)
                        cbp(ig+1) = cmplx(xm, ZERO, kind=wp)

                        cycle main_loop
                    !
                    ! case of a complex zero
                    !
                    end block block_construct

                    it = 0
                    icv = 0
                    cx = cmplx(xm, ZERO, kind=wp)

                    loop_120: do

                        fsg = cmplx(ONE, ZERO, kind=wp)
                        hsg = cmplx(ONE, ZERO, kind=wp)
                        fp = ZERO
                        fpp = ZERO

                        do j = 1, iz
                            dd = ONE/(cx - bh(j))
                            fsg = fsg*a(j)*dd
                            hsg = hsg*c(j)*dd
                            fp = fp + dd
                            fpp = fpp - dd**2
                        end do

                        select case (modiz)
                            case (0)
                                f = cmplx(ONE, ZERO, kind=wp) - fsg - hsg
                            case default
                                f = cmplx(ONE, ZERO, kind=wp) + fsg + hsg
                        end select

                        i3 = 0

                        if (abs(fp) > ZERO) then
                            i3 = 1
                            r3 = -f/fp
                        end if

                        i2 = 0

                        if (abs(fpp) > ZERO) then

                            i2 = 1
                            cdis = sqrt((fp**2) - TWO * f * fpp)
                            r1 = cdis - fp
                            r2 = (-fp) - cdis

                            if (abs(r1) - abs(r2) > ZERO) then
                                r1 = r1/fpp
                            else
                                r1 = r2/fpp
                            end if

                            r2 = ((TWO*f)/fpp)/r1

                            if (abs(r2) < abs(r1)) r1 = r2

                            if (i3 > 0 .and. abs(r3) < abs(r1)) r1 = r3

                        else
                            r1 = r3
                        end if

                        cx = cx + r1
                        it = it + 1

                        if (it > 50) then
                            ierror = 4
                            return
                        end if

                        if (abs(r1) > scnv) cycle loop_120

                        if (icv <= 0) then
                            icv = 1
                            cycle loop_120
                        end if

                        exit loop_120
                    end do loop_120

                    cbp(ig) = cx
                    cbp(ig+1) = conjg(cx)

                end do main_loop

                if (abs(cbp(n)) - abs(cbp(1)) <= ZERO) then
                    if (abs(cbp(n)) - abs(cbp(1)) == ZERO) then
                        ierror = 4
                        return
                    end if

                    nhalf = n/2
                    do j = 1, nhalf
                        nt = n - j
                        cx = cbp(j)
                        cbp(j) = cbp(nt+1)
                        cbp(nt+1) = cx
                    end do
                end if

                ncmplx = 1

                if (any(aimag(cbp(2:iz)) /= ZERO)) return

                ncmplx = 0
                bp(1) = real(cbp(1), kind=wp)
                bp(2:iz) = real(cbp(2:iz))

            end associate

        end associate common_variables

    end subroutine compute_eigen_values

    ! Purpose:
    !
    ! Applies a sequence of matrix operations to the vector x and
    ! stores the result in y
    ! bd, bm1, bm2 are arrays containing roots of certian b polynomials
    ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
    ! aa   array containing scalar multipliers of the vector x
    ! na is the length of the array aa
    ! x, y  the matrix operations are applied to x and the result is y
    ! a, b, c  are arrays which contain the tridiagonal matrix
    ! m  is the order of the matrix
    ! d, w, u are working arrays
    ! is  determines whether or not a change in sign is made
    !
    pure subroutine prod(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, w, u)

        ! Dummy arguments
        integer(ip), intent(in)  :: nd
        integer(ip), intent(in)  :: nm1
        integer(ip), intent(in)  :: nm2
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: bd(nd)
        real(wp),    intent(in)  :: bm1(nm1)
        real(wp),    intent(in)  :: bm2(nm2)
        real(wp),    intent(in)  :: aa(na)
        real(wp),    intent(in)  :: x(m)
        real(wp),    intent(out) :: y(m)
        real(wp),    intent(in)  :: a(m)
        real(wp),    intent(in)  :: b(m)
        real(wp),    intent(in)  :: c(m)
        real(wp),    intent(out) :: d(m)
        real(wp),    intent(out) :: w(m)
        real(wp),    intent(out) :: u(m)

        ! Local variables
        integer(ip) :: j, mm, id, ibr, m1, m2, ia, k
        real(wp)    :: rt, den

        w = x
        y = w
        mm = m - 1
        id = nd
        ibr = 0
        m1 = nm1
        m2 = nm2
        ia = na

        main_loop: do

            if (ia > 0) then

                if (nd == 0) then
                    rt = -aa(ia)
                else
                    rt = aa(ia)
                end if

                ia = ia - 1
                !
                ! scalar multiplication
                !
                y = rt*w
            end if

            if (id <= 0) return

            rt = bd(id)
            id = id - 1

            if (id == 0) ibr = 1

            !
            ! begin solution to system
            !
            d(m) = a(m)/(b(m)-rt)
            w(m) = y(m)/(b(m)-rt)

            do j = 2, mm
                k = m - j
                den = b(k+1) - rt - c(k+1)*d(k+2)
                d(k+1) = a(k+1)/den
                w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
            end do

            den = b(1) - rt - c(1)*d(2)
            w(1) = ONE

            if (den /= ZERO) w(1) = (y(1)-c(1)*w(2))/den

            do j = 2, m
                w(j) = w(j) - d(j)*w(j-1)
            end do

            if (na > 0) cycle main_loop

            if (m1 <= 0) then
                if (m2 <= 0) then
                    y = w
                    ibr = 1
                    cycle main_loop
                end if
            else
                if (.not.(m2 > 0 .and. abs(bm1(m1)) <= abs(bm2(m2)))) then
                    if (ibr <= 0) then
                        if (abs(bm1(m1)-bd(id)) < abs(bm1(m1)-rt)) then
                            y(:m) = w(:m)
                            ibr = 1
                            cycle main_loop
                        end if
                    end if
                    rt = rt - bm1(m1)
                    m1 = m1 - 1
                    y = y + rt*w
                    cycle main_loop
                end if
            end if

            if (ibr <= 0 .and. abs(bm2(m2)-bd(id)) < abs(bm2(m2)-rt)) then
                y = w
                ibr = 1
                cycle main_loop
            end if

            rt = rt - bm2(m2)
            m2 = m2 - 1
            y = y + rt*w

        end do main_loop

    end subroutine prod

    ! Purpose:
    !
    ! Applies a sequence of matrix operations to the vector x and
    ! stores the result in y periodic boundary conditions
    !
    ! bd, bm1, bm2 are arrays containing roots of certain b polynomials
    !
    ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
    !
    ! aa   array containing scalar multipliers of the vector x
    !
    ! na is the length of the array aa
    !
    ! x, y  the matrix operations are applied to x and the result is y
    !
    ! a, b, c  are arrays which contain the tridiagonal matrix
    !
    ! m  is the order of the matrix
    !
    ! d, u, w are working arrays
    !
    ! is  determines whether or not a change in sign is made
    !
    pure subroutine prodp(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, u, w)

        ! Dummy arguments
        integer(ip), intent(in)  :: nd
        integer(ip), intent(in)  :: nm1
        integer(ip), intent(in)  :: nm2
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: bd(nd)
        real(wp),    intent(in)  :: bm1(nm1)
        real(wp),    intent(in)  :: bm2(nm2)
        real(wp),    intent(in)  :: aa(na)
        real(wp),    intent(in)  :: x(m)
        real(wp),    intent(out) :: y(m)
        real(wp),    intent(in)  :: a(m)
        real(wp),    intent(in)  :: b(m)
        real(wp),    intent(in)  :: c(m)
        real(wp),    intent(out) :: d(m)
        real(wp),    intent(out) :: u(m)
        real(wp),    intent(out) :: w(m)

        ! Local variables
        integer(ip) :: j, mm, mm2, id, ibr, m1, m2, ia, k
        real(wp)    :: rt, bh, ym, den, v, am

        y = x
        w = y
        mm = m - 1
        mm2 = m - 2
        id = nd
        ibr = 0
        m1 = nm1
        m2 = nm2
        ia = na

        main_loop: do

            if (ia > 0) then
                rt = aa(ia)

                if (nd == 0) rt = -rt

                ia = ia - 1
                y = rt*w
            end if

            if (id <= 0) return

            rt = bd(id)
            id = id - 1

            if (id == 0) ibr = 1

            !
            ! begin solution to system
            !
            bh = b(m) - rt
            ym = y(m)
            den = b(1) - rt
            d(1) = c(1)/den
            u(1) = a(1)/den
            w(1) = y(1)/den
            v = c(m)

            if (mm2 >= 2) then
                do j = 2, mm2
                    den = b(j) - rt - a(j)*d(j-1)
                    d(j) = c(j)/den
                    u(j) = -a(j)*u(j-1)/den
                    w(j) = (y(j)-a(j)*w(j-1))/den
                    bh = bh - v*u(j-1)
                    ym = ym - v*w(j-1)
                    v = -v*d(j-1)
                end do
            end if

            den = b(m-1) - rt - a(m-1)*d(m-2)
            d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
            w(m-1) = (y(m-1)-a(m-1)*w(m-2))/den
            am = a(m) - v*d(m-2)
            bh = bh - v*u(m-2)
            ym = ym - v*w(m-2)
            den = bh - am*d(m-1)

            if (den /= ZERO) then
                w(m) = (ym - am*w(m-1))/den
            else
                w(m) = ONE
            end if

            w(m-1) = w(m-1) - d(m-1)*w(m)

            do j = 2, mm
                k = m - j
                w(k) = w(k) - d(k)*w(k+1) - u(k)*w(m)
            end do

            if (na <= 0) then
                if (m1 <= 0 .and. m2 <= 0) then
                    y = w
                    ibr = 1
                    cycle main_loop
                else

                    if (m2 > 0 .and. abs(bm1(m1)) - abs(bm2(m2)) <= ZERO) then
                        if (ibr <= 0 .and. abs(bm2(m2)-bd(id)) - abs(bm2(m2)-rt) < ZERO) then
                            y = w
                            ibr = 1
                            cycle main_loop
                        else
                            rt = rt - bm2(m2)
                            m2 = m2 - 1
                            y = y + rt*w
                        end if
                    end if

                    if (ibr <= 0 .and. abs(bm1(m1)-bd(id)) - abs(bm1(m1)-rt) < ZERO) then
                        y = w
                        ibr = 1
                        cycle main_loop
                    else
                        rt = rt - bm1(m1)
                        m1 = m1 - 1
                        y = y + rt*w
                        cycle main_loop
                    end if
                end if

                if (ibr <= 0 .and. abs(bm2(m2)-bd(id)) - abs(bm2(m2)-rt) < ZERO) then
                    y = w
                    ibr = 1
                else
                    rt = rt - bm2(m2)
                    m2 = m2 - 1
                    y = y + rt*w
                end if
            end if
        end do main_loop

    end subroutine prodp

    ! Purpose:
    !
    !     This subroutine is a modification of the eispack subroutine tqlrat
    !     algorithm 464, comm. acm 16, 689(1973) by reinsch.
    !
    !     this subroutine finds the eigenvalues of a symmetric
    !     tridiagonal matrix by the rational ql method.
    !
    !     on input-
    !
    !        d contains the diagonal elements of the input matrix,
    !
    !        e2 contains the subdiagonal elements of the
    !          input matrix in its last n-1 positions.  e2(1) is arbitrary.
    !
    !      on output-
    !
    !        d contains the eigenvalues in ascending order.  if an
    !          error exit is made, the eigenvalues are correct and
    !          ordered for indices 1, 2, ...error_flag-1, but may not be
    !          the smallest eigenvalues,
    !
    !        e2 has been destroyed,
    !
    !        error_flag is set to
    !          zero       for normal return,
    !          j          if the j-th eigenvalue has not been
    !                     determined after 30 iterations.
    !
    !     questions and comments should be directed to b. s. garbow,
    !     applied mathematics division, argonne national laboratory
    !
    !
    !     eps is a machine dependent parameter specifying
    !     the relative precision of floating point arithmetic.
    !
    subroutine tevls(self, diagonal, subdiagonal, error_flag)

        ! Dummy arguments
        class(GeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(out)   :: error_flag
        real(wp),         intent(inout) :: diagonal(:)
        real(wp),         intent(inout) :: subdiagonal(:)

        ! Local variables
        integer(ip) :: i, j, l, m, ii, l1, mml, nhalf, ntop
        real(wp)    :: b, c, f, g, h, p, r, s, dhold

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv &
            )

            error_flag = 0

            associate( &
                n => size(diagonal), &
                d => diagonal, &
                e2 => subdiagonal &
                )

                if (n == 1) return

                e2(:n-1) = e2(2:n)**2
                f = ZERO
                b = ZERO
                e2(n) = ZERO

                main_loop: do l = 1, n
                    j = 0
                    h = MACHINE_EPSILON*(abs(d(l))+sqrt(e2(l)))

                    if (b <= h) then
                        b = h
                        c = b**2
                    end if
                    !
                    !  look for small squared sub-diagonal element
                    !
                    do m = l, n
                        if (e2(m) > c) then
                            cycle
                        end if
                        exit
                    !
                    ! e2(n) is always zero, so there is no exit
                    !    through the bottom of the loop
                    !
                    end do

                    if_construct: if (m /= l) then

                        loop_105: do

                            if (j == 30) then
                                !
                                ! set error: no convergence to an
                                !    eigenvalue after 30 iterations
                                !
                                error_flag = l
                                return
                            end if

                            j = j + 1
                            !
                            ! form shift
                            !
                            l1 = l + 1
                            s = sqrt(e2(l))
                            g = d(l)
                            p = (d(l1)-g)/(2.0*s)
                            r = sqrt(p**2 + ONE)
                            d(l) = s/(p + sign(r, p))
                            h = g - d(l)
                            d(l1:n) = d(l1:n) - h
                            f = f + h
                            !
                            ! rational ql transformation
                            !
                            if (g == ZERO) then
                                g = b
                            else
                                g = d(m)
                            end if

                            h = g
                            s = ZERO
                            mml = m - l
                            !
                            ! for i = m-1 step -1 until l do
                            !
                            do ii = 1, mml
                                i = m - ii
                                p = g*h
                                r = p + e2(i)
                                e2(i+1) = s*r
                                s = e2(i)/r
                                d(i+1) = h + s*(h + d(i))
                                g = d(i) - e2(i)/g

                                if (g == ZERO) g = b

                                h = g*p/r
                            end do
                            !
                            e2(l) = s*g
                            d(l) = h
                            !
                            ! guard against underflowed h
                            !
                            if (h == ZERO) exit if_construct

                            if (abs(e2(l)) <= abs(c/h)) exit if_construct

                            e2(l) = h*e2(l)

                            if (e2(l) == ZERO) exit loop_105

                        end do loop_105
                    end if if_construct

                    p = d(l) + f
                    !
                    ! order eigenvalues
                    !
                    if (l /= 1) then
                        !
                        ! for i=l step -1 until 2 do
                        !
                        do ii = 2, l
                            i = l + 2 - ii
                            if (p >= d(i-1)) then
                                d(i) = p
                                cycle main_loop
                            else
                                d(i) = d(i-1)
                            end if
                        end do
                    end if

                    i = 1
                    d(i) = p

                end do main_loop

                if (abs(d(n)) < abs(d(1))) then
                    nhalf = n/2
                    do i = 1, nhalf
                        ntop = n - i
                        dhold = d(i)
                        d(i) = d(ntop+1)
                        d(ntop+1) = dhold
                    end do
                end if

            end associate
        end associate common_variables

    end subroutine tevls

end module type_GeneralizedCyclicReductionUtility

