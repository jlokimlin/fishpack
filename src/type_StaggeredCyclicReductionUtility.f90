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
module type_StaggeredCyclicReductionUtility

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_CyclicReductionUtility, only: &
        CyclicReductionUtility

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: CyclicReductionUtility

    ! Parameters confined to the module
    real(wp),    parameter :: ZERO = 0.0_wp
    real(wp),    parameter :: HALF = 0.5_wp
    real(wp),    parameter :: ONE = 1.0_wp
    real(wp),    parameter :: TWO = 2.0_wp
    integer(ip), parameter :: SIZE_OF_WORKSPACE_INDICES = 11

    type, public, extends(CyclicReductionUtility) :: StaggeredCyclicReductionUtility
    contains
        ! Public type-bound procedures
        procedure, public :: poistg_lower_routine
        ! Private type-bound procedures
        procedure, private :: poistg_solve_poisson_on_staggered_grid
    end type

contains

    subroutine poistg_lower_routine(self, nperod, n, mperod, m, a, b, c, idimy, y, ierror, w)

        ! Local variables
        class(StaggeredCyclicReductionUtility), intent(inout) :: self
        integer(ip), intent(in)     :: nperod
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: mperod
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: idimy
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: a(m)
        real(wp),    intent(in)     :: b(m)
        real(wp),    intent(in)     :: c(m)
        real(wp),    intent(inout)  :: y(idimy,n)
        real(wp),    intent(inout)  :: w(:)

        ! Local variables
        integer(ip)     :: workspace_indices(SIZE_OF_WORKSPACE_INDICES)
        integer(ip)     :: i, k, j, np, mp
        integer(ip)     :: nby2, mskip, ipstor, irev, mh, mhm1, m_odd
        real(wp)        :: temp

        ! Check validity of calling arguments
        call poistg_check_input_arguments(nperod, n, mperod, m, idimy, ierror, a, b, c)

        ! Check error flag
        if (ierror /= 0) return

        ! Compute workspace indices
        workspace_indices = poistg_get_workspace_indices(n,m)

        associate( &
            iwba => workspace_indices(1), &
            iwbb => workspace_indices(2), &
            iwbc => workspace_indices(3), &
            iwb2 => workspace_indices(4), &
            iwb3 => workspace_indices(5), &
            iww1 => workspace_indices(6), &
            iww2 => workspace_indices(7), &
            iww3 => workspace_indices(8), &
            iwd => workspace_indices(9), &
            iwtcos => workspace_indices(10), &
            iwp => workspace_indices(11) &
            )

            do i = 1, m
                k = iwba + i - 1
                w(k) = -a(i)
                k = iwbc + i - 1
                w(k) = -c(i)
                k = iwbb + i - 1
                w(k) = TWO - b(i)
                y(i, :n) = -y(i, :n)
            end do

            np = nperod
            mp = mperod + 1

            if (mp == 1) then
                mh = (m + 1)/2
                mhm1 = mh - 1

                if (mh*2 == m) then
                    m_odd = 2
                else
                    m_odd = 1
                end if

                do j = 1, n
                    do i = 1, mhm1
                        w(i) = y(mh-i, j) - y(i+mh, j)
                        w(i+mh) = y(mh-i, j) + y(i+mh, j)
                    end do
                    w(mh) = TWO * y(mh, j)
                    select case (m_odd)
                        case (1)
                            y(:m, j) = w(:m)
                        case (2)
                            w(m) = TWO * y(m, j)
                            y(:m, j) = w(:m)
                    end select
                end do

                k = iwbc + mhm1 - 1
                i = iwba + mhm1
                w(k) = ZERO
                w(i) = ZERO
                w(k+1) = TWO * w(k+1)

                select case (m_odd)
                    case (2)
                        w(iwbb-1) = w(k+1)
                    case default
                        k = iwbb + mhm1 - 1
                        w(k) = w(k) - w(i-1)
                        w(iwbc-1) = w(iwbc-1) + w(iwbb-1)
                end select
            end if

            loop_107: do
                loop_108: do
                    if (nperod /= 4) then

                        call self%poistg_solve_poisson_on_staggered_grid(np, n, m, w(iwba:), w(iwbb:), w(iwbc:), idimy, y, &
                            w, w(iwb2:), w(iwb3:), w(iww1:), w(iww2:), w(iww3:), &
                            w(iwd:), w(iwtcos:), w(iwp:))

                        ipstor = int(w(iww1), kind=ip)
                        irev = 2
                        select case (mp)
                            case (1)
                                exit loop_107
                            case (2)
                                w(1) = real(ipstor + iwp - 1, kind=wp)
                                return
                        end select
                        cycle loop_107
                    end if

                    irev = 1
                    nby2 = n/2
                    np = 2

                    do j = 1, nby2
                        mskip = n + 1 - j
                        do i = 1, m
                            temp = y(i, j)
                            y(i, j) = y(i, mskip)
                            y(i, mskip) = temp
                        end do
                    end do

                    select case (irev)
                        case (1)
                            cycle loop_108
                        case (2)
                            select case (mp)
                                case (1)
                                    exit loop_107
                                case (2)
                                    w(1) = real(ipstor + iwp - 1, kind=wp)
                                    return
                            end select
                            cycle loop_107
                    end select

                    exit loop_108
                end do loop_108
                exit loop_107
            end do loop_107

            do j = 1, n
                w(mh-1:mh-mhm1:(-1)) = HALF * (y(mh+1:mhm1+mh, j)+y(:mhm1, j))
                w(mh+1:mhm1+mh) = HALF * (y(mh+1:mhm1+mh, j)-y(:mhm1, j))
                w(mh) = HALF * y(mh, j)
                select case (m_odd)
                    case (1)
                        y(:m, j) = w(:m)
                    case (2)
                        w(m) = HALF * y(m, j)
                        y(:m, j) = w(:m)
                end select
            end do

            w(1) = real(ipstor + iwp - 1, kind=wp)

        end associate

    end subroutine poistg_lower_routine

    pure subroutine poistg_check_input_arguments(nperod, n, mperod, m, idimy, ierror, a, b, c)

        ! Dummy arguments
        integer(ip), intent(in)  :: nperod
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: mperod
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: idimy
        integer(ip), intent(out) :: ierror
        real(wp),    intent(in)  :: a(m)
        real(wp),    intent(in)  :: b(m)
        real(wp),    intent(in)  :: c(m)

        if (3 > m) then
            ierror = 1
        else if (3 > n) then
            ierror = 2
        else if (idimy < m) then
            ierror = 3
        else if (nperod < 1 .or. nperod > 4) then
            ierror = 4
        else if (mperod < 0 .or. mperod > 1) then
            ierror = 5
        else if (mperod /= 1) then
            if (any(a /= c(1)) .or. any(c /= c(1)) .or. any(b /= b(1))) then
                ierror = 6
            end if
        else if (a(1) /= ZERO .or. c(m) /= ZERO) then
            ierror = 7
        else
            ierror = 0
        end if

    end subroutine poistg_check_input_arguments

    function poistg_get_workspace_indices(n, m) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: n
        integer(ip), intent(in) :: m
        integer(ip)             :: return_value(SIZE_OF_WORKSPACE_INDICES)

        integer(ip) :: j !! Counter

        associate( i => return_value)
            i(1) = m + 1
            do j = 2, SIZE_OF_WORKSPACE_INDICES - 1
                i(j) = i(j-1) + m
            end do
            i(SIZE_OF_WORKSPACE_INDICES) = i(SIZE_OF_WORKSPACE_INDICES - 1) + 4 * n
        end associate

    end function poistg_get_workspace_indices

    subroutine poistg_solve_poisson_on_staggered_grid(self, &
        nperod, n, m, a, bb, c, idimq, q, b, b2, b3, w, w2, w3, d, tcos, p)

        ! Dummy arguments
        class(StaggeredCyclicReductionUtility), intent(inout) :: self
        integer(ip), intent(in)    :: nperod
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: idimq
        real(wp),    intent(in)    :: a(m)
        real(wp),    intent(in)    :: bb(m)
        real(wp),    intent(in)    :: c(m)
        real(wp),    intent(inout) :: q(idimq, n)
        real(wp),    intent(inout) :: b(m)
        real(wp),    intent(inout) :: b2(m)
        real(wp),    intent(inout) :: b3(m)
        real(wp),    intent(inout) :: w(m)
        real(wp),    intent(inout) :: w2(m)
        real(wp),    intent(inout) :: w3(m)
        real(wp),    intent(inout) :: d(m)
        real(wp),    intent(inout) :: tcos(m)
        real(wp),    intent(inout) :: p(4*n)

        ! Local variables
        integer(ip)     :: k(4)
        integer(ip)     :: np, mr, ipp, ipstor, i2r, jr, nr, nlast
        integer(ip)     :: kr, lr, nrod, jstart, jstop, i2rby2
        integer(ip)     :: j, ijump, jp1, jp2, jp3
        integer(ip)     :: jm1, jm2, jm3, i, nrodpr, ii, nlastp, jstep
        real(wp)        :: fnum, fnum2, fi, t

        associate( &
            k1 => k(1), &
            k2 => k(2), &
            k3 => k(3), &
            k4 => k(4) &
            )

            np = nperod
            fnum = HALF * real(np/3, kind=wp)
            fnum2 = HALF * real(np/2, kind=wp)
            mr = m
            ipp = -mr
            ipstor = 0
            i2r = 1
            jr = 2
            nr = n
            nlast = n
            kr = 1
            lr = 0

            if (nr > 3) then

                loop_101: do

                    jr = 2*i2r

                    if ((nr/2)*2 == nr) then
                        nrod = 0
                    else
                        nrod = 1
                    end if

                    jstart = 1
                    jstop = nlast - jr

                    if (nrod == 0) jstop = jstop - i2r

                    i2rby2 = i2r/2

                    if (jstop < jstart) then
                        j = jr
                    else
                        ijump = 1
                        do j = jstart, jstop, jr
                            jp1 = j + i2rby2
                            jp2 = j + i2r
                            jp3 = jp2 + i2rby2
                            jm1 = j - i2rby2
                            jm2 = j - i2r
                            jm3 = jm2 - i2rby2

                            if (j == 1) then
                                call self%generate_cosines(i2r, 1, fnum, HALF, tcos)
                                select case (i2r)
                                    case (1)
                                        b(:mr) = q(:mr, 1)
                                        q(:mr, 1) = q(:mr, 2)
                                    case default
                                        b(:mr) = q(:mr, 1) + HALF * (q(:mr, jp2)-q(:mr, jp1)-q(:mr, &
                                            jp3))
                                        q(:mr, 1) = q(:mr, jp2) + q(:mr, 1) - q(:mr, jp1)
                                end select
                            else
                                if (ijump == 1) then
                                    ijump = 2
                                    call self%generate_cosines(i2r, 1, HALF, ZERO, tcos)
                                end if

                                select case (i2r)
                                    case (1)
                                        b(:mr) = TWO*q(:mr, j)
                                        q(:mr, j) = q(:mr, jm2) + q(:mr, jp2)
                                    case default
                                        do i = 1, mr
                                            fi = q(i, j)
                                            q(i, j)=q(i, j)-q(i, jm1)-q(i, jp1)+q(i, jm2)+q(i, jp2)
                                            b(i) = fi + q(i, j) - q(i, jm3) - q(i, jp3)
                                        end do
                                end select
                            end if

                            call self%solve_tridiag(i2r, 0, mr, a, bb, c, b, tcos, d, w)

                            q(:mr, j) = q(:mr, j) + b(:mr)
                            !
                            ! End of reduction for regular unknowns.
                            !
                        end do
                        !
                        ! Begin special reduction for last unknown.
                        !
                        j = jstop + jr
                    end if

                    nlast = j
                    jm1 = j - i2rby2
                    jm2 = j - i2r
                    jm3 = jm2 - i2rby2

                    if (nrod /= 0) then
                        !
                        ! Odd number of unknowns
                        !
                        select case (i2r)
                            case (1)
                                b(:mr) = q(:mr, j)
                                q(:mr, j) = q(:mr, jm2)
                            case default
                                b(:mr)=q(:mr, j)+HALF * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))

                                select case (nrodpr)
                                    case (0)
                                        q(:mr, j) = q(:mr, jm2) + p(ipp+1:mr+ipp)
                                        ipp = ipp - mr
                                    case default
                                        q(:mr, j) = q(:mr, j) - q(:mr, jm1) + q(:mr, jm2)
                                end select

                                if (lr /= 0) call self%generate_cosines(lr, 1, fnum2, HALF, tcos(kr+1:))

                        end select

                        call self%generate_cosines(kr, 1, fnum2, HALF, tcos)
                        call self%solve_tridiag(kr, lr, mr, a, bb, c, b, tcos, d, w)

                        q(:mr, j) = q(:mr, j) + b(:mr)
                        kr = kr + i2r
                    else
                        jp1 = j + i2rby2
                        jp2 = j + i2r

                        select case (i2r)
                            case (1)
                                b(:mr) = q(:mr, j)
                                tcos(1) = ZERO

                                call self%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)

                                ipp = 0
                                ipstor = mr
                                p(:mr) = b(:mr)
                                b(:mr) = b(:mr) + q(:mr, n)
                                tcos(1) = -ONE + TWO*real(np/2, kind=wp)
                                tcos(2) = ZERO

                                call self%solve_tridiag(1, 1, mr, a, bb, c, b, tcos, d, w)

                                q(:mr, j) = q(:mr, jm2) + p(:mr) + b(:mr)

                            case default

                                b(:mr) = q(:mr, j) + HALF * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))

                                select case (nrodpr)
                                    case (0)
                                        b(:mr) = b(:mr) + p(ipp+1:mr+ipp)
                                    case default
                                        b(:mr) = b(:mr) + q(:mr, jp2) - q(:mr, jp1)
                                end select

                                call self%generate_cosines(i2r, 1, HALF, ZERO, tcos)
                                call self%solve_tridiag(i2r, 0, mr, a, bb, c, b, tcos, d, w)

                                ipp = ipp + mr
                                ipstor = max(ipstor, ipp + mr)
                                p(ipp+1:mr+ipp) = b(:mr) + HALF * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
                                b(:mr) = p(ipp+1:mr+ipp) + q(:mr, jp2)

                                if (lr /= 0) then
                                    call self%generate_cosines(lr, 1, fnum2, HALF, tcos(i2r+1:))
                                    call self%merge_cosines(tcos, 0, i2r, i2r, lr, kr)
                                else
                                    do i = 1, i2r
                                        ii = kr + i
                                        tcos(ii) = tcos(i)
                                    end do
                                end if

                                call self%generate_cosines(kr, 1, fnum2, HALF, tcos)
                                call self%solve_tridiag(kr, kr, mr, a, bb, c, b, tcos, d, w)

                                q(:mr, j) = q(:mr, jm2) + p(ipp+1:mr+ipp) + b(:mr)

                        end select
                        lr = kr
                        kr = kr + jr
                    end if

                    nr = (nlast - 1)/jr + 1

                    if (nr <= 3) exit loop_101

                    i2r = jr
                    nrodpr = nrod

                end do loop_101
            end if


            j = 1 + jr
            jm1 = j - i2r
            jp1 = j + i2r
            jm2 = nlast - i2r

            block_construct: block
                if_nr:  if (nr /= 2) then
                    if_lr:  if (lr == 0) then
                        if_n: if (n == 3) then
                            !
                            ! case n = 3.
                            !
                            select case (np)
                                case (1,3)
                                    b(:mr) = q(:mr, 2)
                                    b2(:mr) = q(:mr, 1) + q(:mr, 3)
                                    b3(:mr) = ZERO

                                    select case (np)
                                        case (1:2)
                                            tcos(1) = -TWO
                                            tcos(2) = ONE
                                            tcos(3) = -ONE
                                            k1 = 2
                                        case default
                                            tcos(1) = -ONE
                                            tcos(2) = ONE
                                            k1 = 1
                                    end select

                                    k2 = 1
                                    k3 = 0
                                    k4 = 0
                                case (2)
                                    b(:mr) = q(:mr, 2)
                                    b2(:mr) = q(:mr, 3)
                                    b3(:mr) = q(:mr, 1)

                                    call self%generate_cosines(3, 1, HALF, ZERO, tcos)

                                    tcos(4) = -ONE
                                    tcos(5) = ONE
                                    tcos(6) = -ONE
                                    tcos(7) = ONE

                                    k1 = 3
                                    k2 = 2
                                    k3 = 1
                                    k4 = 1
                            end select


                            call self%solve_tridiag3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)

                            b(:mr) = b(:mr) + b2(:mr) + b3(:mr)

                            if (np == 3) then
                                tcos(1) = TWO
                                call self%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)
                            end if

                            q(:mr, 2) = b(:mr)
                            b(:mr) = q(:mr, 1) + b(:mr)
                            tcos(1) = -ONE + 4.0_wp*fnum

                            call self%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)

                            q(:mr, 1) = b(:mr)
                            jr = 1
                            i2r = 0
                            exit block_construct
                        end if if_n
                        !
                        !     case n = 2**p+1
                        !
                        b(:mr)=q(:mr, j)+q(:mr, 1)-q(:mr, jm1)+q(:mr, nlast)-q(:mr, jm2)

                        select case (np)
                            case (1, 3)
                                b2(:mr) = q(:mr, 1) + q(:mr, nlast) &
                                    + q(:mr, j) - q(:mr, jm1) - q(:mr, jp1)
                                b3(:mr) = ZERO
                                k1 = nlast - 1
                                k2 = nlast + jr - 1

                                call self%generate_cosines(jr - 1, 1, ZERO, ONE, tcos(nlast:))

                                tcos(k2) = TWO*real(np - 2, kind=wp)

                                call self%generate_cosines(jr, 1, HALF - fnum, HALF, tcos(k2+1:))

                                k3 = (3 - np)/2

                                call self%merge_cosines(tcos, k1, jr - k3, k2 - k3, jr + k3, 0)

                                k1 = k1 - 1 + k3

                                call self%generate_cosines(jr, 1, fnum, HALF, tcos(k1+1:))

                                k2 = jr
                                k3 = 0
                                k4 = 0
                            case (2)
                                do i = 1, mr
                                    fi = (q(i, j)-q(i, jm1)-q(i, jp1))/2
                                    b2(i) = q(i, 1) + fi
                                    b3(i) = q(i, nlast) + fi
                                end do

                                k1 = nlast + jr - 1
                                k2 = k1 + jr - 1

                                call self%generate_cosines(jr - 1, 1, ZERO, ONE, tcos(k1+1:))
                                call self%generate_cosines(nlast, 1, HALF, ZERO, tcos(k2+1:))
                                call self%merge_cosines(tcos, k1, jr - 1, k2, nlast, 0)

                                k3 = k1 + nlast - 1
                                k4 = k3 + jr

                                call self%generate_cosines(jr, 1, HALF, HALF, tcos(k3+1:))
                                call self%generate_cosines(jr, 1, ZERO, HALF, tcos(k4+1:))
                                call self%merge_cosines(tcos, k3, jr, k4, jr, k1)

                                k2 = nlast - 1
                                k3 = jr
                                k4 = jr
                        end select

                        call self%solve_tridiag3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
                        b(:mr) = b(:mr) + b2(:mr) + b3(:mr)

                        if (np == 3) then
                            tcos(1) = TWO
                            call self%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)
                        end if

                        q(:mr, j) = b(:mr) + HALF * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
                        b(:mr) = q(:mr, j) + q(:mr, 1)

                        call self%generate_cosines(jr, 1, fnum, HALF, tcos)
                        call self%solve_tridiag(jr, 0, mr, a, bb, c, b, tcos, d, w)

                        q(:mr, 1) = q(:mr, 1) - q(:mr, jm1) + b(:mr)

                        exit block_construct
                    end if if_lr
                    !
                    ! case of general n with nr = 3 .
                    !
                    b(:mr) = q(:mr, 1) - q(:mr, jm1) + q(:mr, j)

                    if (nrod == 0) then
                        b(:mr) = b(:mr) + p(ipp+1:mr+ipp)
                    else
                        b(:mr) = b(:mr) + q(:mr, nlast) - q(:mr, jm2)
                    end if

                    do i = 1, mr
                        t = HALF * (q(i, j)-q(i, jm1)-q(i, jp1))
                        q(i, j) = t
                        b2(i) = q(i, nlast) + t
                        b3(i) = q(i, 1) + t
                    end do

                    k1 = kr + 2*jr

                    call self%generate_cosines(jr - 1, 1, ZERO, ONE, tcos(k1+1:))

                    k2 = k1 + jr
                    tcos(k2) = TWO*real(np - 2, kind=wp)
                    k4 = (np - 1)*(3 - np)
                    k3 = k2 + 1 - k4

                    !generate_cosines(n, ijump, fnum, fden, a)
                    associate( &
                        n_arg => kr+jr+k4, &
                        fnum_arg => real(k4, kind=wp)/2, &
                        fden_arg => ONE-real(k4, kind=wp) &
                        )
                        call self%generate_cosines(n_arg, 1, fnum_arg , fden_arg, tcos(k3:))
                    end associate

                    k4 = 1 - np/3

                    call self%merge_cosines(tcos, k1, jr - k4, k2 - k4, kr + jr + k4, 0)

                    if (np == 3) k1 = k1 - 1

                    k2 = kr + jr
                    k4 = k1 + k2

                    call self%generate_cosines(kr, 1, fnum2, HALF, tcos(k4+1:))

                    k3 = k4 + kr

                    call self%generate_cosines(jr, 1, fnum, HALF, tcos(k3+1:))
                    call self%merge_cosines(tcos, k4, kr, k3, jr, k1)

                    k4 = k3 + jr

                    call self%generate_cosines(lr, 1, fnum2, HALF, tcos(k4+1:))
                    call self%merge_cosines(tcos, k3, jr, k4, lr, k1 + k2)
                    call self%generate_cosines(kr, 1, fnum2, HALF, tcos(k3+1:))

                    k3 = kr
                    k4 = kr

                    call self%solve_tridiag3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)

                    b(:mr) = b(:mr) + b2(:mr) + b3(:mr)

                    if (np == 3) then
                        tcos(1) = TWO
                        call self%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)
                    end if

                    q(:mr, j) = q(:mr, j) + b(:mr)
                    b(:mr) = q(:mr, 1) + q(:mr, j)

                    call self%generate_cosines(jr, 1, fnum, HALF, tcos)
                    call self%solve_tridiag(jr, 0, mr, a, bb, c, b, tcos, d, w)

                    if (jr == 1) then
                        q(:mr, 1) = b(:mr)
                        exit block_construct
                    end if

                    q(:mr, 1) = q(:mr, 1) - q(:mr, jm1) + b(:mr)

                    exit block_construct
                end if if_nr

                b3(:mr) = ZERO
                b(:mr) = q(:mr, 1) + p(ipp+1:mr+ipp)
                q(:mr, 1) = q(:mr, 1) - q(:mr, jm1)
                b2(:mr) = q(:mr, 1) + q(:mr, nlast)
                k1 = kr + jr
                k2 = k1 + jr

                call self%generate_cosines(jr - 1, 1, ZERO, ONE, tcos(k1+1:))

                select case (np)
                    case (1, 3)
                        tcos(k2) = TWO*real(np - 2, kind=wp)
                        call self%generate_cosines(kr, 1, ZERO, ONE, tcos(k2+1:))
                    case (2)
                        call self%generate_cosines(kr + 1, 1, HALF, ZERO, tcos(k2:))
                end select

                k4 = 1 - np/3

                call self%merge_cosines(tcos, k1, jr - k4, k2 - k4, kr + k4, 0)

                if (np == 3) k1 = k1 - 1

                k2 = kr

                call self%generate_cosines(kr, 1, fnum2, HALF, tcos(k1+1:))

                k4 = k1 + kr

                call self%generate_cosines(lr, 1, fnum2, HALF, tcos(k4+1:))

                k3 = lr
                k4 = 0

                call self%solve_tridiag3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)

                b(:mr) = b(:mr) + b2(:mr)

                if (np == 3) then
                    tcos(1) = TWO
                    call self%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)
                end if

                q(:mr, 1) = q(:mr, 1) + b(:mr)

            end block block_construct

            loop_188: do

                j = nlast - jr
                b(:mr) = q(:mr, nlast) + q(:mr, j)
                jm2 = nlast - i2r

                if (jr == 1) then
                    q(:mr, nlast) = ZERO
                else
                    select case (nrod)
                        case (0)
                            q(:mr, nlast) = p(ipp+1:mr+ipp)
                            ipp = ipp - mr
                        case default
                            q(:mr, nlast) = q(:mr, nlast) - q(:mr, jm2)
                    end select
                end if

                call self%generate_cosines(kr, 1, fnum2, HALF, tcos)
                call self%generate_cosines(lr, 1, fnum2, HALF, tcos(kr+1:))
                call self%solve_tridiag(kr, lr, mr, a, bb, c, b, tcos, d, w)

                q(:mr, nlast) = q(:mr, nlast) + b(:mr)
                nlastp = nlast

                loop_197: do

                    jstep = jr
                    jr = i2r
                    i2r = i2r/2

                    if (jr == 0) then
                        w(1) = real(ipstor, kind=wp)
                        return
                    end if

                    jstart = 1 + jr
                    kr = kr - jr

                    if (nlast + jr <= n) then
                        kr = kr - jr
                        nlast = nlast + jr
                        jstop = nlast - jstep
                    else
                        jstop = nlast - jr
                    end if

                    lr = kr - jr

                    call self%generate_cosines(jr, 1, HALF, ZERO, tcos)

                    do j = jstart, jstop, jstep
                        jm2 = j - jr
                        jp2 = j + jr

                        if (j == jr) then
                            b(:mr) = q(:mr, j) + q(:mr, jp2)
                        else
                            b(:mr) = q(:mr, j) + q(:mr, jm2) + q(:mr, jp2)
                        end if

                        if (jr == 1) then
                            q(:mr, j) = ZERO
                        else
                            jm1 = j - i2r
                            jp1 = j + i2r
                            q(:mr, j) = HALF * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
                        end if

                        call self%solve_tridiag(jr, 0, mr, a, bb, c, b, tcos, d, w)

                        q(:mr, j) = q(:mr, j) + b(:mr)
                    end do

                    if (nlast + i2r <= n) then
                        nrod = 0
                    else
                        nrod = 1
                    end if

                    if (nlastp /= nlast) cycle loop_188

                end do loop_197
            end do loop_188

        end associate

    end subroutine poistg_solve_poisson_on_staggered_grid

end module type_StaggeredCyclicReductionUtility
