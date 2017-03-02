module three_dimensional_solvers

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        PI

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_PeriodicFastFourierTransform, only: &
        PeriodicFastFourierTransform

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: pois3d
    public :: hw3crt
    public :: SolverUtility3D

    ! Parameters confined to the module
    real(wp),    parameter :: ZERO = 0.0_wp
    real(wp),    parameter :: HALF = 0.5_wp
    real(wp),    parameter :: ONE = 1.0_wp
    real(wp),    parameter :: TWO = 2.0_wp
    real(wp),    parameter :: FOUR = 4.0_wp
    integer(ip), parameter :: SIZE_OF_WORKSPACE_INDICES = 6

    type, public :: SolverUtility3D
        ! Type components
        integer(ip) :: IIWK = SIZE_OF_WORKSPACE_INDICES
    contains
        ! Type-bound procedures
        procedure, nopass :: pois3dd
        procedure, nopass :: check_input_arguments
        procedure, nopass :: get_workspace_indices
    end type SolverUtility3D

    ! Declare interfaces for submodule implementation
    interface
        module subroutine pois3d(lperod, l, c1, mperod, m, c2, nperod, n, a, b, c, &
            ldimf, mdimf, f, ierror)

            ! Dummy arguments
            integer(ip), intent(in)    :: lperod
            integer(ip), intent(in)    :: l
            integer(ip), intent(in)    :: mperod
            integer(ip), intent(in)    :: m
            integer(ip), intent(in)    :: nperod
            integer(ip), intent(in)    :: n
            integer(ip), intent(in)    :: ldimf
            integer(ip), intent(in)    :: mdimf
            integer(ip), intent(out)   :: ierror
            real(wp),    intent(in)    :: c1
            real(wp),    intent(in)    :: c2
            real(wp),    intent(inout) :: a(:)
            real(wp),    intent(inout) :: b(:)
            real(wp),    intent(inout) :: c(:)
            real(wp),    intent(inout) :: f(:,:,:)
        end subroutine pois3d

        module subroutine hw3crt(xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, mbdcnd, &
            bdys, bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, ldimf, &
            mdimf, f, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)     :: l
            integer(ip), intent(in)     :: lbdcnd
            integer(ip), intent(in)     :: m
            integer(ip), intent(in)     :: mbdcnd
            integer(ip), intent(in)     :: n
            integer(ip), intent(in)     :: nbdcnd
            integer(ip), intent(in)     :: ldimf
            integer(ip), intent(in)     :: mdimf
            integer(ip), intent(out)    :: ierror
            real(wp),    intent(in)     :: xs
            real(wp),    intent(in)     :: xf
            real(wp),    intent(in)     :: ys
            real(wp),    intent(in)     :: yf
            real(wp),    intent(in)     :: zs
            real(wp),    intent(in)     :: zf
            real(wp),    intent(in)     :: elmbda
            real(wp),    intent(out)    :: pertrb
            real(wp),    intent(in)     :: bdxs(:,:)
            real(wp),    intent(in)     :: bdxf(:,:)
            real(wp),    intent(in)     :: bdys(:,:)
            real(wp),    intent(in)     :: bdyf(:,:)
            real(wp),    intent(in)     :: bdzs(:,:)
            real(wp),    intent(in)     :: bdzf(:,:)
            real(wp),    intent(inout)  :: f(:,:,:)
        end subroutine hw3crt
    end interface

contains

    subroutine pois3dd(lperod, l, c1, mperod, m, c2, nperod, n, a, b, &
        c, ldimf, mdimf, f, ierror, w, workspace_indices)

        ! Dummy arguments
        integer(ip), intent(in)    :: lperod
        integer(ip), intent(in)    :: l
        integer(ip), intent(in)    :: mperod
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: nperod
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: ldimf
        integer(ip), intent(in)    :: mdimf
        integer(ip), intent(out)   :: ierror
        integer(ip), intent(in)    :: workspace_indices(:)
        real(wp),    intent(in)    :: c1
        real(wp),    intent(in)    :: c2
        real(wp),    intent(inout) :: a(n)
        real(wp),    intent(inout) :: b(n)
        real(wp),    intent(inout) :: c(n)
        real(wp),    intent(inout) :: f(:,:,:)
        real(wp),    intent(inout) :: w(:)

        ! Local variables
        integer(ip) :: nh, nhm1, nodd, i, j, k
        real(wp)    :: temp_save(SIZE_OF_WORKSPACE_INDICES)

        ! Check input arguments
        call check_input_arguments(lperod, l, mperod, m, nperod, n, &
            a, b, c, ldimf, mdimf, ierror)

        ! Check error flag
        if (ierror /= 0) return

        associate( &
            lp => lperod + 1, &
            mp => mperod + 1, &
            np => nperod + 1, &
            iwyrt => workspace_indices(1), &
            iwt => workspace_indices(2), &
            iwd => workspace_indices(3), &
            iwbb => workspace_indices(4), &
            iwx => workspace_indices(5), &
            iwy => workspace_indices(6) &
            )

            if (np == 1) then
                ! Reorder unknowns when nperod = 0
                nh = (n + 1)/2
                nhm1 = nh - 1

                if (2*nh == n) then
                    nodd = 2
                else
                    nodd = 1
                end if

                do i = 1, l
                    do j = 1, m
                        do k = 1, nhm1
                            w(k) = f(i, j, nh-k) - f(i, j, k+nh)
                            w(k+nh) = f(i, j, nh-k) + f(i, j, k+nh)
                        end do

                        w(nh) = TWO*f(i, j, nh)

                        select case (nodd)
                            case (1)
                                f(i, j,:n) = w(:n)
                            case (2)
                                w(n) = TWO * f(i, j, n)
                                f(i, j,:n) = w(:n)
                        end select

                    end do
                end do

                temp_save(1) = c(nhm1)
                temp_save(2) = a(nh)
                temp_save(3) = c(nh)
                temp_save(4) = b(nhm1)
                temp_save(5) = b(n)
                temp_save(6) = a(n)

                c(nhm1) = ZERO
                a(nh) = ZERO
                c(nh) = TWO*c(nh)

                select case (nodd)
                    case (2)
                        a(n) = c(nh)
                    case default
                        b(nhm1) = b(nhm1) - a(nh-1)
                        b(n) = b(n) + a(n)
                end select
            end if

            call pois3dd_lower_routine(lp, l, mp, m, n, a, b, c, f, &
                w, w(iwyrt:), w(iwt:), w(iwd:), w(iwx:), w(iwy:), c1, c2, w(iwbb:))

            if (np == 1) then
                do i = 1, l
                    do j = 1, m
                        w(nh-1:nh-nhm1:(-1))= &
                            HALF *(f(i, j, nh+1:nhm1+nh)+f(i, j,:nhm1))
                        w(nh+1:nhm1+nh) = &
                            HALF *(f(i, j, nh+1:nhm1+nh)-f(i, j,:nhm1))
                        w(nh) = HALF *f(i, j, nh)
                        select case (nodd)
                            case (1)
                                f(i, j,:n) = w(:n)
                            case (2)
                                w(n) = HALF * f(i, j, n)
                                f(i, j,:n) = w(:n)
                        end select
                    end do
                end do

                c(nhm1) = temp_save(1)
                a(nh) = temp_save(2)
                c(nh) = temp_save(3)
                b(nhm1) = temp_save(4)
                b(n) = temp_save(5)
                a(n) = temp_save(6)
            end if
        end associate

    end subroutine pois3dd

    pure subroutine check_input_arguments(lperod, l, mperod, m, nperod, n, &
        a, b, c, ldimf, mdimf, ierror)

        ! Dummy arguments
        integer(ip), intent(in)   :: lperod
        integer(ip), intent(in)   :: l
        integer(ip), intent(in)   :: mperod
        integer(ip), intent(in)   :: m
        integer(ip), intent(in)   :: nperod
        integer(ip), intent(in)   :: n
        integer(ip), intent(in)   :: ldimf
        integer(ip), intent(in)   :: mdimf
        integer(ip), intent(out)  :: ierror
        real(wp),    intent(in)   :: a(n)
        real(wp),    intent(in)   :: b(n)
        real(wp),    intent(in)   :: c(n)

        associate( &
            lp => lperod + 1, &
            mp => mperod + 1, &
            np => nperod + 1 &
            )

            if (lp < 1 .or. lp > 5) then
                ierror = 1
            else if (l < 3) then
                ierror = 2
            else if (mp < 1 .or. mp > 5) then
                ierror = 3
            else if (m < 3) then
                ierror = 4
            else if (np < 1 .or. np > 2) then
                ierror = 5
            else if (n < 3) then
                ierror = 6
            else if (ldimf < l) then
                ierror = 7
            else if (mdimf < m) then
                ierror = 8
            else if (np == 1) then
                if (any(a /= c(1)) .or. any(c /= c(1)) .or. any(b /= b(1))) then
                    ierror = 9
                end if
            else if (nperod == 1 .and. (a(1) /= ZERO .or. c(n) /= ZERO)) then
                ierror = 10
            else
                ierror = 0
            end if
        end associate

    end subroutine check_input_arguments

    pure function get_workspace_indices(l, m, n) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)     :: l
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: n
        integer(ip)                 :: return_value(SIZE_OF_WORKSPACE_INDICES)

        associate( indx => return_value )
            indx(1) = l + 1
            indx(2) = indx(1) + m
            indx(3) = indx(2) + max(l, m, n) + 1
            indx(4) = indx(3) + n
            indx(5) = indx(4) + n
            indx(6) = indx(5) + 7*((l + 1)/2) + 15
        end associate

    end function get_workspace_indices

    subroutine pois3dd_lower_routine(lp, l, mp, m, n, a, b, c, &
        f, xrt, yrt, t, d, wx, wy, c1, c2, bb)

        ! Dummy arguments
        integer(ip), intent(in)    :: lp
        integer(ip), intent(in)    :: l
        integer(ip), intent(in)    :: mp
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: n
        real(wp),    intent(in)    :: c1
        real(wp),    intent(in)    :: c2
        real(wp),    intent(inout) :: a(n)
        real(wp),    intent(in)    :: b(n)
        real(wp),    intent(inout) :: c(n)
        real(wp),    intent(inout) :: f(:,:,:)
        real(wp),    intent(inout) :: xrt(:)
        real(wp),    intent(inout) :: yrt(:)
        real(wp),    intent(inout) :: t(:)
        real(wp),    intent(inout) :: d(:)
        real(wp),    intent(inout) :: wx(:)
        real(wp),    intent(inout) :: wy(:)
        real(wp),    intent(inout) :: bb(:)

        ! Local variables
        integer(ip) :: lr, mr, nr, lrdel, i, mrdel, j, ifwrd, is, k
        real(wp)    :: scalx, dx, di, scaly, dy, dj
        type(PeriodicFastFourierTransform)  :: fft

        lr = l
        mr = m
        nr = n

        ! generate transform roots
        lrdel = ((lp - 1)*(lp - 3)*(lp - 5))/3
        scalx = lr + lrdel
        dx = PI/(TWO*scalx)

        select case (lp)
            case (1)
                xrt(1) = ZERO
                xrt(lr) = -FOUR*c1

                do i = 3, lr, 2
                    xrt(i-1) = -FOUR*c1*sin(real(i - 1, kind=wp)*dx)**2
                    xrt(i) = xrt(i-1)
                end do

                call fft%rffti(lr, wx)

            case default

                select case (lp)
                    case (2)
                        di = ZERO
                    case (3, 5)
                        di = HALF
                        scalx = TWO*scalx
                    case (4)
                        di = ONE
                end select

                do i = 1, lr
                    xrt(i) = -FOUR*c1*sin((real(i, kind=wp) - di)*dx)**2
                end do

                scalx = TWO*scalx

                if (lp /= 1) then
                    select case (lp)
                        case (2)
                            call fft%sinti(lr, wx)
                        case (3)
                            call fft%sinqi(lr, wx)
                        case (4)
                            call fft%costi(lr, wx)
                        case (5)
                            call fft%cosqi(lr, wx)
                        case default

                            xrt(1) = ZERO
                            xrt(lr) = -FOUR*c1

                            do i = 3, lr, 2
                                xrt(i-1) = -FOUR*c1*sin(real(i - 1)*dx)**2
                                xrt(i) = xrt(i-1)
                            end do

                            call fft%rffti(lr, wx)
                    end select
                end if
        end select

        mrdel = ((mp - 1)*(mp - 3)*(mp - 5))/3
        scaly = mr + mrdel
        dy = PI/(TWO*scaly)

        select case (mp)
            case (1)
                yrt(1) = ZERO
                yrt(mr) = -FOUR*c2

                do j = 3, mr, 2
                    yrt(j-1) = -FOUR*c2*sin(real(j - 1, kind=wp)*dy)**2
                    yrt(j) = yrt(j-1)
                end do

                call fft%rffti(mr, wy)

            case default

                select case (mp)
                    case (2)
                        dj = ZERO
                    case (3, 5)
                        dj = HALF
                        scaly = TWO*scaly
                    case (4)
                        dj = ONE
                end select

                do j = 1, mr
                    yrt(j) = -FOUR*c2*sin((real(j, kind=wp) - dj)*dy)**2
                end do

                scaly = TWO * scaly

                select case (mp)
                    case (1)
                        yrt(1) = ZERO
                        yrt(mr) = -FOUR*c2

                        do j = 3, mr, 2
                            yrt(j-1) = -FOUR*c2*sin(real(j - 1, kind=wp)*dy)**2
                            yrt(j) = yrt(j-1)
                        end do

                        call fft%rffti(mr, wy)

                    case default

                        select case (mp)
                            case (2)
                                call fft%sinti(mr, wy)
                            case (3)
                                call fft%sinqi(mr, wy)
                            case (4)
                                call fft%costi(mr, wy)
                            case (5)
                                call fft%cosqi(mr, wy)
                        end select

                end select
        end select

        ifwrd = 1
        is = 1

        x_transform: do
            do j=1, mr
                do k=1, nr
                    do i=1, lr
                        t(i) = f(i, j, k)
                    end do
                    select case (lp)
                        case (1)
                            select case (ifwrd)
                                case (1)
                                    call fft%rfftf(lr, t, wx)
                                case (2)
                                    call fft%rfftb(lr, t, wx)
                            end select
                        case (2)
                            call fft%sint(lr, t, wx)
                        case (3)
                            select case (ifwrd)
                                case (1)
                                    call fft%sinqf(lr, t, wx)
                                case (2)
                                    call fft%sinqb(lr, t, wx)
                            end select
                        case (4)
                            call fft%cost(lr, t, wx)
                        case (5)
                            select case (ifwrd)
                                case (1)
                                    call fft%cosqf(lr, t, wx)
                                case (2)
                                    call fft%cosqb(lr, t, wx)
                            end select
                    end select

                    do i=1, lr
                        f(i, j, k) = t(i)
                    end do
                end do
            end do

            if (ifwrd == 2) then
                f(:lr,:mr,:nr) = f(:lr,:mr,:nr)/(scalx*scaly)
                return
            end if

            y_transform: do
                do i=1, lr
                    do k=1, nr
                        t(1:mr) = f(i,1:mr, k)
                        select case (mp)
                            case (1)
                                select case (ifwrd)
                                    case (1)
                                        call fft%rfftf(mr, t, wy)
                                    case (2)
                                        call fft%rfftb(mr, t, wy)
                                end select
                            case (2)
                                call fft%sint(mr, t, wy)
                            case (3)
                                select case (ifwrd)
                                    case (1)
                                        call fft%sinqf(mr, t, wy)
                                    case (2)
                                        call fft%sinqb(mr, t, wy)
                                end select
                            case (4)
                                call fft%cost(mr, t, wy)
                            case (5)
                                select case (ifwrd)
                                    case (1)
                                        call fft%cosqf(mr, t, wy)
                                    case (2)
                                        call fft%cosqb(mr, t, wy)
                                end select
                        end select
                        f(i,1:mr, k) = t(1:mr)
                    end do
                end do

                select case (ifwrd)
                    case (1)
                        do i = 1, lr
                            do j = 1, mr
                                bb(:nr) = b(:nr) + xrt(i) + yrt(j)
                                t(:nr) = f(i, j,:nr)

                                call solve_tridiag(nr, a, bb, c, t, d)

                                f(i, j,:nr) = t(:nr)
                            end do
                        end do
                        ifwrd = 2
                        is = -1
                        cycle y_transform
                    case (2)
                        cycle x_transform
                end select

                exit y_transform
            end do y_transform

            exit x_transform
        end do x_transform

        f(:lr,:mr,:nr) = f(:lr,:mr,:nr)/(scalx*scaly)

    end subroutine pois3dd_lower_routine

    subroutine solve_tridiag(m, a, b, c, y, d)

        ! Dummy arguments
        integer(ip), intent(in)    :: m
        real(wp),    intent(in)    :: a(m) ! subdiagonal, i.e., the diagonal below the main diagonal
        real(wp),    intent(in)    :: b(m) ! the main diagonal
        real(wp),    intent(in)    :: c(m) ! supdiagonal, i.e, the diagonal above the main diagonal
        real(wp),    intent(inout) :: y(m) ! the answer
        real(wp),    intent(inout) :: d(m) ! right part

        ! Local variables
        integer(ip) :: i
        real(wp)    :: temp

        ! Initialize
        d(1) = c(1)/b(1)
        y(1) = y(1)/b(1)

        do i = 2, m - 1
            temp = b(i) - d(i - 1) * a(i)
            d(i) = c(i)/temp
            y(i) = (y(i) - a(i) * y(i - 1))/temp
        end do

        temp = b(m) - a(m) * d(m - 1)

        if (temp == ZERO) then
            y(m) = ZERO
        else
            y(m) = (y(m) - a(m) * y(m - 1))/temp
        end if

        do i = m - 1, 1, -1
            y(i) = y(i) - d(i) * y(i + 1)
        end do

    end subroutine solve_tridiag

end module three_dimensional_solvers
