!
!     file pois3d.f90
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
! SUBROUTINE pois3d(lperod, l, c1, mperod, m, c2, nperod, n, &
!            a, b, c, ldimf, mdimf, f, ierror)
!
!
! DIMENSION OF           a(n), b(n), c(n), f(ldimf, mdimf, n)
! ARGUMENTS
!
! LATEST REVISION        April 2016
!
! PURPOSE                Solves the linear system of equations
!                        for unknown x values, where i=1, 2, ..., l,
!                        j=1, 2, ..., m, and k=1, 2, ..., n
!
!                        c1*(x(i-1, j, k) -2.*x(i, j, k) +x(i+1, j, k)) +
!                        c2*(x(i, j-1, k) -2.*x(i, j, k) +x(i, j+1, k)) +
!                        a(k)*x(i, j, k-1) +b(k)*x(i, j, k)+ c(k)*x(i, j, k+1)
!                        = f(i, j, k)
!
!                        the indices k-1 and k+1 are evaluated modulo n,
!                        i.e. x(i, j, 0)=x(i, j, n) and x(i, j, n+1)=x(i, j, 1).
!                        the unknowns
!                        x(0, j, k), x(l+1, j, k), x(i, 0, k), and x(i, m+1, k)
!                        are assumed to take on certain prescribed
!                        values described below.
!
! USAGE                  call pois3d (lperod, l, c1, mperod, m, c2, nperod,
!                        n, a, b, c, ldimf, mdimf, f, ierror)
!
! ARGUMENTS
!
! ON INPUT
!                        lperod
!                          indicates the values that x(0, j, k) and
!                          x(l+1, j, k) are assumed to have.
!                          = 0  x(0, j, k)=x(l, j, k), x(l+1, j, k)=x(1, j, k)
!                          = 1  x(0, j, k) = 0,      x(l+1, j, k) = 0
!                          = 2  x(0, j, k)=0,        x(l+1, j, k)=x(l-1, j, k)
!                          = 3  x(0, j, k)=x(2, j, k), x(l+1, j, k)=x(l-1, j, k)
!                          = 4  x(0, j, k)=x(2, j, k), x(l+1, j, k) = 0.
!
!                        l
!                          the number of unknowns in the i-direction.
!                          l must be at least 3.
!
!                        c1
!                          real constant in the above linear system
!                          of equations to be solved.
!
!                        mperod
!                          indicates the values that x(i, 0, k) and
!                          x(i, m+1, k) are assumed to have.
!                          = 0  x(i, 0, k)=x(i, m, k), x(i, m+1, k)=x(i, 1, k)
!                          = 1  x(i, 0, k)=0,        x(i, m+1, k)=0
!                          = 2  x(i, 0, k)=0,        x(i, m+1, k)=x(i, m-1, k)
!                          = 3  x(i, 0, k)=x(i, 2, k)  x(i, m+1, k)=x(i, m-1, k)
!                          = 4  x(i, 0, k)=x(i, 2, k)  x(i, m+1, k)=0
!
!                        m
!                          the number of unknowns in the j-direction.
!                          m must be at least 3.
!
!                        c2
!                          real constant in the above linear system
!                          of equations to be solved.
!
!                        nperod
!                          = 0  if a(1) and c(n) are not zero.
!                          = 1  if a(1) = c(n) = 0.
!
!                        n
!                          the number of unknowns in the k-direction.
!                          n must be at least 3.
!
!                        a, b, c
!                          one-dimensional arrays of length n that
!                          specify the coefficients in the linear
!                          equations given above.
!
!                          if nperod = 0 the array elements must not
!                          depend upon index k, but must be constant.
!                          specifically, the subroutine checks the
!                          following condition
!                            a(k) = c(1)
!                            c(k) = c(1)
!                            b(k) = b(1)
!                          for k=1, 2, ..., n.
!
!                        ldimf
!                          The row (or first) dimension of the three-
!                          dimensional array f as it appears in the
!                          program calling pois3d.  this parameter is
!                          used to specify the variable dimension
!                          of f.  ldimf must be at least l.
!
!                        mdimf
!                          The column (or second) dimension of the three
!                          dimensional array f as it appears in the
!                          program calling pois3d.  this parameter is
!                          used to specify the variable dimension
!                          of f.  mdimf must be at least m.
!
!                        f
!                          A three-dimensional array that specifies the
!                          values of the right side of the linear system
!                          of equations given above.  f must be
!                          dimensioned at least l x m x n.
!
! ON OUTPUT
!
!                        f
!                          Contains the solution x.
!
!                        ierror
!                          An error flag that indicates invalid input
!                          parameters.  except for number zero, a
!                          solution is not attempted.
!                          = 0  no error
!                          = 1  if lperod < 0 or > 4
!                          = 2  if l < 3
!                          = 3  if mperod < 0 or > 4
!                          = 4  if m < 3
!                          = 5  if nperod < 0 or > 1
!                          = 6  if n < 3
!                          = 7  if ldimf < l
!                          = 8  if mdimf < m
!                          = 9  if a(k) /= c(1) or c(k) /= c(1)
!                               or b(i) /=b(1) for some k=1, 2, ..., n.
!                          = 10 if nperod = 1 and a(1) /= 0
!                               or c(n) /= 0
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
!
!                          Since this is the only means of indicating a
!                          possibly incorrect call to pois3d, the user
!                          should test ierror after the call.
!
! SPECIAL CONDITIONS     None
!
! I/O                    None
!
! PRECISION              64-bit float and 32-bit integer
!
! REQUIRED files         type_FishpackWorkspace.f90, comf.f90, type_FFTpack.f90
!
! STANDARD               Fortran 2008
!
! HISTORY                * Written by Roland Sweet at NCAR in the late
!                          1970's.  released on NCAR's public software
!                          libraries in January, 1980.
!                        * Revised in June 2004 by John Adams using
!                          Fortran 90 dynamically allocated work space.
!                        * Revised in April 2016 using features from
!                          Fortran 2008
!
! ALGORITHM              This subroutine solves three-dimensional block
!                        tridiagonal linear systems arising from finite
!                        difference approximations to three-dimensional
!                        poisson equations using the FFT package
!                        fftpack written by Paul Swarztrauber.
!
! TIMING                 For large l, m and n, the operation count
!                        is roughly proportional to
!                          l*m*n*(log2(l)+log2(m)+5)
!                        but also depends on input parameters lperod
!                        and mperod.
!
! ACCURACY               To measure the accuracy of the algorithm a
!                        uniform random number generator was used to
!                        create a solution array x for the system given
!                        in the 'purpose' section with
!                          a(k) = c(k) = -0.5_wp *b(k) = 1,  k=1, 2, ..., n
!                        and, when nperod = 1
!                          a(1) = c(n) = 0
!                          a(n) = c(1) = 2.
!
!                        The solution x was substituted into the given
!                        system and, using double precision, a right
!                        side y was computed.  using this array y
!                        subroutine pois3d was called to produce an
!                        approximate solution z.  Relative error
!
!                         = max(abs(z(i, j, k)-x(i, j, k)))/max(abs(x(i, j, k
!
!                        was computed, where the two maxima are taken
!                        over i=1, 2, ..., l, j=1, 2, ..., m and k=1, 2, ..., n.
!                        values of e are given in the table below for
!                        some typical values of l, m and n.
!
!                        l(=m=n)   lperod    mperod       e
!                        ------    ------    ------     ------
!
!                          16        0         0        1.e-13
!                          15        1         1        4.e-13
!                          17        3         3        2.e-13
!                          32        0         0        2.e-13
!                          31        1         1        2.e-12
!                          33        3         3        7.e-13
!
! REFERENCES              None
!
module module_pois3d

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_FFTpack, only: &
        FFTpack

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: pois3d
    public :: pois3dd


contains


    subroutine pois3d( lperod, l, c1, mperod, m, c2, nperod, n, a, b, c, &
        ldimf, mdimf, f, ierror)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: lperod
        integer (ip), intent (in)     :: l
        integer (ip), intent (in)     :: mperod
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: nperod
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: ldimf
        integer (ip), intent (in)     :: mdimf
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in)     :: c1
        real (wp),    intent (in)     :: c2
        real (wp),    intent (inout ) :: a(:)
        real (wp),    intent (in out) :: b(:)
        real (wp),    intent (in out) :: c(:)
        real (wp),    intent (in out) :: f(:,:,:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip)             :: irwk, icwk
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------

        !
        !==> Compute required workspace dimensions
        !
        irwk = 30+l+m+2*n + max(l, m, n) + 7 * ( (l+1)/2+(m+1)/2)
        icwk = 0

        !
        !==> Allocate memory
        !
        call workspace%create(irwk, icwk)

        !
        !==> Solve system
        !
        associate( rew => workspace%real_workspace )

            call pois3dd(lperod, l, c1, mperod, m, c2, nperod, n, &
                a, b, c, ldimf, mdimf, f, ierror, rew)

        end associate

        !
        !==> Release memory
        !
        call workspace%destroy()

    end subroutine pois3d


    subroutine pois3dd(lperod, l, c1, mperod, m, c2, nperod, n, a, b, &
        c, ldimf, mdimf, f, ierror, w)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: lperod
        integer (ip), intent (in)     :: l
        integer (ip), intent (in)     :: mperod
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: nperod
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: ldimf
        integer (ip), intent (in)     :: mdimf
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in)     :: c1
        real (wp),    intent (in)     :: c2
        real (wp),    intent (in out) :: a(*)
        real (wp),    intent (in out) :: b(*)
        real (wp),    intent (in out) :: c(*)
        real (wp),    intent (in out) :: f(ldimf, mdimf, *)
        real (wp),    intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: nh, nhm1, nodd, i, j, k
        integer (ip) :: workspace_indices(6)
        real (wp)    :: temp_save(6)
        !-----------------------------------------------

        !
        !==> Check input arguments
        !
        call check_input_arguments(lperod, l, mperod, m, nperod, n, &
            a, b, c, ldimf, mdimf, ierror)

        ! Check error flag
        if (ierror /= 0) then
            return
        end if

        !
        !==> Compute workspace indices
        !
        workspace_indices = get_pois3dd_workspace_indices(l, m, n)

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
                !
                !==> Reorder unknowns when nperod = 0
                !
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

                        w(nh) = 2.0_wp*f(i, j, nh)

                        select case (nodd)
                            case (1)
                                f(i, j, :n) = w(:n)
                            case (2)
                                w(n) = 2.0_wp * f(i, j, n)
                                f(i, j, :n) = w(:n)
                        end select

                    end do
                end do

                temp_save(1) = c(nhm1)
                temp_save(2) = a(nh)
                temp_save(3) = c(nh)
                temp_save(4) = b(nhm1)
                temp_save(5) = b(n)
                temp_save(6) = a(n)

                c(nhm1) = 0.0_wp
                a(nh) = 0.0_wp
                c(nh) = 2.0_wp*c(nh)

                select case (nodd)
                    case default
                        b(nhm1) = b(nhm1) - a(nh-1)
                        b(n) = b(n) + a(n)
                    case (2)
                        a(n) = c(nh)
                end select

            end if

            call pos3d1(lp, l, mp, m, n, a, b, c, ldimf, mdimf, f, &
                w, w(iwyrt), w(iwt), w(iwd), w(iwx), w(iwy), c1, c2, w(iwbb))

            select case (np)
                case (1)
                    do i = 1, l
                        do j = 1, m
                            w(nh-1:nh-nhm1:(-1))= &
                                0.5_wp *(f(i, j, nh+1:nhm1+nh)+f(i, j, :nhm1))
                            w(nh+1:nhm1+nh) = &
                                0.5_wp *(f(i, j, nh+1:nhm1+nh)-f(i, j, :nhm1))
                            w(nh) = 0.5_wp *f(i, j, nh)
                            select case (nodd)
                                case (1)
                                    f(i, j, :n) = w(:n)
                                case (2)
                                    w(n) = 0.5_wp * f(i, j, n)
                                    f(i, j, :n) = w(:n)
                            end select
                        end do
                    end do

                    c(nhm1) = temp_save(1)
                    a(nh) = temp_save(2)
                    c(nh) = temp_save(3)
                    b(nhm1) = temp_save(4)
                    b(n) = temp_save(5)
                    a(n) = temp_save(6)

                case (2)
                    return
            end select

        end associate


    contains


        pure subroutine check_input_arguments(lperod, l, mperod, m, nperod, n, &
            a, b, c, ldimf, mdimf, ierror)
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in)     :: lperod
            integer (ip), intent (in)     :: l
            integer (ip), intent (in)     :: mperod
            integer (ip), intent (in)     :: m
            integer (ip), intent (in)     :: nperod
            integer (ip), intent (in)     :: n
            integer (ip), intent (in)     :: ldimf
            integer (ip), intent (in)     :: mdimf
            integer (ip), intent (out)    :: ierror
            real (wp),    intent (in out) :: a(*)
            real (wp),    intent (in out) :: b(*)
            real (wp),    intent (in out) :: c(*)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip) :: k !! Counter
            !-----------------------------------------------

            associate( &
                lp => lperod + 1, &
                mp => mperod + 1, &
                np => nperod + 1 &
                )

                if (lp < 1 .or. lp > 5) then
                    ierror = 1
                    return
                else if (l < 3) then
                    ierror = 2
                    return
                else if (mp < 1 .or. mp > 5) then
                    ierror = 3
                    return
                else if (m < 3) then
                    ierror = 4
                    return
                else if (np < 1 .or. np > 2) then
                    ierror = 5
                    return
                else if (n < 3) then
                    ierror = 6
                else if (ldimf < l) then
                    ierror = 7
                    return
                else if (mdimf < m) then
                    ierror = 8
                    return
                else if (np == 1) then
                    do k = 1, n
                        if (a(k) /= c(1)) then
                            ierror = 9
                            return
                        end if
                        if (c(k) /= c(1)) then
                            ierror = 9
                            return
                        end if
                        if (b(k) /= b(1)) then
                            ierror = 9
                            return
                        end if
                    end do
                else if (nperod == 1 .and. (a(1) /= 0.0_wp .or. c(n) /= 0.0_wp)) then
                    ierror = 10
                    return
                else
                    ierror = 0
                end if

            end associate

        end subroutine check_input_arguments



        pure function get_pois3dd_workspace_indices(l, m, n) result (return_value)
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            integer (ip), intent (in)     :: l
            integer (ip), intent (in)     :: m
            integer (ip), intent (in)     :: n
            integer (ip)                  :: return_value(6)
            !--------------------------------------------------------------------------------

            associate( i => return_value )

                i(1) = l + 1
                i(2) = i(1) + m
                i(3) = i(2) + max(l, m, n) + 1
                i(4) = i(3) + n
                i(5) = i(4) + n
                i(6) = i(5) + 7*((l + 1)/2) + 15

            end associate

        end function get_pois3dd_workspace_indices

    end subroutine pois3dd


    subroutine pos3d1(lp, l, mp, m, n, a, b, c, ldimf, mdimf, f, xrt, &
        yrt, t, d, wx, wy, c1, c2, bb)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: lp
        integer (ip), intent (in)     :: l
        integer (ip), intent (in)     :: mp
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: ldimf
        integer (ip), intent (in)     :: mdimf
        real (wp),    intent (in)     :: c1
        real (wp),    intent (in)     :: c2
        real (wp),    intent (in out) :: a(*)
        real (wp),    intent (in)     :: b(*)
        real (wp),    intent (in out) :: c(*)
        real (wp),    intent (in out) :: f(ldimf, mdimf, 1)
        real (wp),    intent (in out) :: xrt(*)
        real (wp),    intent (in out) :: yrt(*)
        real (wp),    intent (in out) :: t(*)
        real (wp),    intent (in out) :: d(*)
        real (wp),    intent (in out) :: wx(*)
        real (wp),    intent (in out) :: wy(*)
        real (wp),    intent (in out) :: bb(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip)         :: lr, mr, nr, lrdel, i, mrdel, j, ifwrd, is, k
        real (wp), parameter :: PI = acos(-1.0_wp)
        real (wp)            :: scalx, dx, di, scaly, dy, dj
        type (FFTpack)       :: fft
        !-----------------------------------------------

        lr = l
        mr = m
        nr = n
        !
        !==> generate transform roots
        !
        lrdel = ((lp - 1)*(lp - 3)*(lp - 5))/3
        scalx = lr + lrdel
        dx = PI/(2.0_wp*scalx)

        select case (lp)
            case (1)
                xrt(1) = 0.0_wp
                xrt(lr) = -4.0_wp*c1

                do i = 3, lr, 2
                    xrt(i-1) = -4.0_wp*c1*sin(real(i - 1)*dx)**2
                    xrt(i) = xrt(i-1)
                end do

                call fft%rffti(lr, wx)

            case default

                select case (lp)
                    case (2)
                        di = 0.0_wp
                    case (3, 5)
                        di = 0.5_wp
                        scalx = 2.0_wp*scalx
                    case (4)
                        di = 1.0_wp
                end select

                do i = 1, lr
                    xrt(i) = -4.0_wp*c1*sin((real(i, kind=wp) - di)*dx)**2
                end do

                scalx = 2.0_wp*scalx

                if ( lp /= 1 ) then
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

                            xrt(1) = 0.0_wp
                            xrt(lr) = -4.0_wp*c1

                            do i = 3, lr, 2
                                xrt(i-1) = -4.0_wp*c1*sin(real(i - 1)*dx)**2
                                xrt(i) = xrt(i-1)
                            end do

                            call fft%rffti(lr, wx)
                    end select
                end if
        end select

        mrdel = ((mp - 1)*(mp - 3)*(mp - 5))/3
        scaly = mr + mrdel
        dy = PI/(2.0_wp*scaly)

        select case (mp)
            case (1)
                yrt(1) = 0.0_wp
                yrt(mr) = -4.0_wp*c2

                do j = 3, mr, 2
                    yrt(j-1) = -4.0_wp*c2*sin(real(j - 1, kind=wp)*dy)**2
                    yrt(j) = yrt(j-1)
                end do

                call fft%rffti(mr, wy)

            case default

                select case (mp)
                    case (2)
                        dj = 0.0_wp
                    case (3, 5)
                        dj = 0.5_wp
                        scaly = 2.0_wp*scaly
                    case (4)
                        dj = 1.0_wp
                end select

                do j = 1, mr
                    yrt(j) = -4.0_wp*c2*sin((real(j, kind=wp) - dj)*dy)**2
                end do

                scaly = 2.0_wp * scaly

                select case (mp)
                    case (1)
                        yrt(1) = 0.0_wp
                        yrt(mr) = -4.0_wp*c2

                        do j = 3, mr, 2
                            yrt(j-1) = -4.0_wp*c2*sin(real(j - 1, kind=wp)*dy)**2
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
            !
            !==> Transform x
            !
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
                f(:lr, :mr, :nr) = f(:lr, :mr, :nr)/(scalx*scaly)
                return
            end if

            y_transform: do
                !
                !==> Transform y
                !
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
                                t(:nr) = f(i, j, :nr)

                                call trid(nr, a, bb, c, t, d)

                                f(i, j, :nr) = t(:nr)
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

        f(:lr, :mr, :nr) = f(:lr, :mr, :nr)/(scalx*scaly)

    contains


        subroutine trid(mr, a, b, c, y, d)

            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in)     :: mr
            real (wp),    intent (in)     :: a(*)
            real (wp),    intent (in)     :: b(*)
            real (wp),    intent (in)     :: c(*)
            real (wp),    intent (in out) :: y(*)
            real (wp),    intent (in out) :: d(*)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip) :: m, mm1, i, iip
            real (wp)    :: z
            !-----------------------------------------------

            m = mr
            mm1 = m - 1
            z = 1.0_wp/b(1)
            d(1) = c(1)*z
            y(1) = y(1)*z

            do i = 2, mm1
                z = 1.0_wp/(b(i)-a(i)*d(i-1))
                d(i) = c(i)*z
                y(i) = (y(i)-a(i)*y(i-1))*z
            end do

            z = b(m) - a(m)*d(mm1)

            if (z == 0.0_wp) then
                y(m) = 0.0_wp
            else
                y(m) = (y(m)-a(m)*y(mm1))/z
            end if

            do iip = 1, mm1
                i = m - iip
                y(i) = y(i) - d(i)*y(i+1)
            end do

        end subroutine trid

    end subroutine pos3d1

end module module_pois3d
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
! April     2016    Fortran 2008 changes
!
