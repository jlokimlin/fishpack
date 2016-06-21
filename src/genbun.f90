!
!     file genbun.f90
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
! SUBROUTINE genbun(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
!
!
! DIMENSION OF           a(m), b(m), c(m), y(idimy, n)
! ARGUMENTS
!
! LATEST REVISION        April 2016
!
! PURPOSE                The name of this package is a mnemonic for the
!                        generalized buneman algorithm.
!
!                        It solves the real linear system of equations
!
!                        a(i)*x(i-1, j) + b(i)*x(i, j) + c(i)*x(i+1, j)
!                        + x(i, j-1) - 2.0*x(i, j) + x(i, j+1) = y(i, j)
!
!                        for i = 1, 2, ..., m  and  j = 1, 2, ..., n.
!
!                        indices i+1 and i-1 are evaluated modulo m,
!                        i.e., x(0, j) = x(m, j) and x(m+1, j) = x(1, j),
!                        and x(i, 0) may equal 0, x(i, 2), or x(i, n),
!                        and x(i, n+1) may equal 0, x(i, n-1), or x(i, 1)
!                        depending on an input parameter.
!
! USAGE                  call genbun(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
!
! ARGUMENTS
!
! ON INPUT               nperod
!
!                          Indicates the values that x(i, 0) and
!                          x(i, n+1) are assumed to have.
!
!                          = 0  if x(i, 0) = x(i, n) and x(i, n+1) =
!                               x(i, 1).
!                          = 1  if x(i, 0) = x(i, n+1) = 0  .
!                          = 2  if x(i, 0) = 0 and x(i, n+1) = x(i, n-1).
!                          = 3  if x(i, 0) = x(i, 2) and x(i, n+1) =
!                               x(i, n-1).
!                          = 4  if x(i, 0) = x(i, 2) and x(i, n+1) = 0.
!
!                        n
!                          The number of unknowns in the j-direction.
!                          n must be greater than 2.
!
!                        mperod
!                          = 0 if a(1) and c(m) are not zero
!                          = 1 if a(1) = c(m) = 0
!
!                        m
!                          The number of unknowns in the i-direction.
!                          n must be greater than 2.
!
!                        a, b, c
!                          One-dimensional arrays of length m that
!                          specify the coefficients in the linear
!                          equations given above.  if mperod = 0
!                          the array elements must not depend upon
!                          the index i, but must be constant.
!                          specifically, the subroutine checks the
!                          following condition .
!
!                            a(i) = c(1)
!                            c(i) = c(1)
!                            b(i) = b(1)
!
!                          for i=1, 2, ..., m.
!
!                        idimy
!                          The row (or first) dimension of the
!                          two-dimensional array y as it appears
!                          in the program calling genbun.
!                          this parameter is used to specify the
!                          variable dimension of y.
!                          idimy must be at least m.
!
!                        y
!                          A two-dimensional complex array that
!                          specifies the values of the right side
!                          of the linear system of equations given
!                          above.
!                          y must be dimensioned at least m*n.
!
!
!  ON OUTPUT             y
!
!                          Contains the solution x.
!
!                        ierror
!                          An error flag which indicates invalid
!                          input parameters  except for number
!                          zero, a solution is not attempted.
!
!                          = 0  no error
!                          = 1  m <= 2
!                          = 2  n <= 2
!                          = 3  idimy <= m
!                          = 4  nperod <= 0 or nperod > 4
!                          = 5  mperod < 0 or mperod > 1
!                          = 6  a(i) /= c(1) or c(i) /= c(1) or
!                               b(i) /= b(1) for
!                               some i=1, 2, ..., m.
!                          = 7  a(1) /= 0 or c(m) /= 0 and
!                                 mperod = 1
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
!
!
! SPECIAL CONDITONS      None
!
! I/O                    None
!
! PRECISION              64-bit precision float and 32-bit precision integer
!
! REQUIRED FILES         comf.f90, gnbnaux.f90, type_FishpackWorkspace.f90
!
!
! STANDARD               Fortran 2008
!
! HISTORY                Written in 1979 by Roland Sweet of NCAR'S
!                        scientific computing division. made available
!                        on NCAR'S public libraries in January, 1980.
!
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated work space.
!
! ALGORITHM              The linear system is solved by a cyclic
!                        reduction algorithm described in the
!                        reference.
!
! PORTABILITY            Fortran 2008 --
!                        the machine dependent constant pi is
!                        defined as acos(-1.0_wp) where wp = REAL64
!                        from the intrinsic module iso_fortran_env.
!
! REFERENCES             Sweet, R., "A cyclic reduction algorithm for
!                        solving block tridiagonal systems of arbitrary
!                        dimensions, " SIAM J. On Numer. Anal., 14 (1977)
!                        PP. 706-720.
!
! ACCURACY               This test was performed on a platform with
!                        64 bit floating point arithmetic.
!                        a uniform random number generator was used
!                        to create a solution array x for the system
!                        given in the 'purpose' description above
!                        with
!                          a(i) = c(i) = -0.5*b(i) = 1, i=1, 2, ..., m
!
!                        and, when mperod = 1
!
!                          a(1) = c(m) = 0
!                          a(m) = c(1) = 2.
!
!                        The solution x was substituted into the
!                        given system  and, using double precision
!                        a right side y was computed.
!                        using this array y, subroutine genbun
!                        was called to produce approximate
!                        solution z.  then relative error
!                          e = max(abs(z(i, j)-x(i, j)))/
!                              max(abs(x(i, j)))
!                        was computed, where the two maxima are taken
!                        over i=1, 2, ..., m and j=1, ..., n.
!
!                        The value of e is given in the table
!                        below for some typical values of m and n.
!
!                   m (=n)    mperod    nperod        e
!                   ------    ------    ------      ------
!
!                     31        0         0         6.e-14
!                     31        1         1         4.e-13
!                     31        1         3         3.e-13
!                     32        0         0         9.e-14
!                     32        1         1         3.e-13
!                     32        1         3         1.e-13
!                     33        0         0         9.e-14
!                     33        1         1         4.e-13
!                     33        1         3         1.e-13
!                     63        0         0         1.e-13
!                     63        1         1         1.e-12
!                     63        1         3         2.e-13
!                     64        0         0         1.e-13
!                     64        1         1         1.e-12
!                     64        1         3         6.e-13
!                     65        0         0         2.e-13
!                     65        1         1         1.e-12
!                     65        1         3         4.e-13
!
module module_genbun

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_gnbnaux, only: &
        GenbunAux

    ! Explicit typing only
    implicit None

    ! Everything is private unless stated otherwise
    private
    public :: genbun
    public :: genbunn


contains


    subroutine genbun(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip),          intent (in)     :: nperod
        integer (ip),          intent (in)     :: n
        integer (ip),          intent (in)     :: mperod
        integer (ip),          intent (in)     :: m
        integer (ip),          intent (in)     :: idimy
        integer (ip),          intent (out)    :: ierror
        real (wp), contiguous, intent (in)     :: a(:)
        real (wp), contiguous, intent (in)     :: b(:)
        real (wp), contiguous, intent (in)     :: c(:)
        real (wp),             intent (in out) :: y(idimy,n)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)             :: irwk
        type (FishpackWorkspace) :: workspace
        !--------------------------------------------------------------------------------

        !
        !==> Allocate memory
        !
        associate( icwk => 0 )

            ! Compute real workspace size for genbun algorithm
            call workspace%get_genbun_workspace_dimensions(n, m, irwk)

            call workspace%create(irwk, icwk)

        end associate

        !
        !==> Solve system
        !
        associate( rew => workspace%real_workspace )

            call genbunn(nperod, n, mperod, m, a, b, c, idimy, y, ierror, rew)

        end associate

        !
        !==> Release memory
        !
        call workspace%destroy()

    end subroutine genbun


    subroutine genbunn(nperod, n, mperod, m, a, b, c, idimy, y, ierror, w)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: nperod
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: mperod
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: idimy
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in)     :: a(m)
        real (wp),    intent (in)     :: b(m)
        real (wp),    intent (in)     :: c(m)
        real (wp),    intent (in out) :: y(idimy,n)
        real (wp),    intent (in out) :: w(*)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)  :: workspace_indices(12)
        integer (ip)  :: i, k, j, mp, np, irev, mh, mhm1, modd
        integer (ip)  :: nby2, mskip
        real (wp)     :: temp, ipstor
        !--------------------------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        call check_input_arguments(nperod, n, mperod, m, idimy, ierror, a, b, c)

        ! Check error flag
        if (ierror /= 0) return

        !
        !==> Compute workspace indices
        !
        workspace_indices = get_genbunn_workspace_indices(n, m)

        associate( &
            mp1 => workspace_indices(1), &
            iwba => workspace_indices(2), &
            iwbb => workspace_indices(3), &
            iwbc => workspace_indices(4), &
            iwb2 => workspace_indices(5), &
            iwb3 => workspace_indices(6), &
            iww1 => workspace_indices(7), &
            iww2 => workspace_indices(8), &
            iww3 => workspace_indices(9), &
            iwd => workspace_indices(10), &
            iwtcos => workspace_indices(11), &
            iwp => workspace_indices(12) &
            )

            w(iwba:m-1+iwba) = -a(:m)
            w(iwbc:m-1+iwbc) = -c(:m)
            w(iwbb:m-1+iwbb) = 2.0_wp - b(:m)
            y(:m, :n) = -y(:m, :n)

            mp = mperod + 1
            np = nperod + 1

            select case (mp)
                case (1)
                    mh = (m + 1)/2
                    mhm1 = mh - 1

                    if (mh*2 == m) then
                        modd = 2
                    else
                        modd = 1
                    end if

                    do j = 1, n
                        w(:mhm1) = y(mh-1:mh-mhm1:(-1), j) - y(mh+1:mhm1+mh, j)
                        w(mh+1:mhm1+mh) = y(mh-1:mh-mhm1:(-1), j) + y(mh+1:mhm1+mh, j)
                        w(mh) = 2.0_wp*y(mh, j)
                        select case (modd)
                            case (1)
                                y(:m, j) = w(:m)
                            case (2)
                                w(m) = 2.0_wp*y(m, j)
                                y(:m, j) = w(:m)
                        end select
                    end do

                    k = iwbc + mhm1 - 1
                    i = iwba + mhm1
                    w(k) = 0.0_wp
                    w(i) = 0.0_wp
                    w(k+1) = 2.0_wp*w(k+1)

                    select case (modd)
                        case default
                            k = iwbb + mhm1 - 1
                            w(k) = w(k) - w(i-1)
                            w(iwbc-1) = w(iwbc-1) + w(iwbb-1)
                        case (2)
                            w(iwbb-1) = w(k+1)
                    end select
            end select

            main_loop: do

                select case (np)
                    case (1)
                        call poisp2(m, n, w(iwba), w(iwbb), w(iwbc), &
                            y, idimy, w, w(iwb2), w(iwb3), w(iww1), &
                            w(iww2), w(iww3), w(iwd), w(iwtcos), w(iwp) )
                    case (2)
                        call poisd2(m, n, 1, w(iwba), w(iwbb), w(iwbc), &
                            y, idimy, w, w(iww1), w(iwd), w(iwtcos), w(iwp))
                    case (3)
                        call poisn2(m, n, 1, 2, w(iwba), w(iwbb), w(iwbc), &
                            y, idimy, w, w(iwb2), w(iwb3), w(iww1), w(iww2),  &
                            w(iww3), w(iwd), w(iwtcos), w(iwp))
                    case (4)
                        call poisn2(m, n, 1, 1, w(iwba), w(iwbb), w(iwbc), &
                            y, idimy, w, w(iwb2), w(iwb3), w(iww1), w(iww2), &
                            w(iww3), w(iwd), w(iwtcos), w(iwp))
                end select

                loop_113: do

                    select case (np)
                        case (1:4)
                            select case (mp)
                                case (1)
                                    do j = 1, n
                                        w(mh-1:mh-mhm1:(-1)) = 0.5_wp*(y(mh+1:mhm1+mh, j)+y(:mhm1, j))
                                        w(mh+1:mhm1+mh) = 0.5_wp*(y(mh+1:mhm1+mh, j)-y(:mhm1, j))
                                        w(mh) = 0.5_wp*y(mh, j)
                                        select case (modd)
                                            case (1)
                                                y(:m, j) = w(:m)
                                            case (2)
                                                w(m) = 0.5_wp*y(m, j)
                                                y(:m, j) = w(:m)
                                        end select
                                    end do
                                    w(1) = ipstor + real(iwp - 1, kind=wp)
                                    return
                                case (2)
                                    w(1) = ipstor + real(iwp - 1, kind=wp)
                                    return
                            end select

                            mh = (m + 1)/2
                            mhm1 = mh - 1

                            if (mh*2 == m) then
                                modd = 2
                            else
                                modd = 1
                            end if

                            do j = 1, n
                                w(:mhm1) = y(mh-1:mh-mhm1:(-1), j) - y(mh+1:mhm1+mh, j)
                                w(mh+1:mhm1+mh) = y(mh-1:mh-mhm1:(-1), j) + y(mh+1:mhm1+mh, j)
                                w(mh) = 2.0_wp*y(mh, j)
                                select case (modd)
                                    case (1)
                                        y(:m, j) = w(:m)
                                    case (2)
                                        w(m) = 2.0_wp*y(m, j)
                                        y(:m, j) = w(:m)
                                end select
                            end do

                            k = iwbc + mhm1 - 1
                            i = iwba + mhm1
                            w(k) = 0.0_wp
                            w(i) = 0.0_wp
                            w(k+1) = 2.0_wp*w(k+1)

                            select case (modd)
                                case default
                                    k = iwbb + mhm1 - 1
                                    w(k) = w(k) - w(i-1)
                                    w(iwbc-1) = w(iwbc-1) + w(iwbb-1)
                                case (2)
                                    w(iwbb-1) = w(k+1)
                            end select

                            cycle main_loop
                        !
                        !     reverse columns when nperod = 4.
                        !
                    end select

                    irev = 1
                    nby2 = n/2

                    loop_124: do

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
                                call poisn2(m, n, 1, 2, w(iwba), w(iwbb), w(iwbc), &
                                    y, idimy, w, w(iwb2), w(iwb3), w(iww1), w(iww2),  &
                                    w(iww3), w(iwd), w(iwtcos), w(iwp))

                                ipstor = w(iww1)
                                irev = 2

                                if (nperod == 4) cycle loop_124
                                cycle loop_113

                            case (2)
                                cycle loop_113
                        end select
                        exit loop_124
                    end do loop_124
                    exit loop_113
                end do loop_113
            end do main_loop

            do j = 1, n
                w(mh-1:mh-mhm1:(-1)) = 0.5_wp*(y(mh+1:mhm1+mh, j)+y(:mhm1, j))
                w(mh+1:mhm1+mh) = 0.5_wp*(y(mh+1:mhm1+mh, j)-y(:mhm1, j))
                w(mh) = 0.5_wp*y(mh, j)
                select case (modd)
                    case (1)
                        y(:m, j) = w(:m)
                    case (2)
                        w(m) = 0.5_wp*y(m, j)
                        y(:m, j) = w(:m)
                end select
            end do

            w(1) = ipstor + real(iwp - 1, kind=wp)

        end associate

    contains

        pure subroutine check_input_arguments(nperod, n, mperod, m, idimy, ierror, a, b, c)
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            integer (ip), intent (in)  :: nperod
            integer (ip), intent (in)  :: n
            integer (ip), intent (in)  :: mperod
            integer (ip), intent (in)  :: m
            integer (ip), intent (in)  :: idimy
            integer (ip), intent (out) :: ierror
            real (wp),    intent (in)  :: a(m)
            real (wp),    intent (in)  :: b(m)
            real (wp),    intent (in)  :: c(m)
            !--------------------------------------------------------------------------------

            if (m <= 2) then
                ierror = 1
                return
            else if (n <= 2) then
                ierror = 2
                return
            else if (idimy < m) then
                ierror = 3
                return
            else if (nperod < 0 .or. nperod > 4) then
                ierror = 4
                return
            else if (mperod < 0 .or. mperod > 1) then
                ierror = 5
                return
            else if (mperod /= 1) then
                if (any(a /= c(1)) .or. any(c /= c(1)) .or. any(b /= b(1))) then
                    ierror = 6
                    return
                end if
            else if (a(1) /= 0.0_wp .or. c(m) /= 0.0_wp) then
                ierror = 7
                return
            else
                ierror = 0
            end if

        end subroutine check_input_arguments


        pure function get_genbunn_workspace_indices(n, m) result (return_value)
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            integer (ip), intent (in) :: n
            integer (ip), intent (in) :: m
            integer (ip)              :: return_value(12)
            !--------------------------------------------------------------------------------
            integer (ip) :: j !! Counter
            !--------------------------------------------------------------------------------

            associate( i => return_value)

                i(1:2) = m + 1

                do j = 3, 11
                    i(j) = i(j-1) + m
                end do

                i(12) = i(11) + 4 * n
            end associate

        end function get_genbunn_workspace_indices


    end subroutine genbunn



    subroutine poisd2(mr, nr, istag, ba, bb, bc, q, idimq, b, w, d, tcos, p)
        !
        ! Purpose:
        !
        !     subroutine to solve poisson's equation for dirichlet boundary
        !     conditions.
        !
        !     istag = 1 if the last diagonal block is the matrix a.
        !     istag = 2 if the last diagonal block is the matrix a+i.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: mr
        integer (ip), intent (in)     :: nr
        integer (ip), intent (in)     :: istag
        integer (ip), intent (in)     :: idimq
        real (wp),    intent (in)     :: ba(*)
        real (wp),    intent (in)     :: bb(*)
        real (wp),    intent (in)     :: bc(*)
        real (wp),    intent (in out) :: q(idimq,*)
        real (wp),    intent (in out) :: b(*)
        real (wp),    intent (in out) :: w(*)
        real (wp),    intent (in out) :: d(*)
        real (wp),    intent (in out) :: tcos(*)
        real (wp),    intent (in out) :: p(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip)     :: m, n, jsh, ipp, ipstor, kr
        integer (ip)     :: irreg, jstsav, i, lr, nun, nodd, noddpr
        integer (ip)     :: jst, jsp, l, j, jm1, jp1, jm2, jp2, jm3, jp3
        integer (ip)     :: krpi, ideg, jdeg
        real (wp)        :: fi, t
        type (GenbunAux) :: genbun_aux
        !-----------------------------------------------

        m = mr
        n = nr
        jsh = 0
        fi = 1.0_wp/istag
        ipp = -m
        ipstor = 0

        select case (istag)
            case (2)
                kr = 1
                jstsav = 1
                irreg = 2
                if (n <= 1) then
                    tcos(1) = -1.0_wp
                    b(:m) = q(:m, 1)
                    call genbun_aux%trix(1, 0, m, ba, bb, bc, b, tcos, d, w)
                    q(:m, 1) = b(:m)
                    w(1) = ipstor
                    return
                end if
            case default
                kr = 0
                irreg = 1
                if (n <= 1) then
                    tcos(1) = 0.0_wp
                    b(:m) = q(:m, 1)
                    call genbun_aux%trix(1, 0, m, ba, bb, bc, b, tcos, d, w)
                    q(:m, 1) = b(:m)
                    w(1) = ipstor
                    return
                end if
        end select

        lr = 0
        p(:m) = 0.0_wp
        nun = n
        jst = 1
        jsp = n
        !
        !==> irreg = 1 when no irregularities have occurred, otherwise it is 2.
        !
        loop_108: do

            l = 2*jst
            nodd = 2 - 2*((nun + 1)/2) + nun
            !
            !==> nodd = 1 when nun is odd, otherwise it is 2.
            !
            select case (nodd)
                case (1)
                    jsp = jsp - jst
                    if (irreg /= 1) then
                        jsp = jsp - l
                    end if
                case default
                    jsp = jsp - l
            end select

            call genbun_aux%cosgen(jst, 1, 0.5_wp, 0.0_wp, tcos)

            if (l <= jsp) then
                do j = l, jsp, l
                    jm1 = j - jsh
                    jp1 = j + jsh
                    jm2 = j - jst
                    jp2 = j + jst
                    jm3 = jm2 - jsh
                    jp3 = jp2 + jsh

                    if (jst == 1) then
                        b(:m) = 2.0_wp*q(:m, j)
                        q(:m, j) = q(:m, jm2) + q(:m, jp2)
                    else
                        do i = 1, m
                            t = q(i, j) - q(i, jm1) - q(i, jp1) + q(i, jm2) + q(i, jp2)
                            b(i) = t + q(i, j) - q(i, jm3) - q(i, jp3)
                            q(i, j) = t
                        end do
                    end if

                    call genbun_aux%trix(jst, 0, m, ba, bb, bc, b, tcos, d, w)
                    q(:m, j) = q(:m, j) + b(:m)

                end do
            end if
            !
            !==> reduction for last unknown
            !
            case_block: block
                select case (nodd)
                    !
                    !==> even number of unknowns
                    !
                    case (2)
                        jsp = jsp + l
                        j = jsp
                        jm1 = j - jsh
                        jp1 = j + jsh
                        jm2 = j - jst
                        jp2 = j + jst
                        jm3 = jm2 - jsh

                        select case (irreg)
                            case (2)
                                call genbun_aux%cosgen(kr, jstsav, 0.0, fi, tcos)
                                call genbun_aux%cosgen(lr, jstsav, 0.0, fi, tcos(kr+1))
                                ideg = kr
                                kr = kr + jst
                            case default
                                jstsav = jst
                                ideg = jst
                                kr = l
                        end select

                        if (jst == 1) then
                            irreg = 2
                            b(:m) = q(:m, j)
                            q(:m, j) = q(:m, jm2)
                        else
                            b(:m) = q(:m, j) &
                                + 0.5_wp*(q(:m, jm2)-q(:m, jm1)-q(:m, jm3))
                            select case (irreg)
                                case (2)
                                    select case (noddpr)
                                        case (2)
                                            q(:m, j) = q(:m, jm2) + q(:m, j) - q(:m, jm1)
                                        case default
                                            q(:m, j) = q(:m, jm2) + p(ipp+1:m+ipp)
                                            ipp = ipp - m
                                    end select
                                case default
                                    q(:m, j) = q(:m, jm2) + 0.5*(q(:m, j)-q(:m, jm1)-q(:m, jp1))
                                    irreg = 2
                            end select
                        end if

                        call genbun_aux%trix(ideg, lr, m, ba, bb, bc, b, tcos, d, w)
                        q(:m, j) = q(:m, j) + b(:m)

                    case default
                        !
                        !==> odd number of unknowns
                        !
                        if (irreg == 1) exit case_block

                        jsp = jsp + l
                        j = jsp
                        jm1 = j - jsh
                        jp1 = j + jsh
                        jm2 = j - jst
                        jp2 = j + jst
                        jm3 = jm2 - jsh

                        select case (istag)
                            case (1)
                                select case (noddpr)
                                    case (2)
                                        b(:m) = 0.5_wp*(q(:m, jm2)-q(:m, jm1)-q(:m, jm3)) &
                                            + q(:m, jp2) - q(:m, jp1) + q(:m, j)
                                    case default
                                        b(:m) = 0.5_wp*(q(:m, jm2)-q(:m, jm1)&
                                            -q(:m, jm3)) + p(ipp+1:m+ipp) + q(:m, j)
                                end select

                                q(:m, j) = 0.5*(q(:m, j)-q(:m, jm1)-q(:m, jp1))
                            case (2)
                                if (jst /= 1) then
                                    select case (noddpr)
                                        case (2)
                                            b(:m) = 0.5_wp*(q(:m, jm2)-q(:m, jm1)-q(:m, jm3)) &
                                                + q(:m, jp2) - q(:m, jp1) + q(:m, j)
                                        case default
                                            b(:m) = 0.5_wp*(q(:m, jm2)-q(:m, jm1)&
                                                -q(:m, jm3)) + p(ipp+1:m+ipp) + q(:m, j)
                                    end select
                                    q(:m, j) = 0.5*(q(:m, j)-q(:m, jm1)-q(:m, jp1))
                                else

                                    b(:m) = q(:m, j)
                                    q(:m, j) = 0.0_wp
                                end if
                        end select

                        call genbun_aux%trix(jst, 0, m, ba, bb, bc, b, tcos, d, w)

                        ipp = ipp + m
                        ipstor = max(ipstor, ipp + m)
                        p(ipp+1:m+ipp) = q(:m, j) + b(:m)
                        b(:m) = q(:m, jp2) + p(ipp+1:m+ipp)

                        if (lr == 0) then
                            do i = 1, jst
                                krpi = kr + i
                                tcos(krpi) = tcos(i)
                            end do
                        else
                            call genbun_aux%cosgen(lr, jstsav, 0.0_wp, fi, tcos(jst+1))
                            call genbun_aux%merger(tcos, 0, jst, jst, lr, kr)
                        end if

                        call genbun_aux%cosgen(kr, jstsav, 0.0_wp, fi, tcos)
                        call genbun_aux%trix(kr, kr, m, ba, bb, bc, b, tcos, d, w)

                        q(:m, j) = q(:m, jm2) + b(:m) + p(ipp+1:m+ipp)
                        lr = kr
                        kr = kr + l

                end select
            end block case_block

            nun = nun/2
            noddpr = nodd
            jsh = jst
            jst = 2*jst

            if (nun < 2) exit loop_108
        end do loop_108
        !
        !==> start solution.
        !
        j = jsp
        b(:m) = q(:m, j)

        select case (irreg)
            case (2)
                kr = lr + jst
                call genbun_aux%cosgen(kr, jstsav, 0.0_wp, fi, tcos)
                call genbun_aux%cosgen(lr, jstsav, 0.0_wp, fi, tcos(kr+1))
                ideg = kr
            case default
                call genbun_aux%cosgen(jst, 1, 0.5_wp, 0.0_wp, tcos)
                ideg = jst
        end select

156 continue

    call genbun_aux%trix(ideg, lr, m, ba, bb, bc, b, tcos, d, w)
    jm1 = j - jsh
    jp1 = j + jsh

    select case (irreg)
        case (2)
            select case (noddpr)
                case (2)
                    q(:m, j) = q(:m, j) - q(:m, jm1) + b(:m)
                case default
                    q(:m, j) = p(ipp+1:m+ipp) + b(:m)
                    ipp = ipp - m
            end select
        case default
            q(:m, j) = 0.5*(q(:m, j)-q(:m, jm1)-q(:m, jp1)) + b(:m)
    end select

164 continue

    jst = jst/2
    jsh = jst/2
    nun = 2*nun

    if (nun > n) then
        w(1) = ipstor
        return
    end if

    do j = jst, n, l
        jm1 = j - jsh
        jp1 = j + jsh
        jm2 = j - jst
        jp2 = j + jst

        if (j <= jst) then
            b(:m) = q(:m, j) + q(:m, jp2)
        else
            if (jp2 <= n) then
                goto 168
            end if
            b(:m) = q(:m, j) + q(:m, jm2)

            if (jst < jstsav) then
                irreg = 1
            end if

            select case (irreg)
                case (1)
                    goto 170
                case (2)
                    goto 171
            end select
168     continue
        b(:m) = q(:m, j) + q(:m, jm2) + q(:m, jp2)
    end if

170 continue

    call genbun_aux%cosgen(jst, 1, 0.5_wp, 0.0_wp, tcos)
    ideg = jst
    jdeg = 0
    goto 172

171 continue

    if (j + l > n) then
        lr = lr - jst
    end if

    kr = jst + lr
    call genbun_aux%cosgen(kr, jstsav, 0.0_wp, fi, tcos)
    call genbun_aux%cosgen(lr, jstsav, 0.0_wp, fi, tcos(kr+1))
    ideg = kr
    jdeg = lr

172 continue

    call genbun_aux%trix(ideg, jdeg, m, ba, bb, bc, b, tcos, d, w)

    if (jst <= 1) then
        q(:m, j) = b(:m)
    else
        if (jp2 > n) then
            goto 177
        end if

175 continue

    q(:m, j) = 0.5*(q(:m, j)-q(:m, jm1)-q(:m, jp1)) + b(:m)
    cycle

177 continue

    select case (irreg)
        case (1)
            goto 175
        case (2)
            if (j + jsh <= n) then
                q(:m, j) = b(:m) + p(ipp+1:m+ipp)
                ipp = ipp - m
            else
                q(:m, j) = b(:m) + q(:m, j) - q(:m, jm1)
            end if
    end select

end if
end do

l = l/2

goto 164

w(1) = ipstor

end subroutine poisd2




subroutine poisn2(m, n, istag, mixbnd, a, bb, c, q, idimq, b, b2, &
    b3, w, w2, w3, d, tcos, p)
    !
    ! Purpose
    !
    !     To solve poisson's equation with neumann boundary
    !     conditions.
    !
    !     istag = 1 if the last diagonal block is a.
    !     istag = 2 if the last diagonal block is a-i.
    !     mixbnd = 1 if have neumann boundary conditions at both boundaries.
    !     mixbnd = 2 if have neumann boundary conditions at bottom and
    !     dirichlet condition at top.  (for this case, must have istag = 1.)
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: m
    integer (ip), intent (in)     :: n
    integer (ip), intent (in)     :: istag
    integer (ip), intent (in)     :: mixbnd
    integer (ip), intent (in)     :: idimq
    real (wp),    intent (in)     :: a(*)
    real (wp),    intent (in)     :: bb(*)
    real (wp),    intent (in)     :: c(*)
    real (wp),    intent (in out) :: q(idimq, *)
    real (wp),    intent (in out) :: b(*)
    real (wp),    intent (in out) :: b2(*)
    real (wp),    intent (in out) :: b3(*)
    real (wp),    intent (in out) :: w(*)
    real (wp),    intent (in out) :: w2(*)
    real (wp),    intent (in out) :: w3(*)
    real (wp),    intent (in out) :: d(*)
    real (wp),    intent (in out) :: tcos(*)
    real (wp),    intent (in out) :: p(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip)     :: k(4)
    integer (ip)     :: mr, ipp, ipstor, i2r, jr, nr, nlast, kr
    integer (ip)     :: lr, i, nrod, jstart, jstop, i2rby2, j, jp1, jp2, jp3, jm1
    integer (ip)     :: jm2, jm3, nrodpr, ii, i1, i2, jr2, nlastp, jstep
    real (wp)        :: fistag, fnum, fden, fi, t
    type (GenbunAux) :: genbun_aux
    !-----------------------------------------------

    associate( &
        k1 => k(1), &
        k2 => k(2), &
        k3 => k(3), &
        k4 => k(4) &
        )

        fistag = 3 - istag
        fnum = 1.0_wp/istag
        fden = 0.5_wp * real(istag - 1, kind=wp)
        mr = m
        ipp = -mr
        ipstor = 0
        i2r = 1
        jr = 2
        nr = n
        nlast = n
        kr = 1
        lr = 0

        select case (istag)
            case (1)
                goto 101
            case (2)
                goto 103
        end select

101 continue

    q(:mr, n) = 0.5_wp * q(:mr, n)

    select case (mixbnd)
        case (1)
            goto 103
        case (2)
            goto 104
    end select

103 continue

    if (n <= 3) then
        goto 155
    end if

104 continue

    jr = 2*i2r

    if ((nr/2)*2 == nr) then
        nrod = 0
    else
        nrod = 1
    end if

    select case (mixbnd)
        case default
            jstart = 1
        case (2)
            jstart = jr
            nrod = 1 - nrod
    end select

    jstop = nlast - jr

    if (nrod == 0) then
        jstop = jstop - i2r
    end if

    call genbun_aux%cosgen(i2r, 1, 0.5_wp, 0.0_wp, tcos)

    i2rby2 = i2r/2

    if (jstop < jstart) then
        j = jr
    else
        do j = jstart, jstop, jr
            jp1 = j + i2rby2
            jp2 = j + i2r
            jp3 = jp2 + i2rby2
            jm1 = j - i2rby2
            jm2 = j - i2r
            jm3 = jm2 - i2rby2

            if (j == 1) then
                jm1 = jp1
                jm2 = jp2
                jm3 = jp3
            end if

            if (i2r == 1) then
                if (j == 1) jm2 = jp2
                b(:mr) = 2.0_wp*q(:mr, j)
                q(:mr, j) = q(:mr, jm2) + q(:mr, jp2)
            else
                do i = 1, mr
                    fi = q(i, j)
                    q(i, j)=q(i, j)-q(i, jm1)-q(i, jp1)+q(i, jm2)+q(i, jp2)
                    b(i) = fi + q(i, j) - q(i, jm3) - q(i, jp3)
                end do
            end if

            call genbun_aux%trix(i2r, 0, mr, a, bb, c, b, tcos, d, w)

            q(:mr, j) = q(:mr, j) + b(:mr)
            !
            !==> end of reduction for regular unknowns.
            !
        end do
        !
        !==> begin special reduction for last unknown.
        !
        j = jstop + jr
    end if

    nlast = j
    jm1 = j - i2rby2
    jm2 = j - i2r
    jm3 = jm2 - i2rby2

    if (nrod /= 0) then
        !
        !     odd number of unknowns
        !
        if (i2r == 1) then
            b(:mr) = fistag*q(:mr, j)
            q(:mr, j) = q(:mr, jm2)
        else
            b(:mr) = q(:mr, j) + 0.5_wp*(q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))
            if (nrodpr == 0) then
                q(:mr, j) = q(:mr, jm2) + p(ipp+1:mr+ipp)
                ipp = ipp - mr
            else
                q(:mr, j) = q(:mr, j) - q(:mr, jm1) + q(:mr, jm2)
            end if
            if (lr /= 0) then
                call genbun_aux%cosgen(lr, 1, 0.5_wp, fden, tcos(kr+1))
            else
                b(:mr) = fistag*b(:mr)
            end if
        end if
        call genbun_aux%cosgen(kr, 1, 0.5_wp, fden, tcos)
        call genbun_aux%trix(kr, lr, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
        kr = kr + i2r
    else
        jp1 = j + i2rby2
        jp2 = j + i2r
        if (i2r == 1) then
            b(:mr) = q(:mr, j)
            call genbun_aux%trix(1, 0, mr, a, bb, c, b, tcos, d, w)
            ipp = 0
            ipstor = mr
            select case (istag)
                case default
                    p(:mr) = b(:mr)
                    b(:mr) = b(:mr) + q(:mr, n)
                    tcos(1) = 1.0_wp
                    tcos(2) = 0.0_wp
                    call genbun_aux%trix(1, 1, mr, a, bb, c, b, tcos, d, w)
                    q(:mr, j) = q(:mr, jm2) + p(:mr) + b(:mr)
                    goto 150
                case (1)
                    p(:mr) = b(:mr)
                    q(:mr, j) = q(:mr, jm2) + 2.0_wp*q(:mr, jp2) + 3.0_wp*b(:mr)
                    goto 150
            end select
        end if

        b(:mr) = q(:mr, j) + 0.5*(q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))

        if (nrodpr == 0) then
            b(:mr) = b(:mr) + p(ipp+1:mr+ipp)
        else
            b(:mr) = b(:mr) + q(:mr, jp2) - q(:mr, jp1)
        end if

        call genbun_aux%trix(i2r, 0, mr, a, bb, c, b, tcos, d, w)
        ipp = ipp + mr
        ipstor = max(ipstor, ipp + mr)
        p(ipp+1:mr+ipp) = b(:mr) + 0.5_wp*(q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
        b(:mr) = p(ipp+1:mr+ipp) + q(:mr, jp2)

        if (lr /= 0) then
            call genbun_aux%cosgen(lr, 1, 0.5_wp, fden, tcos(i2r+1))
            call genbun_aux%merger(tcos, 0, i2r, i2r, lr, kr)
        else
            do i = 1, i2r
                ii = kr + i
                tcos(ii) = tcos(i)
            end do
        end if

        call genbun_aux%cosgen(kr, 1, 0.5_wp, fden, tcos)

        if (lr == 0) then
            select case (istag)
                case (1)
                    goto 146
                case (2)
                    goto 145
            end select
        end if

145 continue

    call genbun_aux%trix(kr, kr, mr, a, bb, c, b, tcos, d, w)
    goto 148

146 continue

    b(:mr) = fistag*b(:mr)

148 continue

    q(:mr, j) = q(:mr, jm2) + p(ipp+1:mr+ipp) + b(:mr)

150 continue
    lr = kr
    kr = kr + jr
end if

select case (mixbnd)
    case default

        nr = (nlast - 1)/jr + 1

        if (nr <= 3) then
            goto 155
        end if
    case (2)

        nr = nlast/jr

        if (nr <= 1) then
            goto 192
        end if
end select

i2r = jr
nrodpr = nrod

goto 104

155 continue

    j = 1 + jr
    jm1 = j - i2r
    jp1 = j + i2r
    jm2 = nlast - i2r

    if (nr /= 2) then
        if (lr /= 0) then
            goto 170
        end if

        if (n == 3) then
            !
            !     case n = 3.
            !
            select case (istag)
                case (1)
                    goto 156
                case (2)
                    goto 168
            end select

156     continue

        b(:mr) = q(:mr, 2)
        tcos(1) = 0.0_wp

        call genbun_aux%trix(1, 0, mr, a, bb, c, b, tcos, d, w)

        q(:mr, 2) = b(:mr)
        b(:mr) = 4.0_wp*b(:mr) + q(:mr, 1) + 2.0_wp*q(:mr, 3)
        tcos(1) = -2.0_wp
        tcos(2) = 2.0_wp
        i1 = 2
        i2 = 0

        call genbun_aux%trix(i1, i2, mr, a, bb, c, b, tcos, d, w)

        q(:mr, 2) = q(:mr, 2) + b(:mr)
        b(:mr) = q(:mr, 1) + 2.0_wp*q(:mr, 2)
        tcos(1) = 0.0_wp

        call genbun_aux%trix(1, 0, mr, a, bb, c, b, tcos, d, w)

        q(:mr, 1) = b(:mr)
        jr = 1
        i2r = 0
        goto 194
    end if
    !
    !==> case n = 2**p+1
    !
    select case (istag)
        case (1)
            goto 162
        case (2)
            goto 170
    end select

162 continue

    b(:mr) = &
        q(:mr, j) + 0.5_wp*q(:mr, 1) &
        - q(:mr, jm1) + q(:mr, nlast) - &
        q(:mr, jm2)

    call genbun_aux%cosgen(jr, 1, 0.5_wp, 0.0_wp, tcos)

    call genbun_aux%trix(jr, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, j) = 0.5_wp*(q(:mr, j)-q(:mr, jm1)-q(:mr, jp1)) + b(:mr)

    b(:mr) = q(:mr, 1) + 2.0_wp*q(:mr, nlast) + 4.0_wp*q(:mr, j)

    jr2 = 2*jr

    call genbun_aux%cosgen(jr, 1, 0.0_wp, 0.0_wp, tcos)

    tcos(jr+1:jr*2) = -tcos(jr:1:(-1))

    call genbun_aux%trix(jr2, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, j) = q(:mr, j) + b(:mr)
    b(:mr) = q(:mr, 1) + 2.0_wp*q(:mr, j)

    call genbun_aux%cosgen(jr, 1, 0.5_wp, 0.0_wp, tcos)

    call genbun_aux%trix(jr, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, 1) = 0.5_wp*q(:mr, 1) - q(:mr, jm1) + b(:mr)

    goto 194
!
!==> case of general n with nr = 3 .
!
168 continue

    b(:mr) = q(:mr, 2)
    q(:mr, 2) = 0.0_wp
    b2(:mr) = q(:mr, 3)
    b3(:mr) = q(:mr, 1)
    jr = 1
    i2r = 0
    j = 2

    goto 177

170 continue

    b(:mr) = 0.5_wp*q(:mr, 1) - q(:mr, jm1) + q(:mr, j)

    if (nrod == 0) then
        b(:mr) = b(:mr) + p(ipp+1:mr+ipp)
    else
        b(:mr) = b(:mr) + q(:mr, nlast) - q(:mr, jm2)
    end if

    do i = 1, mr
        t = 0.5_wp*(q(i, j)-q(i, jm1)-q(i, jp1))
        q(i, j) = t
        b2(i) = q(i, nlast) + t
        b3(i) = q(i, 1) + 2.0_wp*t
    end do

177 continue

    k1 = kr + 2*jr - 1
    k2 = kr + jr
    tcos(k1+1) = -2.
    k4 = k1 + 3 - istag

    call genbun_aux%cosgen(k2 + istag - 2, 1, 0.0_wp, fnum, tcos(k4))

    k4 = k1 + k2 + 1

    call genbun_aux%cosgen(jr - 1, 1, 0.0_wp, 1.0_wp, tcos(k4))
    call genbun_aux%merger(tcos, k1, k2, k1 + k2, jr - 1, 0)

    k3 = k1 + k2 + lr

    call genbun_aux%cosgen(jr, 1, 0.5, 0.0, tcos(k3+1))

    k4 = k3 + jr + 1

    call genbun_aux%cosgen(kr, 1, 0.5, fden, tcos(k4))
    call genbun_aux%merger(tcos, k3, jr, k3 + jr, kr, k1)

    if (lr /= 0) then
        call genbun_aux%cosgen(lr, 1, 0.5, fden, tcos(k4))
        call genbun_aux%merger(tcos, k3, jr, k3 + jr, lr, k3 - lr)
        call genbun_aux%cosgen(kr, 1, 0.5, fden, tcos(k4))
    end if

    k3 = kr
    k4 = kr

    call genbun_aux%tri3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)

    b(:mr) = b(:mr) + b2(:mr) + b3(:mr)
    tcos(1) = 2.0_wp

    call genbun_aux%trix(1, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, j) = q(:mr, j) + b(:mr)
    b(:mr) = q(:mr, 1) + 2.0_wp*q(:mr, j)

    call genbun_aux%cosgen(jr, 1, 0.5_wp, 0.0_wp, tcos)
    call genbun_aux%trix(jr, 0, mr, a, bb, c, b, tcos, d, w)

    if (jr == 1) then
        q(:mr, 1) = b(:mr)
        goto 194
    end if

    q(:mr, 1) = 0.5_wp*q(:mr, 1) - q(:mr, jm1) + b(:mr)
    goto 194
end if

if (n == 2) then
    !
    !==> case  n = 2
    !
    b(:mr) = q(:mr, 1)
    tcos(1) = 0.0_wp

    call genbun_aux%trix(1, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, 1) = b(:mr)
    b(:mr) = 2.*(q(:mr, 2)+b(:mr))*fistag
    tcos(1) = -fistag
    tcos(2) = 2.0_wp

    call genbun_aux%trix(2, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, 1) = q(:mr, 1) + b(:mr)
    jr = 1
    i2r = 0
    goto 194
end if

b3(:mr) = 0.0_wp
b(:mr) = q(:mr, 1) + 2.*p(ipp+1:mr+ipp)
q(:mr, 1) = 0.5*q(:mr, 1) - q(:mr, jm1)
b2(:mr) = 2.*(q(:mr, 1)+q(:mr, nlast))
k1 = kr + jr - 1
tcos(k1+1) = -2.
k4 = k1 + 3 - istag

call genbun_aux%cosgen(kr + istag - 2, 1, 0.0, fnum, tcos(k4))

k4 = k1 + kr + 1

call genbun_aux%cosgen(jr - 1, 1, 0.0, 1.0, tcos(k4))
call genbun_aux%merger(tcos, k1, kr, k1 + kr, jr - 1, 0)
call genbun_aux%cosgen(kr, 1, 0.5, fden, tcos(k1+1))

k2 = kr
k4 = k1 + k2 + 1

call genbun_aux%cosgen(lr, 1, 0.5, fden, tcos(k4))

k3 = lr
k4 = 0

call genbun_aux%tri3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)

b(:mr) = b(:mr) + b2(:mr)
tcos(1) = 2.0_wp

call genbun_aux%trix(1, 0, mr, a, bb, c, b, tcos, d, w)

q(:mr, 1) = q(:mr, 1) + b(:mr)
goto 194

192 continue

    b(:mr) = q(:mr, nlast)
    goto 196

194 continue

    j = nlast - jr
    b(:mr) = q(:mr, nlast) + q(:mr, j)

196 continue

    jm2 = nlast - i2r

    if (jr == 1) then
        q(:mr, nlast) = 0.0_wp
    else
        if (nrod == 0) then
            q(:mr, nlast) = p(ipp+1:mr+ipp)
            ipp = ipp - mr
        else
            q(:mr, nlast) = q(:mr, nlast) - q(:mr, jm2)
        end if
    end if
    call genbun_aux%cosgen(kr, 1, 0.5, fden, tcos)
    call genbun_aux%cosgen(lr, 1, 0.5, fden, tcos(kr+1))

    if (lr == 0) then
        b(:mr) = fistag*b(:mr)
    end if

    call genbun_aux%trix(kr, lr, mr, a, bb, c, b, tcos, d, w)

    q(:mr, nlast) = q(:mr, nlast) + b(:mr)
    nlastp = nlast

206 continue

    jstep = jr
    jr = i2r
    i2r = i2r/2

    if (jr == 0) then
        w(1) = ipstor
        return
    end if

    select case (mixbnd)
        case (2)
            jstart = jr
        case default
            jstart = 1 + jr
    end select

    kr = kr - jr

    if (nlast + jr <= n) then
        kr = kr - jr
        nlast = nlast + jr
        jstop = nlast - jstep
    else
        jstop = nlast - jr
    end if

    lr = kr - jr
    call genbun_aux%cosgen(jr, 1, 0.5_wp, 0.0_wp, tcos)

    do j = jstart, jstop, jstep
        jm2 = j - jr
        jp2 = j + jr

        if (j == jr) then
            b(:mr) = q(:mr, j) + q(:mr, jp2)
        else
            b(:mr) = q(:mr, j) + q(:mr, jm2) + q(:mr, jp2)
        end if

        if (jr == 1) then
            q(:mr, j) = 0.0_wp
        else
            jm1 = j - i2r
            jp1 = j + i2r
            q(:mr, j) = 0.5_wp*(q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
        end if

        call genbun_aux%trix(jr, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
    end do

    nrod = 1
    if (nlast + i2r <= n) then
        nrod = 0
    end if

    if (nlastp /= nlast) then
        goto 194
    end if

    goto 206

    w(1) = ipstor

end associate

end subroutine poisn2


subroutine poisp2(m, n, a, bb, c, q, idimq, b, b2, b3, w, w2, w3, d, tcos, p)
    !
    ! Purpose:
    !
    !     subroutine to solve poisson equation with periodic boundary
    !     conditions.
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: m
    integer (ip), intent (in)     :: n
    integer (ip), intent (in)     :: idimq
    real (wp),    intent (in)     :: a(*)
    real (wp),    intent (in)     :: bb(*)
    real (wp),    intent (in)     :: c(*)
    real (wp),    intent (in out) :: q(idimq,*)
    real (wp),    intent (in out) :: b(*)
    real (wp),    intent (in out) :: b2(*)
    real (wp),    intent (in out) :: b3(*)
    real (wp),    intent (in out) :: w(*)
    real (wp),    intent (in out) :: w2(*)
    real (wp),    intent (in out) :: w3(*)
    real (wp),    intent (in out) :: d(*)
    real (wp),    intent (in out) :: tcos(*)
    real (wp),    intent (in out) :: p(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: mr, nr, nrm1, j, nrmj, nrpj, i, lh
    real (wp)    :: ipstor
    real (wp)    :: s, t
    !-----------------------------------------------

    mr = m
    nr = (n + 1)/2
    nrm1 = nr - 1

    if (2*nr == n) then
        !
        !==> even number of unknowns
        !
        do j = 1, nrm1
            nrmj = nr - j
            nrpj = nr + j
            do i = 1, mr
                s = q(i, nrmj) - q(i, nrpj)
                t = q(i, nrmj) + q(i, nrpj)
                q(i, nrmj) = s
                q(i, nrpj) = t
            end do
        end do

        q(:mr, nr) = 2.0_wp * q(:mr, nr)
        q(:mr, n) = 2.0_wp * q(:mr, n)

        call poisd2(mr, nrm1, 1, a, bb, c, q, idimq, b, w, d, tcos, p)

        ipstor = w(1)

        call poisn2(mr, nr + 1, 1, 1, a, bb, c, q(1, nr), idimq, b, b2, &
            b3, w, w2, w3, d, tcos, p)

        ipstor = max(ipstor, w(1))

        do j = 1, nrm1
            nrmj = nr - j
            nrpj = nr + j
            do i = 1, mr
                s = 0.5_wp*(q(i, nrpj)+q(i, nrmj))
                t = 0.5_wp*(q(i, nrpj)-q(i, nrmj))
                q(i, nrmj) = s
                q(i, nrpj) = t
            end do
        end do
        q(:mr, nr) = 0.5_wp*q(:mr, nr)
        q(:mr, n) = 0.5_wp*q(:mr, n)
    else
        do j = 1, nrm1
            nrpj = n + 1 - j
            do i = 1, mr
                s = q(i, j) - q(i, nrpj)
                t = q(i, j) + q(i, nrpj)
                q(i, j) = s
                q(i, nrpj) = t
            end do
        end do

        q(:mr, nr) = 2.0_wp*q(:mr, nr)
        lh = nrm1/2

        do j = 1, lh
            nrmj = nr - j
            do i = 1, mr
                s = q(i, j)
                q(i, j) = q(i, nrmj)
                q(i, nrmj) = s
            end do
        end do

        call poisd2(mr, nrm1, 2, a, bb, c, q, &
            idimq, b, w, d, tcos, p)

        ipstor = w(1)

        call poisn2(mr, nr, 2, 1, a, bb, c, q(1, nr), &
            idimq, b, b2, b3, w, w2, w3, d, tcos, p)

        ipstor = max(ipstor, w(1))

        do j = 1, nrm1
            nrpj = nr + j
            do i = 1, mr
                s = 0.5_wp*(q(i, nrpj)+q(i, j))
                t = 0.5_wp*(q(i, nrpj)-q(i, j))
                q(i, nrpj) = t
                q(i, j) = s
            end do
        end do

        q(:mr, nr) = 0.5_wp*q(:mr, nr)

        do j = 1, lh
            nrmj = nr - j
            do i = 1, mr
                s = q(i, j)
                q(i, j) = q(i, nrmj)
                q(i, nrmj) = s
            end do
        end do
    end if

    w(1) = ipstor

end subroutine poisp2


end module module_genbun
!
! Revision history
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
