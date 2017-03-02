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
module complex_linear_systems_solver

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        PI

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: cmgnbn

    ! Parameters confined to the module
    real(wp),    parameter :: ZERO = 0.0_wp
    real(wp),    parameter :: HALF = 0.5_wp
    real(wp),    parameter :: ONE = 1.0_wp
    real(wp),    parameter :: TWO = 2.0_wp
    real(wp),    parameter :: THREE = 3.0_wp
    real(wp),    parameter :: FOUR = 4.0_wp
    integer(ip), parameter :: IIWK = 11 ! Size of workspace_indices

contains
    !
    !     SUBROUTINE cmgnbn(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
    !
    !
    ! DIMENSION OF           a(m), b(m), c(m), y(idimy, n)
    ! ARGUMENTS
    !
    ! PURPOSE                The name of this package is a mnemonic for the
    !                        complex generalized Buneman algorithm.
    !                        it solves the complex linear system of equation
    !
    !                        a(i)*x(i-1, j) + b(i)*x(i, j) + c(i)*x(i+1, j)
    !                        + x(i, j-1) - TWO * x(i, j) + x(i, j+1) = y(i, j)
    !
    !                        for i = 1, 2, ..., m  and  j = 1, 2, ..., n.
    !
    !                        indices i+1 and i-1 are evaluated modulo m,
    !                        i.e., x(0, j) = x(m, j) and x(m+1, j) = x(1, j),
    !                        and x(i, 0) may equal 0, x(i, 2), or x(i, n),
    !                        and x(i, n+1) may equal 0, x(i, n-1), or x(i, 1)
    !                        depending on an input parameter.
    !
    ! USAGE                  call cmgnbn (nperod, n, mperod, m, a, b, c, idimy, y,
    !                                     ierror)
    !
    ! ARGUMENTS
    !
    ! ON INPUT               nperod
    !
    !                          indicates the values that x(i, 0) and
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
    !                          the number of unknowns in the j-direction.
    !                          n must be greater than 2.
    !
    !                        mperod
    !                          = 0 if a(1) and c(m) are not zero
    !                          = 1 if a(1) = c(m) = 0
    !
    !                        m
    !                          the number of unknowns in the i-direction.
    !                          n must be greater than 2.
    !
    !                        a, b, c
    !                          one-dimensional complex arrays of length m
    !                          that specify the coefficients in the linear
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
    !                          the row (or first) dimension of the
    !                          two-dimensional array y as it appears
    !                          in the program calling cmgnbn.
    !                          this parameter is used to specify the
    !                          variable dimension of y.
    !                          idimy must be at least m.
    !
    !                        y
    !                          a two-dimensional complex array that
    !                          specifies the values of the right side
    !                          of the linear system of equations given
    !                          above.
    !                          y must be dimensioned at least m*n.
    !
    !
    !  ON OUTPUT             y
    !
    !                          contains the solution x.
    !
    !                        ierror
    !                          an error flag which indicates invalid
    !                          input parameters  except for number
    !                          zero, a solution is not attempted.
    !
    !                          = 0  no error.
    !                          = 1  m <= 2  .
    !                          = 2  n <= 2
    !                          = 3  idimy < m
    !                          = 4  nperod < 0 or nperod > 4
    !                          = 5  mperod < 0 or mperod > 1
    !                          = 6  a(i) /= c(1) or c(i) /= c(1) or
    !                               b(i) /= b(1) for
    !                               some i=1, 2, ..., m.
    !                          = 7  a(1) /= 0 or c(m) /= 0 and
    !                                 mperod = 1
    !                          = 20 if the dynamic allocation of real and
    !                               complex workspace required for solution
    !                               fails (for example if n, m are too large
    !                               for your computer)
    !
    ! HISTORY                Written in 1979 by Roland Sweet of NCAR'S
    !                        scientific computing division. Made available
    !                        on NCAR's public libraries in January, 1980.
    !                        Revised in June 2004 by John Adams using
    !                        Fortran 90 dynamically allocated workspace.
    !
    ! ALGORITHM              The linear system is solved by a cyclic
    !                        reduction algorithm described in the
    !                        reference below.
    !
    ! REFERENCES             Sweet, R., 'A cyclic reduction algorithm for
    !                        solving block tridiagonal systems of arbitrary
    !                        dimensions, ' SIAM J. on Numer. Anal.,
    !                          14(Sept., 1977), pp. 706-720.
    !
    ! ACCURACY               this test was performed on a platform with
    !                        64 bit floating point arithmetic.
    !                        a uniform random number generator was used
    !                        to create a solution array x for the system
    !                        given in the 'purpose' description above
    !                        with
    !                          a(i) = c(i) = -0.5 * b(i) = 1, i=1, 2, ..., m
    !
    !                        and, when mperod = 1
    !
    !                          a(1) = c(m) = 0
    !                          a(m) = c(1) = 2.
    !
    !                        the solution x was substituted into the
    !                        given system  and a right side y was
    !                        computed.  using this array y, subroutine
    !                        cmgnbn was called to produce approximate
    !                        solution z.  then relative error
    !                          e = max(abs(z(i, j)-x(i, j)))/
    !                              max(abs(x(i, j)))
    !                        was computed, where the two maxima are taken
    !                        over i=1, 2, ..., m and j=1, ..., n.
    !
    !                        the value of e is given in the table
    !                        below for some typical values of m and n.
    !
    !                   m (=n)    mperod    nperod       e
    !                   ------    ------    ------     ------
    !
    !                     31        0         0        1.e-12
    !                     31        1         1        4.e-13
    !                     31        1         3        2.e-12
    !                     32        0         0        7.e-14
    !                     32        1         1        5.e-13
    !                     32        1         3        2.e-13
    !                     33        0         0        6.e-13
    !                     33        1         1        5.e-13
    !                     33        1         3        3.e-12
    !                     63        0         0        5.e-12
    !                     63        1         1        6.e-13
    !                     63        1         3        1.e-11
    !                     64        0         0        1.e-13
    !                     64        1         1        3.e-12
    !                     64        1         3        3.e-13
    !                     65        0         0        2.e-12
    !                     65        1         1        5.e-13
    !                     65        1         3        1.e-11
    !
    !
    subroutine cmgnbn(nperod, n, mperod, m, a, b, c, idimy, y, ierror)

        ! Dummy arguments
        integer(ip), intent(in)     :: nperod
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: mperod
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: idimy
        integer(ip), intent(out)    :: ierror
        complex(wp), intent(in)     :: a(:)
        complex(wp), intent(in)     :: b(:)
        complex(wp), intent(in)     :: c(:)
        complex(wp), intent(inout)  :: y(:,:)

        ! Local variables
        type(FishpackWorkspace) :: workspace

        ! Check input arguments
        call check_input_arguments(nperod, n, mperod, m, a, b, c, idimy, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Allocate memory
        call initialize_workspace(n, m, workspace)

        ! Solve system
        associate( &
            cxw => workspace%complex_workspace, &
            indx => workspace%workspace_indices &
            )
            call cmgnbn_lower_routine(nperod, n, mperod, m, a, b, c, idimy, y, cxw, indx)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine cmgnbn

    pure subroutine check_input_arguments(nperod, n, mperod, m, a, b, c, idimy, ierror)

        ! Dummy arguments
        integer(ip), intent(in)     :: nperod
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: mperod
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: idimy
        integer(ip), intent(out)    :: ierror
        complex(wp), intent(in)     :: a(:)
        complex(wp), intent(in)     :: b(:)
        complex(wp), intent(in)     :: c(:)

        if (3 > m) then
            ierror = 1
            return
        else if (3 > n) then
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
            if (any(abs(a(2:m)-c(1)) /= ZERO) .or. &
                any(abs(c(2:m)-c(1)) /= ZERO) .or. &
                any(abs(b(2:m)-b(1)) /= ZERO)) then
                ierror = 6
                return
            end if
        else if (abs(a(1)) /= ZERO .and. abs(c(m)) /= ZERO) then
            ierror = 7
            return
        else
            ierror = 0
        end if

    end subroutine check_input_arguments

    subroutine initialize_workspace(n, m, workspace)

        ! Dummy arguments
        integer(ip),              intent(in)  :: n
        integer(ip),              intent(in)  :: m
        class(FishpackWorkspace), intent(out) :: workspace

        ! Local variables
        integer(ip) :: irwk, icwk, j

        ! Compute required workspace sizes
        irwk = 0
        icwk = (10 + int(log(real(n, kind=wp))/log(TWO), kind=ip))*m + 4*n

        ! Allocate memory
        call workspace%create(irwk, icwk, IIWK)

        ! Compute workspace indices
        associate( indx => workspace%workspace_indices )
            indx(1) = m + 1
            do j = 1, IIWK - 2
                indx(j + 1) = indx(j) + m
            end do
            indx(IIWK) = indx(IIWK-1) + 4*n
        end associate

    end subroutine initialize_workspace

    subroutine cmgnbn_lower_routine(nperod, n, mperod, m, a, b, c, idimy, y, w, workspace_indices)

        ! Dummy arguments
        integer(ip), intent(in)     :: nperod
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: mperod
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: idimy
        complex(wp), intent(in)     :: a(:)
        complex(wp), intent(in)     :: b(:)
        complex(wp), intent(in)     :: c(:)
        complex(wp), intent(inout)  :: y(:,:)
        complex(wp), intent(out)    :: w(:)
        integer(ip), intent(in)     :: workspace_indices(:)

        ! Local variables
        integer(ip) :: i, k, j, mp, np, ipstor
        integer(ip) :: irev, mh, mhm1, modd, nby2, mskip
        complex(wp) :: temp_save

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

            mp = mperod + 1
            np = nperod + 1

            select case (mp)
                case (1)
                    goto 114
                case (2)
                    goto 107
            end select
107     continue
        select case (np)
            case (1)
                goto 108
            case (2)
                goto 109
            case (3)
                goto 110
            case (4)
                goto 111
            case (5)
                goto 123
        end select
108 continue
    call solve_poisson_periodic(m, n, w(iwba:), w(iwbb:), w(iwbc:), y, idimy, w, &
        w(iwb2:), w(iwb3:), w(iww1:), w(iww2:), w(iww3:), w(iwd:), w(iwtcos:), w(iwp:))
    goto 112
109 continue
    call solve_poisson_dirichlet(m, n, 1, w(iwba:), w(iwbb:), w(iwbc:), y, idimy, w, &
        w(iww1:), w(iwd:), w(iwtcos:), w(iwp:))
    goto 112
110 continue
    call solve_poisson_neumann(m, n, 1, 2, w(iwba:), w(iwbb:), w(iwbc:), y, idimy, w, &
        w(iwb2:), w(iwb3:), w(iww1:), w(iww2:), w(iww3:), w(iwd:), w(iwtcos:), w(iwp:))
    goto 112
111 continue
    call solve_poisson_neumann(m, n, 1, 1, w(iwba:), w(iwbb:), w(iwbc:), y, idimy, w, &
        w(iwb2:), w(iwb3:), w(iww1:), w(iww2:), w(iww3:), w(iwd:), w(iwtcos:), w(iwp:))
112 continue


    ipstor = int(w(iww1), kind=ip)
    irev = 2

    if (nperod == 4) goto 124

113 continue

    select case (mp)
        case (1)
            goto 127
        case (2)
            w(1) = cmplx(real(ipstor + iwp - 1, kind=wp), ZERO, kind=wp)
            return
    end select
114 continue

    mh = (m + 1)/2
    mhm1 = mh - 1

    if (mh*2 == m) then
        modd = 2
    else
        modd = 1
    end if

    do j = 1, n
        do i = 1, mhm1
            w(i) = y(mh-i, j) - y(i+mh, j)
            w(i+mh) = y(mh-i, j) + y(i+mh, j)
        end do
        w(mh) = TWO * y(mh, j)
        select case (modd)
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

    select case (modd)
        case default
            k = iwbb + mhm1 - 1
            w(k) = w(k) - w(i-1)
            w(iwbc-1) = w(iwbc-1) + w(iwbb-1)
        case (2)
            w(iwbb-1) = w(k+1)
    end select

    goto 107
!
!     reverse columns when nperod = 4
!

123 continue

    irev = 1
    nby2 = n/2

124 continue

    do j = 1, nby2
        mskip = n + 1 - j
        do i = 1, m
            temp_save = y(i, j)
            y(i, j) = y(i, mskip)
            y(i, mskip) = temp_save
        end do
    end do

    select case (irev)
        case (1)
            goto 110
        case (2)
            goto 113
    end select

127 continue

    do j = 1, n
        w(mh-1:mh-mhm1:(-1)) = HALF * (y(mh+1:mhm1+mh, j)+y(:mhm1, j))
        w(mh+1:mhm1+mh) = HALF * (y(mh+1:mhm1+mh, j)-y(:mhm1, j))
        w(mh) = HALF * y(mh, j)
        select case (modd)
            case (1)
                y(:m, j) = w(:m)
            case (2)
                w(m) = HALF * y(m, j)
                y(:m, j) = w(:m)
        end select
    end do

    w(1) = cmplx(real(ipstor + iwp - 1, kind=wp), ZERO, kind=wp)

end associate

end subroutine cmgnbn_lower_routine

! Purpose:
!
!     To solve poisson's equation for dirichlet boundary
!     conditions.
!
!     istag = 1 if the last diagonal block is the matrix a.
!     istag = 2 if the last diagonal block is the matrix a+i.
!
subroutine solve_poisson_dirichlet(mr, nr, istag, ba, bb, bc, q, idimq, b, w, d, tcos, p)

    ! Dummy arguments
    integer(ip), intent(in)     :: mr
    integer(ip), intent(in)     :: nr
    integer(ip), intent(in)     :: istag
    integer(ip), intent(in)     :: idimq
    complex(wp), intent(in)     :: ba(mr)
    complex(wp), intent(in)     :: bb(mr)
    complex(wp), intent(in)     :: bc(mr)
    complex(wp), intent(inout) :: q(idimq,nr)
    complex(wp), intent(inout) :: b(mr)
    complex(wp), intent(inout) :: w(mr)
    complex(wp), intent(inout) :: d(mr)
    complex(wp), intent(inout) :: tcos(mr)
    complex(wp), intent(inout) :: p(nr*4)

    ! Local variables
    integer(ip) :: m, n, iip, ipstor, jsh, kr, irreg, jstsav, i, lr, nun
    integer(ip) :: jst, jsp, l, nodd, j, jm1, jp1, jm2, jp2, jm3, jp3, noddpr
    integer(ip) :: krpi, ideg, jdeg
    real(wp)    :: fi
    complex(wp) :: t

    m = mr
    n = nr
    fi = ONE/istag
    iip = -m
    ipstor = 0
    jsh = 0

    select case (istag)
        case default
            kr = 0
            irreg = 1
            if (n > 1) goto 106
            tcos(1) = ZERO
        case (2)
            kr = 1
            jstsav = 1
            irreg = 2
            if (n > 1) goto 106
            tcos(1) = cmplx(-ONE, ZERO, kind=wp)
    end select

    b(:m) = q(:m, 1)
    call solve_linear_system(1, 0, m, ba, bb, bc, b, tcos, d, w)
    q(:m, 1) = b(:m)
    goto 183
106 continue
    lr = 0
    p(1:m) = ZERO
    nun = n
    jst = 1
    jsp = n
!
!     irreg = 1 when no irregularities have occurred, otherwise it is 2.
!
108 continue
    l = 2*jst
    nodd = 2 - 2*((nun + 1)/2) + nun
    !
    !     nodd = 1 when nun is odd, otherwise it is 2.
    !
    select case (nodd)
        case default
            jsp = jsp - l
        case (1)
            jsp = jsp - jst
            if (irreg /= 1) jsp = jsp - l
    end select

    call generate_cosines(jst, 1, HALF, ZERO, tcos)

    if (l <= jsp) then
        do j = l, jsp, l
            jm1 = j - jsh
            jp1 = j + jsh
            jm2 = j - jst
            jp2 = j + jst
            jm3 = jm2 - jsh
            jp3 = jp2 + jsh
            if (jst == 1) then
                b(:m) = TWO * q(:m, j)
                q(:m, j) = q(:m, jm2) + q(:m, jp2)
            else
                do i = 1, m
                    t = q(i, j) - q(i, jm1) - q(i, jp1) + q(i, jm2) + q(i, jp2)
                    b(i) = t + q(i, j) - q(i, jm3) - q(i, jp3)
                    q(i, j) = t
                end do
            end if
            call solve_linear_system(jst, 0, m, ba, bb, bc, b, tcos, d, w)
            q(:m, j) = q(:m, j) + b(:m)
        end do
    end if
    !
    !     reduction for last unknown
    !
    select case (nodd)
        case default
            select case (irreg)
                case (1)
                    goto 152
                case (2)
                    goto 120
            end select
        !
        !     odd number of unknowns
        !
120     continue
        jsp = jsp + l
        j = jsp
        jm1 = j - jsh
        jp1 = j + jsh
        jm2 = j - jst
        jp2 = j + jst
        jm3 = jm2 - jsh
        select case (istag)
            case (1)
                goto 123
            case (2)
                goto 121
        end select
121 continue
    if (jst /= 1) goto 123
    do i = 1, m
        b(i) = q(i, j)
        q(i, j) = ZERO
    end do
    goto 130
123 continue
    select case (noddpr)
        case default
            b(:m) = HALF * (q(:m, jm2)-q(:m, jm1)-q(:m, jm3)) + p(iip+1:m+iip) &
                + q(:m, j)
        case (2)
            b(:m) = HALF * (q(:m, jm2)-q(:m, jm1)-q(:m, jm3)) + q(:m, jp2) - q( &
                :m, jp1) + q(:m, j)
    end select

    q(:m, j) = HALF * (q(:m, j)-q(:m, jm1)-q(:m, jp1))
130 continue
    call solve_linear_system(jst, 0, m, ba, bb, bc, b, tcos, d, w)
    iip = iip + m
    ipstor = max(ipstor, iip + m)
    p(iip+1:m+iip) = q(:m, j) + b(:m)
    b(:m) = q(:m, jp2) + p(iip+1:m+iip)
    if (lr == 0) then
        do i = 1, jst
            krpi = kr + i
            tcos(krpi) = tcos(i)
        end do
    else
        call generate_cosines(lr, jstsav, ZERO, fi, tcos(jst+1))
        call merge_tcos(tcos, 0, jst, jst, lr, kr)
    end if
    call generate_cosines(kr, jstsav, ZERO, fi, tcos)
    call solve_linear_system(kr, kr, m, ba, bb, bc, b, tcos, d, w)
    q(:m, j) = q(:m, jm2) + b(:m) + p(iip+1:m+iip)
    lr = kr
    kr = kr + l
!
!     even number of unknowns
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
        case default
            jstsav = jst
            ideg = jst
            kr = l
        case (2)
            call generate_cosines(kr, jstsav, ZERO, fi, tcos)
            call generate_cosines(lr, jstsav, ZERO, fi, tcos(kr+1))
            ideg = kr
            kr = kr + jst
    end select

    if (jst == 1) then
        irreg = 2
        b(:m) = q(:m, j)
        q(:m, j) = q(:m, jm2)
    else
        b(:m) = q(:m, j) + HALF * (q(:m, jm2)-q(:m, jm1)-q(:m, jm3))
        select case (irreg)
            case default
                q(:m, j) = q(:m, jm2) + HALF * (q(:m, j)-q(:m, jm1)-q(:m, jp1))
                irreg = 2
            case (2)
                select case (noddpr)
                    case default
                        q(:m, j) = q(:m, jm2) + p(iip+1:m+iip)
                        iip = iip - m
                    case (2)
                        q(:m, j) = q(:m, jm2) + q(:m, j) - q(:m, jm1)
                end select
        end select
    end if

    call solve_linear_system(ideg, lr, m, ba, bb, bc, b, tcos, d, w)
    q(:m, j) = q(:m, j) + b(:m)
end select
152 continue
    nun = nun/2
    noddpr = nodd
    jsh = jst
    jst = 2*jst
    if (nun >= 2) goto 108
    !
    !     start solution.
    !
    j = jsp
    b(:m) = q(:m, j)
    select case (irreg)
        case default
            call generate_cosines(jst, 1, HALF, ZERO, tcos)
            ideg = jst
        case (2)
            kr = lr + jst
            call generate_cosines(kr, jstsav, ZERO, fi, tcos)
            call generate_cosines(lr, jstsav, ZERO, fi, tcos(kr+1))
            ideg = kr
    end select

    call solve_linear_system(ideg, lr, m, ba, bb, bc, b, tcos, d, w)
    jm1 = j - jsh
    jp1 = j + jsh
    select case (irreg)
        case default
            q(:m, j) = HALF * (q(:m, j)-q(:m, jm1)-q(:m, jp1)) + b(:m)
        case (2)
            select case (noddpr)
                case default
                    q(:m, j) = p(iip+1:m+iip) + b(:m)
                    iip = iip - m
                case (2)
                    q(:m, j) = q(:m, j) - q(:m, jm1) + b(:m)
            end select
    end select
164 continue
    jst = jst/2
    jsh = jst/2
    nun = 2*nun
    if (nun > n) goto 183
    do j = jst, n, l
        jm1 = j - jsh
        jp1 = j + jsh
        jm2 = j - jst
        jp2 = j + jst
        if (j <= jst) then
            b(:m) = q(:m, j) + q(:m, jp2)
        else
            if (jp2 <= n) goto 168
            b(:m) = q(:m, j) + q(:m, jm2)
            if (jst < jstsav) irreg = 1
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
    call generate_cosines(jst, 1, HALF, ZERO, tcos)
    ideg = jst
    jdeg = 0
    goto 172
171 continue
    if (j + l > n) lr = lr - jst
    kr = jst + lr
    call generate_cosines(kr, jstsav, ZERO, fi, tcos)
    call generate_cosines(lr, jstsav, ZERO, fi, tcos(kr+1))
    ideg = kr
    jdeg = lr
172 continue
    call solve_linear_system(ideg, jdeg, m, ba, bb, bc, b, tcos, d, w)
    if (jst <= 1) then
        q(:m, j) = b(:m)
    else
        if (jp2 > n) goto 177
175 continue
    q(:m, j) = HALF * (q(:m, j)-q(:m, jm1)-q(:m, jp1)) + b(:m)
    cycle
177 continue
    select case (irreg)
        case (1)
            goto 175
        case (2)
            goto 178
    end select
178 continue
    if (j + jsh <= n) then
        q(:m, j) = b(:m) + p(iip+1:m+iip)
        iip = iip - m
    else
        q(:m, j) = b(:m) + q(:m, j) - q(:m, jm1)
    end if
end if
end do
l = l/2
goto 164
183 continue
    w(1) = cmplx(real(ipstor, kind=wp), ZERO, kind=wp)

end subroutine solve_poisson_dirichlet

subroutine solve_poisson_neumann(m, n, istag, mixbnd, a, bb, c, q, idimq, b, b2, &
    b3, w, w2, w3, d, tcos, p)
    !
    ! Purpose:
    !
    !     subroutine to solve poisson's equation with neumann boundary
    !     conditions.
    !
    !     istag = 1 if the last diagonal block is a.
    !     istag = 2 if the last diagonal block is a-i.
    !     mixbnd = 1 if have neumann boundary conditions at both boundaries.
    !     mixbnd = 2 if have neumann boundary conditions at bottom and
    !     dirichlet condition at top.  (for this case, must have istag = 1.)
    !

    ! Dummy arguments

    integer(ip), intent(in) :: m
    integer(ip), intent(in) :: n
    integer(ip), intent(in) :: istag
    integer(ip), intent(in) :: mixbnd
    integer(ip), intent(in) :: idimq
    complex(wp) :: a(m)
    complex(wp) :: bb(m)
    complex(wp) :: c(m)
    complex(wp), intent(inout) :: q(idimq,n)
    complex(wp) :: b(m)
    complex(wp) :: b2(m)
    complex(wp) :: b3(m)
    complex(wp) :: w(m)
    complex(wp) :: w2(m)
    complex(wp) :: w3(m)
    complex(wp) :: d(m)
    complex(wp) :: tcos(m)
    complex(wp), intent(inout) :: p(n*4)

    ! Local variables

    integer(ip) :: k(4)
    integer(ip) :: mr, iip
    integer(ip) :: ipstor, i2r, jr, nr, nlast, kr
    integer(ip) :: lr, i, nrod, jstart, jstop, i2rby2, j, jp1, jp2, jp3, jm1
    integer(ip) :: jm2, jm3, nrodpr, ii, i1, i2, jr2, nlastp, jstep
    real(wp)    :: fistag, fnum, fden
    complex(wp) :: fi, t


    associate( &
        k1 => k(1), &
        k2 => k(2), &
        k3 => k(3), &
        k4 => k(4) &
        )

        fistag = 3 - istag
        fnum = ONE/istag
        fden = HALF * real(istag - 1, kind=wp)
        mr = m
        iip = -mr
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
    q(:mr, n) = HALF * q(:mr, n)
    select case (mixbnd)
        case (1)
            goto 103
        case (2)
            goto 104
    end select
103 continue
    if (n <= 3) goto 155
104 continue
    jr = 2*i2r
    nrod = 1
    if ((nr/2)*2 == nr) nrod = 0
    select case (mixbnd)
        case default
            jstart = 1
        case (2)
            jstart = jr
            nrod = 1 - nrod
    end select

    jstop = nlast - jr
    if (nrod == 0) jstop = jstop - i2r
    call generate_cosines(i2r, 1, HALF, ZERO, tcos)
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
                b(:mr) = TWO * q(:mr, j)
                q(:mr, j) = q(:mr, jm2) + q(:mr, jp2)
            else
                do i = 1, mr
                    fi = q(i, j)
                    q(i, j)=q(i, j)-q(i, jm1)-q(i, jp1)+q(i, jm2)+q(i, jp2)
                    b(i) = fi + q(i, j) - q(i, jm3) - q(i, jp3)
                end do
            end if
            call solve_linear_system(i2r, 0, mr, a, bb, c, b, tcos, d, w)
            q(:mr, j) = q(:mr, j) + b(:mr)
        !
        !     end of reduction for regular unknowns.
        !
        end do
        !
        !     begin special reduction for last unknown.
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
            b(:mr) = q(:mr, j) + HALF * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))
            if (nrodpr == 0) then
                q(:mr, j) = q(:mr, jm2) + p(iip+1:mr+iip)
                iip = iip - mr
            else
                q(:mr, j) = q(:mr, j) - q(:mr, jm1) + q(:mr, jm2)
            end if
            if (lr /= 0) then
                call generate_cosines(lr, 1, HALF, fden, tcos(kr+1))
            else
                b(:mr) = fistag*b(:mr)
            end if
        end if
        call generate_cosines(kr, 1, HALF, fden, tcos)
        call solve_linear_system(kr, lr, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
        kr = kr + i2r
    else
        jp1 = j + i2rby2
        jp2 = j + i2r
        if (i2r == 1) then
            b(:mr) = q(:mr, j)
            call solve_linear_system(1, 0, mr, a, bb, c, b, tcos, d, w)
            iip = 0
            ipstor = mr
            select case (istag)
                case default
                    p(:mr) = b(:mr)
                    b(:mr) = b(:mr) + q(:mr, n)
                    tcos(1) = cmplx(ONE, ZERO, kind=wp)
                    tcos(2) = ZERO
                    call solve_linear_system(1, 1, mr, a, bb, c, b, tcos, d, w)
                    q(:mr, j) = q(:mr, jm2) + p(:mr) + b(:mr)
                    goto 150
                case (1)
                    p(:mr) = b(:mr)
                    q(:mr, j) = q(:mr, jm2) + TWO * q(:mr, jp2) + THREE*b(:mr)
                    goto 150
            end select
        end if
        b(:mr) = q(:mr, j) + HALF * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))
        if (nrodpr == 0) then
            b(:mr) = b(:mr) + p(iip+1:mr+iip)
        else
            b(:mr) = b(:mr) + q(:mr, jp2) - q(:mr, jp1)
        end if
        call solve_linear_system(i2r, 0, mr, a, bb, c, b, tcos, d, w)
        iip = iip + mr
        ipstor = max(ipstor, iip + mr)
        p(iip+1:mr+iip) = b(:mr) + HALF * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
        b(:mr) = p(iip+1:mr+iip) + q(:mr, jp2)
        if (lr /= 0) then
            call generate_cosines(lr, 1, HALF, fden, tcos(i2r+1))
            call merge_tcos(tcos, 0, i2r, i2r, lr, kr)
        else
            do i = 1, i2r
                ii = kr + i
                tcos(ii) = tcos(i)
            end do
        end if
        call generate_cosines(kr, 1, HALF, fden, tcos)
        if (lr == 0) then
            select case (istag)
                case (1)
                    goto 146
                case (2)
                    goto 145
            end select
        end if
145 continue
    call solve_linear_system(kr, kr, mr, a, bb, c, b, tcos, d, w)
    goto 148
146 continue
    b(:mr) = fistag*b(:mr)
148 continue
    q(:mr, j) = q(:mr, jm2) + p(iip+1:mr+iip) + b(:mr)
150 continue
    lr = kr
    kr = kr + jr
end if
select case (mixbnd)
    case default
        nr = (nlast - 1)/jr + 1
        if (nr <= 3) goto 155
    case (2)
        nr = nlast/jr
        if (nr <= 1) goto 192
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
        if (lr /= 0) goto 170
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
        tcos(1) = ZERO
        call solve_linear_system(1, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 2) = b(:mr)
        b(:mr) = FOUR*b(:mr) + q(:mr, 1) + TWO * q(:mr, 3)
        tcos(1) = cmplx(-TWO, ZERO, kind=wp)
        tcos(2) = cmplx(TWO, ZERO, kind=wp)
        i1 = 2
        i2 = 0
        call solve_linear_system(i1, i2, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 2) = q(:mr, 2) + b(:mr)
        b(:mr) = q(:mr, 1) + TWO * q(:mr, 2)
        tcos(1) = ZERO
        call solve_linear_system(1, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 1) = b(:mr)
        jr = 1
        i2r = 0
        goto 194
    end if
    !
    !     case n = 2**p+1
    !
    select case (istag)
        case (1)
            goto 162
        case (2)
            goto 170
    end select
162 continue
    b(:mr) = q(:mr, j) + HALF * q(:mr, 1) - q(:mr, jm1) + q(:mr, nlast) - &
        q(:mr, jm2)
    call generate_cosines(jr, 1, HALF, ZERO, tcos)
    call solve_linear_system(jr, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = HALF * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1)) + b(:mr)
    b(:mr) = q(:mr, 1) + TWO * q(:mr, nlast) + FOUR*q(:mr, j)
    jr2 = 2*jr
    call generate_cosines(jr, 1, ZERO, ZERO, tcos)
    tcos(jr+1:jr*2) = -tcos(jr:1:(-1))
    call solve_linear_system(jr2, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = q(:mr, j) + b(:mr)
    b(:mr) = q(:mr, 1) + TWO * q(:mr, j)
    call generate_cosines(jr, 1, HALF, ZERO, tcos)
    call solve_linear_system(jr, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = HALF * q(:mr, 1) - q(:mr, jm1) + b(:mr)
    goto 194
!
!     case of general n with nr = 3 .
!
168 continue
    b(:mr) = q(:mr, 2)
    q(:mr, 2) = ZERO
    b2(:mr) = q(:mr, 3)
    b3(:mr) = q(:mr, 1)
    jr = 1
    i2r = 0
    j = 2
    goto 177
170 continue
    b(:mr) = HALF * q(:mr, 1) - q(:mr, jm1) + q(:mr, j)
    if (nrod == 0) then
        b(:mr) = b(:mr) + p(iip+1:mr+iip)
    else
        b(:mr) = b(:mr) + q(:mr, nlast) - q(:mr, jm2)
    end if
    do i = 1, mr
        t = HALF * (q(i, j)-q(i, jm1)-q(i, jp1))
        q(i, j) = t
        b2(i) = q(i, nlast) + t
        b3(i) = q(i, 1) + TWO * t
    end do
177 continue
    k1 = kr + 2*jr - 1
    k2 = kr + jr
    tcos(k1+1) = cmplx(-TWO, ZERO, kind=wp)
    k4 = k1 + 3 - istag
    call generate_cosines(k2 + istag - 2, 1, ZERO, fnum, tcos(k4))
    k4 = k1 + k2 + 1
    call generate_cosines(jr - 1, 1, ZERO, ONE, tcos(k4))
    call merge_tcos(tcos, k1, k2, k1 + k2, jr - 1, 0)
    k3 = k1 + k2 + lr
    call generate_cosines(jr, 1, HALF, ZERO, tcos(k3+1))
    k4 = k3 + jr + 1
    call generate_cosines(kr, 1, HALF, fden, tcos(k4))
    call merge_tcos(tcos, k3, jr, k3 + jr, kr, k1)
    if (lr /= 0) then
        call generate_cosines(lr, 1, HALF, fden, tcos(k4))
        call merge_tcos(tcos, k3, jr, k3 + jr, lr, k3 - lr)
        call generate_cosines(kr, 1, HALF, fden, tcos(k4))
    end if
    k3 = kr
    k4 = kr
    call solve_tridiagonal_system(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
    b(:mr) = b(:mr) + b2(:mr) + b3(:mr)
    tcos(1) = cmplx(TWO, ZERO, kind=wp)
    call solve_linear_system(1, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = q(:mr, j) + b(:mr)
    b(:mr) = q(:mr, 1) + TWO * q(:mr, j)
    call generate_cosines(jr, 1, HALF, ZERO, tcos)
    call solve_linear_system(jr, 0, mr, a, bb, c, b, tcos, d, w)
    if (jr == 1) then
        q(:mr, 1) = b(:mr)
        goto 194
    end if
    q(:mr, 1) = HALF * q(:mr, 1) - q(:mr, jm1) + b(:mr)
    goto 194
end if
if (n == 2) then
    !
    !     case  n = 2
    !
    b(:mr) = q(:mr, 1)
    tcos(1) = ZERO
    call solve_linear_system(1, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = b(:mr)
    b(:mr) = TWO * (q(:mr, 2)+b(:mr))*fistag
    tcos(1) = cmplx((-fistag), ZERO, kind=wp)
    tcos(2) = cmplx(TWO, ZERO, kind=wp)
    call solve_linear_system(2, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = q(:mr, 1) + b(:mr)
    jr = 1
    i2r = 0
    goto 194
end if
b3(:mr) = ZERO
b(:mr) = q(:mr, 1) + TWO * p(iip+1:mr+iip)
q(:mr, 1) = HALF * q(:mr, 1) - q(:mr, jm1)
b2(:mr) = TWO * (q(:mr, 1)+q(:mr, nlast))
k1 = kr + jr - 1
tcos(k1+1) = cmplx(-TWO, ZERO, kind=wp)
k4 = k1 + 3 - istag
call generate_cosines(kr + istag - 2, 1, ZERO, fnum, tcos(k4))
k4 = k1 + kr + 1
call generate_cosines(jr - 1, 1, ZERO, ONE, tcos(k4))
call merge_tcos(tcos, k1, kr, k1 + kr, jr - 1, 0)
call generate_cosines(kr, 1, HALF, fden, tcos(k1+1))
k2 = kr
k4 = k1 + k2 + 1
call generate_cosines(lr, 1, HALF, fden, tcos(k4))
k3 = lr
k4 = 0
call solve_tridiagonal_system(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
b(:mr) = b(:mr) + b2(:mr)
tcos(1) = cmplx(TWO, ZERO, kind=wp)
call solve_linear_system(1, 0, mr, a, bb, c, b, tcos, d, w)
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
        q(:mr, nlast) = ZERO
    else
        if (nrod == 0) then
            q(:mr, nlast) = p(iip+1:mr+iip)
            iip = iip - mr
        else
            q(:mr, nlast) = q(:mr, nlast) - q(:mr, jm2)
        end if
    end if
    call generate_cosines(kr, 1, HALF, fden, tcos)
    call generate_cosines(lr, 1, HALF, fden, tcos(kr+1))
    if (lr == 0) then
        b(:mr) = fistag*b(:mr)
    end if
    call solve_linear_system(kr, lr, mr, a, bb, c, b, tcos, d, w)
    q(:mr, nlast) = q(:mr, nlast) + b(:mr)
    nlastp = nlast
206 continue
    jstep = jr
    jr = i2r
    i2r = i2r/2
    if (jr == 0) goto 222
    select case (mixbnd)
        case default
            jstart = 1 + jr
        case (2)
            jstart = jr
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
    call generate_cosines(jr, 1, HALF, ZERO, tcos)
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
        call solve_linear_system(jr, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
    end do
    nrod = 1
    if (nlast + i2r <= n) nrod = 0
    if (nlastp /= nlast) goto 194
    goto 206
222 continue
    w(1) = cmplx(real(ipstor, kind=wp), ZERO, kind=wp)

end associate

end subroutine solve_poisson_neumann

! Purpose:
!
! To solve poisson equation with periodic boundary conditions.
!
subroutine solve_poisson_periodic(m, n, a, bb, c, q, idimq, b, b2, b3, w, w2, w3, d, tcos, p)

    ! Dummy arguments
    integer(ip), intent(in) :: m
    integer(ip), intent(in) :: n
    integer(ip), intent(in) :: idimq
    complex(wp) :: a(*)
    complex(wp) :: bb(*)
    complex(wp) :: c(*)
    complex(wp) :: q(idimq,n)
    complex(wp) :: b(*)
    complex(wp) :: b2(*)
    complex(wp) :: b3(*)
    complex(wp) :: w(*)
    complex(wp) :: w2(*)
    complex(wp) :: w3(*)
    complex(wp) :: d(*)
    complex(wp) :: tcos(*)
    complex(wp) :: p(n*4)

    ! Local variables
    integer(ip) :: mr, nr, nrm1, j, nrmj, nrpj, i, lh
    real(wp)    :: ipstor
    complex(wp) :: s, t

    mr = m
    nr = (n + 1)/2
    nrm1 = nr - 1

    if ((2*nr) == n) then
        !
        !     even number of unknowns
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

        q(:mr, nr) = TWO * q(:mr, nr)
        q(:mr, n) = TWO * q(:mr, n)

        call solve_poisson_dirichlet (mr, nrm1, 1, a, bb, c, q, idimq, b, w, d, tcos, p)

        ipstor = real(w(1), kind=wp)

        call solve_poisson_neumann(mr, nr + 1, 1, 1, a, bb, c, q(1, nr), idimq, b, b2, &
            b3, w, w2, w3, d, tcos, p)

        ipstor = max(ipstor, real(w(1), kind=wp))

        do j = 1, nrm1
            nrmj = nr - j
            nrpj = nr + j
            do i = 1, mr
                s = HALF * (q(i, nrpj)+q(i, nrmj))
                t = HALF * (q(i, nrpj)-q(i, nrmj))
                q(i, nrmj) = s
                q(i, nrpj) = t
            end do
        end do

        q(:mr, nr) = HALF * q(:mr, nr)
        q(:mr, n) = HALF * q(:mr, n)

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

        q(:mr, nr) = TWO * q(:mr, nr)
        lh = nrm1/2

        do j = 1, lh
            nrmj = nr - j
            do i = 1, mr
                s = q(i, j)
                q(i, j) = q(i, nrmj)
                q(i, nrmj) = s
            end do
        end do

        call solve_poisson_dirichlet(mr, nrm1, 2, a, bb, c, q, idimq, b, w, d, tcos, p)

        ipstor = real(w(1), kind=wp)

        call solve_poisson_neumann(mr, nr, 2, 1, a, bb, c, q(1, nr), idimq, b, b2, b3, &
            w, w2, w3, d, tcos, p)

        ipstor = max(ipstor, real(w(1), kind=wp))

        do j = 1, nrm1
            nrpj = nr + j
            do i = 1, mr
                s = HALF * (q(i, nrpj)+q(i, j))
                t = HALF * (q(i, nrpj)-q(i, j))
                q(i, nrpj) = t
                q(i, j) = s
            end do
        end do

        q(:mr, nr) = HALF * q(:mr, nr)

        do j = 1, lh
            nrmj = nr - j
            do i = 1, mr
                s = q(i, j)
                q(i, j) = q(i, nrmj)
                q(i, nrmj) = s
            end do
        end do
    end if

    w(1) = cmplx(ipstor, ZERO, kind=wp)

end subroutine solve_poisson_periodic

! Purpose:
!
! Computes required cosine values in ascending
! order. When ijump > 1 the routine computes values
!
! 2*cos(j*pi/l) , j=1, 2, ..., l and j /= 0(mod n/ijump+1)
!
! where l = ijump*(n/ijump+1).
!
!
! when ijump = 1 it computes
!
! 2*cos((j-fnum)*pi/(n+fden)) ,  j=1, 2, ... , n
!
! where
!
!        fnum = 0.5, fden = 0.0, for regular reduction values
!        fnum = 0.0, fden = 1.0, for b-r and c-r when istag = 1
!        fnum = 0.0, fden = 0.5, for b-r and c-r when istag = 2
!        fnum = 0.5, fden = 0.5, for b-r and c-r when istag = 2
!                                in solve_poisson_neumann only.
!
!
pure subroutine generate_cosines(n, ijump, fnum, fden, a)

    ! Dummy arguments
    integer(ip), intent(in)  :: n
    integer(ip), intent(in)  :: ijump
    real(wp),    intent(in)  :: fnum
    real(wp),    intent(in)  :: fden
    complex(wp), intent(out) :: a(*)

    ! Local variables
    integer(ip)         :: k3, k4, k, k1, k5, i, k2, np1
    real(wp)            :: pibyn, x, y

    if (n /= 0) then
        if (ijump /= 1) then
            k3 = n/ijump + 1
            k4 = k3 - 1
            pibyn = PI/(n + ijump)
            do k = 1, ijump
                k1 = (k - 1)*k3
                k5 = (k - 1)*k4
                do i = 1, k4
                    x = k1 + i
                    k2 = k5 + i
                    a(k2) = cmplx((-TWO * cos(x*pibyn)), ZERO, kind=wp)
                end do
            end do
        else
            np1 = n + 1
            y = PI/(real(n, kind=wp) + fden)
            do i = 1, n
                x = real(np1 - i, kind=wp) - fnum
                a(i) = cmplx(TWO * cos(x*y), ZERO, kind=wp)
            end do
        end if
    end if

end subroutine generate_cosines

! Purpose:
!
! Merges two ascending strings of numbers in the
! array tcos.  The first string is of length m1 and starts at
! tcos(i1+1).  The second string is of length m2 and starts at
! tcos(i2+1).  The merged string goes into tcos(i3+1).
!
subroutine merge_tcos(tcos, i1, m1, i2, m2, i3)

    ! Dummy arguments
    integer(ip), intent(in)     :: i1
    integer(ip), intent(in)     :: m1
    integer(ip), intent(in)     :: i2
    integer(ip), intent(in)     :: m2
    integer(ip), intent(in)     :: i3
    complex(wp), intent(inout) :: tcos(*)

    ! Local variables
    integer(ip) :: j11, j3, j1, j2, j, l, k, m
    complex(wp) :: x, y

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
                        if (real(x - y, kind=wp) > ZERO) exit block_construct
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


end subroutine merge_tcos

! Purpose:
!
! To solve a system of linear equations where the
! coefficient matrix is a rational function in the matrix given by
! tridiagonal  ( . . . , a(i), b(i), c(i), . . . ).
!
subroutine solve_linear_system(idegbr, idegcr, m, a, b, c, y, tcos, d, w)

    ! Dummy arguments
    integer(ip), intent(in)     :: idegbr
    integer(ip), intent(in)     :: idegcr
    integer(ip), intent(in)     :: m
    complex(wp), intent(in)     :: a(m)
    complex(wp), intent(in)     :: b(m)
    complex(wp), intent(in)     :: c(m)
    complex(wp), intent(inout) :: y(m)
    complex(wp), intent(in)     :: tcos(*)
    complex(wp), intent(out)    :: d(m)
    complex(wp), intent(out)    :: w(m)

    ! Local variables
    integer(ip) :: mm1, ifb, ifc, l, lint, k, i, iip
    complex(wp) :: x, xx, z

    mm1 = m - 1
    ifb = idegbr + 1
    ifc = idegcr + 1
    l = ifb/ifc
    lint = 1

    do k = 1, idegbr
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

        if (abs(z) == ZERO) then
            y(m) = ZERO
        else
            y(m) = (y(m)-a(m)*y(mm1))/z
        end if

        do iip = 1, mm1
            y(m-iip) = y(m-iip) - d(m-iip)*y(m+1-iip)
        end do

        if (k /= l) cycle

        y = y + w
        lint = lint + 1
        l = (lint*ifb)/ifc
    end do

end subroutine solve_linear_system

subroutine solve_tridiagonal_system(m, a, b, c, k, y1, y2, y3, tcos, d, w1, w2, w3)

    ! Dummy arguments
    integer(ip), intent(in)     :: m
    integer(ip), intent(in)     :: k(4)
    complex(wp), intent(in)     :: a(m)
    complex(wp), intent(in)     :: b(m)
    complex(wp), intent(in)     :: c(m)
    complex(wp), intent(inout) :: y1(m)
    complex(wp), intent(inout) :: y2(m)
    complex(wp), intent(inout) :: y3(m)
    complex(wp), intent(in)     :: tcos(*)
    complex(wp), intent(out)    :: d(m)
    complex(wp), intent(out)    :: w1(m)
    complex(wp), intent(out)    :: w2(m)
    complex(wp), intent(out)    :: w3(m)

    ! Local variables
    integer(ip) :: mm1, k1, k2, k3, k4
    integer(ip) :: if1, if2, if3, if4, k2k3k4, l1, l2
    integer(ip) :: l3, lint1, lint2, lint3
    integer(ip) :: kint1, kint2, kint3, n, i, iip
    complex(wp) :: x, xx, z

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

    do n = 1, k1
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

        do iip = 1, mm1
            y1(m-iip) = y1(m-iip) - d(m-iip)*y1(m+1-iip)
            y2(m-iip) = y2(m-iip) - d(m-iip)*y2(m+1-iip)
            y3(m-iip) = y3(m-iip) - d(m-iip)*y3(m+1-iip)
        end do

        if (k2k3k4 == 0) cycle

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

        if (n /= l3) cycle

        i = lint3 + kint3
        xx = x - tcos(i)
        y3 = xx*y3 + w3
        lint3 = lint3 + 1
        l3 = (lint3*if1)/if4
    end do

end subroutine solve_tridiagonal_system

end module complex_linear_systems_solver
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
! May       2016    Fortran 2008 changes
!
