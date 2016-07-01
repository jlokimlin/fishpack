!
!     file cmgnbn.f90
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
!     SUBROUTINE cmgnbn(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
!
!
! DIMENSION OF           a(m), b(m), c(m), y(idimy, n)
! ARGUMENTS
!
! LATEST REVISION        May 2016
!
! PURPOSE                The name of this package is a mnemonic for the
!                        complex generalized Buneman algorithm.
!                        it solves the complex linear system of equation
!
!                        a(i)*x(i-1, j) + b(i)*x(i, j) + c(i)*x(i+1, j)
!                        + x(i, j-1) - 2.0_wp * x(i, j) + x(i, j+1) = y(i, j)
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
!                               complex work space required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
!
! SPECIAL CONDITONS      None
!
! I/O                    None
!
! PRECISION              64-bit double precision
!
! REQUIRED LIBRARY       comf.f90, type_FishpackWorkspace.f90
! FILES
!
! STANDARD               Fortran 2008
!
! HISTORY                Written in 1979 by Roland Sweet of NCAR'S
!                        scientific computing division. Made available
!                        on NCAR's public libraries in January, 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated work space.
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
!                          a(i) = c(i) = -0.5_wp * b(i) = 1, i=1, 2, ..., m
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
module module_cmgnbn

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        PI

    use type_FishpackWorkspace, only: &
        Fish => FishpackWorkspace

    ! Explicit typing only
    implicit None

    ! Everything is private unless stated otherwise
    private
    public :: cmgnbn


contains


    subroutine cmgnbn(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: nperod
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: mperod
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: idimy
        integer (ip), intent (out)    :: ierror
        complex (wp), intent (in)     :: a(:)
        complex (wp), intent (in)     :: b(:)
        complex (wp), intent (in)     :: c(:)
        complex (wp), intent (in out) :: y(:,:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: irwk, icwk
        type (Fish)  :: workspace
        !-----------------------------------------------

        !
        !==> Check validity of input arguments
        !
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
            if (any(abs(a(2:m)-c(1)) /= 0.0_wp) .or. &
                any(abs(c(2:m)-c(1)) /= 0.0_wp) .or. &
                any(abs(b(2:m)-b(1)) /= 0.0_wp)) then
                ierror = 6
                return
            end if
        else if (abs(a(1)) /= 0.0_wp .and. abs(c(m)) /= 0.0_wp) then
            ierror = 7
            return
        else
            ierror = 0
        end if

        !
        !==> Allocate memory
        !
        irwk = 0
        icwk = (10 + int(log(real(n, kind=wp))/log(2.0_wp), kind=ip))*m + 4*n
        call workspace%create(irwk, icwk, ierror)

        ! Check if allocation was succesful
        if (ierror == 20) return

        associate( cxw => workspace%complex_workspace )
            !
            !==> Solve system
            !
            call cmgnbnn(nperod, n, mperod, m, a, b, c, idimy, y, cxw)

        end associate

        !
        !==> Release memory
        !
        call workspace%destroy()

    end subroutine cmgnbn



    subroutine cmgnbnn(nperod, n, mperod, m, a, b, c, idimy, y, w)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: nperod
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: mperod
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: idimy
        complex (wp), intent (in)     :: a(:)
        complex (wp), intent (in)     :: b(:)
        complex (wp), intent (in)     :: c(:)
        complex (wp), intent (in out) :: y(:,:)
        complex (wp), intent (in out) :: w(:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: workspace_indices(11)
        integer (ip) :: i, k, j, mp, np, ipstor
        integer (ip) :: irev, mh, mhm1, modd, nby2, mskip
        complex (wp) :: temp_save
        !-----------------------------------------------

        workspace_indices = get_workspace_indices(n, m)

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
                w(k) = 2.0_wp - b(i)
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
    call cmposp(m, n, w(iwba:), w(iwbb:), w(iwbc:), y, idimy, w, &
        w(iwb2:), w(iwb3:), w(iww1:), w(iww2:), w(iww3:), w(iwd:), w(iwtcos:), w(iwp:))
    goto 112
109 continue
    call cmposd(m, n, 1, w(iwba:), w(iwbb:), w(iwbc:), y, idimy, w, &
        w(iww1:), w(iwd:), w(iwtcos:), w(iwp:))
    goto 112
110 continue
    call cmposn(m, n, 1, 2, w(iwba:), w(iwbb:), w(iwbc:), y, idimy, w, &
        w(iwb2:), w(iwb3:), w(iww1:), w(iww2:), w(iww3:), w(iwd:), w(iwtcos:), w(iwp:))
    goto 112
111 continue
    call cmposn(m, n, 1, 1, w(iwba:), w(iwbb:), w(iwbc:), y, idimy, w, &
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
            w(1) = cmplx(real(ipstor + iwp - 1, kind=wp), 0.0_wp, kind=wp)
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
        w(mh) = 2.0_wp * y(mh, j)
        select case (modd)
            case (1)
                y(:m, j) = w(:m)
            case (2)
                w(m) = 2.0_wp * y(m, j)
                y(:m, j) = w(:m)
        end select
    end do

    k = iwbc + mhm1 - 1
    i = iwba + mhm1
    w(k) = 0.0_wp
    w(i) = 0.0_wp
    w(k+1) = 2.0_wp * w(k+1)

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
        w(mh-1:mh-mhm1:(-1)) = 0.5_wp * (y(mh+1:mhm1+mh, j)+y(:mhm1, j))
        w(mh+1:mhm1+mh) = 0.5_wp * (y(mh+1:mhm1+mh, j)-y(:mhm1, j))
        w(mh) = 0.5_wp * y(mh, j)
        select case (modd)
            case (1)
                y(:m, j) = w(:m)
            case (2)
                w(m) = 0.5_wp * y(m, j)
                y(:m, j) = w(:m)
        end select
    end do

    w(1) = cmplx(real(ipstor + iwp - 1, kind=wp), 0.0_wp, kind=wp)

end associate

end subroutine cmgnbnn



pure function get_workspace_indices(n, m) result (return_value)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in) :: n
    integer (ip), intent (in) :: m
    integer (ip)              :: return_value(11)
    !-----------------------------------------------
    integer (ip) :: j
    !-----------------------------------------------

    associate (i => return_value)
        i(1) = m + 1
        do j = 1, 9
            i(j + 1) = i(j) + m
        end do
        i(11) = i(10) + 4*n
    end associate

end function get_workspace_indices



subroutine cmposd(mr, nr, istag, ba, bb, bc, q, idimq, b, w, d, tcos, p)
    !
    !     subroutine to solve poisson's equation for dirichlet boundary
    !     conditions.
    !
    !     istag = 1 if the last diagonal block is the matrix a.
    !     istag = 2 if the last diagonal block is the matrix a+i.
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: mr
    integer (ip), intent (in)     :: nr
    integer (ip), intent (in)     :: istag
    integer (ip), intent (in)     :: idimq
    complex (wp), intent (in)     :: ba(mr)
    complex (wp), intent (in)     :: bb(mr)
    complex (wp), intent (in)     :: bc(mr)
    complex (wp), intent (in out) :: q(idimq,nr)
    complex (wp), intent (in out) :: b(mr)
    complex (wp), intent (in out) :: w(mr)
    complex (wp), intent (in out) :: d(mr)
    complex (wp), intent (in out) :: tcos(mr)
    complex (wp), intent (in out) :: p(nr*4)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: m, n, iip, ipstor, jsh, kr, irreg, jstsav, i, lr, nun
    integer (ip) :: jst, jsp, l, nodd, j, jm1, jp1, jm2, jp2, jm3, jp3, noddpr
    integer (ip) :: krpi, ideg, jdeg
    real (wp)    :: fi
    complex (wp) :: t
    !-----------------------------------------------

    m = mr
    n = nr
    fi = 1.0_wp/istag
    iip = -m
    ipstor = 0
    jsh = 0

    select case (istag)
        case default
            kr = 0
            irreg = 1
            if (n > 1) goto 106
            tcos(1) = 0.0_wp
        case (2)
            kr = 1
            jstsav = 1
            irreg = 2
            if (n > 1) goto 106
            tcos(1) = cmplx(-1.0_wp, 0.0_wp, kind=wp)
    end select

    b(:m) = q(:m, 1)
    call cmptrx(1, 0, m, ba, bb, bc, b, tcos, d, w)
    q(:m, 1) = b(:m)
    goto 183
106 continue
    lr = 0
    p(1:m) = 0.0_wp
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

    call cmpcsg(jst, 1, 0.5_wp, 0.0_wp, tcos)

    if (l <= jsp) then
        do j = l, jsp, l
            jm1 = j - jsh
            jp1 = j + jsh
            jm2 = j - jst
            jp2 = j + jst
            jm3 = jm2 - jsh
            jp3 = jp2 + jsh
            if (jst == 1) then
                b(:m) = 2.0_wp * q(:m, j)
                q(:m, j) = q(:m, jm2) + q(:m, jp2)
            else
                do i = 1, m
                    t = q(i, j) - q(i, jm1) - q(i, jp1) + q(i, jm2) + q(i, jp2)
                    b(i) = t + q(i, j) - q(i, jm3) - q(i, jp3)
                    q(i, j) = t
                end do
            end if
            call cmptrx(jst, 0, m, ba, bb, bc, b, tcos, d, w)
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
        q(i, j) = 0.0_wp
    end do
    goto 130
123 continue
    select case (noddpr)
        case default
            b(:m) = 0.5_wp * (q(:m, jm2)-q(:m, jm1)-q(:m, jm3)) + p(iip+1:m+iip) &
                + q(:m, j)
        case (2)
            b(:m) = 0.5_wp * (q(:m, jm2)-q(:m, jm1)-q(:m, jm3)) + q(:m, jp2) - q( &
                :m, jp1) + q(:m, j)
    end select

    q(:m, j) = 0.5_wp * (q(:m, j)-q(:m, jm1)-q(:m, jp1))
130 continue
    call cmptrx(jst, 0, m, ba, bb, bc, b, tcos, d, w)
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
        call cmpcsg (lr, jstsav, 0.0_wp, fi, tcos(jst+1))
        call cmpmrg (tcos, 0, jst, jst, lr, kr)
    end if
    call cmpcsg (kr, jstsav, 0.0_wp, fi, tcos)
    call cmptrx(kr, kr, m, ba, bb, bc, b, tcos, d, w)
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
            call cmpcsg (kr, jstsav, 0.0_wp, fi, tcos)
            call cmpcsg (lr, jstsav, 0.0_wp, fi, tcos(kr+1))
            ideg = kr
            kr = kr + jst
    end select

    if (jst == 1) then
        irreg = 2
        b(:m) = q(:m, j)
        q(:m, j) = q(:m, jm2)
    else
        b(:m) = q(:m, j) + 0.5_wp * (q(:m, jm2)-q(:m, jm1)-q(:m, jm3))
        select case (irreg)
            case default
                q(:m, j) = q(:m, jm2) + 0.5_wp * (q(:m, j)-q(:m, jm1)-q(:m, jp1))
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

    call cmptrx(ideg, lr, m, ba, bb, bc, b, tcos, d, w)
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
            call cmpcsg (jst, 1, 0.5_wp, 0.0_wp, tcos)
            ideg = jst
        case (2)
            kr = lr + jst
            call cmpcsg (kr, jstsav, 0.0_wp, fi, tcos)
            call cmpcsg (lr, jstsav, 0.0_wp, fi, tcos(kr+1))
            ideg = kr
    end select

    call cmptrx(ideg, lr, m, ba, bb, bc, b, tcos, d, w)
    jm1 = j - jsh
    jp1 = j + jsh
    select case (irreg)
        case default
            q(:m, j) = 0.5_wp * (q(:m, j)-q(:m, jm1)-q(:m, jp1)) + b(:m)
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
    call cmpcsg (jst, 1, 0.5_wp, 0.0_wp, tcos)
    ideg = jst
    jdeg = 0
    goto 172
171 continue
    if (j + l > n) lr = lr - jst
    kr = jst + lr
    call cmpcsg (kr, jstsav, 0.0_wp, fi, tcos)
    call cmpcsg (lr, jstsav, 0.0_wp, fi, tcos(kr+1))
    ideg = kr
    jdeg = lr
172 continue
    call cmptrx(ideg, jdeg, m, ba, bb, bc, b, tcos, d, w)
    if (jst <= 1) then
        q(:m, j) = b(:m)
    else
        if (jp2 > n) goto 177
175 continue
    q(:m, j) = 0.5_wp * (q(:m, j)-q(:m, jm1)-q(:m, jp1)) + b(:m)
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
    w(1) = cmplx(real(ipstor, kind=wp), 0.0_wp, kind=wp)

end subroutine cmposd


subroutine cmposn(m, n, istag, mixbnd, a, bb, c, q, idimq, b, b2, &
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
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in) :: m
    integer (ip), intent (in) :: n
    integer (ip), intent (in) :: istag
    integer (ip), intent (in) :: mixbnd
    integer (ip), intent (in) :: idimq
    complex (wp) :: a(m)
    complex (wp) :: bb(m)
    complex (wp) :: c(m)
    complex (wp), intent (in out) :: q(idimq,n)
    complex (wp) :: b(m)
    complex (wp) :: b2(m)
    complex (wp) :: b3(m)
    complex (wp) :: w(m)
    complex (wp) :: w2(m)
    complex (wp) :: w3(m)
    complex (wp) :: d(m)
    complex (wp) :: tcos(m)
    complex (wp), intent (in out) :: p(n*4)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: k(4)
    integer (ip) :: mr, iip
    integer (ip) :: ipstor, i2r, jr, nr, nlast, kr
    integer (ip) :: lr, i, nrod, jstart, jstop, i2rby2, j, jp1, jp2, jp3, jm1
    integer (ip) :: jm2, jm3, nrodpr, ii, i1, i2, jr2, nlastp, jstep
    real (wp)    :: fistag, fnum, fden
    complex (wp) :: fi, t
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
    q(:mr, n) = 0.5_wp * q(:mr, n)
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
    call cmpcsg (i2r, 1, 0.5_wp, 0.0_wp, tcos)
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
                b(:mr) = 2.0_wp * q(:mr, j)
                q(:mr, j) = q(:mr, jm2) + q(:mr, jp2)
            else
                do i = 1, mr
                    fi = q(i, j)
                    q(i, j)=q(i, j)-q(i, jm1)-q(i, jp1)+q(i, jm2)+q(i, jp2)
                    b(i) = fi + q(i, j) - q(i, jm3) - q(i, jp3)
                end do
            end if
            call cmptrx(i2r, 0, mr, a, bb, c, b, tcos, d, w)
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
            b(:mr) = q(:mr, j) + 0.5_wp * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))
            if (nrodpr == 0) then
                q(:mr, j) = q(:mr, jm2) + p(iip+1:mr+iip)
                iip = iip - mr
            else
                q(:mr, j) = q(:mr, j) - q(:mr, jm1) + q(:mr, jm2)
            end if
            if (lr /= 0) then
                call cmpcsg (lr, 1, 0.5_wp, fden, tcos(kr+1))
            else
                b(:mr) = fistag*b(:mr)
            end if
        end if
        call cmpcsg (kr, 1, 0.5_wp, fden, tcos)
        call cmptrx(kr, lr, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
        kr = kr + i2r
    else
        jp1 = j + i2rby2
        jp2 = j + i2r
        if (i2r == 1) then
            b(:mr) = q(:mr, j)
            call cmptrx(1, 0, mr, a, bb, c, b, tcos, d, w)
            iip = 0
            ipstor = mr
            select case (istag)
                case default
                    p(:mr) = b(:mr)
                    b(:mr) = b(:mr) + q(:mr, n)
                    tcos(1) = cmplx(1.0_wp, 0.0_wp, kind=wp)
                    tcos(2) = 0.0_wp
                    call cmptrx(1, 1, mr, a, bb, c, b, tcos, d, w)
                    q(:mr, j) = q(:mr, jm2) + p(:mr) + b(:mr)
                    goto 150
                case (1)
                    p(:mr) = b(:mr)
                    q(:mr, j) = q(:mr, jm2) + 2.0_wp * q(:mr, jp2) + 3.0_wp*b(:mr)
                    goto 150
            end select
        end if
        b(:mr) = q(:mr, j) + 0.5_wp * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))
        if (nrodpr == 0) then
            b(:mr) = b(:mr) + p(iip+1:mr+iip)
        else
            b(:mr) = b(:mr) + q(:mr, jp2) - q(:mr, jp1)
        end if
        call cmptrx(i2r, 0, mr, a, bb, c, b, tcos, d, w)
        iip = iip + mr
        ipstor = max(ipstor, iip + mr)
        p(iip+1:mr+iip) = b(:mr) + 0.5_wp * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
        b(:mr) = p(iip+1:mr+iip) + q(:mr, jp2)
        if (lr /= 0) then
            call cmpcsg (lr, 1, 0.5_wp, fden, tcos(i2r+1))
            call cmpmrg (tcos, 0, i2r, i2r, lr, kr)
        else
            do i = 1, i2r
                ii = kr + i
                tcos(ii) = tcos(i)
            end do
        end if
        call cmpcsg (kr, 1, 0.5_wp, fden, tcos)
        if (lr == 0) then
            select case (istag)
                case (1)
                    goto 146
                case (2)
                    goto 145
            end select
        end if
145 continue
    call cmptrx(kr, kr, mr, a, bb, c, b, tcos, d, w)
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
        tcos(1) = 0.0_wp
        call cmptrx(1, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 2) = b(:mr)
        b(:mr) = 4.0_wp*b(:mr) + q(:mr, 1) + 2.0_wp * q(:mr, 3)
        tcos(1) = cmplx(-2.0_wp, 0.0_wp, kind=wp)
        tcos(2) = cmplx(2.0_wp, 0.0_wp, kind=wp)
        i1 = 2
        i2 = 0
        call cmptrx(i1, i2, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 2) = q(:mr, 2) + b(:mr)
        b(:mr) = q(:mr, 1) + 2.0_wp * q(:mr, 2)
        tcos(1) = 0.0_wp
        call cmptrx(1, 0, mr, a, bb, c, b, tcos, d, w)
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
    b(:mr) = q(:mr, j) + 0.5_wp * q(:mr, 1) - q(:mr, jm1) + q(:mr, nlast) - &
        q(:mr, jm2)
    call cmpcsg (jr, 1, 0.5_wp, 0.0_wp, tcos)
    call cmptrx(jr, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = 0.5_wp * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1)) + b(:mr)
    b(:mr) = q(:mr, 1) + 2.0_wp * q(:mr, nlast) + 4.0_wp*q(:mr, j)
    jr2 = 2*jr
    call cmpcsg (jr, 1, 0.0_wp, 0.0_wp, tcos)
    tcos(jr+1:jr*2) = -tcos(jr:1:(-1))
    call cmptrx(jr2, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = q(:mr, j) + b(:mr)
    b(:mr) = q(:mr, 1) + 2.0_wp * q(:mr, j)
    call cmpcsg (jr, 1, 0.5_wp, 0.0_wp, tcos)
    call cmptrx(jr, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = 0.5_wp * q(:mr, 1) - q(:mr, jm1) + b(:mr)
    goto 194
!
!     case of general n with nr = 3 .
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
    b(:mr) = 0.5_wp * q(:mr, 1) - q(:mr, jm1) + q(:mr, j)
    if (nrod == 0) then
        b(:mr) = b(:mr) + p(iip+1:mr+iip)
    else
        b(:mr) = b(:mr) + q(:mr, nlast) - q(:mr, jm2)
    end if
    do i = 1, mr
        t = 0.5_wp * (q(i, j)-q(i, jm1)-q(i, jp1))
        q(i, j) = t
        b2(i) = q(i, nlast) + t
        b3(i) = q(i, 1) + 2.0_wp * t
    end do
177 continue
    k1 = kr + 2*jr - 1
    k2 = kr + jr
    tcos(k1+1) = cmplx(-2.0_wp, 0.0_wp, kind=wp)
    k4 = k1 + 3 - istag
    call cmpcsg (k2 + istag - 2, 1, 0.0_wp, fnum, tcos(k4))
    k4 = k1 + k2 + 1
    call cmpcsg (jr - 1, 1, 0.0_wp, 1.0_wp, tcos(k4))
    call cmpmrg (tcos, k1, k2, k1 + k2, jr - 1, 0)
    k3 = k1 + k2 + lr
    call cmpcsg (jr, 1, 0.5_wp, 0.0_wp, tcos(k3+1))
    k4 = k3 + jr + 1
    call cmpcsg (kr, 1, 0.5_wp, fden, tcos(k4))
    call cmpmrg (tcos, k3, jr, k3 + jr, kr, k1)
    if (lr /= 0) then
        call cmpcsg (lr, 1, 0.5_wp, fden, tcos(k4))
        call cmpmrg (tcos, k3, jr, k3 + jr, lr, k3 - lr)
        call cmpcsg (kr, 1, 0.5_wp, fden, tcos(k4))
    end if
    k3 = kr
    k4 = kr
    call cmptr3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
    b(:mr) = b(:mr) + b2(:mr) + b3(:mr)
    tcos(1) = cmplx(2.0_wp, 0.0_wp, kind=wp)
    call cmptrx(1, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = q(:mr, j) + b(:mr)
    b(:mr) = q(:mr, 1) + 2.0_wp * q(:mr, j)
    call cmpcsg (jr, 1, 0.5_wp, 0.0_wp, tcos)
    call cmptrx(jr, 0, mr, a, bb, c, b, tcos, d, w)
    if (jr == 1) then
        q(:mr, 1) = b(:mr)
        goto 194
    end if
    q(:mr, 1) = 0.5_wp * q(:mr, 1) - q(:mr, jm1) + b(:mr)
    goto 194
end if
if (n == 2) then
    !
    !     case  n = 2
    !
    b(:mr) = q(:mr, 1)
    tcos(1) = 0.0_wp
    call cmptrx(1, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = b(:mr)
    b(:mr) = 2.0_wp * (q(:mr, 2)+b(:mr))*fistag
    tcos(1) = cmplx((-fistag), 0.0_wp, kind=wp)
    tcos(2) = cmplx(2.0_wp, 0.0_wp, kind=wp)
    call cmptrx(2, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = q(:mr, 1) + b(:mr)
    jr = 1
    i2r = 0
    goto 194
end if
b3(:mr) = 0.0_wp
b(:mr) = q(:mr, 1) + 2.0_wp * p(iip+1:mr+iip)
q(:mr, 1) = 0.5_wp * q(:mr, 1) - q(:mr, jm1)
b2(:mr) = 2.0_wp * (q(:mr, 1)+q(:mr, nlast))
k1 = kr + jr - 1
tcos(k1+1) = cmplx(-2.0_wp, 0.0_wp, kind=wp)
k4 = k1 + 3 - istag
call cmpcsg (kr + istag - 2, 1, 0.0_wp, fnum, tcos(k4))
k4 = k1 + kr + 1
call cmpcsg (jr - 1, 1, 0.0_wp, 1.0_wp, tcos(k4))
call cmpmrg (tcos, k1, kr, k1 + kr, jr - 1, 0)
call cmpcsg (kr, 1, 0.5_wp, fden, tcos(k1+1))
k2 = kr
k4 = k1 + k2 + 1
call cmpcsg (lr, 1, 0.5_wp, fden, tcos(k4))
k3 = lr
k4 = 0
call cmptr3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
b(:mr) = b(:mr) + b2(:mr)
tcos(1) = cmplx(2.0_wp, 0.0_wp, kind=wp)
call cmptrx(1, 0, mr, a, bb, c, b, tcos, d, w)
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
            q(:mr, nlast) = p(iip+1:mr+iip)
            iip = iip - mr
        else
            q(:mr, nlast) = q(:mr, nlast) - q(:mr, jm2)
        end if
    end if
    call cmpcsg (kr, 1, 0.5_wp, fden, tcos)
    call cmpcsg (lr, 1, 0.5_wp, fden, tcos(kr+1))
    if (lr == 0) then
        b(:mr) = fistag*b(:mr)
    end if
    call cmptrx(kr, lr, mr, a, bb, c, b, tcos, d, w)
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
    call cmpcsg (jr, 1, 0.5_wp, 0.0_wp, tcos)
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
            q(:mr, j) = 0.5_wp * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
        end if
        call cmptrx(jr, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
    end do
    nrod = 1
    if (nlast + i2r <= n) nrod = 0
    if (nlastp /= nlast) goto 194
    goto 206
222 continue
    w(1) = cmplx(real(ipstor, kind=wp), 0.0_wp, kind=wp)

end associate

end subroutine cmposn


subroutine cmposp(m, n, a, bb, c, q, idimq, b, b2, b3, w, w2, w3, d, tcos, p)
    !
    ! Purpose:
    !
    !     subroutine to solve poisson equation with periodic boundary
    !     conditions.
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in) :: m
    integer (ip), intent (in) :: n
    integer (ip), intent (in) :: idimq
    complex (wp) :: a(*)
    complex (wp) :: bb(*)
    complex (wp) :: c(*)
    complex (wp) :: q(idimq,n)
    complex (wp) :: b(*)
    complex (wp) :: b2(*)
    complex (wp) :: b3(*)
    complex (wp) :: w(*)
    complex (wp) :: w2(*)
    complex (wp) :: w3(*)
    complex (wp) :: d(*)
    complex (wp) :: tcos(*)
    complex (wp) :: p(n*4)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: mr, nr, nrm1, j, nrmj, nrpj, i, lh
    real (wp)    :: ipstor
    complex (wp) :: s, t
    !-----------------------------------------------

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

        q(:mr, nr) = 2.0_wp * q(:mr, nr)
        q(:mr, n) = 2.0_wp * q(:mr, n)

        call cmposd (mr, nrm1, 1, a, bb, c, q, idimq, b, w, d, tcos, p)

        ipstor = real(w(1), kind=wp)

        call cmposn(mr, nr + 1, 1, 1, a, bb, c, q(1, nr), idimq, b, b2, &
            b3, w, w2, w3, d, tcos, p)

        ipstor = max(ipstor, real(w(1), kind=wp))

        do j = 1, nrm1
            nrmj = nr - j
            nrpj = nr + j
            do i = 1, mr
                s = 0.5_wp * (q(i, nrpj)+q(i, nrmj))
                t = 0.5_wp * (q(i, nrpj)-q(i, nrmj))
                q(i, nrmj) = s
                q(i, nrpj) = t
            end do
        end do

        q(:mr, nr) = 0.5_wp * q(:mr, nr)
        q(:mr, n) = 0.5_wp * q(:mr, n)

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

        q(:mr, nr) = 2.0_wp * q(:mr, nr)
        lh = nrm1/2

        do j = 1, lh
            nrmj = nr - j
            do i = 1, mr
                s = q(i, j)
                q(i, j) = q(i, nrmj)
                q(i, nrmj) = s
            end do
        end do

        call cmposd(mr, nrm1, 2, a, bb, c, q, idimq, b, w, d, tcos, p)

        ipstor = real(w(1), kind=wp)

        call cmposn(mr, nr, 2, 1, a, bb, c, q(1, nr), idimq, b, b2, b3, &
            w, w2, w3, d, tcos, p)

        ipstor = max(ipstor, real(w(1), kind=wp))

        do j = 1, nrm1
            nrpj = nr + j
            do i = 1, mr
                s = 0.5_wp * (q(i, nrpj)+q(i, j))
                t = 0.5_wp * (q(i, nrpj)-q(i, j))
                q(i, nrpj) = t
                q(i, j) = s
            end do
        end do

        q(:mr, nr) = 0.5_wp * q(:mr, nr)

        do j = 1, lh
            nrmj = nr - j
            do i = 1, mr
                s = q(i, j)
                q(i, j) = q(i, nrmj)
                q(i, nrmj) = s
            end do
        end do
    end if

    w(1) = cmplx(ipstor, 0.0_wp, kind=wp)

end subroutine cmposp


pure subroutine cmpcsg(n, ijump, fnum, fden, a)
    !
    ! Purpose:
    !
    !     this subroutine computes required cosine values in ascending
    !     order.  when ijump > 1 the routine computes values
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
    !                                in cmposn only.
    !
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)  :: n
    integer (ip), intent (in)  :: ijump
    real (wp),    intent (in)  :: fnum
    real (wp),    intent (in)  :: fden
    complex (wp), intent (out) :: a(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip)         :: k3, k4, k, k1, k5, i, k2, np1
    real (wp)            :: pibyn, x, y
    !-----------------------------------------------

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
                    a(k2) = cmplx((-2.0_wp * cos(x*pibyn)), 0.0_wp, kind=wp)
                end do
            end do
        else
            np1 = n + 1
            y = PI/(real(n, kind=wp) + fden)
            do i = 1, n
                x = real(np1 - i, kind=wp) - fnum
                a(i) = cmplx(2.0_wp * cos(x*y), 0.0_wp, kind=wp)
            end do
        end if
    end if

end subroutine cmpcsg


subroutine cmpmrg(tcos, i1, m1, i2, m2, i3)
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
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: i1
    integer (ip), intent (in)     :: m1
    integer (ip), intent (in)     :: i2
    integer (ip), intent (in)     :: m2
    integer (ip), intent (in)     :: i3
    complex (wp), intent (in out) :: tcos(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: j11, j3, j1, j2, j, l, k, m
    complex (wp) :: x, y
    !-----------------------------------------------

    j1 = 1
    j2 = 1
    j = i3

    if_construct: if (m1 /= 0) then
        if (m2 /= 0) then
            loop_101: do
                j11 = j1
                j3 = max(m1, j11)
                block_construct: block
                    do j1 = j11, j3
                        j = j + 1
                        l = j1 + i1
                        x = tcos(l)
                        l = j2 + i2
                        y = tcos(l)
                        if (real(x - y, kind=wp) > 0.0_wp) exit block_construct
                        tcos(j) = x
                    end do
                    if (j2 > m2) return
                    exit if_construct
                end block block_construct
                tcos(j) = y
                j2 = j2 + 1
                if (j2 > m2) exit loop_101
            end do loop_101
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


end subroutine cmpmrg


subroutine cmptrx(idegbr, idegcr, m, a, b, c, y, tcos, d, w)
    !
    !     Subroutine to solve a system of linear equations where the
    !     coefficient matrix is a rational function in the matrix given by
    !     tridiagonal  ( . . . , a(i), b(i), c(i), . . . ).
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: idegbr
    integer (ip), intent (in)     :: idegcr
    integer (ip), intent (in)     :: m
    complex (wp), intent (in)     :: a(m)
    complex (wp), intent (in)     :: b(m)
    complex (wp), intent (in)     :: c(m)
    complex (wp), intent (in out) :: y(m)
    complex (wp), intent (in)     :: tcos(*)
    complex (wp), intent (out)    :: d(m)
    complex (wp), intent (out)    :: w(m)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: mm1, ifb, ifc, l, lint, k, i, iip
    complex (wp) :: x, xx, z
    !-----------------------------------------------

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

        z = 1.0_wp/(b(1)-x)
        d(1) = c(1)*z
        y(1) = y(1)*z

        do i = 2, mm1
            z = 1.0_wp/(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y(i) = (y(i)-a(i)*y(i-1))*z
        end do

        z = b(m) - x - a(m)*d(mm1)

        if (abs(z) == 0.0_wp) then
            y(m) = 0.0_wp
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

end subroutine cmptrx



subroutine cmptr3(m, a, b, c, k, y1, y2, y3, tcos, d, w1, w2, w3)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: m
    integer (ip), intent (in)     :: k(4)
    complex (wp), intent (in)     :: a(m)
    complex (wp), intent (in)     :: b(m)
    complex (wp), intent (in)     :: c(m)
    complex (wp), intent (in out) :: y1(m)
    complex (wp), intent (in out) :: y2(m)
    complex (wp), intent (in out) :: y3(m)
    complex (wp), intent (in)     :: tcos(*)
    complex (wp), intent (out)    :: d(m)
    complex (wp), intent (out)    :: w1(m)
    complex (wp), intent (out)    :: w2(m)
    complex (wp), intent (out)    :: w3(m)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: mm1, k1, k2, k3, k4
    integer (ip) :: if1, if2, if3, if4, k2k3k4, l1, l2
    integer (ip) :: l3, lint1, lint2, lint3
    integer (ip) :: kint1, kint2, kint3, n, i, iip
    complex (wp) :: x, xx, z
    !-----------------------------------------------
    !
    !     subroutine to solve tridiagonal systems
    !
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

        z = 1.0_wp/(b(1)-x)
        d(1) = c(1)*z
        y1(1) = y1(1)*z
        y2(1) = y2(1)*z
        y3(1) = y3(1)*z

        do i = 2, m
            z = 1.0_wp/(b(i)-x-a(i)*d(i-1))
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

end subroutine cmptr3



end module module_cmgnbn
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
