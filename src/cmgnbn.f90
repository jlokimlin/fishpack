module module_cmgnbn

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: cmgnbn


contains


    subroutine cmgnbn(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
        !
        !     file cmgnbn.f
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
        !     SUBROUTINE cmgnbn (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, ierror)
        !
        !
        ! DIMENSION OF           A(M), B(M), C(M), Y(IDIMY, N)
        ! ARGUMENTS
        !
        ! LATEST REVISION        NOVEMBER 2004
        !
        ! PURPOSE                THE NAME OF THIS PACKAGE IS A MNEMONIC FOR THE
        !                        COMPLEX GENERALIZED BUNEMAN ALGORITHM.
        !                        IT SOLVES THE COMPLEX LINEAR SYSTEM OF EQUATION
        !
        !                        A(I)*X(I-1, J) + B(I)*X(I, J) + C(I)*X(I+1, J)
        !                        + X(I, J-1) - 2.0_wp * X(I, J) + X(I, J+1) = Y(I, J)
        !
        !                        FOR I = 1, 2, ..., M  AND  J = 1, 2, ..., N.
        !
        !                        INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
        !                        I.E., X(0, J) = X(M, J) AND X(M+1, J) = X(1, J),
        !                        AND X(I, 0) MAY EQUAL 0, X(I, 2), OR X(I, N),
        !                        AND X(I, N+1) MAY EQUAL 0, X(I, N-1), OR X(I, 1)
        !                        DEPENDING ON AN INPUT PARAMETER.
        !
        ! USAGE                  CALL cmgnbn (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y,
        !                                     ierror)
        !
        ! ARGUMENTS
        !
        ! ON INPUT               NPEROD
        !
        !                          INDICATES THE VALUES THAT X(I, 0) AND
        !                          X(I, N+1) ARE ASSUMED TO HAVE.
        !
        !                          = 0  IF X(I, 0) = X(I, N) AND X(I, N+1) =
        !                               X(I, 1).
        !                          = 1  IF X(I, 0) = X(I, N+1) = 0  .
        !                          = 2  IF X(I, 0) = 0 AND X(I, N+1) = X(I, N-1).
        !                          = 3  IF X(I, 0) = X(I, 2) AND X(I, N+1) =
        !                               X(I, N-1).
        !                          = 4  IF X(I, 0) = X(I, 2) AND X(I, N+1) = 0.
        !
        !                        N
        !                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
        !                          N MUST BE GREATER THAN 2.
        !
        !                        MPEROD
        !                          = 0 IF A(1) AND C(M) ARE NOT ZERO
        !                          = 1 IF A(1) = C(M) = 0
        !
        !                        M
        !                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
        !                          N MUST BE GREATER THAN 2.
        !
        !                        A, B, C
        !                          ONE-DIMENSIONAL COMPLEX ARRAYS OF LENGTH M
        !                          THAT SPECIFY THE COEFFICIENTS IN THE LINEAR
        !                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0
        !                          THE ARRAY ELEMENTS MUST NOT DEPEND UPON
        !                          THE INDEX I, BUT MUST BE CONSTANT.
        !                          SPECIFICALLY, THE SUBROUTINE CHECKS THE
        !                          FOLLOWING CONDITION .
        !
        !                            A(I) = C(1)
        !                            C(I) = C(1)
        !                            B(I) = B(1)
        !
        !                          FOR I=1, 2, ..., M.
        !
        !                        IDIMY
        !                          THE ROW (OR FIRST) DIMENSION OF THE
        !                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS
        !                          IN THE PROGRAM CALLING cmgnbn.
        !                          THIS PARAMETER IS USED TO SPECIFY THE
        !                          VARIABLE DIMENSION OF Y.
        !                          IDIMY MUST BE AT LEAST M.
        !
        !                        Y
        !                          A TWO-DIMENSIONAL COMPLEX ARRAY THAT
        !                          SPECIFIES THE VALUES OF THE RIGHT SIDE
        !                          OF THE LINEAR SYSTEM OF EQUATIONS GIVEN
        !                          ABOVE.
        !                          Y MUST BE DIMENSIONED AT LEAST M*N.
        !
        !
        !  ON OUTPUT             Y
        !
        !                          CONTAINS THE SOLUTION X.
        !
        !                        ierror
        !                          AN ERROR FLAG WHICH INDICATES INVALID
        !                          INPUT PARAMETERS  EXCEPT FOR NUMBER
        !                          ZERO, A SOLUTION IS NOT ATTEMPTED.
        !
        !                          = 0  NO ERROR.
        !                          = 1  M .LE. 2  .
        !                          = 2  N .LE. 2
        !                          = 3  IDIMY .LT. M
        !                          = 4  NPEROD .LT. 0 OR NPEROD .GT. 4
        !                          = 5  MPEROD .LT. 0 OR MPEROD .GT. 1
        !                          = 6  A(I) .NE. C(1) OR C(I) .NE. C(1) OR
        !                               B(I) .NE. B(1) FOR
        !                               SOME I=1, 2, ..., M.
        !                          = 7  A(1) .NE. 0 OR C(M) .NE. 0 AND
        !                                 MPEROD = 1
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
        !
        ! SPECIAL CONDITONS      NONE
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED LIBRARY       comf.f, fish.f
        ! FILES
        !
        ! LANGUAGE               FORTRAN 90
        !
        ! HISTORY                WRITTEN IN 1979 BY ROLAND SWEET OF NCAR'S
        !                        SCIENTIFIC COMPUTING DIVISION.  MADE AVAILABLE
        !                        ON NCAR'S PUBLIC LIBRARIES IN JANUARY, 1980.
        !                        Revised in June 2004 by John Adams using
        !                        Fortran 90 dynamically allocated work space.
        !
        ! ALGORITHM              THE LINEAR SYSTEM IS SOLVED BY A CYCLIC
        !                        REDUCTION ALGORITHM DESCRIBED IN THE
        !                        REFERENCE BELOW.
        !
        ! PORTABILITY            FORTRAN 90.  ALL MACHINE DEPENDENT CONSTANTS
        !                        ARE DEFINED IN FUNCTION P1MACH.
        !
        ! REFERENCES             SWEET, R., 'A CYCLIC REDUCTION ALGORITHM FOR
        !                        SOLVING BLOCK TRIDIAGONAL SYSTEMS OF ARBITRARY
        !                        DIMENSIONS, ' SIAM J. ON NUMER. ANAL.,
        !                          14(SEPT., 1977), PP. 706-720.
        !
        ! ACCURACY               THIS TEST WAS PERFORMED ON A Platform with
        !                        64 bit floating point arithmetic.
        !                        A UNIFORM RANDOM NUMBER GENERATOR WAS USED
        !                        TO CREATE A SOLUTION ARRAY X FOR THE SYSTEM
        !                        GIVEN IN THE 'PURPOSE' DESCRIPTION ABOVE
        !                        WITH
        !                          A(I) = C(I) = -0.5_wp * B(I) = 1, I=1, 2, ..., M
        !
        !                        AND, WHEN MPEROD = 1
        !
        !                          A(1) = C(M) = 0
        !                          A(M) = C(1) = 2.
        !
        !                        THE SOLUTION X WAS SUBSTITUTED INTO THE
        !                        GIVEN SYSTEM  AND A RIGHT SIDE Y WAS
        !                        COMPUTED.  USING THIS ARRAY Y, SUBROUTINE
        !                        cmgnbn WAS CALLED TO PRODUCE APPROXIMATE
        !                        SOLUTION Z.  THEN RELATIVE ERROR
        !                          E = MAX(abs(Z(I, J)-X(I, J)))/
        !                              MAX(abs(X(I, J)))
        !                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
        !                        OVER I=1, 2, ..., M AND J=1, ..., N.
        !
        !                        THE VALUE OF E IS GIVEN IN THE TABLE
        !                        BELOW FOR SOME TYPICAL VALUES OF M AND N.
        !
        !                   M (=N)    MPEROD    NPEROD       E
        !                   ------    ------    ------     ------
        !
        !                     31        0         0        1.E-12
        !                     31        1         1        4.E-13
        !                     31        1         3        2.E-12
        !                     32        0         0        7.E-14
        !                     32        1         1        5.E-13
        !                     32        1         3        2.E-13
        !                     33        0         0        6.E-13
        !                     33        1         1        5.E-13
        !                     33        1         3        3.E-12
        !                     63        0         0        5.E-12
        !                     63        1         1        6.E-13
        !                     63        1         3        1.E-11
        !                     64        0         0        1.E-13
        !                     64        1         1        3.E-12
        !                     64        1         3        3.E-13
        !                     65        0         0        2.E-12
        !                     65        1         1        5.E-13
        !                     65        1         3        1.E-11
        !
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip),             intent (in)     :: nperod
        integer (ip),             intent (in)     :: n
        integer (ip),             intent (in)     :: mperod
        integer (ip),             intent (in)     :: m
        integer (ip),             intent (in)     :: idimy
        integer (ip),             intent (out)    :: ierror
        complex (wp), contiguous, intent (in)     :: a(:)
        complex (wp), contiguous, intent (in)     :: b(:)
        complex (wp), contiguous, intent (in)     :: c(:)
        complex (wp), contiguous, intent (in out) :: y(:,:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip)             :: i !! counter
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------

        ! initialize error flag
        ierror = 0

        ! Check input arguments
        if (m <= 2) ierror = 1
        if (n <= 2) ierror = 2
        if (idimy < m) ierror = 3
        if (nperod < 0 .or. nperod > 4) ierror = 4
        if (mperod < 0 .or. mperod > 1) ierror = 5

        if (mperod /= 1) then
            do i = 2, m
                if (abs(a(i)-c(1)) /= 0.) then
                    ierror = 6
                    exit
                end if
                if (abs(c(i)-c(1)) /= 0.) then
                    ierror = 6
                    exit
                end if
                if (abs(b(i)-b(1)) /= 0.) then
                    ierror = 6
                    exit
                end if
            end do
        end if

        if (abs(a(1))/=0. .and. abs(c(m))/=0.) then
            ierror = 7
        end if

        if (ierror /= 0) return

        ! allocate required workspace
        associate( &
            icwk => (10 + int(log(real(n))/log(2.0)))*m + 4*n, &
            irwk => 0 &
            )
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! Solve system
        associate( cxw => workspace%complex_workspace )
            call cmgnbnn(nperod, n, mperod, m, a, b, c, idimy, y, cxw)
        end associate

        ! Release memory
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
        complex (wp), intent (in)     :: a(*)
        complex (wp), intent (in)     :: b(*)
        complex (wp), intent (in)     :: c(*)
        complex (wp), intent (in out) :: y(idimy, *)
        complex (wp), intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: iwba, iwbb, iwbc, iwb2, iwb3, iww1, iww2, iww3, iwd, &
            iwtcos, iwp, i, k, j, mp, np, ipstor, irev, mh, mhm1, modd, nby2, mskip
        complex (wp) :: a1
        !-----------------------------------------------
        iwba = m + 1
        iwbb = iwba + m
        iwbc = iwbb + m
        iwb2 = iwbc + m
        iwb3 = iwb2 + m
        iww1 = iwb3 + m
        iww2 = iww1 + m
        iww3 = iww2 + m
        iwd = iww3 + m
        iwtcos = iwd + m
        iwp = iwtcos + 4*n
        do i = 1, m
            k = iwba + i - 1
            w(k) = -a(i)
            k = iwbc + i - 1
            w(k) = -c(i)
            k = iwbb + i - 1
            w(k) = 2. - b(i)
            y(i, :n) = -y(i, :n)
        end do
        mp = mperod + 1
        np = nperod + 1
        go to (114, 107) mp
107 continue
    go to (108, 109, 110, 111, 123) np
108 continue
    call cmposp (m, n, w(iwba), w(iwbb), w(iwbc), y, idimy, w, w(iwb2) &
        , w(iwb3), w(iww1), w(iww2), w(iww3), w(iwd), w(iwtcos), w(iwp) &
        )
    go to 112
109 continue
    call cmposd (m, n, 1, w(iwba), w(iwbb), w(iwbc), y, idimy, w, w( &
        iww1), w(iwd), w(iwtcos), w(iwp))
    go to 112
110 continue
    call cmposn(m, n, 1, 2, w(iwba), w(iwbb), w(iwbc), y, idimy, w, w &
        (iwb2), w(iwb3), w(iww1), w(iww2), w(iww3), w(iwd), w(iwtcos), &
        w(iwp))
    go to 112
111 continue
    call cmposn(m, n, 1, 1, w(iwba), w(iwbb), w(iwbc), y, idimy, w, w &
        (iwb2), w(iwb3), w(iww1), w(iww2), w(iww3), w(iwd), w(iwtcos), &
        w(iwp))
112 continue
    ipstor = real(w(iww1))
    irev = 2
    if (nperod == 4) go to 124
113 continue
    go to (127, 133) mp
114 continue
    mh = (m + 1)/2
    mhm1 = mh - 1
    modd = 1
    if (mh*2 == m) modd = 2
    do j = 1, n
        do i = 1, mhm1
            w(i) = y(mh-i, j) - y(i+mh, j)
            w(i+mh) = y(mh-i, j) + y(i+mh, j)
        end do
        w(mh) = 2.0_wp * y(mh, j)
        go to (117, 116) modd
116 continue
    w(m) = 2.0_wp * y(m, j)
117 continue
    y(:m, j) = w(:m)
end do
k = iwbc + mhm1 - 1
i = iwba + mhm1
w(k) = (0., 0.)
w(i) = (0., 0.)
w(k+1) = 2.0_wp * w(k+1)
select case (modd)
    case default
        k = iwbb + mhm1 - 1
        w(k) = w(k) - w(i-1)
        w(iwbc-1) = w(iwbc-1) + w(iwbb-1)
    case (2)
        w(iwbb-1) = w(k+1)
end select
122 continue
    go to 107
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
            a1 = y(i, j)
            y(i, j) = y(i, mskip)
            y(i, mskip) = a1
        end do
    end do
    go to (110, 113) irev
127 continue
    do j = 1, n
        w(mh-1:mh-mhm1:(-1)) = 0.5_wp * (y(mh+1:mhm1+mh, j)+y(:mhm1, j))
        w(mh+1:mhm1+mh) = 0.5_wp * (y(mh+1:mhm1+mh, j)-y(:mhm1, j))
        w(mh) = 0.5_wp * y(mh, j)
        go to (130, 129) modd
129 continue
    w(m) = 0.5_wp * y(m, j)
130 continue
    y(:m, j) = w(:m)
end do
133 continue
    w(1) = cmplx(real(ipstor + iwp - 1), 0.)
    return
end subroutine cmgnbnn


subroutine cmposd(mr, nr, istag, ba, bb, bc, q, idimq, b, w, d, tcos, p)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: mr
    integer (ip), intent (in)     :: nr
    integer (ip), intent (in)     :: istag
    integer (ip), intent (in)     :: idimq
    complex (wp), intent (in)     :: ba(*)
    complex (wp), intent (in)     :: bb(*)
    complex (wp), intent (in)     :: bc(*)
    complex (wp), intent (in out) :: q(idimq, 1)
    complex (wp), intent (in out) :: b(*)
    complex (wp), intent (in out) :: w(*)
    complex (wp), intent (in out) :: d(*)
    complex (wp), intent (in out) :: tcos(*)
    complex (wp), intent (in out) :: p(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: m, n, ip_rename, ipstor, jsh, kr, irreg, jstsav, i, lr, nun, &
        jst, jsp, l, nodd, j, jm1, jp1, jm2, jp2, jm3, jp3, noddpr &
        , krpi, ideg, jdeg
    real (wp) :: fi
    complex (wp) :: t
    !-----------------------------------------------
    !
    !     SUBROUTINE TO SOLVE POISSON'S EQUATION FOR DIRICHLET BOUNDARY
    !     CONDITIONS.
    !
    !     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A.
    !     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A+I.
    !
    m = mr
    n = nr
    fi = 1./real(istag)
    ip_rename = -m
    ipstor = 0
    jsh = 0
    select case (istag)
        case default
            kr = 0
            irreg = 1
            if (n > 1) go to 106
            tcos(1) = (0., 0.)
        case (2)
            kr = 1
            jstsav = 1
            irreg = 2
            if (n > 1) go to 106
            tcos(1) = CMPLX(-1., 0.)
    end select
103 continue
    b(:m) = Q(:m, 1)
    call cmptrx (1, 0, m, ba, bb, bc, b, tcos, d, w)
    q(:m, 1) = B(:m)
    go to 183
106 continue
    lr = 0
    do i = 1, m
        p(i) = CMPLX(0., 0.)
    end do
    nun = n
    jst = 1
    jsp = n
!
!     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
!
108 continue
    l = 2*jst
    nodd = 2 - 2*((nun + 1)/2) + nun
    !
    !     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
    !
    select case (nodd)
        case default
            jsp = jsp - l
        case (1)
            jsp = jsp - jst
            if (irreg /= 1) jsp = jsp - l
    end select
111 continue
    call cmpcsg (jst, 1, 0.5, 0.0, tcos)
    if (l <= jsp) then
        do j = l, jsp, l
            jm1 = j - jsh
            jp1 = j + jsh
            jm2 = j - jst
            jp2 = j + jst
            jm3 = jm2 - jsh
            jp3 = jp2 + jsh
            if (jst == 1) then
                b(:m) = 2.0_wp * Q(:m, j)
                q(:m, j) = Q(:m, jm2) + Q(:m, jp2)
            else
                do i = 1, m
                    t = Q(i, j) - Q(i, jm1) - Q(i, jp1) + Q(i, jm2) + Q(i, jp2)
                    b(i) = t + Q(i, j) - Q(i, jm3) - Q(i, jp3)
                    q(i, j) = t
                end do
            end if
            call cmptrx (jst, 0, m, ba, bb, bc, b, tcos, d, w)
            q(:m, j) = Q(:m, j) + B(:m)
        end do
    end if
    !
    !     REDUCTION FOR LAST UNKNOWN
    !
    select case (nodd)
        case default
            go to (152, 120) irreg
        !
        !     ODD NUMBER OF UNKNOWNS
        !
120     continue
        jsp = jsp + l
        j = jsp
        jm1 = j - jsh
        jp1 = j + jsh
        jm2 = j - jst
        jp2 = j + jst
        jm3 = jm2 - jsh
        go to (123, 121) istag
121 continue
    if (jst /= 1) go to 123
    do i = 1, m
        b(i) = Q(i, j)
        q(i, j) = CMPLX(0., 0.)
    end do
    go to 130
123 continue
    select case (noddpr)
        case default
            b(:m) = 0.5_wp * (Q(:m, jm2)-Q(:m, jm1)-Q(:m, jm3)) + P(ip_rename+1:m+ip_rename) &
                + Q(:m, j)
        case (2)
            b(:m) = 0.5_wp * (Q(:m, jm2)-Q(:m, jm1)-Q(:m, jm3)) + Q(:m, jp2) - Q( &
                :m, jp1) + Q(:m, j)
    end select
128 continue
    q(:m, j) = 0.5_wp * (Q(:m, j)-Q(:m, jm1)-Q(:m, jp1))
130 continue
    call cmptrx (jst, 0, m, ba, bb, bc, b, tcos, d, w)
    ip_rename = ip_rename + m
    ipstor = max(ipstor, ip_rename + m)
    p(ip_rename+1:m+ip_rename) = Q(:m, j) + B(:m)
    b(:m) = Q(:m, jp2) + P(ip_rename+1:m+ip_rename)
    if (lr == 0) then
        do i = 1, jst
            krpi = kr + i
            tcos(krpi) = TCOS(i)
        end do
    else
        call cmpcsg (lr, jstsav, 0., fi, TCOS(jst+1))
        call cmpmrg (tcos, 0, jst, jst, lr, kr)
    end if
    call cmpcsg (kr, jstsav, 0.0, fi, tcos)
    call cmptrx (kr, kr, m, ba, bb, bc, b, tcos, d, w)
    q(:m, j) = Q(:m, jm2) + B(:m) + P(ip_rename+1:m+ip_rename)
    lr = kr
    kr = kr + l
!
!     EVEN NUMBER OF UNKNOWNS
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
            call cmpcsg (kr, jstsav, 0.0, fi, tcos)
            call cmpcsg (lr, jstsav, 0.0, fi, TCOS(kr+1))
            ideg = kr
            kr = kr + jst
    end select
139 continue
    if (jst == 1) then
        irreg = 2
        b(:m) = Q(:m, j)
        q(:m, j) = Q(:m, jm2)
    else
        b(:m) = Q(:m, j) + 0.5_wp * (Q(:m, jm2)-Q(:m, jm1)-Q(:m, jm3))
        select case (irreg)
            case default
                q(:m, j) = Q(:m, jm2) + 0.5_wp * (Q(:m, j)-Q(:m, jm1)-Q(:m, jp1))
                irreg = 2
            case (2)
                select case (noddpr)
                    case default
                        q(:m, j) = Q(:m, jm2) + P(ip_rename+1:m+ip_rename)
                        ip_rename = ip_rename - m
                    case (2)
                        q(:m, j) = Q(:m, jm2) + Q(:m, j) - Q(:m, jm1)
                end select
        end select
    end if
150 continue
    call cmptrx (ideg, lr, m, ba, bb, bc, b, tcos, d, w)
    q(:m, j) = Q(:m, j) + B(:m)
end select
152 continue
    nun = nun/2
    noddpr = nodd
    jsh = jst
    jst = 2*jst
    if (nun >= 2) go to 108
    !
    !     START SOLUTION.
    !
    j = jsp
    b(:m) = Q(:m, j)
    select case (irreg)
        case default
            call cmpcsg (jst, 1, 0.5, 0.0, tcos)
            ideg = jst
        case (2)
            kr = lr + jst
            call cmpcsg (kr, jstsav, 0.0, fi, tcos)
            call cmpcsg (lr, jstsav, 0.0, fi, TCOS(kr+1))
            ideg = kr
    end select
156 continue
    call cmptrx (ideg, lr, m, ba, bb, bc, b, tcos, d, w)
    jm1 = j - jsh
    jp1 = j + jsh
    select case (irreg)
        case default
            q(:m, j) = 0.5_wp * (Q(:m, j)-Q(:m, jm1)-Q(:m, jp1)) + B(:m)
        case (2)
            select case (noddpr)
                case default
                    q(:m, j) = P(ip_rename+1:m+ip_rename) + B(:m)
                    ip_rename = ip_rename - m
                case (2)
                    q(:m, j) = Q(:m, j) - Q(:m, jm1) + B(:m)
            end select
    end select
164 continue
    jst = jst/2
    jsh = jst/2
    nun = 2*nun
    if (nun > n) go to 183
    do j = jst, n, l
        jm1 = j - jsh
        jp1 = j + jsh
        jm2 = j - jst
        jp2 = j + jst
        if (j <= jst) then
            b(:m) = Q(:m, j) + Q(:m, jp2)
        else
            if (jp2 <= n) go to 168
            b(:m) = Q(:m, j) + Q(:m, jm2)
            if (jst < jstsav) irreg = 1
            go to (170, 171) irreg
168     continue
        b(:m) = Q(:m, j) + Q(:m, jm2) + Q(:m, jp2)
    end if
170 continue
    call cmpcsg (jst, 1, 0.5, 0.0, tcos)
    ideg = jst
    jdeg = 0
    go to 172
171 continue
    if (j + l > n) lr = lr - jst
    kr = jst + lr
    call cmpcsg (kr, jstsav, 0.0, fi, tcos)
    call cmpcsg (lr, jstsav, 0.0, fi, TCOS(kr+1))
    ideg = kr
    jdeg = lr
172 continue
    call cmptrx (ideg, jdeg, m, ba, bb, bc, b, tcos, d, w)
    if (jst <= 1) then
        q(:m, j) = B(:m)
    else
        if (jp2 > n) go to 177
175 continue
    q(:m, j) = 0.5_wp * (Q(:m, j)-Q(:m, jm1)-Q(:m, jp1)) + B(:m)
    cycle
177 continue
    go to (175, 178) irreg
178 continue
    if (j + jsh <= n) then
        q(:m, j) = B(:m) + P(ip_rename+1:m+ip_rename)
        ip_rename = ip_rename - m
    else
        q(:m, j) = B(:m) + Q(:m, j) - Q(:m, jm1)
    end if
end if
end do
l = l/2
go to 164
183 continue
    w(1) = CMPLX(real(ipstor), 0.)
    return
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
    complex (wp) :: a(*)
    complex (wp) :: bb(*)
    complex (wp) :: c(*)
    complex (wp), intent (in out) :: q(idimq, *)
    complex (wp) :: b(*)
    complex (wp) :: b2(*)
    complex (wp) :: b3(*)
    complex (wp) :: w(*)
    complex (wp) :: w2(*)
    complex (wp) :: w3(*)
    complex (wp) :: d(*)
    complex (wp) :: tcos(*)
    complex (wp), intent (in out) :: p(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: k(4)
    integer (ip) :: k1, k2, k3, k4, mr, ip_rename
    integer (ip) :: ipstor, i2r, jr, nr, nlast, kr
    integer (ip) :: lr, i, nrod, jstart, jstop, i2rby2, j, jp1, jp2, jp3, jm1
    integer (ip) :: jm2, jm3, nrodpr, ii, i1, i2, jr2, nlastp, jstep
    real (wp)    :: fistag, fnum, fden
    complex (wp) :: fi, t
    !-----------------------------------------------

    equivalence (k(1), k1), (k(2), k2), (k(3), k3), (k(4), k4)
    fistag = 3 - istag
    fnum = 1./real(istag)
    fden = 0.5_wp * real(istag - 1)
    mr = m
    ip_rename = -mr
    ipstor = 0
    i2r = 1
    jr = 2
    nr = n
    nlast = n
    kr = 1
    lr = 0
    go to (101, 103) istag
101 continue
    q(:mr, n) = 0.5_wp * q(:mr, n)
    go to (103, 104) mixbnd
103 continue
    if (n <= 3) go to 155
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
107 continue
    jstop = nlast - jr
    if (nrod == 0) jstop = jstop - i2r
    call cmpcsg (i2r, 1, 0.5, 0.0, tcos)
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
            call cmptrx (i2r, 0, mr, a, bb, c, b, tcos, d, w)
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
                q(:mr, j) = q(:mr, jm2) + p(ip_rename+1:mr+ip_rename)
                ip_rename = ip_rename - mr
            else
                q(:mr, j) = q(:mr, j) - q(:mr, jm1) + q(:mr, jm2)
            end if
            if (lr /= 0) then
                call cmpcsg (lr, 1, 0.5, fden, tcos(kr+1))
            else
                b(:mr) = fistag*b(:mr)
            end if
        end if
        call cmpcsg (kr, 1, 0.5, fden, tcos)
        call cmptrx (kr, lr, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
        kr = kr + i2r
    else
        jp1 = j + i2rby2
        jp2 = j + i2r
        if (i2r == 1) then
            b(:mr) = q(:mr, j)
            call cmptrx (1, 0, mr, a, bb, c, b, tcos, d, w)
            ip_rename = 0
            ipstor = mr
            select case (istag)
                case default
                    p(:mr) = b(:mr)
                    b(:mr) = b(:mr) + q(:mr, n)
                    tcos(1) = cmplx(1., 0.)
                    tcos(2) = cmplx(0., 0.)
                    call cmptrx (1, 1, mr, a, bb, c, b, tcos, d, w)
                    q(:mr, j) = q(:mr, jm2) + p(:mr) + b(:mr)
                    go to 150
                case (1)
                    p(:mr) = b(:mr)
                    q(:mr, j) = q(:mr, jm2) + 2.0_wp * q(:mr, jp2) + 3.*b(:mr)
                    go to 150
            end select
        end if
        b(:mr) = q(:mr, j) + 0.5_wp * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))
        if (nrodpr == 0) then
            b(:mr) = b(:mr) + p(ip_rename+1:mr+ip_rename)
        else
            b(:mr) = b(:mr) + q(:mr, jp2) - q(:mr, jp1)
        end if
        call cmptrx (i2r, 0, mr, a, bb, c, b, tcos, d, w)
        ip_rename = ip_rename + mr
        ipstor = max(ipstor, ip_rename + mr)
        p(ip_rename+1:mr+ip_rename) = b(:mr) + 0.5_wp * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
        b(:mr) = p(ip_rename+1:mr+ip_rename) + q(:mr, jp2)
        if (lr /= 0) then
            call cmpcsg (lr, 1, 0.5, fden, tcos(i2r+1))
            call cmpmrg (tcos, 0, i2r, i2r, lr, kr)
        else
            do i = 1, i2r
                ii = kr + i
                tcos(ii) = tcos(i)
            end do
        end if
        call cmpcsg (kr, 1, 0.5, fden, tcos)
        if (lr == 0) then
            go to (146, 145) istag
        end if
145 continue
    call cmptrx (kr, kr, mr, a, bb, c, b, tcos, d, w)
    go to 148
146 continue
    b(:mr) = fistag*b(:mr)
148 continue
    q(:mr, j) = q(:mr, jm2) + p(ip_rename+1:mr+ip_rename) + b(:mr)
150 continue
    lr = kr
    kr = kr + jr
end if
select case (mixbnd)
    case default
        nr = (nlast - 1)/jr + 1
        if (nr <= 3) go to 155
    case (2)
        nr = nlast/jr
        if (nr <= 1) go to 192
end select
154 continue
    i2r = jr
    nrodpr = nrod
    go to 104
155 continue
    j = 1 + jr
    jm1 = j - i2r
    jp1 = j + i2r
    jm2 = nlast - i2r
    if (nr /= 2) then
        if (lr /= 0) go to 170
        if (n == 3) then
            !
            !     case n = 3.
            !
            go to (156, 168) istag
156     continue
        b(:mr) = q(:mr, 2)
        tcos(1) = cmplx(0., 0.)
        call cmptrx (1, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 2) = b(:mr)
        b(:mr) = 4.*b(:mr) + q(:mr, 1) + 2.0_wp * q(:mr, 3)
        tcos(1) = cmplx(-2., 0.)
        tcos(2) = cmplx(2., 0.)
        i1 = 2
        i2 = 0
        call cmptrx (i1, i2, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 2) = q(:mr, 2) + b(:mr)
        b(:mr) = q(:mr, 1) + 2.0_wp * q(:mr, 2)
        tcos(1) = (0., 0.)
        call cmptrx (1, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 1) = b(:mr)
        jr = 1
        i2r = 0
        go to 194
    end if
    !
    !     case n = 2**p+1
    !
    go to (162, 170) istag
162 continue
    b(:mr) = q(:mr, j) + 0.5_wp * q(:mr, 1) - q(:mr, jm1) + q(:mr, nlast) - &
        q(:mr, jm2)
    call cmpcsg (jr, 1, 0.5, 0.0, tcos)
    call cmptrx (jr, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = 0.5_wp * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1)) + b(:mr)
    b(:mr) = q(:mr, 1) + 2.0_wp * q(:mr, nlast) + 4.*q(:mr, j)
    jr2 = 2*jr
    call cmpcsg (jr, 1, 0.0, 0.0, tcos)
    tcos(jr+1:jr*2) = -tcos(jr:1:(-1))
    call cmptrx (jr2, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = q(:mr, j) + b(:mr)
    b(:mr) = q(:mr, 1) + 2.0_wp * q(:mr, j)
    call cmpcsg (jr, 1, 0.5, 0.0, tcos)
    call cmptrx (jr, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = 0.5_wp * q(:mr, 1) - q(:mr, jm1) + b(:mr)
    go to 194
!
!     case of general n with nr = 3 .
!
168 continue
    b(:mr) = q(:mr, 2)
    q(:mr, 2) = (0., 0.)
    b2(:mr) = q(:mr, 3)
    b3(:mr) = q(:mr, 1)
    jr = 1
    i2r = 0
    j = 2
    go to 177
170 continue
    b(:mr) = 0.5_wp * q(:mr, 1) - q(:mr, jm1) + q(:mr, j)
    if (nrod == 0) then
        b(:mr) = b(:mr) + p(ip_rename+1:mr+ip_rename)
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
    tcos(k1+1) = (-2., 0.)
    k4 = k1 + 3 - istag
    call cmpcsg (k2 + istag - 2, 1, 0.0, fnum, tcos(k4))
    k4 = k1 + k2 + 1
    call cmpcsg (jr - 1, 1, 0.0, 1.0, tcos(k4))
    call cmpmrg (tcos, k1, k2, k1 + k2, jr - 1, 0)
    k3 = k1 + k2 + lr
    call cmpcsg (jr, 1, 0.5, 0.0, tcos(k3+1))
    k4 = k3 + jr + 1
    call cmpcsg (kr, 1, 0.5, fden, tcos(k4))
    call cmpmrg (tcos, k3, jr, k3 + jr, kr, k1)
    if (lr /= 0) then
        call cmpcsg (lr, 1, 0.5, fden, tcos(k4))
        call cmpmrg (tcos, k3, jr, k3 + jr, lr, k3 - lr)
        call cmpcsg (kr, 1, 0.5, fden, tcos(k4))
    end if
    k3 = kr
    k4 = kr
    call cmptr3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
    b(:mr) = b(:mr) + b2(:mr) + b3(:mr)
    tcos(1) = (2., 0.)
    call cmptrx (1, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = q(:mr, j) + b(:mr)
    b(:mr) = q(:mr, 1) + 2.0_wp * q(:mr, j)
    call cmpcsg (jr, 1, 0.5, 0.0, tcos)
    call cmptrx (jr, 0, mr, a, bb, c, b, tcos, d, w)
    if (jr == 1) then
        q(:mr, 1) = b(:mr)
        go to 194
    end if
    q(:mr, 1) = 0.5_wp * q(:mr, 1) - q(:mr, jm1) + b(:mr)
    go to 194
end if
if (n == 2) then
    !
    !     case  n = 2
    !
    b(:mr) = q(:mr, 1)
    tcos(1) = (0., 0.)
    call cmptrx (1, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = b(:mr)
    b(:mr) = 2.0_wp * (q(:mr, 2)+b(:mr))*fistag
    tcos(1) = cmplx((-fistag), 0.)
    tcos(2) = cmplx(2., 0.)
    call cmptrx (2, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = q(:mr, 1) + b(:mr)
    jr = 1
    i2r = 0
    go to 194
end if
b3(:mr) = (0., 0.)
b(:mr) = q(:mr, 1) + 2.0_wp * p(ip_rename+1:mr+ip_rename)
q(:mr, 1) = 0.5_wp * q(:mr, 1) - q(:mr, jm1)
b2(:mr) = 2.0_wp * (q(:mr, 1)+q(:mr, nlast))
k1 = kr + jr - 1
tcos(k1+1) = (-2., 0.)
k4 = k1 + 3 - istag
call cmpcsg (kr + istag - 2, 1, 0.0, fnum, tcos(k4))
k4 = k1 + kr + 1
call cmpcsg (jr - 1, 1, 0.0, 1.0, tcos(k4))
call cmpmrg (tcos, k1, kr, k1 + kr, jr - 1, 0)
call cmpcsg (kr, 1, 0.5, fden, tcos(k1+1))
k2 = kr
k4 = k1 + k2 + 1
call cmpcsg (lr, 1, 0.5, fden, tcos(k4))
k3 = lr
k4 = 0
call cmptr3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
b(:mr) = b(:mr) + b2(:mr)
tcos(1) = (2., 0.)
call cmptrx (1, 0, mr, a, bb, c, b, tcos, d, w)
q(:mr, 1) = q(:mr, 1) + b(:mr)
go to 194
192 continue
    b(:mr) = q(:mr, nlast)
    go to 196
194 continue
    j = nlast - jr
    b(:mr) = q(:mr, nlast) + q(:mr, j)
196 continue
    jm2 = nlast - i2r
    if (jr == 1) then
        q(:mr, nlast) = (0., 0.)
    else
        if (nrod == 0) then
            q(:mr, nlast) = p(ip_rename+1:mr+ip_rename)
            ip_rename = ip_rename - mr
        else
            q(:mr, nlast) = q(:mr, nlast) - q(:mr, jm2)
        end if
    end if
    call cmpcsg (kr, 1, 0.5, fden, tcos)
    call cmpcsg (lr, 1, 0.5, fden, tcos(kr+1))
    if (lr == 0) then
        b(:mr) = fistag*b(:mr)
    end if
    call cmptrx (kr, lr, mr, a, bb, c, b, tcos, d, w)
    q(:mr, nlast) = q(:mr, nlast) + b(:mr)
    nlastp = nlast
206 continue
    jstep = jr
    jr = i2r
    i2r = i2r/2
    if (jr == 0) go to 222
    select case (mixbnd)
        case default
            jstart = 1 + jr
        case (2)
            jstart = jr
    end select
209 continue
    kr = kr - jr
    if (nlast + jr <= n) then
        kr = kr - jr
        nlast = nlast + jr
        jstop = nlast - jstep
    else
        jstop = nlast - jr
    end if
    lr = kr - jr
    call cmpcsg (jr, 1, 0.5, 0.0, tcos)
    do j = jstart, jstop, jstep
        jm2 = j - jr
        jp2 = j + jr
        if (j == jr) then
            b(:mr) = q(:mr, j) + q(:mr, jp2)
        else
            b(:mr) = q(:mr, j) + q(:mr, jm2) + q(:mr, jp2)
        end if
        if (jr == 1) then
            q(:mr, j) = (0., 0.)
        else
            jm1 = j - i2r
            jp1 = j + i2r
            q(:mr, j) = 0.5_wp * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
        end if
        call cmptrx (jr, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
    end do
    nrod = 1
    if (nlast + i2r <= n) nrod = 0
    if (nlastp /= nlast) go to 194
    go to 206
222 continue
    w(1) = cmplx(real(ipstor), 0.)

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
    integer (ip) :: idimq
    complex (wp) :: a(*)
    complex (wp) :: bb(*)
    complex (wp) :: c(*)
    complex (wp) :: q(idimq, 1)
    complex (wp) :: b(*)
    complex (wp) :: b2(*)
    complex (wp) :: b3(*)
    complex (wp) :: w(*)
    complex (wp) :: w2(*)
    complex (wp) :: w3(*)
    complex (wp) :: d(*)
    complex (wp) :: tcos(*)
    complex (wp) :: p(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: mr, nr, nrm1, j, nrmj, nrpj, i, ipstor, lh
    complex (wp) :: s, t
    !-----------------------------------------------

    mr = m
    nr = (n + 1)/2
    nrm1 = nr - 1
    if (2*nr == n) then
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
        ipstor = real(w(1))
        call cmposn(mr, nr + 1, 1, 1, a, bb, c, q(1, nr), idimq, b, b2 &
            , b3, w, w2, w3, d, tcos, p)
        ipstor = max(ipstor, int(real(w(1))))
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
        call cmposd (mr, nrm1, 2, a, bb, c, q, idimq, b, w, d, tcos, p)
        ipstor = real(w(1))
        call cmposn(mr, nr, 2, 1, a, bb, c, q(1, nr), idimq, b, b2, b3 &
            , w, w2, w3, d, tcos, p)
        ipstor = max(ipstor, int(real(w(1))))
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
    w(1) = cmplx(real(ipstor), 0.)

end subroutine cmposp


pure subroutine cmpcsg(n, ijump, fnum, fden, a)
    !
    ! Purpose:
    !
    !     this subroutine computes required cosine values in ascending
    !     order.  when ijump .gt. 1 the routine computes values
    !
    !        2*cos(j*pi/l) , j=1, 2, ..., l and j .ne. 0(mod n/ijump+1)
    !
    !     where l = ijump*(n/ijump+1).
    !
    !
    !     when ijump = 1 it computes
    !
    !            2*cos((j-fnum)*pi/(n+fden)) ,  j=1, 2, ... , n
    !
    !     where
    !        fnum = 0.5, fden = 0.0,  for regular reduction values
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
    integer (ip) :: k3, k4, k, k1, k5, i, k2, np1
    real (wp), parameter :: pi = acos( -1.0_wp )
    real (wp) :: pibyn, x, y
    !-----------------------------------------------

    if (n /= 0) then
        if (ijump /= 1) then
            k3 = n/ijump + 1
            k4 = k3 - 1
            pibyn = pi/(n + ijump)
            do k = 1, ijump
                k1 = (k - 1)*k3
                k5 = (k - 1)*k4
                do i = 1, k4
                    x = k1 + i
                    k2 = k5 + i
                    a(k2) = cmplx((-2.0_wp * cos(x*pibyn)), 0.)
                end do
            end do
        else
            np1 = n + 1
            y = pi/(real(n) + fden)
            do i = 1, n
                x = real(np1 - i) - fnum
                a(i) = cmplx(2.0_wp * cos(x*y), 0.)
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
    if (m1 == 0) go to 107
    if (m2 == 0) go to 104
101 continue
    j11 = j1
    j3 = max(m1, j11)
    do j1 = j11, j3
        j = j + 1
        l = j1 + i1
        x = tcos(l)
        l = j2 + i2
        y = tcos(l)
        if (real(x - y) > 0.) go to 103
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
        tcos(m) = tcos(l)
    end do
    go to 109
106 continue
    if (j2 > m2) go to 109
107 continue
    k = j - j2 + 1
    do j = j2, m2
        m = k + j
        l = j + i2
        tcos(m) = tcos(l)
    end do
109 continue

end subroutine cmpmrg


subroutine cmptrx(idegbr, idegcr, m, a, b, c, y, tcos, d, w)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in) :: idegbr
    integer (ip), intent (in) :: idegcr
    integer (ip), intent (in) :: m
    complex (wp), intent (in) :: a(*)
    complex (wp), intent (in) :: b(*)
    complex (wp), intent (in) :: c(*)
    complex (wp), intent (in out) :: y(*)
    complex (wp), intent (in) :: tcos(*)
    complex (wp), intent (in out) :: d(*)
    complex (wp), intent (in out) :: w(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: mm1, ifb, ifc, l, lint, k, i, ip_rename
    complex (wp) :: x, xx, z
    !-----------------------------------------------
    !
    !     SUBROUTINE TO SOLVE A SYSTEM OF LINEAR EQUATIONS WHERE THE
    !     COEFFICIENT MATRIX IS A RATIONAL FUNCTION IN THE MATRIX GIVEN BY
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
        if (abs(z) == 0.) then
            y(m) = (0., 0.)
        else
            y(m) = (Y(m)-A(m)*Y(mm1))/z
        end if
        do ip_rename = 1, mm1
            y(m-ip_rename) = Y(m-ip_rename) - D(m-ip_rename)*Y(m+1-ip_rename)
        end do
        if (k /= l) cycle
        y(:m) = Y(:m) + W(:m)
        lint = lint + 1
        l = (lint*ifb)/ifc
    end do

end subroutine cmptrx


subroutine cmptr3(m, a, b, c, k, y1, y2, y3, tcos, d, w1, w2, w3)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in) :: m
    integer (ip), intent (in) :: k(4)
    complex (wp), intent (in) :: a(*)
    complex (wp), intent (in) :: b(*)
    complex (wp), intent (in) :: c(*)
    complex (wp), intent (in out) :: y1(*)
    complex (wp), intent (in out) :: y2(*)
    complex (wp), intent (in out) :: y3(*)
    complex (wp), intent (in) :: tcos(*)
    complex (wp), intent (in out) :: d(*)
    complex (wp), intent (in out) :: w1(*)
    complex (wp), intent (in out) :: w2(*)
    complex (wp), intent (in out) :: w3(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: mm1, k1, k2, k3, k4, if1, if2, if3, if4, k2k3k4, l1, l2
    integer (ip) :: l3, lint1, lint2, lint3, kint1, kint2, kint3, n, i, ip_rename
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
            if (n == l1) then
                w1(:m) = y1(:m)
            end if
            if (n == l2) then
                w2(:m) = y2(:m)
            end if
            if (n == l3) then
                w3(:m) = y3(:m)
            end if
        end if
        z = 1./(b(1)-x)
        d(1) = c(1)*z
        y1(1) = y1(1)*z
        y2(1) = y2(1)*z
        y3(1) = y3(1)*z
        do i = 2, m
            z = 1./(b(i)-x-a(i)*d(i-1))
            d(i) = c(i)*z
            y1(i) = (y1(i)-a(i)*y1(i-1))*z
            y2(i) = (y2(i)-a(i)*y2(i-1))*z
            y3(i) = (y3(i)-a(i)*y3(i-1))*z
        end do
        do ip_rename = 1, mm1
            y1(m-ip_rename) = y1(m-ip_rename) - d(m-ip_rename)*y1(m+1-ip_rename)
            y2(m-ip_rename) = y2(m-ip_rename) - d(m-ip_rename)*y2(m+1-ip_rename)
            y3(m-ip_rename) = y3(m-ip_rename) - d(m-ip_rename)*y3(m+1-ip_rename)
        end do
        if (k2k3k4 == 0) cycle
        if (n == l1) then
            i = lint1 + kint1
            xx = x - tcos(i)
            y1(:m) = xx*y1(:m) + w1(:m)
            lint1 = lint1 + 1
            l1 = (lint1*if1)/if2
        end if
        if (n == l2) then
            i = lint2 + kint2
            xx = x - tcos(i)
            y2(:m) = xx*y2(:m) + w2(:m)
            lint2 = lint2 + 1
            l2 = (lint2*if1)/if3
        end if
        if (n /= l3) cycle
        i = lint3 + kint3
        xx = x - tcos(i)
        y3(:m) = xx*y3(:m) + w3(:m)
        lint3 = lint3 + 1
        l3 = (lint3*if1)/if4
    end do

end subroutine cmptr3

end module module_cmgnbn
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 changes
!-----------------------------------------------------------------------
