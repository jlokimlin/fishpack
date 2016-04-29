module module_genbun

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_gnbnaux, only: &
        cosgen, &
        merge_rename, &
        trix, &
        tri3

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: genbun
    public :: genbunn


contains


    subroutine genbun(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
        !
        !     file genbun.f
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
        !     SUBROUTINE genbun (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, ierror)
        !
        !
        ! DIMENSION OF           A(M), B(M), C(M), Y(IDIMY, N)
        ! ARGUMENTS
        !
        ! LATEST REVISION        JUNE 2004
        !
        ! PURPOSE                THE NAME OF THIS PACKAGE IS A MNEMONIC FOR THE
        !                        GENERALIZED BUNEMAN ALGORITHM.
        !
        !                        IT SOLVES THE REAL LINEAR SYSTEM OF EQUATIONS
        !
        !                        A(I)*X(I-1, J) + B(I)*X(I, J) + C(I)*X(I+1, J)
        !                        + X(I, J-1) - 2.*X(I, J) + X(I, J+1) = Y(I, J)
        !
        !                        FOR I = 1, 2, ..., M  AND  J = 1, 2, ..., N.
        !
        !                        INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
        !                        I.E., X(0, J) = X(M, J) AND X(M+1, J) = X(1, J),
        !                        AND X(I, 0) MAY EQUAL 0, X(I, 2), OR X(I, N),
        !                        AND X(I, N+1) MAY EQUAL 0, X(I, N-1), OR X(I, 1)
        !                        DEPENDING ON AN INPUT PARAMETER.
        !
        ! USAGE                  CALL genbun (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y,
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
        !                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
        !                          SPECIFY THE COEFFICIENTS IN THE LINEAR
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
        !                          IN THE PROGRAM CALLING genbun.
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
        !
        ! SPECIAL CONDITONS      NONE
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED FILES         comf.f, gnbnaux.f, fish.f
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
        !                        REFERENCE.
        !
        ! PORTABILITY            FORTRAN 90 --
        !                        THE MACHINE DEPENDENT CONSTANT PI IS
        !                        DEFINED IN FUNCTION PI_MACH.
        !
        ! REFERENCES             SWEET, R., "A CYCLIC REDUCTION ALGORITHM FOR
        !                        SOLVING BLOCK TRIDIAGONAL SYSTEMS OF ARBITRARY
        !                        DIMENSIONS, " SIAM J. ON NUMER. ANAL., 14 (1977)
        !                        PP. 706-720.
        !
        ! ACCURACY               THIS TEST WAS PERFORMED ON a platform with
        !                        64 bit floating point arithmetic.
        !                        A UNIFORM RANDOM NUMBER GENERATOR WAS USED
        !                        TO CREATE A SOLUTION ARRAY X FOR THE SYSTEM
        !                        GIVEN IN THE 'PURPOSE' DESCRIPTION ABOVE
        !                        WITH
        !                          A(I) = C(I) = -0.5*B(I) = 1, I=1, 2, ..., M
        !
        !                        AND, WHEN MPEROD = 1
        !
        !                          A(1) = C(M) = 0
        !                          A(M) = C(1) = 2.
        !
        !                        THE SOLUTION X WAS SUBSTITUTED INTO THE
        !                        GIVEN SYSTEM  AND, USING DOUBLE PRECISION
        !                        A RIGHT SIDE Y WAS COMPUTED.
        !                        USING THIS ARRAY Y, SUBROUTINE genbun
        !                        WAS CALLED TO PRODUCE APPROXIMATE
        !                        SOLUTION Z.  THEN RELATIVE ERROR
        !                          E = MAX(abs(Z(I, J)-X(I, J)))/
        !                              MAX(abs(X(I, J)))
        !                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
        !                        OVER I=1, 2, ..., M AND J=1, ..., N.
        !
        !                        THE VALUE OF E IS GIVEN IN THE TABLE
        !                        BELOW FOR SOME TYPICAL VALUES OF M AND N.
        !
        !                   M (=N)    MPEROD    NPEROD        E
        !                   ------    ------    ------      ------
        !
        !                     31        0         0         6.E-14
        !                     31        1         1         4.E-13
        !                     31        1         3         3.E-13
        !                     32        0         0         9.E-14
        !                     32        1         1         3.E-13
        !                     32        1         3         1.E-13
        !                     33        0         0         9.E-14
        !                     33        1         1         4.E-13
        !                     33        1         3         1.E-13
        !                     63        0         0         1.E-13
        !                     63        1         1         1.E-12
        !                     63        1         3         2.E-13
        !                     64        0         0         1.E-13
        !                     64        1         1         1.E-12
        !                     64        1         3         6.E-13
        !                     65        0         0         2.E-13
        !                     65        1         1         1.E-12
        !                     65        1         3         4.E-13
        !
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
        real (wp),             intent (in out) :: y(idimy,*)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)             :: irwk
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------

        ! initialize error flag
        ierror = 0

        ! check input arguments: case 1
        if (m <= 2) then
            ierror = 1
            return
        end if

        ! check input arguments: case 2
        if (n <= 2) then
            ierror = 2
            return
        end if

        ! check input arguments: case 3
        if (idimy < m) then
            ierror = 3
            return
        end if

        ! check input arguments: case 4
        if (nperod < 0 .or. nperod > 4) then
            ierror = 4
            return
        end if

        ! check input arguments: case 5
        if (mperod < 0 .or. mperod > 1) then
            ierror = 5
            return
        end if

        ! compute real workspace size for genbun
        call workspace%get_genbun_workspace_dimensions(n, m, irwk)

        ! create real workspace
        associate( icwk => 0 )
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! solve system
        associate( rew => workspace%real_workspace )
            call genbunn(nperod, n, mperod, m, a, b, c, idimy, y, ierror, rew )
        end associate

        ! Release allocated work space
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
        real (wp),    intent (in)     :: a(*)
        real (wp),    intent (in)     :: b(*)
        real (wp),    intent (in)     :: c(*)
        real (wp),    intent (in out) :: y(idimy, *)
        real (wp),    intent (in out) :: w(*)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)    :: i, mp1, iwba, iwbb, iwbc, iwb2, iwb3, iww1, iww2, iww3
        integer (ip)    :: iwd, iwtcos, iwp, k, j, mp, np, irev, mh, mhm1, modd
        integer (ip)    :: nby2, mskip
        real (wp)       :: a1, ipstor
        !-----------------------------------------------

        ! initialize error flag
        ierror = 0

        ! check input arguments: case 6
        if (mperod /= 1) then
            loop: do i = 2, m
                if (a(i) /= c(1)) then
                    ierror = 6
                    exit loop
                end if
                if (c(i) /= c(1)) then
                    ierror = 6
                    exit loop
                end if
                if (b(i) /= b(1)) then
                    ierror = 6
                    exit loop
                end if
            end do loop
        end if

        ! check input arguments: case 7
        if (a(1) /= 0.0_wp .or. c(m) /= 0.0_wp) then
            ierror = 7
        end if

        ! check the error flag
        if (ierror /= 0) then
            return
        end if

        mp1 = m + 1
        iwba = mp1
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

        w(iwba:m-1+iwba) = -a(:m)
        w(iwbc:m-1+iwbc) = -c(:m)
        w(iwbb:m-1+iwbb) = 2. - b(:m)
        y(:m, :n) = -y(:m, :n)

        mp = mperod + 1
        np = nperod + 1

        select case (mp)
            case (1)
                go to 114
            case (2)
                go to 107
        end select

107 continue

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
        case (5)
            go to 123
    end select

112 continue

    ipstor = w(iww1)
    irev = 2

    if (nperod == 4) then
        go to 124
    end if

113 continue

    select case (mp)
        case (1)
            go to 127
        case (2)
            w(1) = ipstor + real(iwp - 1, kind=wp)
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

    go to 107
!
!     reverse columns when nperod = 4.
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

    select case (irev)
        case (1)
            call poisn2(m, n, 1, 2, w(iwba), w(iwbb), w(iwbc), &
                y, idimy, w, w(iwb2), w(iwb3), w(iww1), w(iww2),  &
                w(iww3), w(iwd), w(iwtcos), w(iwp))
            go to 112
        case (2)
            go to 113
    end select

127 continue

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
    real (wp),    intent (in out) :: q(idimq, 1)
    real (wp),    intent (in out) :: b(*)
    real (wp),    intent (in out) :: w(*)
    real (wp),    intent (in out) :: d(*)
    real (wp),    intent (in out) :: tcos(*)
    real (wp),    intent (in out) :: p(*)
    !-----------------------------------------------
    ! dictionary: local variables
    !-----------------------------------------------
    integer (ip)    :: m, n, jsh, ipp, ipstor, kr, irreg, jstsav, i, lr, nun
    integer (ip)    :: jst, jsp, l, nodd, j, jm1, jp1, jm2, jp2, jm3, jp3, noddpr
    integer (ip)    :: krpi, ideg, jdeg
    real (wp)       :: fi, t
    !-----------------------------------------------

    m = mr
    n = nr
    jsh = 0
    fi = 1.0_wp/istag
    ipp = -m
    ipstor = 0

    select case (istag)
        case default
            kr = 0
            irreg = 1
            if (n > 1) then
                go to 106
            end if
            tcos(1) = 0.0_wp
        case (2)
            kr = 1
            jstsav = 1
            irreg = 2
            if (n > 1) then
                go to 106
            end if
            tcos(1) = -1.0_wp
    end select

103 continue

    b(:m) = q(:m, 1)
    call trix (1, 0, m, ba, bb, bc, b, tcos, d, w)
    q(:m, 1) = b(:m)

    w(1) = ipstor
    return

106 continue

    lr = 0
    p(:m) = 0.0_wp
    nun = n
    jst = 1
    jsp = n
!
!==> irreg = 1 when no irregularities have occurred, otherwise it is 2.
!
108 continue

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

111 continue

    call cosgen(jst, 1, 0.5_wp, 0.0_wp, tcos)

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

            call trix (jst, 0, m, ba, bb, bc, b, tcos, d, w)
            q(:m, j) = q(:m, j) + b(:m)

        end do
    end if
    !
    !==> reduction for last unknown
    !
    select case (nodd)
        case default
            select case (irreg)
                case (1)
                    go to 152
                case (2)
                    go to 120
            end select
        !
        !==> odd number of unknowns
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
                go to 123
            case (2)
                go to 121
        end select

121 continue

    if (jst /= 1) then
        go to 123
    end if

    b(:m) = q(:m, j)
    q(:m, j) = 0.0_wp
    go to 130

123 continue

    select case (noddpr)
        case (2)
            b(:m) = 0.5_wp*(q(:m, jm2)-q(:m, jm1)-q(:m, jm3)) &
                + q(:m, jp2) - q(:m, jp1) + q(:m, j)
        case default
            b(:m) = 0.5_wp*(q(:m, jm2)-q(:m, jm1)&
                -q(:m, jm3)) + p(ipp+1:m+ipp) + q(:m, j)
    end select

128 continue

    q(:m, j) = 0.5*(q(:m, j)-q(:m, jm1)-q(:m, jp1))

130 continue

    call trix (jst, 0, m, ba, bb, bc, b, tcos, d, w)
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
        call cosgen(lr, jstsav, 0.0_wp, fi, tcos(jst+1))
        call merge_rename(tcos, 0, jst, jst, lr, kr)
    end if

    call cosgen(kr, jstsav, 0.0_wp, fi, tcos)
    call trix (kr, kr, m, ba, bb, bc, b, tcos, d, w)
    q(:m, j) = q(:m, jm2) + b(:m) + p(ipp+1:m+ipp)
    lr = kr
    kr = kr + l
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
        case default
            jstsav = jst
            ideg = jst
            kr = l
        case (2)
            call cosgen(kr, jstsav, 0.0, fi, tcos)
            call cosgen(lr, jstsav, 0.0, fi, tcos(kr+1))
            ideg = kr
            kr = kr + jst
    end select

139 continue

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

150 continue

    call trix (ideg, lr, m, ba, bb, bc, b, tcos, d, w)
    q(:m, j) = q(:m, j) + b(:m)

end select

152 continue

    nun = nun/2
    noddpr = nodd
    jsh = jst
    jst = 2*jst

    if (nun >= 2) then
        go to 108
    end if
    !
    !==> start solution.
    !
    j = jsp
    b(:m) = q(:m, j)

    select case (irreg)
        case (2)
            kr = lr + jst
            call cosgen(kr, jstsav, 0.0_wp, fi, tcos)
            call cosgen(lr, jstsav, 0.0_wp, fi, tcos(kr+1))
            ideg = kr
        case default
            call cosgen(jst, 1, 0.5_wp, 0.0_wp, tcos)
            ideg = jst
    end select

156 continue

    call trix(ideg, lr, m, ba, bb, bc, b, tcos, d, w)
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
                go to 168
            end if
            b(:m) = q(:m, j) + q(:m, jm2)

            if (jst < jstsav) then
                irreg = 1
            end if

            select case (irreg)
                case (1)
                    go to 170
                case (2)
                    go to 171
            end select
168     continue
        b(:m) = q(:m, j) + q(:m, jm2) + q(:m, jp2)
    end if

170 continue

    call cosgen(jst, 1, 0.5_wp, 0.0_wp, tcos)
    ideg = jst
    jdeg = 0
    go to 172

171 continue

    if (j + l > n) then
        lr = lr - jst
    end if

    kr = jst + lr
    call cosgen(kr, jstsav, 0.0_wp, fi, tcos)
    call cosgen(lr, jstsav, 0.0_wp, fi, tcos(kr+1))
    ideg = kr
    jdeg = lr

172 continue

    call trix (ideg, jdeg, m, ba, bb, bc, b, tcos, d, w)

    if (jst <= 1) then
        q(:m, j) = b(:m)
    else
        if (jp2 > n) then
            go to 177
        end if

175 continue

    q(:m, j) = 0.5*(q(:m, j)-q(:m, jm1)-q(:m, jp1)) + b(:m)
    cycle

177 continue

    select case (irreg)
        case (1)
            go to 175
        case (2)
            go to 178
    end select

178 continue

    if (j + jsh <= n) then
        q(:m, j) = b(:m) + p(ipp+1:m+ipp)
        ipp = ipp - m
    else
        q(:m, j) = b(:m) + q(:m, j) - q(:m, jm1)
    end if

end if
end do

l = l/2
go to 164
w(1) = ipstor

end subroutine poisd2


subroutine poisn2(m, n, istag, mixbnd, a, bb, c, q, idimq, b, b2, &
    b3, w, w2, w3, d, tcos, p)
    !-----------------------------------------------
    ! dictionary: calling arguments
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
    ! dictionary: local variables
    !-----------------------------------------------
    integer (ip)    :: k(4)
    integer (ip)    :: k1, k2, k3, k4, mr, ipp, ipstor, i2r, jr, nr, nlast, kr
    integer (ip)    :: lr, i, nrod, jstart, jstop, i2rby2, j, jp1, jp2, jp3, jm1
    integer (ip)    :: jm2, jm3, nrodpr, ii, i1, i2, jr2, nlastp, jstep
    real (wp)       :: fistag, fnum, fden, fi, t
    !-----------------------------------------------
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
    equivalence (k(1), k1), (k(2), k2), (k(3), k3), (k(4), k4)

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
            go to 101
        case (2)
            go to 103
    end select

101 continue

    q(:mr, n) = 0.5_wp * q(:mr, n)

    select case (mixbnd)
        case (1)
            go to 103
        case (2)
            go to 104
    end select

103 continue

    if (n <= 3) then
        go to 155
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

107 continue

    jstop = nlast - jr

    if (nrod == 0) then
        jstop = jstop - i2r
    end if

    call cosgen(i2r, 1, 0.5_wp, 0.0_wp, tcos)

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

            call trix (i2r, 0, mr, a, bb, c, b, tcos, d, w)

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
                call cosgen(lr, 1, 0.5_wp, fden, tcos(kr+1))
            else
                b(:mr) = fistag*b(:mr)
            end if
        end if
        call cosgen(kr, 1, 0.5_wp, fden, tcos)
        call trix (kr, lr, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
        kr = kr + i2r
    else
        jp1 = j + i2rby2
        jp2 = j + i2r
        if (i2r == 1) then
            b(:mr) = q(:mr, j)
            call trix (1, 0, mr, a, bb, c, b, tcos, d, w)
            ipp = 0
            ipstor = mr
            select case (istag)
                case default
                    p(:mr) = b(:mr)
                    b(:mr) = b(:mr) + q(:mr, n)
                    tcos(1) = 1.0_wp
                    tcos(2) = 0.0_wp
                    call trix (1, 1, mr, a, bb, c, b, tcos, d, w)
                    q(:mr, j) = q(:mr, jm2) + p(:mr) + b(:mr)
                    go to 150
                case (1)
                    p(:mr) = b(:mr)
                    q(:mr, j) = q(:mr, jm2) + 2.0_wp*q(:mr, jp2) + 3.0_wp*b(:mr)
                    go to 150
            end select
        end if

        b(:mr) = q(:mr, j) + 0.5*(q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))

        if (nrodpr == 0) then
            b(:mr) = b(:mr) + p(ipp+1:mr+ipp)
        else
            b(:mr) = b(:mr) + q(:mr, jp2) - q(:mr, jp1)
        end if

        call trix (i2r, 0, mr, a, bb, c, b, tcos, d, w)
        ipp = ipp + mr
        ipstor = max(ipstor, ipp + mr)
        p(ipp+1:mr+ipp) = b(:mr) + 0.5_wp*(q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
        b(:mr) = p(ipp+1:mr+ipp) + q(:mr, jp2)

        if (lr /= 0) then
            call cosgen(lr, 1, 0.5_wp, fden, tcos(i2r+1))
            call merge_rename(tcos, 0, i2r, i2r, lr, kr)
        else
            do i = 1, i2r
                ii = kr + i
                tcos(ii) = tcos(i)
            end do
        end if

        call cosgen(kr, 1, 0.5_wp, fden, tcos)

        if (lr == 0) then
            select case (istag)
                case (1)
                    go to 146
                case (2)
                    go to 145
            end select
        end if

145 continue
    call trix (kr, kr, mr, a, bb, c, b, tcos, d, w)
    go to 148
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
            go to 155
        end if
    case (2)

        nr = nlast/jr

        if (nr <= 1) then
            go to 192
        end if
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
        if (lr /= 0) then
            go to 170
        end if

        if (n == 3) then
            !
            !     case n = 3.
            !
            select case (istag)
                case (1)
                    go to 156
                case (2)
                    go to 168
            end select

156     continue

        b(:mr) = q(:mr, 2)
        tcos(1) = 0.0_wp

        call trix (1, 0, mr, a, bb, c, b, tcos, d, w)

        q(:mr, 2) = b(:mr)
        b(:mr) = 4.0_wp*b(:mr) + q(:mr, 1) + 2.0_wp*q(:mr, 3)
        tcos(1) = -2.0_wp
        tcos(2) = 2.0_wp
        i1 = 2
        i2 = 0

        call trix (i1, i2, mr, a, bb, c, b, tcos, d, w)

        q(:mr, 2) = q(:mr, 2) + b(:mr)
        b(:mr) = q(:mr, 1) + 2.0_wp*q(:mr, 2)
        tcos(1) = 0.0_wp

        call trix (1, 0, mr, a, bb, c, b, tcos, d, w)

        q(:mr, 1) = b(:mr)
        jr = 1
        i2r = 0
        go to 194
    end if
    !
    !==> case n = 2**p+1
    !
    select case (istag)
        case (1)
            go to 162
        case (2)
            go to 170
    end select

162 continue

    b(:mr) = &
        q(:mr, j) + 0.5_wp*q(:mr, 1) &
        - q(:mr, jm1) + q(:mr, nlast) - &
        q(:mr, jm2)

    call cosgen(jr, 1, 0.5_wp, 0.0_wp, tcos)

    call trix (jr, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, j) = 0.5_wp*(q(:mr, j)-q(:mr, jm1)-q(:mr, jp1)) + b(:mr)

    b(:mr) = q(:mr, 1) + 2.0_wp*q(:mr, nlast) + 4.0_wp*q(:mr, j)

    jr2 = 2*jr

    call cosgen(jr, 1, 0.0_wp, 0.0_wp, tcos)

    tcos(jr+1:jr*2) = -tcos(jr:1:(-1))

    call trix (jr2, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, j) = q(:mr, j) + b(:mr)
    b(:mr) = q(:mr, 1) + 2.0_wp*q(:mr, j)

    call cosgen(jr, 1, 0.5_wp, 0.0_wp, tcos)

    call trix (jr, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, 1) = 0.5_wp*q(:mr, 1) - q(:mr, jm1) + b(:mr)

    go to 194
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

    go to 177

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

    call cosgen(k2 + istag - 2, 1, 0.0_wp, fnum, tcos(k4))

    k4 = k1 + k2 + 1

    call cosgen(jr - 1, 1, 0.0_wp, 1.0_wp, tcos(k4))
    call merge_rename(tcos, k1, k2, k1 + k2, jr - 1, 0)

    k3 = k1 + k2 + lr

    call cosgen(jr, 1, 0.5, 0.0, tcos(k3+1))

    k4 = k3 + jr + 1

    call cosgen(kr, 1, 0.5, fden, tcos(k4))
    call merge_rename(tcos, k3, jr, k3 + jr, kr, k1)

    if (lr /= 0) then
        call cosgen(lr, 1, 0.5, fden, tcos(k4))
        call merge_rename(tcos, k3, jr, k3 + jr, lr, k3 - lr)
        call cosgen(kr, 1, 0.5, fden, tcos(k4))
    end if

    k3 = kr
    k4 = kr

    call tri3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)

    b(:mr) = b(:mr) + b2(:mr) + b3(:mr)
    tcos(1) = 2.0_wp

    call trix (1, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, j) = q(:mr, j) + b(:mr)
    b(:mr) = q(:mr, 1) + 2.0_wp*q(:mr, j)

    call cosgen(jr, 1, 0.5_wp, 0.0_wp, tcos)
    call trix (jr, 0, mr, a, bb, c, b, tcos, d, w)

    if (jr == 1) then
        q(:mr, 1) = b(:mr)
        go to 194
    end if

    q(:mr, 1) = 0.5_wp*q(:mr, 1) - q(:mr, jm1) + b(:mr)
    go to 194
end if

if (n == 2) then
    !
    !==> case  n = 2
    !
    b(:mr) = q(:mr, 1)
    tcos(1) = 0.0_wp

    call trix (1, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, 1) = b(:mr)
    b(:mr) = 2.*(q(:mr, 2)+b(:mr))*fistag
    tcos(1) = -fistag
    tcos(2) = 2.0_wp

    call trix (2, 0, mr, a, bb, c, b, tcos, d, w)

    q(:mr, 1) = q(:mr, 1) + b(:mr)
    jr = 1
    i2r = 0
    go to 194
end if

b3(:mr) = 0.0_wp
b(:mr) = q(:mr, 1) + 2.*p(ipp+1:mr+ipp)
q(:mr, 1) = 0.5*q(:mr, 1) - q(:mr, jm1)
b2(:mr) = 2.*(q(:mr, 1)+q(:mr, nlast))
k1 = kr + jr - 1
tcos(k1+1) = -2.
k4 = k1 + 3 - istag

call cosgen(kr + istag - 2, 1, 0.0, fnum, tcos(k4))

k4 = k1 + kr + 1

call cosgen(jr - 1, 1, 0.0, 1.0, tcos(k4))
call merge_rename(tcos, k1, kr, k1 + kr, jr - 1, 0)
call cosgen(kr, 1, 0.5, fden, tcos(k1+1))

k2 = kr
k4 = k1 + k2 + 1

call cosgen(lr, 1, 0.5, fden, tcos(k4))

k3 = lr
k4 = 0

call tri3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)

b(:mr) = b(:mr) + b2(:mr)
tcos(1) = 2.0_wp

call trix (1, 0, mr, a, bb, c, b, tcos, d, w)

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
        q(:mr, nlast) = 0.0_wp
    else
        if (nrod == 0) then
            q(:mr, nlast) = p(ipp+1:mr+ipp)
            ipp = ipp - mr
        else
            q(:mr, nlast) = q(:mr, nlast) - q(:mr, jm2)
        end if
    end if
    call cosgen(kr, 1, 0.5, fden, tcos)
    call cosgen(lr, 1, 0.5, fden, tcos(kr+1))

    if (lr == 0) then
        b(:mr) = fistag*b(:mr)
    end if

    call trix (kr, lr, mr, a, bb, c, b, tcos, d, w)

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
    call cosgen(jr, 1, 0.5_wp, 0.0_wp, tcos)

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

        call trix (jr, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
    end do

    nrod = 1
    if (nlast + i2r <= n) then
        nrod = 0
    end if

    if (nlastp /= nlast) then
        go to 194
    end if

    go to 206

    w(1) = ipstor

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
    real (wp),    intent (in out) :: q(idimq, 1)
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
    ! dictionary: local variables
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
