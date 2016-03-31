module module_poistg

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
    public :: poistg
    public :: poistgg


contains


    subroutine poistg( nperod, n, mperod, m, a, b, c, idimy, y, ierror )
        !
        !     file poistg.f
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
        !     SUBROUTINE poistg (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, ierror)
        !
        !
        ! DIMENSION OF           A(M),  B(M),  C(M),  Y(IDIMY, N)
        ! ARGUMENTS
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
        !                        FOR UNKNOWN X VALUES, WHERE I=1, 2, ..., M
        !                        AND J=1, 2, ..., N
        !
        !                        A(I)*X(I-1, J) + B(I)*X(I, J) + C(I)*X(I+1, J)
        !                        + X(I, J-1) - 2.*X(I, J) + X(I, J+1)
        !                        = Y(I, J)
        !
        !                        THE INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
        !                        I.E. X(0, J) = X(M, J) AND X(M+1, J) = X(1, J), AND
        !                        X(I, 0) MAY BE EQUAL TO X(I, 1) OR -X(I, 1), AND
        !                        X(I, N+1) MAY BE EQUAL TO X(I, N) OR -X(I, N),
        !                        DEPENDING ON AN INPUT PARAMETER.
        !
        ! USAGE                  CALL poistg (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y,
        !                                     ierror)
        !
        ! ARGUMENTS
        !
        ! ON INPUT
        !
        !                        NPEROD
        !                          INDICATES VALUES WHICH X(I, 0) AND X(I, N+1)
        !                          ARE ASSUMED TO HAVE.
        !                          = 1 IF X(I, 0) = -X(I, 1) AND X(I, N+1) = -X(I, N
        !                          = 2 IF X(I, 0) = -X(I, 1) AND X(I, N+1) =  X(I, N
        !                          = 3 IF X(I, 0) =  X(I, 1) AND X(I, N+1) =  X(I, N
        !                          = 4 IF X(I, 0) =  X(I, 1) AND X(I, N+1) = -X(I, N
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
        !                          M MUST BE GREATER THAN 2.
        !
        !                        A, B, C
        !                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
        !                          SPECIFY THE COEFFICIENTS IN THE LINEAR
        !                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0 THE
        !                          ARRAY ELEMENTS MUST NOT DEPEND ON INDEX I,
        !                          BUT MUST BE CONSTANT.  SPECIFICALLY, THE
        !                          SUBROUTINE CHECKS THE FOLLOWING CONDITION
        !                            A(I) = C(1)
        !                            B(I) = B(1)
        !                            C(I) = C(1)
        !                          FOR I = 1, 2, ..., M.
        !
        !                        IDIMY
        !                          THE ROW (OR FIRST) DIMENSION OF THE TWO-
        !                          DIMENSIONAL ARRAY Y AS IT APPEARS IN THE
        !                          PROGRAM CALLING poistg.  THIS PARAMETER IS
        !                          USED TO SPECIFY THE VARIABLE DIMENSION OF Y.
        !                          IDIMY MUST BE AT LEAST M.
        !
        !                        Y
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
        !                          OF EQUATIONS GIVEN ABOVE.
        !                          Y MUST BE DIMENSIONED AT LEAST M X N.
        !
        ! ON OUTPUT
        !
        !                        Y
        !                          CONTAINS THE SOLUTION X.
        !
        !                        ierror
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
        !                          SOLUTION IS NOT ATTEMPTED.
        !                          = 0  NO ERROR
        !                          = 1  IF M .LE. 2
        !                          = 2  IF N .LE. 2
        !                          = 3  IDIMY .LT. M
        !                          = 4  IF NPEROD .LT. 1 OR NPEROD .GT. 4
        !                          = 5  IF MPEROD .LT. 0 OR MPEROD .GT. 1
        !                          = 6  IF MPEROD = 0 AND A(I) .NE. C(1)
        !                               OR B(I) .NE. B(1) OR C(I) .NE. C(1)
        !                               FOR SOME I = 1, 2, ..., M.
        !                          = 7  IF MPEROD .EQ. 1 .AND.
        !                               (A(1).NE.0 .OR. C(M).NE.0)
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
        !
        !                          SINCE THIS IS THE ONLY MEANS OF INDICATING A
        !                          POSSIBLY INCORRECT CALL TO POIS3D, THE USER
        !                          SHOULD TEST ierror AFTER THE CALL.
        !
        !
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED files         fish.f, gnbnaux.f, comf.f
        !
        ! LANGUAGE               FORTRAN 90
        !
        ! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
        !                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
        !                        LIBRARIES IN JANUARY, 1980.
        !                        Revised in June 2004 by John Adams using
        !                        Fortran 90 dynamically allocated work space.
        !
        ! PORTABILITY            FORTRAN 90
        !
        ! ALGORITHM              THIS SUBROUTINE IS AN IMPLEMENTATION OF THE
        !                        ALGORITHM PRESENTED IN THE REFERENCE BELOW.
        !
        ! TIMING                 FOR LARGE M AND N, THE EXECUTION TIME IS
        !                        ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
        !
        ! ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
        !                        UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
        !                        CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
        !                        IN THE 'PURPOSE' SECTION ABOVE, WITH
        !                          A(I) = C(I) = -0.5_wp * B(I) = 1,    I=1, 2, ..., M
        !                        AND, WHEN MPEROD = 1
        !                          A(1) = C(M) = 0
        !                          B(1) = B(M) =-1.
        !
        !                        THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
        !                        SYSTEM AND, USING DOUBLE PRECISION, A RIGHT SID
        !                        Y WAS COMPUTED.  USING THIS ARRAY Y SUBROUTINE
        !                        poistg WAS CALLED TO PRODUCE AN APPROXIMATE
        !                        SOLUTION Z.  THEN THE RELATIVE ERROR, DEFINED A
        !                          E = MAX(abs(Z(I, J)-X(I, J)))/MAX(abs(X(I, J)))
        !                        WHERE THE TWO MAXIMA ARE TAKEN OVER I=1, 2, ..., M
        !                        AND J=1, 2, ..., N, WAS COMPUTED.  VALUES OF E ARE
        !                        GIVEN IN THE TABLE BELOW FOR SOME TYPICAL VALUE
        !                        OF M AND N.
        !
        !                        M (=N)    MPEROD    NPEROD      E
        !                        ------    ------    ------    ------
        !
        !                          31        0-1       1-4     9.E-13
        !                          31        1         1       4.E-13
        !                          31        1         3       3.E-13
        !                          32        0-1       1-4     3.E-12
        !                          32        1         1       3.E-13
        !                          32        1         3       1.E-13
        !                          33        0-1       1-4     1.E-12
        !                          33        1         1       4.E-13
        !                          33        1         3       1.E-13
        !                          63        0-1       1-4     3.E-12
        !                          63        1         1       1.E-12
        !                          63        1         3       2.E-13
        !                          64        0-1       1-4     4.E-12
        !                          64        1         1       1.E-12
        !                          64        1         3       6.E-13
        !                          65        0-1       1-4     2.E-13
        !                          65        1         1       1.E-11
        !                          65        1         3       4.E-13
        !
        ! REFERENCES             SCHUMANN, U. AND R. SWEET, "A DIRECT METHOD
        !                        FOR THE SOLUTION OF POISSON"S EQUATION WITH
        !                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
        !                        GRID OF ARBITRARY SIZE, " J. COMP. PHYS.
        !                        20(1976), PP. 171-182.
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
        real (wp), contiguous, intent (in out) :: y(:,:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)             :: i !! Counter
        type (FishpackWorkspace) :: workspace
        !--------------------------------------------------------------------------------

        ! initialize error flag
        ierror = 0

        ! check input arguments: case 1
        if (m <= 2) then
            ierror = 1
            return
        end if

        ! check input arguments: case 1
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
        if (nperod < 1 .or. nperod > 4) then
            ierror = 4
            return
        end if

        ! check input arguments: case 5
        if (mperod < 0 .or. mperod > 1) then
            ierror = 5
            return
        end if

        ! check input arguments: case 6
        if (mperod /= 1) then
            do i = 1, m
                if (a(i) /= c(1)) then
                    ierror = 6
                    exit
                end if
                if (c(i) /= c(1)) then
                    ierror = 6
                    exit
                end if
                if (b(i) /= b(1)) then
                    ierror = 6
                    exit
                end if
            end do
        end if

       ! check input arguments: case 7
        if (a(1) /= 0.0_wp .or. c(m) /= 0.0_wp) then
            ierror = 7
            return
        end if

        ! final sanity check
        if (ierror /= 0) then
            return
        end if

        ! compute and allocate real work space for poistg
        associate( &
            irwk => m * (9 + int(log(real(n, kind=wp))/log(2.0_wp),kind=ip)) + 4 * n, &
            icwk => 0 &
            )
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! solve system
        associate( rew => workspace%real_workspace )
            call  poistgg( nperod, n, mperod, m, a, b, c, idimy, y, ierror, rew )
        end associate

        ! release work space
        call workspace%destroy()

    end subroutine poistg


    subroutine poistgg(nperod, n, mperod, m, a, b, c, idimy, y, ierror, w)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: nperod
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: mperod
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: idimy
        integer (ip), intent (in)     :: ierror
        real (wp),    intent (in)     :: a(*)
        real (wp),    intent (in)     :: b(*)
        real (wp),    intent (in)     :: c(*)
        real (wp),    intent (in out) :: y(idimy, *)
        real (wp),    intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: iwba, iwbb, iwbc, iwb2, iwb3, iww1, iww2, iww3, iwd
        integer (ip) :: iwtcos, iwp, i, k, j, np, mp, ipstor, irev, mh, mhm1, modd
        integer (ip) :: nby2, mskip
        real (wp)    :: a1
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
        np = nperod
        mp = mperod + 1
        go to (110, 107) mp
107 continue
    go to (108, 108, 108, 119) nperod
108 continue
    call postg2 (np, n, m, w(iwba), w(iwbb), w(iwbc), idimy, y, w, w( &
        iwb2), w(iwb3), w(iww1), w(iww2), w(iww3), w(iwd), w(iwtcos), w &
        (iwp))
    ipstor = w(iww1)
    irev = 2
    if (nperod == 4) go to 120
109 continue
    go to (123, 129) mp
110 continue
    mh = (m + 1)/2
    mhm1 = mh - 1
    modd = 1
    if (mh*2 == m) modd = 2
    do j = 1, n
        do i = 1, mhm1
            w(i) = y(mh-i, j) - y(i+mh, j)
            w(i+mh) = y(mh-i, j) + y(i+mh, j)
        end do
        w(mh) = 2.*y(mh, j)
        go to (113, 112) modd
112 continue
    w(m) = 2.*y(m, j)
113 continue
    y(:m, j) = w(:m)
end do
k = iwbc + mhm1 - 1
i = iwba + mhm1
w(k) = 0.
w(i) = 0.
w(k+1) = 2.*w(k+1)
select case (modd) 
    case default
        k = iwbb + mhm1 - 1
        w(k) = w(k) - w(i-1)
        w(iwbc-1) = w(iwbc-1) + w(iwbb-1)
    case (2)
        w(iwbb-1) = w(k+1)
end select
118 continue
    go to 107
119 continue
    irev = 1
    nby2 = n/2
    np = 2
120 continue
    do j = 1, nby2
        mskip = n + 1 - j
        do i = 1, m
            a1 = y(i, j)
            y(i, j) = y(i, mskip)
            y(i, mskip) = a1
        end do
    end do
    go to (108, 109) irev
123 continue
    do j = 1, n
        w(mh-1:mh-mhm1:(-1)) = 0.5_wp * (y(mh+1:mhm1+mh, j)+y(:mhm1, j))
        w(mh+1:mhm1+mh) = 0.5_wp * (y(mh+1:mhm1+mh, j)-y(:mhm1, j))
        w(mh) = 0.5_wp * y(mh, j)
        go to (126, 125) modd
125 continue
    w(m) = 0.5_wp * y(m, j)
126 continue
    y(:m, j) = w(:m)
end do
129 continue
    w(1) = ipstor + iwp - 1

end subroutine poistgg


subroutine postg2(nperod, n, m, a, bb, c, idimq, q, b, b2, b3, w, &
    w2, w3, d, tcos, p)
    !
    ! Purpose:
    !
    ! Subroutine to solve poisson's equation on a staggered grid.
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: nperod
    integer (ip), intent (in)     :: n
    integer (ip), intent (in)     :: m
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
    integer (ip) :: k(4)
    integer (ip) :: k1, k2, k3, k4, np, mr, ipp, ipstor, i2r, jr, nr, nlast
    integer (ip) :: kr, lr, nrod, jstart, jstop, i2rby2, j, ijump, jp1, jp2, jp3
    integer (ip) :: jm1, jm2, jm3, i, nrodpr, ii, nlastp, jstep
    real (wp)    :: fnum, fnum2, fi, t
    !-----------------------------------------------

    !  specifies that variables share the same memory
    equivalence (k(1), k1), (k(2), k2), (k(3), k3), (k(4), k4)

    np = nperod
    fnum = 0.5_wp * real(np/3)
    fnum2 = 0.5_wp * real(np/2)
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
101 continue
    jr = 2*i2r
    nrod = 1
    if ((nr/2)*2 == nr) nrod = 0
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
                call cosgen(i2r, 1, fnum, 0.5, tcos)
                if (i2r == 1) then
                    b(:mr) = q(:mr, 1)
                    q(:mr, 1) = q(:mr, 2)
                    go to 112
                end if
                b(:mr) = q(:mr, 1) + 0.5_wp * (q(:mr, jp2)-q(:mr, jp1)-q(:mr, &
                    jp3))
                q(:mr, 1) = q(:mr, jp2) + q(:mr, 1) - q(:mr, jp1)
                go to 112
            end if
            go to (107, 108) ijump
107     continue
        ijump = 2
        call cosgen(i2r, 1, 0.5, 0.0, tcos)
108 continue
    if (i2r == 1) then
        b(:mr) = 2.*q(:mr, j)
        q(:mr, j) = q(:mr, jm2) + q(:mr, jp2)
    else
        do i = 1, mr
            fi = q(i, j)
            q(i, j)=q(i, j)-q(i, jm1)-q(i, jp1)+q(i, jm2)+q(i, jp2)
            b(i) = fi + q(i, j) - q(i, jm3) - q(i, jp3)
        end do
    end if
112 continue
    call trix(i2r, 0, mr, a, bb, c, b, tcos, d, w)
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
        b(:mr) = q(:mr, j)
        q(:mr, j) = q(:mr, jm2)
    else
        b(:mr)=q(:mr, j)+0.5_wp * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))
        if (nrodpr == 0) then
            q(:mr, j) = q(:mr, jm2) + p(ipp+1:mr+ipp)
            ipp = ipp - mr
        else
            q(:mr, j) = q(:mr, j) - q(:mr, jm1) + q(:mr, jm2)
        end if
        if (lr /= 0) call cosgen(lr, 1, fnum2, 0.5, tcos(kr+1))
    end if
    call cosgen(kr, 1, fnum2, 0.5, tcos)
    call trix(kr, lr, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = q(:mr, j) + b(:mr)
    kr = kr + i2r
else
    jp1 = j + i2rby2
    jp2 = j + i2r
    if (i2r == 1) then
        b(:mr) = q(:mr, j)
        tcos(1) = 0.
        call trix(1, 0, mr, a, bb, c, b, tcos, d, w)
        ipp = 0
        ipstor = mr
        p(:mr) = b(:mr)
        b(:mr) = b(:mr) + q(:mr, n)
        tcos(1) = (-1.) + 2.*real(np/2)
        tcos(2) = 0.
        call trix(1, 1, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, jm2) + p(:mr) + b(:mr)
    else
        b(:mr)=q(:mr, j)+0.5_wp * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))
        if (nrodpr == 0) then
            b(:mr) = b(:mr) + p(ipp+1:mr+ipp)
        else
            b(:mr) = b(:mr) + q(:mr, jp2) - q(:mr, jp1)
        end if
        call cosgen(i2r, 1, 0.5, 0.0, tcos)
        call trix(i2r, 0, mr, a, bb, c, b, tcos, d, w)
        ipp = ipp + mr
        ipstor = max(ipstor, ipp + mr)
        p(ipp+1:mr+ipp) = b(:mr) + 0.5_wp * (q(:mr, j)-q(:mr, jm1)-q(:mr, &
            jp1))
        b(:mr) = p(ipp+1:mr+ipp) + q(:mr, jp2)
        if (lr /= 0) then
            call cosgen(lr, 1, fnum2, 0.5, tcos(i2r+1))
            call merge_rename(tcos, 0, i2r, i2r, lr, kr)
        else
            do i = 1, i2r
                ii = kr + i
                tcos(ii) = tcos(i)
            end do
        end if
        call cosgen(kr, 1, fnum2, 0.5, tcos)
        call trix(kr, kr, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, jm2) + p(ipp+1:mr+ipp) + b(:mr)
    end if
    lr = kr
    kr = kr + jr
end if
nr = (nlast - 1)/jr + 1
if (nr <= 3) go to 142
i2r = jr
nrodpr = nrod
go to 101
end if
142 continue
    j = 1 + jr
    jm1 = j - i2r
    jp1 = j + i2r
    jm2 = nlast - i2r
    if (nr /= 2) then
        if (lr == 0) then
            if (n == 3) then
                !
                !     case n = 3.
                !
                go to (143, 148, 143) np
143         continue
            b(:mr) = q(:mr, 2)
            b2(:mr) = q(:mr, 1) + q(:mr, 3)
            b3(:mr) = 0.
            select case (np)
                case default
                    tcos(1) = -1.
                    tcos(2) = 1.
                    k1 = 1
                case (1:2)
                    tcos(1) = -2.
                    tcos(2) = 1.
                    tcos(3) = -1.
                    k1 = 2
            end select
147     continue
        k2 = 1
        k3 = 0
        k4 = 0
        go to 150
148 continue
    b(:mr) = q(:mr, 2)
    b2(:mr) = q(:mr, 3)
    b3(:mr) = q(:mr, 1)
    call cosgen(3, 1, 0.5, 0.0, tcos)
    tcos(4) = -1.
    tcos(5) = 1.
    tcos(6) = -1.
    tcos(7) = 1.
    k1 = 3
    k2 = 2
    k3 = 1
    k4 = 1
150 continue
    call tri3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
    b(:mr) = b(:mr) + b2(:mr) + b3(:mr)
    go to (153, 153, 152) np
152 continue
    tcos(1) = 2.
    call trix(1, 0, mr, a, bb, c, b, tcos, d, w)
153 continue
    q(:mr, 2) = b(:mr)
    b(:mr) = q(:mr, 1) + b(:mr)
    tcos(1) = (-1.) + 4.*fnum
    call trix(1, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = b(:mr)
    jr = 1
    i2r = 0
    go to 188
end if
!
!     case n = 2**p+1
!
b(:mr)=q(:mr, j)+q(:mr, 1)-q(:mr, jm1)+q(:mr, nlast)-q(:mr, jm2)
go to (158, 160, 158) np
158 continue
    b2(:mr) = q(:mr, 1) + q(:mr, nlast) + q(:mr, j) - q(:mr, jm1) - &
        q(:mr, jp1)
    b3(:mr) = 0.
    k1 = nlast - 1
    k2 = nlast + jr - 1
    call cosgen(jr - 1, 1, 0.0, 1.0, tcos(nlast))
    tcos(k2) = 2.*real(np - 2)
    call cosgen(jr, 1, 0.5 - fnum, 0.5, tcos(k2+1))
    k3 = (3 - np)/2
    call merge_rename(tcos, k1, jr - k3, k2 - k3, jr + k3, 0)
    k1 = k1 - 1 + k3
    call cosgen(jr, 1, fnum, 0.5, tcos(k1+1))
    k2 = jr
    k3 = 0
    k4 = 0
    go to 162
160 continue
    do i = 1, mr
        fi = (q(i, j)-q(i, jm1)-q(i, jp1))/2.
        b2(i) = q(i, 1) + fi
        b3(i) = q(i, nlast) + fi
    end do
    k1 = nlast + jr - 1
    k2 = k1 + jr - 1
    call cosgen(jr - 1, 1, 0.0, 1.0, tcos(k1+1))
    call cosgen(nlast, 1, 0.5, 0.0, tcos(k2+1))
    call merge_rename(tcos, k1, jr - 1, k2, nlast, 0)
    k3 = k1 + nlast - 1
    k4 = k3 + jr
    call cosgen(jr, 1, 0.5, 0.5, tcos(k3+1))
    call cosgen(jr, 1, 0.0, 0.5, tcos(k4+1))
    call merge_rename(tcos, k3, jr, k4, jr, k1)
    k2 = nlast - 1
    k3 = jr
    k4 = jr
162 continue
    call tri3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
    b(:mr) = b(:mr) + b2(:mr) + b3(:mr)
    if (np == 3) then
        tcos(1) = 2.
        call trix(1, 0, mr, a, bb, c, b, tcos, d, w)
    end if
    q(:mr, j) = b(:mr) + 0.5_wp * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
    b(:mr) = q(:mr, j) + q(:mr, 1)
    call cosgen(jr, 1, fnum, 0.5, tcos)
    call trix(jr, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = q(:mr, 1) - q(:mr, jm1) + b(:mr)
    go to 188
end if
!
!     case of general n with nr = 3 .
!
b(:mr) = q(:mr, 1) - q(:mr, jm1) + q(:mr, j)
if (nrod == 0) then
    b(:mr) = b(:mr) + p(ipp+1:mr+ipp)
else
    b(:mr) = b(:mr) + q(:mr, nlast) - q(:mr, jm2)
end if
do i = 1, mr
    t = 0.5_wp * (q(i, j)-q(i, jm1)-q(i, jp1))
    q(i, j) = t
    b2(i) = q(i, nlast) + t
    b3(i) = q(i, 1) + t
end do
k1 = kr + 2*jr
call cosgen(jr - 1, 1, 0.0, 1.0, tcos(k1+1))
k2 = k1 + jr
tcos(k2) = 2.*real(np - 2)
k4 = (np - 1)*(3 - np)
k3 = k2 + 1 - k4
call cosgen(kr+jr+k4, 1, real(k4)/2., 1.-real(k4), tcos(k3))
k4 = 1 - np/3
call merge_rename(tcos, k1, jr - k4, k2 - k4, kr + jr + k4, 0)
if (np == 3) k1 = k1 - 1
k2 = kr + jr
k4 = k1 + k2
call cosgen(kr, 1, fnum2, 0.5, tcos(k4+1))
k3 = k4 + kr
call cosgen(jr, 1, fnum, 0.5, tcos(k3+1))
call merge_rename(tcos, k4, kr, k3, jr, k1)
k4 = k3 + jr
call cosgen(lr, 1, fnum2, 0.5, tcos(k4+1))
call merge_rename(tcos, k3, jr, k4, lr, k1 + k2)
call cosgen(kr, 1, fnum2, 0.5, tcos(k3+1))
k3 = kr
k4 = kr
call tri3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
b(:mr) = b(:mr) + b2(:mr) + b3(:mr)
if (np == 3) then
    tcos(1) = 2.
    call trix(1, 0, mr, a, bb, c, b, tcos, d, w)
end if
q(:mr, j) = q(:mr, j) + b(:mr)
b(:mr) = q(:mr, 1) + q(:mr, j)
call cosgen(jr, 1, fnum, 0.5, tcos)
call trix(jr, 0, mr, a, bb, c, b, tcos, d, w)
if (jr == 1) then
    q(:mr, 1) = b(:mr)
    go to 188
end if
q(:mr, 1) = q(:mr, 1) - q(:mr, jm1) + b(:mr)
go to 188
end if
b3(:mr) = 0.
b(:mr) = q(:mr, 1) + p(ipp+1:mr+ipp)
q(:mr, 1) = q(:mr, 1) - q(:mr, jm1)
b2(:mr) = q(:mr, 1) + q(:mr, nlast)
k1 = kr + jr
k2 = k1 + jr
call cosgen(jr - 1, 1, 0.0, 1.0, tcos(k1+1))
go to (182, 183, 182) np
182 continue
    tcos(k2) = 2.*real(np - 2)
    call cosgen(kr, 1, 0.0, 1.0, tcos(k2+1))
    go to 184
183 continue
    call cosgen(kr + 1, 1, 0.5, 0.0, tcos(k2))
184 continue
    k4 = 1 - np/3
    call merge_rename(tcos, k1, jr - k4, k2 - k4, kr + k4, 0)
    if (np == 3) k1 = k1 - 1
    k2 = kr
    call cosgen(kr, 1, fnum2, 0.5, tcos(k1+1))
    k4 = k1 + kr
    call cosgen(lr, 1, fnum2, 0.5, tcos(k4+1))
    k3 = lr
    k4 = 0
    call tri3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
    b(:mr) = b(:mr) + b2(:mr)
    if (np == 3) then
        tcos(1) = 2.
        call trix(1, 0, mr, a, bb, c, b, tcos, d, w)
    end if
    q(:mr, 1) = q(:mr, 1) + b(:mr)
188 continue
    j = nlast - jr
    b(:mr) = q(:mr, nlast) + q(:mr, j)
    jm2 = nlast - i2r
    if (jr == 1) then
        q(:mr, nlast) = 0.
    else
        if (nrod == 0) then
            q(:mr, nlast) = p(ipp+1:mr+ipp)
            ipp = ipp - mr
        else
            q(:mr, nlast) = q(:mr, nlast) - q(:mr, jm2)
        end if
    end if
    call cosgen(kr, 1, fnum2, 0.5, tcos)
    call cosgen(lr, 1, fnum2, 0.5, tcos(kr+1))
    call trix(kr, lr, mr, a, bb, c, b, tcos, d, w)
    q(:mr, nlast) = q(:mr, nlast) + b(:mr)
    nlastp = nlast
197 continue
    jstep = jr
    jr = i2r
    i2r = i2r/2
    if (jr == 0) go to 210
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
    call cosgen(jr, 1, 0.5, 0.0, tcos)
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
        call trix(jr, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = q(:mr, j) + b(:mr)
    end do
    nrod = 1
    if (nlast + i2r <= n) nrod = 0
    if (nlastp /= nlast) go to 188
    go to 197
210 continue
    w(1) = ipstor

end subroutine postg2


end module module_poistg
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 Changes
!-----------------------------------------------------------------------
