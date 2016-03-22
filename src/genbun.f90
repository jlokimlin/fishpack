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
    public :: test_genbun
    public :: genbun
    public :: genbunn

contains


    subroutine test_genbun()
        !     file tgenbun.f
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
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: idimf, m, mp1, mperod, n, nperod, i, j, ierror
        real (wp) , dimension(22, 40) :: f
        real (wp), dimension(20) :: a, b, c
        real (wp), dimension(21) :: x
        real (wp), dimension(41) :: y
        real (wp)       :: dx, pi, dy, s, t, tsq, t4, discretization_error
        !-----------------------------------------------
        !
        !     FROM THE DIMENSION STATEMENT WE GET THAT IDIMF = 22
        !
        idimf = 22
        m = 20
        mp1 = m + 1
        mperod = 1
        dx = 0.05
        n = 40
        nperod = 0
        pi = acos( -1.0 )
        dy = pi/20.
        !
        !     GENERATE GRID POINTS FOR LATER USE.
        !
        do i = 1, mp1
            x(i) = real(i - 1)*dx
        end do
        do j = 1, n
            y(j) = (-pi) + real(j - 1)*dy
        end do
        !
        !     GENERATE COEFFICIENTS.
        !
        s = (dy/dx)**2
        do i = 2, 19
            t = 1. + X(i)
            tsq = t**2
            a(i) = (tsq + t*dx)*s
            b(i) = -2.*tsq*s
            c(i) = (tsq - t*dx)*s
        end do
        a(1) = 0.
        b(1) = -2.*s
        c(1) = -B(1)
        b(20) = -2.*s*(1. + X(20))**2
        a(20) = (-B(20)/2.) + (1. + X(20))*dx*s
        c(20) = 0.
        !
        !     GENERATE RIGHT SIDE.
        !
        do i = 2, 19
            do j = 1, n
                f(i, j) = 3.*(1. + X(i))**4*dy**2*SIN(Y(j))
            end do
        end do
        t = 1. + X(20)
        tsq = t**2
        t4 = tsq**2
        do j = 1, n
            f(1, j) = (11. + 8./dx)*dy**2*SIN(Y(j))
            f(20, j) = (3.*t4*dy**2 - 16.*tsq*s + 16.*t*s*dx)*SIN(Y(j))
        end do

        call genbun(nperod, n, mperod, m, a, b, c, idimf, f, ierror)
        !
        !     compute discretization error.  the exact solution is
        !
        !            u(x, y) = (1+x)**4*sin(y) .
        !
        discretization_error = 0.0_wp
        do i = 1, m
            do j = 1, n
                t = abs(f(i,j)-(1.+x(i))**4*sin(y(j)))
                discretization_error = max(t, discretization_error)
            end do
        end do
        !     Print earlier output from platforms with 32 and 64 bit floating point
        !     arithemtic followed by the output from this computer
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     genbun *** TEST RUN *** '
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') '     ierror = 0,  discretization error = 9.6406E-3'
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)') &
            '     ierror =', ierror, ' discretization error = ', discretization_error

    end subroutine test_genbun


    subroutine genbun( nperod, n, mperod, m, a, b, c, idimy, y, ierror)
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
        integer (ip)             :: irwk, icwk
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
        integer (ip)    :: iwd, iwtcos, iwp, k, j, mp, np, ipstor, irev, mh, mhm1, modd
        integer (ip)    :: mhpi, mhmi, nby2, mskip
        real (wp)       :: a1
        !-----------------------------------------------

        ! initialize error flag
        ierror = 0

        ! check input arguments: case 6
        if (mperod /= 1) then
            do i = 2, m
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
        end if

        ! check the error flag
        if (ierror /= 0) return

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
        go to (114, 107) mp
107 continue
    go to (108, 109, 110, 111, 123) np
108 continue
    call poisp2 (m, n, w(iwba), w(iwbb), w(iwbc), y, idimy, w, w(iwb2) &
        , w(iwb3), w(iww1), w(iww2), w(iww3), w(iwd), w(iwtcos), w(iwp) &
        )
    go to 112
109 continue
    call poisd2 (m, n, 1, w(iwba), w(iwbb), w(iwbc), y, idimy, w, w( &
        iww1), w(iwd), w(iwtcos), w(iwp))
    go to 112
110 continue
    call poisn2 (m, n, 1, 2, w(iwba), w(iwbb), w(iwbc), y, idimy, w, w &
        (iwb2), w(iwb3), w(iww1), w(iww2), w(iww3), w(iwd), w(iwtcos), &
        w(iwp))
    go to 112
111 continue
    call poisn2 (m, n, 1, 1, w(iwba), w(iwbb), w(iwbc), y, idimy, w, w &
        (iwb2), w(iwb3), w(iww1), w(iww2), w(iww3), w(iwd), w(iwtcos), &
        w(iwp))
112 continue
    ipstor = w(iww1)
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
        w(:mhm1) = y(mh-1:mh-mhm1:(-1), j) - y(mh+1:mhm1+mh, j)
        w(mh+1:mhm1+mh) = y(mh-1:mh-mhm1:(-1), j) + y(mh+1:mhm1+mh, j)
        w(mh) = 2.*y(mh, j)
        go to (117, 116) modd
116 continue
    w(m) = 2.*y(m, j)
117 continue
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
122 continue
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
    go to (110, 113) irev
127 continue
    do j = 1, n
        w(mh-1:mh-mhm1:(-1)) = 0.5*(y(mh+1:mhm1+mh, j)+y(:mhm1, j))
        w(mh+1:mhm1+mh) = 0.5*(y(mh+1:mhm1+mh, j)-y(:mhm1, j))
        w(mh) = 0.5*y(mh, j)
        go to (130, 129) modd
129 continue
    w(m) = 0.5*y(m, j)
130 continue
    y(:m, j) = w(:m)
end do
133 continue
    w(1) = ipstor + iwp - 1

end subroutine genbunn



subroutine POISD2(mr, nr, istag, ba, bb, bc, q, idimq, b, w, d, tcos, p)
    !--------------------------------------------------------------------------------
    ! Dictionary: calling arguments
    !--------------------------------------------------------------------------------
    integer (ip), intent (in) :: mr
    integer (ip), intent (in) :: nr
    integer (ip), intent (in) :: istag
    integer (ip), intent (in) :: idimq
    real (wp) :: ba(*)
    real (wp) :: bb(*)
    real (wp) :: bc(*)
    real (wp), intent (in out) :: q(idimq, 1)
    real (wp) :: b(*)
    real (wp), intent (in out) :: w(*)
    real (wp) :: d(*)
    real (wp) :: tcos(*)
    real (wp), intent (in out) :: p(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip)    :: m, n, jsh, ipp, ipstor, kr, irreg, jstsav, i, lr, nun
    integer (ip)    :: jst, jsp, l, nodd, j, jm1, jp1, jm2, jp2, jm3, jp3, noddpr, ip1
    integer (ip)    :: krpi, ideg, jdeg
    real (wp)       :: fi, t
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
    jsh = 0
    fi = 1./real(istag)
    ipp = -m
    ipstor = 0
    select case (istag)
        case default
            kr = 0
            irreg = 1
            if (n > 1) go to 106
            tcos(1) = 0.
        case (2)
            kr = 1
            jstsav = 1
            irreg = 2
            if (n > 1) go to 106
            tcos(1) = -1.
    end select
103 continue
    b(:m) = Q(:m, 1)
    call TRIX (1, 0, m, ba, bb, bc, b, tcos, d, w)
    q(:m, 1) = B(:m)
    go to 183
106 continue
    lr = 0
    p(:m) = 0.
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
    call cosgen(jst, 1, 0.5, 0.0, tcos)
    if (l <= jsp) then
        do j = l, jsp, l
            jm1 = j - jsh
            jp1 = j + jsh
            jm2 = j - jst
            jp2 = j + jst
            jm3 = jm2 - jsh
            jp3 = jp2 + jsh
            if (jst == 1) then
                b(:m) = 2.*Q(:m, j)
                q(:m, j) = Q(:m, jm2) + Q(:m, jp2)
            else
                do i = 1, m
                    t = Q(i, j) - Q(i, jm1) - Q(i, jp1) + Q(i, jm2) + Q(i, jp2)
                    b(i) = t + Q(i, j) - Q(i, jm3) - Q(i, jp3)
                    q(i, j) = t
                end do
            end if
            call TRIX (jst, 0, m, ba, bb, bc, b, tcos, d, w)
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
    b(:m) = Q(:m, j)
    q(:m, j) = 0.
    go to 130
123 continue
    select case (noddpr)
        case default
            b(:m) = 0.5*(Q(:m, jm2)-Q(:m, jm1)-Q(:m, jm3)) + P(ipp+1:m+ipp) &
                + Q(:m, j)
        case (2)
            b(:m) = 0.5*(Q(:m, jm2)-Q(:m, jm1)-Q(:m, jm3)) + Q(:m, jp2) - Q( &
                :m, jp1) + Q(:m, j)
    end select
128 continue
    q(:m, j) = 0.5*(Q(:m, j)-Q(:m, jm1)-Q(:m, jp1))
130 continue
    call TRIX (jst, 0, m, ba, bb, bc, b, tcos, d, w)
    ipp = ipp + m
    ipstor = max(ipstor, ipp + m)
    p(ipp+1:m+ipp) = Q(:m, j) + B(:m)
    b(:m) = Q(:m, jp2) + P(ipp+1:m+ipp)
    if (lr == 0) then
        do i = 1, jst
            krpi = kr + i
            tcos(krpi) = TCOS(i)
        end do
    else
        call cosgen(lr, jstsav, 0., fi, TCOS(jst+1))
        call merge_rename(tcos, 0, jst, jst, lr, kr)
    end if
    call cosgen(kr, jstsav, 0.0, fi, tcos)
    call TRIX (kr, kr, m, ba, bb, bc, b, tcos, d, w)
    q(:m, j) = Q(:m, jm2) + B(:m) + P(ipp+1:m+ipp)
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
            call cosgen(kr, jstsav, 0.0, fi, tcos)
            call cosgen(lr, jstsav, 0.0, fi, TCOS(kr+1))
            ideg = kr
            kr = kr + jst
    end select
139 continue
    if (jst == 1) then
        irreg = 2
        b(:m) = Q(:m, j)
        q(:m, j) = Q(:m, jm2)
    else
        b(:m) = Q(:m, j) + 0.5*(Q(:m, jm2)-Q(:m, jm1)-Q(:m, jm3))
        select case (irreg)
            case default
                q(:m, j) = Q(:m, jm2) + 0.5*(Q(:m, j)-Q(:m, jm1)-Q(:m, jp1))
                irreg = 2
            case (2)
                select case (noddpr)
                    case default
                        q(:m, j) = Q(:m, jm2) + P(ipp+1:m+ipp)
                        ipp = ipp - m
                    case (2)
                        q(:m, j) = Q(:m, jm2) + Q(:m, j) - Q(:m, jm1)
                end select
        end select
    end if
150 continue
    call TRIX (ideg, lr, m, ba, bb, bc, b, tcos, d, w)
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
            call cosgen(jst, 1, 0.5, 0.0, tcos)
            ideg = jst
        case (2)
            kr = lr + jst
            call cosgen(kr, jstsav, 0.0, fi, tcos)
            call cosgen(lr, jstsav, 0.0, fi, TCOS(kr+1))
            ideg = kr
    end select
156 continue
    call TRIX (ideg, lr, m, ba, bb, bc, b, tcos, d, w)
    jm1 = j - jsh
    jp1 = j + jsh
    select case (irreg)
        case default
            q(:m, j) = 0.5*(Q(:m, j)-Q(:m, jm1)-Q(:m, jp1)) + B(:m)
        case (2)
            select case (noddpr)
                case default
                    q(:m, j) = P(ipp+1:m+ipp) + B(:m)
                    ipp = ipp - m
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
    call cosgen(jst, 1, 0.5, 0.0, tcos)
    ideg = jst
    jdeg = 0
    go to 172
171 continue
    if (j + l > n) lr = lr - jst
    kr = jst + lr
    call cosgen(kr, jstsav, 0.0, fi, tcos)
    call cosgen(lr, jstsav, 0.0, fi, TCOS(kr+1))
    ideg = kr
    jdeg = lr
172 continue
    call TRIX (ideg, jdeg, m, ba, bb, bc, b, tcos, d, w)
    if (jst <= 1) then
        q(:m, j) = B(:m)
    else
        if (jp2 > n) go to 177
175 continue
    q(:m, j) = 0.5*(Q(:m, j)-Q(:m, jm1)-Q(:m, jp1)) + B(:m)
    cycle
177 continue
    go to (175, 178) irreg
178 continue
    if (j + jsh <= n) then
        q(:m, j) = B(:m) + P(ipp+1:m+ipp)
        ipp = ipp - m
    else
        q(:m, j) = B(:m) + Q(:m, j) - Q(:m, jm1)
    end if
end if
end do
l = l/2
go to 164
183 continue
    w(1) = ipstor
    return
end subroutine POISD2


subroutine POISN2(m, n, istag, mixbnd, a, bb, c, q, idimq, b, b2, &
    b3, w, w2, w3, d, tcos, p)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in) :: m
    integer (ip), intent (in) :: n
    integer (ip), intent (in) :: istag
    integer (ip), intent (in) :: mixbnd
    integer (ip), intent (in) :: idimq
    real (wp) :: a(*)
    real (wp) :: bb(*)
    real (wp) :: c(*)
    real (wp), intent (in out) :: q(idimq, *)
    real (wp) :: b(*)
    real (wp) :: b2(*)
    real (wp) :: b3(*)
    real (wp), intent (in out) :: w(*)
    real (wp) :: w2(*)
    real (wp) :: w3(*)
    real (wp) :: d(*)
    real (wp) :: tcos(*)
    real (wp), intent (in out) :: p(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip)    :: k(4)
    integer (ip)    :: k1, k2, k3, k4, mr, ipp, ipstor, i2r, jr, nr, nlast, kr
    integer (ip)    :: lr, i, nrod, jstart, jstop, i2rby2, j, jp1, jp2, jp3, jm1
    integer (ip)    :: jm2, jm3, nrodpr, ii, i1, i2, jr2, nlastp, jstep
    real (wp)       :: fistag, fnum, fden, fi, t
    !-----------------------------------------------
    !
    !     SUBROUTINE TO SOLVE POISSON'S EQUATION WITH NEUMANN BOUNDARY
    !     CONDITIONS.
    !
    !     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS A.
    !     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS A-I.
    !     MIXBND = 1 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTH BOUNDARIES.
    !     MIXBND = 2 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTTOM AND
    !     DIRICHLET CONDITION AT TOP.  (FOR THIS CASE, MUST HAVE ISTAG = 1.)
    !
    equivalence (K(1), K1), (K(2), K2), (K(3), K3), (K(4), K4)
    fistag = 3 - istag
    fnum = 1.0_wp/real(istag, kind=wp)
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
    go to (101, 103) istag
101 continue
    q(:mr, n) = 0.5_wp * Q(:mr, n)
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
    call cosgen(i2r, 1, 0.5, 0.0, tcos)
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
                b(:mr) = 2.*Q(:mr, j)
                q(:mr, j) = Q(:mr, jm2) + Q(:mr, jp2)
            else
                do i = 1, mr
                    fi = Q(i, j)
                    q(i, j)=Q(i, j)-Q(i, jm1)-Q(i, jp1)+Q(i, jm2)+Q(i, jp2)
                    b(i) = fi + Q(i, j) - Q(i, jm3) - Q(i, jp3)
                end do
            end if
            call TRIX (i2r, 0, mr, a, bb, c, b, tcos, d, w)
            q(:mr, j) = Q(:mr, j) + B(:mr)
        !
        !     END OF REDUCTION FOR REGULAR UNKNOWNS.
        !
        end do
        !
        !     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
        !
        j = jstop + jr
    end if
    nlast = j
    jm1 = j - i2rby2
    jm2 = j - i2r
    jm3 = jm2 - i2rby2
    if (nrod /= 0) then
        !
        !     ODD NUMBER OF UNKNOWNS
        !
        if (i2r == 1) then
            b(:mr) = fistag*Q(:mr, j)
            q(:mr, j) = Q(:mr, jm2)
        else
            b(:mr) = Q(:mr, j) + 0.5*(Q(:mr, jm2)-Q(:mr, jm1)-Q(:mr, jm3))
            if (nrodpr == 0) then
                q(:mr, j) = Q(:mr, jm2) + P(ipp+1:mr+ipp)
                ipp = ipp - mr
            else
                q(:mr, j) = Q(:mr, j) - Q(:mr, jm1) + Q(:mr, jm2)
            end if
            if (lr /= 0) then
                call cosgen(lr, 1, 0.5, fden, TCOS(kr+1))
            else
                b(:mr) = fistag*B(:mr)
            end if
        end if
        call cosgen(kr, 1, 0.5, fden, tcos)
        call TRIX (kr, lr, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = Q(:mr, j) + B(:mr)
        kr = kr + i2r
    else
        jp1 = j + i2rby2
        jp2 = j + i2r
        if (i2r == 1) then
            b(:mr) = Q(:mr, j)
            call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
            ipp = 0
            ipstor = mr
            select case (istag)
                case default
                    p(:mr) = B(:mr)
                    b(:mr) = B(:mr) + Q(:mr, n)
                    tcos(1) = 1.
                    tcos(2) = 0.
                    call TRIX (1, 1, mr, a, bb, c, b, tcos, d, w)
                    q(:mr, j) = Q(:mr, jm2) + P(:mr) + B(:mr)
                    go to 150
                case (1)
                    p(:mr) = B(:mr)
                    q(:mr, j) = Q(:mr, jm2) + 2.*Q(:mr, jp2) + 3.*B(:mr)
                    go to 150
            end select
        end if
        b(:mr) = Q(:mr, j) + 0.5*(Q(:mr, jm2)-Q(:mr, jm1)-Q(:mr, jm3))
        if (nrodpr == 0) then
            b(:mr) = B(:mr) + P(ipp+1:mr+ipp)
        else
            b(:mr) = B(:mr) + Q(:mr, jp2) - Q(:mr, jp1)
        end if
        call TRIX (i2r, 0, mr, a, bb, c, b, tcos, d, w)
        ipp = ipp + mr
        ipstor = max(ipstor, ipp + mr)
        p(ipp+1:mr+ipp) = B(:mr) + 0.5*(Q(:mr, j)-Q(:mr, jm1)-Q(:mr, jp1))
        b(:mr) = P(ipp+1:mr+ipp) + Q(:mr, jp2)
        if (lr /= 0) then
            call cosgen(lr, 1, 0.5, fden, TCOS(i2r+1))
            call merge_rename(tcos, 0, i2r, i2r, lr, kr)
        else
            do i = 1, i2r
                ii = kr + i
                tcos(ii) = TCOS(i)
            end do
        end if
        call cosgen(kr, 1, 0.5, fden, tcos)
        if (lr == 0) then
            go to (146, 145) istag
        end if
145 continue
    call TRIX (kr, kr, mr, a, bb, c, b, tcos, d, w)
    go to 148
146 continue
    b(:mr) = fistag*B(:mr)
148 continue
    q(:mr, j) = Q(:mr, jm2) + P(ipp+1:mr+ipp) + B(:mr)
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
            !     CASE N = 3.
            !
            go to (156, 168) istag
156     continue
        b(:mr) = Q(:mr, 2)
        tcos(1) = 0.
        call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 2) = B(:mr)
        b(:mr) = 4.*B(:mr) + Q(:mr, 1) + 2.*Q(:mr, 3)
        tcos(1) = -2.
        tcos(2) = 2.
        i1 = 2
        i2 = 0
        call TRIX (i1, i2, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 2) = Q(:mr, 2) + B(:mr)
        b(:mr) = Q(:mr, 1) + 2.*Q(:mr, 2)
        tcos(1) = 0.
        call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, 1) = B(:mr)
        jr = 1
        i2r = 0
        go to 194
    end if
    !
    !     CASE N = 2**P+1
    !
    go to (162, 170) istag
162 continue
    b(:mr) = Q(:mr, j) + 0.5*Q(:mr, 1) - Q(:mr, jm1) + Q(:mr, nlast) - &
        Q(:mr, jm2)
    call cosgen(jr, 1, 0.5, 0.0, tcos)
    call TRIX (jr, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = 0.5*(Q(:mr, j)-Q(:mr, jm1)-Q(:mr, jp1)) + B(:mr)
    b(:mr) = Q(:mr, 1) + 2.*Q(:mr, nlast) + 4.*Q(:mr, j)
    jr2 = 2*jr
    call cosgen(jr, 1, 0.0, 0.0, tcos)
    tcos(jr+1:jr*2) = -TCOS(jr:1:(-1))
    call TRIX (jr2, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = Q(:mr, j) + B(:mr)
    b(:mr) = Q(:mr, 1) + 2.*Q(:mr, j)
    call cosgen(jr, 1, 0.5, 0.0, tcos)
    call TRIX (jr, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = 0.5*Q(:mr, 1) - Q(:mr, jm1) + B(:mr)
    go to 194
!
!     CASE OF GENERAL N WITH NR = 3 .
!
168 continue
    b(:mr) = Q(:mr, 2)
    q(:mr, 2) = 0.
    b2(:mr) = Q(:mr, 3)
    b3(:mr) = Q(:mr, 1)
    jr = 1
    i2r = 0
    j = 2
    go to 177
170 continue
    b(:mr) = 0.5*Q(:mr, 1) - Q(:mr, jm1) + Q(:mr, j)
    if (nrod == 0) then
        b(:mr) = B(:mr) + P(ipp+1:mr+ipp)
    else
        b(:mr) = B(:mr) + Q(:mr, nlast) - Q(:mr, jm2)
    end if
    do i = 1, mr
        t = 0.5*(Q(i, j)-Q(i, jm1)-Q(i, jp1))
        q(i, j) = t
        b2(i) = Q(i, nlast) + t
        b3(i) = Q(i, 1) + 2.*t
    end do
177 continue
    k1 = kr + 2*jr - 1
    k2 = kr + jr
    tcos(k1+1) = -2.
    k4 = k1 + 3 - istag
    call cosgen(k2 + istag - 2, 1, 0.0, fnum, TCOS(k4))
    k4 = k1 + k2 + 1
    call cosgen(jr - 1, 1, 0.0, 1.0, TCOS(k4))
    call merge_rename(tcos, k1, k2, k1 + k2, jr - 1, 0)
    k3 = k1 + k2 + lr
    call cosgen(jr, 1, 0.5, 0.0, TCOS(k3+1))
    k4 = k3 + jr + 1
    call cosgen(kr, 1, 0.5, fden, TCOS(k4))
    call merge_rename(tcos, k3, jr, k3 + jr, kr, k1)
    if (lr /= 0) then
        call cosgen(lr, 1, 0.5, fden, TCOS(k4))
        call merge_rename(tcos, k3, jr, k3 + jr, lr, k3 - lr)
        call cosgen(kr, 1, 0.5, fden, TCOS(k4))
    end if
    k3 = kr
    k4 = kr
    call TRI3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
    b(:mr) = B(:mr) + B2(:mr) + B3(:mr)
    tcos(1) = 2.
    call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = Q(:mr, j) + B(:mr)
    b(:mr) = Q(:mr, 1) + 2.*Q(:mr, j)
    call cosgen(jr, 1, 0.5, 0.0, tcos)
    call TRIX (jr, 0, mr, a, bb, c, b, tcos, d, w)
    if (jr == 1) then
        q(:mr, 1) = B(:mr)
        go to 194
    end if
    q(:mr, 1) = 0.5*Q(:mr, 1) - Q(:mr, jm1) + B(:mr)
    go to 194
end if
if (n == 2) then
    !
    !     CASE  N = 2
    !
    b(:mr) = Q(:mr, 1)
    tcos(1) = 0.
    call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = B(:mr)
    b(:mr) = 2.*(Q(:mr, 2)+B(:mr))*fistag
    tcos(1) = -fistag
    tcos(2) = 2.
    call TRIX (2, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = Q(:mr, 1) + B(:mr)
    jr = 1
    i2r = 0
    go to 194
end if
b3(:mr) = 0.
b(:mr) = Q(:mr, 1) + 2.*P(ipp+1:mr+ipp)
q(:mr, 1) = 0.5*Q(:mr, 1) - Q(:mr, jm1)
b2(:mr) = 2.*(Q(:mr, 1)+Q(:mr, nlast))
k1 = kr + jr - 1
tcos(k1+1) = -2.
k4 = k1 + 3 - istag
call cosgen(kr + istag - 2, 1, 0.0, fnum, TCOS(k4))
k4 = k1 + kr + 1
call cosgen(jr - 1, 1, 0.0, 1.0, TCOS(k4))
call merge_rename(tcos, k1, kr, k1 + kr, jr - 1, 0)
call cosgen(kr, 1, 0.5, fden, TCOS(k1+1))
k2 = kr
k4 = k1 + k2 + 1
call cosgen(lr, 1, 0.5, fden, TCOS(k4))
k3 = lr
k4 = 0
call TRI3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
b(:mr) = B(:mr) + B2(:mr)
tcos(1) = 2.
call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
q(:mr, 1) = Q(:mr, 1) + B(:mr)
go to 194
192 continue
    b(:mr) = Q(:mr, nlast)
    go to 196
194 continue
    j = nlast - jr
    b(:mr) = Q(:mr, nlast) + Q(:mr, j)
196 continue
    jm2 = nlast - i2r
    if (jr == 1) then
        q(:mr, nlast) = 0.
    else
        if (nrod == 0) then
            q(:mr, nlast) = P(ipp+1:mr+ipp)
            ipp = ipp - mr
        else
            q(:mr, nlast) = Q(:mr, nlast) - Q(:mr, jm2)
        end if
    end if
    call cosgen(kr, 1, 0.5, fden, tcos)
    call cosgen(lr, 1, 0.5, fden, TCOS(kr+1))
    if (lr == 0) then
        b(:mr) = fistag*B(:mr)
    end if
    call TRIX (kr, lr, mr, a, bb, c, b, tcos, d, w)
    q(:mr, nlast) = Q(:mr, nlast) + B(:mr)
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
    call cosgen(jr, 1, 0.5, 0.0, tcos)
    do j = jstart, jstop, jstep
        jm2 = j - jr
        jp2 = j + jr
        if (j == jr) then
            b(:mr) = Q(:mr, j) + Q(:mr, jp2)
        else
            b(:mr) = Q(:mr, j) + Q(:mr, jm2) + Q(:mr, jp2)
        end if
        if (jr == 1) then
            q(:mr, j) = 0.
        else
            jm1 = j - i2r
            jp1 = j + i2r
            q(:mr, j) = 0.5*(Q(:mr, j)-Q(:mr, jm1)-Q(:mr, jp1))
        end if
        call TRIX (jr, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = Q(:mr, j) + B(:mr)
    end do
    nrod = 1
    if (nlast + i2r <= n) nrod = 0
    if (nlastp /= nlast) go to 194
    go to 206
222 continue
    w(1) = ipstor

end subroutine POISN2
    !
    !*****************************************************************************************
    !
subroutine POISP2(m, n, a, bb, c, q, idimq, b, b2, b3, w, w2, w3, d, tcos, p)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in) :: m
    integer (ip), intent (in) :: n
    integer (ip) :: idimq
    real (wp) :: a(*)
    real (wp) :: bb(*)
    real (wp) :: c(*)
    real (wp) :: q(idimq, 1)
    real (wp) :: b(*)
    real (wp) :: b2(*)
    real (wp) :: b3(*)
    real (wp), intent (in out) :: w(*)
    real (wp) :: w2(*)
    real (wp) :: w3(*)
    real (wp) :: d(*)
    real (wp) :: tcos(*)
    real (wp) :: p(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: mr, nr, nrm1, j, nrmj, nrpj, i, ipstor, lh
    real (wp)    ::  s, t
    !real, save ::  alll
    !-----------------------------------------------
    !
    !     SUBROUTINE TO SOLVE POISSON EQUATION WITH PERIODIC BOUNDARY
    !     CONDITIONS.
    !
    mr = m
    nr = (n + 1)/2
    nrm1 = nr - 1
    if (2*nr == n) then
        !
        !     EVEN NUMBER OF UNKNOWNS
        !
        do j = 1, nrm1
            nrmj = nr - j
            nrpj = nr + j
            do i = 1, mr
                s = Q(i, nrmj) - Q(i, nrpj)
                t = Q(i, nrmj) + Q(i, nrpj)
                q(i, nrmj) = s
                q(i, nrpj) = t
            end do
        end do
        q(:mr, nr) = 2.*Q(:mr, nr)
        q(:mr, n) = 2.*Q(:mr, n)
        call POISD2 (mr, nrm1, 1, a, bb, c, q, idimq, b, w, d, tcos, p)
        ipstor = W(1)
        call POISN2 (mr, nr + 1, 1, 1, a, bb, c, Q(1, nr), idimq, b, b2 &
            , b3, w, w2, w3, d, tcos, p)
        ipstor = max(ipstor, INT(W(1)))
        do j = 1, nrm1
            nrmj = nr - j
            nrpj = nr + j
            do i = 1, mr
                s = 0.5*(Q(i, nrpj)+Q(i, nrmj))
                t = 0.5*(Q(i, nrpj)-Q(i, nrmj))
                q(i, nrmj) = s
                q(i, nrpj) = t
            end do
        end do
        q(:mr, nr) = 0.5*Q(:mr, nr)
        q(:mr, n) = 0.5*Q(:mr, n)
    else
        do j = 1, nrm1
            nrpj = n + 1 - j
            do i = 1, mr
                s = Q(i, j) - Q(i, nrpj)
                t = Q(i, j) + Q(i, nrpj)
                q(i, j) = s
                q(i, nrpj) = t
            end do
        end do
        q(:mr, nr) = 2.*Q(:mr, nr)
        lh = nrm1/2
        do j = 1, lh
            nrmj = nr - j
            do i = 1, mr
                s = Q(i, j)
                q(i, j) = Q(i, nrmj)
                q(i, nrmj) = s
            end do
        end do
        call POISD2 (mr, nrm1, 2, a, bb, c, q, idimq, b, w, d, tcos, p)
        ipstor = W(1)
        call POISN2 (mr, nr, 2, 1, a, bb, c, Q(1, nr), idimq, b, b2, b3 &
            , w, w2, w3, d, tcos, p)
        ipstor = max(ipstor, INT(W(1)))
        do j = 1, nrm1
            nrpj = nr + j
            do i = 1, mr
                s = 0.5*(Q(i, nrpj)+Q(i, j))
                t = 0.5*(Q(i, nrpj)-Q(i, j))
                q(i, nrpj) = t
                q(i, j) = s
            end do
        end do
        q(:mr, nr) = 0.5*Q(:mr, nr)
        do j = 1, lh
            nrmj = nr - j
            do i = 1, mr
                s = Q(i, j)
                q(i, j) = Q(i, nrmj)
                q(i, nrmj) = s
            end do
        end do
    end if
    w(1) = ipstor
    return

end subroutine POISP2
    !
    !*****************************************************************************************
    !
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
