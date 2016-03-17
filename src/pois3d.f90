module module_pois3d

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_fftpack, only: &
        rffti,&
        rfftf,&
        rfftb,&
        ezffti,&
        ezfftf,&
        ezfftb,&
        sint,&
        sinti, &
        costi,&
        cost,&
        sinqi,&
        sinqf,&
        sinqb,&
        cosqi,&
        cosqf,&
        cosqb,&
        cffti,&
        cfftf, &
        cfftb

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: pois3d
    public :: pois3dd
    public :: pois3d_unit_test

contains

    subroutine pois3d_unit_test()
        !
        !     file tpois3d.f
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
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer::ldimf, mdimf, lperod, l, mperod, m, nperod, n, i, j, k, ierror
        real , dimension(32, 33, 10) :: f
        real , dimension(10) :: a, b, c
        real , dimension(30) :: x, y
        real , dimension(10) :: z
        real :: pi, dx, c1, dy, c2, dz, dzsq, t, err
        !-----------------------------------------------
        !
        !     FROM THE DIMENSION STATEMENT WE GET THAT LDIMF = 32, MDIMF = 33,
        !
        ldimf = 32
        mdimf = 33
        pi = acos( -1.0 )
        lperod = 0
        l = 30
        dx = 2.*pi/real(l)
        c1 = 1./dx**2
        mperod = 0
        m = 30
        dy = 2.*pi/real(m)
        c2 = 1./dy**2
        nperod = 1
        n = 10
        dz = 1./real(n)
        dzsq = 1./dz**2
        !
        !     GENERATE GRID POINTS FOR LATER USE.
        !
        do i = 1, l
            x(i) = (-pi) + real(i - 1)*dx
        end do
        do j = 1, m
            y(j) = (-pi) + real(j - 1)*dy
        end do
        !
        !     GENERATE COEFFICIENTS
        !
        a(1) = 0.
        b(1) = -2.*dzsq
        c(1) = -B(1)
        z(1) = 0.
        do k = 2, n
            z(k) = real(k - 1)*dz
            t = 1. + Z(k)
            a(k) = t**2*dzsq + t/dz
            b(k) = -2.*t**2*dzsq
            c(k) = t**2*dzsq - t/dz
        end do
        !
        !     GENERATE RIGHT SIDE OF EQUATION
        !
        do i = 1, l
            do j = 1, m
                do k = 2, n
                    f(i, j, k) = 2.*SIN(X(i))*SIN(Y(j))*(1. + Z(k))**4
                end do
            end do
        end do
        do i = 1, l
            do j = 1, l
                f(i, j, 1) = (10. + 8./dz)*SIN(X(i))*SIN(Y(j))
                f(i, j, n) = F(i, j, n) - C(n)*16.*SIN(X(i))*SIN(Y(j))
            end do
        end do
        c(n) = 0.
        !
        !     CALL POIS3D TO SOLVE EQUATIONS.
        !
        call POIS3D(lperod, l, c1, mperod, m, c2, nperod, n, a, b, c, &
            ldimf, mdimf, f, ierror)
        !
        !     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
        !
        !              U(X, Y, Z) = SIN(X)*SIN(Y)*(1+Z)**4
        !
        err = 0.
        do i = 1, l
            do j = 1, m
                do k = 1, n
                    t = abs(F(i, j, k)-SIN(X(i))*SIN(Y(j))*(1.+Z(k))**4)
                    err = max(t, err)
                end do
            end do
        end do
        !     Print earlier output from platforms with 32 and 64 bit floating point
        !     arithemtic followed by the output from this computer
        write( *, *) ''
        write( *, *) '    POIS3D TEST RUN *** '
        write( *, *) &
            '    Previous 64 bit floating point arithmetic result '
        write( *, *) '    IERROR = 0,  Discretization Error = 2.93277E-2'

        write( *, *) '    The output from your computer is: '
        write( *, *) '    IERROR =', ierror, ' Discretization Error = ', &
            err

    end subroutine pois3d_unit_test

    subroutine POIS3D( lperod, l, c1, mperod, m, c2, nperod, n, a, b, c, &
        ldimf, mdimf, f, ierror)
            !
        !     file pois3d.f
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
        !     SUBROUTINE POIS3D (LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, C, LDIMF,
        !    +                   MDIMF, F, IERROR)
        !
        !
        ! DIMENSION OF           A(N), B(N), C(N), F(LDIMF, MDIMF, N)
        ! ARGUMENTS
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
        !                        FOR UNKNOWN X VALUES, WHERE I=1, 2, ..., L,
        !                        J=1, 2, ..., M, AND K=1, 2, ..., N
        !
        !                        C1*(X(I-1, J, K) -2.*X(I, J, K) +X(I+1, J, K)) +
        !                        C2*(X(I, J-1, K) -2.*X(I, J, K) +X(I, J+1, K)) +
        !                        A(K)*X(I, J, K-1) +B(K)*X(I, J, K)+ C(K)*X(I, J, K+1)
        !                        = F(I, J, K)
        !
        !                        THE INDICES K-1 AND K+1 ARE EVALUATED MODULO N,
        !                        I.E. X(I, J, 0)=X(I, J, N) AND X(I, J, N+1)=X(I, J, 1).
        !                        THE UNKNOWNS
        !                        X(0, J, K), X(L+1, J, K), X(I, 0, K), AND X(I, M+1, K)
        !                        ARE ASSUMED TO TAKE ON CERTAIN PRESCRIBED
        !                        VALUES DESCRIBED BELOW.
        !
        ! USAGE                  CALL POIS3D (LPEROD, L, C1, MPEROD, M, C2, NPEROD,
        !                        N, A, B, C, LDIMF, MDIMF, F, IERROR)
        !
        ! ARGUMENTS
        !
        ! ON INPUT
        !                        LPEROD
        !                          INDICATES THE VALUES THAT X(0, J, K) AND
        !                          X(L+1, J, K) ARE ASSUMED TO HAVE.
        !                          = 0  X(0, J, K)=X(L, J, K), X(L+1, J, K)=X(1, J, K)
        !                          = 1  X(0, J, K) = 0,      X(L+1, J, K) = 0
        !                          = 2  X(0, J, K)=0,        X(L+1, J, K)=X(L-1, J, K)
        !                          = 3  X(0, J, K)=X(2, J, K), X(L+1, J, K)=X(L-1, J, K)
        !                          = 4  X(0, J, K)=X(2, J, K), X(L+1, J, K) = 0.
        !
        !                        L
        !                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
        !                          L MUST BE AT LEAST 3.
        !
        !                        C1
        !                          REAL CONSTANT IN THE ABOVE LINEAR SYSTEM
        !                          OF EQUATIONS TO BE SOLVED.
        !
        !                        MPEROD
        !                          INDICATES THE VALUES THAT X(I, 0, K) AND
        !                          X(I, M+1, K) ARE ASSUMED TO HAVE.
        !                          = 0  X(I, 0, K)=X(I, M, K), X(I, M+1, K)=X(I, 1, K)
        !                          = 1  X(I, 0, K)=0,        X(I, M+1, K)=0
        !                          = 2  X(I, 0, K)=0,        X(I, M+1, K)=X(I, M-1, K)
        !                          = 3  X(I, 0, K)=X(I, 2, K)  X(I, M+1, K)=X(I, M-1, K)
        !                          = 4  X(I, 0, K)=X(I, 2, K)  X(I, M+1, K)=0
        !
        !                        M
        !                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
        !                          M MUST BE AT LEAST 3.
        !
        !                        C2
        !                          REAL CONSTANT IN THE ABOVE LINEAR SYSTEM
        !                          OF EQUATIONS TO BE SOLVED.
        !
        !                        NPEROD
        !                          = 0  IF A(1) AND C(N) ARE NOT ZERO.
        !                          = 1  IF A(1) = C(N) = 0.
        !
        !                        N
        !                          THE NUMBER OF UNKNOWNS IN THE K-DIRECTION.
        !                          N MUST BE AT LEAST 3.
        !
        !                        A, B, C
        !                          ONE-DIMENSIONAL ARRAYS OF LENGTH N THAT
        !                          SPECIFY THE COEFFICIENTS IN THE LINEAR
        !                          EQUATIONS GIVEN ABOVE.
        !
        !                          IF NPEROD = 0 THE ARRAY ELEMENTS MUST NOT
        !                          DEPEND UPON INDEX K, BUT MUST BE CONSTANT.
        !                          SPECIFICALLY, THE SUBROUTINE CHECKS THE
        !                          FOLLOWING CONDITION
        !                            A(K) = C(1)
        !                            C(K) = C(1)
        !                            B(K) = B(1)
        !                          FOR K=1, 2, ..., N.
        !
        !                        LDIMF
        !                          THE ROW (OR FIRST) DIMENSION OF THE THREE-
        !                          DIMENSIONAL ARRAY F AS IT APPEARS IN THE
        !                          PROGRAM CALLING POIS3D.  THIS PARAMETER IS
        !                          USED TO SPECIFY THE VARIABLE DIMENSION
        !                          OF F.  LDIMF MUST BE AT LEAST L.
        !
        !                        MDIMF
        !                          THE COLUMN (OR SECOND) DIMENSION OF THE THREE
        !                          DIMENSIONAL ARRAY F AS IT APPEARS IN THE
        !                          PROGRAM CALLING POIS3D.  THIS PARAMETER IS
        !                          USED TO SPECIFY THE VARIABLE DIMENSION
        !                          OF F.  MDIMF MUST BE AT LEAST M.
        !
        !                        F
        !                          A THREE-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
        !                          OF EQUATIONS GIVEN ABOVE.  F MUST BE
        !                          DIMENSIONED AT LEAST L X M X N.
        !
        ! ON OUTPUT
        !
        !                        F
        !                          CONTAINS THE SOLUTION X.
        !
        !                        IERROR
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
        !                          SOLUTION IS NOT ATTEMPTED.
        !                          = 0  NO ERROR
        !                          = 1  IF LPEROD .LT. 0 OR .GT. 4
        !                          = 2  IF L .LT. 3
        !                          = 3  IF MPEROD .LT. 0 OR .GT. 4
        !                          = 4  IF M .LT. 3
        !                          = 5  IF NPEROD .LT. 0 OR .GT. 1
        !                          = 6  IF N .LT. 3
        !                          = 7  IF LDIMF .LT. L
        !                          = 8  IF MDIMF .LT. M
        !                          = 9  IF A(K) .NE. C(1) OR C(K) .NE. C(1)
        !                               OR B(I) .NE.B(1) FOR SOME K=1, 2, ..., N.
        !                          = 10 IF NPEROD = 1 AND A(1) .NE. 0
        !                               OR C(N) .NE. 0
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
        !
        !                          SINCE THIS IS THE ONLY MEANS OF INDICATING A
        !                          POSSIBLY INCORRECT CALL TO POIS3D, THE USER
        !                          SHOULD TEST IERROR AFTER THE CALL.
        !
        ! SPECIAL CONDITIONS     NONE
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED files         fish.f, comf.f, fftpack.f
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
        ! ALGORITHM              THIS SUBROUTINE SOLVES THREE-DIMENSIONAL BLOCK
        !                        TRIDIAGONAL LINEAR SYSTEMS ARISING FROM FINITE
        !                        DIFFERENCE APPROXIMATIONS TO THREE-DIMENSIONAL
        !                        POISSON EQUATIONS USING THE FFT PACKAGE
        !                        FFTPACK WRITTEN BY PAUL SWARZTRAUBER.
        !
        ! TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
        !                        IS ROUGHLY PROPORTIONAL TO
        !                          L*M*N*(LOG2(L)+LOG2(M)+5)
        !                        BUT ALSO DEPENDS ON INPUT PARAMETERS LPEROD
        !                        AND MPEROD.
        !
        ! ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
        !                        UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
        !                        CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
        !                        IN THE 'PURPOSE' SECTION WITH
        !                          A(K) = C(K) = -0.5*B(K) = 1,  K=1, 2, ..., N
        !                        AND, WHEN NPEROD = 1
        !                          A(1) = C(N) = 0
        !                          A(N) = C(1) = 2.
        !
        !                        THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
        !                        SYSTEM AND, USING DOUBLE PRECISION, A RIGHT
        !                        SIDE Y WAS COMPUTED.  USING THIS ARRAY Y
        !                        SUBROUTINE POIS3D WAS CALLED TO PRODUCE AN
        !                        APPROXIMATE SOLUTION Z.  RELATIVE ERROR
        !
        !                        E = MAX(abs(Z(I, J, K)-X(I, J, K)))/MAX(abs(X(I, J, K
        !
        !                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
        !                        OVER I=1, 2, ..., L, J=1, 2, ..., M AND K=1, 2, ..., N.
        !                        VALUES OF E ARE GIVEN IN THE TABLE BELOW FOR
        !                        SOME TYPICAL VALUES OF L, M AND N.
        !
        !                        L(=M=N)   LPEROD    MPEROD       E
        !                        ------    ------    ------     ------
        !
        !                          16        0         0        1.E-13
        !                          15        1         1        4.E-13
        !                          17        3         3        2.E-13
        !                          32        0         0        2.E-13
        !                          31        1         1        2.E-12
        !                          33        3         3        7.E-13
        !
        ! REFERENCES              NONE
        ! ********************************************************************

        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer  :: lperod
        integer  :: l
        integer  :: mperod
        integer  :: m
        integer  :: nperod
        integer  :: n
        integer  :: ldimf
        integer  :: mdimf
        integer  :: ierror
        real  :: c1
        real  :: c2
        real  :: a(*)
        real  :: b(*)
        real  :: c(*)
        real  :: f(ldimf, mdimf, *)
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: lp, mp, np, k, irwk, icwk
        
        !-----------------------------------------------
        lp = lperod + 1
        mp = mperod + 1
        np = nperod + 1
        !
        !     CHECK FOR INVALID INPUT.
        !
        ierror = 0
        if (lp<1 .or. lp>5) ierror = 1
        if (l < 3) ierror = 2
        if (mp<1 .or. mp>5) ierror = 3
        if (m < 3) ierror = 4
        if (np<1 .or. np>2) ierror = 5
        if (n < 3) ierror = 6
        if (ldimf < l) ierror = 7
        if (mdimf < m) ierror = 8
        if (np == 1) then
            do k = 1, n
                if (A(k) /= C(1)) go to 102
                if (C(k) /= C(1)) go to 102
                if (B(k) /= B(1)) go to 102
            end do
            go to 104
102     continue
        ierror = 9
    end if
    if (nperod==1 .and. (A(1)/=0. .or. C(n)/=0.)) ierror = 10
! 104 IF (IERROR .NE. 0) GO TO 122
104 continue
    if (ierror /= 0) return
    !     allocate required work space length (generous estimate)
    irwk=30+l+m+2*n+max(l, m, n)+7*(INT((l+1)/2)+INT((m+1)/2))
    icwk = 0
    call workspace%create( irwk, icwk, ierror )
    !     check that allocation was successful
    if (ierror == 20) return
    call pois3dd(lperod, l, c1, mperod, m, c2, nperod, n, a, b, c, ldimf, &
        mdimf, f, ierror, workspace%rew)
    !     release work space
    call workspace%destroy()

end subroutine POIS3D

subroutine POIS3DD(lperod, l, c1, mperod, m, c2, nperod, n, a, b, &
    c, ldimf, mdimf, f, ierror, w)

    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: lperod
    integer  :: l
    integer , intent (in) :: mperod
    integer  :: m
    integer , intent (in) :: nperod
    integer  :: n
    integer  :: ldimf
    integer  :: mdimf
    integer  :: ierror
    real  :: c1
    real  :: c2
    real  :: a(*)
    real  :: b(*)
    real  :: c(*)
    real  :: f(ldimf, mdimf, *)
    real  :: w(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: lp, mp, np, iwyrt, iwt, iwd, iwbb, iwx, iwy, nh, nhm1, &
        nodd, i, j, k
    real, dimension(6) :: save
    !-----------------------------------------------
    lp = lperod + 1
    mp = mperod + 1
    np = nperod + 1
    iwyrt = l + 1
    iwt = iwyrt + m
    iwd = iwt + max(l, m, n) + 1
    iwbb = iwd + n
    iwx = iwbb + n
    iwy = iwx + 7*((l + 1)/2) + 15
    go to (105, 114) np
!
!     REORDER UNKNOWNS WHEN NPEROD = 0.
!
105 continue
    nh = (n + 1)/2
    nhm1 = nh - 1
    nodd = 1
    if (2*nh == n) nodd = 2
    do i = 1, l
        do j = 1, m
            do k = 1, nhm1
                w(k) = F(i, j, nh-k) - F(i, j, k+nh)
                w(k+nh) = F(i, j, nh-k) + F(i, j, k+nh)
            end do
            w(nh) = 2.*F(i, j, nh)
            go to (108, 107) nodd
107     continue
        w(n) = 2.*F(i, j, n)
108 continue
    f(i, j, :n) = W(:n)
end do
end do
save(1) = C(nhm1)
save(2) = A(nh)
save(3) = C(nh)
save(4) = B(nhm1)
save(5) = B(n)
save(6) = A(n)
c(nhm1) = 0.
a(nh) = 0.
c(nh) = 2.*C(nh)
select case (nodd)
    case default
        b(nhm1) = B(nhm1) - A(nh-1)
        b(n) = B(n) + A(n)
    case (2)
        a(n) = C(nh)
end select
114 continue
    call POS3D1 (lp, l, mp, m, n, a, b, c, ldimf, mdimf, f, w, W(iwyrt &
        ), W(iwt), W(iwd), W(iwx), W(iwy), c1, c2, W(iwbb))
    go to (115, 122) np
115 continue
    do i = 1, l
        do j = 1, m
            w(nh-1:nh-nhm1:(-1))=0.5*(F(i, j, nh+1:nhm1+nh)+F(i, j, :nhm1))
            w(nh+1:nhm1+nh) = 0.5*(F(i, j, nh+1:nhm1+nh)-F(i, j, :nhm1))
            w(nh) = 0.5*F(i, j, nh)
            go to (118, 117) nodd
117     continue
        w(n) = 0.5*F(i, j, n)
118 continue
    f(i, j, :n) = W(:n)
end do
end do
c(nhm1) = SAVE(1)
a(nh) = SAVE(2)
c(nh) = SAVE(3)
b(nhm1) = SAVE(4)
b(n) = SAVE(5)
a(n) = SAVE(6)
122 continue

end subroutine POIS3DD

subroutine POS3D1(lp, l, mp, m, n, a, b, c, ldimf, mdimf, f, xrt, &
    yrt, t, d, wx, wy, c1, c2, bb)

    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: lp
    integer , intent (in) :: l
    integer , intent (in) :: mp
    integer , intent (in) :: m
    integer , intent (in) :: n
    integer , intent (in) :: ldimf
    integer , intent (in) :: mdimf
    real , intent (in) :: c1
    real , intent (in) :: c2
    real  :: a(*)
    real , intent (in) :: b(*)
    real  :: c(*)
    real , intent (in out) :: f(ldimf, mdimf, 1)
    real , intent (in out) :: xrt(*)
    real , intent (in out) :: yrt(*)
    real  :: t(*)
    real  :: d(*)
    real  :: wx(*)
    real  :: wy(*)
    real  :: bb(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: lr, mr, nr, lrdel, i, mrdel, j, ifwrd, is, k
    real :: pi, scalx, dx, di, scaly, dy, dj
    !-----------------------------------------------
    pi = acos( -1.0 )
    lr = l
    mr = m
    nr = n
    !
    !     GENERATE TRANSFORM ROOTS
    !
    lrdel = ((lp - 1)*(lp - 3)*(lp - 5))/3
    scalx = lr + lrdel
    dx = pi/(2.*scalx)
    go to (108, 103, 101, 102, 101) lp
101 continue
    di = 0.5
    scalx = 2.*scalx
    go to 104
102 continue
    di = 1.0
    go to 104
103 continue
    di = 0.0
104 continue
    do i = 1, lr
        xrt(i) = -4.*c1*SIN((real(i) - di)*dx)**2
    end do
    scalx = 2.*scalx
    go to (112, 106, 110, 107, 111) lp
106 continue
    call SINTI (lr, wx)
    go to 112
107 continue
    call COSTI (lr, wx)
    go to 112
108 continue
    xrt(1) = 0.
    xrt(lr) = -4.*c1
    do i = 3, lr, 2
        xrt(i-1) = -4.*c1*SIN(real(i - 1)*dx)**2
        xrt(i) = XRT(i-1)
    end do
    call RFFTI (lr, wx)
    go to 112
110 continue
    call SINQI (lr, wx)
    go to 112
111 continue
    call COSQI (lr, wx)
112 continue
    mrdel = ((mp - 1)*(mp - 3)*(mp - 5))/3
    scaly = mr + mrdel
    dy = pi/(2.*scaly)
    go to (120, 115, 113, 114, 113) mp
113 continue
    dj = 0.5
    scaly = 2.*scaly
    go to 116
114 continue
    dj = 1.0
    go to 116
115 continue
    dj = 0.0
116 continue
    do j = 1, mr
        yrt(j) = -4.*c2*SIN((real(j) - dj)*dy)**2
    end do
    scaly = 2.*scaly
    go to (124, 118, 122, 119, 123) mp
118 continue
    call SINTI (mr, wy)
    go to 124
119 continue
    call COSTI (mr, wy)
    go to 124
120 continue
    yrt(1) = 0.
    yrt(mr) = -4.*c2
    do j = 3, mr, 2
        yrt(j-1) = -4.*c2*SIN(real(j - 1)*dy)**2
        yrt(j) = YRT(j-1)
    end do
    call RFFTI (mr, wy)
    go to 124
122 continue
    call SINQI (mr, wy)
    go to 124
123 continue
    call COSQI (mr, wy)
124 continue
    ifwrd = 1
    is = 1
125 continue
    !
    !     TRANSFORM X
    !
    do j=1, mr
        do k=1, nr
            do i=1, lr
                t(i) = F(i, j, k)
            end do
            go to (127, 130, 131, 134, 135), lp
127         go to (128, 129), ifwrd
128         call RFFTF (lr, t, wx)
            go to 138
129         call RFFTB (lr, t, wx)
            go to 138
130         call SINT (lr, t, wx)
            go to 138
131         go to (132, 133), ifwrd
132         call SINQF (lr, t, wx)
            go to 138
133         call SINQB (lr, t, wx)
            go to 138
134         call COST (lr, t, wx)
            go to 138
135         go to (136, 137), ifwrd
136         call COSQF (lr, t, wx)
            go to 138
137         call COSQB (lr, t, wx)
138     continue
        do i=1, lr
            f(i, j, k) = T(i)
        end do
    end do
end do
go to (142, 164) ifwrd
!
!     TRANSFORM Y
!
142 continue
    do i=1, lr
        do k=1, nr
            do j=1, mr
                t(j) = F(i, j, k)
            end do
            go to (144, 147, 148, 151, 152), mp
144         go to (145, 146), ifwrd
145         call RFFTF (mr, t, wy)
            go to 155
146         call RFFTB (mr, t, wy)
            go to 155
147         call SINT (mr, t, wy)
            go to 155
148         go to (149, 150), ifwrd
149         call SINQF (mr, t, wy)
            go to 155
150         call SINQB (mr, t, wy)
            go to 155
151         call COST (mr, t, wy)
            go to 155
152         go to (153, 154), ifwrd
153         call COSQF (mr, t, wy)
            go to 155
154         call COSQB (mr, t, wy)
155     continue
        do j=1, mr
            f(i, j, k) = T(j)
        end do
    end do
end do
go to (159, 125) ifwrd
159 continue
    do i = 1, lr
        do j = 1, mr
            bb(:nr) = B(:nr) + XRT(i) + YRT(j)
            t(:nr) = F(i, j, :nr)
            call TRID (nr, a, bb, c, t, d)
            f(i, j, :nr) = T(:nr)
        end do
    end do
    ifwrd = 2
    is = -1
    go to 142
164 continue
    f(:lr, :mr, :nr) = F(:lr, :mr, :nr)/(scalx*scaly)

end subroutine POS3D1

subroutine TRID(mr, a, b, c, y, d)

    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer , intent (in) :: mr
    real , intent (in) :: a(*)
    real , intent (in) :: b(*)
    real , intent (in) :: c(*)
    real , intent (in out) :: y(*)
    real , intent (in out) :: d(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: m, mm1, i, ip
    real :: z
    !-----------------------------------------------
    m = mr
    mm1 = m - 1
    z = 1./B(1)
    d(1) = C(1)*z
    y(1) = Y(1)*z
    do i = 2, mm1
        z = 1./(B(i)-A(i)*D(i-1))
        d(i) = C(i)*z
        y(i) = (Y(i)-A(i)*Y(i-1))*z
    end do
    z = B(m) - A(m)*D(mm1)
    if (z == 0.) then
        y(m) = 0.
    else
        y(m) = (Y(m)-A(m)*Y(mm1))/z
    end if
    do ip = 1, mm1
        i = m - ip
        y(i) = Y(i) - D(i)*Y(i+1)
    end do

end subroutine TRID

end module module_pois3d
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
