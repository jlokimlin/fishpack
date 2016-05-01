module module_hwscyl

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_genbun, only: &
        genbunn

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hwscyl


contains


    subroutine hwscyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror)
        !
        !     file hwscyl.f
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
        !     SUBROUTINE hwscyl (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD,
        !    +                   ELMBDA, F, IDIMF, PERTRB, ierror)
        !
        !
        ! DIMENSION OF           BDA(N), BDB(N), BDC(M), BDD(M), F(IDIMF, N+1)
        ! ARGUMENTS
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION
        !                        TO THE HELMHOLTZ EQUATION IN CYLINDRICAL
        !                        COORDINATES.  THIS MODIFIED HELMHOLTZ EQUATION
        !
        !                          (1/R)(D/DR)(R(DU/DR)) + (D/DZ)(DU/DZ)
        !
        !                          + (LAMBDA/R**2)U = F(R, Z)
        !
        !                        RESULTS FROM THE FOURIER TRANSFORM OF THE
        !                        THREE-DIMENSIONAL POISSON EQUATION.
        !
        ! USAGE                  CALL hwscyl (A, B, M, MBDCND, BDA, BDB, C, D, N,
        !                                     NBDCND, BDC, BDD, ELMBDA, F, IDIMF,
        !                                     PERTRB, ierror, W)
        !
        ! ARGUMENTS
        ! ON INPUT               A, B
        !                          THE RANGE OF R, I.E., A .LE. R .LE. B.
        !                          A MUST BE LESS THAN B AND A MUST BE
        !                          NON-NEGATIVE.
        !
        !                        M
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL (A, B) IS SUBDIVIDED.  HENCE,
        !                          THERE WILL BE M+1 GRID POINTS IN THE
        !                          R-DIRECTION GIVEN BY R(I) = A+(I-1)DR,
        !                          FOR I = 1, 2, ..., M+1, WHERE DR = (B-A)/M
        !                          IS THE PANEL WIDTH. M MUST BE GREATER
        !                          THAN 3.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT R = A AND R = B.
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               R = A AND R = B.
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               R = A AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO R IS
        !                               SPECIFIED AT R = B.
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO R IS SPECIFIED AT
        !                               R = A (SEE NOTE BELOW) AND R = B.
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO R IS SPECIFIED AT
        !                               R = A (SEE NOTE BELOW) AND THE
        !                               SOLUTION IS SPECIFIED AT R = B.
        !                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
        !                               R = A = 0 AND THE SOLUTION IS
        !                               SPECIFIED AT R = B.
        !                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
        !                               R = A = 0 AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO R IS SPECIFIED
        !                               AT R = B.
        !
        !                          IF A = 0, DO NOT USE MBDCND = 3 OR 4,
        !                          BUT INSTEAD USE MBDCND = 1, 2, 5, OR 6  .
        !
        !                        BDA
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
        !                          THE SOLUTION WITH RESPECT TO R AT R = A.
        !
        !                          WHEN MBDCND = 3 OR 4,
        !                            BDA(J) = (D/DR)U(A, Z(J)), J = 1, 2, ..., N+1.
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
        !                          A DUMMY VARIABLE.
        !
        !                        BDB
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE
        !                          OF THE SOLUTION WITH RESPECT TO R AT R = B.
        !
        !                          WHEN MBDCND = 2, 3, OR 6,
        !                            BDB(J) = (D/DR)U(B, Z(J)), J = 1, 2, ..., N+1.
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
        !                          A DUMMY VARIABLE.
        !
        !                        C, D
        !                          THE RANGE OF Z, I.E., C .LE. Z .LE. D.
        !                          C MUST BE LESS THAN D.
        !
        !                        N
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL (C, D) IS SUBDIVIDED.  HENCE,
        !                          THERE WILL BE N+1 GRID POINTS IN THE
        !                          Z-DIRECTION GIVEN BY Z(J) = C+(J-1)DZ,
        !                          FOR J = 1, 2, ..., N+1,
        !                          WHERE DZ = (D-C)/N IS THE PANEL WIDTH.
        !                          N MUST BE GREATER THAN 3.
        !
        !                        NBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT Z = C AND Z = D.
        !
        !                          = 0  IF THE SOLUTION IS PERIODIC IN Z,
        !                               I.E., U(I, 1) = U(I, N+1).
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               Z = C AND Z = D.
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               Z = C AND THE DERIVATIVE OF
        !                               THE SOLUTION WITH RESPECT TO Z IS
        !                               SPECIFIED AT Z = D.
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO Z IS
        !                               SPECIFIED AT Z = C AND Z = D.
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO Z IS SPECIFIED AT
        !                               Z = C AND THE SOLUTION IS SPECIFIED
        !                               AT Z = D.
        !
        !                        BDC
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE
        !                          OF THE SOLUTION WITH RESPECT TO Z AT Z = C.
        !
        !                          WHEN NBDCND = 3 OR 4,
        !                            BDC(I) = (D/DZ)U(R(I), C), I = 1, 2, ..., M+1.
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
        !                          A DUMMY VARIABLE.
        !
        !                        BDD
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
        !                          THE SOLUTION WITH RESPECT TO Z AT Z = D.
        !
        !                          WHEN NBDCND = 2 OR 3,
        !                            BDD(I) = (D/DZ)U(R(I), D), I = 1, 2, ..., M+1
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
        !                          A DUMMY VARIABLE.
        !
        !                        ELMBDA
        !                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
        !                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
        !                          MAY NOT EXIST.  HOWEVER, hwscyl WILL
        !                          ATTEMPT TO FIND A SOLUTION.  LAMBDA MUST
        !                          BE ZERO WHEN MBDCND = 5 OR 6  .
        !
        !                        F
        !                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
        !                          LEAST (M+1)*(N+1), SPECIFYING VALUES
        !                          OF THE RIGHT SIDE OF THE HELMHOLTZ
        !                          EQUATION AND BOUNDARY DATA (IF ANY).
        !
        !                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
        !                          FOR I = 2, 3, ..., M AND J = 2, 3, ..., N
        !                          F(I, J) = F(R(I), Z(J)).
        !
        !                          ON THE BOUNDARIES F IS DEFINED AS FOLLOWS:
        !                          FOR J = 1, 2, ..., N+1 AND I = 1, 2, ..., M+1
        !
        !                          MBDCND   F(1, J)            F(M+1, J)
        !                          ------   ---------         ---------
        !
        !                            1      U(A, Z(J))         U(B, Z(J))
        !                            2      U(A, Z(J))         F(B, Z(J))
        !                            3      F(A, Z(J))         F(B, Z(J))
        !                            4      F(A, Z(J))         U(B, Z(J))
        !                            5      F(0, Z(J))         U(B, Z(J))
        !                            6      F(0, Z(J))         F(B, Z(J))
        !
        !                          NBDCND   F(I, 1)            F(I, N+1)
        !                          ------   ---------         ---------
        !
        !                            0      F(R(I), C)         F(R(I), C)
        !                            1      U(R(I), C)         U(R(I), D)
        !                            2      U(R(I), C)         F(R(I), D)
        !                            3      F(R(I), C)         F(R(I), D)
        !                            4      F(R(I), C)         U(R(I), D)
        !
        !                          NOTE:
        !                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
        !                          U AND THE RIGHT SIDE F AT A CORNER THEN
        !                          THE SOLUTION MUST BE SPECIFIED.
        !
        !                        IDIMF
        !                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
        !                          F AS IT APPEARS IN THE PROGRAM CALLING
        !                          hwscyl.  THIS PARAMETER IS USED TO SPECIFY
        !                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
        !                          BE AT LEAST M+1  .
        !
        !
        ! ON OUTPUT              F
        !                          CONTAINS THE SOLUTION U(I, J) OF THE FINITE
        !                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
        !                          (R(I), Z(J)), I =1, 2, ..., M+1, J =1, 2, ..., N+1.
        !
        !                        PERTRB
        !                          IF ONE SPECIFIES A COMBINATION OF PERIODIC,
        !                          DERIVATIVE, AND UNSPECIFIED BOUNDARY
        !                          CONDITIONS FOR A POISSON EQUATION
        !                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
        !                          PERTRB IS A CONSTANT, CALCULATED AND
        !                          SUBTRACTED FROM F, WHICH ENSURES THAT A
        !                          SOLUTION EXISTS.  hwscyl THEN COMPUTES
        !                          THIS SOLUTION, WHICH IS A LEAST SQUARES
        !                          SOLUTION TO THE ORIGINAL APPROXIMATION.
        !                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
        !                          A SOLUTION.  HENCE, THE SOLUTION IS NOT
        !                          UNIQUE. THE VALUE OF PERTRB SHOULD BE
        !                          SMALL COMPARED TO THE RIGHT SIDE F.
        !                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
        !                          ESSENTIALLY DIFFERENT PROBLEM.  THIS
        !                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
        !                          THAT A MEANINGFUL SOLUTION HAS BEEN OBTAINED.
        !
        !                        ierror
        !                          AN ERROR FLAG WHICH INDICATES INVALID INPUT
        !                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 11,
        !                          A SOLUTION IS NOT ATTEMPTED.
        !
        !                          =  0  NO ERROR.
        !                          =  1  A .LT. 0  .
        !                          =  2  A .GE. B.
        !                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6  .
        !                          =  4  C .GE. D.
        !                          =  5  N .LE. 3
        !                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4  .
        !                          =  7  A = 0, MBDCND = 3 OR 4  .
        !                          =  8  A .GT. 0, MBDCND .GE. 5  .
        !                          =  9  A = 0, LAMBDA .NE. 0, MBDCND .GE. 5  .
        !                          = 10  IDIMF .LT. M+1  .
        !                          = 11  LAMBDA .GT. 0  .
        !                          = 12  M .LE. 3
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
        !
        !                          SINCE THIS IS THE ONLY MEANS OF INDICATING
        !                          A POSSIBLY INCORRECT CALL TO hwscyl, THE
        !                          USER SHOULD TEST ierror AFTER THE CALL.
        !
        !
        ! SPECIAL CONDITIONS     NONE
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED files         fish.f, genbun.f, gnbnaux.f, comf.f
        !
        ! LANGUAGE               FORTRAN 90
        !
        ! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
        !                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
        !                        LIBRARIES IN January 1980.
        !                        Revised in June 2004 by John Adams using
        !                        Fortran 90 dynamically allocated work space.
        !
        !
        ! PORTABILITY            FORTRAN 90
        !
        ! ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
        !                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
        !                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
        !                        AND THEN CALLS genbun TO SOLVE THE SYSTEM.
        !
        ! TIMING                 FOR LARGE  M AND N, THE OPERATION COUNT
        !                        IS ROUGHLY PROPORTIONAL TO
        !                          M*N*(LOG2(N)
        !                        BUT ALSO DEPENDS ON INPUT PARAMETERS NBDCND
        !                        AND MBDCND.
        !
        ! ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A LOSS
        !                        OF NO MORE THAN THREE SIGNIFICANT DIGITS FOR N
        !                        AND M AS LARGE AS 64.  MORE DETAILS ABOUT
        !                        ACCURACY CAN BE FOUND IN THE DOCUMENTATION FOR
        !                        SUBROUTINE genbun WHICH IS THE ROUTINE THAT
        !                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
        !
        ! REFERENCES             SWARZTRAUBER, P. AND R. SWEET, "EFFICIENT
        !                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
        !                        ELLIPTIC EQUATIONS"
        !                          NCAR TN/IA-109, JULY, 1975, 138 PP.
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip),          intent (in)     :: m
        integer (ip),          intent (in)     :: mbdcnd
        integer (ip),          intent (in)     :: n
        integer (ip),          intent (in)     :: nbdcnd
        integer (ip),          intent (in)     :: idimf
        integer (ip),          intent (out)    :: ierror
        real (wp),             intent (in)     :: a
        real (wp),             intent (in)     :: b
        real (wp),             intent (in)     :: c
        real (wp),             intent (in)     :: d
        real (wp),             intent (in)     :: elmbda
        real (wp),             intent (out)    :: pertrb
        real (wp), contiguous, intent (in)     :: bda(:)
        real (wp), contiguous, intent (in)     :: bdb(:)
        real (wp), contiguous, intent (in)     :: bdc(:)
        real (wp), contiguous, intent (in)     :: bdd(:)
        real (wp), contiguous, intent (in out) :: f(:,:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        type (FishpackWorkspace) :: workspace
        integer (ip)             :: irwk
        !-----------------------------------------------

        ! initialize error flag
        ierror = 0

        ! check if input values are valid
        if (a < 0.) ierror = 1
        if (a >= b) ierror = 2
        if (mbdcnd<=0 .or. mbdcnd>=7) ierror = 3
        if (c >= d) ierror = 4
        if (n <= 3) ierror = 5
        if (nbdcnd<=(-1) .or. nbdcnd>=5) ierror = 6
        if (a==0. .and. (mbdcnd==3 .or. mbdcnd==4)) ierror = 7
        if (a>0. .and. mbdcnd>=5) ierror = 8
        if (a==0. .and. elmbda/=0. .and. mbdcnd>=5) ierror = 9
        if (idimf < m + 1) ierror = 10
        if (m <= 3) ierror = 12
        if (ierror /= 0) return

        ! compute real workspace size
        call workspace%get_genbun_workspace_dimensions(n, m, irwk)

        ! Allocate workspace
        irwk = irwk + 3*m
        associate( icwk => 0 )
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! Check if allocation was succcessful
        if (ierror == 20) return

        ! solve system
        associate( rew => workspace%real_workspace )
            call hwscyll(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hwscyl

    subroutine hwscylL(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, w)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: mbdcnd
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: nbdcnd
        integer (ip), intent (in)     :: idimf
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in)     :: a
        real (wp),    intent (in)     :: b
        real (wp),    intent (in)     :: c
        real (wp),    intent (in)     :: d
        real (wp),    intent (in)     :: elmbda
        real (wp),    intent (out)    :: pertrb
        real (wp),    intent (in)     :: bda(*)
        real (wp),    intent (in)     :: bdb(*)
        real (wp),    intent (in)     :: bdc(*)
        real (wp),    intent (in)     :: bdd(*)
        real (wp),    intent (in out) :: f(idimf,*)
        real (wp),    intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: mp1, np1, np, mstart, mstop, munk, nstart, nstop, nunk &
            , id2, id3, id4, id5, id6, istart, ij, i, j, k, l, nsp1, nstm1 &
            , ierr1, i1
        real::deltar, dlrby2, dlrsq, deltht, dlthsq, a1, r, a2, s, s1, s2
        !-----------------------------------------------
        mp1 = m + 1
        deltar = (b - a)/real(m)
        dlrby2 = deltar/2.
        dlrsq = deltar**2
        np1 = n + 1
        deltht = (d - c)/real(n)
        dlthsq = deltht**2
        np = nbdcnd + 1
        !
        !     DEFINE RANGE OF INDICES I AND J FOR UNKNOWNS U(I, J).
        !
        mstart = 2
        mstop = m
        go to (104, 103, 102, 101, 101, 102) mbdcnd
101 continue
    mstart = 1
    go to 104
102 continue
    mstart = 1
103 continue
    mstop = mp1
104 continue
    munk = mstop - mstart + 1
    nstart = 1
    nstop = n
    go to (108, 105, 106, 107, 108) np
105 continue
    nstart = 2
    go to 108
106 continue
    nstart = 2
107 continue
    nstop = np1
108 continue
    nunk = nstop - nstart + 1
    !
    !     DEFINE A, B, C COEFFICIENTS IN W-ARRAY.
    !
    id2 = munk
    id3 = id2 + munk
    id4 = id3 + munk
    id5 = id4 + munk
    id6 = id5 + munk
    istart = 1
    a1 = 2./dlrsq
    ij = 0
    if (mbdcnd==3 .or. mbdcnd==4) ij = 1
    if (mbdcnd > 4) then
        w(1) = 0.
        w(id2+1) = -2.*a1
        w(id3+1) = 2.*a1
        istart = 2
        ij = 1
    end if
    do i = istart, munk
        r = a + real(i - ij)*deltar
        j = id5 + i
        w(j) = r
        j = id6 + i
        w(j) = 1./r**2
        w(i) = (r - dlrby2)/(r*dlrsq)
        j = id3 + i
        w(j) = (r + dlrby2)/(r*dlrsq)
        k = id6 + i
        j = id2 + i
        w(j) = (-a1) + elmbda*W(k)
    end do
    go to (114, 111, 112, 113, 114, 112) mbdcnd
111 continue
    w(id2) = a1
    go to 114
112 continue
    w(id2) = a1
113 continue
    w(id3+1) = a1*real(istart)
114 continue
    go to (115, 115, 117, 117, 119, 119) mbdcnd
115 continue
    a1 = W(1)
    f(2, nstart:nstop) = F(2, nstart:nstop) - a1*F(1, nstart:nstop)
    go to 119
117 continue
    a1 = 2.*deltar*W(1)
    f(1, nstart:nstop) = F(1, nstart:nstop) + a1*BDA(nstart:nstop)
119 continue
    go to (120, 122, 122, 120, 120, 122) mbdcnd
120 continue
    a1 = W(id4)
    f(m, nstart:nstop) = F(m, nstart:nstop) - a1*F(mp1, nstart:nstop)
    go to 124
122 continue
    a1 = 2.*deltar*W(id4)
    f(mp1, nstart:nstop) = F(mp1, nstart:nstop) - a1*BDB(nstart:nstop)
!
!     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
!
124 continue
    a1 = 1./dlthsq
    l = id5 - mstart + 1
    go to (134, 125, 125, 127, 127) np
125 continue
    f(mstart:mstop, 2) = F(mstart:mstop, 2) - a1*F(mstart:mstop, 1)
    go to 129
127 continue
    a1 = 2./deltht
    f(mstart:mstop, 1) = F(mstart:mstop, 1) + a1*BDC(mstart:mstop)
129 continue
    a1 = 1./dlthsq
    go to (134, 130, 132, 132, 130) np
130 continue
    f(mstart:mstop, n) = F(mstart:mstop, n) - a1*F(mstart:mstop, np1)
    go to 134
132 continue
    a1 = 2./deltht
    f(mstart:mstop, np1) = F(mstart:mstop, np1) - a1*BDD(mstart:mstop)
134 continue
    pertrb = 0.
    if (elmbda >= 0.) then
        if (elmbda /= 0.) then
            ierror = 11
        else
            w(id5+1) = 0.5*(W(id5+2)-dlrby2)
            go to (146, 146, 138, 146, 146, 137) mbdcnd
137     continue
        w(id5+1) = 0.5*W(id5+1)
138 continue
    go to (140, 146, 146, 139, 146) np
139 continue
    a2 = 2.
    go to 141
140 continue
    a2 = 1.
141 continue
    k = id5 + munk
    w(k) = 0.5*(W(k-1)+dlrby2)
    s = 0.
    do i = mstart, mstop
        s1 = 0.
        nsp1 = nstart + 1
        nstm1 = nstop - 1
        s1 = SUM(F(i, nsp1:nstm1))
        k = i + l
        s = s + (a2*s1 + F(i, nstart)+F(i, nstop))*W(k)
    end do
    s2 = real(m)*a + (0.75 + real((m - 1)*(m + 1)))*dlrby2
    if (mbdcnd == 3) s2 = s2 + 0.25*dlrby2
    s1 = (2. + a2*real(nunk - 2))*s2
    pertrb = s/s1
    f(mstart:mstop, nstart:nstop) = F(mstart:mstop, nstart:nstop) &
        - pertrb
end if
end if
146 continue
    w(:mstop-mstart+1) = W(:mstop-mstart+1)*dlthsq
    w(id2+1:mstop-mstart+1+id2) = W(id2+1:mstop-mstart+1+id2)*dlthsq
    w(id3+1:mstop-mstart+1+id3) = W(id3+1:mstop-mstart+1+id3)*dlthsq
    f(mstart:mstop, nstart:nstop) = F(mstart:mstop, nstart:nstop)*dlthsq
    w(1) = 0.
    w(id4) = 0.
    !
    !     SOLVE THE SYSTEM OF EQUATIONS.
    !
    ierr1 = 0
    i1 = 1
    call genbunn(nbdcnd, nunk, i1, munk, W(1), W(id2+1), W(id3+1), &
        idimf, F(mstart, nstart), ierr1, W(id4+1))
    if (nbdcnd == 0) then
        f(mstart:mstop, np1) = F(mstart:mstop, 1)
    end if

end subroutine hwscylL

end module module_hwscyl
!
! REVISION HISTORY
!
! September 1973    Version 1
! April     1976    Version 2
! January   1978    Version 3
! December  1979    Version 3.1
! February  1985    Documentation upgrade
! November  1988    Version 3.2, FORTRAN 77 changes
! June      2004    Version 5.0, Fortran 90 Changes
!-----------------------------------------------------------------------
