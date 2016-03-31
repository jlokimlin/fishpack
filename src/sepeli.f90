module module_sepeli

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_blktri, only: &
        blktrii

    use module_sepaux, only: &
        seport, &
        sepmin, &
        septri, &
        sepdx, &
        sepdy, &
        kswx, kswy, k, l, mit, nit, is, ms, js, ns, & ! saved integer constants
        ait, bit, cit, dit, dlx, dly, tdlx3, tdly3, dlx4, dly4, & ! saved real constants
        get_coefficients

    ! Explicit typing only!
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: sepeli


contains


    subroutine sepeli(intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, &
        beta, c, d, n, nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, grhs, &
        usol, idmn, w, pertrb, ierror)
        !
        !     file sepeli.f
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
        !     SUBROUTINE sepeli (INTL, IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA, C,
        !    +                   D, N, NBDCND, BDC, GAMA, BDD, XNU, COFX, COFY, GRHS,
        !    +                   USOL, IDMN, W, PERTRB, ierror)
        !
        ! DIMENSION OF           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
        ! ARGUMENTS              USOL(IDMN, N+1), GRHS(IDMN, N+1),
        !
        ! LATEST REVISION        JUNE 2004
        !
        ! PURPOSE                sepeli SOLVES FOR EITHER THE SECOND-ORDER
        !                        FINITE DIFFERENCE APPROXIMATION OR A
        !                        FOURTH-ORDER APPROXIMATION TO A SEPARABLE
        !                        ELLIPTIC EQUATION
        !
        !                                 2    2
        !                          AF(X)*D U/DX + BF(X)*DU/DX  + CF(X)*U +
        !                                 2    2
        !                          DF(Y)*D U/DY  + EF(Y)*DU/DY + FF(Y)*U
        !
        !                          = G(X, Y)
        !
        !                        ON A RECTANGLE (X GREATER THAN OR EQUAL TO A
        !                        AND LESS THAN OR EQUAL TO B; Y GREATER THAN
        !                        OR EQUAL TO C AND LESS THAN OR EQUAL TO D).
        !                        ANY COMBINATION OF PERIODIC OR MIXED BOUNDARY
        !                        CONDITIONS IS ALLOWED.
        !
        !                        THE POSSIBLE BOUNDARY CONDITIONS ARE:
        !                        IN THE X-DIRECTION:
        !                        (0) PERIODIC, U(X+B-A, Y)=U(X, Y) FOR ALL
        !                            Y, X (1) U(A, Y), U(B, Y) ARE SPECIFIED FOR
        !                            ALL Y
        !                        (2) U(A, Y), DU(B, Y)/DX+BETA*U(B, Y) ARE
        !                            SPECIFIED FOR ALL Y
        !                        (3) DU(A, Y)/DX+ALPHA*U(A, Y), DU(B, Y)/DX+
        !                            BETA*U(B, Y) ARE SPECIFIED FOR ALL Y
        !                        (4) DU(A, Y)/DX+ALPHA*U(A, Y), U(B, Y) ARE
        !                            SPECIFIED FOR ALL Y
        !
        !                        IN THE Y-DIRECTION:
        !                        (0) PERIODIC, U(X, Y+D-C)=U(X, Y) FOR ALL X, Y
        !                        (1) U(X, C), U(X, D) ARE SPECIFIED FOR ALL X
        !                        (2) U(X, C), DU(X, D)/DY+XNU*U(X, D) ARE
        !                            SPECIFIED FOR ALL X
        !                        (3) DU(X, C)/DY+GAMA*U(X, C), DU(X, D)/DY+
        !                            XNU*U(X, D) ARE SPECIFIED FOR ALL X
        !                        (4) DU(X, C)/DY+GAMA*U(X, C), U(X, D) ARE
        !                            SPECIFIED FOR ALL X
        !
        ! USAGE                  CALL sepeli (INTL, IORDER, A, B, M, MBDCND, BDA,
        !                                     ALPHA, BDB, BETA, C, D, N, NBDCND, BDC,
        !                                     GAMA, BDD, XNU, COFX, COFY, GRHS, USOL,
        !                                     IDMN, W, PERTRB, ierror)
        !
        ! ARGUMENTS
        ! ON INPUT               INTL
        !                          = 0 ON INITIAL ENTRY TO sepeli OR IF ANY
        !                              OF THE ARGUMENTS C, D, N, NBDCND, COFY
        !                              ARE CHANGED FROM A PREVIOUS CALL
        !                          = 1 IF C, D, N, NBDCND, COFY ARE UNCHANGED
        !                              FROM THE PREVIOUS CALL.
        !
        !                        IORDER
        !                          = 2 IF A SECOND-ORDER APPROXIMATION
        !                              IS SOUGHT
        !                          = 4 IF A FOURTH-ORDER APPROXIMATION
        !                              IS SOUGHT
        !
        !                        A, B
        !                          THE RANGE OF THE X-INDEPENDENT VARIABLE,
        !                          I.E., X IS GREATER THAN OR EQUAL TO A
        !                          AND LESS THAN OR EQUAL TO B.  A MUST BE
        !                          LESS THAN B.
        !
        !                        M
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL [A, B] IS SUBDIVIDED. HENCE,
        !                          THERE WILL BE M+1 GRID POINTS IN THE X-
        !                          DIRECTION GIVEN BY XI=A+(I-1)*DLX
        !                          FOR I=1, 2, ..., M+1 WHERE DLX=(B-A)/M IS
        !                          THE PANEL WIDTH.  M MUST BE LESS THAN
        !                          IDMN AND GREATER THAN 5.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITION
        !                          AT X=A AND X=B
        !
        !                          = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.,
        !                              U(X+B-A, Y)=U(X, Y)  FOR ALL Y, X
        !                          = 1 IF THE SOLUTION IS SPECIFIED AT X=A
        !                              AND X=B, I.E., U(A, Y) AND U(B, Y) ARE
        !                              SPECIFIED FOR ALL Y
        !                          = 2 IF THE SOLUTION IS SPECIFIED AT X=A AND
        !                              THE BOUNDARY CONDITION IS MIXED AT X=B,
        !                              I.E., U(A, Y) AND DU(B, Y)/DX+BETA*U(B, Y)
        !                              ARE SPECIFIED FOR ALL Y
        !                          = 3 IF THE BOUNDARY CONDITIONS AT X=A AND
        !                              X=B ARE MIXED, I.E.,
        !                              DU(A, Y)/DX+ALPHA*U(A, Y) AND
        !                              DU(B, Y)/DX+BETA*U(B, Y) ARE SPECIFIED
        !                              FOR ALL Y
        !                          = 4 IF THE BOUNDARY CONDITION AT X=A IS
        !                              MIXED AND THE SOLUTION IS SPECIFIED
        !                              AT X=B, I.E., DU(A, Y)/DX+ALPHA*U(A, Y)
        !                              AND U(B, Y) ARE SPECIFIED FOR ALL Y
        !
        !                        BDA
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
        !                          THAT SPECIFIES THE VALueS OF
        !                          DU(A, Y)/DX+ ALPHA*U(A, Y) AT X=A, WHEN
        !                          MBDCND=3 OR 4.
        !                          BDA(J) = DU(A, YJ)/DX+ALPHA*U(A, YJ),
        !                          J=1, 2, ..., N+1. WHEN MBDCND HAS ANY OTHER
        !                          OTHER VALue, BDA IS A DUMMY PARAMETER.
        !
        !                        ALPHA
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN
        !                          CASE OF A MIXED BOUNDARY CONDITION AT X=A
        !                          (SEE ARGUMENT BDA).  IF MBDCND IS NOT
        !                          EQUAL TO 3 OR 4 THEN ALPHA IS A DUMMY
        !                          PARAMETER.
        !
        !                        BDB
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
        !                          THAT SPECIFIES THE VALueS OF
        !                          DU(B, Y)/DX+ BETA*U(B, Y) AT X=B.
        !                          WHEN MBDCND=2 OR 3
        !                          BDB(J) = DU(B, YJ)/DX+BETA*U(B, YJ),
        !                          J=1, 2, ..., N+1. WHEN MBDCND HAS ANY OTHER
        !                          OTHER VALue, BDB IS A DUMMY PARAMETER.
        !
        !                        BETA
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN
        !                          CASE OF A MIXED BOUNDARY CONDITION AT
        !                          X=B (SEE ARGUMENT BDB).  IF MBDCND IS
        !                          NOT EQUAL TO 2 OR 3 THEN BETA IS A DUMMY
        !                          PARAMETER.
        !
        !                        C, D
        !                          THE RANGE OF THE Y-INDEPENDENT VARIABLE,
        !                          I.E., Y IS GREATER THAN OR EQUAL TO C
        !                          AND LESS THAN OR EQUAL TO D.  C MUST BE
        !                          LESS THAN D.
        !
        !                        N
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL [C, D] IS SUBDIVIDED.
        !                          HENCE, THERE WILL BE N+1 GRID POINTS
        !                          IN THE Y-DIRECTION GIVEN BY
        !                          YJ=C+(J-1)*DLY FOR J=1, 2, ..., N+1 WHERE
        !                          DLY=(D-C)/N IS THE PANEL WIDTH.
        !                          IN ADDITION, N MUST BE GREATER THAN 4.
        !
        !                        NBDCND
        !                          INDICATES THE TYPES OF BOUNDARY CONDITIONS
        !                          AT Y=C AND Y=D
        !
        !                          = 0 IF THE SOLUTION IS PERIODIC IN Y,
        !                              I.E., U(X, Y+D-C)=U(X, Y)  FOR ALL X, Y
        !                          = 1 IF THE SOLUTION IS SPECIFIED AT Y=C
        !                              AND Y = D, I.E., U(X, C) AND U(X, D)
        !                              ARE SPECIFIED FOR ALL X
        !                          = 2 IF THE SOLUTION IS SPECIFIED AT Y=C
        !                              AND THE BOUNDARY CONDITION IS MIXED
        !                              AT Y=D, I.E., U(X, C) AND
        !                              DU(X, D)/DY+XNU*U(X, D) ARE SPECIFIED
        !                              FOR ALL X
        !                          = 3 IF THE BOUNDARY CONDITIONS ARE MIXED
        !                              AT Y=C AND Y=D, I.E.,
        !                              DU(X, D)/DY+GAMA*U(X, C) AND
        !                              DU(X, D)/DY+XNU*U(X, D) ARE SPECIFIED
        !                              FOR ALL X
        !                          = 4 IF THE BOUNDARY CONDITION IS MIXED
        !                              AT Y=C AND THE SOLUTION IS SPECIFIED
        !                              AT Y=D, I.E. DU(X, C)/DY+GAMA*U(X, C)
        !                              AND U(X, D) ARE SPECIFIED FOR ALL X
        !
        !                        BDC
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
        !                          THAT SPECIFIES THE VALue OF
        !                          DU(X, C)/DY+GAMA*U(X, C) AT Y=C.
        !                          WHEN NBDCND=3 OR 4 BDC(I) = DU(XI, C)/DY +
        !                          GAMA*U(XI, C), I=1, 2, ..., M+1.
        !                          WHEN NBDCND HAS ANY OTHER VALue, BDC
        !                          IS A DUMMY PARAMETER.
        !
        !                        GAMA
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN
        !                          CASE OF A MIXED BOUNDARY CONDITION AT
        !                          Y=C (SEE ARGUMENT BDC).  IF NBDCND IS
        !                          NOT EQUAL TO 3 OR 4 THEN GAMA IS A DUMMY
        !                          PARAMETER.
        !
        !                        BDD
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
        !                          THAT SPECIFIES THE VALue OF
        !                          DU(X, D)/DY + XNU*U(X, D) AT Y=C.
        !                          WHEN NBDCND=2 OR 3 BDD(I) = DU(XI, D)/DY +
        !                          XNU*U(XI, D), I=1, 2, ..., M+1.
        !                          WHEN NBDCND HAS ANY OTHER VALue, BDD
        !                          IS A DUMMY PARAMETER.
        !
        !                        XNU
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN
        !                          CASE OF A MIXED BOUNDARY CONDITION AT
        !                          Y=D (SEE ARGUMENT BDD).  IF NBDCND IS
        !                          NOT EQUAL TO 2 OR 3 THEN XNU IS A
        !                          DUMMY PARAMETER.
        !
        !                        COFX
        !                          A USER-SUPPLIED SUBPROGRAM WITH
        !                          PARAMETERS X, AFUN, BFUN, CFUN WHICH
        !                          RETURNS THE VALueS OF THE X-DEPENDENT
        !                          COEFFICIENTS AF(X), BF(X), CF(X) IN THE
        !                          ELLIPTIC EQUATION AT X.
        !
        !                        COFY
        !                          A USER-SUPPLIED SUBPROGRAM WITH PARAMETERS
        !                          Y, DFUN, EFUN, FFUN WHICH RETURNS THE
        !                          VALueS OF THE Y-DEPENDENT COEFFICIENTS
        !                          DF(Y), EF(Y), FF(Y) IN THE ELLIPTIC
        !                          EQUATION AT Y.
        !
        !                          NOTE:  COFX AND COFY MUST BE DECLARED
        !                          EXTERNAL IN THE CALLING ROUTINE.
        !                          THE VALueS RETURNED IN AFUN AND DFUN
        !                          MUST SATISFY AFUN*DFUN GREATER THAN 0
        !                          FOR A LESS THAN X LESS THAN B, C LESS
        !                          THAN Y LESS THAN D (SEE ierror=10).
        !                          THE COEFFICIENTS PROVIDED MAY LEAD TO A
        !                          MATRIX EQUATION WHICH IS NOT DIAGONALLY
        !                          DOMINANT IN WHICH CASE SOLUTION MAY FAIL
        !                          (SEE ierror=4).
        !
        !                        GRHS
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALueS OF THE RIGHT-HAND SIDE OF THE
        !                          ELLIPTIC EQUATION, I.E.,
        !                          GRHS(I, J)=G(XI, YI), FOR I=2, ..., M,
        !                          J=2, ..., N.  AT THE BOUNDARIES, GRHS IS
        !                          DEFINED BY
        !
        !                          MBDCND   GRHS(1, J)   GRHS(M+1, J)
        !                          ------   ---------   -----------
        !                            0      G(A, YJ)     G(B, YJ)
        !                            1         *           *
        !                            2         *        G(B, YJ)  J=1, 2, ..., N+1
        !                            3      G(A, YJ)     G(B, YJ)
        !                            4      G(A, YJ)        *
        !
        !                          NBDCND   GRHS(I, 1)   GRHS(I, N+1)
        !                          ------   ---------   -----------
        !                            0      G(XI, C)     G(XI, D)
        !                            1         *           *
        !                            2         *        G(XI, D)  I=1, 2, ..., M+1
        !                            3      G(XI, C)     G(XI, D)
        !                            4      G(XI, C)        *
        !
        !                          WHERE * MEANS THESE QUANTITIES ARE NOT USED.
        !                          GRHS SHOULD BE DIMENSIONED IDMN BY AT LEAST
        !                          N+1 IN THE CALLING ROUTINE.
        !
        !                        USOL
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALueS OF THE SOLUTION ALONG THE BOUNDARIES.
        !                          AT THE BOUNDARIES, USOL IS DEFINED BY
        !
        !                          MBDCND   USOL(1, J)   USOL(M+1, J)
        !                          ------   ---------   -----------
        !                            0         *           *
        !                            1      U(A, YJ)     U(B, YJ)
        !                            2      U(A, YJ)        *     J=1, 2, ..., N+1
        !                            3         *           *
        !                            4         *        U(B, YJ)
        !
        !                          NBDCND   USOL(I, 1)   USOL(I, N+1)
        !                          ------   ---------   -----------
        !                            0         *           *
        !                            1      U(XI, C)     U(XI, D)
        !                            2      U(XI, C)        *     I=1, 2, ..., M+1
        !                            3         *           *
        !                            4         *        U(XI, D)
        !
        !                          WHERE * MEANS THE QUANTITIES ARE NOT USED
        !                          IN THE SOLUTION.
        !
        !                          IF IORDER=2, THE USER MAY EQUIVALENCE GRHS
        !                          AND USOL TO SAVE SPACE.  NOTE THAT IN THIS
        !                          CASE THE TABLES SPECIFYING THE BOUNDARIES
        !                          OF THE GRHS AND USOL ARRAYS DETERMINE THE
        !                          BOUNDARIES UNIQueLY EXCEPT AT THE CORNERS.
        !                          IF THE TABLES CALL FOR BOTH G(X, Y) AND
        !                          U(X, Y) AT A CORNER THEN THE SOLUTION MUST
        !                          BE CHOSEN.  FOR EXAMPLE, IF MBDCND=2 AND
        !                          NBDCND=4, THEN U(A, C), U(A, D), U(B, D) MUST
        !                          BE CHOSEN AT THE CORNERS IN ADDITION
        !                          TO G(B, C).
        !
        !                          IF IORDER=4, THEN THE TWO ARRAYS, USOL AND
        !                          GRHS, MUST BE DISTINCT.
        !
        !                          USOL SHOULD BE DIMENSIONED IDMN BY AT LEAST
        !                          N+1 IN THE CALLING ROUTINE.
        !
        !                        IDMN
        !                          THE ROW (OR FIRST) DIMENSION OF THE ARRAYS
        !                          GRHS AND USOL AS IT APPEARS IN THE PROGRAM
        !                          CALLING sepeli.  THIS PARAMETER IS USED
        !                          TO SPECIFY THE VARIABLE DIMENSION OF GRHS
        !                          AND USOL.  IDMN MUST BE AT LEAST 7 AND
        !                          GREATER THAN OR EQUAL TO M+1.
        !
        !                        W
        !                          A fortran 90 derived TYPE (FishpackWorkspace) variable
        !                          that must be declared by the user.  The first
        !                          two declarative statements in the user program
        !                          calling sepeli must be:
        !
        !                               use type_FishpackWorkspace
        !                               TYPE (FishpackWorkspace) :: W
        !
        !                          The first statement makes the fishpack module
        !                          defined in the file "fish.f" available to the
        !                          user program calling sepeli.  The second statement
        !                          declares a derived type variable (defined in
        !                          the module "fish.f") which is used internally
        !                          in sepeli to dynamically allocate real and complex
        !                          work space used in solution.  An error flag
        !                          (ierror = 20) is set if the required work space
        !                          allocation fails (for example if N, M are too large)
        !                          Real and complex values are set in the components
        !                          of W on a initial (INTL=0) call to sepeli.  These
        !                          must be preserved on non-initial calls (INTL=1)
        !                          to sepeli.  This eliminates redundant calculations
        !                          and saves compute time.
        !               ****       IMPORTANT!  The user program calling sepeli should
        !                          include the statement:
        !
        !                               CALL FISHFIN(W)
        !
        !                          after the final approximation is generated by
        !                          sepeli.  The will deallocate the real and complex
        !                          work space of W.  Failure to include this statement
        !                          could result in serious memory leakage.
        !
        ! ON OUTPUT              USOL
        !                          CONTAINS THE APPROXIMATE SOLUTION TO THE
        !                          ELLIPTIC EQUATION.
        !                          USOL(I, J) IS THE APPROXIMATION TO U(XI, YJ)
        !                          FOR I=1, 2..., M+1 AND J=1, 2, ..., N+1.
        !                          THE APPROXIMATION HAS ERROR
        !                          O(DLX**2+DLY**2) IF CALLED WITH IORDER=2
        !                          AND O(DLX**4+DLY**4) IF CALLED WITH
        !                          IORDER=4.
        !
        !                        W
        !                          The derived type (FishpackWorkspace) variable W
        !                          contains real and complex values that must not
        !                          be destroyed if sepeli is called again with
        !                          INTL=1.
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
        !                          BOUNDARY CONDITIONS
        !                          (I.E., ALPHA=BETA=0 IF MBDCND=3;
        !                          GAMA=XNU=0 IF NBDCND=3) IS SPECIFIED
        !                          AND IF THE COEFFICIENTS OF U(X, Y) IN THE
        !                          SEPARABLE ELLIPTIC EQUATION ARE ZERO
        !                          (I.E., CF(X)=0 FOR X GREATER THAN OR EQUAL
        !                          TO A AND LESS THAN OR EQUAL TO B;
        !                          FF(Y)=0 FOR Y GREATER THAN OR EQUAL TO C
        !                          AND LESS THAN OR EQUAL TO D) THEN A
        !                          SOLUTION MAY NOT EXIST.  PERTRB IS A
        !                          CONSTANT CALCULATED AND SUBTRACTED FROM
        !                          THE RIGHT-HAND SIDE OF THE MATRIX EQUATIONS
        !                          GENERATED BY sepeli WHICH INSURES THAT A
        !                          SOLUTION EXISTS. sepeli THEN COMPUTES THIS
        !                          SOLUTION WHICH IS A WEIGHTED MINIMAL LEAST
        !                          SQUARES SOLUTION TO THE ORIGINAL PROBLEM.
        !
        !                        ierror
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS OR FAILURE TO FIND A SOLUTION
        !                          = 0 NO ERROR
        !                          = 1 IF A GREATER THAN B OR C GREATER THAN D
        !                          = 2 IF MBDCND LESS THAN 0 OR MBDCND GREATER
        !                              THAN 4
        !                          = 3 IF NBDCND LESS THAN 0 OR NBDCND GREATER
        !                              THAN 4
        !                          = 4 IF ATTEMPT TO FIND A SOLUTION FAILS.
        !                              (THE LINEAR SYSTEM GENERATED IS NOT
        !                              DIAGONALLY DOMINANT.)
        !                          = 5 IF IDMN IS TOO SMALL
        !                              (SEE DISCUSSION OF IDMN)
        !                          = 6 IF M IS TOO SMALL OR TOO LARGE
        !                              (SEE DISCUSSION OF M)
        !                          = 7 IF N IS TOO SMALL (SEE DISCUSSION OF N)
        !                          = 8 IF IORDER IS NOT 2 OR 4
        !                          = 9 IF INTL IS NOT 0 OR 1
        !                          = 10 IF AFUN*DFUN LESS THAN OR EQUAL TO 0
        !                               FOR SOME INTERIOR MESH POINT (XI, YJ)
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space in the derived type
        !                               (FishpackWorkspace) variable W fails (e.g.,
        !                               if N, M are too large for the platform used)
        !
        !                          NOTE (CONCERNING ierror=4):  FOR THE
        !                          COEFFICIENTS INPUT THROUGH COFX, COFY,
        !                          THE DISCRETIZATION MAY LEAD TO A BLOCK
        !                          TRIDIAGONAL LINEAR SYSTEM WHICH IS NOT
        !                          DIAGONALLY DOMINANT (FOR EXAMPLE, THIS
        !                          HAPPENS IF CFUN=0 AND BFUN/(2.*DLX) GREATER
        !                          THAN AFUN/DLX**2).  IN THIS CASE SOLUTION
        !                          MAY FAIL.  THIS CANNOT HAPPEN IN THE LIMIT
        !                          AS DLX, DLY APPROACH ZERO.  HENCE, THE
        !                          CONDITION MAY BE REMEDIED BY TAKING LARGER
        !                          VALueS FOR M OR N.
        !
        ! SPECIAL CONDITIONS     SEE COFX, COFY ARGUMENT DESCRIPTIONS ABOVE.
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED FILES         blktri.f, comf.f, sepaux.f, fish.f
        !
        ! LANGUAGE               Fortran 90
        !
        ! HISTORY                DEVELOPED AT NCAR DURING 1975-76 BY
        !                        JOHN C. ADAMS OF THE SCIENTIFIC COMPUTING
        !                        DIVISION.  RELEASED ON NCAR'S PUBLIC SOFTWARE
        !                        LIBRARIES IN JANUARY 1980. Revised in June
        !                        2004 using Fortan 90 dynamically allocated work
        !                        space and derived data types to eliminate mixed
        !                        mode conflicts in the earlier versions. All
        !                        statement labels, arithmetic if statements and
        !                        computed GO TO statements have been removed from
        !                        the current version of sepeli.
        !
        ! ALGORITHM              sepeli AUTOMATICALLY DISCRETIZES THE
        !                        SEPARABLE ELLIPTIC EQUATION WHICH IS THEN
        !                        SOLVED BY A GENERALIZED CYCLIC REDUCTION
        !                        ALGORITHM IN THE SUBROUTINE, blktri.  THE
        !                        FOURTH-ORDER SOLUTION IS OBTAINED USING
        !                        'deferRED CORRECTIONS' WHICH IS DESCRIBED
        !                        AND REFERENCED IN SECTIONS, REFERENCES AND
        !                        METHOD.
        !
        ! TIMING                 THE OPERATIONAL COUNT IS PROPORTIONAL TO
        !                        M*N*LOG2(N).
        !
        ! ACCURACY               THE FOLLOWING ACCURACY RESULTS WERE OBTAINED
        !                        using 64 bit floating point arithmetic.  Note
        !                        THAT THE FOURTH-ORDER accuracy is not realized
        !                        UNTIL THE MESH IS sufficiently refined.
        !
        !                                     SECOND-ORDER  FOURTH-ORDER
        !                            M    N     ERROR         ERROR
        !
        !                             6    6    6.8E-1        1.2E0
        !                            14   14    1.4E-1        1.8E-1
        !                            30   30    3.2E-2        9.7E-3
        !                            62   62    7.5E-3        3.0E-4
        !                           126  126    1.8E-3        3.5E-6
        !
        !
        ! REFERENCES             KELLER, H.B., NUMERICAL METHODS FOR TWO-POINT
        !                        BOUNDARY-VALue PROBLEMS, BLAISDEL (1968),
        !                        WALTHAM, MASS.
        !
        !                        SWARZTRAUBER, P., AND R. SWEET (1975):
        !                        EFFICIENT FORTRAN SUBPROGRAMS FOR THE
        !                        SOLUTION OF ELLIPTIC PARTIAL DIFFERENTIAL
        !                        EQUATIONS.  NCAR TECHNICAL NOTE
        !                        NCAR-TN/IA-109, PP. 135-137.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip),              intent (in)     :: intl
        integer (ip),              intent (in)     :: iorder
        integer (ip),              intent (in)     :: m
        integer (ip),              intent (in)     :: mbdcnd
        integer (ip),              intent (in)     :: n
        integer (ip),              intent (in)     :: nbdcnd
        integer (ip),              intent (in)     :: idmn
        integer (ip),              intent (out)    :: ierror
        real (wp),                 intent (in)     :: a
        real (wp),                 intent (in)     :: b
        real (wp),                 intent (in)     :: alpha
        real (wp),                 intent (in)     :: beta
        real (wp),                 intent (in)     :: c
        real (wp),                 intent (in)     :: d
        real (wp),                 intent (in)     :: gama
        real (wp),                 intent (in)     :: xnu
        real (wp),                 intent (out)    :: pertrb
        real (wp), contiguous,     intent (in)     :: bda(:)
        real (wp), contiguous,     intent (in)     :: bdb(:)
        real (wp), contiguous,     intent (in)     :: bdc(:)
        real (wp), contiguous,     intent (in)     :: bdd(:)
        real (wp)                                  :: grhs(idmn, *)
        real (wp)                                  :: usol(idmn, *)
        class (FishpackWorkspace)                  :: w
        procedure (get_coefficients)               :: cofx
        procedure (get_coefficients)               :: cofy
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)       :: k, l, np, irwk, icwk
        integer (ip), save :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12
        !-----------------------------------------------

        !     save local variable work space pointers for noninitial call
        !     check input arguments
        call chkprm (intl, iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, cofx, &
            cofy, idmn, ierror)

        if (ierror /= 0) return
        if (intl == 0) then
            !     allocate space and set work space indices on initial call only
            k = m + 1
            l = n + 1
            !          compute required blktri work space lengths
            np = nbdcnd

            call w%get_block_tridiagonal_workpace_dimensions (n, m, irwk, icwk)
            !
            !     SET WORK SPACE INDICES
            !
            i1 = irwk + 1
            i2 = i1 + l
            i3 = i2 + l
            i4 = i3 + l
            i5 = i4 + l
            i6 = i5 + l
            i7 = i6 + l
            i8 = i7 + k
            i9 = i8 + k
            i10 = i9 + k
            i11 = i10 + k
            i12 = i11 + k
            !          set sepeli work space requirements
            irwk = i12 + k
            icwk = icwk + 3*(m + 1)
            !          allocate required real and complex work space
            call w%create( irwk, icwk, ierror )
            !          return if allocation failure
            if (ierror == 20) return
        end if
        ierror = 0
        !     compute second or fourth order solution
        associate( &
            rew => w%real_workspace, &
            cxw => w%complex_workspace &
            )
            call spelip(intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, d, n, &
                nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, rew(i1), rew(i2), rew(i3), &
                rew(i4), rew(i5), rew(i6), rew(i7), rew(i8), rew(i9), &
                rew(i10), rew(i11), rew(i12), grhs, usol, idmn, rew, cxw, &
                pertrb, ierror)
        end associate

    end subroutine sepeli


    subroutine spelip( intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, &
        beta, c, d, n, nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, an, bn, &
        cn, dn, un, zn, am, bm, cm, dm, um, zm, grhs, usol, idmn, w, &
        wc, pertrb, ierror )
        !
        ! Purpose:
        !
        !     spelip sets up vectors and arrays for input to blktri
        !     and computes a second order solution in usol.  a return jump to
        !     sepeli occurrs if iorder=2.  if iorder=4 a fourth order
        !     solution is generated in usol.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: intl
        integer (ip), intent (in)     :: iorder
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: mbdcnd
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: nbdcnd
        integer (ip), intent (in)     :: idmn
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in)     :: a
        real (wp),    intent (in)     :: b
        real (wp),    intent (in)     :: alpha
        real (wp),    intent (in)     :: beta
        real (wp),    intent (in)     :: c
        real (wp),    intent (in)     :: d
        real (wp),    intent (in)     :: gama
        real (wp),    intent (in)     :: xnu
        real (wp),    intent (out)    :: pertrb
        real (wp),    intent (in)     :: bda(*)
        real (wp),    intent (in)     :: bdb(*)
        real (wp),    intent (in)     :: bdc(*)
        real (wp),    intent (in)     :: bdd(*)
        real (wp),    intent (in out) :: an(*)
        real (wp),    intent (in out) :: bn(*)
        real (wp),    intent (in out) :: cn(*)
        real (wp),    intent (in out) :: dn(*)
        real (wp),    intent (in out) :: un(*)
        real (wp),    intent (in out) :: zn(*)
        real (wp),    intent (in out) :: am(*)
        real (wp),    intent (in out) :: bm(*)
        real (wp),    intent (in out) :: cm(*)
        real (wp),    intent (in out) :: dm(*)
        real (wp),    intent (in out) :: um(*)
        real (wp),    intent (in out) :: zm(*)
        real (wp),    intent (in out) :: grhs(idmn, *)
        real (wp),    intent (in out) :: usol(idmn, *)
        real (wp),    intent (in out) :: w(*)
        complex (wp)                  :: wc(*)
        procedure (get_coefficients)  :: cofx
        procedure (get_coefficients)  :: cofy
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: i, j, i1, mp, np
        real (wp)    :: xi, ai, bi, ci, axi, bxi, cxi
        real (wp)    :: yj, dj, ej, fj, dyj, eyj
        real (wp)    :: fyj, ax1, cxm, dy1, fyn, prtrb
        logical      :: singlr
        !--------------------------------------------------------------------------------

        !     set parameters internally
        !
        kswx = mbdcnd + 1
        kswy = nbdcnd + 1
        k = m + 1
        l = n + 1
        ait = a
        bit = b
        cit = c
        dit = d
        !
        !     set right hand side values from grhs in usol on the interior
        !     and non-specified boundaries.
        !
        usol(2:m, 2:n) = grhs(2:m, 2:n)
        if (kswx/=2 .and. kswx/=3) then
            usol(1, 2:n) = grhs(1, 2:n)
        end if
        if (kswx/=2 .and. kswx/=5) then
            usol(k, 2:n) = grhs(k, 2:n)
        end if
        if (kswy/=2 .and. kswy/=3) then
            usol(2:m, 1) = grhs(2:m, 1)
        end if
        if (kswy/=2 .and. kswy/=5) then
            usol(2:m, l) = grhs(2:m, l)
        end if
        if (kswx/=2 .and. kswx/=3 .and. kswy/=2 .and. kswy/=3) usol(1, 1) &
            = grhs(1, 1)
        if (kswx/=2 .and. kswx/=5 .and. kswy/=2 .and. kswy/=3) usol(k, 1) &
            = grhs(k, 1)
        if (kswx/=2 .and. kswx/=3 .and. kswy/=2 .and. kswy/=5) usol(1, l) &
            = grhs(1, l)
        if (kswx/=2 .and. kswx/=5 .and. kswy/=2 .and. kswy/=5) usol(k, l) &
            = grhs(k, l)
        i1 = 1
        !
        !     set switches for periodic or non-periodic boundaries
        !
        mp = 1
        np = 1
        if (kswx == 1) mp = 0
        if (kswy == 1) np = 0
        !
        !     set dlx, dly and size of block tri-diagonal system generated
        !     in nint, mint
        !
        dlx = (bit - ait)/real(m)
        mit = k - 1
        if (kswx == 2) mit = k - 2
        if (kswx == 4) mit = k
        dly = (dit - cit)/real(n)
        nit = l - 1
        if (kswy == 2) nit = l - 2
        if (kswy == 4) nit = l
        tdlx3 = 2.0_wp * dlx**3
        dlx4 = dlx**4
        tdly3 = 2.0_wp * dly**3
        dly4 = dly**4
        !
        !     set subscript limits for portion of array to input to blktri
        !
        is = 1
        js = 1
        if (kswx==2 .or. kswx==3) is = 2
        if (kswy==2 .or. kswy==3) js = 2
        ns = nit + js - 1
        ms = mit + is - 1
        !
        !     set x - direction
        !
        do i = 1, mit
            xi = ait + real(is + i - 2)*dlx
            call cofx(xi, ai, bi, ci)
            axi = (ai/dlx - 0.5*bi)/dlx
            bxi = (-2.*ai/dlx**2) + ci
            cxi = (ai/dlx + 0.5*bi)/dlx
            am(i) = axi
            bm(i) = bxi
            cm(i) = cxi
        end do
        !
        !     set y direction
        !
        do j = 1, nit
            yj = cit + real(js + j - 2)*dly
            call cofy (yj, dj, ej, fj)
            dyj = (dj/dly - 0.5*ej)/dly
            eyj = (-2.*dj/dly**2) + fj
            fyj = (dj/dly + 0.5*ej)/dly
            an(j) = dyj
            bn(j) = eyj
            cn(j) = fyj
        end do
        !
        !     adjust edges in x direction unless periodic
        !
        ax1 = am(1)
        cxm = cm(mit)
        select case (kswx)
            case (2)
                !
                !     dirichlet-dirichlet in x direction
                !
                am(1) = 0.0
                cm(mit) = 0.0
            case (5)
                !
                !     mixed-dirichlet in x direction
                !
                am(1) = 0.0
                bm(1) = bm(1) + 2.*alpha*dlx*ax1
                cm(1) = cm(1) + ax1
                cm(mit) = 0.0
            case (3)
                !
                !     dirichlet-mixed in x direction
                !
                am(1) = 0.0
                am(mit) = am(mit) + cxm
                bm(mit) = bm(mit) - 2.*beta*dlx*cxm
                cm(mit) = 0.0
            !
            !     mixed - mixed in x direction
            !
            case (4)
                am(1) = 0.0
                bm(1) = bm(1) + 2.*dlx*alpha*ax1
                cm(1) = cm(1) + ax1
                am(mit) = am(mit) + cxm
                bm(mit) = bm(mit) - 2.*dlx*beta*cxm
                cm(mit) = 0.0
        end select
        !
        !     adjust in y direction unless periodic
        !
        dy1 = an(1)
        fyn = cn(nit)
        select case (kswy)
            case (2)
                !
                !     dirichlet-dirichlet in y direction
                !
                an(1) = 0.0
                cn(nit) = 0.0
            case (5)
                !
                !     mixed-dirichlet in y direction
                !
                an(1) = 0.0
                bn(1) = bn(1) + 2.*dly*gama*dy1
                cn(1) = cn(1) + dy1
                cn(nit) = 0.0
            case (3)
                !
                !     dirichlet-mixed in y direction
                !
                an(1) = 0.0
                an(nit) = an(nit) + fyn
                bn(nit) = bn(nit) - 2.*dly*xnu*fyn
                cn(nit) = 0.0
            case (4)
                !
                !     mixed - mixed direction in y direction
                !
                an(1) = 0.0
                bn(1) = bn(1) + 2.*dly*gama*dy1
                cn(1) = cn(1) + dy1
                an(nit) = an(nit) + fyn
                bn(nit) = bn(nit) - 2.0_wp * dly*xnu*fyn
                cn(nit) = 0.0
        end select
        if (kswx /= 1) then
            !
            !     adjust usol along x edge
            !
            if (kswx==2 .or. kswx==3) then
                if (kswx==2 .or. kswx==5) then
                    usol(is, js:ns) = usol(is, js:ns) - ax1*usol(1, js:ns)
                    usol(ms, js:ns) = usol(ms, js:ns) - cxm*usol(k, js:ns)
                else
                    usol(is, js:ns) = usol(is, js:ns) - ax1*usol(1, js:ns)
                    usol(ms, js:ns) = usol(ms, js:ns) - 2.0_wp * dlx*cxm*bdb(js:ns)
                end if
            else
                if (kswx==2 .or. kswx==5) then
                    usol(is, js:ns) = usol(is, js:ns) + 2.0_wp * dlx*ax1*bda(js:ns)
                    usol(ms, js:ns) = usol(ms, js:ns) - cxm*usol(k, js:ns)
                else
                    usol(is, js:ns) = usol(is, js:ns) + 2.0_wp * dlx*ax1*bda(js:ns)
                    usol(ms, js:ns) = usol(ms, js:ns) - 2.0_wp * dlx*cxm*bdb(js:ns)
                end if
            end if
        end if
        if (kswy /= 1) then
            !
            !     adjust usol along y edge
            !
            if (kswy==2 .or. kswy==3) then
                if (kswy==2 .or. kswy==5) then
                    usol(is:ms, js) = usol(is:ms, js) - dy1*usol(is:ms, 1)
                    usol(is:ms, ns) = usol(is:ms, ns) - fyn*usol(is:ms, l)
                else
                    usol(is:ms, js) = usol(is:ms, js) - dy1*usol(is:ms, 1)
                    usol(is:ms, ns) = usol(is:ms, ns) - 2.0_wp * dly*fyn*bdd(is:ms)
                end if
            else
                if (kswy==2 .or. kswy==5) then
                    usol(is:ms, js) = usol(is:ms, js) + 2.0_wp * dly*dy1*bdc(is:ms)
                    usol(is:ms, ns) = usol(is:ms, ns) - fyn*usol(is:ms, l)
                else
                    usol(is:ms, js) = usol(is:ms, js) + 2.0_wp * dly*dy1*bdc(is:ms)
                    usol(is:ms, ns) = usol(is:ms, ns) - 2.0_wp * dly*fyn*bdd(is:ms)
                end if
            end if
        end if
        !
        !     save adjusted edges in grhs if iorder=4
        !
        if (iorder == 4) then
            grhs(is, js:ns) = usol(is, js:ns)
            grhs(ms, js:ns) = usol(ms, js:ns)
            grhs(is:ms, js) = usol(is:ms, js)
            grhs(is:ms, ns) = usol(is:ms, ns)
        end if
        pertrb = 0.0
        !
        !     check if operator is singular
        !
        call chksng(mbdcnd, nbdcnd, alpha, beta, gama, xnu, cofx, cofy, singlr)
        !
        !     compute non-zero eigenvector in null space of transpose
        !     if singular
        !
        if (singlr) call septri (mit, am, bm, cm, dm, um, zm)
        if (singlr) call septri (nit, an, bn, cn, dn, un, zn)
        !
        !     make initialization call to blktrii
        !
        if (intl == 0) then
            call blktrii (intl, np, nit, an, bn, cn, mp, mit, am, bm, cm, &
                idmn, usol(is, js), ierror, w, wc)
            if (ierror /= 0) return
        end if
        !
        !     adjust right hand side if necessary
        !
        if (singlr) call seport (usol, idmn, zn, zm, pertrb)
        !
        !     compute solution
        !
        call blktrii (i1, np, nit, an, bn, cn, mp, mit, am, bm, cm, idmn, &
            usol(is, js), ierror, w, wc)
        if (ierror /= 0) return
        !
        !     set periodic boundaries if necessary
        !
        if (kswx == 1) then
            usol(k, :l) = usol(1, :l)
        end if
        if (kswy == 1) then
            usol(:k, l) = usol(:k, 1)
        end if
        !
        !     minimize solution with respect to weighted least squares
        !     norm if operator is singular
        !
        if (singlr) call sepmin (usol, idmn, zn, zm, prtrb)
        !
        !     return if deferred corrections and a fourth order solution are
        !     not flagged
        !
        if (iorder == 2) return
        !
        !     compute new right hand side for fourth order solution
        !
        call defer (cofx, cofy, idmn, usol, grhs)
        if (singlr) call seport (usol, idmn, zn, zm, pertrb)
        !
        !     compute fourth order solution
        !
        call blktrii (i1, np, nit, an, bn, cn, mp, mit, am, bm, cm, idmn, &
            usol(is, js), ierror, w, wc)
        if (ierror /= 0) return
        !
        !     set periodic boundaries if necessary
        !
        if (kswx == 1) then
            usol(k, :l) = usol(1, :l)
        end if
        if (kswy == 1) then
            usol(:k, l) = usol(:k, 1)
        end if
        !
        !     minimize solution with respect to weighted least squares
        !     norm if operator is singular
        !
        if (singlr) call sepmin (usol, idmn, zn, zm, prtrb)

    end subroutine spelip


    subroutine chkprm( intl, iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, &
        cofx, cofy, idmn, ierror )
        !
        ! Purpose:
        !
        ! This program checks the input arguments for errors
        !
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)    :: intl
        integer (ip), intent (in)    :: iorder
        integer (ip), intent (in)    :: m
        integer (ip), intent (in)    :: mbdcnd
        integer (ip), intent (in)    :: n
        integer (ip), intent (in)    :: nbdcnd
        integer (ip), intent (in)    :: idmn
        integer (ip), intent (out)   :: ierror
        real (wp),    intent (in)    :: a
        real (wp),    intent (in)    :: b
        real (wp),    intent (in)    :: c
        real (wp),    intent (in)    :: d
        procedure (get_coefficients) :: cofx
        procedure (get_coefficients) :: cofy
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip) :: i, j
        real (wp)    :: dlx, dly, xi, ai, bi, ci, yj, dj, ej, fj
        !-----------------------------------------------


        !     check definition of solution region
        !
        if (a>=b .or. c>=d) then
            ierror = 1
            return
        end if
        !
        !     check boundary condition arguments
        !
        if (mbdcnd<0 .or. mbdcnd>4) then
            ierror = 2
            return
        end if
        if (nbdcnd<0 .or. nbdcnd>4) then
            ierror = 3
            return
        end if
        !
        !     check first dimension in calling routine
        !
        if (idmn < 7) then
            ierror = 5
            return
        end if
        !
        !     check m, n
        !
        if (m>idmn - 1 .or. m<6) then
            ierror = 6
            return
        end if
        if (n < 5) then
            ierror = 7
            return
        end if
        !
        !     check iorder
        !
        if (iorder/=2 .and. iorder/=4) then
            ierror = 8
            return
        end if
        !
        !     check intl
        !
        if (intl/=0 .and. intl/=1) then
            ierror = 9
            return
        end if
        !
        !     check that equation is elliptic (only on initial call)
        !
        if (intl == 0) then
            dlx = (b - a)/real(m)
            dly = (d - c)/real(n)
            do i = 2, m
                xi = a + real(i - 1)*dlx
                call cofx(xi, ai, bi, ci)
                do j = 2, n
                    yj = c + real(j - 1)*dly
                    call cofy (yj, dj, ej, fj)
                    if (ai*dj > 0.0) cycle
                    ierror = 10
                    return
                end do
            end do
        end if
        !
        !     no error found
        !
        ierror = 0

    end subroutine chkprm


    subroutine chksng( mbdcnd, nbdcnd, alpha, beta, gama, xnu, cofx, &
        cofy, singlr )
        !
        ! Purpose:
        !
        !     this subroutine checks if the pde   sepeli
        !     must solve is a singular operator
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: mbdcnd
        integer (ip), intent (in)     :: nbdcnd
        real (wp),    intent (in)     :: alpha
        real (wp),    intent (in)     :: beta
        real (wp),    intent (in)     :: gama
        real (wp),    intent (in)     :: xnu
        logical ,     intent (out)    :: singlr
        procedure (get_coefficients)  :: cofx
        procedure (get_coefficients)  :: cofy
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: i, j
        real :: xi, ai, bi, ci, yj, dj, ej, fj
        !-----------------------------------------------


        ! Initialize flag
        singlr = .false.
        !
        !     check if the boundary conditions are
        !     entirely periodic and/or mixed
        !
        if(mbdcnd/=0.and.mbdcnd/=3.or.nbdcnd/=0.and.nbdcnd/=3)return
        !
        !     check that mixed conditions are pure neuman
        !
        if (mbdcnd == 3) then
            if (alpha/=0.0 .or. beta/=0.0) return
        end if

        if (nbdcnd == 3) then
            if (gama/=0.0 .or. xnu/=0.0) return
        end if
        !
        !     check that non-derivative coefficient functions
        !     are zero
        !
        do i = is, ms
            xi = ait + real(i - 1)*dlx
            call cofx(xi, ai, bi, ci)
            if (ci == 0.0) cycle
            return
        end do

        do j = js, ns
            yj = cit + real(j - 1)*dly
            call cofy (yj, dj, ej, fj)
            if (fj == 0.0) cycle
            return
        end do
        !
        !     the operator must be singular if this point is reached
        !
        singlr = .true.

    end subroutine chksng


    subroutine defer( cofx, cofy, idmn, usol, grhs )
        !
        ! Purpose:
        !
        !     this subroutine first approximates the truncation error given by
        !     trun1(x, y)=dlx**2*tx+dly**2*ty where
        !     tx=afun(x)*uxxxx/12.0+bfun(x)*uxxx/6.0 on the interior and
        !     at the boundaries if periodic(here uxxx, uxxxx are the third
        !     and fourth partial derivatives of u with respect to x).
        !     tx is of the form afun(x)/3.0_wp * (uxxxx/4.0+uxxx/dlx)
        !     at x=a or x=b if the boundary condition there is mixed.
        !     tx=0.0 along specified boundaries.  ty has symmetric form
        !     in y with x, afun(x), bfun(x) replaced by y, dfun(y), efun(y).
        !     the second order solution in usol is used to approximate
        !     (via second order finite differencing) the truncation error
        !     and the result is added to the right hand side in grhs
        !     and then transferred to usol to be used as a new right
        !     hand side when calling blktri for a fourth order solution.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: idmn
        real (wp),    intent (in out) :: usol(idmn, *)
        real (wp),    intent (in out) :: grhs(idmn, *)
        procedure (get_coefficients)  :: cofx
        procedure (get_coefficients)  :: cofy
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: j, i
        real (wp)    :: yj, dj, ej, fj, xi, ai, bi, ci
        real (wp)    :: uxxx, uxxxx, uyyy, uyyyy, tx, ty
        !-----------------------------------------------

        !
        !     compute truncation error approximation over the entire mesh
        !
        do j = js, ns
            yj = cit + real(j - 1)*dly
            call cofy (yj, dj, ej, fj)
            do i = is, ms
                xi = ait + real(i - 1)*dlx
                call cofx(xi, ai, bi, ci)
                !
                !     compute partial derivative approximations at (xi, yj)
                !
                call sepdx (usol, idmn, i, j, uxxx, uxxxx)
                call sepdy (usol, idmn, i, j, uyyy, uyyyy)
                tx = ai*uxxxx/12.0 + bi*uxxx/6.0
                ty = dj*uyyyy/12.0 + ej*uyyy/6.0
                !
                !     reset form of truncation if at boundary which is non-periodic
                !
                if (kswx/=1 .and. (i==1 .or. i==k)) tx = ai/3.0_wp * (uxxxx/4.0 &
                    + uxxx/dlx)
                if (kswy/=1 .and. (j==1 .or. j==l)) ty = dj/3.0_wp * (uyyyy/4.0 &
                    + uyyy/dly)
                grhs(i, j) = grhs(i, j) + dlx**2*tx + dly**2*ty
            end do
        end do
        !
        !     reset the right hand side in usol
        !
        usol(is:ms, js:ns) = grhs(is:ms, js:ns)

    end subroutine defer

end module module_sepeli
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    version 5.0, fortran 90 changes
!-----------------------------------------------------------------------
