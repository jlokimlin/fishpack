module module_sepx4

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_genbun, only: &
        genbun

    use module_sepaux, only: &
        SepAux, &
        get_coefficients

    ! Explicit typing only!
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: sepx4

contains


    subroutine sepx4(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, &
        d, n, nbdcnd, bdc, bdd, cofx, grhs, usol, idmn, pertrb, &
        ierror)
        !
        !     file sepx4.f
        !
        !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        !     *                                                               *
        !     *                  copyright(c) 2005 by UCAR                   *
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
        !     *                Boulder, Colorado (80307)  U.S.A.             *
        !     *                                                               *
        !     *                   which is sponsored by                       *
        !     *                                                               *
        !     *              the National Science Foundation                  *
        !     *                                                               *
        !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        !
        !     SUBROUTINE sepx4(IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA, C, D, N,
        !    +                  NBDCND, BDC, BDD, COFX, GRHS, USOL, IDMN, PERTRB,
        !    +                  ierror)
        !
        !
        !
        ! DIMENSION OF           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
        ! ARGUMENTS              USOL(IDMN, N+1),     GRHS(IDMN, N+1),
        !
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                sepx4 SOLVES FOR EITHER THE SECOND-ORDER
        !                        FINITE DIFFERENCE APPROXIMATION OR A
        !                        FOURTH-ORDER APPROXIMATION TO A SEPARABLE
        !                        ELLIPTIC EQUATION
        !
        !                          AF(X)*UXX+BF(X)*UX+CF(X)*U+UYY = G(X, Y)
        !
        !                        ON A RECTANGLE(X GREATER THAN OR EQUAL TO
        !                        A AND LESS THAN OR EQUAL TO B, Y GREATER THAN
        !                        OR EQUAL TO C AND LESS THAN OR EQUAL TO D).
        !                        ANY COMBINATION OF PERIODIC OR MIXED BOUNDARY
        !                        CONDITIONS IS ALLOWED.  IF BOUNDARY
        !                        CONDITIONS IN THE X DIRECTION ARE PERIODIC
        !                       (SEE MBDCND=0 BELOW) THEN THE COEFFICIENTS
        !                        MUST SATISFY
        !
        !                          AF(X)=C1, BF(X)=0, CF(X)=C2 FOR ALL X.
        !
        !                        HERE C1, C2 ARE CONSTANTS, C1.GT.0.
        !
        !                        THE POSSIBLE BOUNDARY CONDITIONS ARE:
        !                        IN THE X-DIRECTION:
        !                         (0) PERIODIC, U(X+B-A, Y)=U(X, Y) FOR
        !                              ALL Y, X
        !                         (1) U(A, Y), U(B, Y) ARE SPECIFIED FOR ALL Y
        !                         (2) U(A, Y), DU(B, Y)/DX+BETA*U(B, Y) ARE
        !                              SPECIFIED FOR ALL Y
        !                         (3) DU(A, Y)/DX+ALPHA*U(A, Y), DU(B, Y)/DX+
        !                              BETA*U(B, Y) ARE SPECIFIED FOR ALL Y
        !                         (4) DU(A, Y)/DX+ALPHA*U(A, Y), U(B, Y) ARE
        !                              SPECIFIED FOR ALL Y
        !
        !                        IN THE Y-DIRECTION:
        !                         (0) PERIODIC, U(X, Y+D-C)=U(X, Y) FOR ALL X, Y
        !                         (1) U(X, C), U(X, D) ARE SPECIFIED FOR ALL X
        !                         (2) U(X, C), DU(X, D)/DY ARE SPECIFIED FOR
        !                              ALL X
        !                         (3) DU(X, C)/DY, DU(X, D)/DY ARE SPECIFIED FOR
        !                              ALL X
        !                         (4) DU(X, C)/DY, U(X, D) ARE SPECIFIED FOR
        !                              ALL X
        !
        ! USAGE                  CALL sepx4(IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB,
        !                                   BETA, C, D, N, NBDCND, BDC, BDD, COFX,
        !                                   GRHS, USOL, IDMN, W, PERTRB, ierror)
        !
        ! ARGUMENTS
        ! ON INPUT               IORDER
        !                          = 2 IF A SECOND-ORDER APPROXIMATION IS
        !                              SOUGHT
        !                          = 4 IF A FOURTH-ORDER APPROXIMATION IS
        !                              SOUGHT
        !
        ! *** caution ***          GRHS SHOULD BE RESET IF sepx4 WAS FIRST CALLED
        !                          WITH IORDER=2 AND WILL BE CALLED AGAIN WITH
        !                          IORDER=4.  VALueS IN GRHS ARE DESTROYED BY THE
        !                          IORDER=2 CALL.
        !
        !
        !                        A, B
        !                          THE RANGE OF THE X-INDEPENDENT VARIABLE,
        !                          I.E., X IS GREATER THAN OR EQUAL TO A
        !                          AND LESS THAN OR EQUAL TO B.  A MUST BE
        !                          LESS THAN B.
        !
        !                        M
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL(A, B) IS SUBDIVIDED.  HENCE,
        !                          THERE WILL BE M+1 GRID POINTS IN THE X-
        !                          DIRECTION GIVEN BY XI=A+(I-1)*DLX
        !                          FOR I=1, 2, ..., M+1 WHERE DLX=(B-A)/M IS
        !                          THE PANEL WIDTH.  M MUST BE LESS THAN
        !                          IDMN AND GREATER THAN 5.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITION
        !                          AT X=A AND X=B
        !                          = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.,
        !                              U(X+B-A, Y)=U(X, Y)  FOR ALL Y, X
        !                          = 1 IF THE SOLUTION IS SPECIFIED AT X=A
        !                              AND X=B, I.E., U(A, Y) AND U(B, Y) ARE
        !                              SPECIFIED FOR ALL Y
        !                          = 2 IF THE SOLUTION IS SPECIFIED AT X=A
        !                              AND THE BOUNDARY CONDITION IS MIXED AT
        !                              X=B, I.E., U(A, Y) AND
        !                              DU(B, Y)/DX+BETA*U(B, Y) ARE SPECIFIED
        !                              FOR ALL Y
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
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALueS OF
        !                          DU(A, Y)/DX+ ALPHA*U(A, Y) AT X=A, WHEN
        !                          MBDCND=3 OR 4.
        !                          BDA(J) = DU(A, YJ)/DX+ALPHA*U(A, YJ),
        !                          J=1, 2, ..., N+1
        !                          WHEN MBDCND HAS ANY OTHER VALue, BDA IS
        !                          A DUMMY PARAMETER.
        !
        !                        ALPHA
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN CASE
        !                          OF A MIXED BOUNDARY CONDITION AT X=A
        !                         (SEE ARGUMENT BDA).  IF MBDCND IS NOT EQUAL
        !                          TO EITHER 3 OR 4, THEN ALPHA IS A DUMMY
        !                          PARAMETER.
        !
        !                        BDB
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALueS OF
        !                          DU(B, Y)/DX+ BETA*U(B, Y) AT X=B.
        !                          WHEN MBDCND=2 OR 3
        !                          BDB(J) = DU(B, YJ)/DX+BETA*U(B, YJ),
        !                          J=1, 2, ..., N+1
        !                          WHEN MBDCND HAS ANY OTHER VALue, BDB IS
        !                          A DUMMY PARAMETER.
        !
        !                        BETA
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN
        !                          CASE OF A MIXED BOUNDARY CONDITION AT X=B
        !                         (SEE ARGUMENT BDB).  IF MBDCND IS NOT EQUAL
        !                          TO 2 OR 3, THEN BETA IS A DUMMY PARAMETER.
        !
        !                        C, D
        !                          THE RANGE OF THE Y-INDEPENDENT VARIABLE,
        !                          I.E., Y IS GREATER THAN OR EQUAL TO C AND
        !                          LESS THAN OR EQUAL TO D.  C MUST BE LESS
        !                          THAN D.
        !
        !                        N
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL(C, D) IS SUBDIVIDED.  HENCE,
        !                          THERE WILL BE N+1 GRID POINTS IN THE Y-
        !                          DIRECTION GIVEN BY YJ=C+(J-1)*DLY FOR
        !                          J=1, 2, ..., N+1 WHERE DLY=(D-C)/N IS THE
        !                          PANEL WIDTH.  IN ADDITION, N MUST BE
        !                          GREATER THAN 4.
        !
        !                        NBDCND
        !                          INDICATES THE TYPES OF BOUNDARY CONDITIONS
        !                          AT Y=C AND Y=D
        !                          = 0 IF THE SOLUTION IS PERIODIC IN Y,
        !                              I.E., U(X, Y+D-C)=U(X, Y) FOR ALL X, Y
        !                          = 1 IF THE SOLUTION IS SPECIFIED AT Y=C
        !                              AND Y = D, I.E., U(X, C)  AND U(X, D)
        !                              ARE SPECIFIED FOR ALL X
        !                          = 2 IF THE SOLUTION IS SPECIFIED AT Y=C
        !                              AND THE BOUNDARY CONDITION IS MIXED
        !                              AT Y=D, I.E., DU(X, C)/DY AND U(X, D)
        !                              ARE SPECIFIED FOR ALL X
        !                          = 3 IF THE BOUNDARY CONDITIONS ARE MIXED
        !                              AT Y=CAND Y=D I.E.,
        !                              DU(X, D)/DY AND DU(X, D)/DY ARE
        !                              SPECIFIED FOR ALL X
        !                          = 4 IF THE BOUNDARY CONDITION IS MIXED
        !                              AT Y=C AND THE SOLUTION IS SPECIFIED
        !                              AT Y=D, I.E. DU(X, C)/DY+GAMA*U(X, C)
        !                              AND U(X, D) ARE SPECIFIED FOR ALL X
        !
        !                        BDC
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIES THE VALue DU(X, C)/DY AT Y=C.
        !
        !                          WHEN NBDCND=3 OR 4
        !                            BDC(I) = DU(XI, C)/DY I=1, 2, ..., M+1.
        !
        !                          WHEN NBDCND HAS ANY OTHER VALue, BDC IS
        !                          A DUMMY PARAMETER.
        !
        !                        BDD
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIED THE VALue OF DU(X, D)/DY AT Y=D.
        !
        !                          WHEN NBDCND=2 OR 3
        !                            BDD(I)=DU(XI, D)/DY I=1, 2, ..., M+1.
        !
        !                          WHEN NBDCND HAS ANY OTHER VALue, BDD IS
        !                          A DUMMY PARAMETER.
        !
        !                        COFX
        !                          A USER-SUPPLIED SUBPROGRAM WITH PARAMETERS
        !                          X, AFUN, BFUN, CFUN WHICH RETURNS THE
        !                          VALueS OF THE X-DEPENDENT COEFFICIENTS
        !                          AF(X), BF(X), CF(X) IN THE ELLIPTIC
        !                          EQUATION AT X.  IF BOUNDARY CONDITIONS IN
        !                          THE X DIRECTION ARE PERIODIC THEN THE
        !                          COEFFICIENTS MUST SATISFY AF(X)=C1, BF(X)=0,
        !                          CF(X)=C2 FOR ALL X.  HERE C1.GT.0
        !                          AND C2 ARE CONSTANTS.
        !
        !                          NOTE THAT COFX MUST BE DECLARED EXTERNAL
        !                          IN THE CALLING ROUTINE.
        !
        !                        GRHS
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALueS OF THE RIGHT-HAND SIDE OF THE
        !                          ELLIPTIC EQUATION, I.E., GRHS(I, J)=G(XI, YI),
        !                          FOR I=2, ..., M, J=2, ..., N.  AT THE
        !                          BOUNDARIES, GRHS IS DEFINED BY
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
        !                          WHERE * MEANS THESE QUANTITES ARE NOT USED.
        !                          GRHS SHOULD BE DIMENSIONED IDMN BY AT LEAST
        !                          N+1 IN THE CALLING ROUTINE.
        !
        ! *** caution              GRHS SHOULD BE RESET IF sepx4 WAS FIRST CALLED
        !                          WITH IORDER=2 AND WILL BE CALLED AGAIN WITH
        !                          IORDER=4.  VALueS IN GRHS ARE DESTROYED BY THE
        !                          IORDER=2 CALL.
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
        !                          WHERE * MEANS THE QUANTITES ARE NOT USED
        !                          IN THE SOLUTION.
        !
        !                          IF IORDER=2, THE USER MAY EQUIVALENCE GRHS
        !                          AND USOL TO SAVE SPACE.  NOTE THAT IN THIS
        !                          CASE THE TABLES SPECIFYING THE BOUNDARIES
        !                          OF THE GRHS AND USOL ARRAYS DETERMINE THE
        !                          BOUNDARIES UNIQueLY EXCEPT AT THE CORNERS.
        !                          IF THE TABLES CALL FOR BOTH G(X, Y) AND
        !                          U(X, Y) AT A CORNER THEN THE SOLUTION MUST
        !                          BE CHOSEN.
        !                          FOR EXAMPLE, IF MBDCND=2 AND NBDCND=4,
        !                          THEN U(A, C), U(A, D), U(B, D) MUST BE CHOSEN
        !                          AT THE CORNERS IN ADDITION TO G(B, C).
        !
        !                          IF IORDER=4, THEN THE TWO ARRAYS, USOL AND
        !                          GRHS, MUST BE DISTINCT.
        !
        !                          USOL SHOULD BE DIMENSIONED IDMN BY AT LEAST
        !                          N+1 IN THE CALLING ROUTINE.
        !
        !                        IDMN
        !                          THE ROW(OR FIRST) DIMENSION OF THE ARRAYS
        !                          GRHS AND USOL AS IT APPEARS IN THE PROGRAM
        !                          CALLING SEPELI.  THIS PARAMETER IS USED
        !                          TO SPECIFY THE VARIABLE DIMENSION OF GRHS
        !                          AND USOL.  IDMN MUST BE AT LEAST 7 AND
        !                          GREATER THAN OR EQUAL TO M+1.
        !
        !
        ! ON OUTPUT              USOL
        !                          CONTAINS THE APPROXIMATE SOLUTION TO THE
        !                          ELLIPTIC EQUATION. USOL(I, J) IS THE
        !                          APPROXIMATION TO U(XI, YJ) FOR I=1, 2..., M+1
        !                          AND J=1, 2, ..., N+1.  THE APPROXIMATION HAS
        !                          ERROR O(DLX**2+DLY**2) IF CALLED WITH
        !                          IORDER=2 AND O(DLX**4+DLY**4) IF CALLED
        !                          WITH IORDER=4.
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
        !                          BOUNDARY CONDITIONS(I.E., ALPHA=BETA=0 IF
        !                          MBDCND=3) IS SPECIFIED AND IF CF(X)=0 FOR
        !                          ALL X THEN A SOLUTION TO THE DISCRETIZED
        !                          MATRIX EQUATION MAY NOT EXIST
        !                         (REFLECTING THE NON-UNIQueNESS OF SOLUTIONS
        !                          TO THE PDE).
        !                          PERTRB IS A CONSTANT CALCULATED AND
        !                          SUBTRACTED FROM THE RIGHT HAND SIDE OF THE
        !                          MATRIX EQUATION INSURING THE EXISTENCE OF A
        !                          SOLUTION.  sepx4 COMPUTES THIS SOLUTION
        !                          WHICH IS A WEIGHTED MINIMAL LEAST SQUARES
        !                          SOLUTION TO THE ORIGINAL PROBLEM.  IF
        !                          SINGULARITY IS NOT DETECTED PERTRB=0.0 IS
        !                          RETURNED BY sepx4.
        !
        !                        ierror
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS OR FAILURE TO FIND A SOLUTION
        !
        !                          =  0 NO ERROR
        !                          =  1 IF A GREATER THAN B OR C GREATER
        !                               THAN D
        !                          =  2 IF MBDCND LESS THAN 0 OR MBDCND
        !                               GREATER THAN 4
        !                          =  3 IF NBDCND LESS THAN 0 OR NBDCND
        !                               GREATER THAN 4
        !                          =  4 IF ATTEMPT TO FIND A SOLUTION FAILS.
        !                              (THE LINEAR SYSTEM GENERATED IS NOT
        !                               DIAGONALLY DOMINANT.)
        !                          =  5 IF IDMN IS TOO SMALL(SEE DISCUSSION
        !                               OF IDMN)
        !                          =  6 IF M IS TOO SMALL OR TOO LARGE
        !                              (SEE DISCUSSION OF M)
        !                          =  7 IF N IS TOO SMALL(SEE DISCUSSION OF N)
        !                          =  8 IF IORDER IS NOT 2 OR 4
        !                          =  9 IF INTL IS NOT 0 OR 1
        !                          = 10 IF AFUN IS LESS THAN OR EQUAL TO ZERO
        !                               FOR SOME INTERIOR MESH POINT XI SOME
        !                               INTERIOR MESH POINT(XI, YJ)
        !                          = 12 IF MBDCND=0 AND AF(X)=CF(X)=CONSTANT
        !                               OR BF(X)=0 FOR ALL X IS NOT TRue.
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails(for example if N, M are too large
        !                               for your computer)
        !
        ! SPECIAL CONDITIONS     NONE
        !
        ! I/O                    NONE
        !
        ! REQUIRED files         fish.f, comf.f, genbun.f, gnbnaux.f, sepaux.f
        !
        !
        ! PRECISION              SINGLE
        !
        !
        ! LANGUAGE               FORTRAN 90
        !
        ! HISTORY                sepx4 WAS DEVELOPED AT NCAR BY JOHN C.
        !                        ADAMS OF THE SCIENTIFIC COMPUTING DIVISION
        !                        IN OCTOBER 1978.  THE BASIS OF THIS CODE IS
        !                        NCAR ROUTINE SEPELI.  BOTH PACKAGES WERE
        !                        RELEASED ON NCAR'S PUBLIC LIBRARIES IN
        !                        JANUARY 1980. sepx4 was modified in June 2004
        !                        incorporating fortran 90 dynamical storage
        !                        allocation for work space requirements
        !
        ! PORTABILITY            FORTRAN 90
        !
        ! ALGORITHM              sepx4 AUTOMATICALLY DISCRETIZES THE SEPARABLE
        !                        ELLIPTIC EQUATION WHICH IS THEN SOLVED BY A
        !                        GENERALIZED CYCLIC REDUCTION ALGORITHM IN THE
        !                        SUBROUTINE POIS.  THE FOURTH ORDER SOLUTION
        !                        IS OBTAINED USING THE TECHNIQue OF DEFFERRED
        !                        CORRECTIONS REFERENCED BELOW.
        !
        ! TIMING                 WHEN POSSIBLE, sepx4 SHOULD BE USED INSTEAD
        !                        OF PACKAGE SEPELI.  THE INCREASE IN SPEED
        !                        IS AT LEAST A FACTOR OF THREE.
        !
        ! REFERENCES             KELLER, H.B., NUMERICAL METHODS FOR TWO-POINT
        !                        BOUNDARY-VALue PROBLEMS, BLAISDEL(1968),
        !                        WALTHAM, MASS.
        !
        !                        SWARZTRAUBER, P., AND R. SWEET(1975):
        !                        EFFICIENT FORTRAN SUBPROGRAMS FOR THE
        !                        SOLUTION OF ELLIPTIC PARTIAL DIFFERENTIAL
        !                        EQUATIONS.  NCAR TECHNICAL NOTE
        !                          NCAR-TN/IA-109, PP. 135-137.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip),          intent (in)     :: iorder
        integer (ip),          intent (in)     :: m
        integer (ip),          intent (in)     :: mbdcnd
        integer (ip),          intent (in)     :: n
        integer (ip),          intent (in)     :: nbdcnd
        integer (ip),          intent (in)     :: idmn
        integer (ip),          intent (out)    :: ierror
        real (wp),             intent (in)     :: a
        real (wp),             intent (in)     :: b
        real (wp),             intent (in)     :: alpha
        real (wp),             intent (in)     :: beta
        real (wp),             intent (in)     :: c
        real (wp),             intent (in)     :: d
        real (wp),             intent (out)    :: pertrb
        real (wp), contiguous, intent (in)     :: bda(:)
        real (wp), contiguous, intent (in)     :: bdb(:)
        real (wp), contiguous, intent (in)     :: bdc(:)
        real (wp), contiguous, intent (in)     :: bdd(:)
        real (wp), contiguous, intent (in out) :: grhs(:,:)
        real (wp), contiguous, intent (out)    :: usol(:,:)
        procedure (get_coefficients)           :: cofx
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)             :: l, k, length, irwk, icwk
        integer (ip)             :: sepx4_workspace_indices(13)
        type (FishpackWorkspace) :: workspace
        type (SepAux)            :: sep_aux
        !-----------------------------------------------

        !
        !==> Check input parameters
        !
        call check_input_parameters(iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, cofx, idmn, ierror)

        if (ierror /= 0) then
            return
        end if

        !
        !==> compute minimum work space and check work space length input
        !
        select case (nbdcnd)
            case (0)
                l = n
                k = m + 1
            case default
                l = n + 1
                k = m + 1
        end select

        !
        !==> set real and complex workspace sizes
        !
        call compute_sepx4_workspace_dimensions(n, m, l, k, length, irwk, icwk)
        !
        !==> Allocate memory
        !
        call workspace%create(irwk, icwk, ierror)

        !
        !==> set sepx4 workspace indices
        !
        sepx4_workspace_indices=get_sepx4_workspace_indices(length, l, k)

        !
        !==> Solve system
        !
        associate( &
            i => sepx4_workspace_indices, &
            rew => workspace%real_workspace &
            )

            call s4elip(sep_aux, iorder, a, b, m, mbdcnd, &
                bda, alpha, bdb, beta, c, d, n,  nbdcnd, bdc, bdd, cofx, &
                rew(i(1)), rew(i(2)), rew(i(3)),  rew(i(4)), &
                rew(i(5)), rew(i(6)), rew(i(7):i(7)), rew(i(8):i(8)), &
                rew(i(9):i(9)), rew(i(10)), rew(i(11)), rew(i(12)), &
                grhs, usol, idmn, rew(i(13)), pertrb, ierror)

        end associate

        !
        !==> Release memory
        !
        call workspace%destroy()

    end subroutine sepx4



    pure subroutine compute_sepx4_workspace_dimensions(n, m, l, k, length, irwk, icwk)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)  :: n, m, l, k
        integer (ip), intent (out) :: length
        integer (ip), intent (out) :: irwk
        integer (ip), intent (out) :: icwk
        !--------------------------------------------------------------------------------
        integer (ip) :: log2n
        !--------------------------------------------------------------------------------

        log2n = int(log(real(n + 1, kind=wp))/log(2.0_wp) + 0.5_wp, kind=ip)
        length = 4*(n + 1) +(10 + log2n) * k

        ! set real and complex workspace sizes
        irwk = length + 6 * (k + l) + 1
        icwk = 0

    end subroutine compute_sepx4_workspace_dimensions


    pure function get_sepx4_workspace_indices(length, l, k) result (return_value)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent (in) :: length
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: k
        integer (ip)              :: return_value(13)
        !--------------------------------------------------------------------------------
        integer (ip) :: j !! Counter
        !--------------------------------------------------------------------------------

        associate( i => return_value)
            i(1) = length + 1

            do j = 1, 6
                i(j+1) = i(j) + l
            end do

            do j = 7, 11
                i(j+1) = i(j) + k
            end do

            i(13) = 1
        end associate

    end function get_sepx4_workspace_indices



    subroutine s4elip(sep_aux, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, &
        c, d, n, nbdcnd, bdc, bdd, cofx, an, bn, cn, dn, un, zn, am, bm, &
        cm, dm, um, zm, grhs, usol, idmn, w, pertrb, ierror)
        !
        ! Purpose:
        !
        !     s4elip sets up vectors and arrays for input to blktri
        !     and computes a second order solution in usol.  a return jump to
        !     sepeli occurrs if iorder=2.  if iorder=4 a fourth order
        !     solution is generated in usol.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SepAux),        intent (in out) :: sep_aux
        integer (ip),          intent (in)     :: iorder
        integer (ip),          intent (in)     :: m
        integer (ip),          intent (in)     :: mbdcnd
        integer (ip),          intent (in)     :: n
        integer (ip),          intent (in)     :: nbdcnd
        integer (ip),          intent (in)     :: idmn
        integer (ip),          intent (in out) :: ierror
        real (wp),             intent (in)     :: a
        real (wp),             intent (in)     :: b
        real (wp),             intent (in)     :: alpha
        real (wp),             intent (in)     :: beta
        real (wp),             intent (in)     :: c
        real (wp),             intent (in)     :: d
        real (wp),             intent (out)    :: pertrb
        real (wp), contiguous, intent (in)     :: bda(:)
        real (wp), contiguous, intent (in)     :: bdb(:)
        real (wp), contiguous, intent (in)     :: bdc(:)
        real (wp), contiguous, intent (in)     :: bdd(:)
        real (wp),             intent (in out) :: an(*)
        real (wp),             intent (in out) :: bn(*)
        real (wp),             intent (in out) :: cn(*)
        real (wp),             intent (in out) :: dn(*)
        real (wp),             intent (in out) :: un(*)
        real (wp),             intent (in out) :: zn(*)
        real (wp), contiguous, intent (in out) :: am(:)
        real (wp), contiguous, intent (in out) :: bm(:)
        real (wp), contiguous, intent (in out) :: cm(:)
        real (wp),             intent (in out) :: dm(*)
        real (wp),             intent (in out) :: um(*)
        real (wp),             intent (in out) :: zm(*)
        real (wp),             intent (in out) :: grhs(idmn,*)
        real (wp),             intent (in out) :: usol(idmn,*)
        real (wp),             intent (in out) :: w(*)
        procedure (get_coefficients)           :: cofx
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: i, j, i1, mp, np, ieror
        real (wp)    :: xi, ai, bi, ci, axi, bxi, cxi
        real (wp)    :: dyj, eyj, fyj, ax1, cxm
        real (wp)    :: dy1, fyn, gama, xnu, prtrb
        logical      :: singular
        !-----------------------------------------------

        ! Associate various quantities
        associate( &
            kswx => sep_aux%kswx, &
            kswy => sep_aux%kswy, &
            k => sep_aux%k, &
            l=>sep_aux%l, &
            mit=>sep_aux%mit, &
            nit=> sep_aux%nit, &
            is=> sep_aux%is, &
            ms=> sep_aux%ms, &
            js=> sep_aux%js, &
            ns=> sep_aux%ns, &
            ait => sep_aux%ait, &
            bit => sep_aux%bit, &
            cit => sep_aux%cit, &
            dit => sep_aux%dit, &
            dlx => sep_aux%dlx, &
            dly => sep_aux%dly, &
            tdlx3 => sep_aux%tdlx3, &
            tdly3 => sep_aux%tdly3, &
            dlx4 => sep_aux%dlx4, &
            dly4 => sep_aux%dly4 &
            )

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
            dly =(dit - cit)/n
            !
            !==> set right hand side values from grhs in usol on the interior
            !     and non-specified boundaries.
            !
            usol(2:m, 2:n) = (dly**2) * grhs(2:m, 2:n)

            if (kswx /= 2 .and. kswx /= 3) then
                usol(1, 2:n) = (dly**2) * grhs(1, 2:n)
            end if

            if (kswx /= 2 .and. kswx /= 5) then
                usol(k, 2:n) = (dly**2) * grhs(k, 2:n)
            end if

            if (kswy /= 2 .and. kswy /= 3) then
                usol(2:m, 1) = (dly**2) * grhs(2:m, 1)
            end if

            if (kswy /= 2 .and. kswy /= 5) then
                usol(2:m, l) = (dly**2) * grhs(2:m, l)
            end if

            if (kswx /= 2 .and. kswx /= 3 .and. kswy /= 2 .and. kswy /= 3) then
                usol(1, 1) = (dly**2) * grhs(1, 1)
            end if

            if (kswx /= 2 .and. kswx /= 5 .and. kswy /= 2 .and. kswy /= 3) then
                usol(k, 1) = (dly**2) * grhs(k, 1)
            end if

            if (kswx /= 2 .and. kswx /= 3 .and. kswy /= 2 .and. kswy /= 5) then
                usol(1, l) = (dly**2) * grhs(1, l)
            end if

            if (kswx /= 2 .and. kswx /= 5 .and. kswy /= 2 .and. kswy /= 5) then
                usol(k, l) = (dly**2) * grhs(k, l)
            end if

            i1 = 1
            !
            !==> set switches for periodic or non-periodic boundaries
            !
            if (kswx == 1) then
                mp = 0
            else
                mp = 1
            end if

            np = nbdcnd
            !
            !     set dlx, dly and size of block tri-diagonal system generated
            !     in nint, mint
            !
            dlx =(bit - ait)/m
            mit = k - 1

            if (kswx == 2) then
                mit = k - 2
            end if

            if (kswx == 4) then
                mit = k
            end if

            dly =(dit - cit)/n
            nit = l - 1

            if (kswy == 2) then
                nit = l - 2
            end if

            if (kswy == 4) then
                nit = l
            end if

            tdlx3 = 2.0_wp * (dlx**3)
            dlx4 = dlx**4
            tdly3 = 2.0_wp * (dly**3)
            dly4 = dly**4
            !
            !==> set subscript limits for portion of array to input to blktri
            !
            if (kswx==2 .or. kswx==3) then
                is = 2
            else
                is = 1
            end if

            if (kswy==2 .or. kswy==3) then
                js = 2
            else
                js = 1
            end if

            ns = nit + js - 1
            ms = mit + is - 1
            !
            !     set x - direction
            !
            do i = 1, mit
                xi = ait + real(is + i - 2)*dlx
                call cofx(xi, ai, bi, ci)
                axi =(ai/dlx - 0.5*bi)/dlx
                bxi =(-2.0 * ai/dlx**2) + ci
                cxi =(ai/dlx + 0.5*bi)/dlx
                am(i) = (dly**2) * axi
                bm(i) = dly**2*bxi
                cm(i) = dly**2*cxi
            end do
            !
            !     set y direction
            !
            dyj = 1.0_wp
            eyj = -2.0_wp
            fyj = 1.0_wp
            an(:nit) = dyj
            bn(:nit) = eyj
            cn(:nit) = fyj
            !
            !     adjust edges in x direction unless periodic
            !
            ax1 = am(1)
            cxm = cm(mit)
            select case(kswx)
                case(2)
                    !
                    !     dirichlet-dirichlet in x direction
                    !
                    am(1) = 0.0
                    cm(mit) = 0.0
                case(3)
                    !
                    !     dirichlet-mixed in x direction
                    !
                    am(1) = 0.0
                    am(mit) = am(mit) + cxm
                    bm(mit) = bm(mit) - 2.0 * beta*dlx*cxm
                    cm(mit) = 0.0
                case(4)
                    !
                    !     mixed - mixed in x direction
                    !
                    am(1) = 0.0_wp
                    bm(1) = bm(1) + 2.0_wp * dlx * alpha * ax1
                    cm(1) = cm(1) + ax1
                    am(mit) = am(mit) + cxm
                    bm(mit) = bm(mit) - 2.0_wp * dlx * beta * cxm
                    cm(mit) = 0.0_wp
                case(5)
                    !
                    !     mixed-dirichlet in x direction
                    !
                    am(1) = 0.0_wp
                    bm(1) = bm(1) + 2.0_wp * alpha * dlx * ax1
                    cm(1) = cm(1) + ax1
                    cm(mit) = 0.0_wp
            end select
            !
            !     adjust in y direction unless periodic
            !
            dy1 = an(1)
            fyn = cn(nit)
            gama = 0.0_wp
            xnu = 0.0_wp
            select case(kswy)
                case(2)
                    !
                    !     dirichlet-dirichlet in y direction
                    !
                    an(1) = 0.0_wp
                    cn(nit) = 0.0_wp
                case(3)
                    !
                    !     dirichlet-mixed in y direction
                    !
                    an(1) = 0.0_wp
                    an(nit) = an(nit) + fyn
                    bn(nit) = bn(nit) - 2.0_wp * dly * xnu * fyn
                    cn(nit) = 0.0_wp
                case(4)
                    !
                    !     mixed - mixed direction in y direction
                    !
                    an(1) = 0.0_wp
                    bn(1) = bn(1) + 2.0_wp * dly * gama * dy1
                    cn(1) = cn(1) + dy1
                    an(nit) = an(nit) + fyn
                    bn(nit) = bn(nit) - 2.0_wp * dly * xnu * fyn
                    cn(nit) = 0.0_wp
                case(5)
                    !
                    !     mixed-dirichlet in y direction
                    !
                    an(1) = 0.0_wp
                    bn(1) = bn(1) + 2.0_wp * dly * gama * dy1
                    cn(1) = cn(1) + dy1
                    cn(nit) = 0.0_wp
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
            !==> save adjusted edges in grhs if iorder=4
            !
            if (iorder == 4) then
                grhs(is, js:ns) = usol(is, js:ns)
                grhs(ms, js:ns) = usol(ms, js:ns)
                grhs(is:ms, js) = usol(is:ms, js)
                grhs(is:ms, ns) = usol(is:ms, ns)
            end if

            pertrb = 0.0_wp
            !
            !==> check if operator is singular
            !
            call check_if_singular(sep_aux, mbdcnd, nbdcnd, alpha, beta, cofx, singular)
            !
            !     compute non-zero eigenvector in null space of transpose
            !     if singular
            !
            if (singular .eqv. .true.) then
                call sep_aux%septri(mit, am, bm, cm, dm, um, zm)
            end if

            if (singular .eqv. .true.) then
                call sep_aux%septri(nit, an, bn, cn, dn, un, zn)
            end if
            !
            !     adjust right hand side if necessary
            !
            if (singular .eqv. .true.) then
                call sep_aux%seport(usol, idmn, zn, zm, pertrb)
            end if

            !
            !==> compute solution
            !
            !     save adjusted right hand side in grhs
            grhs(is:ms, js:ns) = usol(is:ms, js:ns)

            call genbun(np, nit, mp, mit, am, bm, cm, idmn, usol(is, js), ieror)
            !
            !==>  Check if error detected in pois
            !     this can only correspond to ierror=12
            if (ieror /= 0) then
                !       set error flag if improper coefficients input to pois
                ierror = 12
                return
            end if

            if (ierror /= 0) then
                return
            end if
            !
            !==> set periodic boundaries if necessary
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
            if (singular .eqv. .true.) then
                call sep_aux%sepmin(usol, idmn, zn, zm, prtrb)
            end if
            !
            !     return if deferred corrections and a fourth order solution are
            !     not flagged
            !
            if (iorder == 2) then
                return
            end if
            !
            !==> compute new right hand side for fourth order solution
            !
            call d4fer(sep_aux, cofx, idmn, usol, grhs)

            if (singular .eqv. .true.) then
                call sep_aux%seport(usol, idmn, zn, zm, pertrb)
            end if

            !
            !==> compute solution
            !
            !     save adjusted right hand side in grhs
            grhs(is:ms, js:ns) = usol(is:ms, js:ns)

            call genbun(np, nit, mp, mit, am, bm, cm, idmn, usol(is, js), ieror)

            !
            !==> check if error detected in pois
            !    this can only correspond to ierror=12
            !
            if (ieror /= 0) then
                !       set error flag if improper coefficients input to pois
                ierror = 12
                return
            end if

            if (ierror /= 0) then
                return
            end if
            !
            !==> set periodic boundaries if necessary
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
            if (singular .eqv. .true.) then
                call sep_aux%sepmin(usol, idmn, zn, zm, prtrb)
            end if

        end associate

    end subroutine s4elip


    subroutine check_input_parameters(iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, cofx, &
        idmn, ierror)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
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
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: i
        real (wp)    :: xi, ai, bi, ci
        real (wp)    :: dlx
        !-----------------------------------------------


        !
        !     this program checks the input parameters for errors
        !
        !
        !
        !     check definition of solution region
        !
        if (a >= b .or. c >= d) then
            ierror = 1
            return
        end if
        !
        !     check boundary switches
        !
        if (mbdcnd < 0 .or. mbdcnd > 4) then
            ierror = 2
            return
        end if

        if (nbdcnd < 0 .or. nbdcnd > 4) then
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
        !     check m
        !
        if (m > idmn - 1 .or. m < 6) then
            ierror = 6
            return
        end if
        !
        !     check n
        !
        if (n < 5) then
            ierror = 7
            return
        end if
        !
        !     check iorder
        !
        if (iorder /= 2 .and. iorder /= 4) then
            ierror = 8
            return
        end if
        !
        !     check that equation is elliptic
        !
        dlx =(b - a)/m
        do i = 2, m
            xi = a + real(i - 1, kind=wp ) * dlx
            call cofx(xi, ai, bi, ci)
            if (ai > 0.0_wp) then
                cycle
            end if
            ierror = 10
            return
        end do
        !
        !==> no error found
        !
        ierror = 0


    end subroutine check_input_parameters


    subroutine check_if_singular(sep_aux, mbdcnd, nbdcnd, alpha, beta, cofx, singlr)
        !
        ! Purpose:
        !
        !     this subroutine checks if the pde sepeli
        !     must solve is a singular operator
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SepAux), intent (in out) :: sep_aux
        integer (ip), intent (in)    :: mbdcnd
        integer (ip), intent (in)    :: nbdcnd
        real (wp),    intent (in)    :: alpha
        real (wp),    intent (in)    :: beta
        logical ,     intent (out)   :: singlr
        procedure (get_coefficients) :: cofx
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: i
        real (wp)    :: xi, ai, bi, ci
        !-----------------------------------------------

        ! Associate various quantities
        associate( &
            kswx => sep_aux%kswx, &
            kswy => sep_aux%kswy, &
            k => sep_aux%k, &
            l=>sep_aux%l, &
            mit=>sep_aux%mit, &
            nit=> sep_aux%nit, &
            is=> sep_aux%is, &
            ms=> sep_aux%ms, &
            js=> sep_aux%js, &
            ns=> sep_aux%ns, &
            ait => sep_aux%ait, &
            bit => sep_aux%bit, &
            cit => sep_aux%cit, &
            dit => sep_aux%dit, &
            dlx => sep_aux%dlx, &
            dly => sep_aux%dly, &
            tdlx3 => sep_aux%tdlx3, &
            tdly3 => sep_aux%tdly3, &
            dlx4 => sep_aux%dlx4, &
            dly4 => sep_aux%dly4 &
            )

            singlr = .false.
            !
            !     check if the boundary conditions are
            !     entirely periodic and/or mixed
            !
            if (mbdcnd /=0 .and. mbdcnd /=3 .or. nbdcnd /=0 .and. nbdcnd /= 3) then
                return
            end if
            !
            !     check that mixed conditions are pure neuman
            !
            if (mbdcnd == 3) then
                if (alpha /= 0.0_wp .or. beta /= 0.0_wp) then
                    return
                end if
            end if
            !
            !     check that non-derivative coefficient functions
            !     are zero
            !
            do i = is, ms
                xi = ait + real(i - 1, kind=wp)*dlx
                call cofx(xi, ai, bi, ci)
                if (ci == 0.0_wp) then
                    cycle
                end if
                return
            end do
            !
            !     the operator must be singular if this point is reached
            !
            singlr = .true.

        end associate

    end subroutine check_if_singular



    subroutine d4fer(sep_aux, cofx, idmn, usol, grhs)
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
        !    (via second order finite differencing) the truncation error
        !     and the result is added to the right hand side in grhs
        !     and then transferred to usol to be used as a new right
        !     hand side when calling blktri for a fourth order solution.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SepAux), intent (in out) :: sep_aux
        integer (ip), intent (in)     :: idmn
        real (wp),    intent (in out) :: usol(idmn, *)
        real (wp),    intent (in out) :: grhs(idmn, *)
        procedure (get_coefficients)  :: cofx
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: i, j
        real (wp)    :: xi, ai, bi, ci, uxxx, uxxxx, uyyy, uyyyy, tx, ty
        !-----------------------------------------------

        ! Associate various quantities
        associate( &
            kswx => sep_aux%kswx, &
            kswy => sep_aux%kswy, &
            k => sep_aux%k, &
            l=>sep_aux%l, &
            mit=>sep_aux%mit, &
            nit=> sep_aux%nit, &
            is=> sep_aux%is, &
            ms=> sep_aux%ms, &
            js=> sep_aux%js, &
            ns=> sep_aux%ns, &
            ait => sep_aux%ait, &
            bit => sep_aux%bit, &
            cit => sep_aux%cit, &
            dit => sep_aux%dit, &
            dlx => sep_aux%dlx, &
            dly => sep_aux%dly, &
            tdlx3 => sep_aux%tdlx3, &
            tdly3 => sep_aux%tdly3, &
            dlx4 => sep_aux%dlx4, &
            dly4 => sep_aux%dly4 &
            )

            !     compute truncation error approximation over the entire mesh
            !
            do i = is, ms
                xi = ait + real(i - 1)*dlx
                call cofx(xi, ai, bi, ci)
                do j = js, ns
                    !
                    !     compute partial derivative approximations at(xi, yj)
                    !
                    call sep_aux%sepdx(usol, idmn, i, j, uxxx, uxxxx)
                    call sep_aux%sepdy(usol, idmn, i, j, uyyy, uyyyy)
                    tx = ai*(uxxxx/12) + bi*(uxxx/6)
                    ty = uyyyy/12
                    !
                    !     reset form of truncation if at boundary which is non-periodic
                    !
                    if (kswx /= 1 .and. (i==1 .or. i==k)) then
                        tx = (ai/3) * ((uxxxx/4) + uxxx/dlx)
                    end if

                    if (kswy /= 1.and.(j==1.or.j==l)) then
                        ty=((uyyyy/4)+uyyy/dly)/3
                    end if

                    grhs(i, j) = grhs(i, j) + dly**2*(dlx**2*tx + dly**2*ty)
                end do
            end do
            !
            !     reset the right hand side in usol
            !
            usol(is:ms, js:ns) = grhs(is:ms, js:ns)

        end associate

    end subroutine d4fer


end module module_sepx4
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
