module module_hstcsp

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_blktri, only:&
        blktrii

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hstcsp


contains


    subroutine hstcsp(intl, a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, &
        bdc, bdd, elmbda, f, idimf, pertrb, ierror, w)

        !     file hstcsp.f
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
        !     SUBROUTINE hstcsp (INTL, A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC,
        !    +                   BDD, ELMBDA, F, IDIMF, PERTRB, ierror, W)
        !
        !
        ! DIMENSION OF           BDA(N), BDB(N), BDC(M), BDD(M), F(IDIMF, N)
        ! ARGUMENTS
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
        !                        DIFFERENCE APPROXIMATION ON A STAGGERED
        !                        GRID TO THE MODIFIED HELMHOLTZ EQUATION IN
        !                        SPHERICAL COORDINATES ASSUMING AXISYMMETRY
        !                        (NO DEPENDENCE ON LONGITUDE).
        !
        !                        THE EQUATION IS
        !
        !                           (1/R**2)(D/DR)(R**2(DU/DR)) +
        !                           1/(R**2*SIN(THETA))(D/DTHETA)
        !                           (SIN(THETA)(DU/DTHETA)) +
        !                           (LAMBDA/(R*SIN(THETA))**2)U  =  F(THETA, R)
        !
        !                        WHERE THETA IS COLATITUDE AND R IS THE
        !                        RADIAL COORDINATE. THIS TWO-DIMENSIONAL
        !                        MODIFIED HELMHOLTZ EQUATION RESULTS FROM
        !                        THE FOURIER TRANSFORM OF THE THREE-
        !                        DIMENSIONAL POISSON EQUATION.
        !
        !
        ! USAGE                  CALL hstcsp (INTL, A, B, M, MBDCND, BDA, BDB, C, D, N,
        !                                     NBDCND, BDC, BDD, ELMBDA, F, IDIMF,
        !                                     PERTRB, ierror, W)
        !
        ! ARGUMENTS
        !  ON INPUT              INTL
        !
        !                          = 0  ON INITIAL ENTRY TO hstcsp OR IF ANY
        !                               OF THE ARGUMENTS C, D, N, OR NBDCND
        !                               ARE CHANGED FROM A PREVIOUS CALL
        !
        !                          = 1  IF C, D, N, AND NBDCND ARE ALL
        !                               UNCHANGED FROM PREVIOUS CALL TO hstcsp
        !
        !                          NOTE:
        !                          A CALL WITH INTL = 0 TAKES APPROXIMATELY
        !                          1.5 TIMES AS MUCH TIME AS A CALL WITH
        !                          INTL = 1.  ONCE A CALL WITH INTL = 0
        !                          HAS BEEN MADE THEN SUBSEQUENT SOLUTIONS
        !                          CORRESPONDING TO DIFFERENT F, BDA, BDB,
        !                          BDC, AND BDD CAN BE OBTAINED FASTER WITH
        !                          INTL = 1 SINCE INITIALIZATION IS NOT
        !                          REPEATED.
        !
        !                        A, B
        !                          THE RANGE OF THETA (COLATITUDE),
        !                          I.E. A .LE. THETA .LE. B.  A
        !                          MUST BE LESS THAN B AND A MUST BE
        !                          NON-NEGATIVE.  A AND B ARE IN RADIANS.
        !                          A = 0 CORRESPONDS TO THE NORTH POLE AND
        !                          B = PI CORRESPONDS TO THE SOUTH POLE.
        !
        !                          * * *  IMPORTANT  * * *
        !
        !                          IF B IS EQUAL TO PI, THEN B MUST BE
        !                          COMPUTED USING THE STATEMENT
        !                              B = PI_MACH(DUM)
        !                          THIS INSURES THAT B IN THE USER'S PROGRAM
        !                          IS EQUAL TO PI IN THIS PROGRAM, PERMITTING
        !                          SEVERAL TESTS OF THE INPUT PARAMETERS THAT
        !                          OTHERWISE WOULD NOT BE POSSIBLE.
        !
        !                          * * * * * * * * * * * *
        !
        !                        M
        !                          THE NUMBER OF GRID POINTS IN THE INTERVAL
        !                          (A, B).  THE GRID POINTS IN THE THETA-
        !                          DIRECTION ARE GIVEN BY
        !                            THETA(I) = A + (I-0.5)DTHETA
        !                          FOR I=1, 2, ..., M WHERE DTHETA =(B-A)/M.
        !                          M MUST BE GREATER THAN 4.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT THETA = A AND THETA = B.
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = A AND THETA = B.
        !                               (SEE NOTES 1, 2 BELOW)
        !
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = A AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO THETA IS
        !                               SPECIFIED AT THETA = B
        !                               (SEE NOTES 1, 2 BELOW).
        !
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               AT THETA = A (SEE NOTES 1, 2 BELOW)
        !                               AND THETA = B.
        !
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED AT
        !                               THETA = A (SEE NOTES 1, 2 BELOW) AND
        !                               THE SOLUTION IS SPECIFIED AT THETA = B.
        !
        !                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = A = 0 AND THE SOLUTION IS
        !                               SPECIFIED AT THETA = B.
        !                               (SEE NOTE 2 BELOW)
        !
        !                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = A = 0 AND THE DERIVATIVE OF
        !                               THE SOLUTION WITH RESPECT TO THETA IS
        !                               SPECIFIED AT THETA = B
        !                               (SEE NOTE 2 BELOW).
        !
        !                          = 7  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = A AND THE SOLUTION IS
        !                               UNSPECIFIED AT THETA = B = PI.
        !
        !                          = 8  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED AT
        !                               THETA = A (SEE NOTE 1 BELOW)
        !                               AND THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = B = PI.
        !
        !                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
        !                                THETA = A = 0 AND THETA = B = PI.
        !
        !                          NOTE 1:
        !                          IF A = 0, DO NOT USE MBDCND = 1, 2, 3, 4, 7
        !                          OR 8, BUT INSTEAD USE MBDCND = 5, 6, OR 9.
        !
        !                          NOTE 2:
        !                          IF B = PI, DO NOT USE MBDCND = 1, 2, 3, 4, 5,
        !                          OR 6, BUT INSTEAD USE MBDCND = 7, 8, OR 9.
        !
        !                          NOTE 3:
        !                          WHEN A = 0  AND/OR B = PI THE ONLY
        !                          MEANINGFUL BOUNDARY CONDITION IS
        !                          DU/DTHETA = 0.   SEE D. GREENSPAN,
        !                          'NUMERICAL ANALYSIS OF ELLIPTIC
        !                           BOUNDARY VALUE PROBLEMS, '
        !                          HARPER AND ROW, 1965, CHAPTER 5.)
        !
        !                        BDA
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
        !                          SPECIFIES THE BOUNDARY VALUES (IF ANY) OF
        !                          THE SOLUTION AT THETA = A.
        !
        !                          WHEN  MBDCND = 1, 2, OR 7,
        !                            BDA(J) = U(A, R(J)),   J=1, 2, ..., N.
        !
        !                          WHEN MBDCND = 3, 4, OR 8,
        !                            BDA(J) = (D/DTHETA)U(A, R(J)), J=1, 2, ..., N.
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS A
        !                          DUMMY VARIABLE.
        !
        !                        BDB
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
        !                          SPECIFIES THE BOUNDARY VALUES OF THE
        !                          SOLUTION AT THETA = B.
        !
        !                          WHEN MBDCND = 1, 4, OR 5,
        !                            BDB(J) = U(B, R(J)),     J=1, 2, ..., N.
        !
        !                          WHEN MBDCND = 2, 3, OR 6,
        !                            BDB(J) = (D/DTHETA)U(B, R(J)), J=1, 2, ..., N.
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
        !                          A DUMMY VARIABLE.
        !
        !                        C, D
        !                          THE RANGE OF R , I.E. C .LE. R .LE. D.
        !                          C MUST BE LESS THAN D AND NON-NEGATIVE.
        !
        !                        N
        !                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
        !                          (C, D).  THE UNKNOWNS IN THE R-DIRECTION
        !                          ARE GIVEN BY R(J) = C + (J-0.5)DR,
        !                          J=1, 2, ..., N, WHERE DR = (D-C)/N.
        !                          N MUST BE GREATER THAN 4.
        !
        !                        NBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT R = C AND R = D.
        !
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               R = C AND R = D.
        !
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               R = C AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO R IS
        !                               SPECIFIED AT R = D. (SEE NOTE 1 BELOW)
        !
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO R IS SPECIFIED AT
        !                               R = C AND R = D.
        !
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO R IS
        !                               SPECIFIED AT R = C AND THE SOLUTION
        !                               IS SPECIFIED AT R = D.
        !
        !                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
        !                               R = C = 0 (SEE NOTE 2 BELOW) AND THE
        !                               SOLUTION IS SPECIFIED AT R = D.
        !
        !                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
        !                               R = C = 0 (SEE NOTE 2 BELOW)
        !                               AND THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO R IS SPECIFIED AT
        !                               R = D.
        !
        !                          NOTE 1:
        !                          IF C = 0 AND MBDCND = 3, 6, 8 OR 9, THE
        !                          SYSTEM OF EQUATIONS TO BE SOLVED IS
        !                          SINGULAR.  THE UNIQUE SOLUTION IS
        !                          DETERMINED BY EXTRAPOLATION TO THE
        !                          SPECIFICATION OF U(THETA(1), C).
        !                          BUT IN THESE CASES THE RIGHT SIDE OF THE
        !                          SYSTEM WILL BE PERTURBED BY THE CONSTANT
        !                          PERTRB.
        !
        !                          NOTE 2:
        !                          NBDCND = 5 OR 6 CANNOT BE USED WITH
        !                          MBDCND =1, 2, 4, 5, OR 7
        !                          (THE FORMER INDICATES THAT THE SOLUTION IS
        !                          UNSPECIFIED AT R = 0; THE LATTER INDICATES
        !                          SOLUTION IS SPECIFIED).
        !                          USE INSTEAD NBDCND = 1 OR 2.
        !
        !                        BDC
        !                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
        !                          SPECIFIES THE BOUNDARY VALUES OF THE
        !                          SOLUTION AT R = C.  WHEN NBDCND = 1 OR 2,
        !                            BDC(I) = U(THETA(I), C),    I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 3 OR 4,
        !                            BDC(I) = (D/DR)U(THETA(I), C), I=1, 2, ..., M.
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
        !                          A DUMMY VARIABLE.
        !
        !                        BDD
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
        !                          SPECIFIES THE BOUNDARY VALUES OF THE
        !                          SOLUTION AT R = D.  WHEN NBDCND = 1 OR 4,
        !                            BDD(I) = U(THETA(I), D) ,    I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 2 OR 3,
        !                            BDD(I) = (D/DR)U(THETA(I), D), I=1, 2, ..., M.
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
        !                          A DUMMY VARIABLE.
        !
        !                        ELMBDA
        !                          THE CONSTANT LAMBDA IN THE MODIFIED
        !                          HELMHOLTZ EQUATION.  IF LAMBDA IS GREATER
        !                          THAN 0, A SOLUTION MAY NOT EXIST.
        !                          HOWEVER, hstcsp WILL ATTEMPT TO FIND A
        !                          SOLUTION.
        !
        !                        F
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUES OF THE RIGHT SIDE OF THE MODIFIED
        !                          HELMHOLTZ EQUATION.  FOR I=1, 2, ..., M AND
        !                          J=1, 2, ..., N
        !
        !                                F(I, J) = F(THETA(I), R(J)) .
        !
        !                          F MUST BE DIMENSIONED AT LEAST M X N.
        !
        !                        IDIMF
        !                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
        !                          F AS IT APPEARS IN THE PROGRAM CALLING
        !                          hstcsp.  THIS PARAMETER IS USED TO SPECIFY
        !                          THE VARIABLE DIMENSION OF F.
        !                          IDIMF MUST BE AT LEAST M.
        !
        !                        W
        !                          A fortran 90 derived TYPE (FishpackWorkspace) variable
        !                          that must be declared by the user.  The first
        !                          two declarative statements in the user program
        !                          calling hstcsp must be:
        !
        !                               use type_FishpackWorkspace
        !                               TYPE (FishpackWorkspace) :: W
        !
        !                          The first statement makes the fishpack module
        !                          defined in the file "fish.f" available to the
        !                          user program calling hstcsp.  The second statement
        !                          declares a derived type variable (defined in
        !                          the module "fish.f") which is used internally
        !                          in blktri to dynamically allocate real and complex
        !                          work space used in solution.  An error flag
        !                          (ierror = 20) is set if the required work space
        !                          allocation fails (for example if N, M are too large)
        !                          Real and complex values are set in the components
        !                          of W on a initial (IFLG=0) call to hstcsp.  These
        !                          must be preserved on non-initial calls (INTL=1)
        !                          to hstcsp.  This eliminates redundant calculations
        !                          and saves compute time.
        !               ****       IMPORTANT!  The user program calling hstcsp should
        !                          include the statement:
        !
        !                               CALL FISHFIN(W)
        !
        !                          after the final approximation is generated by
        !                          hstcsp.  The will deallocate the real and complex
        !                          work space of W.  Failure to include this statement
        !                          could result in serious memory leakage.
        !
        !
        !
        ! ON OUTPUT              F
        !                          CONTAINS THE SOLUTION U(I, J) OF THE FINITE
        !                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
        !                          (THETA(I), R(J)) FOR I=1, 2, .., M, J=1, 2, ..., N.
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
        !                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
        !                          SPECIFIED FOR A POISSON EQUATION
        !                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
        !                          PERTRB IS A CONSTANT, CALCULATED AND
        !                          SUBTRACTED FROM F, WHICH ENSURES THAT A
        !                          SOLUTION EXISTS.  hstcsp THEN COMPUTES THIS
        !                          SOLUTION, WHICH IS A LEAST SQUARES SOLUTION
        !                          TO THE ORIGINAL APPROXIMATION.
        !                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
        !                          A SOLUTION; HENCE, THE SOLUTION IS NOT
        !                          UNIQUE.  THE VALUE OF PERTRB SHOULD BE
        !                          SMALL COMPARED TO THE RIGHT SIDE F.
        !                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
        !                          ESSENTIALLY DIFFERENT PROBLEM.
        !                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
        !                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
        !                          OBTAINED.
        !
        !                        ierror
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS. EXCEPT FOR NUMBERS 0 AND 10,
        !                          A SOLUTION IS NOT ATTEMPTED.
        !
        !                          =  0  NO ERROR
        !
        !                          =  1  A .LT. 0 OR B .GT. PI
        !
        !                          =  2  A .GE. B
        !
        !                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 9
        !
        !                          =  4  C .LT. 0
        !
        !                          =  5  C .GE. D
        !
        !                          =  6  NBDCND .LT. 1 OR NBDCND .GT. 6
        !
        !                          =  7  N .LT. 5
        !
        !                          =  8  NBDCND = 5 OR 6 AND
        !                                MBDCND = 1, 2, 4, 5, OR 7
        !
        !                          =  9  C .GT. 0 AND NBDCND .GE. 5
        !
        !                          = 10  ELMBDA .GT. 0
        !
        !                          = 11  IDIMF .LT. M
        !
        !                          = 12  M .LT. 5
        !
        !                          = 13  A = 0 AND MBDCND =1, 2, 3, 4, 7 OR 8
        !
        !                          = 14  B = PI AND MBDCND .LE. 6
        !
        !                          = 15  A .GT. 0 AND MBDCND = 5, 6, OR 9
        !
        !                          = 16  B .LT. PI AND MBDCND .GE. 7
        !
        !                          = 17  LAMBDA .NE. 0 AND NBDCND .GE. 5
        !
        !                          SINCE THIS IS THE ONLY MEANS OF INDICATING
        !                          A POSSIBLY INCORRECT CALL TO hstcsp,
        !                          THE USER SHOULD TEST ierror AFTER THE CALL.
        !
        !                        = 20 If the dynamic allocation of real and
        !                             complex work space in the derived type
        !                             (FishpackWorkspace) variable W fails (e.g.,
        !                             if N, M are too large for the platform used)
        !
        !                        W
        !                             The derived type (FishpackWorkspace) variable W
        !                             contains real and complex values that must not
        !                             be destroyed if hstcsp is called again with
        !                             IFLG=1.
        !
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED LIBRARY       fish.f, blktri.f, comf.f
        ! FILES
        !
        ! LANGUAGE               FORTRAN 90
        !
        ! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
        !                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
        !                        IN January 1980. Revised by John Adams in June
        !                        2004 using Fortan 90 dynamically allocated work
        !                        space and derived data types to eliminate mixed
        !                        mode conflicts in the earlier versions.
        !
        ! PORTABILITY            FORTRAN 90
        !
        ! ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
        !                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
        !                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR
        !                        AND CALLS blktri WHICH SOLVES THE LINEAR
        !                        SYSTEM OF EQUATIONS.
        !
        !
        ! TIMING                 FOR LARGE M AND N, THE OPERATION COUNT IS
        !                        ROUGHLY PROPORTIONAL TO M*N*LOG2(N).  THE
        !                        TIMING ALSO DEPENDS ON INPUT PARAMETER INTL.
        !
        ! ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
        !                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
        !                        DIGITS FOR N AND M AS LARGE AS 64.
        !                        MORE DETAILED INFORMATION ABOUT ACCURACY
        !                        CAN BE FOUND IN THE DOCUMENTATION FOR
        !                        SUBROUTINE blktri WHICH IS THE ROUTINE
        !                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
        !
        ! REFERENCES             P.N. SWARZTRAUBER, "A DIRECT METHOD FOR
        !                        THE DISCRETE SOLUTION OF SEPARABLE ELLIPTIC
        !                        EQUATIONS", 
        !                        SIAM J. NUMER. ANAL. 11(1974), PP. 1136-1150.
        !
        !                        U. SCHUMANN AND R. SWEET, "A DIRECT METHOD FOR
        !                        THE SOLUTION OF POISSON'S EQUATION WITH NEUMANN
        !                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
        !                        ARBITRARY SIZE, " J. COMP. PHYS. 20(1976), 
        !                        PP. 171-182.
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip) :: intl
        integer (ip) :: m
        integer (ip) :: mbdcnd
        integer (ip) :: n
        integer (ip) :: nbdcnd
        integer (ip) :: idimf
        integer (ip) :: ierror
        real (wp) :: a
        real (wp) :: b
        real (wp) :: c
        real (wp) :: d
        real (wp) :: elmbda
        real (wp) :: pertrb
        real (wp) :: bda(*)
        real (wp) :: bdb(*)
        real (wp) :: bdc(*)
        real (wp) :: bdd(*)
        real (wp), intent (in out) :: f(idimf,*)
        class (FishpackWorkspace) :: w
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip)             :: k, l, np, irwk, icwk, ierr1
        real (wp)                :: pi
        integer (ip), save       :: iw1, iwbm, iwcm, iwan, iwbn, iwcn, iwsnth, iwrsq
        !-----------------------------------------------

        pi = acos( -1.0 )
        !
        !     CHECK FOR INVALID INPUT PARAMETERS
        !
        ierror = 0
        if (a<0. .or. b>pi) ierror = 1
        if (a >= b) ierror = 2
        if (mbdcnd<1 .or. mbdcnd>9) ierror = 3
        if (c < 0.) ierror = 4
        if (c >= d) ierror = 5
        if (nbdcnd<1 .or. nbdcnd>6) ierror = 6
        if (n < 5) ierror = 7
        if ((nbdcnd==5 .or. nbdcnd==6) .and. (mbdcnd==1 .or. mbdcnd==2 &
            .or. mbdcnd==4 .or. mbdcnd==5 .or. mbdcnd==7)) ierror = 8
        if (c>0. .and. nbdcnd>=5) ierror = 9
        if (idimf < m) ierror = 11
        if (m < 5) ierror = 12
        if(a==0..and.mbdcnd/=5.and.mbdcnd/=6.and.mbdcnd/=9)ierror=13
        if (b==pi .and. mbdcnd<=6) ierror = 14
        if(a>0..and.(mbdcnd==5.or.mbdcnd==6.or.mbdcnd==9))ierror=15
        if (b<pi .and. mbdcnd>=7) ierror = 16
        if (elmbda/=0. .and. nbdcnd>=5) ierror = 17
        if (ierror == 0) then
            if (intl == 0) then
                !     allocate required work space
                k = m + 1
                l = n + 1
                np = nbdcnd
                !          compute blktri requirements in irwk, icwk
                call w%get_block_tridiagonal_workpace_dimensions (n, m, irwk, icwk)
                !     set work space indices
                iw1 = irwk + 1
                iwbm = iw1 + m
                iwcm = iwbm + m
                iwan = iwcm + m
                iwbn = iwan + n
                iwcn = iwbn + n
                iwsnth = iwcn + n
                iwrsq = iwsnth + m
                !     allocate hstcsp required work spac
                irwk = iwrsq + n
                icwk = icwk + 3*k
                call w%create( irwk, icwk, ierror )
                if (ierror == 20) return
            end if
            ierr1 = 0

            associate( &
                rew => w%real_workspace, &
                cxw => w%complex_workspace &
                )
                call hstcs1(intl, a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                    elmbda, f, idimf, pertrb, ierr1, rew(iw1), rew(iwbm), rew(iwcm), &
                    rew(iwan), rew(iwbn), rew(iwcn), rew(iwsnth), rew(iwrsq), &
                    rew, cxw)
            end associate
            ierror = ierr1
        end if

    end subroutine hstcsp


    subroutine hstcs1( intl, a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, &
        bdc, bdd, elmbda, f, idimf, pertrb, ierr1, am, bm, cm, an, bn, &
        cn, snth, rsq, w, wc)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer , intent (in) :: intl
        integer (ip) :: m
        integer , intent (in) :: mbdcnd
        integer (ip) :: n
        integer , intent (in) :: nbdcnd
        integer (ip) :: idimf
        integer (ip) :: ierr1
        real (wp), intent (in) :: a
        real (wp), intent (in) :: b
        real (wp), intent (in) :: c
        real (wp), intent (in) :: d
        real (wp), intent (in) :: elmbda
        real (wp), intent (out) :: pertrb
        real (wp), intent (in) :: bda(*)
        real (wp), intent (in) :: bdb(*)
        real (wp), intent (in) :: bdc(*)
        real (wp), intent (in) :: bdd(*)
        real (wp), intent (in out) :: f(idimf,*)
        real (wp) :: am(*)
        real (wp) :: bm(*)
        real (wp) :: cm(*)
        real (wp) :: an(*)
        real (wp) :: bn(*)
        real (wp) :: cn(*)
        real (wp), intent (in out) :: snth(*)
        real (wp), intent (in out) :: rsq(*)
        real (wp), intent (in out) :: w(*)
        complex  :: wc(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: i, j, isw, nb
        real (wp)    :: dth, dthsq, dr, x, y, a2, a1, a3
        !-----------------------------------------------

        dth = (b - a)/m
        dthsq = dth**2

        do i = 1, m
            snth(i) = sin(a + (real(i) - 0.5_wp)*dth)
        end do

        dr = (d - c)/n

        do j = 1, n
            rsq(j) = (c + (real(j) - 0.5_wp)*dr)**2
        end do
        !
        !     multiply right side by r(j)**2
        !
        do j = 1, n
            x = rsq(j)
            f(:m, j) = x*f(:m, j)
        end do
        !
        !      define coefficients am, bm, cm
        !
        x = 1./(2.0_wp*cos(dth/2.))
        am(2:m) = (snth(:m-1)+snth(2:m))*x
        cm(:m-1) = am(2:m)
        am(1) = sin(a)
        cm(m) = sin(b)
        do i = 1, m
            x = 1./snth(i)
            y = x/dthsq
            am(i) = am(i)*y
            cm(i) = cm(i)*y
            bm(i) = elmbda*x*x - am(i) - cm(i)
        end do
        !
        !==> Define coefficients an, bn, cn
        !
        x = c/dr
        do j = 1, n
            an(j) = (x + real(j - 1, kind=wp))**2
            cn(j) = (x + real(j, kind=wp))**2
            bn(j) = -(an(j)+cn(j))
        end do
        isw = 1
        nb = nbdcnd

        if (c == 0. .and. nb == 2) nb = 6
        !
        !==> Enter data on theta boundaries
        !
        select case (mbdcnd)
            case (1:2, 7)
                bm(1) = bm(1) - am(1)
                x = 2.0_wp*am(1)
                f(1, :n) = f(1, :n) - x*bda(:n)
            case (3:4, 8)
                bm(1) = bm(1) + am(1)
                x = dth*am(1)
                f(1, :n) = f(1, :n) + x*bda(:n)
        end select

        select case (mbdcnd)
            case (1, 4:5)
                bm(m) = bm(m) - cm(m)
                x = 2.0_wp*cm(m)
                f(m, :n) = f(m, :n) - x*bdb(:n)
            case (2:3, 6)
                bm(m) = bm(m) + cm(m)
                x = dth*cm(m)
                f(m, :n) = f(m, :n) - x*bdb(:n)
        end select

        select case (nb)
            case (1:2)
                bn(1) = bn(1) - an(1)
                x = 2.0_wp*an(1)
                f(:m, 1) = f(:m, 1) - x*bdc(:m)
            case (3:4)
                bn(1) = bn(1) + an(1)
                x = dr*an(1)
                f(:m, 1) = f(:m, 1) + x*bdc(:m)
        end select

        select case (nb)
            case (1, 4:5)
                bn(n) = bn(n) - cn(n)
                x = 2.0_wp*cn(n)
                f(:m, n) = f(:m, n) - x*bdd(:m)
            case (2:3, 6)
                bn(n) = bn(n) + cn(n)
                x = dr*cn(n)
                f(:m, n) = f(:m, n) - x*bdd(:m)
        end select

        pertrb = 0.0_wp

        case_block: block
            select case (mbdcnd)
                case (1:2, 4:5, 7)
                    exit case_block
                case (3, 6, 8:9)
                    select case (nb)
                        case (1:2, 4:5)
                            exit case_block
                        case (3, 6)
                            if (elmbda >= 0.0_wp) then
                                if (elmbda /= 0.0_wp) then
                                    ierr1 = 10
                                else
                                    isw = 2
                                    do i = 1, m
                                        x = 0.0_wp
                                        x = sum(f(i, :n))
                                        pertrb = pertrb + x*snth(i)
                                    end do
                                    x = 0.0_wp
                                    x = sum(rsq(:n))
                                    pertrb = 2.0_wp*(pertrb*sin(dth/2.))/(x*(cos(a) - cos(b)))
                                    do j = 1, n
                                        x = rsq(j)*pertrb
                                        f(:m, j) = f(:m, j) - x
                                    end do
                                end if
                            end if
                    end select
            end select
        end block case_block

        a2 = sum(f(:m, 1))
        a2 = a2/rsq(1)
        !
        !     initialize blktri
        !
        ierr1 = 0
        if (intl == 0) call blktrii(0, 1, n, an, bn, cn, 1, m, am, bm, cm, idimf, f, ierr1, w, wc)

        call blktrii(1, 1, n, an, bn, cn, 1, m, am, bm, cm, idimf, f, ierr1, w, wc)

        if (.not.(isw /=2 .or. c /= 0. .or. nbdcnd /= 2)) then
            a3 = 0.0_wp
            a1 = dot_product(snth(:m), f(:m, 1))
            a3 = sum(snth(:m))
            a1 = a1 + rsq(1)*a2/2

            if (mbdcnd == 3) a1=a1+(sin(b)*bdb(1)-sin(a)*bda(1))/(2.0_wp*(b-a))

            a1 = a1/a3
            a1 = bdc(1) - a1
            f(:m, :n) = f(:m, :n) + a1
        end if

    end subroutine hstcs1

end module module_hstcsp
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
!-----------------------------------------------------------------------
