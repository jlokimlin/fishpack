module module_hwscsp

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_blktri, only: &
        blktrii

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hwscsp


contains


    subroutine hwscsp(intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, &
        nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, ierror, w)
        !
        !     file hwscsp.f
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
        !     SUBROUTINE hwscsp (INTL, TS, TF, M, MBDCND, BDTS, BDTF, RS, RF, N, NBDCND,
        !    +                   BDRS, BDRF, ELMBDA, F, IDIMF, PERTRB, ierror, W)
        !
        !
        ! DIMENSION OF           BDTS(N+1),     BDTF(N+1), BDRS(M+1), BDRF(M+1),
        ! ARGUMENTS              F(IDIMF, N+1)
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION
        !                        TO THE MODIFIED HELMHOLTZ EQUATION IN
        !                        SPHERICAL COORDINATES ASSUMING AXISYMMETRY
        !                        (NO DEPENDENCE ON LONGITUDE).  THE EQUATION
        !                        IS
        !
        !                          (1/R**2)(D/DR)((R**2)(D/DR)U) +
        !
        !                          (1/(R**2)SIN(THETA))(D/DTHETA)
        !
        !                          (SIN(THETA)(D/DTHETA)U) +
        !
        !                          (LAMBDA/(RSIN(THETA))**2)U = F(THETA, R).
        !
        !                        THIS TWO DIMENSIONAL MODIFIED HELMHOLTZ
        !                        EQUATION RESULTS FROM THE FOURIER TRANSFORM
        !                        OF THE THREE DIMENSIONAL POISSON EQUATION.
        !
        ! USAGE                  CALL hwscsp (INTL, TS, TF, M, MBDCND, BDTS, BDTF,
        !                                     RS, RF, N, NBDCND, BDRS, BDRF, ELMBDA,
        !                                     F, IDIMF, PERTRB, ierror, W)
        !
        ! ARGUMENTS
        ! ON INPUT               INTL
        !                          = 0  ON INITIAL ENTRY TO hwscsp OR IF ANY
        !                               OF THE ARGUMENTS RS, RF, N, NBDCND
        !                               ARE CHANGED FROM A PREVIOUS CALL.
        !                          = 1  IF RS, RF, N, NBDCND ARE ALL UNCHANGED
        !                               FROM PREVIOUS CALL TO hwscsp.
        !
        !                          NOTE:
        !                          A CALL WITH INTL=0 TAKES APPROXIMATELY
        !                          1.5 TIMES AS MUCH TIME AS A CALL WITH
        !                          INTL = 1  .  ONCE A CALL WITH INTL = 0
        !                          HAS BEEN MADE THEN SUBSEQUENT SOLUTIONS
        !                          CORRESPONDING TO DIFFERENT F, BDTS, BDTF,
        !                          BDRS, BDRF CAN BE OBTAINED FASTER WITH
        !                          INTL = 1 SINCE INITIALIZATION IS NOT
        !                          REPEATED.
        !
        !                        TS, TF
        !                          THE RANGE OF THETA (COLATITUDE), I.E.,
        !                          TS .LE. THETA .LE. TF. TS MUST BE LESS
        !                          THAN TF.  TS AND TF ARE IN RADIANS. A TS OF
        !                          ZERO CORRESPONDS TO THE NORTH POLE AND A
        !                          TF OF PI CORRESPONDS TO THE SOUTH POLE.
        !
        !                          **** IMPORTANT ****
        !
        !                          IF TF IS EQUAL TO PI THEN IT MUST BE
        !                          COMPUTED USING THE STATEMENT
        !                          TF = PI_MACH(DUM). THIS INSURES THAT TF
        !                          IN THE USER'S PROGRAM IS EQUAL TO PI IN
        !                          THIS PROGRAM WHICH PERMITS SEVERAL TESTS
        !                          OF THE  INPUT PARAMETERS THAT OTHERWISE
        !                          WOULD NOT BE POSSIBLE.
        !
        !                        M
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL (TS, TF) IS SUBDIVIDED.
        !                          HENCE, THERE WILL BE M+1 GRID POINTS
        !                          IN THE THETA-DIRECTION GIVEN BY
        !                          THETA(K) = (I-1)DTHETA+TS FOR
        !                          I = 1, 2, ..., M+1, WHERE DTHETA = (TF-TS)/M
        !                          IS THE PANEL WIDTH.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITION
        !                          AT THETA = TS AND  THETA = TF.
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = TS AND THETA = TF.
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = TS AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO THETA IS
        !                               SPECIFIED AT THETA = TF
        !                               (SEE NOTE 2 BELOW).
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               AT THETA = TS AND THETA = TF
        !                               (SEE NOTES 1, 2 BELOW).
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               AT THETA = TS (SEE NOTE 1 BELOW) AND
        !                               SOLUTION IS SPECIFIED AT THETA = TF.
        !                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = TS = 0 AND THE SOLUTION IS
        !                                SPECIFIED AT THETA = TF.
        !                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = TS = 0 AND THE DERIVATIVE
        !                               OF THE SOLUTION WITH RESPECT TO THETA
        !                               IS SPECIFIED AT THETA = TF
        !                               (SEE NOTE 2 BELOW).
        !                          = 7  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = TS AND THE SOLUTION IS
        !                                UNSPECIFIED AT THETA = TF = PI.
        !                          = 8  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               AT THETA = TS (SEE NOTE 1 BELOW)
        !                               AND THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = TF = PI.
        !                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = TS = 0 AND THETA = TF = PI.
        !
        !                          NOTE 1:
        !                          IF TS = 0, DO NOT USE MBDCND = 3, 4, OR 8,
        !                          BUT INSTEAD USE MBDCND = 5, 6, OR 9  .
        !
        !                          NOTE 2:
        !                          IF TF = PI, DO NOT USE MBDCND = 2, 3, OR 6,
        !                          BUT INSTEAD USE MBDCND = 7, 8, OR 9  .
        !
        !                        BDTS
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
        !                          THE SOLUTION WITH RESPECT TO THETA AT
        !                          THETA = TS.  WHEN MBDCND = 3, 4, OR 8,
        !
        !                            BDTS(J) = (D/DTHETA)U(TS, R(J)),
        !                            J = 1, 2, ..., N+1  .
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDTS IS
        !                          A DUMMY VARIABLE.
        !
        !                        BDTF
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
        !                          THE SOLUTION WITH RESPECT TO THETA AT
        !                          THETA = TF.  WHEN MBDCND = 2, 3, OR 6,
        !
        !                          BDTF(J) = (D/DTHETA)U(TF, R(J)),
        !                          J = 1, 2, ..., N+1  .
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDTF IS
        !                          A DUMMY VARIABLE.
        !
        !                        RS, RF
        !                          THE RANGE OF R, I.E., RS .LE. R .LT. RF.
        !                          RS MUST BE LESS THAN RF.  RS MUST BE
        !                          NON-NEGATIVE.
        !
        !                        N
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL (RS, RF) IS SUBDIVIDED.
        !                          HENCE, THERE WILL BE N+1 GRID POINTS IN THE
        !                          R-DIRECTION GIVEN BY R(J) = (J-1)DR+RS
        !                          FOR J = 1, 2, ..., N+1, WHERE DR = (RF-RS)/N
        !                          IS THE PANEL WIDTH.
        !                          N MUST BE GREATER THAN 2
        !
        !                        NBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITION
        !                          AT R = RS AND R = RF.
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               R = RS AND R = RF.
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               R = RS AND THE DERIVATIVE
        !                               OF THE SOLUTION WITH RESPECT TO R
        !                               IS SPECIFIED AT R = RF.
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO R IS SPECIFIED AT
        !                               R = RS AND R = RF.
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO R IS SPECIFIED AT
        !                               RS AND THE SOLUTION IS SPECIFIED AT
        !                               R = RF.
        !                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
        !                               R = RS = 0 (SEE NOTE BELOW)  AND THE
        !                               SOLUTION IS SPECIFIED AT R = RF.
        !                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
        !                               R = RS = 0 (SEE NOTE BELOW) AND THE
        !                               DERIVATIVE OF THE SOLUTION WITH
        !                               RESPECT TO R IS SPECIFIED AT R = RF.
        !
        !                          NOTE:
        !                          NBDCND = 5 OR 6 CANNOT BE USED WITH
        !                          MBDCND = 1, 2, 4, 5, OR 7.  THE FORMER
        !                          INDICATES THAT THE SOLUTION IS UNSPECIFIED
        !                          AT R = 0, THE LATTER INDICATES THAT THE
        !                          SOLUTION IS SPECIFIED).
        !                          USE INSTEAD   NBDCND = 1 OR 2  .
        !
        !                        BDRS
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
        !                          THE SOLUTION WITH RESPECT TO R AT R = RS.
        !
        !                          WHEN NBDCND = 3 OR 4,
        !                            BDRS(I) = (D/DR)U(THETA(I), RS),
        !                            I = 1, 2, ..., M+1  .
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDRS IS
        !                          A DUMMY VARIABLE.
        !
        !                        BDRF
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
        !                          THAT SPECIFIES THE VALUES OF THE
        !                          DERIVATIVE OF THE SOLUTION WITH RESPECT
        !                          TO R AT R = RF.
        !
        !                          WHEN NBDCND = 2, 3, OR 6,
        !                            BDRF(I) = (D/DR)U(THETA(I), RF),
        !                            I = 1, 2, ..., M+1  .
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDRF IS
        !                          A DUMMY VARIABLE.
        !
        !                        ELMBDA
        !                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
        !                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
        !                          MAY NOT EXIST.  HOWEVER, hwscsp WILL
        !                          ATTEMPT TO FIND A SOLUTION.  IF NBDCND = 5
        !                          OR 6 OR  MBDCND = 5, 6, 7, 8, OR 9, ELMBDA
        !                          MUST BE ZERO.
        !
        !                        F
        !                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
        !                          LEAST (M+1)*(N+1), SPECIFYING VALUES OF THE
        !                          RIGHT SIDE OF THE HELMHOLTZ EQUATION AND
        !                          BOUNDARY VALUES (IF ANY).
        !
        !                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
        !                          FOR I = 2, 3, ..., M AND J = 2, 3, ..., N
        !                          F(I, J) = F(THETA(I), R(J)).
        !
        !                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
        !                          FOR J=1, 2, ..., N+1,  I=1, 2, ..., M+1,
        !
        !                          MBDCND   F(1, J)            F(M+1, J)
        !                          ------   ----------        ----------
        !
        !                            1      U(TS, R(J))        U(TF, R(J))
        !                            2      U(TS, R(J))        F(TF, R(J))
        !                            3      F(TS, R(J))        F(TF, R(J))
        !                            4      F(TS, R(J))        U(TF, R(J))
        !                            5      F(0, R(J))         U(TF, R(J))
        !                            6      F(0, R(J))         F(TF, R(J))
        !                            7      U(TS, R(J))        F(PI, R(J))
        !                            8      F(TS, R(J))        F(PI, R(J))
        !                            9      F(0, R(J))         F(PI, R(J))
        !
        !                            NBDCND   F(I, 1)            F(I, N+1)
        !                            ------   --------------    --------------
        !
        !                              1      U(THETA(I), RS)    U(THETA(I), RF)
        !                              2      U(THETA(I), RS)    F(THETA(I), RF)
        !                              3      F(THETA(I), RS)    F(THETA(I), RF)
        !                              4      F(THETA(I), RS)    U(THETA(I), RF)
        !                              5      F(TS, 0)           U(THETA(I), RF)
        !                              6      F(TS, 0)           F(THETA(I), RF)
        !
        !                          NOTE:
        !                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
        !                          U AND THE RIGHT SIDE F AT A CORNER THEN
        !                          THE SOLUTION MUST BE SPECIFIED.
        !
        !                        IDIMF
        !                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
        !                          F AS IT APPEARS IN THE PROGRAM CALLING
        !                          hwscsp.  THIS PARAMETER IS USED TO SPECIFY
        !                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
        !                          BE AT LEAST M+1  .
        !
        !                        W
        !                          A fortran 90 derived TYPE (FishpackWorkspace) variable
        !                          that must be declared by the user.  The first
        !                          two declarative statements in the user program
        !                          calling SEPELI must be:
        !
        !                               use type_FishpackWorkspace
        !                               TYPE (FishpackWorkspace) :: W
        !
        !                          The first statement makes the fishpack module
        !                          defined in the file "fish.f" available to the
        !                          user program calling hwscsp.  The second statement
        !                          declares a derived type variable (defined in
        !                          the module "fish.f") which is used internally
        !                          in hwscsp to dynamically allocate real and complex
        !                          work space used in solution.  An error flag
        !                          (ierror = 20) is set if the required work space
        !                          allocation fails (for example if N, M are too large)
        !                          Real and complex values are set in the components
        !                          of W on a initial (INTL=0) call to hwscsp.  These
        !                          must be preserved on non-initial calls (INTL=1)
        !                          to hwscsp.  This eliminates redundant calculations
        !                          and saves compute time.
        !               ****       IMPORTANT!  The user program calling hwscsp should
        !                          include the statement:
        !
        !                               CALL FISHFIN(W)
        !
        !                          after the final approximation is generated by
        !                          hwscsp.  The will deallocate the real and complex
        !                          work space of W.  Failure to include this statement
        !                          could result in serious memory leakage.
        !
        !
        ! ON OUTPUT              F
        !                          CONTAINS THE SOLUTION U(I, J) OF THE FINITE
        !                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
        !                          (THETA(I), R(J)),  I = 1, 2, ..., M+1,
        !                                            J = 1, 2, ..., N+1  .
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
        !                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
        !                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
        !                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
        !                          CALCULATED AND SUBTRACTED FROM F, WHICH
        !                          ENSURES THAT A SOLUTION EXISTS.  hwscsp
        !                          THEN COMPUTES THIS SOLUTION, WHICH IS A
        !                          LEAST SQUARES SOLUTION TO THE ORIGINAL
        !                          APPROXIMATION. THIS SOLUTION IS NOT UNIQUE
        !                          AND IS UNNORMALIZED. THE VALUE OF PERTRB
        !                          SHOULD BE SMALL COMPARED TO THE RIGHT SIDE
        !                          F. OTHERWISE , A SOLUTION IS OBTAINED TO
        !                          AN ESSENTIALLY DIFFERENT PROBLEM. THIS
        !                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
        !                          THAT A MEANINGFUL SOLUTION HAS BEEN OBTAINED.
        !
        !                        ierror
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 10,
        !                          A SOLUTION IS NOT ATTEMPTED.
        !
        !                          = 1  TS.LT.0. OR TF.GT.PI
        !                          = 2  TS.GE.TF
        !                          = 3  M.LT.5
        !                          = 4  MBDCND.LT.1 OR MBDCND.GT.9
        !                          = 5  RS.LT.0
        !                          = 6  RS.GE.RF
        !                          = 7  N.LT.5
        !                          = 8  NBDCND.LT.1 OR NBDCND.GT.6
        !                          = 9  ELMBDA.GT.0
        !                          = 10 IDIMF.LT.M+1
        !                          = 11 ELMBDA.NE.0 AND MBDCND.GE.5
        !                          = 12 ELMBDA.NE.0 AND NBDCND EQUALS 5 OR 6
        !                          = 13 MBDCND EQUALS 5, 6 OR 9 AND TS.NE.0
        !                          = 14 MBDCND.GE.7 AND TF.NE.PI
        !                          = 15 TS.EQ.0 AND MBDCND EQUALS 3, 4 OR 8
        !                          = 16 TF.EQ.PI AND MBDCND EQUALS 2, 3 OR 6
        !                          = 17 NBDCND.GE.5 AND RS.NE.0
        !                          = 18 NBDCND.GE.5 AND MBDCND EQUALS 1, 2, 4, 5 OR
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space in the derived type
        !                               (FishpackWorkspace) variable W fails (e.g.,
        !                               if N, M are too large for the platform used)
        !
        !                          SINCE THIS IS THE ONLY MEANS OF INDICATING
        !                          A POSSLIBY INCORRECT CALL TO hwscsp, THE
        !                          USER SHOULD TEST ierror AFTER A CALL.
        !
        !                        W
        !                          The derived type (FishpackWorkspace) variable W
        !                          contains real and complex values that must not
        !                          be destroyed if hwscsp is called again with
        !                          INTL=1.
        !
        ! SPECIAL CONDITIONS     NONE
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED files         fish.f, blktri.f, comf.f
        !
        ! LANGUAGE               FORTRAN 90
        !
        ! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
        !                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
        !                        LIBRARIES IN January 1980. Revised by John
        !                        Adams in June 2004 using Fortran 90 dynamically
        !                        allocated work space and derived datat types
        !                        to eliminate mixed mode conflicts in the earlier
        !                        versions.
        !
        ! PORTABILITY            FORTRAN 90
        !
        ! ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
        !                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
        !                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
        !                        AND THEN CALLS blktri TO SOLVE THE SYSTEM.
        !
        ! REFERENCES             SWARZTRAUBER, P. AND R. SWEET, "EFFICIENT
        !                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
        !                        ELLIPTIC EQUATIONS"
        !                          NCAR TN/IA-109, JULY, 1975, 138 PP.
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: intl
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: mbdcnd
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: nbdcnd
        integer (ip), intent (in)     :: idimf
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in)     :: ts
        real (wp),    intent (in)     :: tf
        real (wp),    intent (in)     :: rs
        real (wp),    intent (in)     :: rf
        real (wp),    intent (in)     :: elmbda
        real (wp),    intent (out)    :: pertrb
        real (wp),    intent (in)     :: bdts(*)
        real (wp),    intent (in)     :: bdtf(*)
        real (wp),    intent (in)     :: bdrs(*)
        real (wp),    intent (in)     :: bdrf(*)
        real (wp),    intent (in out) :: f(idimf, *)
        class (fishpackworkspace)     :: w
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip), save   :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10
        integer (ip)         :: nck, l, k, np, irwk, icwk, np1, mp1
        real (wp), parameter :: PI = acos(-1.0_wp)
        !-----------------------------------------------

        ! Initialize error flag
        ierror = 0

        if (ts < 0.0_wp .or. tf > PI) then
            ierror = 1
        end if

        if (ts >= tf) then
            ierror = 2
        end if

        if (m < 5) then
            ierror = 3
        end if

        if (mbdcnd < 1 .or. mbdcnd > 9) then
            ierror = 4
        end if

        if (rs < 0.0_wp) then
            ierror = 5
        end if

        if (rs >= rf) then
            ierror = 6
        end if

        if (n < 5) then
            ierror = 7
        end if

        if (nbdcnd < 1 .or. nbdcnd > 6) then
            ierror = 8
        end if

        if (elmbda > 0.0_wp) then
            ierror = 9
        end if

        if (idimf < m + 1) then
            ierror = 10
        end if

        if (elmbda /= 0.0_wp .and. mbdcnd >= 5) then
            ierror = 11
        end if

        if (elmbda /= 0.0_wp .and. (nbdcnd == 5 .or. nbdcnd == 6)) then
            ierror = 12
        end if

        if ((mbdcnd == 5 .or. mbdcnd == 6 .or. mbdcnd == 9) .and. ts /= 0.0_wp) then
            ierror=13
        end if

        if (mbdcnd >= 7 .and. tf /= PI) then
            ierror = 14
        end if

        if (ts == 0.0_wp .and. (mbdcnd == 4.or.mbdcnd == 8 .or. mbdcnd == 3)) then
            ierror = 15
        end if

        if (tf == PI .and. (mbdcnd == 2 .or. mbdcnd == 3 .or. mbdcnd == 6)) then
            ierror = 16
        end if

        if (nbdcnd >= 5 .and. rs /= 0.0_wp) then
            ierror = 17
        end if

        if ( &
            (nbdcnd >= 5) &
            .and. &
            (mbdcnd == 1 .or. mbdcnd == 2 .or. mbdcnd == 5 .or. mbdcnd == 7) &
            ) then
            ierror = 18
        end if

        if (ierror /= 0 .and. ierror /= 9) then
            return
        end if

        nck = n

        select case (nbdcnd)
            case (1, 5)
                nck = nck - 1
            case (3)
                nck = nck + 1
        end select

        l = 2
        k = 1
        l = l + l
        k = k + 1

        do while (nck - l > 0)
            l = l + l
            k = k + 1
        end do

        l = l + l

        if (intl == 0) then
            !
            !==> Compute blktri work space lengths
            !
            np = nbdcnd

            call w%get_block_tridiagonal_workpace_dimensions (N, M, IRWK, ICWK)

            np1 = n + 1
            mp1 = m + 1
            i1 = (k - 2)*l + k + max(2*n, 6*m) + 13
            i2 = i1 + np1
            i3 = i2 + np1
            i4 = i3 + np1
            i5 = i4 + np1
            i6 = i5 + np1
            i7 = i6 + mp1
            i8 = i7 + mp1
            i9 = i8 + mp1
            i10 = i9 + mp1
            !
            !==> Set real and complex work space requirements
            !
            irwk = i10 + mp1
            icwk = icwk + 3*m
            !
            !==> Allocate memory
            !
            call w%create(irwk, icwk, ierror)
            !
            !==> return if allocation fails
            !
            if (ierror == 20) then
                return
            end if
        end if

        associate( &
            rew => w%real_workspace, &
            cxw => w%complex_workspace &
            )
            call hwscs1(intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, nbdcnd, bdrs, &
                bdrf, elmbda, f, idimf, pertrb, rew, cxw, rew(i1), rew(i2), &
                rew(i3), rew(i4), rew(i5), rew(i6), rew(i7), rew(i8), &
                rew(i9), rew(i10), ierror)
        end associate

    end subroutine hwscsp


    subroutine hwscs1(intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, &
        nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, w, wc, s, an, bn &
        , cn, r, am, bm, cm, sint, bmh, ierror)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: intl
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: mbdcnd
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: nbdcnd
        integer (ip), intent (in)     :: idimf
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in)     :: ts
        real (wp),    intent (in)     :: tf
        real (wp),    intent (in)     :: rs
        real (wp),    intent (in)     :: rf
        real (wp),    intent (in)     :: elmbda
        real (wp),    intent (out)    :: pertrb
        real (wp),    intent (in)     :: bdts(*)
        real (wp),    intent (in)     :: bdtf(*)
        real (wp),    intent (in)     :: bdrs(*)
        real (wp),    intent (in)     :: bdrf(*)
        real (wp),    intent (in out) :: f(idimf, *)
        real (wp)                     :: w(*)
        real (wp),    intent (in out) :: s(*)
        real (wp)                     :: an(*)
        real (wp)                     :: bn(*)
        real (wp)                     :: cn(*)
        real (wp),    intent (in out) :: r(*)
        real (wp)                     :: am(*)
        real (wp)                     :: bm(*)
        real (wp)                     :: cm(*)
        real (wp),    intent (in out) :: sint(*)
        real (wp),    intent (in out) :: bmh(*)
        complex (wp)                  :: wc(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip)         :: hack_counter, mp1, i, np1, j, mp, np
        integer (ip)         :: its, itf, itsp, itfm, ictr, jrs
        integer (ip)         :: l, jrf, jrsp, jrfm, munk, nunk, ising, iflg
        real (wp)            :: dum, eps, dth, tdt, hdth, sdts
        real (wp)            :: theta, t1, dr, hdr
        real (wp)            :: tdr, dr2, czr, at, ct, wts, wtf
        real (wp)            :: ar, wtnm, yps, cr, wrs, wrf
        real (wp)            :: wrz, summation, r2, hne, yhld
        real (wp)            :: rs2, rf2, rsq, xp, yph, xps
        !-----------------------------------------------

        eps = epsilon(dum)
        mp1 = m + 1
        dth = (tf - ts)/m
        tdt = dth + dth
        hdth = dth/2
        sdts = 1.0_wp/(dth**2)

        do i = 1, mp1
            theta = ts + real(i - 1, kind=wp)*dth
            sint(i) = sin(theta)
            if (sint(i) == 0.0_wp) then
                cycle
            end if
            t1 = sdts/sint(i)
            am(i) = t1*sin(theta - hdth)
            cm(i) = t1*sin(theta + hdth)
            bm(i) = -(am(i)+cm(i))
        end do

        np1 = n + 1
        dr = (rf - rs)/n
        hdr = dr/2
        tdr = dr + dr
        dr2 = dr**2
        czr = 6.0_wp*dth/(dr2*(cos(ts) - cos(tf)))

        do j = 1, np1
            r(j) = rs + real(j - 1, kind=wp)*dr
            an(j) = (r(j)-hdr)**2/dr2
            cn(j) = (r(j)+hdr)**2/dr2
            bn(j) = -(an(j)+cn(j))
        end do

        mp = 1
        np = 1

        !
        !==> Boundary condition at phi=ps
        !
        select case (mbdcnd)
            case (1:2,7)
                at = am(2)
                its = 2
            case (3:4,8)
                at = am(1)
                its = 1
                cm(1) = cm(1) + am(1)
            case (5:6,9)
                its = 1
                bm(1) = -4.0_wp*sdts
                cm(1) = -bm(1)
        end select

        !
        !==> Boundary condition at phi=pf
        !
        select case (mbdcnd)
            case (1,4:5)
                ct = cm(m)
                itf = m
            case (2:3,6)
                ct = cm(m+1)
                am(m+1) = am(m+1) + cm(m+1)
                itf = m + 1
            case (7:9)
                itf = m + 1
                am(m+1) = 4.0_wp*sdts
                bm(m+1) = -am(m+1)
        end select

        wts = sint(its+1)*am(its+1)/cm(its)
        wtf = sint(itf-1)*cm(itf-1)/am(itf)
        itsp = its + 1
        itfm = itf - 1
        !
        ! Boundary condition at r=rs
        !
        ictr = 0
        select case (nbdcnd)
            case default
                ar = an(2)
                jrs = 2
            case (3:4)
                ar = an(1)
                jrs = 1
                cn(1) = cn(1) + an(1)
            case (5:6)
                jrs = 2
                ictr = 1
                s(n) = an(n)/bn(n)
                do j = 3, n
                    l = n - j + 2
                    s(l) = an(l)/(bn(l)-cn(l)*s(l+1))
                end do
                s(2) = -s(2)
                do j = 3, n
                    s(j) = -s(j)*s(j-1)
                end do
                wtnm = wts + wtf
                do i = itsp, itfm
                    wtnm = wtnm + sint(i)
                end do
                yps = czr*wtnm*(s(2)-1.0_wp)
        end select

        !
        !==> Boundary condition at r=rf
        !
        select case (nbdcnd)
            case (1,4:5)
                cr = cn(n)
                jrf = n
            case (2:3,6)
                cr = cn(n+1)
                an(n+1) = an(n+1) + cn(n+1)
                jrf = n + 1
        end select

        wrs = an(jrs+1)*r(jrs)**2/cn(jrs)
        wrf = cn(jrf-1)*r(jrf)**2/an(jrf)
        wrz = an(jrs)/czr
        jrsp = jrs + 1
        jrfm = jrf - 1
        munk = itf - its + 1
        nunk = jrf - jrs + 1
        bmh(its:itf) = bm(its:itf)
        ising = 0

        if  (nbdcnd == 3 .or. nbdcnd == 6 ) then
            select case (mbdcnd)
                case (3,6,8:9)
                    if (elmbda >= 0.0_wp) then
                        ising = 1
                        summation = wts*wrs + wts*wrf + wtf*wrs + wtf*wrf
                        if (ictr /= 0) then
                            summation = summation + wrz
                        end if
                        do j = jrsp, jrfm
                            r2 = r(j)**2
                            do i = itsp, itfm
                                summation = summation + r2*sint(i)
                            end do
                        end do
                        do j = jrsp, jrfm
                            summation = summation + (wts + wtf)*r(j)**2
                        end do
                        do i = itsp, itfm
                            summation = summation + (wrs + wrf)*sint(i)
                        end do
                        hne = summation
                    end if
            end select
        end if

        select case (mbdcnd)
            case (1:4, 7:8)
                bm(its) = bmh(its) + elmbda/sint(its)**2
        end select

        select case (mbdcnd)
            case (1:6)
                bm(itf) = bmh(itf) + elmbda/sint(itf)**2
        end select

        bm(itsp:itfm) = bmh(itsp:itfm) + elmbda/sint(itsp:itfm)**2

        select case (mbdcnd)
            case (1:2, 7)
                f(2, jrs:jrf) = f(2, jrs:jrf) - at*f(1, jrs:jrf)/r(jrs:jrf)**2
            case (3:4, 8)
                f(1, jrs:jrf) = f(1, jrs:jrf) + tdt*bdts(jrs:jrf)*at/r(jrs:jrf)**2
        end select

        select case (mbdcnd)
            case (1, 4:5)
                f(m, jrs:jrf) = f(m, jrs:jrf) - ct*f(m+1, jrs:jrf)/r(jrs:jrf)**2
            case (2:3, 6)
                f(m+1, jrs:jrf)=f(m+1, jrs:jrf)-tdt*bdtf(jrs:jrf)*ct/r(jrs:jrf)**2
        end select

        !
        !==> GCC 5.1 doesn't support exit for the select case construct yet
        !    only do loops are fully supported
        !
        dumb_hack: do hack_counter = 1, 1
            select case (nbdcnd)
                case default
                    if (mbdcnd - 3 /= 0) then
                        exit dumb_hack
                    end if
                    yhld = f(its, 1) - czr/tdt*(sin(tf)*bdtf(2)-sin(ts)*bdts(2))
                    f(:mp1, 1) = yhld
                case (1:2)
                    rs2 = (rs + dr)**2
                    f(its:itf, 2) = f(its:itf, 2) - ar*f(its:itf, 1)/rs2
                case (3:4)
                    f(its:itf, 1) = f(its:itf, 1) + tdr*bdrs(its:itf)*ar/rs**2
            end select
        end do dumb_hack

        select case (nbdcnd)
            case (1, 4:5)
                rf2 = (rf - dr)**2
                f(its:itf, n) = f(its:itf, n) - cr*f(its:itf, n+1)/rf2
            case (2:3, 6)
                f(its:itf, n+1) = f(its:itf, n+1) - tdr*bdrf(its:itf)*cr/rf**2
        end select

        pertrb = 0.0_wp

        if (ising /= 0) then
            summation = wts*wrs*f(its, jrs) + wts*wrf*f(its, jrf) + wtf*wrs*f(itf, &
                jrs) + wtf*wrf*f(itf, jrf)

            if (ictr /= 0) then
                summation = summation + wrz*f(its, 1)
            end if

            do j = jrsp, jrfm
                r2 = r(j)**2
                do i = itsp, itfm
                    summation = summation + r2*sint(i)*f(i, j)
                end do
            end do

            summation = summation &
                + dot_product(r(jrsp:jrfm)**2, wts*f(its, jrsp:jrfm)+ &
                wtf*f(itf, jrsp:jrfm))

            summation = summation &
                + dot_product(sint(itsp:itfm), wrs*f(itsp:itfm, jrs)+ &
                wrf*f(itsp:itfm, jrf))

            pertrb = summation/hne

            f(:mp1, :np1) = f(:mp1, :np1) - pertrb

        end if

        do j = jrs, jrf
            rsq = r(j)**2
            f(its:itf, j) = rsq*f(its:itf, j)
        end do

        iflg = intl

        call blktrii(iflg, np, nunk, an(jrs), bn(jrs), cn(jrs), mp, munk, &
            am(its), bm(its), cm(its), idimf, f(its, jrs), ierror, w, wc)

        iflg = iflg + 1

        do while(iflg - 1 == 0)
            call blktrii(iflg, np, nunk, an(jrs), bn(jrs), cn(jrs), mp, &
                munk, am(its), bm(its), cm(its), idimf, f(its, jrs), ierror, &
                w, wc)
            iflg = iflg + 1
        end do

        if (nbdcnd == 0) then
            f(:mp1, jrf+1) = f(:mp1, jrs)
        end if

        if (mbdcnd == 0) then
            f(itf+1, :np1) = f(its, :np1)
        end if

        xp = 0.0_wp

        if (ictr /= 0) then
            if (ising == 0) then
                summation = wts*f(its, 2) + wtf*f(itf, 2)
                summation = summation + dot_product(sint(itsp:itfm), f(itsp:itfm, 2))
                yph = czr*summation
                xp = (f(its, 1)-yph)/yps
                do j = jrs, jrf
                    xps = xp*s(j)
                    f(its:itf, j) = f(its:itf, j) + xps
                end do
            end if
            f(:mp1, 1) = xp
        end if

    end subroutine hwscs1

end module module_hwscsp
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
