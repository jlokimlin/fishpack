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
    public :: test_hwscsp

contains

    subroutine test_hwscsp()
        !
        !     file thwscsp.f
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
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer::intl, m, mbdcnd, n, nbdcnd, idimf, mp1, i, np1, j, ierror
        real , dimension(48, 33) :: f
        real , dimension(33) :: bdtf, bdts, bdrs, bdrf
        real , dimension(48) :: theta
        real , dimension(33) :: r
        real :: pi, dum, ts, tf, rs, rf, elmbda, dtheta, dr, ci4, pertrb, discretization_error, z, dphi, si
        !-----------------------------------------------
        !
        !     PROGRAM TO ILLUSTRATE THE USE OF hwscsp
        !
        !
        pi = acos( -1.0 )
        intl = 0
        ts = 0.
        tf = pi/2.
        m = 36
        mbdcnd = 6
        rs = 0.
        rf = 1.
        n = 32
        nbdcnd = 5
        elmbda = 0.
        idimf = 48
        !
        !     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
        !     BOUNDARY DATA AND THE RIGHT SIDE OF THE EQUATION.
        !
        mp1 = m + 1
        dtheta = tf/real(m)
        do i = 1, mp1
            theta(i) = real(i - 1)*dtheta
        end do
        np1 = n + 1
        dr = 1./real(n)
        do j = 1, np1
            r(j) = real(j - 1)*dr
        end do
        !
        !     GENERATE NORMAL DERIVATIVE DATA AT EQUATOR
        !
        bdtf(:np1) = 0.
        !
        !     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
        !
        do i = 1, mp1
            f(i, n+1) = COS(THETA(i))**4
        end do
        !
        !     COMPUTE RIGHT SIDE OF EQUATION
        !
        do i = 1, mp1
            ci4 = 12.*COS(THETA(i))**2
            f(i, :n) = ci4*R(:n)**2
        end do

        call hwscsp (intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, &
            nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, ierror, workspace)
        !
        !     COMPUTE DISCRETIZATION ERROR
        !
        discretization_error = 0.
        do i = 1, mp1
            ci4 = COS(THETA(i))**4
            do j = 1, n
                z = abs(F(i, j)-ci4*R(j)**4)
                discretization_error = max(z, discretization_error)
            end do
        end do
        !     Print earlier output from platforms with 32 and 64 bit floating point
        !     arithemtic followed by the output from this computer
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     hwscsp *** TEST RUN, EXAMPLE 1 *** '
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') '     ierror = 0,  discretization error = 7.9984E-4 '
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)') &
            '     ierror =', ierror, '     discretization error =', discretization_error
        !
        !     THE FOLLOWING PROGRAM ILLUSTRATES THE USE OF hwscsp TO SOLVE
        !     A THREE DIMENSIONAL PROBLEM WHICH HAS LONGITUDNAL DEPENDENCE
        !
        mbdcnd = 2
        nbdcnd = 1
        dphi = pi/72.
        elmbda = -2.*(1. - COS(dphi))/dphi**2
        !
        !     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
        !
        do i = 1, mp1
            f(i, n+1) = SIN(THETA(i))
        end do
        !
        !     COMPUTE RIGHT SIDE OF THE EQUATION
        !
        f(:mp1, :n) = 0.
        call hwscsp (intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, &
            nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, ierror, workspace)
        !
        !     COMPUTE DISCRETIZATION ERROR   (FOURIER COEFFICIENTS)
        !
        discretization_error = 0
        do i = 1, mp1
            si = SIN(THETA(i))
            do j = 1, np1
                z = abs(F(i, j)-R(j)*si)
                discretization_error = max(z, discretization_error)
            end do
        end do
        !
        !     Print earlier output from platforms with 32 and 64 bit floating point
        !     arithemtic followed by the output from this computer

        write( stdout, '(A)') ''
        write( stdout, '(A)') '     hwscsp *** TEST RUN, EXAMPLE 2 *** '
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') '     ierror = 0, discretization error = 5.8682E-5 '
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)') &
            '     ierror =', ierror, '     discretization error =', discretization_error


        ! release dynamically allocated workspace arrays
        call workspace%destroy()

    end subroutine test_hwscsp


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
        !                        LIBRARIES IN JANUARY 1980. Revised by John
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
        !***********************************************************************
        type (FishpackWorkspace) :: w
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer  :: INTL
        integer  :: M
        integer  :: MBDCND
        integer  :: N
        integer  :: NBDCND
        integer  :: IDIMF
        integer  :: ierror
        real  :: TS
        real  :: TF
        real  :: RS
        real  :: RF
        real  :: ELMBDA
        real  :: PERTRB
        real  :: BDTS(*)
        real  :: BDTF(*)
        real  :: BDRS(*)
        real  :: BDRF(*)
        real  :: F(IDIMF, *)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer :: I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, NCK, L, K, NP &
            , IRWK, ICWK, NP1, MP1
        real :: PI

        save I1, I2, I3, I4, I5, I6, I7, I8, I9, I10
        !-----------------------------------------------
        pi = acos( -1.0 )
        ierror = 0
        if (TS<0. .or. TF>PI) ierror = 1
        if (TS >= TF) ierror = 2
        if (M < 5) ierror = 3
        if (MBDCND<1 .or. MBDCND>9) ierror = 4
        if (RS < 0.) ierror = 5
        if (RS >= RF) ierror = 6
        if (N < 5) ierror = 7
        if (NBDCND<1 .or. NBDCND>6) ierror = 8
        if (ELMBDA > 0.) ierror = 9
        if (IDIMF < M + 1) ierror = 10
        if (ELMBDA/=0. .and. MBDCND>=5) ierror = 11
        if (ELMBDA/=0. .and. (NBDCND==5 .or. NBDCND==6)) ierror = 12
        if((MBDCND==5.or.MBDCND==6.or.MBDCND==9).and.TS/=0.)ierror=13
        if (MBDCND>=7 .and. TF/=PI) ierror = 14
        if(TS==0..and.(MBDCND==4.or.MBDCND==8.or.MBDCND==3))ierror=15
        if(TF==PI.and.(MBDCND==2.or.MBDCND==3.or.MBDCND==6))ierror=16
        if (NBDCND>=5 .and. RS/=0.) ierror = 17
        if (NBDCND>=5 .and. (MBDCND==1 .or. MBDCND==2 .or. MBDCND==5 .or. &
            MBDCND==7)) ierror = 18
        if (ierror/=0 .and. ierror/=9) return
        NCK = N
        go to (101, 103, 102, 103, 101, 103) NBDCND
101 continue
    NCK = NCK - 1
    go to 103
102 continue
    NCK = NCK + 1
103 continue
    L = 2
    K = 1
    L = L + L
    K = K + 1
    do while(NCK - L > 0)
        L = L + L
        K = K + 1
    end do
    L = L + L

    if (INTL == 0) then
        !          compute blktri work space lengths
        NP = NBDCND
        call w%get_block_tridiagonal_workpace_dimensions (N, M, IRWK, ICWK)
        NP1 = N + 1
        MP1 = M + 1
        I1 = (K - 2)*L + K + max(2*N, 6*M) + 13
        I2 = I1 + NP1
        I3 = I2 + NP1
        I4 = I3 + NP1
        I5 = I4 + NP1
        I6 = I5 + NP1
        I7 = I6 + MP1
        I8 = I7 + MP1
        I9 = I8 + MP1
        I10 = I9 + MP1
        !          set real and complex work space requirements
        IRWK = I10 + MP1
        ICWK = ICWK + 3*M
        !          allocate work space
        call w%create( irwk, icwk, ierror )
        !          return if allocation fails
        if (ierror == 20) return
    end if

    associate( &
        rew => w%rew, &
        cxw => w%cxw &
        )
        call hwscs1 (intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, nbdcnd, bdrs, &
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
    integer , intent (in) :: INTL
    integer , intent (in) :: M
    integer , intent (in) :: MBDCND
    integer , intent (in) :: N
    integer , intent (in) :: NBDCND
    integer  :: IDIMF
    integer  :: ierror
    real , intent (in) :: TS
    real , intent (in) :: TF
    real , intent (in) :: RS
    real , intent (in) :: RF
    real , intent (in) :: ELMBDA
    real , intent (out) :: PERTRB
    real , intent (in) :: BDTS(*)
    real , intent (in) :: BDTF(*)
    real , intent (in) :: BDRS(*)
    real , intent (in) :: BDRF(*)
    real  :: F(IDIMF, *)
    real  :: W(*)
    real , intent (in out) :: S(*)
    real  :: AN(*)
    real  :: BN(*)
    real  :: CN(*)
    real , intent (in out) :: R(*)
    real  :: AM(*)
    real  :: BM(*)
    real  :: CM(*)
    real , intent (in out) :: SINT(*)
    real , intent (in out) :: BMH(*)
    complex  :: WC(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer :: MP1, I, NP1, J, MP, NP, ITS, ITF, ITSP, ITFM, ICTR, JRS &
        , L, JRF, JRSP, JRFM, MUNK, NUNK, ISING, IFLG
    real :: PI, DUM, EPS, DTH, TDT, HDTH, SDTS, THETA, T1, DR, HDR, &
        TDR, DR2, CZR, AT, CT, WTS, WTF, AR, WTNM, YPS, CR, WRS, WRF, &
        WRZ, summation, R2, HNE, YHLD, RS2, RF2, RSQ, XP, YPH, XPS

    !-----------------------------------------------
    pi = acos( -1.0 )
    EPS = epsilon(DUM)
    MP1 = M + 1
    DTH = (TF - TS)/real(M)
    TDT = DTH + DTH
    HDTH = DTH/2.
    SDTS = 1./(DTH*DTH)
    do I = 1, MP1
        THETA = TS + real(I - 1)*DTH
        SINT(I) = SIN(THETA)
        if (SINT(I) == 0.) cycle
        T1 = SDTS/SINT(I)
        AM(I) = T1*SIN(THETA - HDTH)
        CM(I) = T1*SIN(THETA + HDTH)
        BM(I) = -(AM(I)+CM(I))
    end do
    NP1 = N + 1
    DR = (RF - RS)/real(N)
    HDR = DR/2.
    TDR = DR + DR
    DR2 = DR*DR
    CZR = 6.*DTH/(DR2*(COS(TS) - COS(TF)))
    do J = 1, NP1
        R(J) = RS + real(J - 1)*DR
        AN(J) = (R(J)-HDR)**2/DR2
        CN(J) = (R(J)+HDR)**2/DR2
        BN(J) = -(AN(J)+CN(J))
    end do
    MP = 1
    NP = 1
    !
    ! BOUNDARY CONDITION AT PHI=PS
    !
    go to (104, 104, 105, 105, 106, 106, 104, 105, 106) MBDCND
104 continue
    AT = AM(2)
    ITS = 2
    go to 107
105 continue
    AT = AM(1)
    ITS = 1
    CM(1) = CM(1) + AM(1)
    go to 107
106 continue
    ITS = 1
    BM(1) = -4.*SDTS
    CM(1) = -BM(1)
!
! BOUNDARY CONDITION AT PHI=PF
!
107 continue
    go to (108, 109, 109, 108, 108, 109, 110, 110, 110) MBDCND
108 continue
    CT = CM(M)
    ITF = M
    go to 111
109 continue
    CT = CM(M+1)
    AM(M+1) = AM(M+1) + CM(M+1)
    ITF = M + 1
    go to 111
110 continue
    ITF = M + 1
    AM(M+1) = 4.*SDTS
    BM(M+1) = -AM(M+1)
111 continue
    WTS = SINT(ITS+1)*AM(ITS+1)/CM(ITS)
    WTF = SINT(ITF-1)*CM(ITF-1)/AM(ITF)
    ITSP = ITS + 1
    ITFM = ITF - 1
    !
    ! BOUNDARY CONDITION AT R=RS
    !
    ICTR = 0
    select case (NBDCND)
        case default
            AR = AN(2)
            JRS = 2
        case (3:4)
            AR = AN(1)
            JRS = 1
            CN(1) = CN(1) + AN(1)
        case (5:6)
            JRS = 2
            ICTR = 1
            S(N) = AN(N)/BN(N)
            do J = 3, N
                L = N - J + 2
                S(L) = AN(L)/(BN(L)-CN(L)*S(L+1))
            end do
            S(2) = -S(2)
            do J = 3, N
                S(J) = -S(J)*S(J-1)
            end do
            WTNM = WTS + WTF
            do I = ITSP, ITFM
                WTNM = WTNM + SINT(I)
            end do
            YPS = CZR*WTNM*(S(2)-1.)
    end select
!
! BOUNDARY CONDITION AT R=RF
!
118 continue
    go to (119, 120, 120, 119, 119, 120) NBDCND
119 continue
    CR = CN(N)
    JRF = N
    go to 121
120 continue
    CR = CN(N+1)
    AN(N+1) = AN(N+1) + CN(N+1)
    JRF = N + 1
121 continue
    WRS = AN(JRS+1)*R(JRS)**2/CN(JRS)
    WRF = CN(JRF-1)*R(JRF)**2/AN(JRF)
    WRZ = AN(JRS)/CZR
    JRSP = JRS + 1
    JRFM = JRF - 1
    MUNK = ITF - ITS + 1
    NUNK = JRF - JRS + 1
    BMH(ITS:ITF) = BM(ITS:ITF)
    ISING = 0
    go to (132, 132, 123, 132, 132, 123) NBDCND
123 continue
    go to (132, 132, 124, 132, 132, 124, 132, 124, 124) MBDCND
124 continue
    if (ELMBDA >= 0.) then
        ISING = 1
        summation = WTS*WRS + WTS*WRF + WTF*WRS + WTF*WRF
        if (ICTR /= 0) then
            summation = summation + WRZ
        end if
        do J = JRSP, JRFM
            R2 = R(J)**2
            do I = ITSP, ITFM
                summation = summation + R2*SINT(I)
            end do
        end do
        do J = JRSP, JRFM
            summation = summation + (WTS + WTF)*R(J)**2
        end do
        do I = ITSP, ITFM
            summation = summation + (WRS + WRF)*SINT(I)
        end do
        HNE = summation
    end if
132 continue
    go to (133, 133, 133, 133, 134, 134, 133, 133, 134) MBDCND
133 continue
    BM(ITS) = BMH(ITS) + ELMBDA/SINT(ITS)**2
134 continue
    go to (135, 135, 135, 135, 135, 135, 136, 136, 136) MBDCND
135 continue
    BM(ITF) = BMH(ITF) + ELMBDA/SINT(ITF)**2
136 continue
    BM(ITSP:ITFM) = BMH(ITSP:ITFM) + ELMBDA/SINT(ITSP:ITFM)**2
    go to (138, 138, 140, 140, 142, 142, 138, 140, 142) MBDCND
138 continue
    F(2, JRS:JRF) = F(2, JRS:JRF) - AT*F(1, JRS:JRF)/R(JRS:JRF)**2
    go to 142
140 continue
    F(1, JRS:JRF) = F(1, JRS:JRF) + TDT*BDTS(JRS:JRF)*AT/R(JRS:JRF)**2
142 continue
    go to (143, 145, 145, 143, 143, 145, 147, 147, 147) MBDCND
143 continue
    F(M, JRS:JRF) = F(M, JRS:JRF) - CT*F(M+1, JRS:JRF)/R(JRS:JRF)**2
    go to 147
145 continue
    F(M+1, JRS:JRF)=F(M+1, JRS:JRF)-TDT*BDTF(JRS:JRF)*CT/R(JRS:JRF)**2
147 continue
    select case (NBDCND)
        case default
            if (MBDCND - 3 /= 0) go to 155
            YHLD = F(ITS, 1) - CZR/TDT*(SIN(TF)*BDTF(2)-SIN(TS)*BDTS(2))
            F(:MP1, 1) = YHLD
        case (1:2)
            RS2 = (RS + DR)**2
            F(ITS:ITF, 2) = F(ITS:ITF, 2) - AR*F(ITS:ITF, 1)/RS2
        case (3:4)
            F(ITS:ITF, 1) = F(ITS:ITF, 1) + TDR*BDRS(ITS:ITF)*AR/RS**2
    end select
155 continue
    go to (156, 158, 158, 156, 156, 158) NBDCND
156 continue
    RF2 = (RF - DR)**2
    F(ITS:ITF, N) = F(ITS:ITF, N) - CR*F(ITS:ITF, N+1)/RF2
    go to 160
158 continue
    F(ITS:ITF, N+1) = F(ITS:ITF, N+1) - TDR*BDRF(ITS:ITF)*CR/RF**2
160 continue
    PERTRB = 0.
    if (ISING /= 0) then
        summation = WTS*WRS*F(ITS, JRS) + WTS*WRF*F(ITS, JRF) + WTF*WRS*F(ITF, &
            JRS) + WTF*WRF*F(ITF, JRF)
        if (ICTR /= 0) then
            summation = summation + WRZ*F(ITS, 1)
        end if
        do J = JRSP, JRFM
            R2 = R(J)**2
            do I = ITSP, ITFM
                summation = summation + R2*SINT(I)*F(I, J)
            end do
        end do
        summation = summation + DOT_PRODUCT(R(JRSP:JRFM)**2, WTS*F(ITS, JRSP:JRFM)+ &
            WTF*F(ITF, JRSP:JRFM))
        summation = summation + DOT_PRODUCT(SINT(ITSP:ITFM), WRS*F(ITSP:ITFM, JRS)+ &
            WRF*F(ITSP:ITFM, JRF))
        PERTRB = summation/HNE
        F(:MP1, :NP1) = F(:MP1, :NP1) - PERTRB
    end if
    do J = JRS, JRF
        RSQ = R(J)**2
        F(ITS:ITF, J) = RSQ*F(ITS:ITF, J)
    end do
    IFLG = INTL
    call blktriI (IFLG, NP, NUNK, AN(JRS), BN(JRS), CN(JRS), MP, MUNK &
        , AM(ITS), BM(ITS), CM(ITS), IDIMF, F(ITS, JRS), ierror, W, WC)
    IFLG = IFLG + 1
    do while(IFLG - 1 == 0)
        call blktriI (IFLG, NP, NUNK, AN(JRS), BN(JRS), CN(JRS), MP, &
            MUNK, AM(ITS), BM(ITS), CM(ITS), IDIMF, F(ITS, JRS), ierror, &
            W, WC)
        IFLG = IFLG + 1
    end do
    if (NBDCND == 0) then
        F(:MP1, JRF+1) = F(:MP1, JRS)
    end if
    if (MBDCND == 0) then
        F(ITF+1, :NP1) = F(ITS, :NP1)
    end if
    XP = 0.
    if (ICTR /= 0) then
        if (ISING == 0) then
            summation = WTS*F(ITS, 2) + WTF*F(ITF, 2)
            summation = summation + DOT_PRODUCT(SINT(ITSP:ITFM), F(ITSP:ITFM, 2))
            YPH = CZR*summation
            XP = (F(ITS, 1)-YPH)/YPS
            do J = JRS, JRF
                XPS = XP*S(J)
                F(ITS:ITF, J) = F(ITS:ITF, J) + XPS
            end do
        end if
        F(:MP1, 1) = XP
    end if

end subroutine HWSCS1

end module module_hwscsp
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
