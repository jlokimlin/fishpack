module module_hwsssp

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_genbun, only: &
        genbunn

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hwsssp
    public :: hwsssp_unit_test

contains

    subroutine hwsssp_unit_test()
        !
        !     file thwsssp.f
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
        !     PROGRAM TO ILLUSTRATE THE USE OF HWSSSP
        !
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: m, mbdcnd, n, nbdcnd, idimf, mp1, i, np1, j, ierror
        real , dimension(19, 73) :: f
        real , dimension(73) :: bdtf
        real , dimension(19) :: sint
        real , dimension(73) :: sinp
        real :: pi, ts, tf, ps, pf, elmbda, dtheta, dphi, bdts(1), bdps(1), bdpf(1) &
            , pertrb, err, z

        !        integer :: m, mbdcnd, n, nbdcnd, idimf, mp1, i, np1, j, ierror
        !        real , dimension(19, 73) :: f
        !        real , dimension(73) :: bdtf, bdts, bdps, bdpf
        !        real , dimension(19) :: sint
        !        real , dimension(73) :: sinp
        !        real :: pi, ts, tf, ps, pf, elmbda, dtheta, dphi, pertrb, err, z
        !-----------------------------------------------
        pi = acos( -1.0 )
        ts = 0
        tf = pi/2.
        m = 18
        mbdcnd = 6
        ps = 0
        pf = pi + pi
        n = 72
        nbdcnd = 0
        elmbda = 0.
        idimf = 19
        !
        !     GENERATE SINES FOR USE IN SUBSEQUENT COMPUTATIONS
        !
        dtheta = tf/real(m)
        mp1 = m + 1
        do i = 1, mp1
            sint(i) = SIN(real(i - 1)*dtheta)
        end do
        dphi = (pi + pi)/real(n)
        np1 = n + 1
        do j = 1, np1
            sinp(j) = SIN(real(j - 1)*dphi)
        end do
        !
        !     COMPUTE RIGHT SIDE OF EQUATION AND STORE IN F
        !
        do j = 1, np1
            f(:mp1, j) = 2. - 6.*(SINT(:mp1)*SINP(j))**2
        end do
        !
        !     STORE DERIVATIVE DATA AT THE EQUATOR
        !
        bdtf(:np1) = 0.
        !
        call HWSSSP(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
            bdps, bdpf, elmbda, f, idimf, pertrb, ierror)
        !
        !     COMPUTE DISCRETIZATION ERROR. SINCE PROBLEM IS SINGULAR, THE
        !     SOLUTION MUST BE NORMALIZED.
        !
        err = 0.0
        do j = 1, np1
            do i = 1, mp1
                z = abs(F(i, j)-(SINT(i)*SINP(j))**2-F(1, 1))
                err = max(z, err)
            end do
        end do

        write( *, *) ''
        write( *, *) '    HWSSSP TEST RUN *** '
        write( *, *) &
            '    Previous 64 bit floating point arithmetic result '
        write( *, *) '    IERROR = 0,  Discretization Error = 3.38107E-3'

        write( *, *) '    The output from your computer is: '
        write( *, *) '    IERROR =', ierror, ' Discretization Error = ', &
            err

    end subroutine hwsssp_unit_test

    subroutine HWSSSP(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
        bdps, bdpf, elmbda, f, idimf, pertrb, ierror)
        !
        !     file hwsssp.f
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
        !     SUBROUTINE HWSSSP (TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, NBDCND, BDPS,
        !    +                   BDPF, ELMBDA, F, IDIMF, PERTRB, IERROR)
        !
        ! DIMENSION OF           BDTS(N+1),    BDTF(N+1), BDPS(M+1), BDPF(M+1),
        ! ARGUMENTS              F(IDIMF, N+1)
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION TO
        !                        THE HELMHOLTZ EQUATION IN SPHERICAL
        !                        COORDINATES AND ON THE SURFACE OF THE UNIT
        !                        SPHERE (RADIUS OF 1).  THE EQUATION IS
        !
        !                          (1/SIN(THETA))(D/DTHETA)(SIN(THETA)
        !                          (DU/DTHETA)) + (1/SIN(THETA)**2)(D/DPHI)
        !                          (DU/DPHI)  + LAMBDA*U = F(THETA, PHI)
        !
        !                        WHERE THETA IS COLATITUDE AND PHI IS
        !                        LONGITUDE.
        !
        ! USAGE                  CALL HWSSSP (TS, TF, M, MBDCND, BDTS, BDTF, PS, PF,
        !                                     N, NBDCND, BDPS, BDPF, ELMBDA, F,
        !                                     IDIMF, PERTRB, IERROR, W)
        !
        ! ARGUMENTS
        ! ON INPUT               TS, TF
        !
        !                          THE RANGE OF THETA (COLATITUDE), I.E.,
        !                          TS .LE. THETA .LE. TF. TS MUST BE LESS
        !                          THAN TF.  TS AND TF ARE IN RADIANS.
        !                          A TS OF ZERO CORRESPONDS TO THE NORTH
        !                          POLE AND A TF OF PI CORRESPONDS TO
        !                          THE SOUTH POLE.
        !
        !                          * * * IMPORTANT * * *
        !
        !                          IF TF IS EQUAL TO PI THEN IT MUST BE
        !                          COMPUTED USING THE STATEMENT
        !                          TF = PI_MACH(DUM). THIS INSURES THAT TF
        !                          IN THE USER'S PROGRAM IS EQUAL TO PI IN
        !                          THIS PROGRAM WHICH PERMITS SEVERAL TESTS
        !                          OF THE INPUT PARAMETERS THAT OTHERWISE
        !                          WOULD NOT BE POSSIBLE.
        !
        !
        !                        M
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL (TS, TF) IS SUBDIVIDED.
        !                          HENCE, THERE WILL BE M+1 GRID POINTS IN THE
        !                          THETA-DIRECTION GIVEN BY
        !                          THETA(I) = (I-1)DTHETA+TS FOR
        !                          I = 1, 2, ..., M+1, WHERE
        !                          DTHETA = (TF-TS)/M IS THE PANEL WIDTH.
        !                          M MUST BE GREATER THAN 5
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITION
        !                          AT THETA = TS AND THETA = TF.
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = TS AND THETA = TF.
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = TS AND THE DERIVATIVE OF
        !                               THE SOLUTION WITH RESPECT TO THETA IS
        !                               SPECIFIED AT THETA = TF
        !                               (SEE NOTE 2 BELOW).
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               SPECIFIED AT THETA = TS AND
        !                               THETA = TF (SEE NOTES 1, 2 BELOW).
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               AT THETA = TS (SEE NOTE 1 BELOW)
        !                               AND THE SOLUTION IS SPECIFIED AT
        !                               THETA = TF.
        !                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = TS = 0 AND THE SOLUTION
        !                               IS SPECIFIED AT THETA = TF.
        !                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = TS = 0 AND THE DERIVATIVE
        !                               OF THE SOLUTION WITH RESPECT TO THETA
        !                               IS SPECIFIED AT THETA = TF
        !                               (SEE NOTE 2 BELOW).
        !                          = 7  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = TS AND THE SOLUTION IS
        !                               IS UNSPECIFIED AT THETA = TF = PI.
        !                          = 8  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               AT THETA = TS (SEE NOTE 1 BELOW) AND
        !                               THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = TF = PI.
        !                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = TS = 0 AND THETA = TF = PI.
        !
        !                          NOTES:
        !                          IF TS = 0, DO NOT USE MBDCND = 3, 4, OR 8,
        !                          BUT INSTEAD USE MBDCND = 5, 6, OR 9  .
        !
        !                          IF TF = PI, DO NOT USE MBDCND = 2, 3, OR 6,
        !                          BUT INSTEAD USE MBDCND = 7, 8, OR 9  .
        !
        !                        BDTS
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
        !                          THE SOLUTION WITH RESPECT TO THETA AT
        !                          THETA = TS.  WHEN MBDCND = 3, 4, OR 8,
        !
        !                          BDTS(J) = (D/DTHETA)U(TS, PHI(J)),
        !                          J = 1, 2, ..., N+1  .
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDTS IS
        !                          A DUMMY VARIABLE.
        !
        !                        BDTF
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
        !                          THAT SPECIFIES THE VALUES OF THE DERIVATIVE
        !                          OF THE SOLUTION WITH RESPECT TO THETA AT
        !                          THETA = TF.  WHEN MBDCND = 2, 3, OR 6,
        !
        !                          BDTF(J) = (D/DTHETA)U(TF, PHI(J)),
        !                          J = 1, 2, ..., N+1  .
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDTF IS
        !                          A DUMMY VARIABLE.
        !
        !                        PS, PF
        !                          THE RANGE OF PHI (LONGITUDE), I.E.,
        !                          PS .LE. PHI .LE. PF.  PS MUST BE LESS
        !                          THAN PF.  PS AND PF ARE IN RADIANS.
        !                          IF PS = 0 AND PF = 2*PI, PERIODIC
        !                          BOUNDARY CONDITIONS ARE USUALLY PRESCRIBED.
        !
        !                          * * * IMPORTANT * * *
        !
        !                          IF PF IS EQUAL TO 2*PI THEN IT MUST BE
        !                          COMPUTED USING THE STATEMENT
        !                          PF = 2.*PI_MACH(DUM). THIS INSURES THAT
        !                          PF IN THE USERS PROGRAM IS EQUAL TO
        !                          2*PI IN THIS PROGRAM WHICH PERMITS TESTS
        !                          OF THE INPUT PARAMETERS THAT OTHERWISE
        !                          WOULD NOT BE POSSIBLE.
        !
        !                        N
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL (PS, PF) IS SUBDIVIDED.
        !                          HENCE, THERE WILL BE N+1 GRID POINTS
        !                          IN THE PHI-DIRECTION GIVEN BY
        !                          PHI(J) = (J-1)DPHI+PS  FOR
        !                          J = 1, 2, ..., N+1, WHERE
        !                          DPHI = (PF-PS)/N IS THE PANEL WIDTH.
        !                          N MUST BE GREATER THAN 4
        !
        !                        NBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITION
        !                          AT PHI = PS AND PHI = PF.
        !
        !                          = 0  IF THE SOLUTION IS PERIODIC IN PHI,
        !                               I.U., U(I, J) = U(I, N+J).
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               PHI = PS AND PHI = PF
        !                               (SEE NOTE BELOW).
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               PHI = PS (SEE NOTE BELOW)
        !                               AND THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO PHI IS SPECIFIED
        !                               AT PHI = PF.
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO PHI IS SPECIFIED
        !                               AT PHI = PS AND PHI = PF.
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO PHI IS SPECIFIED
        !                               AT PS AND THE SOLUTION IS SPECIFIED
        !                               AT PHI = PF
        !
        !                          NOTE:
        !                          NBDCND = 1, 2, OR 4 CANNOT BE USED WITH
        !                          MBDCND = 5, 6, 7, 8, OR 9.  THE FORMER INDICATES
        !                          THAT THE SOLUTION IS SPECIFIED AT A POLE, THE
        !                          LATTER INDICATES THAT THE SOLUTION IS NOT
        !                          SPECIFIED.  USE INSTEAD  MBDCND = 1 OR 2.
        !
        !                        BDPS
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE
        !                          OF THE SOLUTION WITH RESPECT TO PHI AT
        !                          PHI = PS.  WHEN NBDCND = 3 OR 4,
        !
        !                            BDPS(I) = (D/DPHI)U(THETA(I), PS),
        !                            I = 1, 2, ..., M+1  .
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDPS IS
        !                          A DUMMY VARIABLE.
        !
        !                        BDPF
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE
        !                          OF THE SOLUTION WITH RESPECT TO PHI AT
        !                          PHI = PF.  WHEN NBDCND = 2 OR 3,
        !
        !                            BDPF(I) = (D/DPHI)U(THETA(I), PF),
        !                            I = 1, 2, ..., M+1  .
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDPF IS
        !                          A DUMMY VARIABLE.
        !
        !                        ELMBDA
        !                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
        !                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
        !                          MAY NOT EXIST.  HOWEVER, HWSSSP WILL
        !                          ATTEMPT TO FIND A SOLUTION.
        !
        !                        F
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUE OF THE RIGHT SIDE OF THE HELMHOLTZ
        !                          EQUATION AND BOUNDARY VALUES (IF ANY).
        !                          F MUST BE DIMENSIONED AT LEAST (M+1)*(N+1).
        !
        !                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
        !                          FOR I = 2, 3, ..., M AND J = 2, 3, ..., N
        !                          F(I, J) = F(THETA(I), PHI(J)).
        !
        !                          ON THE BOUNDARIES F IS DEFINED AS FOLLOWS:
        !                          FOR J = 1, 2, ..., N+1 AND I = 1, 2, ..., M+1
        !
        !                          MBDCND   F(1, J)            F(M+1, J)
        !                          ------   ------------      ------------
        !
        !                            1      U(TS, PHI(J))      U(TF, PHI(J))
        !                            2      U(TS, PHI(J))      F(TF, PHI(J))
        !                            3      F(TS, PHI(J))      F(TF, PHI(J))
        !                            4      F(TS, PHI(J))      U(TF, PHI(J))
        !                            5      F(0, PS)           U(TF, PHI(J))
        !                            6      F(0, PS)           F(TF, PHI(J))
        !                            7      U(TS, PHI(J))      F(PI, PS)
        !                            8      F(TS, PHI(J))      F(PI, PS)
        !                            9      F(0, PS)           F(PI, PS)
        !
        !                          NBDCND   F(I, 1)            F(I, N+1)
        !                          ------   --------------    --------------
        !
        !                            0      F(THETA(I), PS)    F(THETA(I), PS)
        !                            1      U(THETA(I), PS)    U(THETA(I), PF)
        !                            2      U(THETA(I), PS)    F(THETA(I), PF)
        !                            3      F(THETA(I), PS)    F(THETA(I), PF)
        !                            4      F(THETA(I), PS)    U(THETA(I), PF)
        !
        !                          NOTE:
        !                          IF THE TABLE CALLS FOR BOTH THE SOLUTION U
        !                          AND THE RIGHT SIDE F AT A CORNER THEN THE
        !                          SOLUTION MUST BE SPECIFIED.
        !
        !                        IDIMF
        !                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
        !                          F AS IT APPEARS IN THE PROGRAM CALLING
        !                          HWSSSP.  THIS PARAMETER IS USED TO SPECIFY
        !                          THE VARIABLE DIMENSION OF F. IDIMF MUST BE
        !                          AT LEAST M+1  .
        !
        !
        ! ON OUTPUT              F
        !                          CONTAINS THE SOLUTION U(I, J) OF THE FINITE
        !                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
        !                          (THETA(I), PHI(J)),  I = 1, 2, ..., M+1  AND
        !                          J = 1, 2, ..., N+1  .
        !
        !                        PERTRB
        !                          IF ONE SPECIFIES A COMBINATION OF PERIODIC,
        !                          DERIVATIVE OR UNSPECIFIED BOUNDARY
        !                          CONDITIONS FOR A POISSON EQUATION
        !                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
        !                          PERTRB IS A CONSTANT, CALCULATED AND
        !                          SUBTRACTED FROM F, WHICH ENSURES THAT A
        !                          SOLUTION EXISTS.  HWSSSP THEN COMPUTES
        !                          THIS SOLUTION, WHICH IS A LEAST SQUARES
        !                          SOLUTION TO THE ORIGINAL APPROXIMATION.
        !                          THIS SOLUTION IS NOT UNIQUE AND IS
        !                          UNNORMALIZED. THE VALUE OF PERTRB SHOULD
        !                          BE SMALL COMPARED TO THE RIGHT SIDE F.
        !                          OTHERWISE , A SOLUTION IS OBTAINED TO AN
        !                          ESSENTIALLY DIFFERENT PROBLEM. THIS
        !                          COMPARISON SHOULD ALWAYS BE MADE TO INSURE
        !                          THAT A MEANINGFUL SOLUTION HAS BEEN
        !                          OBTAINED
        !
        !                        IERROR
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 8,
        !                          A SOLUTION IS NOT ATTEMPTED.
        !
        !                          = 0  NO ERROR
        !                          = 1  TS.LT.0 OR TF.GT.PI
        !                          = 2  TS.GE.TF
        !                          = 3  MBDCND.LT.1 OR MBDCND.GT.9
        !                          = 4  PS.LT.0 OR PS.GT.PI+PI
        !                          = 5  PS.GE.PF
        !                          = 6  N.LT.5
        !                          = 7  M.LT.5
        !                          = 8  NBDCND.LT.0 OR NBDCND.GT.4
        !                          = 9  ELMBDA.GT.0
        !                          = 10 IDIMF.LT.M+1
        !                          = 11 NBDCND EQUALS 1, 2 OR 4 AND MBDCND.GE.5
        !                          = 12 TS.EQ.0 AND MBDCND EQUALS 3, 4 OR 8
        !                          = 13 TF.EQ.PI AND MBDCND EQUALS 2, 3 OR 6
        !                          = 14 MBDCND EQUALS 5, 6 OR 9 AND TS.NE.0
        !                          = 15 MBDCND.GE.7 AND TF.NE.PI
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
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
        !                        LIBRARIES IN JANUARY 1980.
        !                        Revised in June 2004 by John Adams using
        !                        Fortran 90 dynamically allocated work space.
        !
        ! PORTABILITY            FORTRAN 90
        !
        ! ALGORITHM              THE ROUTINE DEFINES THE FINITE DIFFERENCE
        !                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
        !                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS
        !                        AND THEN CALLS GENBUN TO SOLVE THE SYSTEM.
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
        !                        SUBROUTINE GENBUN WHICH IS THE ROUTINE THAT
        !                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
        !
        ! REFERENCES             P. N. SWARZTRAUBER, "THE DIRECT SOLUTION OF
        !                        THE DISCRETE POISSON EQUATION ON THE SURFACE OF
        !                        A SPHERE", S.I.A.M. J. NUMER. ANAL., 15(1974), 
        !                        PP 212-215.
        !
        !                        SWARZTRAUBER, P. AND R. SWEET, "EFFICIENT
        !                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
        !                        ELLIPTIC EQUATIONS", NCAR TN/IA-109, JULY, 
        !                        1975, 138 PP.
        !***********************************************************************
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer  :: m
        integer  :: mbdcnd
        integer  :: n
        integer  :: nbdcnd
        integer  :: idimf
        integer  :: ierror
        real  :: ts
        real  :: tf
        real  :: ps
        real  :: pf
        real  :: elmbda
        real  :: pertrb
        real  :: bdts(*)
        real  :: bdtf(*)
        real  :: bdps(*)
        real  :: bdpf(*)
        real  :: f(idimf, 1)
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: nbr, irwk, icwk
        real :: pi, dum, tpi
        !-----------------------------------------------
        !
        nbr = nbdcnd + 1
        pi = acos( -1.0 )
        tpi = 2.*pi
        ierror = 0
        if (ts<0. .or. tf>pi) ierror = 1
        if (ts >= tf) ierror = 2
        if (mbdcnd<1 .or. mbdcnd>9) ierror = 3
        if (ps<0. .or. pf>tpi) ierror = 4
        if (ps >= pf) ierror = 5
        if (n < 5) ierror = 6
        if (m < 5) ierror = 7
        if (nbdcnd<0 .or. nbdcnd>4) ierror = 8
        if (elmbda > 0.) ierror = 9
        if (idimf < m + 1) ierror = 10
        if ((nbdcnd==1 .or. nbdcnd==2 .or. nbdcnd==4) .and. mbdcnd>=5) ierror = 11
        if(ts==0..and.(mbdcnd==3.or.mbdcnd==4.or.mbdcnd==8))ierror=12
        if(tf==pi.and.(mbdcnd==2.or.mbdcnd==3.or.mbdcnd==6))ierror=13
        if((mbdcnd==5.or.mbdcnd==6.or.mbdcnd==9).and.ts/=0.)ierror=14
        if (mbdcnd>=7 .and. tf/=pi) ierror = 15
        if (ierror/=0 .and. ierror/=9) return

        !     allocate generous work space estimate
        irwk=4*(n+1)+(16+INT(log(real(n+1))/log(2.0)))*(m+1)
        icwk = 0

        call workspace%create( irwk, icwk, ierror )
        !     check that allocation was successful
        if (ierror == 20) return

        call hwssspp(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, bdps, &
            bdpf, elmbda, f, idimf, pertrb, ierror, workspace%rew)
        !     release dynamically allocated work space
        call workspace%destroy()

    end subroutine HWSSSP

    subroutine HWSSSPP(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, &
        nbdcnd, bdps, bdpf, elmbda, f, idimf, pertrb, ierror, w)
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer  :: m
        integer  :: mbdcnd
        integer  :: n
        integer  :: nbdcnd
        integer  :: idimf
        integer  :: ierror
        real  :: ts
        real  :: tf
        real  :: ps
        real  :: pf
        real  :: elmbda
        real  :: pertrb
        real  :: bdts(*)
        real  :: bdtf(*)
        real  :: bdps(*)
        real  :: bdpf(*)
        real  :: f(idimf, *)
        real  :: w(*)
        !-----------------------------------------------
        call HWSSS1 (ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
            bdps, bdpf, elmbda, f, idimf, pertrb, w, W(m+2), W(2*m+3), W(3* &
            m+4), W(4*m+5), W(5*m+6), W(6*m+7))

    end subroutine HWSSSPP



    subroutine HWSSS1(TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, NBDCND, &
        BDPS, BDPF, ELMBDA, F, IDIMF, PERTRB, AM, BM, CM, SN, SS, &
        SINT, D)
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: M
        integer , intent (in) :: MBDCND
        integer , intent (in) :: N
        integer  :: NBDCND
        integer  :: IDIMF
        real , intent (in) :: TS
        real , intent (in) :: TF
        real , intent (in) :: PS
        real , intent (in) :: PF
        real , intent (in) :: ELMBDA
        real , intent (out) :: PERTRB
        real , intent (in) :: BDTS(*)
        real , intent (in) :: BDTF(*)
        real , intent (in) :: BDPS(*)
        real , intent (in) :: BDPF(*)
        real  :: F(IDIMF,*)
        real  :: AM(*)
        real  :: BM(*)
        real  :: CM(*)
        real , intent (inout) :: SN(*)
        real , intent (inout) :: SS(*)
        real , intent (inout) :: SINT(*)
        real  :: D(*)
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: MP1, NP1, I, INP, ISP, MBR, ITS, ITF, ITSP, ITFM, MUNK &
            , IID, II, NBR, JPS, JPF, JPSP, JPFM, NUNK, ISING, J, IERROR
        real :: PI, DUM, TPI, HPI, FN, FM, DTH, HDTH, TDT, DPHI, TDP, &
            DPHI2, EDP2, DTH2, CP, WP, FIM1, THETA, T1, AT, CT, WTS, WTF, &
            WPS, WPF, FJJ, CF, SUM, SUM1, HNE, YHLD, SUM2, DFN, DNN, DSN, &
            CNP, HLD, DFS, DSS, DNS, CSP, RTN, RTS, DEN
        !-----------------------------------------------
        !
        PI = 4.0*ATAN(1.0)
        TPI = PI + PI
        HPI = PI/2.
        MP1 = M + 1
        NP1 = N + 1
        FN = N
        FM = M
        DTH = (TF - TS)/FM
        HDTH = DTH/2.
        TDT = DTH + DTH
        DPHI = (PF - PS)/FN
        TDP = DPHI + DPHI
        DPHI2 = DPHI*DPHI
        EDP2 = ELMBDA*DPHI2
        DTH2 = DTH*DTH
        CP = 4./(FN*DTH2)
        WP = FN*SIN(HDTH)/4.
        do I = 1, MP1
            FIM1 = I - 1
            THETA = FIM1*DTH + TS
            SINT(I) = SIN(THETA)
            if (SINT(I) == 0.) cycle
            T1 = 1./(DTH2*SINT(I))
            AM(I) = T1*SIN(THETA - HDTH)
            CM(I) = T1*SIN(THETA + HDTH)
            BM(I) = (-AM(I)) - CM(I) + ELMBDA
        end do
        INP = 0
        ISP = 0
        !
        ! BOUNDARY CONDITION AT THETA=TS
        !
        MBR = MBDCND + 1
        go to (103,104,104,105,105,106,106,104,105,106) MBR
103 continue
    ITS = 1
    go to 107
104 continue
    AT = AM(2)
    ITS = 2
    go to 107
105 continue
    AT = AM(1)
    ITS = 1
    CM(1) = AM(1) + CM(1)
    go to 107
106 continue
    AT = AM(2)
    INP = 1
    ITS = 2
!
! BOUNDARY CONDITION THETA=TF
!
107 continue
    go to (108,109,110,110,109,109,110,111,111,111) MBR
108 continue
    ITF = M
    go to 112
109 continue
    CT = CM(M)
    ITF = M
    go to 112
110 continue
    CT = CM(M+1)
    AM(M+1) = AM(M+1) + CM(M+1)
    ITF = M + 1
    go to 112
111 continue
    ITF = M
    ISP = 1
    CT = CM(M)
!
! COMPUTE HOMOGENEOUS SOLUTION WITH SOLUTION AT POLE EQUAL TO ONE
!
112 continue
    ITSP = ITS + 1
    ITFM = ITF - 1
    WTS = SINT(ITS+1)*AM(ITS+1)/CM(ITS)
    WTF = SINT(ITF-1)*CM(ITF-1)/AM(ITF)
    MUNK = ITF - ITS + 1
    if (ISP > 0) then
        D(ITS) = CM(ITS)/BM(ITS)
        do I = ITSP, M
            D(I) = CM(I)/(BM(I)-AM(I)*D(I-1))
        end do
        SS(M) = -D(M)
        IID = M - ITS
        do II = 1, IID
            I = M - II
            SS(I) = -D(I)*SS(I+1)
        end do
        SS(M+1) = 1.
    endif
    if (INP > 0) then
        SN(1) = 1.
        D(ITF) = AM(ITF)/BM(ITF)
        IID = ITF - 2
        do II = 1, IID
            I = ITF - II
            D(I) = AM(I)/(BM(I)-CM(I)*D(I+1))
        end do
        SN(2) = -D(2)
        do I = 3, ITF
            SN(I) = -D(I)*SN(I-1)
        end do
    endif
    !
    ! BOUNDARY CONDITIONS AT PHI=PS
    !
    NBR = NBDCND + 1
    WPS = 1.
    WPF = 1.
    select case (NBR)
        case default
            JPS = 1
        case (2:3)
            JPS = 2
        case (4:5)
            JPS = 1
            WPS = 0.5
    end select
!
! BOUNDARY CONDITION AT PHI=PF
!
124 continue
    go to (125,126,127,127,126) NBR
125 continue
    JPF = N
    go to 128
126 continue
    JPF = N
    go to 128
127 continue
    WPF = 0.5
    JPF = N + 1
128 continue
    JPSP = JPS + 1
    JPFM = JPF - 1
    NUNK = JPF - JPS + 1
    FJJ = JPFM - JPSP + 1
    !
    ! SCALE COEFFICIENTS FOR SUBROUTINE GENBUN
    !
    do I = ITS, ITF
        CF = DPHI2*SINT(I)*SINT(I)
        AM(I) = CF*AM(I)
        BM(I) = CF*BM(I)
        CM(I) = CF*CM(I)
    end do
    AM(ITS) = 0.
    CM(ITF) = 0.
    ISING = 0
    go to (130,138,138,130,138,138,130,138,130,130) MBR
130 continue
    go to (131,138,138,131,138) NBR
131 continue
    if (ELMBDA >= 0.) then
        ISING = 1
        SUM = WTS*WPS + WTS*WPF + WTF*WPS + WTF*WPF
        if (INP > 0) then
            SUM = SUM + WP
        endif
        if (ISP > 0) then
            SUM = SUM + WP
        endif
        SUM1 = 0.
        do I = ITSP, ITFM
            SUM1 = SUM1 + SINT(I)
        end do
        SUM = SUM + FJJ*(SUM1 + WTS + WTF)
        SUM = SUM + (WPS + WPF)*SUM1
        HNE = SUM
    endif
138 continue
    go to (146,142,142,144,144,139,139,142,144,139) MBR
139 continue
    if (NBDCND - 3 /= 0) go to 146
    YHLD = F(1,JPS) - 4./(FN*DPHI*DTH2)*(BDPF(2)-BDPS(2))
    F(1,:NP1) = YHLD
    go to 146
142 continue
    F(2,JPS:JPF) = F(2,JPS:JPF) - AT*F(1,JPS:JPF)
    go to 146
144 continue
    F(1,JPS:JPF) = F(1,JPS:JPF) + TDT*BDTS(JPS:JPF)*AT
146 continue
    go to (154,150,152,152,150,150,152,147,147,147) MBR
147 continue
    if (NBDCND - 3 /= 0) go to 154
    YHLD = F(M+1,JPS) - 4./(FN*DPHI*DTH2)*(BDPF(M)-BDPS(M))
    F(M+1,:NP1) = YHLD
    go to 154
150 continue
    F(M,JPS:JPF) = F(M,JPS:JPF) - CT*F(M+1,JPS:JPF)
    go to 154
152 continue
    F(M+1,JPS:JPF) = F(M+1,JPS:JPF) - TDT*BDTF(JPS:JPF)*CT
154 continue
    go to (159,155,155,157,157) NBR
155 continue
    F(ITS:ITF,2) = F(ITS:ITF,2) - F(ITS:ITF,1)/(DPHI2*SINT(ITS:ITF)* &
        SINT(ITS:ITF))
    go to 159
157 continue
    F(ITS:ITF,1) = F(ITS:ITF,1) + TDP*BDPS(ITS:ITF)/(DPHI2*SINT(ITS: &
        ITF)*SINT(ITS:ITF))
159 continue
    go to (164,160,162,162,160) NBR
160 continue
    F(ITS:ITF,N) = F(ITS:ITF,N) - F(ITS:ITF,N+1)/(DPHI2*SINT(ITS:ITF)* &
        SINT(ITS:ITF))
    go to 164
162 continue
    F(ITS:ITF,N+1) = F(ITS:ITF,N+1) - TDP*BDPF(ITS:ITF)/(DPHI2*SINT( &
        ITS:ITF)*SINT(ITS:ITF))
164 continue
    PERTRB = 0.
    if (ISING /= 0) then
        SUM = WTS*WPS*F(ITS,JPS) + WTS*WPF*F(ITS,JPF) + WTF*WPS*F(ITF, &
            JPS) + WTF*WPF*F(ITF,JPF)
        if (INP > 0) then
            SUM = SUM + WP*F(1,JPS)
        endif
        if (ISP > 0) then
            SUM = SUM + WP*F(M+1,JPS)
        endif
        do I = ITSP, ITFM
            SUM1 = 0.
            do J = JPSP, JPFM
                SUM1 = SUM1 + F(I,J)
            end do
            SUM = SUM + SINT(I)*SUM1
        end do
        SUM1 = 0.
        SUM2 = 0.
        do J = JPSP, JPFM
            SUM1 = SUM1 + F(ITS,J)
            SUM2 = SUM2 + F(ITF,J)
        end do
        SUM = SUM + WTS*SUM1 + WTF*SUM2
        SUM1 = 0.
        SUM2 = 0.
        SUM1 = DOT_PRODUCT(SINT(ITSP:ITFM),F(ITSP:ITFM,JPS))
        SUM2 = DOT_PRODUCT(SINT(ITSP:ITFM),F(ITSP:ITFM,JPF))
        SUM = SUM + WPS*SUM1 + WPF*SUM2
        PERTRB = SUM/HNE
        F(:MP1,:NP1) = F(:MP1,:NP1) - PERTRB
    endif
    !
    ! SCALE RIGHT SIDE FOR SUBROUTINE GENBUN
    !
    do I = ITS, ITF
        CF = DPHI2*SINT(I)*SINT(I)
        F(I,JPS:JPF) = CF*F(I,JPS:JPF)
    end do
    call GENBUNN (NBDCND, NUNK, 1, MUNK, AM(ITS), BM(ITS), CM(ITS), &
        IDIMF, F(ITS,JPS), IERROR, D)
    if (ISING <= 0) go to 186
    if (INP > 0) then
        if (ISP > 0) go to 186
        F(1,:NP1) = 0.
        go to 209
    endif
    if (ISP <= 0) go to 186
    F(M+1,:NP1) = 0.
    go to 209
186 continue
    if (INP > 0) then
        SUM = WPS*F(ITS,JPS) + WPF*F(ITS,JPF)
        do J = JPSP, JPFM
            SUM = SUM + F(ITS,J)
        end do
        DFN = CP*SUM
        DNN = CP*((WPS + WPF + FJJ)*(SN(2)-1.)) + ELMBDA
        DSN = CP*(WPS + WPF + FJJ)*SN(M)
        if (ISP > 0) go to 194
        CNP = (F(1,1)-DFN)/DNN
        do I = ITS, ITF
            HLD = CNP*SN(I)
            F(I,JPS:JPF) = F(I,JPS:JPF) + HLD
        end do
        F(1,:NP1) = CNP
        go to 209
    endif
    if (ISP <= 0) go to 209
194 continue
    SUM = WPS*F(ITF,JPS) + WPF*F(ITF,JPF)
    do J = JPSP, JPFM
        SUM = SUM + F(ITF,J)
    end do
    DFS = CP*SUM
    DSS = CP*((WPS + WPF + FJJ)*(SS(M)-1.)) + ELMBDA
    DNS = CP*(WPS + WPF + FJJ)*SS(2)
    if (INP <= 0) then
        CSP = (F(M+1,1)-DFS)/DSS
        do I = ITS, ITF
            HLD = CSP*SS(I)
            F(I,JPS:JPF) = F(I,JPS:JPF) + HLD
        end do
        F(M+1,:NP1) = CSP
    else
        RTN = F(1,1) - DFN
        RTS = F(M+1,1) - DFS
        if (ISING > 0) then
            CSP = 0.
            CNP = RTN/DNN
        else
            if (ABS(DNN) - ABS(DSN) > 0.) then
                DEN = DSS - DNS*DSN/DNN
                RTS = RTS - RTN*DSN/DNN
                CSP = RTS/DEN
                CNP = (RTN - CSP*DNS)/DNN
            else
                DEN = DNS - DSS*DNN/DSN
                RTN = RTN - RTS*DNN/DSN
                CSP = RTN/DEN
                CNP = (RTS - DSS*CSP)/DSN
            endif
        endif
        do I = ITS, ITF
            HLD = CNP*SN(I) + CSP*SS(I)
            F(I,JPS:JPF) = F(I,JPS:JPF) + HLD
        end do
        F(1,:NP1) = CNP
        F(M+1,:NP1) = CSP
    endif
209 continue
    if (NBDCND == 0) then
        F(:MP1,JPF+1) = F(:MP1,JPS)
    endif

end subroutine HWSSS1

end module module_hwsssp
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
