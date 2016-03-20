module module_hwsssp

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_genbun, only: &
        genbunn

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hwsssp
    public :: test_hwsssp

contains

    subroutine test_hwsssp()
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
        !     PROGRAM TO ILLUSTRATE THE USE OF hwsssp
        !
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer :: m, mbdcnd, n, nbdcnd, idimf, mp1, i, np1, j, ierror
        real (wp), dimension(19, 73) :: f
        real (wp), dimension(73) :: bdtf, bdts, bdps, bdpf
        real (wp), dimension(19) :: sint
        real (wp), dimension(73) :: sinp
        real :: pi, ts, tf, ps, pf, elmbda, dtheta, dphi, pertrb, discretization_error, z
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
        !     generate sines for use in subsequent computations
        !
        dtheta = tf/real(m)
        mp1 = m + 1
        do i = 1, mp1
            sint(i) = sin(real(i - 1)*dtheta)
        end do
        dphi = (pi + pi)/real(n)
        np1 = n + 1
        do j = 1, np1
            sinp(j) = sin(real(j - 1)*dphi)
        end do
        !
        !     compute right side of equation and store in f
        !
        do j = 1, np1
            f(:mp1, j) = 2. - 6.*(sint(:mp1)*sinp(j))**2
        end do
        !
        !     store derivative data at the equator
        !
        bdtf(:np1) = 0.
        !
        call hwsssp(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
            bdps, bdpf, elmbda, f, idimf, pertrb, ierror)
        !
        !     compute discretization error. since problem is singular, the
        !     solution must be normalized.
        !
        discretization_error = 0.0
        do j = 1, np1
            do i = 1, mp1
                z = abs(f(i, j)-(sint(i)*sinp(j))**2-f(1, 1))
                discretization_error = max(z, discretization_error)
            end do
        end do

        write( stdout, '(A)') ''
        write( stdout, '(A)') '     hwsssp *** TEST RUN *** '
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') '     ierror = 0,  discretization error = 3.38107E-3'
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)') &
            '      ierror =', ierror, ' discretization error = ', discretization_error

    end subroutine test_hwsssp


    subroutine hwsssp(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
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
        !     SUBROUTINE hwsssp (TS, TF, M, MBDCND, BDTS, BDTF, PS, PF, N, NBDCND, BDPS,
        !    +                   BDPF, ELMBDA, F, IDIMF, PERTRB, ierror)
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
        ! USAGE                  CALL hwsssp (TS, TF, M, MBDCND, BDTS, BDTF, PS, PF,
        !                                     N, NBDCND, BDPS, BDPF, ELMBDA, F,
        !                                     IDIMF, PERTRB, ierror, W)
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
        !                          MAY NOT EXIST.  HOWEVER, hwsssp WILL
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
        !                          hwsssp.  THIS PARAMETER IS USED TO SPECIFY
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
        !                          SOLUTION EXISTS.  hwsssp THEN COMPUTES
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
        !                        ierror
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
        ! REFERENCES             P. N. SWARZTRAUBER, "THE DIRECT SOLUTION OF
        !                        THE DISCRETE POISSON EQUATION ON THE SURFACE OF
        !                        A SPHERE", S.I.A.M. J. NUMER. ANAL., 15(1974), 
        !                        PP 212-215.
        !
        !                        SWARZTRAUBER, P. AND R. SWEET, "EFFICIENT
        !                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
        !                        ELLIPTIC EQUATIONS", NCAR TN/IA-109, JULY, 
        !                        1975, 138 PP.
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip),          intent (in)     :: m
        integer (ip),          intent (in)     :: mbdcnd
        integer (ip),          intent (in)     :: n
        integer (ip),          intent (in)     :: nbdcnd
        integer (ip),          intent (in)     :: idimf
        integer (ip),          intent (out)    :: ierror
        real (wp),             intent (in)     :: ts
        real (wp),             intent (in)     :: tf
        real (wp),             intent (in)     :: ps
        real (wp),             intent (in)     :: pf
        real (wp),             intent (in)     :: elmbda
        real (wp),             intent (out)    :: pertrb
        real (wp), contiguous, intent (in)     :: bdts(:)
        real (wp), contiguous, intent (in)     :: bdtf(:)
        real (wp), contiguous, intent (in)     :: bdps(:)
        real (wp), contiguous, intent (in)     :: bdpf(:)
        real (wp), contiguous, intent (in out) :: f(:,:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------

        ! initialize error flag
        ierror = 0

        ! Check if input values are valid
        associate( pi => acos( -1.0_wp ) )
            associate( two_pi => 2.0_wp * pi)
                if (ts<0. .or. tf>pi) ierror = 1
                if (ts >= tf) ierror = 2
                if (mbdcnd<1 .or. mbdcnd>9) ierror = 3
                if (ps<0. .or. pf>two_pi) ierror = 4
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
            end associate
        end associate

        ! allocate generous work space estimate
        associate( &
            irwk => 4 * (n + 1) + (16 + int(log(real(n+1,kind=wp))/log(2.0_wp), kind=ip)) * (m + 1), &
            icwk => 0 &
            )
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! check that allocation was successful
        if (ierror == 20) return

        ! solve system
        associate( rew => workspace%real_workspace )
            call hwssspp(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, bdps, &
                bdpf, elmbda, f, idimf, pertrb, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hwsssp


    subroutine hwssspp(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, &
        nbdcnd, bdps, bdpf, elmbda, f, idimf, pertrb, ierror, w)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: mbdcnd
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: nbdcnd
        integer (ip), intent (in)     :: idimf
        integer (ip), intent (in)     :: ierror
        real (wp),    intent (in)     :: ts
        real (wp),    intent (in)     :: tf
        real (wp),    intent (in)     :: ps
        real (wp),    intent (in)     :: pf
        real (wp),    intent (in)     :: elmbda
        real (wp),    intent (out)    :: pertrb
        real (wp),    intent (in)     :: bdts(*)
        real (wp),    intent (in)     :: bdtf(*)
        real (wp),    intent (in)     :: bdps(*)
        real (wp),    intent (in)     :: bdpf(*)
        real (wp),    intent (in out) :: f(idimf,*)
        real (wp),    intent (in out) :: w(*)
        !-----------------------------------------------

        call hwsss1(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
            bdps, bdpf, elmbda, f, idimf, pertrb, w, w(m+2), w(2*m+3), &
            w(3*m+4), w(4*m+5), w(5*m+6), w(6*m+7))

    end subroutine hwssspp


    subroutine hwsss1(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
        bdps, bdpf, elmbda, f, idimf, pertrb, am, bm, cm, sn, ss, &
        sint, d)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: mbdcnd
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: nbdcnd
        integer (ip), intent (in)     :: idimf
        real (wp),    intent (in)     :: ts
        real (wp),    intent (in)     :: tf
        real (wp),    intent (in)     :: ps
        real (wp),    intent (in)     :: pf
        real (wp),    intent (in)     :: elmbda
        real (wp),    intent (out)    :: pertrb
        real (wp),    intent (in)     :: bdts(*)
        real (wp),    intent (in)     :: bdtf(*)
        real (wp),    intent (in)     :: bdps(*)
        real (wp),    intent (in)     :: bdpf(*)
        real (wp),    intent (in out) :: f(idimf,*)
        real (wp),    intent (in out) :: am(*)
        real (wp),    intent (in out) :: bm(*)
        real (wp),    intent (in out) :: cm(*)
        real (wp),    intent (in out) :: sn(*)
        real (wp),    intent (in out) :: ss(*)
        real (wp),    intent (in out) :: sint(*)
        real (wp),    intent (in out) :: d(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer :: mp1, np1, i, inp, isp, mbr, its, itf, itsp, itfm, munk &
            , iid, ii, nbr, jps, jpf, jpsp, jpfm, nunk, ising, j, ierror
        real :: pi, dum, tpi, hpi, fn, fm, dth, hdth, tdt, dphi, tdp, &
            dphi2, edp2, dth2, cp, wp_rename, fim1, theta, t1, at, ct, wts, wtf, &
            wps, wpf, fjj, cf, summation, sum1, hne, yhld, sum2, dfn, dnn, dsn, &
            cnp, hld, dfs, dss, dns, csp, rtn, rts, den
        !-----------------------------------------------
        !
        pi = acos( -1.0 )
        tpi = pi + pi
        hpi = pi/2.
        mp1 = m + 1
        np1 = n + 1
        fn = n
        fm = m
        dth = (tf - ts)/fm
        hdth = dth/2.
        tdt = dth + dth
        dphi = (pf - ps)/fn
        tdp = dphi + dphi
        dphi2 = dphi*dphi
        edp2 = elmbda*dphi2
        dth2 = dth*dth
        cp = 4./(fn*dth2)
        wp_rename = fn*SIN(hdth)/4.
        do i = 1, mp1
            fim1 = i - 1
            theta = fim1*dth + ts
            sint(i) = SIN(theta)
            if (SINT(i) == 0.) cycle
            t1 = 1./(dth2*SINT(i))
            am(i) = t1*SIN(theta - hdth)
            cm(i) = t1*SIN(theta + hdth)
            bm(i) = (-AM(i)) - CM(i) + elmbda
        end do
        inp = 0
        isp = 0
        !
        ! BOUNDARY CONDITION AT THETA=TS
        !
        mbr = mbdcnd + 1
        go to (103,104,104,105,105,106,106,104,105,106) mbr
103 continue
    its = 1
    go to 107
104 continue
    at = AM(2)
    its = 2
    go to 107
105 continue
    at = AM(1)
    its = 1
    cm(1) = AM(1) + CM(1)
    go to 107
106 continue
    at = AM(2)
    inp = 1
    its = 2
!
! BOUNDARY CONDITION THETA=TF
!
107 continue
    go to (108,109,110,110,109,109,110,111,111,111) mbr
108 continue
    itf = m
    go to 112
109 continue
    ct = CM(m)
    itf = m
    go to 112
110 continue
    ct = CM(m+1)
    am(m+1) = AM(m+1) + CM(m+1)
    itf = m + 1
    go to 112
111 continue
    itf = m
    isp = 1
    ct = CM(m)
!
! COMPUTE HOMOGENEOUS SOLUTION WITH SOLUTION AT POLE EQUAL TO ONE
!
112 continue
    itsp = its + 1
    itfm = itf - 1
    wts = SINT(its+1)*AM(its+1)/CM(its)
    wtf = SINT(itf-1)*CM(itf-1)/AM(itf)
    munk = itf - its + 1
    if (isp > 0) then
        d(its) = CM(its)/BM(its)
        do i = itsp, m
            d(i) = CM(i)/(BM(i)-AM(i)*D(i-1))
        end do
        ss(m) = -D(m)
        iid = m - its
        do ii = 1, iid
            i = m - ii
            ss(i) = -D(i)*SS(i+1)
        end do
        ss(m+1) = 1.
    end if
    if (inp > 0) then
        sn(1) = 1.
        d(itf) = AM(itf)/BM(itf)
        iid = itf - 2
        do ii = 1, iid
            i = itf - ii
            d(i) = AM(i)/(BM(i)-CM(i)*D(i+1))
        end do
        sn(2) = -D(2)
        do i = 3, itf
            sn(i) = -D(i)*SN(i-1)
        end do
    end if
    !
    ! BOUNDARY CONDITIONS AT PHI=PS
    !
    nbr = nbdcnd + 1
    wps = 1.
    wpf = 1.
    select case (nbr)
        case default
            jps = 1
        case (2:3)
            jps = 2
        case (4:5)
            jps = 1
            wps = 0.5
    end select
!
! BOUNDARY CONDITION AT PHI=PF
!
124 continue
    go to (125,126,127,127,126) nbr
125 continue
    jpf = n
    go to 128
126 continue
    jpf = n
    go to 128
127 continue
    wpf = 0.5
    jpf = n + 1
128 continue
    jpsp = jps + 1
    jpfm = jpf - 1
    nunk = jpf - jps + 1
    fjj = jpfm - jpsp + 1
    !
    ! SCALE COEFFICIENTS FOR SUBROUTINE genbun
    !
    do i = its, itf
        cf = dphi2*SINT(i)*SINT(i)
        am(i) = cf*AM(i)
        bm(i) = cf*BM(i)
        cm(i) = cf*CM(i)
    end do
    am(its) = 0.
    cm(itf) = 0.
    ising = 0
    go to (130,138,138,130,138,138,130,138,130,130) mbr
130 continue
    go to (131,138,138,131,138) nbr
131 continue
    if (elmbda >= 0.) then
        ising = 1
        summation = wts*wps + wts*wpf + wtf*wps + wtf*wpf
        if (inp > 0) then
            summation = summation + wp_rename
        end if
        if (isp > 0) then
            summation = summation + wp_rename
        end if
        sum1 = 0.
        do i = itsp, itfm
            sum1 = sum1 + SINT(i)
        end do
        summation = summation + fjj*(sum1 + wts + wtf)
        summation = summation + (wps + wpf)*sum1
        hne = summation
    end if
138 continue
    go to (146,142,142,144,144,139,139,142,144,139) mbr
139 continue
    if (nbdcnd - 3 /= 0) go to 146
    yhld = F(1,jps) - 4./(fn*dphi*dth2)*(BDPF(2)-BDPS(2))
    f(1,:np1) = yhld
    go to 146
142 continue
    f(2,jps:jpf) = F(2,jps:jpf) - at*F(1,jps:jpf)
    go to 146
144 continue
    f(1,jps:jpf) = F(1,jps:jpf) + tdt*BDTS(jps:jpf)*at
146 continue
    go to (154,150,152,152,150,150,152,147,147,147) mbr
147 continue
    if (nbdcnd - 3 /= 0) go to 154
    yhld = F(m+1,jps) - 4./(fn*dphi*dth2)*(BDPF(m)-BDPS(m))
    f(m+1,:np1) = yhld
    go to 154
150 continue
    f(m,jps:jpf) = F(m,jps:jpf) - ct*F(m+1,jps:jpf)
    go to 154
152 continue
    f(m+1,jps:jpf) = F(m+1,jps:jpf) - tdt*BDTF(jps:jpf)*ct
154 continue
    go to (159,155,155,157,157) nbr
155 continue
    f(its:itf,2) = F(its:itf,2) - F(its:itf,1)/(dphi2*SINT(its:itf)* &
        SINT(its:itf))
    go to 159
157 continue
    f(its:itf,1) = F(its:itf,1) + tdp*BDPS(its:itf)/(dphi2*SINT(its: &
        itf)*SINT(its:itf))
159 continue
    go to (164,160,162,162,160) nbr
160 continue
    f(its:itf,n) = F(its:itf,n) - F(its:itf,n+1)/(dphi2*SINT(its:itf)* &
        SINT(its:itf))
    go to 164
162 continue
    f(its:itf,n+1) = F(its:itf,n+1) - tdp*BDPF(its:itf)/(dphi2*SINT( &
        its:itf)*SINT(its:itf))
164 continue
    pertrb = 0.
    if (ising /= 0) then
        summation = wts*wps*F(its,jps) + wts*wpf*F(its,jpf) + wtf*wps*F(itf, &
            jps) + wtf*wpf*F(itf,jpf)
        if (inp > 0) then
            summation = summation + wp_rename*F(1,jps)
        end if
        if (isp > 0) then
            summation = summation + wp_rename*F(m+1,jps)
        end if
        do i = itsp, itfm
            sum1 = 0.
            do j = jpsp, jpfm
                sum1 = sum1 + F(i,j)
            end do
            summation = summation + SINT(i)*sum1
        end do
        sum1 = 0.
        sum2 = 0.
        do j = jpsp, jpfm
            sum1 = sum1 + F(its,j)
            sum2 = sum2 + F(itf,j)
        end do
        summation = summation + wts*sum1 + wtf*sum2
        sum1 = 0.
        sum2 = 0.
        sum1 = DOT_PRODUCT(SINT(itsp:itfm),F(itsp:itfm,jps))
        sum2 = DOT_PRODUCT(SINT(itsp:itfm),F(itsp:itfm,jpf))
        summation = summation + wps*sum1 + wpf*sum2
        pertrb = summation/hne
        f(:mp1,:np1) = F(:mp1,:np1) - pertrb
    end if
    !
    ! SCALE RIGHT SIDE FOR SUBROUTINE genbun
    !
    do i = its, itf
        cf = dphi2*SINT(i)*SINT(i)
        f(i,jps:jpf) = cf*F(i,jps:jpf)
    end do

    call genbunn(nbdcnd, nunk, 1, munk, AM(its), BM(its), CM(its), &
        idimf, F(its,jps), ierror, d)

    if (ising <= 0) go to 186
    if (inp > 0) then
        if (isp > 0) go to 186
        f(1,:np1) = 0.
        go to 209
    end if
    if (isp <= 0) go to 186
    f(m+1,:np1) = 0.
    go to 209
186 continue
    if (inp > 0) then
        summation = wps*F(its,jps) + wpf*F(its,jpf)
        do j = jpsp, jpfm
            summation = summation + F(its,j)
        end do
        dfn = cp*summation
        dnn = cp*((wps + wpf + fjj)*(SN(2)-1.)) + elmbda
        dsn = cp*(wps + wpf + fjj)*SN(m)
        if (isp > 0) go to 194
        cnp = (F(1,1)-dfn)/dnn
        do i = its, itf
            hld = cnp*SN(i)
            f(i,jps:jpf) = F(i,jps:jpf) + hld
        end do
        f(1,:np1) = cnp
        go to 209
    end if
    if (isp <= 0) go to 209
194 continue
    summation = wps*F(itf,jps) + wpf*F(itf,jpf)
    do j = jpsp, jpfm
        summation = summation + F(itf,j)
    end do
    dfs = cp*summation
    dss = cp*((wps + wpf + fjj)*(SS(m)-1.)) + elmbda
    dns = cp*(wps + wpf + fjj)*SS(2)
    if (inp <= 0) then
        csp = (F(m+1,1)-dfs)/dss
        do i = its, itf
            hld = csp*SS(i)
            f(i,jps:jpf) = F(i,jps:jpf) + hld
        end do
        f(m+1,:np1) = csp
    else
        rtn = F(1,1) - dfn
        rts = F(m+1,1) - dfs
        if (ising > 0) then
            csp = 0.
            cnp = rtn/dnn
        else
            if (ABS(dnn) - ABS(dsn) > 0.) then
                den = dss - dns*dsn/dnn
                rts = rts - rtn*dsn/dnn
                csp = rts/den
                cnp = (rtn - csp*dns)/dnn
            else
                den = dns - dss*dnn/dsn
                rtn = rtn - rts*dnn/dsn
                csp = rtn/den
                cnp = (rts - dss*csp)/dsn
            end if
        end if
        do i = its, itf
            hld = cnp*SN(i) + csp*SS(i)
            f(i,jps:jpf) = F(i,jps:jpf) + hld
        end do
        f(1,:np1) = cnp
        f(m+1,:np1) = csp
    end if
209 continue
    if (nbdcnd == 0) then
        f(:mp1,jpf+1) = F(:mp1,jps)
    end if

end subroutine hwsss1

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
