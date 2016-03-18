module module_hwsplr

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_genbun, only: &
        genbunn

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hwsplr
    public :: test_hwsplr

contains

    subroutine test_hwsplr()
        !
        !     file thwsplr.f
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
        ! Dictionary: local variables
        !-----------------------------------------------
        integer :: idimf, m, mbdcnd, n, nbdcnd, mp1, np1, i, j, ierror
        real , dimension(100, 50) :: f
        real , dimension(51) :: bdc, bdd, r, bda, bdb
        real , dimension(49) :: theta
        real :: a, b, c, pi, dum, d, elmbda, pertrb, err, z
        !-----------------------------------------------
        !
        !          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE hwsplr TO SOLVE
        !     THE EQUATION
        !
        !     (1/R)(D/DR)(R*(DU/DR)) + (1/R**2)(D/DTHETA)(DU/DTHETA) = 16*R**2
        !
        !     ON THE QUARTER-DISK 0 .LT. R .LT. 1, 0 .LT. THETA .LT. PI/2 WITH
        !     WITH THE BOUNDARY CONDITIONS
        !
        !     U(1, THETA) = 1 - COS(4*THETA), 0 .LE. THETA .LE. 1
        !
        !     AND
        !
        !     (DU/DTHETA)(R, 0) = (DU/DTHETA)(R, PI/2) = 0,  0 .LE. R .LE. 1.
        !
        !     (NOTE THAT THE SOLUTION U IS UNSPECIFIED AT R = 0.)
        !          THE R-INTERVAL WILL BE DIVIDED INTO 50 PANELS AND THE
        !     THETA-INTERVAL WILL BE DIVIDED INTO 48 PANELS.
        !
        !
        !     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
        !
        idimf = 100
        a = 0.
        b = 1.
        m = 50
        mbdcnd = 5
        c = 0.
        pi = acos( -1.0 )
        d = pi/2.
        n = 48
        nbdcnd = 3
        elmbda = 0.
        !
        !     AUXILIARY QUANTITIES.
        !
        mp1 = m + 1
        np1 = n + 1
        !
        !     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
        !     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
        !
        do i = 1, mp1
            r(i) = real(i - 1)/50.
        end do
        do j = 1, np1
            theta(j) = real(j - 1)*pi/96.
        end do
        !
        !     GENERATE BOUNDARY DATA.
        !
        bdc(:mp1) = 0.
        bdd(:mp1) = 0.
        !
        !     BDA AND BDB ARE DUMMY VARIABLES.
        !
        do j = 1, np1
            f(mp1, j) = 1. - COS(4.*THETA(j))
        end do
        !
        !     GENERATE RIGHT SIDE OF EQUATION.
        !
        do i = 1, m
            f(i, :np1) = 16.*R(i)**2
        end do

        call hwsplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
            elmbda, f, idimf, pertrb, ierror)
        !
        !     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
        !                U(R, THETA) = R**4*(1 - COS(4*THETA))
        !
        err = 0.
        do i = 1, mp1
            do j = 1, np1
                z = abs(F(i, j)-R(i)**4*(1.-COS(4.*THETA(j))))
                err = max(z, err)
            end do
        end do

        write( *, *) ''
        write( *, *) '    hwsplr *** TEST RUN *** '
        write( *, *) &
            '    Previous 64 bit floating point arithmetic result '
        write( *, *) '    ierror = 0,  discretization error = 6.19134E-4'

        write( *, *) '    The output from your computer is: '
        write( *, *) '    ierror =', ierror, ' discretization error = ', &
            err

    end subroutine test_hwsplr

    subroutine hwsplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror)
        !     file hwsplr.f
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
        !     SUBROUTINE hwsplr (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD,
        !    +                   ELMBDA, F, IDIMF, PERTRB, ierror)
        !
        !
        ! DIMENSION OF           BDA(N), BDB(N), BDC(M), BDD(M), F(IDIMF, N+1)
        ! ARGUMENTS
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES A FINITE DIFFERENCE APPROXIMATION TO
        !                        THE HELMHOLTZ EQUATION IN POLAR COORDINATES.
        !                        THE EQUATION IS
        !
        !                            (1/R)(D/DR)(R(DU/DR)) +
        !                            (1/R**2)(D/DTHETA)(DU/DTHETA) +
        !                            LAMBDA*U = F(R, THETA).
        !
        ! USAGE                  CALL hwsplr (A, B, M, MBDCND, BDA, BDB, C, D, N,
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
        !                          FOR I = 1, 2, ..., M+1,
        !                          WHERE DR = (B-A)/M IS THE PANEL WIDTH.
        !                          M MUST BE GREATER THAN 3.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITION
        !                          AT R = A AND R = B.
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               R = A AND R = B.
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               R = A AND THE DERIVATIVE OF
        !                               THE SOLUTION WITH RESPECT TO R IS
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
        !                          NOTE:
        !                          IF A = 0, DO NOT USE MBDCND = 3 OR 4, BUT
        !                          INSTEAD USE MBDCND = 1, 2, 5, OR 6  .
        !
        !                        BDA
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
        !                          THE SOLUTION WITH RESPECT TO R AT R = A.
        !
        !                          WHEN MBDCND = 3 OR 4,
        !                            BDA(J) = (D/DR)U(A, THETA(J)),
        !                            J = 1, 2, ..., N+1  .
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
        !                          A DUMMY VARIABLE.
        !
        !                        BDB
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE OF
        !                          THE SOLUTION WITH RESPECT TO R AT R = B.
        !
        !                          WHEN MBDCND = 2, 3, OR 6,
        !                            BDB(J) = (D/DR)U(B, THETA(J)),
        !                            J = 1, 2, ..., N+1  .
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
        !                          A DUMMY VARIABLE.
        !
        !                        C, D
        !                          THE RANGE OF THETA, I.E., C .LE.
        !                          THETA .LE. D.  C MUST BE LESS THAN D.
        !
        !                        N
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL (C, D) IS SUBDIVIDED.  HENCE,
        !                          THERE WILL BE N+1 GRID POINTS IN THE
        !                          THETA-DIRECTION GIVEN BY
        !                          THETA(J) = C+(J-1)DTHETA FOR
        !                          J = 1, 2, ..., N+1, WHERE
        !                          DTHETA = (D-C)/N IS THE PANEL WIDTH.
        !                          N MUST BE GREATER THAN 3.
        !
        !                        NBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT THETA = C AND AT THETA = D.
        !
        !                          = 0  IF THE SOLUTION IS PERIODIC IN THETA,
        !                               I.E., U(I, J) = U(I, N+J).
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = C AND THETA = D
        !                               (SEE NOTE BELOW).
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = C AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO THETA IS
        !                               SPECIFIED AT THETA = D
        !                               (SEE NOTE BELOW).
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               AT THETA = C AND THE SOLUTION IS
        !                               SPECIFIED AT THETA = D
        !                               (SEE NOTE BELOW).
        !
        !                          NOTE:
        !                          WHEN NBDCND = 1, 2, OR 4, DO NOT USE
        !                          MBDCND = 5 OR 6
        !                          (THE FORMER INDICATES THAT THE SOLUTION
        !                          IS SPECIFIED AT R = 0, THE LATTER INDICATES
        !                          THE SOLUTION IS UNSPECIFIED AT R = 0).
        !                          USE INSTEAD MBDCND = 1 OR 2  .
        !
        !                        BDC
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE
        !                          OF THE SOLUTION WITH RESPECT TO THETA AT
        !                          THETA = C.  WHEN NBDCND = 3 OR 4,
        !
        !                            BDC(I) = (D/DTHETA)U(R(I), C),
        !                            I = 1, 2, ..., M+1  .
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
        !                          A DUMMY VARIABLE.
        !
        !                        BDD
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIES THE VALUES OF THE DERIVATIVE
        !                          OF THE SOLUTION WITH RESPECT TO THETA AT
        !                          THETA = D.  WHEN NBDCND = 2 OR 3,
        !
        !                            BDD(I) = (D/DTHETA)U(R(I), D),
        !                            I = 1, 2, ..., M+1  .
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
        !                          A DUMMY VARIABLE.
        !
        !                        ELMBDA
        !                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
        !                          EQUATION.  IF LAMBDA .LT. 0, A SOLUTION
        !                          MAY NOT EXIST.  HOWEVER, hwsplr WILL
        !                          ATTEMPT TO FIND A SOLUTION.
        !
        !                        F
        !                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
        !                          LEAST (M+1)*(N+1), SPECIFYING VALUES
        !                          OF THE RIGHT SIDE OF THE HELMHOLTZ
        !                          EQUATION AND BOUNDARY DATA (IF ANY).
        !
        !                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
        !                          FOR I = 2, 3, ..., M AND J = 2, 3, ..., N
        !                          F(I, J) = F(R(I), THETA(J)).
        !
        !                          ON THE BOUNDARIES F IS DEFINED AS FOLLOWS:
        !                          FOR J = 1, 2, ..., N+1 AND I = 1, 2, ..., M+1
        !
        !                          MBDCND   F(1, J)            F(M+1, J)
        !                          ------   -------------     -------------
        !
        !                            1      U(A, THETA(J))     U(B, THETA(J))
        !                            2      U(A, THETA(J))     F(B, THETA(J))
        !                            3      F(A, THETA(J))     F(B, THETA(J))
        !                            4      F(A, THETA(J))     U(B, THETA(J))
        !                            5      F(0, 0)            U(B, THETA(J))
        !                            6      F(0, 0)            F(B, THETA(J))
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
        !                          THEN THE SOLUTION MUST BE SPECIFIED.
        !
        !                        IDIMF
        !                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
        !                          F AS IT APPEARS IN THE PROGRAM CALLING
        !                          hwsplr.  THIS PARAMETER IS USED TO SPECIFY
        !                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
        !                          BE AT LEAST M+1.
        !
        ! ON OUTPUT              F
        !                          CONTAINS THE SOLUTION U(I, J) OF THE FINITE
        !                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
        !                          (R(I), THETA(J)),
        !                          I = 1, 2, ..., M+1, J = 1, 2, ..., N+1  .
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
        !                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
        !                          SPECIFIED FOR A POISSON EQUATION
        !                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
        !                          PERTRB IS A CONSTANT, CALCULATED AND
        !                          SUBTRACTED FROM F, WHICH ENSURES THAT A
        !                          SOLUTION EXISTS.  hwsplr THEN COMPUTES
        !                          THIS SOLUTION, WHICH IS A LEAST SQUARES
        !                          SOLUTION TO THE ORIGINAL APPROXIMATION.
        !                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
        !                          A SOLUTION.  HENCE, THE SOLUTION IS NOT
        !                          UNIQUE.  PERTRB SHOULD BE SMALL COMPARED
        !                          TO THE RIGHT SIDE. OTHERWISE, A SOLUTION
        !                          IS OBTAINED TO AN ESSENTIALLY DIFFERENT
        !                          PROBLEM.  THIS COMPARISON SHOULD ALWAYS
        !                          BE MADE TO INSURE THAT A MEANINGFUL
        !                          SOLUTION HAS BEEN OBTAINED.
        !
        !                        ierror
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 11,
        !                          A SOLUTION IS NOT ATTEMPTED.
        !
        !                          =  0  NO ERROR.
        !                          =  1  A .LT. 0  .
        !                          =  2  A .GE. B.
        !                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6  .
        !                          =  4  C .GE. D.
        !                          =  5  N .LE. 3
        !                          =  6  NBDCND .LT. 0 OR .GT. 4  .
        !                          =  7  A = 0, MBDCND = 3 OR 4  .
        !                          =  8  A .GT. 0, MBDCND .GE. 5  .
        !                          =  9  MBDCND .GE. 5, NBDCND .NE. 0
        !                                AND NBDCND .NE. 3  .
        !                          = 10  IDIMF .LT. M+1  .
        !                          = 11  LAMBDA .GT. 0  .
        !                          = 12  M .LE. 3
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
        !
        !                          SINCE THIS IS THE ONLY MEANS OF INDICATING
        !                          A POSSIBLY INCORRECT CALL TO hwsplr, THE
        !                          USER SHOULD TEST ierror AFTER THE CALL.
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
        ! REFERENCES             SWARZTRAUBER, P. AND R. SWEET, "EFFICIENT
        !                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
        !                        ELLIPTIC EQUATIONS"
        !                          NCAR TN/IA-109, JULY, 1975, 138 PP.
        !***********************************************************************
        type (FishpackWorkspace) :: w
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer  :: m
        integer  :: mbdcnd
        integer  :: n
        integer  :: nbdcnd
        integer  :: idimf
        integer  :: ierror
        real  :: a
        real  :: b
        real  :: c
        real  :: d
        real  :: elmbda
        real  :: pertrb
        real  :: bda(*)
        real  :: bdb(*)
        real  :: bdc(*)
        real  :: bdd(*)
        real  :: f(idimf, *)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer :: irwk, icwk
        !-----------------------------------------------
        !
        !     CHECK FOR INVALID PARAMETERS.
        !
        ierror = 0
        if (a < 0.) ierror = 1
        if (a >= b) ierror = 2
        if (mbdcnd<=0 .or. mbdcnd>=7) ierror = 3
        if (c >= d) ierror = 4
        if (n <= 3) ierror = 5
        if (nbdcnd<=(-1) .or. nbdcnd>=5) ierror = 6
        if (a==0. .and. (mbdcnd==3 .or. mbdcnd==4)) ierror = 7
        if (a>0. .and. mbdcnd>=5) ierror = 8
        if (mbdcnd>=5 .and. nbdcnd/=0 .and. nbdcnd/=3) ierror = 9
        if (idimf < m + 1) ierror = 10
        if (m <= 3) ierror = 12
        if (ierror /= 0) return
        !     compute and allocate required work space
        irwk=4*(n+1)+(m+1)*(13+INT(log(real(n+1))/log(2.0)))
        icwk = 0
        call w%create( irwk, icwk, ierror )
        !     check that allocation was successful
        if (ierror == 20) return
        call hwsplrr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
            elmbda, f, idimf, pertrb, ierror, w%rew)
        !     release allocated work space
        call w%destroy()

    end subroutine hwsplr

    subroutine hwsplrR(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc &
        , bdd, elmbda, f, idimf, pertrb, ierror, w)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer , intent (in) :: m
        integer , intent (in) :: mbdcnd
        integer , intent (in) :: n
        integer  :: nbdcnd
        integer  :: idimf
        integer , intent (out) :: ierror
        real , intent (in) :: a
        real , intent (in) :: b
        real , intent (in) :: c
        real , intent (in) :: d
        real , intent (in) :: elmbda
        real , intent (out) :: pertrb
        real , intent (in) :: bda(*)
        real , intent (in) :: bdb(*)
        real , intent (in) :: bdc(*)
        real , intent (in) :: bdd(*)
        real  :: f(idimf, *)
        real  :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer :: mp1, np1, np, mstart, mstop, munk, nstart, nstop, nunk &
            , id2, id3, id4, id5, id6, ij, i, j, l, lp, k, i1, ierr1, ip
        real :: deltar, dlrby2, dlrsq, deltht, dlthsq, a1, r, s2, a2, &
             s, s1, ypole

        !real, save ::  alll
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
        mstop = mp1
        go to (101, 105, 102, 103, 104, 105) mbdcnd
101 continue
    mstop = m
    go to 105
102 continue
    mstart = 1
    go to 105
103 continue
    mstart = 1
104 continue
    mstop = m
105 continue
    munk = mstop - mstart + 1
    nstart = 1
    nstop = n
    go to (109, 106, 107, 108, 109) np
106 continue
    nstart = 2
    go to 109
107 continue
    nstart = 2
108 continue
    nstop = np1
109 continue
    nunk = nstop - nstart + 1
    !
    !     DEFINE A, B, C COEFFICIENTS IN W-ARRAY.
    !
    id2 = munk
    id3 = id2 + munk
    id4 = id3 + munk
    id5 = id4 + munk
    id6 = id5 + munk
    a1 = 2./dlrsq
    ij = 0
    if (mbdcnd==3 .or. mbdcnd==4) ij = 1
    do i = 1, munk
        r = a + real(i - ij)*deltar
        j = id5 + i
        w(j) = r
        j = id6 + i
        w(j) = 1./r**2
        w(i) = (r - dlrby2)/(r*dlrsq)
        j = id3 + i
        w(j) = (r + dlrby2)/(r*dlrsq)
        j = id2 + i
        w(j) = (-a1) + elmbda
    end do
    go to (114, 111, 112, 113, 114, 111) mbdcnd
111 continue
    w(id2) = a1
    go to 114
112 continue
    w(id2) = a1
113 continue
    w(id3+1) = a1
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
!     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
!
124 continue
    a1 = 1./dlthsq
    l = id5 - mstart + 1
    lp = id6 - mstart + 1
    go to (134, 125, 125, 127, 127) np
125 continue
    f(mstart:mstop, 2) = F(mstart:mstop, 2) - a1*W(mstart+lp:mstop+lp)*F &
        (mstart:mstop, 1)
    go to 129
127 continue
    a1 = 2./deltht
    f(mstart:mstop, 1) = F(mstart:mstop, 1) + a1*W(mstart+lp:mstop+lp)* &
        BDC(mstart:mstop)
129 continue
    a1 = 1./dlthsq
    go to (134, 130, 132, 132, 130) np
130 continue
    f(mstart:mstop, n) = F(mstart:mstop, n) - a1*W(mstart+lp:mstop+lp)*F &
        (mstart:mstop, np1)
    go to 134
132 continue
    a1 = 2./deltht
    f(mstart:mstop, np1) = F(mstart:mstop, np1) - a1*W(mstart+lp:mstop+ &
        lp)*BDD(mstart:mstop)
134 continue
    if (mbdcnd>=5 .and. nbdcnd==3) f(1, 1) = F(1, 1) - (BDD(2)-BDC(2))* &
        4./(real(n)*deltht*dlrsq)
    !
    !     ADJUST RIGHT SIDE OF SINGULAR PROBLEMS TO INSURE EXISTENCE OF A
    !     SOLUTION.
    !
    pertrb = 0.
    if (elmbda >= 0.) then
        if (elmbda /= 0.) then
            ierror = 11
        else
            if (nbdcnd==0 .or. nbdcnd==3) then
                s2 = 0.
                go to (144, 144, 137, 144, 144, 138) mbdcnd
137         continue
            w(id5+1) = 0.5*(W(id5+2)-dlrby2)
            s2 = 0.25*deltar
138     continue
        a2 = 2.
        if (nbdcnd == 0) a2 = 1.
        j = id5 + munk
        w(j) = 0.5*(W(j-1)+dlrby2)
        s = 0.
        do i = mstart, mstop
            s1 = 0.
            ij = nstart + 1
            k = nstop - 1
            s1 = SUM(F(i, ij:k))
            j = i + l
            s = s + (a2*s1 + F(i, nstart)+F(i, nstop))*W(j)
        end do
        s2=real(m)*a+deltar*(real((m-1)*(m+1))*0.5+0.25)+s2
        s1 = (2. + a2*real(nunk - 2))*s2
        if (mbdcnd /= 3) then
            s2 = real(n)*a2*deltar/8.
            s = s + F(1, 1)*s2
            s1 = s1 + s2
        end if
        pertrb = s/s1
        f(mstart:mstop, nstart:nstop) = F(mstart:mstop, nstart: &
            nstop) - pertrb
    end if
end if
end if
144 continue
    do i = mstart, mstop
        k = i - mstart + 1
        j = i + lp
        a1 = dlthsq/W(j)
        w(k) = a1*W(k)
        j = id2 + k
        w(j) = a1*W(j)
        j = id3 + k
        w(j) = a1*W(j)
        f(i, nstart:nstop) = a1*F(i, nstart:nstop)
    end do
    w(1) = 0.
    w(id4) = 0.
    !
    !     SOLVE THE SYSTEM OF EQUATIONS.
    !
    i1 = 1
    ierr1 = 0
    call genbunn (nbdcnd, nunk, i1, munk, W(1), W(id2+1), W(id3+1), &
        idimf, F(mstart, nstart), ierr1, W(id4+1))
    go to (157, 157, 157, 157, 148, 147) mbdcnd
!
!     ADJUST THE SOLUTION AS NECESSARY FOR THE PROBLEMS WHERE A = 0.
!
147 continue
    if (elmbda /= 0.) go to 148
    ypole = 0.
    go to 155
148 continue
    j = id5 + munk
    w(j) = W(id2)/W(id3)
    do ip = 3, munk
        i = munk - ip + 2
        j = id5 + i
        lp = id2 + i
        k = id3 + i
        w(j) = W(i)/(W(lp)-W(k)*W(j+1))
    end do
    w(id5+1) = -0.5*dlthsq/(W(id2+1)-W(id3+1)*W(id5+2))
    do i = 2, munk
        j = id5 + i
        w(j) = -W(j)*W(j-1)
    end do
    s = 0.
    s = SUM(F(2, nstart:nstop))
    a2 = nunk
    if (nbdcnd /= 0) then
        s = s - 0.5*(F(2, nstart)+F(2, nstop))
        a2 = a2 - 1.
    end if
    ypole = (0.25*dlrsq*F(1, 1)-s/a2)/(W(id5+1)-1.+elmbda*dlrsq*0.25)
    do i = mstart, mstop
        k = l + i
        f(i, nstart:nstop) = F(i, nstart:nstop) + ypole*W(k)
    end do
155 continue
    f(1, :np1) = ypole
157 continue
    if (nbdcnd == 0) then
        f(mstart:mstop, np1) = F(mstart:mstop, 1)
    end if

end subroutine hwsplrR

end module module_hwsplr
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
