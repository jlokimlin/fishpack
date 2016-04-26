module module_hstssp

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_genbun, only: &
        genbunn

    use module_poistg, only: &
        poistgg

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hstssp


contains


    subroutine hstssp( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror )
        !
        !     file hstssp.f
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
        !     SUBROUTINE hstssp (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD,
        !    +                   ELMBDA, F, IDIMF, PERTRB, ierror)
        !
        !
        ! DIMENSION OF           BDA(N), BDB(N), BDC(M), BDD(M), F(IDIMF, N)
        ! ARGUMENTS
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
        !                        DIFFERENCE APPROXIMATION ON A STAGGERED GRID
        !                        TO THE HELMHOLTZ EQUATION IN SPHERICAL
        !                        COORDINATES AND ON THE SURFACE OF THE UNIT
        !                        SPHERE (RADIUS OF 1).  THE EQUATION IS
        !
        !                          (1/SIN(THETA))(D/DTHETA)(SIN(THETA)
        !                          (DU/DTHETA)) + (1/SIN(THETA)**2)
        !                          (D/DPHI)(DU/DPHI) + LAMBDA*U = F(THETA, PHI)
        !
        !                        WHERE THETA IS COLATITUDE AND PHI IS
        !                        LONGITUDE.
        !
        ! USAGE                  CALL hstssp (A, B, M, MBDCND, BDA, BDB, C, D, N,
        !                                     NBDCND, BDC, BDD, ELMBDA, F, IDIMF,
        !                                     PERTRB, ierror)
        !
        !
        ! ARGUMENTS
        ! ON INPUT
        !
        !                        A, B
        !                          THE RANGE OF THETA (COLATITUDE),
        !                          I.E. A .LE. THETA .LE. B.
        !                          A MUST BE LESS THAN B AND A MUST BE
        !                          NON-NEGATIVE.  A AND B ARE IN RADIANS.
        !                          A = 0 CORRESPONDS TO THE NORTH POLE AND
        !                          B = PI CORRESPONDS TO THE SOUTH POLE.
        !
        !
        !                            * * *  IMPORTANT  * * *
        !
        !                          IF B IS EQUAL TO PI, THEN B MUST BE
        !                          COMPUTED USING THE STATEMENT
        !                            B = PI_MACH(DUM)
        !
        !                          THIS INSURES THAT B IN THE USER"S PROGRAM
        !                          IS EQUAL TO PI IN THIS PROGRAM WHICH
        !                          PERMITS SEVERAL TESTS OF THE INPUT
        !                          PARAMETERS THAT OTHERWISE WOULD NOT BE
        !                          POSSIBLE.
        !
        !                            * * * * * * * * * * * *
        !                        M
        !                          THE NUMBER OF GRID POINTS IN THE INTERVAL
        !                          (A, B).  THE GRID POINTS IN THE THETA
        !                          DIRECTION ARE GIVEN BY
        !                          THETA(I) = A + (I-0.5)DTHETA
        !                          FOR I=1, 2, ..., M WHERE DTHETA =(B-A)/M.
        !                          M MUST BE GREATER THAN 2.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT THETA = A AND THETA = B.
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = A AND THETA = B.
        !                               (SEE NOTE 3 BELOW)
        !
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = A AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO THETA IS
        !                               SPECIFIED AT THETA = B
        !                               (SEE NOTES 2 AND 3 BELOW).
        !
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                             WITH RESPECT TO THETA IS SPECIFIED
        !                             AT THETA = A
        !                             (SEE NOTES 1, 2 BELOW) AND THETA = B.
        !
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               AT THETA = A
        !                               (SEE NOTES 1 AND 2 BELOW) AND THE
        !                               SOLUTION IS SPECIFIED AT THETA = B.
        !
        !                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = A = 0 AND THE SOLUTION IS
        !                               SPECIFIED AT THETA = B.
        !                               (SEE NOTE 3 BELOW)
        !
        !                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = A = 0 AND THE DERIVATIVE
        !                               OF THE SOLUTION WITH RESPECT TO THETA
        !                               IS SPECIFIED AT THETA = B
        !                               (SEE NOTE 2 BELOW).
        !
        !                          = 7  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = A AND THE SOLUTION IS
        !                               UNSPECIFIED AT THETA = B = PI.
        !                               (SEE NOTE 3 BELOW)
        !
        !                          = 8  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED AT
        !                               THETA = A (SEE NOTE 1 BELOW)
        !                               AND THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = B = PI.
        !
        !                          = 9  IF THE SOLUTION IS UNSPECIFIED AT
        !                               THETA = A = 0 AND THETA = B = PI.
        !
        !                          NOTE 1:
        !                          IF A = 0, DO NOT USE MBDCND = 3, 4, OR 8,
        !                          BUT INSTEAD USE MBDCND = 5, 6, OR 9.
        !
        !                          NOTE 2:
        !                          IF B = PI, DO NOT USE MBDCND = 2, 3, OR 6,
        !                          BUT INSTEAD USE MBDCND = 7, 8, OR 9.
        !
        !                          NOTE 3:
        !                          WHEN THE SOLUTION IS SPECIFIED AT
        !                          THETA = 0 AND/OR THETA = PI AND THE OTHER
        !                          BOUNDARY CONDITIONS ARE COMBINATIONS
        !                          OF UNSPECIFIED, NORMAL DERIVATIVE, OR
        !                          PERIODICITY A SINGULAR SYSTEM RESULTS.
        !                          THE UNIQUE SOLUTION IS DETERMINED BY
        !                          EXTRAPOLATION TO THE SPECIFICATION OF THE
        !                          SOLUTION AT EITHER THETA = 0 OR THETA = PI.
        !                          BUT IN THESE CASES THE RIGHT SIDE OF THE
        !                          SYSTEM  WILL BE PERTURBED BY THE CONSTANT
        !                          PERTRB.
        !
        !                        BDA
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
        !                          SPECIFIES THE BOUNDARY VALUES (IF ANY) OF
        !                          THE SOLUTION AT THETA = A.
        !
        !                          WHEN MBDCND = 1, 2, OR 7,
        !                            BDA(J) = U(A, PHI(J)) ,      J=1, 2, ..., N.
        !
        !                          WHEN MBDCND = 3, 4, OR 8,
        !                            BDA(J) = (D/DTHETA)U(A, PHI(J)) ,
        !                            J=1, 2, ..., N.
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE,
        !                          BDA IS A DUMMY VARIABLE.
        !
        !                        BDB
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
        !                          SPECIFIES THE BOUNDARY VALUES OF THE
        !                          SOLUTION AT THETA = B.
        !
        !                          WHEN MBDCND = 1, 4, OR 5,
        !                            BDB(J) = U(B, PHI(J)) ,       J=1, 2, ..., N.
        !
        !                          WHEN MBDCND = 2, 3, OR 6,
        !                            BDB(J) = (D/DTHETA)U(B, PHI(J)) ,
        !                            J=1, 2, ..., N.
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
        !                          A DUMMY VARIABLE.
        !
        !                        C, D
        !                          THE RANGE OF PHI (LONGITUDE),
        !                          I.E. C .LE. PHI .LE. D.
        !                          C MUST BE LESS THAN D.  IF D-C = 2*PI,
        !                          PERIODIC BOUNDARY CONDITIONS ARE USUALLY
        !                          USUALLY PRESCRIBED.
        !
        !                        N
        !                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
        !                          (C, D).  THE UNKNOWNS IN THE PHI-DIRECTION
        !                          ARE GIVEN BY PHI(J) = C + (J-0.5)DPHI,
        !                          J=1, 2, ..., N, WHERE DPHI = (D-C)/N.
        !                          N MUST BE GREATER THAN 2.
        !
        !                        NBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT PHI = C  AND PHI = D.
        !
        !                          = 0  IF THE SOLUTION IS PERIODIC IN PHI,
        !                               I.E.  U(I, J) = U(I, N+J).
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               PHI = C AND PHI = D
        !                               (SEE NOTE BELOW).
        !
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               PHI = C AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO PHI IS
        !                               SPECIFIED AT PHI = D
        !                               (SEE NOTE BELOW).
        !
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO PHI IS SPECIFIED
        !                               AT PHI = C AND PHI = D.
        !
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO PHI IS SPECIFIED
        !                               AT PHI = C AND THE SOLUTION IS
        !                               SPECIFIED AT PHI = D
        !                               (SEE NOTE BELOW).
        !
        !                          NOTE:
        !                          WHEN NBDCND = 1, 2, OR 4, DO NOT USE
        !                          MBDCND = 5, 6, 7, 8, OR 9
        !                          (THE FORMER INDICATES THAT THE SOLUTION
        !                          IS SPECIFIED AT A POLE; THE LATTER
        !                          INDICATES THE SOLUTION IS UNSPECIFIED).
        !                          USE INSTEAD MBDCND = 1 OR 2.
        !
        !                        BDC
        !                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
        !                          SPECIFIES THE BOUNDARY VALUES OF THE
        !                          SOLUTION AT PHI = C.
        !
        !                          WHEN NBDCND = 1 OR 2,
        !                            BDC(I) = U(THETA(I), C) ,     I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 3 OR 4,
        !                            BDC(I) = (D/DPHI)U(THETA(I), C),
        !                            I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
        !
        !                        BDD
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
        !                          SPECIFIES THE BOUNDARY VALUES OF THE
        !                          SOLUTION AT PHI = D.
        !
        !                          WHEN NBDCND = 1 OR 4,
        !                            BDD(I) = U(THETA(I), D) ,     I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 2 OR 3,
        !                            BDD(I) = (D/DPHI)U(THETA(I), D) ,
        !                            I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
        !
        !                        ELMBDA
        !                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
        !                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
        !                          A SOLUTION MAY NOT EXIST.  HOWEVER,
        !                          hstssp WILL ATTEMPT TO FIND A SOLUTION.
        !
        !                        F
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
        !                          THE VALUES OF THE RIGHT SIDE OF THE
        !                          HELMHOLTZ EQUATION.
        !                          FOR I=1, 2, ..., M AND J=1, 2, ..., N
        !
        !                            F(I, J) = F(THETA(I), PHI(J)) .
        !
        !                          F MUST BE DIMENSIONED AT LEAST M X N.
        !
        !                        IDIMF
        !                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
        !                          F AS IT APPEARS IN THE PROGRAM CALLING
        !                          hstssp.  THIS PARAMETER IS USED TO SPECIFY
        !                          THE VARIABLE DIMENSION OF F.
        !                          IDIMF MUST BE AT LEAST M.
        !
        !
        ! ON OUTPUT              F
        !                          CONTAINS THE SOLUTION U(I, J) OF THE FINITE
        !                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
        !                          (THETA(I), PHI(J)) FOR
        !                          I=1, 2, ..., M, J=1, 2, ..., N.
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
        !                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
        !                          SPECIFIED FOR A POISSON EQUATION
        !                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
        !                          PERTRB IS A CONSTANT, CALCULATED AND
        !                          SUBTRACTED FROM F, WHICH ENSURES THAT A
        !                          SOLUTION EXISTS.  hstssp THEN COMPUTES
        !                          THIS SOLUTION, WHICH IS A LEAST SQUARES
        !                          SOLUTION TO THE ORIGINAL APPROXIMATION.
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
        !                          PARAMETERS. EXCEPT TO NUMBERS 0 AND 14,
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
        !                          =  4  C .GE. D
        !
        !                          =  5  N .LE. 2
        !
        !                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4
        !
        !                          =  7  A .GT. 0 AND MBDCND = 5, 6, OR 9
        !
        !                          =  8  A = 0 AND MBDCND = 3, 4, OR 8
        !
        !                          =  9  B .LT. PI AND MBDCND .GE. 7
        !
        !                          = 10  B = PI AND MBDCND = 2, 3, OR 6
        !
        !                          = 11  MBDCND .GE. 5 AND NDBCND = 1, 2, OR 4
        !
        !                          = 12  IDIMF .LT. M
        !
        !                          = 13  M .LE. 2
        !
        !                          = 14  LAMBDA .GT. 0
        !
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
        !
        !                          SINCE THIS IS THE ONLY MEANS OF INDICATING
        !                          A POSSIBLY INCORRECT CALL TO hstssp, THE
        !                          USER SHOULD TEST ierror AFTER THE CALL.
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED FILES         fish.f, comf.f, genbun.f, gnbnaux.f, poistg.f
        !
        ! LANGUAGE               FORTRAN 90
        !
        ! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
        !                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
        !                        IN JANUARY 1980.
        !                        Revised in June 2004 by John Adams using
        !                        Fortran 90 dynamically allocated work space.
        !
        ! PORTABILITY            FORTRAN 90.
        !
        ! ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-
        !                        DIFFERENCE EQUATIONS, INCORPORATES BOUNDARY
        !                        DATA, ADJUSTS THE RIGHT SIDE WHEN THE SYSTEM
        !                        IS SINGULAR AND CALLS EITHER POISTG OR genbun
        !                        WHICH SOLVES THE LINEAR SYSTEM OF EQUATIONS.
        !
        ! TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
        !                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
        !
        ! ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
        !                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
        !                        DIGITS FOR N AND M AS LARGE AS 64.
        !                        MORE DETAILED INFORMATION ABOUT ACCURACY
        !                        CAN BE FOUND IN THE DOCUMENTATION FOR
        !                        ROUTINE POISTG WHICH IS THE ROUTINE THAT
        !                        ACTUALLY SOLVES THE FINITE DIFFERENCE
        !                        EQUATIONS.
        !
        ! REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD
        !                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
        !                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
        !                        GRID OF ARBITRARY SIZE, " J. COMP. PHYS.
        !                        20(1976), PP. 171-182.
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
        integer (ip)             :: irwk
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------

        ! initialize error flag
        ierror = 0

        ! Check input arguments
        associate( pi => acos(-1.0_wp))
            if (a<0. .or. b>pi) ierror = 1
            if (a >= b) ierror = 2
            if (mbdcnd<=0 .or. mbdcnd>9) ierror = 3
            if (c >= d) ierror = 4
            if (n <= 2) ierror = 5
            if (nbdcnd<0 .or. nbdcnd>=5) ierror = 6
            if(a>0..and.(mbdcnd==5.or.mbdcnd==6.or.mbdcnd==9))ierror=7
            if(a==0..and.(mbdcnd==3.or.mbdcnd==4.or.mbdcnd==8))ierror=8
            if (b<pi .and. mbdcnd>=7) ierror = 9
            if(b==pi.and.(mbdcnd==2.or.mbdcnd==3.or.mbdcnd==6))ierror=10
            if (mbdcnd>=5 .and. (nbdcnd==1 .or. nbdcnd==2 .or. nbdcnd==4)) &
                ierror = 11
            if (idimf < m) ierror = 12
            if (m <= 2) ierror = 13
            if (ierror /= 0) return
        end associate

        ! compute required real workspace size
        call workspace%get_genbun_workspace_dimensions( n, m, irwk )
        irwk = irwk + 3 * m

        ! allocate real workspace
        associate( icwk => 0 )
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! Check that allocation was successful
        if (ierror == 20) then
            return
        end if

        ! solve system
        associate( rew => workspace%real_workspace )
            call hstsspp(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hstssp


    subroutine hstsspp( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, w )
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip) :: m
        integer , intent (in) :: mbdcnd
        integer (ip) :: n
        integer , intent (in) :: nbdcnd
        integer (ip) :: idimf
        integer , intent (out) :: ierror
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
        real (wp) :: f(idimf, 1)
        real (wp), intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer::np, isw, jsw, mb, iwb, iwc, iwr, iws, i, j, mm1, k, lp, ierr1, i1
        real :: deltar, dlrsq, deltht, dlthsq, pi, a1, a2, a3
        !-----------------------------------------------

        deltar = (b - a)/real(m)
        dlrsq = deltar**2
        deltht = (d - c)/real(n)
        dlthsq = deltht**2
        np = nbdcnd + 1
        isw = 1
        jsw = 1
        mb = mbdcnd
        if (elmbda == 0.) then
            go to (101, 102, 105, 103, 101, 105, 101, 105, 105) mbdcnd
101     continue
        if (a/=0. .or. b/=pi) go to 105
        mb = 9
        go to 104
102 continue
    if (a /= 0.) go to 105
    mb = 6
    go to 104
103 continue
    if (b /= pi) go to 105
    mb = 8
104 continue
    jsw = 2
end if
105 continue
    iwb = m
    iwc = iwb + m
    iwr = iwc + m
    iws = iwr + m
    do i = 1, m
        j = iwr + i
        w(j) = SIN(a + (real(i) - 0.5)*deltar)
        w(i) = SIN(a + real(i - 1)*deltar)/dlrsq
    end do
    mm1 = m - 1
    w(iwc+1:mm1+iwc) = W(2:mm1+1)
    w(iwb+1:mm1+iwb) = elmbda*W(iwr+1:mm1+iwr) - (W(:mm1)+W(2:mm1+1))
    w(iwr) = SIN(b)/dlrsq
    w(iwc) = elmbda*W(iws) - (W(m)+W(iwr))
    do i = 1, m
        j = iwr + i
        a1 = W(j)
        f(i, :n) = a1*F(i, :n)
    end do
    !
    !     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
    !
    go to (110, 110, 112, 112, 114, 114, 110, 112, 114) mb
110 continue
    a1 = 2.*W(1)
    w(iwb+1) = W(iwb+1) - W(1)
    f(1, :n) = F(1, :n) - a1*BDA(:n)
    go to 114
112 continue
    a1 = deltar*W(1)
    w(iwb+1) = W(iwb+1) + W(1)
    f(1, :n) = F(1, :n) + a1*BDA(:n)
114 continue
    go to (115, 117, 117, 115, 115, 117, 119, 119, 119) mb
115 continue
    a1 = 2.*W(iwr)
    w(iwc) = W(iwc) - W(iwr)
    f(m, :n) = F(m, :n) - a1*BDB(:n)
    go to 119
117 continue
    a1 = deltar*W(iwr)
    w(iwc) = W(iwc) + W(iwr)
    f(m, :n) = F(m, :n) - a1*BDB(:n)
!
!     ENTER BOUNDARY DATA FOR PHI-BOUNDARIES.
!
119 continue
    a1 = 2./dlthsq
    go to (129, 120, 120, 122, 122) np
120 continue
    f(:m, 1) = F(:m, 1) - a1*BDC(:m)/W(iwr+1:m+iwr)
    go to 124
122 continue
    a1 = 1./deltht
    f(:m, 1) = F(:m, 1) + a1*BDC(:m)/W(iwr+1:m+iwr)
124 continue
    a1 = 2./dlthsq
    go to (129, 125, 127, 127, 125) np
125 continue
    f(:m, n) = F(:m, n) - a1*BDD(:m)/W(iwr+1:m+iwr)
    go to 129
127 continue
    a1 = 1./deltht
    f(:m, n) = F(:m, n) - a1*BDD(:m)/W(iwr+1:m+iwr)
129 continue
    pertrb = 0.
    if (elmbda >= 0.) then
        if (elmbda /= 0.) then
            ierror = 14
        else
            go to (139, 139, 132, 139, 139, 132, 139, 132, 132) mb
132     continue
        go to (133, 139, 139, 133, 139) np
133 continue
    isw = 2
    do j = 1, n
        pertrb = pertrb + sum(F(:m, j))
    end do
    a1 = real(n)*(COS(a) - COS(b))/(2.*SIN(0.5*deltar))
    pertrb = pertrb/a1
    do i = 1, m
        j = iwr + i
        a1 = pertrb*W(j)
        f(i, :n) = F(i, :n) - a1
    end do
    a2 = sum(F(1, :n))
    a3 = sum(F(m, :n))
    a2 = a2/W(iwr+1)
    a3 = a3/W(iws)
end if
end if
139 continue
    do i = 1, m
        j = iwr + i
        a1 = dlthsq*W(j)
        w(i) = a1*W(i)
        j = iwc + i
        w(j) = a1*W(j)
        j = iwb + i
        w(j) = a1*W(j)
        f(i, :n) = a1*F(i, :n)
    end do
    lp = nbdcnd
    w(1) = 0.
    w(iwr) = 0.
    !
    !     CALL POISTG OR genbun TO SOLVE THE SYSTEM OF EQUATIONS.
    !
    ierr1 = 0
    i1 = 1
    if (nbdcnd /= 0) then
        call POISTGG(lp, n, i1, m, w, W(iwb+1), W(iwc+1), idimf, f, &
            ierr1, W(iwr+1))
    else
        call genbunn(lp, n, i1, m, w, W(iwb+1), W(iwc+1), idimf, f, &
            ierr1, W(iwr+1))
    end if
    if (isw==2 .and. jsw==2) then
        if (mb == 8) then
            a1 = sum(F(m, :n))
            a1 = (a1 - dlrsq*a3/16.)/real(n)
            if (nbdcnd == 3) a1 = a1 + (BDD(m)-BDC(m))/(d - c)
            a1 = BDB(1) - a1
        else
            a1 = sum(F(1, :n))
            a1 = (a1 - dlrsq*a2/16.)/real(n)
            if (nbdcnd == 3) a1 = a1 + (BDD(1)-BDC(1))/(d - c)
            a1 = BDA(1) - a1
        end if
        f(:m, :n) = F(:m, :n) + a1
    end if

end subroutine hstsspP
    !
    !*****************************************************************************************
    !
end module module_hstssp
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 changes
!-----------------------------------------------------------------------

