!
!     file hstcyl.f
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
!     SUBROUTINE hstcyl (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD, 
!    +                   ELMBDA, F, IDIMF, PERTRB, ierror)
!
! DIMENSION OF           BDA(N), BDB(N), BDC(M), BDD(M), F(IDIMF, N)
! ARGUMENTS
!
! LATEST REVISION        June 2004
!
! PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
!                        DIFFERENCE APPROXIMATION ON A STAGGERED
!                        GRID TO THE MODIFIED HELMHOLTZ EQUATION
!                        IN CYLINDRICAL COORDINATES. THIS EQUATION
!
!                          (1/R)(D/DR)(R(DU/DR)) + (D/DZ)(DU/DZ)
!
!                          + LAMBDA*(1/R**2)*U = F(R, Z)
!
!                        IS A TWO-DIMENSIONAL MODIFIED HELMHOLTZ
!                        EQUATION RESULTING FROM THE FOURIER TRANSFORM
!                        OF A THREE-DIMENSIONAL POISSON EQUATION.
!
! USAGE                  CALL hstcyl (A, B, M, MBDCND, BDA, BDB, C, D, N, 
!                                     NBDCND, BDC, BDD, ELMBDA, F, IDIMF, 
!                                     PERTRB, ierror)
!
! ARGUMENTS
! ON INPUT               A, B
!
!                          THE RANGE OF R, I.E. A .LE. R .LE. B.
!                          A MUST BE LESS THAN B AND A MUST BE
!                          BE NON-NEGATIVE.
!
!                        M
!                          THE NUMBER OF GRID POINTS IN THE INTERVAL
!                          (A, B).  THE GRID POINTS IN THE R-DIRECTION
!                          R-DIRECTION ARE GIVEN BY
!                          R(I) = A + (I-0.5)DR FOR I=1, 2, ..., M
!                          WHERE DR =(B-A)/M.
!                          M MUST BE GREATER THAN 2.
!
!                        MBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT R = A AND R = B.
!
!                          = 1  IF THE SOLUTION IS SPECIFIED AT R = A
!                               (SEE NOTE BELOW) AND R = B.
!
!                          = 2  IF THE SOLUTION IS SPECIFIED AT R = A
!                               (SEE NOTE BELOW) AND THE DERIVATIVE
!                               OF THE SOLUTION WITH RESPECT TO R IS
!                               SPECIFIED AT R = B.
!
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO R IS SPECIFIED AT
!                               R = A (SEE NOTE BELOW) AND R = B.
!
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO R IS SPECIFIED AT
!                               R = A (SEE NOTE BELOW) AND THE
!                               SOLUTION IS SPECIFIED AT R = B.
!
!                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
!                               R = A = 0 AND THE SOLUTION IS
!                               SPECIFIED AT R = B.
!
!                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
!                               R = A = 0 AND THE DERIVATIVE OF THE
!                               SOLUTION WITH RESPECT TO R IS SPECIFIED
!                               AT R = B.
!
!                          NOTE:
!                          IF A = 0, DO NOT USE MBDCND = 1, 2, 3, OR 4, 
!                          BUT INSTEAD USE MBDCND = 5 OR 6.
!                          THE RESULTING APPROXIMATION GIVES THE ONLY
!                          MEANINGFUL BOUNDARY CONDITION, 
!                          I.E. DU/DR = 0.
!                          (SEE D. GREENSPAN, 'INTRODUCTORY NUMERICAL
!                          ANALYSIS OF ELLIPTIC BOUNDARY VALUE
!                          PROBLEMS, ' HARPER AND ROW, 1965, CHAPTER 5.)
!
!                        BDA
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
!                          SPECIFIES THE BOUNDARY VALUES (IF ANY)
!                          OF THE SOLUTION AT R = A.
!
!                          WHEN MBDCND = 1 OR 2, 
!                            BDA(J) = U(A, Z(J)) ,       J=1, 2, ..., N.
!
!                          WHEN MBDCND = 3 OR 4, 
!                            BDA(J) = (D/DR)U(A, Z(J)) ,   J=1, 2, ..., N.
!
!                          WHEN MBDCND = 5 OR 6, BDA IS A DUMMY
!                          VARIABLE.
!
!                        BDB
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT R = B.
!
!                          WHEN MBDCND = 1, 4, OR 5, 
!                            BDB(J) = U(B, Z(J)) ,        J=1, 2, ..., N.
!
!                          WHEN MBDCND = 2, 3, OR 6, 
!                            BDB(J) = (D/DR)U(B, Z(J)) ,   J=1, 2, ..., N.
!
!                        C, D
!                          THE RANGE OF Z, I.E. C .LE. Z .LE. D.
!                          C MUST BE LESS THAN D.
!
!                        N
!                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
!                          (C, D).  THE UNKNOWNS IN THE Z-DIRECTION
!                          ARE GIVEN BY Z(J) = C + (J-0.5)DZ, 
!                          J=1, 2, ..., N, WHERE DZ = (D-C)/N.
!                          N MUST BE GREATER THAN 2.
!
!                        NBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT Z = C  AND Z = D.
!
!                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
!                               U(I, J) = U(I, N+J).
!
!                          = 1  IF THE SOLUTION IS SPECIFIED AT Z = C
!                               AND Z = D.
!
!                          = 2  IF THE SOLUTION IS SPECIFIED AT Z = C
!                               AND THE DERIVATIVE OF THE SOLUTION WITH
!                               RESPECT TO Z IS SPECIFIED AT Z = D.
!
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION WITH
!                               RESPECT TO Z IS SPECIFIED AT Z = C
!                               AND Z = D.
!
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION WITH
!                               RESPECT TO Z IS SPECIFIED AT Z = C AND
!                               THE SOLUTION IS SPECIFIED AT Z = D.
!
!                        BDC
!                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT Z = C.
!
!                          WHEN NBDCND = 1 OR 2, 
!                            BDC(I) = U(R(I), C) ,        I=1, 2, ..., M.
!
!                          WHEN NBDCND = 3 OR 4, 
!                            BDC(I) = (D/DZ)U(R(I), C),    I=1, 2, ..., M.
!
!                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
!
!                        BDD
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT Z = D.
!
!                          WHEN NBDCND = 1 OR 4, 
!                            BDD(I) = U(R(I), D) ,       I=1, 2, ..., M.
!
!                          WHEN NBDCND = 2 OR 3, 
!                            BDD(I) = (D/DZ)U(R(I), D) ,   I=1, 2, ..., M.
!
!                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
!
!                        ELMBDA
!                          THE CONSTANT LAMBDA IN THE MODIFIED
!                          HELMHOLTZ EQUATION.  IF LAMBDA IS GREATER
!                          THAN 0, A SOLUTION MAY NOT EXIST.
!                          HOWEVER, hstcyl WILL ATTEMPT TO FIND A
!                          SOLUTION.  LAMBDA MUST BE ZERO WHEN
!                          MBDCND = 5 OR 6.
!
!                        F
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
!                          THE VALUES OF THE RIGHT SIDE OF THE
!                          MODIFIED HELMHOLTZ EQUATION.
!                          FOR I=1, 2, ..., M   AND J=1, 2, ..., N
!                            F(I, J) = F(R(I), Z(J)) .
!                          F MUST BE DIMENSIONED AT LEAST M X N.
!
!                        IDIMF
!                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
!                          F AS IT APPEARS IN THE PROGRAM CALLING
!                          hstcyl.  THIS PARAMETER IS USED TO SPECIFY
!                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
!                          BE AT LEAST M.
!
! ON OUTPUT
!
!                        F
!                          CONTAINS THE SOLUTION U(I, J) OF THE FINITE
!                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
!                          (R(I), Z(J)) FOR  I=1, 2, ..., M, J=1, 2, ..., N.
!
!                        PERTRB
!                          IF A COMBINATION OF PERIODIC, DERIVATIVE, 
!                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
!                          SPECIFIED FOR A POISSON EQUATION
!                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
!                          PERTRB IS A CONSTANT, CALCULATED AND
!                          SUBTRACTED FROM F, WHICH ENSURES THAT A
!                          SOLUTION EXISTS.  hstcyl THEN COMPUTES
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
!                          PARAMETERS. EXCEPT TO NUMBERS 0 AND 11, 
!                          A SOLUTION IS NOT ATTEMPTED.
!
!                          =  0  NO ERROR
!
!                          =  1  A .LT. 0
!
!                          =  2  A .GE. B
!
!                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6
!
!                          =  4  C .GE. D
!
!                          =  5  N .LE. 2
!
!                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4
!
!                          =  7  A = 0 AND MBDCND = 1, 2, 3, OR 4
!
!                          =  8  A .GT. 0 AND MBDCND .GE. 5
!
!                          =  9  M .LE. 2
!
!                          = 10  IDIMF .LT. M
!
!                          = 11  LAMBDA .GT. 0
!
!                          = 12  A=0, MBDCND .GE. 5, ELMBDA .NE. 0
!
!                          SINCE THIS IS THE ONLY MEANS OF INDICATING
!                          A POSSIBLY INCORRECT CALL TO hstcyl, THE
!                          USER SHOULD TEST ierror AFTER THE CALL.
!
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if N, M are too large
!                               for your computer)
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED LIBRARY       fish.f, comf.f, genbun.f, gnbnaux.f, poistg.f
! FILES
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
!                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
!                        IN JANUARY 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated work space.
!
! PORTABILITY            FORTRAN 90
!
! ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
!                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
!                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR AND
!                        CALLS EITHER POISTG OR genbun WHICH SOLVES THE
!                        LINEAR SYSTEM OF EQUATIONS.
!
! TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
!                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
!
! ACCURACY               THE SOLUTION PROCESS RESULTS IN A LOSS
!                        OF NO MORE THAN FOUR SIGNIFICANT DIGITS
!                        FOR N AND M AS LARGE AS 64.
!                        MORE DETAILED INFORMATION ABOUT ACCURACY
!                        CAN BE FOUND IN THE DOCUMENTATION FOR
!                        SUBROUTINE POISTG WHICH IS THE ROUTINE THAT
!                        ACTUALLY SOLVES THE FINITE DIFFERENCE
!                        EQUATIONS.
!
! REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD FOR
!                        THE SOLUTION OF POISSON'S EQUATION WITH NEUMANN
!                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
!                        ARBITRARY SIZE, " J. COMP. PHYS. 20(1976), 
!                        PP. 171-182.
!***********************************************************************
!
module module_hstcyl

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
    public :: hstcyl
    public :: test_hstcyl

contains
     !
     !*****************************************************************************************
     !
    subroutine test_hstcyl()
        !     file thstcyl.f
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
        !
        !     PROGRAM TO ILLUSTRATE THE USE OF hstcyl TO SOLVE THE EQUATION
        !
        !    (1/R)(D/DR)(R*DU/DR) + (D/DZ)(DU/DZ) = (2*R*Z)**2*(4*Z**2 + 3*R**2)
        !
        !     ON THE RECTANGLE 0 .LT. R .LT. 1 , 0 .LT. Z .LT. 1 WITH THE
        !     BOUNDARY CONDITIONS
        !
        !     (DU/DR)(1, Z) = 4*Z**2  FOR  0 .LE. Z .LE. 1
        !
        !     AND
        !
        !     (DU/DZ)(R, 0) = 0 AND (DU/DZ)(R, 1) = 4*R**2  FOR  0 .LE. R .LE. 1 .
        !
        !     THE SOLUTION TO THIS PROBLEM IS NOT UNIQUE.  IT IS A
        !     ONE-PARAMETER FAMILY OF SOLUTIONS GIVEN BY
        !
        !            U(R, Z) = (R*Z)**4 + ARBITRARY CONSTANT .
        !
        !     THE R-INTERVAL WILL CONTAIN 50 UNKNOWNS AND THE Z-INTERVAL WILL
        !     CONTAIN 52 UNKNOWNS.
        !
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip)                 :: idimf, m, mbdcnd, n, nbdcnd, i, j, ierror
        real (wp), dimension(51, 52) :: f
        real (wp), dimension(52)     :: bda, bdb
        real (wp), dimension(50)     :: bdc, bdd, r
        real (wp), dimension(52)     :: z
        real (wp)                    :: a, b, c, d, elmbda, pertrb, x, discretization_error
        !-----------------------------------------------

        !     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
        !
        idimf = 51
        a = 0.0_wp
        b = 1.0_wp
        m = 50
        mbdcnd = 6
        c = 0.0_wp
        d = 1.0_wp
        n = 52
        nbdcnd = 3
        elmbda = 0.0_wp
        !
        !     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
        !     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
        !
        do i = 1, m
            r(i) = (real(i, kind=wp) - 0.5_wp)/50
        end do

        do j = 1, n
            z(j) = (real(j, kind=wp) - 0.5_wp)/52
        end do
        !
        !     Generate boundary data.
        !
        bdb(:n) = 4.0_wp * ( z(:n)**4 )
        !
        !     Generate boundary data.
        !
        bdc(:m) = 0.0_wp
        bdd(:m) = 4.0_wp * ( r(:m)**4 )
        !
        !     BDA IS A DUMMY VARIABLE.
        !
        !     GENERATE RIGHT SIDE OF EQUATION.
        !
        do i = 1, m
            f(i, :n) = 4.0_wp*(r(i)**2) * (z(:n)**2) * (4.0_wp * (z(:n)**2) &
                + 3.0_wp * (r(i)**2) )
        end do

        call hstcyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
            elmbda, f, idimf, pertrb, ierror)
        !
        !     COMPUTE DISCRETIZATION ERROR BY MINIMIZING OVER ALL A THE FUNCTION
        !     NORM(F(I, J) - A*1 - U(R(I), Z(J))).  THE EXACT SOLUTION IS
        !                U(R, Z) = (R*Z)**4 + ARBITRARY CONSTANT.
        !
        x = 0.0_wp
        do i = 1, m
            x = x + sum(f(i, :n)-(r(i)*z(:n))**4)
        end do
        x = x/(m*n)
        f(:m, :n) = f(:m, :n) - x
        discretization_error = 0.0_wp
        do i = 1, m
            do j = 1, n
                x = abs(f(i, j)-(r(i)*z(j))**4)
                discretization_error = max(x, discretization_error)
            end do
        end do
        !     Print earlier output from platforms with 32 and 64 bit floating point
        !     arithemtic followed by the output from this computer
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     hstcyl *** TEST RUN *** '
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') '     ierror = 0,  PERTRB =-4.4311E-4'
        write( stdout, '(A)') '     discretization error = 7.5280E-5 '
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)') '     ierror =', ierror, ' PERTRB = ', pertrb
        write( stdout, '(A,1pe15.6)') '     discretization error = ', discretization_error

    end subroutine test_hstcyl
    !
    !*****************************************************************************************
    !
    subroutine hstcyl( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror )
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
        integer                  :: irwk
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------

        ! initialize error flag
        ierror = 0

        ! Check input arguments
        if (a < 0.0_wp) ierror = 1
        if (a >= b) ierror = 2
        if (mbdcnd <= 0 .or. mbdcnd >= 7) ierror = 3
        if (c >= d) ierror = 4
        if (n <= 2) ierror = 5
        if (nbdcnd < 0 .or. nbdcnd >= 5) ierror = 6
        if (a == 0.0_wp .and. mbdcnd /= 5 .and. mbdcnd /= 6) ierror = 7
        if (a > 0.0_wp .and. mbdcnd >= 5) ierror = 8
        if (idimf < m) ierror = 10
        if (m <= 2) ierror = 9
        if (a == 0.0_wp .and. mbdcnd >= 5 .and. elmbda /= 0.0_wp) ierror = 12
        if (ierror /= 0) return

        ! calculate required real work space size
        call workspace%get_genbun_workspace_dimensions( n, m, irwk )
        irwk = irwk + 3 * m

        ! Allocate real workspace array
        associate( icwk => 0)
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! if workspace allocation was succesful
        if (ierror == 20) return

        ! solve system
        associate( rew => workspace%rew )
            call hstcyll(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror, rew )
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hstcyl


    subroutine hstcyll( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, w )
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in) :: m
        integer (ip), intent (in) :: mbdcnd
        integer (ip), intent (in) :: n
        integer (ip), intent (in) :: nbdcnd
        integer (ip), intent (in) :: idimf
        integer (ip), intent (out) :: ierror
        real (wp),    intent (in) :: a
        real (wp),    intent (in) :: b
        real (wp),    intent (in) :: c
        real (wp),    intent (in) :: d
        real (wp),    intent (in) :: elmbda
        real (wp),    intent (out) :: pertrb
        real (wp),    intent (in) :: bda(*)
        real (wp),    intent (in) :: bdb(*)
        real (wp),    intent (in) :: bdc(*)
        real (wp),    intent (in) :: bdd(*)
        real (wp),    intent (in out) :: f(idimf,*)
        real (wp),    intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: np, iwb, iwc, iwr, i, j, k, lp, ierr1
        real (wp)    :: deltar, dlrsq, deltht, dlthsq, a1
        !-----------------------------------------------

        deltar = (b - a)/m
        dlrsq = deltar**2
        deltht = (d - c)/n
        dlthsq = deltht**2
        np = nbdcnd + 1
        !
        !     Define a, b, c coefficients in w-array.
        !
        iwb = m
        iwc = iwb + m
        iwr = iwc + m
        do i = 1, m
            j = iwr + i
            w(j) = a + (real(i, kind=wp) - 0.5_wp)*deltar
            w(i) = (a + real(i - 1, kind=wp)*deltar)/(dlrsq*w(j))
            k = iwc + i
            w(k) = (a + real(i, kind=wp)*deltar)/(dlrsq*w(j))
            k = iwb + i
            w(k) = elmbda/w(j)**2 - 2.0_wp/dlrsq
        end do
        !
        !     Enter boundary data for r-boundaries.
        !
        go to (102, 102, 104, 104, 106, 106) mbdcnd
102 continue
    a1 = 2.*w(1)
    w(iwb+1) = w(iwb+1) - w(1)
    f(1, :n) = f(1, :n) - a1*bda(:n)
    go to 106
104 continue
    a1 = deltar*w(1)
    w(iwb+1) = w(iwb+1) + w(1)
    f(1, :n) = f(1, :n) + a1 * bda(:n)
106 continue
    go to (107, 109, 109, 107, 107, 109) mbdcnd
107 continue
    w(iwc) = w(iwc) - w(iwr)
    a1 = 2.0_wp * w(iwr)
    f(m, :n) = f(m, :n) - a1 * bdb(:n)
    go to 111
109 continue
    w(iwc) = w(iwc) + w(iwr)
    a1 = deltar * w(iwr)
    f(m, :n) = f(m, :n) - a1 * bdb(:n)
!
!     Enter boundary data for theta-boundaries.
!
111 continue
    a1 = 2.0_wp/dlthsq
    go to (121, 112, 112, 114, 114) np
112 continue
    f(:m, 1) = f(:m, 1) - a1*bdc(:m)
    go to 116
114 continue
    a1 = 1.0_wp/deltht
    f(:m, 1) = f(:m, 1) + a1*bdc(:m)
116 continue
    a1 = 2.0_wp/dlthsq
    go to (121, 117, 119, 119, 117) np
117 continue
    f(:m, n) = f(:m, n) - a1*bdd(:m)
    go to 121
119 continue
    a1 = 1./deltht
    f(:m, n) = f(:m, n) - a1*bdd(:m)
121 continue
    pertrb = 0.0_wp
    if (elmbda >= 0.0_wp ) then
        if (elmbda /= 0.0_wp ) then
            ierror = 11
        else
            go to (130, 130, 124, 130, 130, 124) mbdcnd
124     continue
        go to (125, 130, 130, 125, 130) np
125 continue
    do i = 1, m
        a1 = 0.0_wp
        a1 = sum(f(i, :n))
        j = iwr + i
        pertrb = pertrb + a1 * w(j)
    end do
    pertrb = pertrb/(real( m * n, kind=wp) * 0.5_wp * (a + b))
    f(:m, :n) = f(:m, :n) - pertrb
end if
end if
130 continue
    w(:m) = w(:m) * dlthsq
    w(iwc+1:m+iwc) = w(iwc+1:m+iwc) * dlthsq
    w(iwb+1:m+iwb) = w(iwb+1:m+iwb) * dlthsq
    f(:m, :n) = f(:m, :n)*dlthsq
    lp = nbdcnd
    w(1) = 0.0_wp
    w(iwr) = 0.0_wp
    !
    !     SOLVE THE SYSTEM OF EQUATIONS.
    !
    ierr1 = 0
    if (nbdcnd /= 0) then
        call poistgg (lp, n, 1, m, w, w(iwb+1), w(iwc+1), idimf, f, &
            ierr1, w(iwr+1))
    else
        call genbunn (lp, n, 1, m, w, w(iwb+1), w(iwc+1), idimf, f, &
            ierr1, w(iwr+1))
    end if

end subroutine hstcyll
    !
    !*****************************************************************************************
    !
end module module_hstcyl
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
