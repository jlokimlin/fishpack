module module_hstplr

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
    public :: hstplr


contains


    subroutine hstplr( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror )
        !
        !     file hstplr.f
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
        !     SUBROUTINE hstplr (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD,
        !    +                   ELMBDA, F, IDIMF, PERTRB, ierror)
        !
        ! DIMENSION OF           BDA(N), BDB(N), BDC(M), BDD(M), F(IDIMF, N)
        ! ARGUMENTS
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
        !                        DIFFERENCE APPROXIMATION ON A STAGGERED
        !                        GRID TO THE HELMHOLTZ EQUATION IN POLAR
        !                        COORDINATES.  THE EQUATION IS
        !
        !                           (1/R)(D/DR)(R(DU/DR)) +
        !                           (1/R**2)(D/DTHETA)(DU/DTHETA) +
        !                           LAMBDA*U = F(R, THETA)
        !
        ! USAGE                  CALL hstplr (A, B, M, MBDCND, BDA, BDB, C, D, N,
        !                                     NBDCND, BDC, BDD, ELMBDA, F,
        !                                     IDIMF, PERTRB, ierror)
        !
        ! ARGUMENTS
        ! ON INPUT               A, B
        !
        !                          THE RANGE OF R, I.E. A .LE. R .LE. B.
        !                          A MUST BE LESS THAN B AND A MUST BE
        !                          NON-NEGATIVE.
        !
        !                        M
        !                          THE NUMBER OF GRID POINTS IN THE INTERVAL
        !                          (A, B).  THE GRID POINTS IN THE R-DIRECTION
        !                          ARE GIVEN BY R(I) = A + (I-0.5)DR FOR
        !                          I=1, 2, ..., M WHERE DR =(B-A)/M.
        !                          M MUST BE GREATER THAN 2.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT R = A AND R = B.
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT R = A
        !                               AND R = B.
        !
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT R = A
        !                               AND THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO R IS SPECIFIED AT R = B.
        !                               (SEE NOTE 1 BELOW)
        !
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO R IS SPECIFIED AT
        !                               R = A (SEE NOTE 2 BELOW) AND R = B.
        !
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO R IS SPECIFIED AT
        !                               SPECIFIED AT R = A (SEE NOTE 2 BELOW)
        !                               AND THE SOLUTION IS SPECIFIED AT R = B.
        !
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
        !                          NOTE 1:
        !                          IF A = 0, MBDCND = 2, AND NBDCND = 0 OR 3,
        !                          THE SYSTEM OF EQUATIONS TO BE SOLVED IS
        !                          SINGULAR.  THE UNIQUE SOLUTION IS
        !                          IS DETERMINED BY EXTRAPOLATION TO THE
        !                          SPECIFICATION OF U(0, THETA(1)).
        !                          BUT IN THIS CASE THE RIGHT SIDE OF THE
        !                          SYSTEM WILL BE PERTURBED BY THE CONSTANT
        !                          PERTRB.
        !
        !                          NOTE 2:
        !                          IF A = 0, DO NOT USE MBDCND = 3 OR 4,
        !                          BUT INSTEAD USE MBDCND = 1, 2, 5, OR 6.
        !
        !                        BDA
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
        !                          SPECIFIES THE BOUNDARY VALUES (IF ANY) OF
        !                          THE SOLUTION AT R = A.
        !
        !                          WHEN MBDCND = 1 OR 2,
        !                            BDA(J) = U(A, THETA(J)) ,     J=1, 2, ..., N.
        !
        !                          WHEN MBDCND = 3 OR 4,
        !                            BDA(J) = (D/DR)U(A, THETA(J)) ,
        !                            J=1, 2, ..., N.
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
        !                            BDB(J) = U(B, THETA(J)) ,     J=1, 2, ..., N.
        !
        !                          WHEN MBDCND = 2, 3, OR 6,
        !                            BDB(J) = (D/DR)U(B, THETA(J)) ,
        !                            J=1, 2, ..., N.
        !
        !                        C, D
        !                          THE RANGE OF THETA, I.E. C .LE. THETA .LE. D.
        !                          C MUST BE LESS THAN D.
        !
        !                        N
        !                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
        !                          (C, D).  THE UNKNOWNS IN THE THETA-
        !                          DIRECTION ARE GIVEN BY THETA(J) = C +
        !                          (J-0.5)DT,   J=1, 2, ..., N, WHERE
        !                          DT = (D-C)/N.  N MUST BE GREATER THAN 2.
        !
        !                        NBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT THETA = C  AND THETA = D.
        !
        !                          = 0  IF THE SOLUTION IS PERIODIC IN THETA,
        !                               I.E. U(I, J) = U(I, N+J).
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = C AND THETA = D
        !                               (SEE NOTE BELOW).
        !
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               THETA = C AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO THETA IS
        !                               SPECIFIED AT THETA = D
        !                               (SEE NOTE BELOW).
        !
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               AT THETA = C AND THETA = D.
        !
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO THETA IS SPECIFIED
        !                               AT THETA = C AND THE SOLUTION IS
        !                               SPECIFIED AT THETA = D
        !                               (SEE NOTE BELOW).
        !
        !                          NOTE:
        !                          WHEN NBDCND = 1, 2, OR 4, DO NOT USE
        !                          MBDCND = 5 OR 6 (THE FORMER INDICATES THAT
        !                          THE SOLUTION IS SPECIFIED AT R =  0; THE
        !                          LATTER INDICATES THE SOLUTION IS UNSPECIFIED
        !                          AT R = 0).  USE INSTEAD MBDCND = 1 OR 2.
        !
        !                        BDC
        !                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
        !                          SPECIFIES THE BOUNDARY VALUES OF THE
        !                          SOLUTION AT THETA = C.
        !
        !                          WHEN NBDCND = 1 OR 2,
        !                            BDC(I) = U(R(I), C) ,        I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 3 OR 4,
        !                            BDC(I) = (D/DTHETA)U(R(I), C),
        !                            I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
        !
        !                        BDD
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
        !                          SPECIFIES THE BOUNDARY VALUES OF THE
        !                          SOLUTION AT THETA = D.
        !
        !                          WHEN NBDCND = 1 OR 4,
        !                            BDD(I) = U(R(I), D) ,         I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 2 OR 3,
        !                            BDD(I) =(D/DTHETA)U(R(I), D), I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
        !
        !                        ELMBDA
        !                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
        !                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
        !                          A SOLUTION MAY NOT EXIST.  HOWEVER, hstplr
        !                          WILL ATTEMPT TO FIND A SOLUTION.
        !
        !                        F
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUES OF THE RIGHT SIDE OF THE HELMHOLTZ
        !                          EQUATION.
        !
        !                          FOR I=1, 2, ..., M AND J=1, 2, ..., N
        !                            F(I, J) = F(R(I), THETA(J)) .
        !
        !                          F MUST BE DIMENSIONED AT LEAST M X N.
        !
        !                        IDIMF
        !                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
        !                          F AS IT APPEARS IN THE PROGRAM CALLING
        !                          hstplr.  THIS PARAMETER IS USED TO SPECIFY
        !                          THE VARIABLE DIMENSION OF F.
        !                          IDIMF MUST BE AT LEAST M.
        !
        !
        ! ON OUTPUT
        !
        !                        F
        !                          CONTAINS THE SOLUTION U(I, J) OF THE FINITE
        !                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
        !                          (R(I), THETA(J)) FOR I=1, 2, ..., M,
        !                          J=1, 2, ..., N.
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
        !                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
        !                          SPECIFIED FOR A POISSON EQUATION
        !                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
        !                          PERTRB IS A CONSTANT CALCULATED AND
        !                          SUBTRACTED FROM F, WHICH ENSURES THAT A
        !                          SOLUTION EXISTS.  hstplr THEN COMPUTES THIS
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
        !                          =  7  A = 0 AND MBDCND = 3 OR 4
        !
        !                          =  8  A .GT. 0 AND MBDCND .GE. 5
        !
        !                          =  9  MBDCND .GE. 5 AND NBDCND .NE. 0 OR 3
        !
        !                          = 10  IDIMF .LT. M
        !
        !                          = 11  LAMBDA .GT. 0
        !
        !                          = 12  M .LE. 2
        !
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
        !
        !                          SINCE THIS IS THE ONLY MEANS OF INDICATING
        !                          A POSSIBLY INCORRECT CALL TO hstplr, THE
        !                          USER SHOULD TEST ierror AFTER THE CALL.
        !
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
        ! PORTABILITY            FORTRAN 90
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
        integer                  :: irwk
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------

        ierror = 0
        if (a < 0.) ierror = 1
        if (a >= b) ierror = 2
        if (mbdcnd<=0 .or. mbdcnd>=7) ierror = 3
        if (c >= d) ierror = 4
        if (n <= 2) ierror = 5
        if (nbdcnd<0 .or. nbdcnd>=5) ierror = 6
        if (a==0. .and. (mbdcnd==3 .or. mbdcnd==4)) ierror = 7
        if (a>0. .and. mbdcnd>=5) ierror = 8
        if (mbdcnd>=5 .and. nbdcnd/=0 .and. nbdcnd/=3) ierror = 9
        if (idimf < m) ierror = 10
        if (m <= 2) ierror = 12
        if (ierror /= 0) return

        ! compute required real workspace size
        call workspace%get_genbun_workspace_dimensions( n, m, irwk )
        irwk = irwk + 3 * m

        ! allocate real workspace array
        associate( icwk => 0 )
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! Check if allocation was successful
        if (ierror == 20) return

        associate( rew => workspace%real_workspace )
            call hstplrr( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror, rew)
        end associate

        ! release memory
        call workspace%destroy()

    end subroutine hstplr


    subroutine hstplrr( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, w )
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: mbdcnd
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: nbdcnd
        integer (ip), intent (in)     :: idimf
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in)     :: a
        real (wp),    intent (in)     :: b
        real (wp),    intent (in)     :: c
        real (wp),    intent (in)     :: d
        real (wp),    intent (in)     :: elmbda
        real (wp),    intent (out)    :: pertrb
        real (wp),    intent (in)     :: bda(*)
        real (wp),    intent (in)     :: bdb(*)
        real (wp),    intent (in)     :: bdc(*)
        real (wp),    intent (in)     :: bdd(*)
        real (wp),    intent (in out) :: f(idimf,*)
        real (wp),    intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer :: np, isw, mb, iwb, iwc, iwr, i, j, k, lp, ierr1
        real :: deltar, dlrsq, deltht, dlthsq, a1, a2
        !-----------------------------------------------
        deltar = (b - a)/m
        dlrsq = deltar**2
        deltht = (d - c)/n
        dlthsq = deltht**2
        np = nbdcnd + 1
        isw = 1
        mb = mbdcnd
        if (a==0. .and. mbdcnd==2) mb = 6
        !
        !     define a, b, c coefficients in w-array.
        !
        iwb = m
        iwc = iwb + m
        iwr = iwc + m
        do i = 1, m
            j = iwr + i
            w(j) = a + (real(i) - 0.5)*deltar
            w(i) = (a + real(i - 1)*deltar)/dlrsq
            k = iwc + i
            w(k) = (a + real(i)*deltar)/dlrsq
            k = iwb + i
            w(k) = (elmbda - 2./dlrsq)*w(j)
        end do
        do i = 1, m
            j = iwr + i
            a1 = w(j)
            f(i, :n) = a1*f(i, :n)
        end do
        !
        !     enter boundary data for r-boundaries.
        !
        go to (104, 104, 106, 106, 108, 108) mb
104 continue
    a1 = 2.*w(1)
    w(iwb+1) = w(iwb+1) - w(1)
    f(1, :n) = f(1, :n) - a1*bda(:n)
    go to 108
106 continue
    a1 = deltar*w(1)
    w(iwb+1) = w(iwb+1) + w(1)
    f(1, :n) = f(1, :n) + a1*bda(:n)
108 continue
    go to (109, 111, 111, 109, 109, 111) mb
109 continue
    a1 = 2.*w(iwr)
    w(iwc) = w(iwc) - w(iwr)
    f(m, :n) = f(m, :n) - a1*bdb(:n)
    go to 113
111 continue
    a1 = deltar*w(iwr)
    w(iwc) = w(iwc) + w(iwr)
    f(m, :n) = f(m, :n) - a1*bdb(:n)
!
!     enter boundary data for theta-boundaries.
!
113 continue
    a1 = 2./dlthsq
    go to (123, 114, 114, 116, 116) np
114 continue
    f(:m, 1) = f(:m, 1) - a1*bdc(:m)/w(iwr+1:m+iwr)
    go to 118
116 continue
    a1 = 1./deltht
    f(:m, 1) = f(:m, 1) + a1*bdc(:m)/w(iwr+1:m+iwr)
118 continue
    a1 = 2./dlthsq
    go to (123, 119, 121, 121, 119) np
119 continue
    f(:m, n) = f(:m, n) - a1*bdd(:m)/w(iwr+1:m+iwr)
    go to 123
121 continue
    a1 = 1./deltht
    f(:m, n) = f(:m, n) - a1*bdd(:m)/w(iwr+1:m+iwr)
123 continue
    pertrb = 0.
    if (elmbda >= 0.) then
        if (elmbda /= 0.) then
            ierror = 11
        else
            go to (133, 133, 126, 133, 133, 126) mb
126     continue
        go to (127, 133, 133, 127, 133) np
127 continue
    isw = 2
    do j = 1, n
        pertrb = pertrb + sum(f(:m, j))
    end do
    pertrb = pertrb/(real(m*n)*0.5*(a + b))
    do i = 1, m
        j = iwr + i
        a1 = pertrb*w(j)
        f(i, :n) = f(i, :n) - a1
    end do
    a2 = sum(f(1, :n))
    a2 = a2/w(iwr+1)
end if
end if
133 continue
    do i = 1, m
        j = iwr + i
        a1 = dlthsq*w(j)
        w(i) = a1*w(i)
        j = iwc + i
        w(j) = a1*w(j)
        j = iwb + i
        w(j) = a1*w(j)
        f(i, :n) = a1*f(i, :n)
    end do
    lp = nbdcnd
    w(1) = 0.
    w(iwr) = 0.
    !
    !     to solve the system of equations.
    !
    ierr1 = 0
    if (lp /= 0) then
        call poistgg(lp, n, 1, m, w, w(iwb+1), w(iwc+1), idimf, f, &
            ierr1, w(iwr+1))
    else
        call genbunn(lp, n, 1, m, w, w(iwb+1), w(iwc+1), idimf, f, &
            ierr1, w(iwr+1))
    end if
    if (.not.(a/=0. .or. mbdcnd/=2 .or. isw/=2)) then
        a1 = sum(f(1, :n))
        a1 = (a1 - dlrsq*a2/16.)/real(n)
        if (nbdcnd == 3) a1 = a1 + (bdd(1)-bdc(1))/(d - c)
        a1 = bda(1) - a1
        f(:m, :n) = f(:m, :n) + a1
    end if

end subroutine hstplrr


end module module_hstplr
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
