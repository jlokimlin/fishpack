module module_hstcrt

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
    public :: hstcrt


contains


    subroutine hstcrt( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror )
        !
        !     file hstcrt.f
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
        !     SUBROUTINE hstcrt (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD,
        !    +                   ELMBDA, F, IDIMF, PERTRB, ierror)
        !
        ! DIMENSION OF           BDA(N), BDB(N), BDC(M), BDD(M), F(IDIMF, N)
        ! ARGUMENTS
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                 SOLVES THE STANDARD FIVE-POINT FINITE
        !                         DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
        !                         EQUATION
        !                           (D/DX)(DU/DX) + (D/DY)(DU/DY) + LAMBDA*U
        !                           = F(X, Y)
        !                         ON A STAGGERED GRID IN CARTESIAN COORDINATES.
        !
        ! USAGE                   CALL hstcrt (A, B, M, MBDCND, BDA, BDB, C, D
        !                                      N, NBDCND, BDC, BDD, ELMBDA,
        !                                      F, IDIMF, PERTRB, ierror)
        !
        ! ARGUMENTS
        ! ON INPUT
        !
        !                        A, B
        !                          THE RANGE OF X, I.E. A .LE. X .LE. B.
        !                          A MUST BE LESS THAN B.
        !
        !                        M
        !                          THE NUMBER OF GRID POINTS IN THE
        !                          INTERVAL (A, B).  THE GRID POINTS
        !                          IN THE X-DIRECTION ARE GIVEN BY
        !                          X(I) = A + (I-0.5)DX FOR I=1, 2, ..., M
        !                          WHERE DX =(B-A)/M.  M MUST BE GREATER
        !                          THAN 2.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT X = A AND X = B.
        !
        !                          = 0  IF THE SOLUTION IS PERIODIC IN X,
        !                               U(M+I, J) = U(I, J).
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               X = A AND X = B.
        !
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               X = A AND THE DERIVATIVE
        !                               OF THE SOLUTION WITH RESPECT TO X
        !                               IS SPECIFIED AT X = B.
        !
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO X IS SPECIFIED
        !                               AT X = A  AND X = B.
        !
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO X IS SPECIFIED
        !                               AT X = A  AND THE SOLUTION IS
        !                               SPECIFIED AT X = B.
        !
        !                        BDA
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N
        !                          THAT SPECIFIES THE BOUNDARY VALUES
        !                          (IF ANY) OF THE SOLUTION AT X = A.
        !
        !                          WHEN MBDCND = 1 OR 2,
        !                            BDA(J) = U(A, Y(J)) ,         J=1, 2, ..., N.
        !
        !                          WHEN MBDCND = 3 OR 4,
        !                            BDA(J) = (D/DX)U(A, Y(J)) ,   J=1, 2, ..., N.
        !
        !                        BDB
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N
        !                          THAT SPECIFIES THE BOUNDARY VALUES
        !                          OF THE SOLUTION AT X = B.
        !
        !                          WHEN MBDCND = 1 OR 4
        !                            BDB(J) = U(B, Y(J)) ,        J=1, 2, ..., N.
        !
        !                          WHEN MBDCND = 2 OR 3
        !                            BDB(J) = (D/DX)U(B, Y(J)) ,  J=1, 2, ..., N.
        !
        !                        C, D
        !                          THE RANGE OF Y, I.E. C .LE. Y .LE. D.
        !                          C MUST BE LESS THAN D.
        !
        !
        !                        N
        !                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
        !                          (C, D).  THE UNKNOWNS IN THE Y-DIRECTION
        !                          ARE GIVEN BY Y(J) = C + (J-0.5)DY,
        !                          J=1, 2, ..., N, WHERE DY = (D-C)/N.
        !                          N MUST BE GREATER THAN 2.
        !
        !                        NBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT Y = C   AND Y = D.
        !
        !
        !                          = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.
        !                               U(I, J) = U(I, N+J).
        !
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT Y = C
        !                               AND Y = D.
        !
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT Y = C
        !                               AND THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO Y IS SPECIFIED AT
        !                               Y = D.
        !
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO Y IS SPECIFIED AT
        !                               Y = C AND Y = D.
        !
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO Y IS SPECIFIED AT
        !                               Y = C AND THE SOLUTION IS SPECIFIED
        !                               AT Y = D.
        !
        !                        BDC
        !                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
        !                          SPECIFIES THE BOUNDARY VALUES OF THE
        !                          SOLUTION AT Y = C.
        !
        !                          WHEN NBDCND = 1 OR 2,
        !                            BDC(I) = U(X(I), C) ,        I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 3 OR 4,
        !                            BDC(I) = (D/DY)U(X(I), C),   I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
        !
        !                        BDD
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
        !                          SPECIFIES THE BOUNDARY VALUES OF THE
        !                          SOLUTION AT Y = D.
        !
        !                          WHEN NBDCND = 1 OR 4,
        !                            BDD(I) = U(X(I), D) ,        I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 2 OR 3,
        !                            BDD(I) = (D/DY)U(X(I), D) ,  I=1, 2, ..., M.
        !
        !                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
        !
        !                        ELMBDA
        !                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
        !                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
        !                          A SOLUTION MAY NOT EXIST. HOWEVER,
        !                          hstcrt WILL  ATTEMPT TO FIND A SOLUTION.
        !
        !                        F
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
        !                          THE VALUES OF THE RIGHT SIDE OF THE
        !                          HELMHOLTZ EQUATION.  FOR I=1, 2, ..., M
        !                          AND J=1, 2, ..., N
        !
        !                            F(I, J) = F(X(I), Y(J)) .
        !
        !                          F MUST BE DIMENSIONED AT LEAST M X N.
        !
        !                        IDIMF
        !                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
        !                          F AS IT APPEARS IN THE PROGRAM CALLING
        !                          hstcrt.  THIS PARAMETER IS USED TO SPECIFY
        !                          THE VARIABLE DIMENSION OF F.
        !                          IDIMF MUST BE AT LEAST M.
        !
        !
        ! ON OUTPUT              F
        !                          CONTAINS THE SOLUTION U(I, J) OF THE FINITE
        !                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
        !                          (X(I), Y(J)) FOR  I=1, 2, ..., M, J=1, 2, ..., N.
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
        !                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
        !                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
        !                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
        !                          CALCULATED AND SUBTRACTED FROM F, WHICH
        !                          ENSURES THAT A SOLUTION EXISTS.  hstcrt
        !                          THEN COMPUTES THIS SOLUTION, WHICH IS A
        !                          LEAST SQUARES SOLUTION TO THE ORIGINAL
        !                          APPROXIMATION.  THIS SOLUTION PLUS ANY
        !                          CONSTANT IS ALSO A SOLUTION; HENCE, THE
        !                          SOLUTION IS NOT UNIQUE.  THE VALUE OF
        !                          PERTRB SHOULD BE SMALL COMPARED TO THE
        !                          RIGHT SIDE F.  OTHERWISE, A SOLUTION IS
        !                          OBTAINED TO AN ESSENTIALLY DIFFERENT PROBLEM.
        !                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
        !                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
        !                          OBTAINED.
        !
        !                        ierror
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS.  EXCEPT TO NUMBERS 0 AND  6,
        !                          A SOLUTION IS NOT ATTEMPTED.
        !
        !                          =  0  NO ERROR
        !
        !                          =  1  A .GE. B
        !
        !                          =  2  MBDCND .LT. 0 OR MBDCND .GT. 4
        !
        !                          =  3  C .GE. D
        !
        !                          =  4  N .LE. 2
        !
        !                         =  5  NBDCND .LT. 0 OR NBDCND .GT. 4
        !
        !                         =  6  LAMBDA .GT. 0
        !
        !                         =  7  IDIMF .LT. M
        !
        !                         =  8  M .LE. 2
        !
        !                         SINCE THIS IS THE ONLY MEANS OF INDICATING
        !                         A POSSIBLY INCORRECT CALL TO hstcrt, THE
        !                         USER SHOULD TEST ierror AFTER THE CALL.
        !
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
        !
        !
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
        !                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR
        !                        AND CALLS EITHER POISTG OR genbun WHICH SOLVES
        !                        THE LINEAR SYSTEM OF EQUATIONS.
        !
        ! TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
        !                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
        !
        ! ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A
        !                        LOSS OF NO MORE THAN FOUR SIGNIFICANT DIGITS
        !                        FOR N AND M AS LARGE AS 64.  MORE DETAILED
        !                        INFORMATION ABOUT ACCURACY CAN BE FOUND IN
        !                        THE DOCUMENTATION FOR PACKAGE POISTG WHICH
        !                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
        !
        ! REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD
        !                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
        !                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
        !                        ARBITRARY SIZE, " J. COMP. PHYS. 20(1976), 
        !                        PP. 171-182.
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
        integer (ip)             :: irwk, icwk
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------
        !
        !     check for invalid parameters.
        !
        ierror = 0
        if (a >= b) ierror = 1
        if (mbdcnd<0 .or. mbdcnd>4) ierror = 2
        if (c >= d) ierror = 3
        if (n <= 2) ierror = 4
        if (nbdcnd < 0 .or. nbdcnd > 4) ierror = 5
        if (idimf < m) ierror = 7
        if (m <= 2) ierror = 8
        if (ierror /= 0) return

        ! compute required real work space
        call workspace%get_genbun_workspace_dimensions(n, m, irwk)
        irwk = irwk + 3 * m

        ! Create real workspace
        associate( icwk => 0 )
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! Solve system
        associate( rew => workspace%real_workspace )
            call hstcrtt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hstcrt


    subroutine hstcrtt( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, w )
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
        real (wp),             intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: nperod, mperod, np, mp, id2, id3, id4, j, ierr1
        real (wp)    :: deltax, twdelx, delxsq, deltay, twdely, delysq, twdysq, s, st2
        !-----------------------------------------------

        nperod = nbdcnd
        mperod = 0
        if (mbdcnd > 0) mperod = 1
        deltax = (b - a)/m
        twdelx = 1./deltax
        delxsq = 2./deltax**2
        deltay = (d - c)/n
        twdely = 1./deltay
        delysq = deltay**2
        twdysq = 2./delysq
        np = nbdcnd + 1
        mp = mbdcnd + 1
        !
        !     define the a, b, c coefficients in w-array.
        !
        id2 = m
        id3 = id2 + m
        id4 = id3 + m
        s = (deltay/deltax)**2
        st2 = 2.*s
        w(:m) = s
        w(id2+1:m+id2) = (-st2) + elmbda*delysq
        w(id3+1:m+id3) = s
        !
        !     enter boundary data for x-boundaries.
        !
        go to (111, 102, 102, 104, 104) mp
102 continue
    f(1, :n) = f(1, :n) - bda(:n)*delxsq
    w(id2+1) = w(id2+1) - w(1)
    go to 106
104 continue
    f(1, :n) = f(1, :n) + bda(:n)*twdelx
    w(id2+1) = w(id2+1) + w(1)
106 continue
    go to (111, 107, 109, 109, 107) mp
107 continue
    f(m, :n) = f(m, :n) - bdb(:n)*delxsq
    w(id3) = w(id3) - w(1)
    go to 111
109 continue
    f(m, :n) = f(m, :n) - bdb(:n)*twdelx
    w(id3) = w(id3) + w(1)
111 continue
    go to (121, 112, 112, 114, 114) np
112 continue
    f(:m, 1) = f(:m, 1) - bdc(:m)*twdysq
    go to 116
114 continue
    f(:m, 1) = f(:m, 1) + bdc(:m)*twdely
116 continue
    go to (121, 117, 119, 119, 117) np
117 continue
    f(:m, n) = f(:m, n) - bdd(:m)*twdysq
    go to 121
119 continue
    f(:m, n) = f(:m, n) - bdd(:m)*twdely
121 continue
    f(:m, :n) = f(:m, :n)*delysq
    if (mperod /= 0) then
        w(1) = 0.
        w(id4) = 0.
    end if
    pertrb = 0.
    if (elmbda >= 0.) then
        if (elmbda /= 0.) then
            ierror = 6
        else
            go to (127, 133, 133, 127, 133) mp
127     continue
        go to (128, 133, 133, 128, 133) np
    !
    !     for singular problems must adjust data to insure that a solution
    !     will exist.
    !
128 continue
    s = 0.
    do j = 1, n
        s = s + sum(f(:m, j))
    end do
    pertrb = s/real(m*n)
    f(:m, :n) = f(:m, :n) - pertrb
    pertrb = pertrb/delysq
!
!     solve the equation.
!
end if
end if
133 continue
    ierr1 = 0
    if (nperod /= 0) then
        call poistgg(nperod, n, mperod, m, w(1), w(id2+1), w(id3+1), &
            idimf, f, ierr1, w(id4+1))
    else
        call genbunn(nperod, n, mperod, m, w(1), w(id2+1), w(id3+1), &
            idimf, f, ierr1, w(id4+1))
    end if

end subroutine hstcrtt


end module module_hstcrt
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
