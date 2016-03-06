!
!     file hwscrt.f
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
!     SUBROUTINE HWSCRT (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD, 
!    +                   ELMBDA, F, IDIMF, PERTRB, IERROR)
!
! DIMENSION OF           BDA(N),      BDB(N),   BDC(M), BDD(M), 
! ARGUMENTS              F(IDIMF, N)
!
! LATEST REVISION        June 2004
!
! PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
!                        DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
!                        EQUATION IN CARTESIAN COORDINATES.  THIS
!                        EQUATION IS
!
!                          (D/DX)(DU/DX) + (D/DY)(DU/DY)
!                          + LAMBDA*U = F(X, Y).
!
! USAGE                  CALL HWSCRT (A, B, M, MBDCND, BDA, BDB, C, D, N, 
!                                     NBDCND, BDC, BDD, ELMBDA, F, IDIMF, 
!                                     PERTRB, IERROR)
!
! ARGUMENTS
! ON INPUT               A, B
!
!                          THE RANGE OF X, I.E., A .LE. X .LE. B.
!                          A MUST BE LESS THAN B.
!
!                        M
!                          THE NUMBER OF PANELS INTO WHICH THE
!                          INTERVAL (A, B) IS SUBDIVIDED.
!                          HENCE, THERE WILL BE M+1 GRID POINTS
!                          IN THE X-DIRECTION GIVEN BY
!                          X(I) = A+(I-1)DX FOR I = 1, 2, ..., M+1, 
!                          WHERE DX = (B-A)/M IS THE PANEL WIDTH.
!                          M MUST BE GREATER THAN 3.
!
!                        MBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT X = A AND X = B.
!
!                          = 0  IF THE SOLUTION IS PERIODIC IN X, 
!                               I.E., U(I, J) = U(M+I, J).
!                          = 1  IF THE SOLUTION IS SPECIFIED AT
!                               X = A AND X = B.
!                          = 2  IF THE SOLUTION IS SPECIFIED AT
!                               X = A AND THE DERIVATIVE OF THE
!                               SOLUTION WITH RESPECT TO X IS
!                               SPECIFIED AT X = B.
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO X IS SPECIFIED AT
!                               AT X = A AND X = B.
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO X IS SPECIFIED AT
!                               X = A AND THE SOLUTION IS SPECIFIED
!                               AT X = B.
!
!                        BDA
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
!                          SPECIFIES THE VALUES OF THE DERIVATIVE
!                          OF THE SOLUTION WITH RESPECT TO X AT X = A.
!
!                          WHEN MBDCND = 3 OR 4, 
!
!                            BDA(J) = (D/DX)U(A, Y(J)), J = 1, 2, ..., N+1.
!
!                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
!                          A DUMMY VARIABLE.
!
!                        BDB
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
!                          THAT SPECIFIES THE VALUES OF THE DERIVATIVE
!                          OF THE SOLUTION WITH RESPECT TO X AT X = B.
!
!                          WHEN MBDCND = 2 OR 3, 
!
!                            BDB(J) = (D/DX)U(B, Y(J)), J = 1, 2, ..., N+1
!
!                          WHEN MBDCND HAS ANY OTHER VALUE BDB IS A
!                          DUMMY VARIABLE.
!
!                        C, D
!                          THE RANGE OF Y, I.E., C .LE. Y .LE. D.
!                          C MUST BE LESS THAN D.
!
!                        N
!                          THE NUMBER OF PANELS INTO WHICH THE
!                          INTERVAL (C, D) IS SUBDIVIDED.  HENCE, 
!                          THERE WILL BE N+1 GRID POINTS IN THE
!                          Y-DIRECTION GIVEN BY Y(J) = C+(J-1)DY
!                          FOR J = 1, 2, ..., N+1, WHERE
!                          DY = (D-C)/N IS THE PANEL WIDTH.
!                          N MUST BE GREATER THAN 3.
!
!                        NBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS AT
!                          Y = C AND Y = D.
!
!                          = 0  IF THE SOLUTION IS PERIODIC IN Y, 
!                               I.E., U(I, J) = U(I, N+J).
!                          = 1  IF THE SOLUTION IS SPECIFIED AT
!                               Y = C AND Y = D.
!                          = 2  IF THE SOLUTION IS SPECIFIED AT
!                               Y = C AND THE DERIVATIVE OF THE
!                               SOLUTION WITH RESPECT TO Y IS
!                               SPECIFIED AT Y = D.
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO Y IS SPECIFIED AT
!                               Y = C AND Y = D.
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO Y IS SPECIFIED AT
!                               Y = C AND THE SOLUTION IS SPECIFIED
!                               AT Y = D.
!
!                        BDC
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
!                          SPECIFIES THE VALUES OF THE DERIVATIVE
!                          OF THE SOLUTION WITH RESPECT TO Y AT Y = C.
!
!                          WHEN NBDCND = 3 OR 4, 
!
!                            BDC(I) = (D/DY)U(X(I), C), I = 1, 2, ..., M+1
!
!                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
!                          A DUMMY VARIABLE.
!
!                        BDD
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
!                          SPECIFIES THE VALUES OF THE DERIVATIVE
!                          OF THE SOLUTION WITH RESPECT TO Y AT Y = D.
!
!                          WHEN NBDCND = 2 OR 3, 
!
!                            BDD(I) = (D/DY)U(X(I), D), I = 1, 2, ..., M+1
!
!                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
!                          A DUMMY VARIABLE.
!
!                        ELMBDA
!                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
!                          EQUATION.  IF LAMBDA .GT. 0, A SOLUTION
!                          MAY NOT EXIST.  HOWEVER, HWSCRT WILL
!                          ATTEMPT TO FIND A SOLUTION.
!
!                        F
!                          A TWO-DIMENSIONAL ARRAY, OF DIMENSION AT
!                          LEAST (M+1)*(N+1), SPECIFYING VALUES OF THE
!                          RIGHT SIDE OF THE HELMHOLTZ  EQUATION AND
!                          BOUNDARY VALUES (IF ANY).
!
!                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
!                          FOR I = 2, 3, ..., M AND J = 2, 3, ..., N
!                          F(I, J) = F(X(I), Y(J)).
!
!                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
!                          FOR J=1, 2, ..., N+1,  I=1, 2, ..., M+1, 
!
!                          MBDCND     F(1, J)        F(M+1, J)
!                          ------     ---------     --------
!
!                            0        F(A, Y(J))     F(A, Y(J))
!                            1        U(A, Y(J))     U(B, Y(J))
!                            2        U(A, Y(J))     F(B, Y(J))
!                            3        F(A, Y(J))     F(B, Y(J))
!                            4        F(A, Y(J))     U(B, Y(J))
!
!
!                          NBDCND     F(I, 1)        F(I, N+1)
!                          ------     ---------     --------
!
!                            0        F(X(I), C)     F(X(I), C)
!                            1        U(X(I), C)     U(X(I), D)
!                            2        U(X(I), C)     F(X(I), D)
!                            3        F(X(I), C)     F(X(I), D)
!                            4        F(X(I), C)     U(X(I), D)
!
!                          NOTE:
!                          IF THE TABLE CALLS FOR BOTH THE SOLUTION U
!                          AND THE RIGHT SIDE F AT A CORNER THEN THE
!                          SOLUTION MUST BE SPECIFIED.
!
!                        IDIMF
!                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
!                          F AS IT APPEARS IN THE PROGRAM CALLING
!                          HWSCRT.  THIS PARAMETER IS USED TO SPECIFY
!                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
!                          BE AT LEAST M+1  .
!
!
! ON OUTPUT              F
!                          CONTAINS THE SOLUTION U(I, J) OF THE FINITE
!                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
!                          (X(I), Y(J)), I = 1, 2, ..., M+1, 
!                          J = 1, 2, ..., N+1  .
!
!                        PERTRB
!                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
!                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
!                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
!                          MAY NOT EXIST.  PERTRB IS A CONSTANT, 
!                          CALCULATED AND SUBTRACTED FROM F, WHICH
!                          ENSURES THAT A SOLUTION EXISTS.  HWSCRT
!                          THEN COMPUTES THIS SOLUTION, WHICH IS A
!                          LEAST SQUARES SOLUTION TO THE ORIGINAL
!                          APPROXIMATION.  THIS SOLUTION PLUS ANY
!                          CONSTANT IS ALSO A SOLUTION.  HENCE, THE
!                          SOLUTION IS NOT UNIQUE.  THE VALUE OF
!                          PERTRB SHOULD BE SMALL COMPARED TO THE
!                          RIGHT SIDE F.  OTHERWISE, A SOLUTION IS
!                          OBTAINED TO AN ESSENTIALLY DIFFERENT
!                          PROBLEM. THIS COMPARISON SHOULD ALWAYS
!                          BE MADE TO INSURE THAT A MEANINGFUL
!                          SOLUTION HAS BEEN OBTAINED.
!
!                        IERROR
!                          AN ERROR FLAG THAT INDICATES INVALID INPUT
!                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 6, 
!                          A SOLUTION IS NOT ATTEMPTED.
!
!                          = 0  NO ERROR
!                          = 1  A .GE. B
!                          = 2  MBDCND .LT. 0 OR MBDCND .GT. 4
!                          = 3  C .GE. D
!                          = 4  N .LE. 3
!                          = 5  NBDCND .LT. 0 OR NBDCND .GT. 4
!                          = 6  LAMBDA .GT. 0
!                          = 7  IDIMF .LT. M+1
!                          = 8  M .LE. 3
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if N, M are too large
!                               for your computer)
!
!                          SINCE THIS IS THE ONLY MEANS OF INDICATING
!                          A POSSIBLY INCORRECT CALL TO HWSCRT, THE
!                          USER SHOULD TEST IERROR AFTER THE CALL.
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
! REFERENCES             SWARZTRAUBER, P. AND R. SWEET, "EFFICIENT
!                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
!                        ELLIPTIC EQUATIONS"
!                          NCAR TN/IA-109, JULY, 1975, 138 PP.
!***********************************************************************
!
module module_hwscrt

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_genbun, only: &
        genbunn

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hwscrt
    public :: hwscrt_unit_test

contains
    !
    !*****************************************************************************************
    !
    subroutine hwscrt_unit_test()
        !     file thwscrt.f
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
        !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer :: idimf, m, mbdcnd, n, nbdcnd, mp1, np1, i, j, ierror
        real , dimension(45, 82) :: f
        real , dimension(81) :: bdb, bda, bdc, bdd, y
        real , dimension(41) :: x
        real::a, b, c, d, elmbda, pi, dum, piby2, pisq, pertrb, discretization_error, z
        !-----------------------------------------------
        !
        !     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
        !
        idimf = 45
        a = 0.
        b = 2.
        m = 40
        mbdcnd = 2
        c = -1.
        d = 3.
        n = 80
        nbdcnd = 0
        elmbda = -4.
        !
        !     AUXILIARY QUANTITIES.
        !
        pi = acos( -1.0 )
        piby2 = pi/2.
        pisq = pi**2
        mp1 = m + 1
        np1 = n + 1
        !
        !     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
        !     BOUNDARY DATA AND THE RIGHT SIDE OF THE HELMHOLTZ EQUATION.
        !
        do i = 1, mp1
            x(i) = real(i - 1)/20.
        end do
        do j = 1, np1
            y(j) = (-1.) + real(j - 1)/20.
        end do
        !
        !     GENERATE BOUNDARY DATA.
        !
        do j = 1, np1
            bdb(j) = 4.*COS((Y(j)+1.)*piby2)
        end do
        !
        !     BDA, BDC, AND BDD ARE DUMMY VARIABLES.
        !
        f(1, :np1) = 0.
        !
        !     GENERATE RIGHT SIDE OF EQUATION.
        !
        do i = 2, mp1
            do j = 1, np1
                f(i, j) = (2. - (4. + pisq/4.)*X(i)**2)*COS((Y(j)+1.)*piby2)
            end do
        end do

        call hwscrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
            elmbda, f, idimf, pertrb, ierror)
        !
        !     Compute discretization error.  the exact solution is
        !                u(x, y) = x**2*cos((y+1)*piby2)
        !
        discretization_error = 0.0_wp
        do i = 1, mp1
            do j = 1, np1
                z = abs(f(i, j)-(x(i)**2) * cos((y(j)+1.0_wp) * piby2))
                discretization_error = max(z, discretization_error)
            end do
        end do
        !     Print earlier output from platforms with 32 and 64 bit floating point
        !     arithemtic followed by the output from this computer
        write( stdout, '(A)') ''
        write( stdout, '(A)') '    HWSCRT UNIT TEST '
        write( stdout, '(A)') '    Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') '    IERROR = 0,  Discretization Error = 5.36508-4'
        write( stdout, '(A)') ''
        write( stdout, '(A)') '    The output from your computer is: '
        write( stdout, '(A,I3,A, E23.15E3)')&
            '    IERROR =', ierror, ' Discretization Error = ', discretization_error

    end subroutine hwscrt_unit_test
    !
    !*****************************************************************************************
    !
    subroutine hwscrt( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror )
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
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
        real (wp) :: f(idimf, *)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip)             :: irwk, icwk
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------

        !! Check for invalid parameters.
        ierror = 0
        if (a >= b) ierror = 1
        if (mbdcnd<0 .or. mbdcnd>4) ierror = 2
        if (c >= d) ierror = 3
        if (n <= 3) ierror = 4
        if (nbdcnd<0 .or. nbdcnd>4) ierror = 5
        if (idimf < m + 1) ierror = 7
        if (m <= 3) ierror = 8
        if (ierror /= 0) return

        !! allocate required work space length (generous estimate)
        associate( int_arg => real( n + 1, kind=wp ) / log(2.0_wp) )
            irwk = 4 * ( n + 1) + (13 + int( int_arg, kind=ip )) *(m + 1)
        end associate
        icwk = 0

        !! create workspace arrays
        call workspace%create( irwk, icwk, ierror )

        !! check that allocation was successful
        if (ierror == 20) return

        call hwscrtt( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
            elmbda, f, idimf, pertrb, ierror, workspace%rew )

        !! release dynamically allocated work space
        call workspace%destroy()

    end subroutine hwscrt
    !
    !*****************************************************************************************
    !
    subroutine hwscrtt( a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, w )
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in) :: m
        integer (ip), intent (in) :: mbdcnd
        integer (ip), intent (in) :: n
        integer (ip), intent (in) :: nbdcnd
        integer (ip) :: idimf
        integer (ip), intent (out) :: ierror
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
        real (wp) :: f(idimf, *)
        real (wp) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: nperod, mperod, np, np1, mp, mp1, nstart, nstop, nskip
        integer (ip) :: nunk, mstart, mstop, mskip, j, munk, id2, id3, id4, msp1
        integer (ip) :: mstm1, nsp1, nstm1, ierr1
        real (wp)    :: deltax, twdelx, delxsq, deltay, twdely, delysq, s, st2, a1, a2, s1
        !-----------------------------------------------

        nperod = nbdcnd
        mperod = 0
        if (mbdcnd > 0) mperod = 1
        deltax = (b - a)/m
        twdelx = 2.0_wp/deltax
        delxsq = 1.0_wp/deltax**2
        deltay = (d - c)/n
        twdely = 2.0_wp/deltay
        delysq = 1.0_wp/deltay**2
        np = nbdcnd + 1
        np1 = n + 1
        mp = mbdcnd + 1
        mp1 = m + 1
        nstart = 1
        nstop = n
        nskip = 1
        go to (104, 101, 102, 103, 104) np
101 continue
    nstart = 2
    go to 104
102 continue
    nstart = 2
103 continue
    nstop = np1
    nskip = 2
104 continue
    nunk = nstop - nstart + 1
    !
    !     Enter boundary data for x-boundaries.
    !
    mstart = 1
    mstop = m
    mskip = 1
    go to (117, 105, 106, 109, 110) mp
105 continue
    mstart = 2
    go to 107
106 continue
    mstart = 2
    mstop = mp1
    mskip = 2
107 continue
    f(2, nstart:nstop) = f(2, nstart:nstop) - f(1, nstart:nstop)*delxsq
    go to 112
109 continue
    mstop = mp1
    mskip = 2
110 continue
    f(1, nstart:nstop) = f(1, nstart:nstop) + bda(nstart:nstop)*twdelx
112 continue
    select case (mskip)
        case default
            f(m, nstart:nstop) = f(m, nstart:nstop) - f(mp1, nstart:nstop)* &
                delxsq
        case (2)
            f(mp1, nstart:nstop) = f(mp1, nstart:nstop) - bdb(nstart:nstop)* &
                twdelx
    end select
117 continue
    munk = mstop - mstart + 1
    !
    !     Enter boundary data for y-boundaries.
    !
    go to (127, 118, 118, 120, 120) np
118 continue
    f(mstart:mstop, 2) = f(mstart:mstop, 2) - f(mstart:mstop, 1)*delysq
    go to 122
120 continue
    f(mstart:mstop, 1) = f(mstart:mstop, 1) + bdc(mstart:mstop)*twdely
122 continue
    select case (nskip)
        case default
            f(mstart:mstop, n) = f(mstart:mstop, n) - f(mstart:mstop, np1)* &
                delysq
        case (2)
            f(mstart:mstop, np1) = f(mstart:mstop, np1) - bdd(mstart:mstop)* &
                twdely
    end select
!
!    Multiply right side by deltay**2.
!
127 continue
    delysq = deltay**2
    f(mstart:mstop, nstart:nstop) = f(mstart:mstop, nstart:nstop)*delysq
    !
    !     Define the a, b, c coefficients in w-array.
    !
    id2 = munk
    id3 = id2 + munk
    id4 = id3 + munk
    s = delysq*delxsq
    st2 = 2.0_wp*s
    w(:munk) = s
    w(id2+1:munk+id2) = (-st2) + elmbda*delysq
    w(id3+1:munk+id3) = s

    if (mp /= 1) then
        w(1) = 0.0_wp
        w(id4) = 0.0_wp
    end if

    go to (135, 135, 132, 133, 134) mp
132 continue
    w(id2) = st2
    go to 135
133 continue
    w(id2) = st2
134 continue
    w(id3+1) = st2
135 continue
    pertrb = 0.0_wp
    if ( elmbda >= 0.0_wp ) then
        if ( elmbda /= 0.0_wp ) then
            ierror = 6
        else
            if ((nbdcnd==0 .or. nbdcnd==3).and.(mbdcnd==0 .or. mbdcnd ==3)) then
                !
                !     For singular problems must adjust data to insure that a solution
                !     will exist.
                !
                a1 = 1.0_wp
                a2 = 1.0_wp
                if (nbdcnd == 3) then
                    a2 = 2.0_wp
                end if
                if (mbdcnd == 3) then
                    a1 = 2.0_wp
                end if
                s1 = 0.0_wp
                msp1 = mstart + 1
                mstm1 = mstop - 1
                nsp1 = nstart + 1
                nstm1 = nstop - 1
                do j = nsp1, nstm1
                    s = 0.0_wp
                    s = sum(f(msp1:mstm1, j))
                    s1 = s1 + s*a1 + f(mstart, j) + f(mstop, j)
                end do
                s1 = a2*s1
                s = 0.0_wp
                s = sum(f(msp1:mstm1, nstart) + f(msp1:mstm1, nstop))
                s1 = s1 + s*a1 + f(mstart, nstart) + f(mstart, nstop) &
                    + f( mstop, nstart) + f(mstop, nstop)
                s = (2.0_wp + real(nunk - 2, kind = wp)*a2) * (2.0_wp + real(munk - 2, kind=wp)*a1)
                pertrb = s1/s
                f(mstart:mstop, nstart:nstop) = f(mstart:mstop, nstart: nstop) - pertrb
                pertrb = pertrb/delysq
            !
            !     Solve the equation.
            !
            end if
        end if
    end if

    ierr1 = 0
    call genbunn( nperod, nunk, mperod, munk, w(1), w(id2+1), w(id3+1), &
        idimf, f(mstart, nstart), ierr1, w(id4+1))

    !
    !     Fill in identical values when have periodic boundary conditions.
    !
    if (nbdcnd == 0) then
        f(mstart:mstop, np1) = f(mstart:mstop, 1)
    end if

    if (mbdcnd == 0) then
        f(mp1, nstart:nstop) = f(1, nstart:nstop)
        if (nbdcnd == 0) then
            f(mp1, np1) = f(1, np1)
        end if
    end if

end subroutine hwscrtt
    !
    !*****************************************************************************************
    !
end module module_hwscrt
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
