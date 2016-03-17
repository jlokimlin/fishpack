module module_poistg

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_gnbnaux, only: &
        cosgen, &
        merge_rename, &
        trix, &
        tri3

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: poistg_unit_test
    public :: poistg
    public :: poistgg

contains

    subroutine poistg_unit_test()
        !     file tpoistg.f
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
        !     PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE poistg TO
        !     SOLVE THE EQUATION
        !
        !     (1/COS(X))(D/DX)(COS(X)(DU/DX)) + (D/DY)(DU/DY) =
        !
        !           2*Y**2*(6-Y**2)*SIN(X)
        !
        !     ON THE RECTANGLE -PI/2 .LT. X .LT. PI/2 AND
        !     0 .LT. Y .LT. 1 WITH THE BOUNDARY CONDITIONS
        !
        !     (DU/DX) (-PI/2, Y) = (DU/DX)(PI/2, Y) = 0 , 0 .LE. Y .LE. 1  (2)
        !
        !     U(X, 0) = 0                                           (3)
        !                                 -PI/2 .LE. X .LE. PI/2
        !     (DU/DY)(X, 1) = 4SIN(X)                               (4)
        !
        !     USING FINITE DIFFERENCES ON A STAGGERED GRID WITH
        !     DELTAX (= DX) = PI/40 AND DELTAY (= DY) = 1/20 .
        !        TO SET UP THE FINITE DIFFERENCE EQUATIONS WE DEFINE
        !     THE GRID POINTS
        !
        !     X(I) = -PI/2 + (I-0.5)DX            I=1, 2, ..., 40
        !
        !     Y(J) = (J-O.5)DY                    J=1, 2, ..., 20
        !
        !     AND LET V(I, J) BE AN APPROXIMATION TO U(X(I), Y(J)).
        !     NUMBERING THE GRID POINTS IN THIS FASHION GIVES THE SET
        !     OF UNKNOWNS AS V(I, J) FOR I=1, 2, ..., 40 AND J=1, 2, ..., 20.
        !     HENCE, IN THE PROGRAM M = 40 AND N = 20.  AT THE INTERIOR
        !     GRID POINT (X(I), Y(J)), WE REPLACE ALL DERIVATIVES IN
        !     EQUATION (1) BY SECOND ORDER CENTRAL FINITE DIFFERENCES,
        !     MULTIPLY BY DY**2, AND COLLECT COEFFICIENTS OF V(I, J) TO
        !     GET THE FINITE DIFFERENCE EQUATION
        !
        !     A(I)V(I-1, J) + B(I)V(I, J) + C(I)V(I+1, J)
        !
        !     + V(I, J-1) - 2V(I, J) + V(I, J+1) = F(I, J)            (5)
        !
        !     WHERE S = (DY/DX)**2, AND FOR I=2, 3, ..., 39
        !
        !     A(I) = S*COS(X(I)-DX/2)
        !
        !     B(I) = -S*(COS(X(I)-DX/2)+COS(X(I)+DX/2))
        !
        !     C(I) = S*COS(X(I)+DX/2)
        !
        !     F(I, J) = 2DY**2*Y(J)**2*(6-Y(J)**2)*SIN(X(I)) , J=1, 2, ..., 19.
        !
        !        TO OBTAIN EQUATIONS FOR I = 1, WE REPLACE EQUATION (2)
        !     BY THE SECOND ORDER APPROXIMATION
        !
        !     (V(1, J)-V(0, J))/DX = 0
        !
        !     AND USE THIS EQUATION TO ELIMINATE V(0, J) IN EQUATION (5)
        !     TO ARRIVE AT THE EQUATION
        !
        !     B(1)V(1, J) + C(1)V(2, J) + V(1, J-1) - 2V(1, J) + V(1, J+1)
        !
        !                       = F(1, J)
        !
        !     WHERE
        !
        !     B(1) = -S*(COS(X(1)-DX/2)+COS(X(1)+DX/2))
        !
        !     C(1) = -B(1)
        !
        !     FOR COMPLETENESS, WE SET A(1) = 0.
        !        TO OBTAIN EQUATIONS FOR I = 40, WE REPLACE THE DERIVATIVE
        !     IN EQUATION (2) AT X=PI/2 IN A SIMILAR FASHION, USE THIS
        !     EQUATION TO ELIMINATE THE VIRTUAL UNKNOWN V(41, J) IN EQUATION
        !     (5) AND ARRIVE AT THE EQUATION
        !
        !     A(40)V(39, J) + B(40)V(40, J)
        !
        !     + V(40, J-1) - 2V(40, J) + V(40, J+1) = F(40, J)
        !
        !     WHERE
        !
        !     A(40) = -B(40) = -S*(COS(X(40)-DX/2)+COS(X(40)+DX/2))
        !
        !     FOR COMPLETENESS, WE SET C(40) = 0.  HENCE, IN THE
        !     PROGRAM MPEROD = 1.
        !        FOR J = 1, WE REPLACE EQUATION (3) BY THE SECOND ORDER
        !     APPROXIMATION
        !
        !                (V(I, 0) + V(I, 1))/2 = 0
        !
        !     TO ARRIVE AT THE CONDITION
        !
        !                V(I, 0) = -V(I, 1) .
        !
        !     FOR J = 20, WE REPLACE EQUATION (4) BY THE SECOND ORDER
        !     APPROXIMATION
        !
        !                (V(I, 21) - V(I, 20))/DY = 4*SIN(X)
        !
        !     AND COMBINE THIS EQUATION WITH EQUATION (5) TO ARRIVE AT
        !     THE EQUATION
        !
        !     A(I)V(I-1, 20) + B(I)V(I, 20) + C(I)V(I+1, 20)
        !
        !     + V(I, 19) - 2V(I, 20) + V(I, 21) = F(I, 20)
        !
        !     WHERE
        !
        !     V(I, 21) = V(I, 20)  AND
        !
        !     F(I, 20) = 2*DY**2*Y(J)**2*(6-Y(J)**2)*SIN(X(I)) - 4*DY*SIN(X(I))
        !
        !     HENCE, IN THE PROGRAM NPEROD = 2 .
        !        THE EXACT SOLUTION TO THIS PROBLEM IS
        !
        !        U(X, Y) = Y**4*COS(X) .
        !
        !
        !     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF = 42
        !
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer :: idimf, mperod, m, nperod, n, i, j, ierror
        real (wp), dimension(42, 20) :: f
        real (wp), dimension(40) :: a, b, c, x
        real (wp), dimension(20) :: y
        real :: pi, dx, dy, s, discretization_error, t
        !-----------------------------------------------

        idimf = 42
        mperod = 1
        m = 40
        pi = acos( -1.0 )
        dx = pi/real(m)
        nperod = 2
        n = 20
        dy = 1./real(n)
        !
        !     GENERATE AND STORE GRID POINTS FOR COMPUTATION.
        !
        do i = 1, m
            x(i) = (-pi/2.) + (real(i) - 0.5)*dx
        end do
        do j = 1, n
            y(j) = (real(j) - 0.5)*dy
        end do
        !
        !     GENERATE COEFFICIENTS .
        !
        s = (dy/dx)**2
        a(1) = 0.
        b(1) = -s*COS((-pi/2.) + dx)/COS(X(1))
        c(1) = -B(1)
        do i = 2, m
            a(i) = s*COS(X(i)-dx/2.)/COS(X(i))
            c(i) = s*COS(X(i)+dx/2.)/COS(X(i))
            b(i) = -(A(i)+C(i))
        end do
        a(40) = -B(40)
        c(40) = 0.
        !
        !     GENERATE RIGHT SIDE OF EQUATION.
        !
        do i = 1, m
            do j = 1, n
                f(i, j) = 2.*dy**2*Y(j)**2*(6. - Y(j)**2)*SIN(X(i))
            end do
        end do
        do i = 1, m
            f(i, n) = F(i, n) - 4.*dy*SIN(X(i))
        end do

        call poistg( nperod, n, mperod, m, a, b, c, idimf, f, ierror)
        !
        !     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
        !
        !     U(X, Y) = Y**4*SIN(X)
        !
        discretization_error = 0.0_wp
        do i = 1, m
            do j = 1, n
                associate( local_error => abs(f(i, j)-y(j)**4*sin(x(i))) )
                    discretization_error = max(local_error, discretization_error)
                end associate
            end do
        end do

        write( stdout, '(A)') ''
        write( stdout, '(A)') '    poistg UNIT TEST *** '
        write( stdout, '(A)') '    Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') '    IERROR = 0,  Discretization Error = 5.6417E-4'
        write( stdout, '(A)') ''
        write( stdout, '(A)') '    The output from your computer is: '
        write( stdout, *)&
            '    IERROR =', ierror, ' Discretization Error = ', discretization_error

    end subroutine poistg_unit_test
    !
    !*****************************************************************************************
    !
    subroutine poistg( nperod, n, mperod, m, a, b, c, idimy, y, ierror )
        !
        !     file poistg.f
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
        !     SUBROUTINE poistg (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, IERROR)
        !
        !
        ! DIMENSION OF           A(M),  B(M),  C(M),  Y(IDIMY, N)
        ! ARGUMENTS
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
        !                        FOR UNKNOWN X VALUES, WHERE I=1, 2, ..., M
        !                        AND J=1, 2, ..., N
        !
        !                        A(I)*X(I-1, J) + B(I)*X(I, J) + C(I)*X(I+1, J)
        !                        + X(I, J-1) - 2.*X(I, J) + X(I, J+1)
        !                        = Y(I, J)
        !
        !                        THE INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
        !                        I.E. X(0, J) = X(M, J) AND X(M+1, J) = X(1, J), AND
        !                        X(I, 0) MAY BE EQUAL TO X(I, 1) OR -X(I, 1), AND
        !                        X(I, N+1) MAY BE EQUAL TO X(I, N) OR -X(I, N),
        !                        DEPENDING ON AN INPUT PARAMETER.
        !
        ! USAGE                  CALL poistg (NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y,
        !                                     IERROR)
        !
        ! ARGUMENTS
        !
        ! ON INPUT
        !
        !                        NPEROD
        !                          INDICATES VALUES WHICH X(I, 0) AND X(I, N+1)
        !                          ARE ASSUMED TO HAVE.
        !                          = 1 IF X(I, 0) = -X(I, 1) AND X(I, N+1) = -X(I, N
        !                          = 2 IF X(I, 0) = -X(I, 1) AND X(I, N+1) =  X(I, N
        !                          = 3 IF X(I, 0) =  X(I, 1) AND X(I, N+1) =  X(I, N
        !                          = 4 IF X(I, 0) =  X(I, 1) AND X(I, N+1) = -X(I, N
        !
        !                        N
        !                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
        !                          N MUST BE GREATER THAN 2.
        !
        !                        MPEROD
        !                          = 0 IF A(1) AND C(M) ARE NOT ZERO
        !                          = 1 IF A(1) = C(M) = 0
        !
        !                        M
        !                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
        !                          M MUST BE GREATER THAN 2.
        !
        !                        A, B, C
        !                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
        !                          SPECIFY THE COEFFICIENTS IN THE LINEAR
        !                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0 THE
        !                          ARRAY ELEMENTS MUST NOT DEPEND ON INDEX I,
        !                          BUT MUST BE CONSTANT.  SPECIFICALLY, THE
        !                          SUBROUTINE CHECKS THE FOLLOWING CONDITION
        !                            A(I) = C(1)
        !                            B(I) = B(1)
        !                            C(I) = C(1)
        !                          FOR I = 1, 2, ..., M.
        !
        !                        IDIMY
        !                          THE ROW (OR FIRST) DIMENSION OF THE TWO-
        !                          DIMENSIONAL ARRAY Y AS IT APPEARS IN THE
        !                          PROGRAM CALLING poistg.  THIS PARAMETER IS
        !                          USED TO SPECIFY THE VARIABLE DIMENSION OF Y.
        !                          IDIMY MUST BE AT LEAST M.
        !
        !                        Y
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
        !                          OF EQUATIONS GIVEN ABOVE.
        !                          Y MUST BE DIMENSIONED AT LEAST M X N.
        !
        ! ON OUTPUT
        !
        !                        Y
        !                          CONTAINS THE SOLUTION X.
        !
        !                        IERROR
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
        !                          SOLUTION IS NOT ATTEMPTED.
        !                          = 0  NO ERROR
        !                          = 1  IF M .LE. 2
        !                          = 2  IF N .LE. 2
        !                          = 3  IDIMY .LT. M
        !                          = 4  IF NPEROD .LT. 1 OR NPEROD .GT. 4
        !                          = 5  IF MPEROD .LT. 0 OR MPEROD .GT. 1
        !                          = 6  IF MPEROD = 0 AND A(I) .NE. C(1)
        !                               OR B(I) .NE. B(1) OR C(I) .NE. C(1)
        !                               FOR SOME I = 1, 2, ..., M.
        !                          = 7  IF MPEROD .EQ. 1 .AND.
        !                               (A(1).NE.0 .OR. C(M).NE.0)
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
        !
        !                          SINCE THIS IS THE ONLY MEANS OF INDICATING A
        !                          POSSIBLY INCORRECT CALL TO POIS3D, THE USER
        !                          SHOULD TEST IERROR AFTER THE CALL.
        !
        !
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED files         fish.f, gnbnaux.f, comf.f
        !
        ! LANGUAGE               FORTRAN 90
        !
        ! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
        !                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
        !                        LIBRARIES IN JANUARY, 1980.
        !                        Revised in June 2004 by John Adams using
        !                        Fortran 90 dynamically allocated work space.
        !
        ! PORTABILITY            FORTRAN 90
        !
        ! ALGORITHM              THIS SUBROUTINE IS AN IMPLEMENTATION OF THE
        !                        ALGORITHM PRESENTED IN THE REFERENCE BELOW.
        !
        ! TIMING                 FOR LARGE M AND N, THE EXECUTION TIME IS
        !                        ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
        !
        ! ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
        !                        UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
        !                        CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
        !                        IN THE 'PURPOSE' SECTION ABOVE, WITH
        !                          A(I) = C(I) = -0.5*B(I) = 1,    I=1, 2, ..., M
        !                        AND, WHEN MPEROD = 1
        !                          A(1) = C(M) = 0
        !                          B(1) = B(M) =-1.
        !
        !                        THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
        !                        SYSTEM AND, USING DOUBLE PRECISION, A RIGHT SID
        !                        Y WAS COMPUTED.  USING THIS ARRAY Y SUBROUTINE
        !                        poistg WAS CALLED TO PRODUCE AN APPROXIMATE
        !                        SOLUTION Z.  THEN THE RELATIVE ERROR, DEFINED A
        !                          E = MAX(abs(Z(I, J)-X(I, J)))/MAX(abs(X(I, J)))
        !                        WHERE THE TWO MAXIMA ARE TAKEN OVER I=1, 2, ..., M
        !                        AND J=1, 2, ..., N, WAS COMPUTED.  VALUES OF E ARE
        !                        GIVEN IN THE TABLE BELOW FOR SOME TYPICAL VALUE
        !                        OF M AND N.
        !
        !                        M (=N)    MPEROD    NPEROD      E
        !                        ------    ------    ------    ------
        !
        !                          31        0-1       1-4     9.E-13
        !                          31        1         1       4.E-13
        !                          31        1         3       3.E-13
        !                          32        0-1       1-4     3.E-12
        !                          32        1         1       3.E-13
        !                          32        1         3       1.E-13
        !                          33        0-1       1-4     1.E-12
        !                          33        1         1       4.E-13
        !                          33        1         3       1.E-13
        !                          63        0-1       1-4     3.E-12
        !                          63        1         1       1.E-12
        !                          63        1         3       2.E-13
        !                          64        0-1       1-4     4.E-12
        !                          64        1         1       1.E-12
        !                          64        1         3       6.E-13
        !                          65        0-1       1-4     2.E-13
        !                          65        1         1       1.E-11
        !                          65        1         3       4.E-13
        !
        ! REFERENCES             SCHUMANN, U. AND R. SWEET, "A DIRECT METHOD
        !                        FOR THE SOLUTION OF POISSON"S EQUATION WITH
        !                        NEUMANN BOUNDARY CONDITIONS ON A STAGGERED
        !                        GRID OF ARBITRARY SIZE, " J. COMP. PHYS.
        !                        20(1976), PP. 171-182.
        ! *********************************************************************
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip),          intent (in)     :: nperod
        integer (ip),          intent (in)     :: n
        integer (ip),          intent (in)     :: mperod
        integer (ip),          intent (in)     :: m
        integer (ip),          intent (in)     :: idimy
        integer (ip),          intent (out)    :: ierror
        real (wp), contiguous, intent (in)     :: a(:)
        real (wp), contiguous, intent (in)     :: b(:)
        real (wp), contiguous, intent (in)     :: c(:)
        real (wp), contiguous, intent (in out) :: y(:,:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)             :: i !! Counter
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------

        ! Initialize error flag
        ierror = 0

        ! Perform sanity check: case 1
        if (m <= 2) then
            ierror = 1
            return
        end if

        ! Perform sanity check: case 2
        if (n <= 2) then
            ierror = 2
            return
        end if

        ! Perform sanity check: case 3
        if (idimy < m) then
            ierror = 3
            return
        end if

        ! Perform sanity check: case 4
        if (nperod < 1 .or. nperod > 4) then
            ierror = 4
            return
        end if

        ! Perform sanity check: case 5
        if (mperod < 0 .or. mperod > 1) then
            ierror = 5
            return
        end if

        ! Perform sanity check: case 6
        if (mperod /= 1) then
            do i = 1, m
                if (A(i) /= C(1)) then
                    ierror = 6
                    exit
                end if
                if (C(i) /= C(1)) then
                    ierror = 6
                    exit
                end if
                if (B(i) /= B(1)) then
                    ierror = 6
                    exit
                end if
            end do
        end if

        ! Perform sanity check: case 7
        if (A(1)/= 0.0_wp .or. C(m)/= 0.0_wp) then
            ierror = 7
            return
        end if

        ! Final sanity check
        if (ierror /= 0) then
            return
        end if

        ! Compute and allocate real work space for poistg

        associate( &
            irwk => m * (9 + int(log(real(n, kind=wp))/log(2.0_wp),kind=ip)) + 4 * n, &
            icwk => 0 &
            )
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! return if allocation failed (e.g., if n, m are too large)
        if (ierror == 20) return

        ! Solve system
        call  poistgg( nperod, n, mperod, m, a, b, c, idimy, y, ierror, workspace%rew )

        ! release work space
        call workspace%destroy()

    end subroutine poistg


    subroutine poistgg(nperod, n, mperod, m, a, b, c, idimy, y, ierror, w)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip), intent (in)     :: nperod
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: mperod
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: idimy
        integer (ip), intent (in)     :: ierror
        real (wp),    intent (in)     :: a(*)
        real (wp),    intent (in)     :: b(*)
        real (wp),    intent (in)     :: c(*)
        real (wp),    intent (in out) :: y(idimy, *)
        real (wp),    intent (in out) :: w(*)
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer (ip) :: iwba, iwbb, iwbc, iwb2, iwb3, iww1, iww2, iww3, iwd
        integer (ip) :: iwtcos, iwp, i, k, j, np, mp, ipstor, irev, mh, mhm1, modd
        integer (ip) :: mhpi, mhmi, nby2, mskip
        real (wp)    :: a1
        !-----------------------------------------------
        iwba = m + 1
        iwbb = iwba + m
        iwbc = iwbb + m
        iwb2 = iwbc + m
        iwb3 = iwb2 + m
        iww1 = iwb3 + m
        iww2 = iww1 + m
        iww3 = iww2 + m
        iwd = iww3 + m
        iwtcos = iwd + m
        iwp = iwtcos + 4*n
        do i = 1, m
            k = iwba + i - 1
            w(k) = -A(i)
            k = iwbc + i - 1
            w(k) = -C(i)
            k = iwbb + i - 1
            w(k) = 2. - B(i)
            y(i, :n) = -Y(i, :n)
        end do
        np = nperod
        mp = mperod + 1
        go to (110, 107) mp
107 continue
    go to (108, 108, 108, 119) nperod
108 continue
    call postg2 (np, n, m, W(iwba), W(iwbb), W(iwbc), idimy, y, w, W( &
        iwb2), W(iwb3), W(iww1), W(iww2), W(iww3), W(iwd), W(iwtcos), W &
        (iwp))
    ipstor = W(iww1)
    irev = 2
    if (nperod == 4) go to 120
109 continue
    go to (123, 129) mp
110 continue
    mh = (m + 1)/2
    mhm1 = mh - 1
    modd = 1
    if (mh*2 == m) modd = 2
    do j = 1, n
        do i = 1, mhm1
            w(i) = Y(mh-i, j) - Y(i+mh, j)
            w(i+mh) = Y(mh-i, j) + Y(i+mh, j)
        end do
        w(mh) = 2.*Y(mh, j)
        go to (113, 112) modd
112 continue
    w(m) = 2.*Y(m, j)
113 continue
    y(:m, j) = W(:m)
end do
k = iwbc + mhm1 - 1
i = iwba + mhm1
w(k) = 0.
w(i) = 0.
w(k+1) = 2.*W(k+1)
select case (modd) 
    case default
        k = iwbb + mhm1 - 1
        w(k) = W(k) - W(i-1)
        w(iwbc-1) = W(iwbc-1) + W(iwbb-1)
    case (2)
        w(iwbb-1) = W(k+1)
end select
118 continue
    go to 107
119 continue
    irev = 1
    nby2 = n/2
    np = 2
120 continue
    do j = 1, nby2
        mskip = n + 1 - j
        do i = 1, m
            a1 = Y(i, j)
            y(i, j) = Y(i, mskip)
            y(i, mskip) = a1
        end do
    end do
    go to (108, 109) irev
123 continue
    do j = 1, n
        w(mh-1:mh-mhm1:(-1)) = 0.5*(Y(mh+1:mhm1+mh, j)+Y(:mhm1, j))
        w(mh+1:mhm1+mh) = 0.5*(Y(mh+1:mhm1+mh, j)-Y(:mhm1, j))
        w(mh) = 0.5*Y(mh, j)
        go to (126, 125) modd
125 continue
    w(m) = 0.5*Y(m, j)
126 continue
    y(:m, j) = W(:m)
end do
129 continue
    w(1) = ipstor + iwp - 1

end subroutine poistgg


subroutine postg2(nperod, n, m, a, bb, c, idimq, q, b, b2, b3, w, &
    w2, w3, d, tcos, p)
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip), intent (in) :: nperod
    integer (ip), intent (in) :: n
    integer (ip), intent (in) :: m
    integer (ip), intent (in) :: idimq
    real (wp) :: a(*)
    real (wp) :: bb(*)
    real (wp) :: c(*)
    real (wp), intent (in out) :: q(idimq, *)
    real (wp) :: b(*)
    real (wp) :: b2(*)
    real (wp) :: b3(*)
    real (wp) :: w(*)
    real (wp) :: w2(*)
    real (wp) :: w3(*)
    real (wp) :: d(*)
    real (wp) :: tcos(*)
    real (wp), intent (in out) :: p(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer (ip) :: k(4)
    integer (ip) :: k1, k2, k3, k4, np, mr, ipp, ipstor, i2r, jr, nr, nlast
    integer (ip) :: kr, lr, nrod, jstart, jstop, i2rby2, j, ijump, jp1, jp2, jp3
    integer (ip) :: jm1, jm2, jm3, i, nrodpr, ii, nlastp, jstep
    real (wp)    :: fnum, fnum2, fi, t
    !-----------------------------------------------
    !
    !     SUBROUTINE TO SOLVE POISSON'S EQUATION ON A STAGGERED GRID.
    !
    !
    equivalence (K(1), K1), (K(2), K2), (K(3), K3), (K(4), K4)
    np = nperod
    fnum = 0.5*real(np/3)
    fnum2 = 0.5*real(np/2)
    mr = m
    ipp = -mr
    ipstor = 0
    i2r = 1
    jr = 2
    nr = n
    nlast = n
    kr = 1
    lr = 0
    if (nr > 3) then
101 continue
    jr = 2*i2r
    nrod = 1
    if ((nr/2)*2 == nr) nrod = 0
    jstart = 1
    jstop = nlast - jr
    if (nrod == 0) jstop = jstop - i2r
    i2rby2 = i2r/2
    if (jstop < jstart) then
        j = jr
    else
        ijump = 1
        do j = jstart, jstop, jr
            jp1 = j + i2rby2
            jp2 = j + i2r
            jp3 = jp2 + i2rby2
            jm1 = j - i2rby2
            jm2 = j - i2r
            jm3 = jm2 - i2rby2
            if (j == 1) then
                call COSGEN (i2r, 1, fnum, 0.5, tcos)
                if (i2r == 1) then
                    b(:mr) = Q(:mr, 1)
                    q(:mr, 1) = Q(:mr, 2)
                    go to 112
                end if
                b(:mr) = Q(:mr, 1) + 0.5*(Q(:mr, jp2)-Q(:mr, jp1)-Q(:mr, &
                    jp3))
                q(:mr, 1) = Q(:mr, jp2) + Q(:mr, 1) - Q(:mr, jp1)
                go to 112
            end if
            go to (107, 108) ijump
107     continue
        ijump = 2
        call COSGEN (i2r, 1, 0.5, 0.0, tcos)
108 continue
    if (i2r == 1) then
        b(:mr) = 2.*Q(:mr, j)
        q(:mr, j) = Q(:mr, jm2) + Q(:mr, jp2)
    else
        do i = 1, mr
            fi = Q(i, j)
            q(i, j)=Q(i, j)-Q(i, jm1)-Q(i, jp1)+Q(i, jm2)+Q(i, jp2)
            b(i) = fi + Q(i, j) - Q(i, jm3) - Q(i, jp3)
        end do
    end if
112 continue
    call TRIX (i2r, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = Q(:mr, j) + B(:mr)
!
!     END OF REDUCTION FOR REGULAR UNKNOWNS.
!
end do
!
!     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
!
j = jstop + jr
end if
nlast = j
jm1 = j - i2rby2
jm2 = j - i2r
jm3 = jm2 - i2rby2
if (nrod /= 0) then
    !
    !     ODD NUMBER OF UNKNOWNS
    !
    if (i2r == 1) then
        b(:mr) = Q(:mr, j)
        q(:mr, j) = Q(:mr, jm2)
    else
        b(:mr)=Q(:mr, j)+0.5*(Q(:mr, jm2)-Q(:mr, jm1)-Q(:mr, jm3))
        if (nrodpr == 0) then
            q(:mr, j) = Q(:mr, jm2) + P(ipp+1:mr+ipp)
            ipp = ipp - mr
        else
            q(:mr, j) = Q(:mr, j) - Q(:mr, jm1) + Q(:mr, jm2)
        end if
        if (lr /= 0) call COSGEN (lr, 1, fnum2, 0.5, TCOS(kr+1))
    end if
    call COSGEN (kr, 1, fnum2, 0.5, tcos)
    call TRIX (kr, lr, mr, a, bb, c, b, tcos, d, w)
    q(:mr, j) = Q(:mr, j) + B(:mr)
    kr = kr + i2r
else
    jp1 = j + i2rby2
    jp2 = j + i2r
    if (i2r == 1) then
        b(:mr) = Q(:mr, j)
        tcos(1) = 0.
        call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
        ipp = 0
        ipstor = mr
        p(:mr) = B(:mr)
        b(:mr) = B(:mr) + Q(:mr, n)
        tcos(1) = (-1.) + 2.*real(np/2)
        tcos(2) = 0.
        call TRIX (1, 1, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = Q(:mr, jm2) + P(:mr) + B(:mr)
    else
        b(:mr)=Q(:mr, j)+0.5*(Q(:mr, jm2)-Q(:mr, jm1)-Q(:mr, jm3))
        if (nrodpr == 0) then
            b(:mr) = B(:mr) + P(ipp+1:mr+ipp)
        else
            b(:mr) = B(:mr) + Q(:mr, jp2) - Q(:mr, jp1)
        end if
        call COSGEN (i2r, 1, 0.5, 0.0, tcos)
        call TRIX (i2r, 0, mr, a, bb, c, b, tcos, d, w)
        ipp = ipp + mr
        ipstor = max(ipstor, ipp + mr)
        p(ipp+1:mr+ipp) = B(:mr) + 0.5*(Q(:mr, j)-Q(:mr, jm1)-Q(:mr, &
            jp1))
        b(:mr) = P(ipp+1:mr+ipp) + Q(:mr, jp2)
        if (lr /= 0) then
            call COSGEN (lr, 1, fnum2, 0.5, TCOS(i2r+1))
            call merge_rename(tcos, 0, i2r, i2r, lr, kr)
        else
            do i = 1, i2r
                ii = kr + i
                tcos(ii) = TCOS(i)
            end do
        end if
        call COSGEN (kr, 1, fnum2, 0.5, tcos)
        call TRIX (kr, kr, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = Q(:mr, jm2) + P(ipp+1:mr+ipp) + B(:mr)
    end if
    lr = kr
    kr = kr + jr
end if
nr = (nlast - 1)/jr + 1
if (nr <= 3) go to 142
i2r = jr
nrodpr = nrod
go to 101
end if
142 continue
    j = 1 + jr
    jm1 = j - i2r
    jp1 = j + i2r
    jm2 = nlast - i2r
    if (nr /= 2) then
        if (lr == 0) then
            if (n == 3) then
                !
                !     CASE N = 3.
                !
                go to (143, 148, 143) np
143         continue
            b(:mr) = Q(:mr, 2)
            b2(:mr) = Q(:mr, 1) + Q(:mr, 3)
            b3(:mr) = 0.
            select case (np)
                case default
                    tcos(1) = -1.
                    tcos(2) = 1.
                    k1 = 1
                case (1:2)
                    tcos(1) = -2.
                    tcos(2) = 1.
                    tcos(3) = -1.
                    k1 = 2
            end select
147     continue
        k2 = 1
        k3 = 0
        k4 = 0
        go to 150
148 continue
    b(:mr) = Q(:mr, 2)
    b2(:mr) = Q(:mr, 3)
    b3(:mr) = Q(:mr, 1)
    call COSGEN (3, 1, 0.5, 0.0, tcos)
    tcos(4) = -1.
    tcos(5) = 1.
    tcos(6) = -1.
    tcos(7) = 1.
    k1 = 3
    k2 = 2
    k3 = 1
    k4 = 1
150 continue
    call TRI3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
    b(:mr) = B(:mr) + B2(:mr) + B3(:mr)
    go to (153, 153, 152) np
152 continue
    tcos(1) = 2.
    call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
153 continue
    q(:mr, 2) = B(:mr)
    b(:mr) = Q(:mr, 1) + B(:mr)
    tcos(1) = (-1.) + 4.*fnum
    call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = B(:mr)
    jr = 1
    i2r = 0
    go to 188
end if
!
!     CASE N = 2**P+1
!
b(:mr)=Q(:mr, j)+Q(:mr, 1)-Q(:mr, jm1)+Q(:mr, nlast)-Q(:mr, jm2)
go to (158, 160, 158) np
158 continue
    b2(:mr) = Q(:mr, 1) + Q(:mr, nlast) + Q(:mr, j) - Q(:mr, jm1) - &
        Q(:mr, jp1)
    b3(:mr) = 0.
    k1 = nlast - 1
    k2 = nlast + jr - 1
    call COSGEN (jr - 1, 1, 0.0, 1.0, TCOS(nlast))
    tcos(k2) = 2.*real(np - 2)
    call COSGEN (jr, 1, 0.5 - fnum, 0.5, TCOS(k2+1))
    k3 = (3 - np)/2
    call merge_rename(tcos, k1, jr - k3, k2 - k3, jr + k3, 0)
    k1 = k1 - 1 + k3
    call COSGEN (jr, 1, fnum, 0.5, TCOS(k1+1))
    k2 = jr
    k3 = 0
    k4 = 0
    go to 162
160 continue
    do i = 1, mr
        fi = (Q(i, j)-Q(i, jm1)-Q(i, jp1))/2.
        b2(i) = Q(i, 1) + fi
        b3(i) = Q(i, nlast) + fi
    end do
    k1 = nlast + jr - 1
    k2 = k1 + jr - 1
    call COSGEN (jr - 1, 1, 0.0, 1.0, TCOS(k1+1))
    call COSGEN (nlast, 1, 0.5, 0.0, TCOS(k2+1))
    call merge_rename(tcos, k1, jr - 1, k2, nlast, 0)
    k3 = k1 + nlast - 1
    k4 = k3 + jr
    call COSGEN (jr, 1, 0.5, 0.5, TCOS(k3+1))
    call COSGEN (jr, 1, 0.0, 0.5, TCOS(k4+1))
    call merge_rename(tcos, k3, jr, k4, jr, k1)
    k2 = nlast - 1
    k3 = jr
    k4 = jr
162 continue
    call TRI3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
    b(:mr) = B(:mr) + B2(:mr) + B3(:mr)
    if (np == 3) then
        tcos(1) = 2.
        call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
    end if
    q(:mr, j) = B(:mr) + 0.5*(Q(:mr, j)-Q(:mr, jm1)-Q(:mr, jp1))
    b(:mr) = Q(:mr, j) + Q(:mr, 1)
    call COSGEN (jr, 1, fnum, 0.5, tcos)
    call TRIX (jr, 0, mr, a, bb, c, b, tcos, d, w)
    q(:mr, 1) = Q(:mr, 1) - Q(:mr, jm1) + B(:mr)
    go to 188
end if
!
!     CASE OF GENERAL N WITH NR = 3 .
!
b(:mr) = Q(:mr, 1) - Q(:mr, jm1) + Q(:mr, j)
if (nrod == 0) then
    b(:mr) = B(:mr) + P(ipp+1:mr+ipp)
else
    b(:mr) = B(:mr) + Q(:mr, nlast) - Q(:mr, jm2)
end if
do i = 1, mr
    t = 0.5*(Q(i, j)-Q(i, jm1)-Q(i, jp1))
    q(i, j) = t
    b2(i) = Q(i, nlast) + t
    b3(i) = Q(i, 1) + t
end do
k1 = kr + 2*jr
call COSGEN (jr - 1, 1, 0.0, 1.0, TCOS(k1+1))
k2 = k1 + jr
tcos(k2) = 2.*real(np - 2)
k4 = (np - 1)*(3 - np)
k3 = k2 + 1 - k4
call COSGEN(kr+jr+k4, 1, real(k4)/2., 1.-real(k4), TCOS(k3))
k4 = 1 - np/3
call merge_rename(tcos, k1, jr - k4, k2 - k4, kr + jr + k4, 0)
if (np == 3) k1 = k1 - 1
k2 = kr + jr
k4 = k1 + k2
call COSGEN (kr, 1, fnum2, 0.5, TCOS(k4+1))
k3 = k4 + kr
call COSGEN (jr, 1, fnum, 0.5, TCOS(k3+1))
call merge_rename(tcos, k4, kr, k3, jr, k1)
k4 = k3 + jr
call COSGEN (lr, 1, fnum2, 0.5, TCOS(k4+1))
call merge_rename(tcos, k3, jr, k4, lr, k1 + k2)
call COSGEN (kr, 1, fnum2, 0.5, TCOS(k3+1))
k3 = kr
k4 = kr
call TRI3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
b(:mr) = B(:mr) + B2(:mr) + B3(:mr)
if (np == 3) then
    tcos(1) = 2.
    call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
end if
q(:mr, j) = Q(:mr, j) + B(:mr)
b(:mr) = Q(:mr, 1) + Q(:mr, j)
call COSGEN (jr, 1, fnum, 0.5, tcos)
call TRIX (jr, 0, mr, a, bb, c, b, tcos, d, w)
if (jr == 1) then
    q(:mr, 1) = B(:mr)
    go to 188
end if
q(:mr, 1) = Q(:mr, 1) - Q(:mr, jm1) + B(:mr)
go to 188
end if
b3(:mr) = 0.
b(:mr) = Q(:mr, 1) + P(ipp+1:mr+ipp)
q(:mr, 1) = Q(:mr, 1) - Q(:mr, jm1)
b2(:mr) = Q(:mr, 1) + Q(:mr, nlast)
k1 = kr + jr
k2 = k1 + jr
call COSGEN (jr - 1, 1, 0.0, 1.0, TCOS(k1+1))
go to (182, 183, 182) np
182 continue
    tcos(k2) = 2.*real(np - 2)
    call COSGEN (kr, 1, 0.0, 1.0, TCOS(k2+1))
    go to 184
183 continue
    call COSGEN (kr + 1, 1, 0.5, 0.0, TCOS(k2))
184 continue
    k4 = 1 - np/3
    call merge_rename(tcos, k1, jr - k4, k2 - k4, kr + k4, 0)
    if (np == 3) k1 = k1 - 1
    k2 = kr
    call COSGEN (kr, 1, fnum2, 0.5, TCOS(k1+1))
    k4 = k1 + kr
    call COSGEN (lr, 1, fnum2, 0.5, TCOS(k4+1))
    k3 = lr
    k4 = 0
    call TRI3 (mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
    b(:mr) = B(:mr) + B2(:mr)
    if (np == 3) then
        tcos(1) = 2.
        call TRIX (1, 0, mr, a, bb, c, b, tcos, d, w)
    end if
    q(:mr, 1) = Q(:mr, 1) + B(:mr)
188 continue
    j = nlast - jr
    b(:mr) = Q(:mr, nlast) + Q(:mr, j)
    jm2 = nlast - i2r
    if (jr == 1) then
        q(:mr, nlast) = 0.
    else
        if (nrod == 0) then
            q(:mr, nlast) = P(ipp+1:mr+ipp)
            ipp = ipp - mr
        else
            q(:mr, nlast) = Q(:mr, nlast) - Q(:mr, jm2)
        end if
    end if
    call COSGEN (kr, 1, fnum2, 0.5, tcos)
    call COSGEN (lr, 1, fnum2, 0.5, TCOS(kr+1))
    call TRIX (kr, lr, mr, a, bb, c, b, tcos, d, w)
    q(:mr, nlast) = Q(:mr, nlast) + B(:mr)
    nlastp = nlast
197 continue
    jstep = jr
    jr = i2r
    i2r = i2r/2
    if (jr == 0) go to 210
    jstart = 1 + jr
    kr = kr - jr
    if (nlast + jr <= n) then
        kr = kr - jr
        nlast = nlast + jr
        jstop = nlast - jstep
    else
        jstop = nlast - jr
    end if
    lr = kr - jr
    call COSGEN (jr, 1, 0.5, 0.0, tcos)
    do j = jstart, jstop, jstep
        jm2 = j - jr
        jp2 = j + jr
        if (j == jr) then
            b(:mr) = Q(:mr, j) + Q(:mr, jp2)
        else
            b(:mr) = Q(:mr, j) + Q(:mr, jm2) + Q(:mr, jp2)
        end if
        if (jr == 1) then
            q(:mr, j) = 0.
        else
            jm1 = j - i2r
            jp1 = j + i2r
            q(:mr, j) = 0.5*(Q(:mr, j)-Q(:mr, jm1)-Q(:mr, jp1))
        end if
        call TRIX (jr, 0, mr, a, bb, c, b, tcos, d, w)
        q(:mr, j) = Q(:mr, j) + B(:mr)
    end do
    nrod = 1
    if (nlast + i2r <= n) nrod = 0
    if (nlastp /= nlast) go to 188
    go to 197
210 continue
    w(1) = ipstor

end subroutine postg2


end module module_poistg
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
