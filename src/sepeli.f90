module module_sepeli


    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_blktri, only: &
        blktrii

    use module_sepaux, only: &
        seport, &
        sepmin, &
        septri, &
        sepdx, &
        sepdy

    ! Explicit typing only!
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: sepeli
    public :: sepeli_unit_test

contains

    subroutine sepeli_unit_test()
        !
        !     file tsepeli.f
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
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer::m, n, nx, ny, i, j, mbdcnd, nbdcnd, idmn, intl, iorder, ierror
        real , dimension(33, 33) :: usol, grhs
        real , dimension(33) :: bda, bdb
        real :: a, b, c, d, dlx, dly, x, af, bf, cf, y, df, ef, ff, alpha &
            , beta, dum(1), pertrb, err, err2, err4
        !-----------------------------------------------

        !     DEFINE ARITHMETIC FUNCTIONS GIVING EXACT SOLUTION
        !
        !
        !     SET LIMITS ON REGION
        !
        a = 0.0
        b = 1.0
        c = 0.0
        d = 1.0
        !
        !     SET GRID SIZE
        !
        m = 32
        n = 32
        dlx = (b - a)/real(m)
        dly = (d - c)/real(n)
        nx = m + 1
        ny = n + 1
        do i = 1, nx
            x = a + real(i - 1)*dlx
            !
            !     SET SPECIFIED BOUNDARY CONDITIONS AT Y=C, D
            !
            usol(i, 1) = UE(x, c)
            usol(i, ny) = UE(x, d)
            call get_coefficients_in_x_direction (x, af, bf, cf)
            do j = 1, ny
                y = c + real(j - 1)*dly
                call get_coefficients_in_y_direction (y, df, ef, ff)
                !
                !     SET RIGHT HAND SIDE
                !
                grhs(i, j) = af*UXXE(x, y) + bf*UXE(x, y) + cf*UE(x, y) + df* &
                    UYYE(x, y) + ef*UYE(x, y) + ff*UE(x, y)
            end do
        end do
        !
        !     SET MIXED BOUNDARY CONDITIONS AT X=A, B
        !
        alpha = 1.0
        beta = 1.0
        do j = 1, ny
            y = c + real(j - 1)*dly
            bda(j) = UXE(a, y) + alpha*UE(a, y)
            bdb(j) = UXE(b, y) + beta*UE(b, y)
        end do
        !
        !     SET BOUNDARY SWITHCES
        !
        mbdcnd = 3
        nbdcnd = 1
        !
        !     SET FIRST DIMENSION OF USOL, GRHS
        !
        idmn = 33
        !     set for initialization of sepeli
        intl = 0
        !
        !     OBTAIN SECOND ORDER APPROXIMATION
        !
        iorder = 2
        call SEPELI(intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta &
            , c, d, n, nbdcnd, dum(1:1), dum(1), dum(1:1), dum(1), &
            get_coefficients_in_x_direction, get_coefficients_in_y_direction, grhs, usol, &
            idmn, workspace, pertrb, ierror)
        err = 0.0
        do i = 1, nx
            x = a + real(i - 1)*dlx
            do j = 1, ny
                y = c + real(j - 1)*dly
                err = max(err, abs((USOL(i, j)-UE(x, y))/UE(x, y)))
            end do
        end do
        err2 = err
        !
        !     OBTAIN FOURTH ORDER APPROXIMATION
        !
        iorder = 4
        !
        !     NON-INITIAL CALL
        !
        intl = 1
        call SEPELI (intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta &
            , c, d, n, nbdcnd, dum(1:1), dum(1), dum(1:1), dum(1), &
            get_coefficients_in_x_direction, get_coefficients_in_y_direction, grhs, usol, &
            idmn, workspace, pertrb, ierror)
        !
        !     COMPUTE DISCRETIZATION ERROR
        !
        err = 0.0
        do j = 1, ny
            y = c + real(j - 1)*dly
            do i = 1, nx
                x = a + real(i - 1)*dlx
                err = max(err, abs((USOL(i, j)-UE(x, y))/UE(x, y)))
            end do
        end do
        err4 = err
        !     Print earlier output from platforms with 32 and 64 bit floating point
        !     arithemtic followed by the output from this computer
        write( *, *) ''
        write( *, *) '    SEPELI TEST RUN *** '
        write( *, *) &
            '    Previous 64 bit floating point arithmetic result '
        write( *, *) '    IERROR = 0'
        write( *, *) '    Second Order Discretization Error = 9.7891E-5'
        write( *, *) '    Fourth Order Discretization Error = 1.4735E-6'

        write( *, *) '    The output from your computer is: '
        write( *, *) '    IERROR =', ierror
        write( *, *) '    Second Order Discretization Error =', err2
        write( *, *) '    Fourth Order Discretization Error =', err4
        !     release dynamically allocated real and complex work space
        call workspace%destroy()

    contains

        real function UE (s, t)
            real, intent (in) :: s
            real, intent (in) :: t
            UE = (s*t)**3 + 1.0
            return
        end function UE

        real function UXE (s, t)
            real, intent (in) :: s
            real, intent (in) :: t
            UXE = 3.0*s**2*t**3
            return
        end function UXE

        real function UXXE (s, t)
            real, intent (in) :: s
            real, intent (in) :: t
            UXXE = 6.0*s*t**3
            return
        end function UXXE


        real function UYE (s, t)
            real, intent (in) :: s
            real, intent (in) :: t
            UYE = 3.0*s**3*t**2
            return
        end function UYE


        real function UYYE (s, t)
            real, intent (in) :: s
            real, intent (in) :: t
            UYYE = 6.0*s**3*t
            return
        end function UYYE

        subroutine get_coefficients_in_x_direction(x, af, bf, cf)
            !-----------------------------------------------
            !   D u m m y   A r g u m e n t s
            !-----------------------------------------------
            real , intent (in) :: x
            real , intent (out) :: af
            real , intent (out) :: bf
            real , intent (out) :: cf
            !-----------------------------------------------
            !
            !     SET COEFFICIENTS IN THE X-DIRECTION.
            !
            af = (x + 1.)**2
            bf = 2.0*(x + 1.)
            cf = -x

        end subroutine get_coefficients_in_x_direction


        subroutine get_coefficients_in_y_direction(y, df, ef, ff)
            !-----------------------------------------------
            !   D u m m y   A r g u m e n t s
            !-----------------------------------------------
            real , intent (in) :: y
            real , intent (out) :: df
            real , intent (out) :: ef
            real , intent (out) :: ff
            !-----------------------------------------------
            !
            !     SET COEFFICIENTS IN Y DIRECTION
            !
            df = exp(y)
            ef = 0.0
            ff = -y

        end subroutine get_coefficients_in_y_direction

    end subroutine sepeli_unit_test

    subroutine SEPELI(intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, &
        beta, c, d, n, nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, grhs, &
        usol, idmn, w, pertrb, ierror)
        !
        !     file sepeli.f
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
        !     SUBROUTINE SEPELI (INTL, IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA, C,
        !    +                   D, N, NBDCND, BDC, GAMA, BDD, XNU, COFX, COFY, GRHS,
        !    +                   USOL, IDMN, W, PERTRB, IERROR)
        !
        ! DIMENSION OF           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
        ! ARGUMENTS              USOL(IDMN, N+1), GRHS(IDMN, N+1),
        !
        ! LATEST REVISION        JUNE 2004
        !
        ! PURPOSE                SEPELI SOLVES FOR EITHER THE SECOND-ORDER
        !                        FINITE DIFFERENCE APPROXIMATION OR A
        !                        FOURTH-ORDER APPROXIMATION TO A SEPARABLE
        !                        ELLIPTIC EQUATION
        !
        !                                 2    2
        !                          AF(X)*D U/DX + BF(X)*DU/DX  + CF(X)*U +
        !                                 2    2
        !                          DF(Y)*D U/DY  + EF(Y)*DU/DY + FF(Y)*U
        !
        !                          = G(X, Y)
        !
        !                        ON A RECTANGLE (X GREATER THAN OR EQUAL TO A
        !                        AND LESS THAN OR EQUAL TO B; Y GREATER THAN
        !                        OR EQUAL TO C AND LESS THAN OR EQUAL TO D).
        !                        ANY COMBINATION OF PERIODIC OR MIXED BOUNDARY
        !                        CONDITIONS IS ALLOWED.
        !
        !                        THE POSSIBLE BOUNDARY CONDITIONS ARE:
        !                        IN THE X-DIRECTION:
        !                        (0) PERIODIC, U(X+B-A, Y)=U(X, Y) FOR ALL
        !                            Y, X (1) U(A, Y), U(B, Y) ARE SPECIFIED FOR
        !                            ALL Y
        !                        (2) U(A, Y), DU(B, Y)/DX+BETA*U(B, Y) ARE
        !                            SPECIFIED FOR ALL Y
        !                        (3) DU(A, Y)/DX+ALPHA*U(A, Y), DU(B, Y)/DX+
        !                            BETA*U(B, Y) ARE SPECIFIED FOR ALL Y
        !                        (4) DU(A, Y)/DX+ALPHA*U(A, Y), U(B, Y) ARE
        !                            SPECIFIED FOR ALL Y
        !
        !                        IN THE Y-DIRECTION:
        !                        (0) PERIODIC, U(X, Y+D-C)=U(X, Y) FOR ALL X, Y
        !                        (1) U(X, C), U(X, D) ARE SPECIFIED FOR ALL X
        !                        (2) U(X, C), DU(X, D)/DY+XNU*U(X, D) ARE
        !                            SPECIFIED FOR ALL X
        !                        (3) DU(X, C)/DY+GAMA*U(X, C), DU(X, D)/DY+
        !                            XNU*U(X, D) ARE SPECIFIED FOR ALL X
        !                        (4) DU(X, C)/DY+GAMA*U(X, C), U(X, D) ARE
        !                            SPECIFIED FOR ALL X
        !
        ! USAGE                  CALL SEPELI (INTL, IORDER, A, B, M, MBDCND, BDA,
        !                                     ALPHA, BDB, BETA, C, D, N, NBDCND, BDC,
        !                                     GAMA, BDD, XNU, COFX, COFY, GRHS, USOL,
        !                                     IDMN, W, PERTRB, IERROR)
        !
        ! ARGUMENTS
        ! ON INPUT               INTL
        !                          = 0 ON INITIAL ENTRY TO SEPELI OR IF ANY
        !                              OF THE ARGUMENTS C, D, N, NBDCND, COFY
        !                              ARE CHANGED FROM A PREVIOUS CALL
        !                          = 1 IF C, D, N, NBDCND, COFY ARE UNCHANGED
        !                              FROM THE PREVIOUS CALL.
        !
        !                        IORDER
        !                          = 2 IF A SECOND-ORDER APPROXIMATION
        !                              IS SOUGHT
        !                          = 4 IF A FOURTH-ORDER APPROXIMATION
        !                              IS SOUGHT
        !
        !                        A, B
        !                          THE RANGE OF THE X-INDEPENDENT VARIABLE,
        !                          I.E., X IS GREATER THAN OR EQUAL TO A
        !                          AND LESS THAN OR EQUAL TO B.  A MUST BE
        !                          LESS THAN B.
        !
        !                        M
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL [A, B] IS SUBDIVIDED. HENCE,
        !                          THERE WILL BE M+1 GRID POINTS IN THE X-
        !                          DIRECTION GIVEN BY XI=A+(I-1)*DLX
        !                          FOR I=1, 2, ..., M+1 WHERE DLX=(B-A)/M IS
        !                          THE PANEL WIDTH.  M MUST BE LESS THAN
        !                          IDMN AND GREATER THAN 5.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITION
        !                          AT X=A AND X=B
        !
        !                          = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.,
        !                              U(X+B-A, Y)=U(X, Y)  FOR ALL Y, X
        !                          = 1 IF THE SOLUTION IS SPECIFIED AT X=A
        !                              AND X=B, I.E., U(A, Y) AND U(B, Y) ARE
        !                              SPECIFIED FOR ALL Y
        !                          = 2 IF THE SOLUTION IS SPECIFIED AT X=A AND
        !                              THE BOUNDARY CONDITION IS MIXED AT X=B,
        !                              I.E., U(A, Y) AND DU(B, Y)/DX+BETA*U(B, Y)
        !                              ARE SPECIFIED FOR ALL Y
        !                          = 3 IF THE BOUNDARY CONDITIONS AT X=A AND
        !                              X=B ARE MIXED, I.E.,
        !                              DU(A, Y)/DX+ALPHA*U(A, Y) AND
        !                              DU(B, Y)/DX+BETA*U(B, Y) ARE SPECIFIED
        !                              FOR ALL Y
        !                          = 4 IF THE BOUNDARY CONDITION AT X=A IS
        !                              MIXED AND THE SOLUTION IS SPECIFIED
        !                              AT X=B, I.E., DU(A, Y)/DX+ALPHA*U(A, Y)
        !                              AND U(B, Y) ARE SPECIFIED FOR ALL Y
        !
        !                        BDA
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
        !                          THAT SPECIFIES THE VALUES OF
        !                          DU(A, Y)/DX+ ALPHA*U(A, Y) AT X=A, WHEN
        !                          MBDCND=3 OR 4.
        !                          BDA(J) = DU(A, YJ)/DX+ALPHA*U(A, YJ),
        !                          J=1, 2, ..., N+1. WHEN MBDCND HAS ANY OTHER
        !                          OTHER VALUE, BDA IS A DUMMY PARAMETER.
        !
        !                        ALPHA
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN
        !                          CASE OF A MIXED BOUNDARY CONDITION AT X=A
        !                          (SEE ARGUMENT BDA).  IF MBDCND IS NOT
        !                          EQUAL TO 3 OR 4 THEN ALPHA IS A DUMMY
        !                          PARAMETER.
        !
        !                        BDB
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1
        !                          THAT SPECIFIES THE VALUES OF
        !                          DU(B, Y)/DX+ BETA*U(B, Y) AT X=B.
        !                          WHEN MBDCND=2 OR 3
        !                          BDB(J) = DU(B, YJ)/DX+BETA*U(B, YJ),
        !                          J=1, 2, ..., N+1. WHEN MBDCND HAS ANY OTHER
        !                          OTHER VALUE, BDB IS A DUMMY PARAMETER.
        !
        !                        BETA
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN
        !                          CASE OF A MIXED BOUNDARY CONDITION AT
        !                          X=B (SEE ARGUMENT BDB).  IF MBDCND IS
        !                          NOT EQUAL TO 2 OR 3 THEN BETA IS A DUMMY
        !                          PARAMETER.
        !
        !                        C, D
        !                          THE RANGE OF THE Y-INDEPENDENT VARIABLE,
        !                          I.E., Y IS GREATER THAN OR EQUAL TO C
        !                          AND LESS THAN OR EQUAL TO D.  C MUST BE
        !                          LESS THAN D.
        !
        !                        N
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL [C, D] IS SUBDIVIDED.
        !                          HENCE, THERE WILL BE N+1 GRID POINTS
        !                          IN THE Y-DIRECTION GIVEN BY
        !                          YJ=C+(J-1)*DLY FOR J=1, 2, ..., N+1 WHERE
        !                          DLY=(D-C)/N IS THE PANEL WIDTH.
        !                          IN ADDITION, N MUST BE GREATER THAN 4.
        !
        !                        NBDCND
        !                          INDICATES THE TYPES OF BOUNDARY CONDITIONS
        !                          AT Y=C AND Y=D
        !
        !                          = 0 IF THE SOLUTION IS PERIODIC IN Y,
        !                              I.E., U(X, Y+D-C)=U(X, Y)  FOR ALL X, Y
        !                          = 1 IF THE SOLUTION IS SPECIFIED AT Y=C
        !                              AND Y = D, I.E., U(X, C) AND U(X, D)
        !                              ARE SPECIFIED FOR ALL X
        !                          = 2 IF THE SOLUTION IS SPECIFIED AT Y=C
        !                              AND THE BOUNDARY CONDITION IS MIXED
        !                              AT Y=D, I.E., U(X, C) AND
        !                              DU(X, D)/DY+XNU*U(X, D) ARE SPECIFIED
        !                              FOR ALL X
        !                          = 3 IF THE BOUNDARY CONDITIONS ARE MIXED
        !                              AT Y=C AND Y=D, I.E.,
        !                              DU(X, D)/DY+GAMA*U(X, C) AND
        !                              DU(X, D)/DY+XNU*U(X, D) ARE SPECIFIED
        !                              FOR ALL X
        !                          = 4 IF THE BOUNDARY CONDITION IS MIXED
        !                              AT Y=C AND THE SOLUTION IS SPECIFIED
        !                              AT Y=D, I.E. DU(X, C)/DY+GAMA*U(X, C)
        !                              AND U(X, D) ARE SPECIFIED FOR ALL X
        !
        !                        BDC
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
        !                          THAT SPECIFIES THE VALUE OF
        !                          DU(X, C)/DY+GAMA*U(X, C) AT Y=C.
        !                          WHEN NBDCND=3 OR 4 BDC(I) = DU(XI, C)/DY +
        !                          GAMA*U(XI, C), I=1, 2, ..., M+1.
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDC
        !                          IS A DUMMY PARAMETER.
        !
        !                        GAMA
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN
        !                          CASE OF A MIXED BOUNDARY CONDITION AT
        !                          Y=C (SEE ARGUMENT BDC).  IF NBDCND IS
        !                          NOT EQUAL TO 3 OR 4 THEN GAMA IS A DUMMY
        !                          PARAMETER.
        !
        !                        BDD
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1
        !                          THAT SPECIFIES THE VALUE OF
        !                          DU(X, D)/DY + XNU*U(X, D) AT Y=C.
        !                          WHEN NBDCND=2 OR 3 BDD(I) = DU(XI, D)/DY +
        !                          XNU*U(XI, D), I=1, 2, ..., M+1.
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDD
        !                          IS A DUMMY PARAMETER.
        !
        !                        XNU
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN
        !                          CASE OF A MIXED BOUNDARY CONDITION AT
        !                          Y=D (SEE ARGUMENT BDD).  IF NBDCND IS
        !                          NOT EQUAL TO 2 OR 3 THEN XNU IS A
        !                          DUMMY PARAMETER.
        !
        !                        COFX
        !                          A USER-SUPPLIED SUBPROGRAM WITH
        !                          PARAMETERS X, AFUN, BFUN, CFUN WHICH
        !                          RETURNS THE VALUES OF THE X-DEPENDENT
        !                          COEFFICIENTS AF(X), BF(X), CF(X) IN THE
        !                          ELLIPTIC EQUATION AT X.
        !
        !                        COFY
        !                          A USER-SUPPLIED SUBPROGRAM WITH PARAMETERS
        !                          Y, DFUN, EFUN, FFUN WHICH RETURNS THE
        !                          VALUES OF THE Y-DEPENDENT COEFFICIENTS
        !                          DF(Y), EF(Y), FF(Y) IN THE ELLIPTIC
        !                          EQUATION AT Y.
        !
        !                          NOTE:  COFX AND COFY MUST BE DECLARED
        !                          EXTERNAL IN THE CALLING ROUTINE.
        !                          THE VALUES RETURNED IN AFUN AND DFUN
        !                          MUST SATISFY AFUN*DFUN GREATER THAN 0
        !                          FOR A LESS THAN X LESS THAN B, C LESS
        !                          THAN Y LESS THAN D (SEE IERROR=10).
        !                          THE COEFFICIENTS PROVIDED MAY LEAD TO A
        !                          MATRIX EQUATION WHICH IS NOT DIAGONALLY
        !                          DOMINANT IN WHICH CASE SOLUTION MAY FAIL
        !                          (SEE IERROR=4).
        !
        !                        GRHS
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUES OF THE RIGHT-HAND SIDE OF THE
        !                          ELLIPTIC EQUATION, I.E.,
        !                          GRHS(I, J)=G(XI, YI), FOR I=2, ..., M,
        !                          J=2, ..., N.  AT THE BOUNDARIES, GRHS IS
        !                          DEFINED BY
        !
        !                          MBDCND   GRHS(1, J)   GRHS(M+1, J)
        !                          ------   ---------   -----------
        !                            0      G(A, YJ)     G(B, YJ)
        !                            1         *           *
        !                            2         *        G(B, YJ)  J=1, 2, ..., N+1
        !                            3      G(A, YJ)     G(B, YJ)
        !                            4      G(A, YJ)        *
        !
        !                          NBDCND   GRHS(I, 1)   GRHS(I, N+1)
        !                          ------   ---------   -----------
        !                            0      G(XI, C)     G(XI, D)
        !                            1         *           *
        !                            2         *        G(XI, D)  I=1, 2, ..., M+1
        !                            3      G(XI, C)     G(XI, D)
        !                            4      G(XI, C)        *
        !
        !                          WHERE * MEANS THESE QUANTITIES ARE NOT USED.
        !                          GRHS SHOULD BE DIMENSIONED IDMN BY AT LEAST
        !                          N+1 IN THE CALLING ROUTINE.
        !
        !                        USOL
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUES OF THE SOLUTION ALONG THE BOUNDARIES.
        !                          AT THE BOUNDARIES, USOL IS DEFINED BY
        !
        !                          MBDCND   USOL(1, J)   USOL(M+1, J)
        !                          ------   ---------   -----------
        !                            0         *           *
        !                            1      U(A, YJ)     U(B, YJ)
        !                            2      U(A, YJ)        *     J=1, 2, ..., N+1
        !                            3         *           *
        !                            4         *        U(B, YJ)
        !
        !                          NBDCND   USOL(I, 1)   USOL(I, N+1)
        !                          ------   ---------   -----------
        !                            0         *           *
        !                            1      U(XI, C)     U(XI, D)
        !                            2      U(XI, C)        *     I=1, 2, ..., M+1
        !                            3         *           *
        !                            4         *        U(XI, D)
        !
        !                          WHERE * MEANS THE QUANTITIES ARE NOT USED
        !                          IN THE SOLUTION.
        !
        !                          IF IORDER=2, THE USER MAY EQUIVALENCE GRHS
        !                          AND USOL TO SAVE SPACE.  NOTE THAT IN THIS
        !                          CASE THE TABLES SPECIFYING THE BOUNDARIES
        !                          OF THE GRHS AND USOL ARRAYS DETERMINE THE
        !                          BOUNDARIES UNIQUELY EXCEPT AT THE CORNERS.
        !                          IF THE TABLES CALL FOR BOTH G(X, Y) AND
        !                          U(X, Y) AT A CORNER THEN THE SOLUTION MUST
        !                          BE CHOSEN.  FOR EXAMPLE, IF MBDCND=2 AND
        !                          NBDCND=4, THEN U(A, C), U(A, D), U(B, D) MUST
        !                          BE CHOSEN AT THE CORNERS IN ADDITION
        !                          TO G(B, C).
        !
        !                          IF IORDER=4, THEN THE TWO ARRAYS, USOL AND
        !                          GRHS, MUST BE DISTINCT.
        !
        !                          USOL SHOULD BE DIMENSIONED IDMN BY AT LEAST
        !                          N+1 IN THE CALLING ROUTINE.
        !
        !                        IDMN
        !                          THE ROW (OR FIRST) DIMENSION OF THE ARRAYS
        !                          GRHS AND USOL AS IT APPEARS IN THE PROGRAM
        !                          CALLING SEPELI.  THIS PARAMETER IS USED
        !                          TO SPECIFY THE VARIABLE DIMENSION OF GRHS
        !                          AND USOL.  IDMN MUST BE AT LEAST 7 AND
        !                          GREATER THAN OR EQUAL TO M+1.
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
        !                          user program calling SEPELI.  The second statement
        !                          declares a derived type variable (defined in
        !                          the module "fish.f") which is used internally
        !                          in SEPELI to dynamically allocate real and complex
        !                          work space used in solution.  An error flag
        !                          (IERROR = 20) is set if the required work space
        !                          allocation fails (for example if N, M are too large)
        !                          Real and complex values are set in the components
        !                          of W on a initial (INTL=0) call to SEPELI.  These
        !                          must be preserved on non-initial calls (INTL=1)
        !                          to SEPELI.  This eliminates redundant calculations
        !                          and saves compute time.
        !               ****       IMPORTANT!  The user program calling SEPELI should
        !                          include the statement:
        !
        !                               CALL FISHFIN(W)
        !
        !                          after the final approximation is generated by
        !                          SEPELI.  The will deallocate the real and complex
        !                          work space of W.  Failure to include this statement
        !                          could result in serious memory leakage.
        !
        ! ON OUTPUT              USOL
        !                          CONTAINS THE APPROXIMATE SOLUTION TO THE
        !                          ELLIPTIC EQUATION.
        !                          USOL(I, J) IS THE APPROXIMATION TO U(XI, YJ)
        !                          FOR I=1, 2..., M+1 AND J=1, 2, ..., N+1.
        !                          THE APPROXIMATION HAS ERROR
        !                          O(DLX**2+DLY**2) IF CALLED WITH IORDER=2
        !                          AND O(DLX**4+DLY**4) IF CALLED WITH
        !                          IORDER=4.
        !
        !                        W
        !                          The derived type (FishpackWorkspace) variable W
        !                          contains real and complex values that must not
        !                          be destroyed if SEPELI is called again with
        !                          INTL=1.
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
        !                          BOUNDARY CONDITIONS
        !                          (I.E., ALPHA=BETA=0 IF MBDCND=3;
        !                          GAMA=XNU=0 IF NBDCND=3) IS SPECIFIED
        !                          AND IF THE COEFFICIENTS OF U(X, Y) IN THE
        !                          SEPARABLE ELLIPTIC EQUATION ARE ZERO
        !                          (I.E., CF(X)=0 FOR X GREATER THAN OR EQUAL
        !                          TO A AND LESS THAN OR EQUAL TO B;
        !                          FF(Y)=0 FOR Y GREATER THAN OR EQUAL TO C
        !                          AND LESS THAN OR EQUAL TO D) THEN A
        !                          SOLUTION MAY NOT EXIST.  PERTRB IS A
        !                          CONSTANT CALCULATED AND SUBTRACTED FROM
        !                          THE RIGHT-HAND SIDE OF THE MATRIX EQUATIONS
        !                          GENERATED BY SEPELI WHICH INSURES THAT A
        !                          SOLUTION EXISTS. SEPELI THEN COMPUTES THIS
        !                          SOLUTION WHICH IS A WEIGHTED MINIMAL LEAST
        !                          SQUARES SOLUTION TO THE ORIGINAL PROBLEM.
        !
        !                        IERROR
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS OR FAILURE TO FIND A SOLUTION
        !                          = 0 NO ERROR
        !                          = 1 IF A GREATER THAN B OR C GREATER THAN D
        !                          = 2 IF MBDCND LESS THAN 0 OR MBDCND GREATER
        !                              THAN 4
        !                          = 3 IF NBDCND LESS THAN 0 OR NBDCND GREATER
        !                              THAN 4
        !                          = 4 IF ATTEMPT TO FIND A SOLUTION FAILS.
        !                              (THE LINEAR SYSTEM GENERATED IS NOT
        !                              DIAGONALLY DOMINANT.)
        !                          = 5 IF IDMN IS TOO SMALL
        !                              (SEE DISCUSSION OF IDMN)
        !                          = 6 IF M IS TOO SMALL OR TOO LARGE
        !                              (SEE DISCUSSION OF M)
        !                          = 7 IF N IS TOO SMALL (SEE DISCUSSION OF N)
        !                          = 8 IF IORDER IS NOT 2 OR 4
        !                          = 9 IF INTL IS NOT 0 OR 1
        !                          = 10 IF AFUN*DFUN LESS THAN OR EQUAL TO 0
        !                               FOR SOME INTERIOR MESH POINT (XI, YJ)
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space in the derived type
        !                               (FishpackWorkspace) variable W fails (e.g.,
        !                               if N, M are too large for the platform used)
        !
        !                          NOTE (CONCERNING IERROR=4):  FOR THE
        !                          COEFFICIENTS INPUT THROUGH COFX, COFY,
        !                          THE DISCRETIZATION MAY LEAD TO A BLOCK
        !                          TRIDIAGONAL LINEAR SYSTEM WHICH IS NOT
        !                          DIAGONALLY DOMINANT (FOR EXAMPLE, THIS
        !                          HAPPENS IF CFUN=0 AND BFUN/(2.*DLX) GREATER
        !                          THAN AFUN/DLX**2).  IN THIS CASE SOLUTION
        !                          MAY FAIL.  THIS CANNOT HAPPEN IN THE LIMIT
        !                          AS DLX, DLY APPROACH ZERO.  HENCE, THE
        !                          CONDITION MAY BE REMEDIED BY TAKING LARGER
        !                          VALUES FOR M OR N.
        !
        ! SPECIAL CONDITIONS     SEE COFX, COFY ARGUMENT DESCRIPTIONS ABOVE.
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED FILES         blktri.f, comf.f, sepaux.f, fish.f
        !
        ! LANGUAGE               Fortran 90
        !
        ! HISTORY                DEVELOPED AT NCAR DURING 1975-76 BY
        !                        JOHN C. ADAMS OF THE SCIENTIFIC COMPUTING
        !                        DIVISION.  RELEASED ON NCAR'S PUBLIC SOFTWARE
        !                        LIBRARIES IN JANUARY 1980. Revised in June
        !                        2004 using Fortan 90 dynamically allocated work
        !                        space and derived data types to eliminate mixed
        !                        mode conflicts in the earlier versions. All
        !                        statement labels, arithmetic if statements and
        !                        computed GO TO statements have been removed from
        !                        the current version of SEPELI.
        !
        ! ALGORITHM              SEPELI AUTOMATICALLY DISCRETIZES THE
        !                        SEPARABLE ELLIPTIC EQUATION WHICH IS THEN
        !                        SOLVED BY A GENERALIZED CYCLIC REDUCTION
        !                        ALGORITHM IN THE SUBROUTINE, BLKTRI.  THE
        !                        FOURTH-ORDER SOLUTION IS OBTAINED USING
        !                        'DEFERRED CORRECTIONS' WHICH IS DESCRIBED
        !                        AND REFERENCED IN SECTIONS, REFERENCES AND
        !                        METHOD.
        !
        ! TIMING                 THE OPERATIONAL COUNT IS PROPORTIONAL TO
        !                        M*N*LOG2(N).
        !
        ! ACCURACY               THE FOLLOWING ACCURACY RESULTS WERE OBTAINED
        !                        using 64 bit floating point arithmetic.  Note
        !                        THAT THE FOURTH-ORDER accuracy is not realized
        !                        UNTIL THE MESH IS sufficiently refined.
        !
        !                                     SECOND-ORDER  FOURTH-ORDER
        !                            M    N     ERROR         ERROR
        !
        !                             6    6    6.8E-1        1.2E0
        !                            14   14    1.4E-1        1.8E-1
        !                            30   30    3.2E-2        9.7E-3
        !                            62   62    7.5E-3        3.0E-4
        !                           126  126    1.8E-3        3.5E-6
        !
        !
        ! REFERENCES             KELLER, H.B., NUMERICAL METHODS FOR TWO-POINT
        !                        BOUNDARY-VALUE PROBLEMS, BLAISDEL (1968),
        !                        WALTHAM, MASS.
        !
        !                        SWARZTRAUBER, P., AND R. SWEET (1975):
        !                        EFFICIENT FORTRAN SUBPROGRAMS FOR THE
        !                        SOLUTION OF ELLIPTIC PARTIAL DIFFERENTIAL
        !                        EQUATIONS.  NCAR TECHNICAL NOTE
        !                        NCAR-TN/IA-109, PP. 135-137.
        !***********************************************************************

        type (FishpackWorkspace) :: w
        external cofx, cofy
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer  :: intl
        integer  :: iorder
        integer  :: m
        integer  :: mbdcnd
        integer  :: n
        integer  :: nbdcnd
        integer  :: idmn
        integer  :: ierror
        real  :: a
        real  :: b
        real  :: alpha
        real  :: beta
        real  :: c
        real  :: d
        real  :: gama
        real  :: xnu
        real  :: pertrb
        real  :: bda(:)
        real  :: bdb(:)
        real  :: bdc(:)
        real  :: bdd(:)
        real  :: grhs(idmn, *)
        real  :: usol(idmn, *)
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer::i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, k, l, np, irwk, icwk

        save I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12
        !-----------------------------------------------
        !   E x t e r n a l   F u n c t i o n s
        !-----------------------------------------------
        !-----------------------------------------------
        !     save local variable work space pointers for noninitial call
        !     check input arguments
        call CHKPRM (intl, iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, cofx &
            , cofy, idmn, ierror)
        if (ierror /= 0) return
        if (intl == 0) then
            !     allocate space and set work space indices on initial call only
            k = m + 1
            l = n + 1
            !          compute required blktri work space lengths
            np = nbdcnd

            call w%get_block_tridiagonal_workpace_dimensions (n, m, irwk, icwk)
            !
            !     SET WORK SPACE INDICES
            !
            i1 = irwk + 1
            i2 = i1 + l
            i3 = i2 + l
            i4 = i3 + l
            i5 = i4 + l
            i6 = i5 + l
            i7 = i6 + l
            i8 = i7 + k
            i9 = i8 + k
            i10 = i9 + k
            i11 = i10 + k
            i12 = i11 + k
            !          set sepeli work space requirements
            irwk = i12 + k
            icwk = icwk + 3*(m + 1)
            !          allocate required real and complex work space
            call w%create( irwk, icwk, ierror )
            !          return if allocation failure
            if (ierror == 20) return
        end if
        ierror = 0
        !     compute second or fourth order solution
        call SPELIP (intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, d, n, &
            nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, w%rew(i1), w%rew(i2), w%rew(i3), &
            w%rew(i4), w%rew(i5), w%rew(i6), w%rew(i7), w%rew(i8), w%rew(i9), &
            w%rew(i10), w%rew(i11), w%rew(i12), grhs, usol, idmn, w%rew, w%cxw, &
            pertrb, ierror)

    end subroutine SEPELI


    subroutine SPELIP(intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, &
        beta, c, d, n, nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, an, bn &
        , cn, dn, un, zn, am, bm, cm, dm, um, zm, grhs, usol, idmn, w, &
        wc, pertrb, ierror)

        external cofx, cofy
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer  :: intl
        integer , intent (in) :: iorder
        integer , intent (in) :: m
        integer  :: mbdcnd
        integer , intent (in) :: n
        integer  :: nbdcnd
        integer  :: idmn
        integer  :: ierror
        real , intent (in) :: a
        real , intent (in) :: b
        real  :: alpha
        real  :: beta
        real , intent (in) :: c
        real , intent (in) :: d
        real  :: gama
        real  :: xnu
        real  :: pertrb
        real , intent (in) :: bda(*)
        real , intent (in) :: bdb(*)
        real , intent (in) :: bdc(*)
        real , intent (in) :: bdd(*)
        real  :: an(*)
        real  :: bn(*)
        real  :: cn(*)
        real  :: dn(*)
        real  :: un(*)
        real  :: zn(*)
        real  :: am(*)
        real  :: bm(*)
        real  :: cm(*)
        real  :: dm(*)
        real  :: um(*)
        real  :: zm(*)
        real  :: grhs(idmn, *)
        real  :: usol(idmn, *)
        real  :: w(*)
        complex  :: wc(*)
        !-----------------------------------------------
        !   C o m m o n   B l o c k s
        !-----------------------------------------------
        !...  /SPLP/
        common /SPLP/ kswx, kswy, k, l, ait, bit, cit, dit, mit, nit, is, &
            ms, js, ns, dlx, dly, tdlx3, tdly3, dlx4, dly4
        integer   kswx, kswy, k, l, mit, nit, is, ms, js, ns
        real   ait, bit, cit, dit, dlx, dly, tdlx3, tdly3, dlx4, dly4
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: i, j, i1, mp, np
        real :: xi, ai, bi, ci, axi, bxi, cxi, yj, dj, ej, fj, dyj, eyj, &
            fyj, ax1, cxm, dy1, fyn, prtrb
        logical :: singlr
        !-----------------------------------------------
        !   E x t e r n a l   F u n c t i o n s
        !-----------------------------------------------
        !-----------------------------------------------
        !
        !     SPELIP SETS UP VECTORS AND ARRAYS FOR INPUT TO BLKTRI
        !     AND COMPUTES A SECOND ORDER SOLUTION IN USOL.  A RETURN JUMP TO
        !     SEPELI OCCURRS IF IORDER=2.  IF IORDER=4 A FOURTH ORDER
        !     SOLUTION IS GENERATED IN USOL.
        !
        !
        !     SET PARAMETERS INTERNALLY
        !
        kswx = mbdcnd + 1
        kswy = nbdcnd + 1
        k = m + 1
        l = n + 1
        ait = a
        bit = b
        cit = c
        dit = d
        !
        !     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
        !     AND NON-SPECIFIED BOUNDARIES.
        !
        usol(2:m, 2:n) = GRHS(2:m, 2:n)
        if (kswx/=2 .and. kswx/=3) then
            usol(1, 2:n) = GRHS(1, 2:n)
        end if
        if (kswx/=2 .and. kswx/=5) then
            usol(k, 2:n) = GRHS(k, 2:n)
        end if
        if (kswy/=2 .and. kswy/=3) then
            usol(2:m, 1) = GRHS(2:m, 1)
        end if
        if (kswy/=2 .and. kswy/=5) then
            usol(2:m, l) = GRHS(2:m, l)
        end if
        if (kswx/=2 .and. kswx/=3 .and. kswy/=2 .and. kswy/=3) usol(1, 1) &
            = GRHS(1, 1)
        if (kswx/=2 .and. kswx/=5 .and. kswy/=2 .and. kswy/=3) usol(k, 1) &
            = GRHS(k, 1)
        if (kswx/=2 .and. kswx/=3 .and. kswy/=2 .and. kswy/=5) usol(1, l) &
            = GRHS(1, l)
        if (kswx/=2 .and. kswx/=5 .and. kswy/=2 .and. kswy/=5) usol(k, l) &
            = GRHS(k, l)
        i1 = 1
        !
        !     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
        !
        mp = 1
        np = 1
        if (kswx == 1) mp = 0
        if (kswy == 1) np = 0
        !
        !     SET DLX, DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
        !     IN NINT, MINT
        !
        dlx = (bit - ait)/real(m)
        mit = k - 1
        if (kswx == 2) mit = k - 2
        if (kswx == 4) mit = k
        dly = (dit - cit)/real(n)
        nit = l - 1
        if (kswy == 2) nit = l - 2
        if (kswy == 4) nit = l
        tdlx3 = 2.0*dlx**3
        dlx4 = dlx**4
        tdly3 = 2.0*dly**3
        dly4 = dly**4
        !
        !     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
        !
        is = 1
        js = 1
        if (kswx==2 .or. kswx==3) is = 2
        if (kswy==2 .or. kswy==3) js = 2
        ns = nit + js - 1
        ms = mit + is - 1
        !
        !     SET X - DIRECTION
        !
        do i = 1, mit
            xi = ait + real(is + i - 2)*dlx
            call COFX (xi, ai, bi, ci)
            axi = (ai/dlx - 0.5*bi)/dlx
            bxi = (-2.*ai/dlx**2) + ci
            cxi = (ai/dlx + 0.5*bi)/dlx
            am(i) = axi
            bm(i) = bxi
            cm(i) = cxi
        end do
        !
        !     SET Y DIRECTION
        !
        do j = 1, nit
            yj = cit + real(js + j - 2)*dly
            call COFY (yj, dj, ej, fj)
            dyj = (dj/dly - 0.5*ej)/dly
            eyj = (-2.*dj/dly**2) + fj
            fyj = (dj/dly + 0.5*ej)/dly
            an(j) = dyj
            bn(j) = eyj
            cn(j) = fyj
        end do
        !
        !     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
        !
        ax1 = AM(1)
        cxm = CM(mit)
        select case (kswx)
            case (2)
                !
                !     DIRICHLET-DIRICHLET IN X DIRECTION
                !
                am(1) = 0.0
                cm(mit) = 0.0
            case (5)
                !
                !     MIXED-DIRICHLET IN X DIRECTION
                !
                am(1) = 0.0
                bm(1) = BM(1) + 2.*alpha*dlx*ax1
                cm(1) = CM(1) + ax1
                cm(mit) = 0.0
            case (3)
                !
                !     DIRICHLET-MIXED IN X DIRECTION
                !
                am(1) = 0.0
                am(mit) = AM(mit) + cxm
                bm(mit) = BM(mit) - 2.*beta*dlx*cxm
                cm(mit) = 0.0
            !
            !     MIXED - MIXED IN X DIRECTION
            !
            case (4)
                am(1) = 0.0
                bm(1) = BM(1) + 2.*dlx*alpha*ax1
                cm(1) = CM(1) + ax1
                am(mit) = AM(mit) + cxm
                bm(mit) = BM(mit) - 2.*dlx*beta*cxm
                cm(mit) = 0.0
        end select
        !
        !     ADJUST IN Y DIRECTION UNLESS PERIODIC
        !
        dy1 = AN(1)
        fyn = CN(nit)
        select case (kswy)
            case (2)
                !
                !     DIRICHLET-DIRICHLET IN Y DIRECTION
                !
                an(1) = 0.0
                cn(nit) = 0.0
            case (5)
                !
                !     MIXED-DIRICHLET IN Y DIRECTION
                !
                an(1) = 0.0
                bn(1) = BN(1) + 2.*dly*gama*dy1
                cn(1) = CN(1) + dy1
                cn(nit) = 0.0
            case (3)
                !
                !     DIRICHLET-MIXED IN Y DIRECTION
                !
                an(1) = 0.0
                an(nit) = AN(nit) + fyn
                bn(nit) = BN(nit) - 2.*dly*xnu*fyn
                cn(nit) = 0.0
            case (4)
                !
                !     MIXED - MIXED DIRECTION IN Y DIRECTION
                !
                an(1) = 0.0
                bn(1) = BN(1) + 2.*dly*gama*dy1
                cn(1) = CN(1) + dy1
                an(nit) = AN(nit) + fyn
                bn(nit) = BN(nit) - 2.0*dly*xnu*fyn
                cn(nit) = 0.0
        end select
        if (kswx /= 1) then
            !
            !     ADJUST USOL ALONG X EDGE
            !
            if (kswx==2 .or. kswx==3) then
                if (kswx==2 .or. kswx==5) then
                    usol(is, js:ns) = USOL(is, js:ns) - ax1*USOL(1, js:ns)
                    usol(ms, js:ns) = USOL(ms, js:ns) - cxm*USOL(k, js:ns)
                else
                    usol(is, js:ns) = USOL(is, js:ns) - ax1*USOL(1, js:ns)
                    usol(ms, js:ns) = USOL(ms, js:ns) - 2.0*dlx*cxm*BDB(js:ns)
                end if
            else
                if (kswx==2 .or. kswx==5) then
                    usol(is, js:ns) = USOL(is, js:ns) + 2.0*dlx*ax1*BDA(js:ns)
                    usol(ms, js:ns) = USOL(ms, js:ns) - cxm*USOL(k, js:ns)
                else
                    usol(is, js:ns) = USOL(is, js:ns) + 2.0*dlx*ax1*BDA(js:ns)
                    usol(ms, js:ns) = USOL(ms, js:ns) - 2.0*dlx*cxm*BDB(js:ns)
                end if
            end if
        end if
        if (kswy /= 1) then
            !
            !     ADJUST USOL ALONG Y EDGE
            !
            if (kswy==2 .or. kswy==3) then
                if (kswy==2 .or. kswy==5) then
                    usol(is:ms, js) = USOL(is:ms, js) - dy1*USOL(is:ms, 1)
                    usol(is:ms, ns) = USOL(is:ms, ns) - fyn*USOL(is:ms, l)
                else
                    usol(is:ms, js) = USOL(is:ms, js) - dy1*USOL(is:ms, 1)
                    usol(is:ms, ns) = USOL(is:ms, ns) - 2.0*dly*fyn*BDD(is:ms)
                end if
            else
                if (kswy==2 .or. kswy==5) then
                    usol(is:ms, js) = USOL(is:ms, js) + 2.0*dly*dy1*BDC(is:ms)
                    usol(is:ms, ns) = USOL(is:ms, ns) - fyn*USOL(is:ms, l)
                else
                    usol(is:ms, js) = USOL(is:ms, js) + 2.0*dly*dy1*BDC(is:ms)
                    usol(is:ms, ns) = USOL(is:ms, ns) - 2.0*dly*fyn*BDD(is:ms)
                end if
            end if
        end if
        !
        !     SAVE ADJUSTED EDGES IN GRHS IF IORDER=4
        !
        if (iorder == 4) then
            grhs(is, js:ns) = USOL(is, js:ns)
            grhs(ms, js:ns) = USOL(ms, js:ns)
            grhs(is:ms, js) = USOL(is:ms, js)
            grhs(is:ms, ns) = USOL(is:ms, ns)
        end if
        pertrb = 0.0
        !
        !     CHECK IF OPERATOR IS SINGULAR
        !
        call CHKSNG(mbdcnd, nbdcnd, alpha, beta, gama, xnu, cofx, cofy, singlr)
        !
        !     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
        !     IF SINGULAR
        !
        if (singlr) call SEPTRI (mit, am, bm, cm, dm, um, zm)
        if (singlr) call SEPTRI (nit, an, bn, cn, dn, un, zn)
        !
        !     MAKE INITIALIZATION CALL TO blktrii
        !
        if (intl == 0) then
            call BLKTRII (intl, np, nit, an, bn, cn, mp, mit, am, bm, cm, &
                idmn, USOL(is, js), ierror, w, wc)
            if (ierror /= 0) return
        end if
        !
        !     ADJUST RIGHT HAND SIDE IF NECESSARY
        !
        if (singlr) call SEPORT (usol, idmn, zn, zm, pertrb)
        !
        !     COMPUTE SOLUTION
        !
        call BLKTRII (i1, np, nit, an, bn, cn, mp, mit, am, bm, cm, idmn, &
            USOL(is, js), ierror, w, wc)
        if (ierror /= 0) return
        !
        !     SET PERIODIC BOUNDARIES IF NECESSARY
        !
        if (kswx == 1) then
            usol(k, :l) = USOL(1, :l)
        end if
        if (kswy == 1) then
            usol(:k, l) = USOL(:k, 1)
        end if
        !
        !     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
        !     NORM IF OPERATOR IS SINGULAR
        !
        if (singlr) call SEPMIN (usol, idmn, zn, zm, prtrb)
        !
        !     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
        !     NOT FLAGGED
        !
        if (iorder == 2) return
        !
        !     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
        !
        call DEFER (cofx, cofy, idmn, usol, grhs)
        if (singlr) call SEPORT (usol, idmn, zn, zm, pertrb)
        !
        !     COMPUTE fourth order SOLUTION
        !
        call BLKTRII (i1, np, nit, an, bn, cn, mp, mit, am, bm, cm, idmn, &
            USOL(is, js), ierror, w, wc)
        if (ierror /= 0) return
        !
        !     SET PERIODIC BOUNDARIES IF NECESSARY
        !
        if (kswx == 1) then
            usol(k, :l) = USOL(1, :l)
        end if
        if (kswy == 1) then
            usol(:k, l) = USOL(:k, 1)
        end if
        !
        !     MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
        !     NORM IF OPERATOR IS SINGULAR
        !
        if (singlr) call SEPMIN (usol, idmn, zn, zm, prtrb)

    end subroutine SPELIP

    subroutine CHKPRM(intl, iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, &
        cofx, cofy, idmn, ierror)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: intl
        integer , intent (in) :: iorder
        integer , intent (in) :: m
        integer , intent (in) :: mbdcnd
        integer , intent (in) :: n
        integer , intent (in) :: nbdcnd
        integer , intent (in) :: idmn
        integer , intent (out) :: ierror
        real , intent (in) :: a
        real , intent (in) :: b
        real , intent (in) :: c
        real , intent (in) :: d
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: i, j
        real :: dlx, dly, xi, ai, bi, ci, yj, dj, ej, fj
        !-----------------------------------------------
        !   E x t e r n a l   F u n c t i o n s
        !-----------------------------------------------
        !-----------------------------------------------
        !
        !     THIS PROGRAM CHECKS THE INPUT arguments FOR ERRORS
        !
        !
        !     CHECK DEFINITION OF SOLUTION REGION
        !
        if (a>=b .or. c>=d) then
            ierror = 1
            return
        end if
        !
        !     check boundary condition arguments
        !
        if (mbdcnd<0 .or. mbdcnd>4) then
            ierror = 2
            return
        end if
        if (nbdcnd<0 .or. nbdcnd>4) then
            ierror = 3
            return
        end if
        !
        !     CHECK FIRST DIMENSION IN CALLING ROUTINE
        !
        if (idmn < 7) then
            ierror = 5
            return
        end if
        !
        !     CHECK M, N
        !
        if (m>idmn - 1 .or. m<6) then
            ierror = 6
            return
        end if
        if (n < 5) then
            ierror = 7
            return
        end if
        !
        !     CHECK IORDER
        !
        if (iorder/=2 .and. iorder/=4) then
            ierror = 8
            return
        end if
        !
        !     CHECK INTL
        !
        if (intl/=0 .and. intl/=1) then
            ierror = 9
            return
        end if
        !
        !     CHECK THAT EQUATION IS ELLIPTIC (only on initial call)
        !
        if (intl == 0) then
            dlx = (b - a)/real(m)
            dly = (d - c)/real(n)
            do i = 2, m
                xi = a + real(i - 1)*dlx
                call COFX (xi, ai, bi, ci)
                do j = 2, n
                    yj = c + real(j - 1)*dly
                    call COFY (yj, dj, ej, fj)
                    if (ai*dj > 0.0) cycle
                    ierror = 10
                    return
                end do
            end do
        end if
        !
        !     NO ERROR FOUND
        !
        ierror = 0

    end subroutine CHKPRM

    subroutine CHKSNG(mbdcnd, nbdcnd, alpha, beta, gama, xnu, cofx, &
        cofy, singlr)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: mbdcnd
        integer , intent (in) :: nbdcnd
        real , intent (in) :: alpha
        real , intent (in) :: beta
        real , intent (in) :: gama
        real , intent (in) :: xnu
        logical , intent (out) :: singlr
        !-----------------------------------------------
        !   C o m m o n   B l o c k s
        !-----------------------------------------------
        !...  /SPLP/
        common /SPLP/ kswx, kswy, k, l, ait, bit, cit, dit, mit, nit, is, &
            ms, js, ns, dlx, dly, tdlx3, tdly3, dlx4, dly4
        integer   kswx, kswy, k, l, mit, nit, is, ms, js, ns
        real   ait, bit, cit, dit, dlx, dly, tdlx3, tdly3, dlx4, dly4
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: i, j
        real :: xi, ai, bi, ci, yj, dj, ej, fj
        !-----------------------------------------------
        !   E x t e r n a l   F u n c t i o n s
        !-----------------------------------------------
        !-----------------------------------------------
        !
        !     THIS SUBROUTINE CHECKS IF THE PDE   SEPELI
        !     MUST SOLVE IS A SINGULAR OPERATOR
        !
        singlr = .false.
        !
        !     CHECK IF THE BOUNDARY CONDITIONS ARE
        !     ENTIRELY PERIODIC AND/OR MIXED
        !
        if(mbdcnd/=0.and.mbdcnd/=3.or.nbdcnd/=0.and.nbdcnd/=3)return
        !
        !     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
        !
        if (mbdcnd == 3) then
            if (alpha/=0.0 .or. beta/=0.0) return
        end if

        if (nbdcnd == 3) then
            if (gama/=0.0 .or. xnu/=0.0) return
        end if
        !
        !     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
        !     ARE ZERO
        !
        do i = is, ms
            xi = ait + real(i - 1)*dlx
            call COFX (xi, ai, bi, ci)
            if (ci == 0.0) cycle
            return
        end do
        do j = js, ns
            yj = cit + real(j - 1)*dly
            call COFY (yj, dj, ej, fj)
            if (fj == 0.0) cycle
            return
        end do
        !
        !     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
        !
        singlr = .true.

    end subroutine CHKSNG

    subroutine DEFER(cofx, cofy, idmn, usol, grhs)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer  :: idmn
        real  :: usol(idmn, *)
        real , intent (in out) :: grhs(idmn, *)
        !-----------------------------------------------
        !   C o m m o n   B l o c k s
        !-----------------------------------------------
        !...  /SPLP/
        common /SPLP/ kswx, kswy, k, l, ait, bit, cit, dit, mit, nit, is, &
            ms, js, ns, dlx, dly, tdlx3, tdly3, dlx4, dly4
        integer   kswx, kswy, k, l, mit, nit, is, ms, js, ns
        real   ait, bit, cit, dit, dlx, dly, tdlx3, tdly3, dlx4, dly4
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: j, i
        real::yj, dj, ej, fj, xi, ai, bi, ci, uxxx, uxxxx, uyyy, uyyyy, tx, ty
        !-----------------------------------------------
        !   E x t e r n a l   F u n c t i o n s
        !-----------------------------------------------
        !-----------------------------------------------
        !
        !     THIS SUBROUTINE FIRST APPROXIMATES THE TRUNCATION ERROR GIVEN BY
        !     TRUN1(X, Y)=DLX**2*TX+DLY**2*TY WHERE
        !     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 ON THE INTERIOR AND
        !     AT THE BOUNDARIES IF PERIODIC(HERE UXXX, UXXXX ARE THE THIRD
        !     AND FOURTH PARTIAL DERIVATIVES OF U WITH RESPECT TO X).
        !     TX IS OF THE FORM AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
        !     AT X=A OR X=B IF THE BOUNDARY CONDITION THERE IS MIXED.
        !     TX=0.0 ALONG SPECIFIED BOUNDARIES.  TY HAS SYMMETRIC FORM
        !     IN Y WITH X, AFUN(X), BFUN(X) REPLACED BY Y, DFUN(Y), EFUN(Y).
        !     THE SECOND ORDER SOLUTION IN USOL IS USED TO APPROXIMATE
        !     (VIA SECOND ORDER FINITE DIFFERENCING) THE TRUNCATION ERROR
        !     AND THE RESULT IS ADDED TO THE RIGHT HAND SIDE IN GRHS
        !     AND THEN TRANSFERRED TO USOL TO BE USED AS A NEW RIGHT
        !     HAND SIDE WHEN CALLING BLKTRI FOR A FOURTH ORDER SOLUTION.
        !
        !
        !     COMPUTE TRUNCATION ERROR APPROXIMATION OVER THE ENTIRE MESH
        !
        do j = js, ns
            yj = cit + real(j - 1)*dly
            call COFY (yj, dj, ej, fj)
            do i = is, ms
                xi = ait + real(i - 1)*dlx
                call COFX (xi, ai, bi, ci)
                !
                !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI, YJ)
                !
                call SEPDX (usol, idmn, i, j, uxxx, uxxxx)
                call SEPDY (usol, idmn, i, j, uyyy, uyyyy)
                tx = ai*uxxxx/12.0 + bi*uxxx/6.0
                ty = dj*uyyyy/12.0 + ej*uyyy/6.0
                !
                !     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
                !
                if (kswx/=1 .and. (i==1 .or. i==k)) tx = ai/3.0*(uxxxx/4.0 &
                    + uxxx/dlx)
                if (kswy/=1 .and. (j==1 .or. j==l)) ty = dj/3.0*(uyyyy/4.0 &
                    + uyyy/dly)
                grhs(i, j) = GRHS(i, j) + dlx**2*tx + dly**2*ty
            end do
        end do
        !
        !     RESET THE RIGHT HAND SIDE IN USOL
        !
        usol(is:ms, js:ns) = GRHS(is:ms, js:ns)

    end subroutine DEFER

end module module_sepeli
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    version 5.0, fortran 90 changes
!-----------------------------------------------------------------------
