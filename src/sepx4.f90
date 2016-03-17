module module_sepx4

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_genbun, only: &
        genbun

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
    public :: sepx4
    public :: sepx4_unit_test

contains

    subroutine sepx4_unit_test()
        !
        !     file tsepx4.f
        !
        !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        !     *                                                               *
        !     *                  copyright(c) 2005 by UCAR                   *
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
        !     *                Boulder, Colorado (80307)  U.S.A.             *
        !     *                                                               *
        !     *                   which is sponsored by                       *
        !     *                                                               *
        !     *              the National Science Foundation                  *
        !     *                                                               *
        !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        !
        !--------------------------------------------------------------------------------
        ! Dictionary
        !--------------------------------------------------------------------------------
        integer(ip) ::m, n, nx, ny, i, j, mbdcnd, nbdcnd, idmn, iorder, ierror
        real (wp), dimension(33, 33) :: usol, grhs
        real (wp), dimension(33) :: bda, bdb
        real (wp) :: a, b, c, d, dlx, dly, x, af, bf, cf, y, alpha, beta, dum(1), &
            pertrb, err, err2, err4
        !-----------------------------------------------

        !
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
        dlx =(b - a)/real(m)
        dly =(d - c)/real(n)
        nx = m + 1
        ny = n + 1
        do i = 1, nx
            x = a + real(i - 1)*dlx
            !
            !     SET SPECIFIED BOUNDARY CONDITIONS AT Y=C, D
            !
            usol(i, 1) = UE(x, c)
            usol(i, ny) = UE(x, d)
            call get_coefficients_in_x_direction(x, af, bf, cf)
            do j = 1, ny
                y = c + real(j - 1)*dly
                !
                !     SET RIGHT HAND SIDE
                !
                grhs(i, j)=af*UXXE(x, y)+bf*UXE(x, y)+cf*UE(x, y)+UYYE(x, y)
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
        !     SET FIRST DIMENSION OF USOL, GRHS AND WORK SPACE LENGTH
        !
        idmn = 33
        !
        !     OBTAIN SECOND ORDER APPROXIMATION
        !
        iorder = 2

        call SEPX4(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, d, &
            n, nbdcnd, dum, dum, get_coefficients_in_x_direction, grhs, usol, idmn, pertrb, ierror)
        !
        !     COMPUTE SECOND ORDER DISCRETIZATION ERROR(RELATIVE)
        !     ALSO RESET SPECIFIED BOUNDARIES AND RIGHT HAND SIDE.
        !
        err = 0.0
        do i = 1, nx
            x = a + real(i - 1)*dlx
            usol(i, 1) = UE(x, c)
            usol(i, ny) = UE(x, d)
            call get_coefficients_in_x_direction(x, af, bf, cf)
            do j = 1, ny
                y = c + real(j - 1)*dly
                err = max(err, abs((USOL(i, j)-UE(x, y))/UE(x, y)))
                !
                !     RESET RIGHT HAND SIDE IN GRHS FOR FOURTH ORDER APPROXIMATION CALL
                !
                grhs(i, j)=af*UXXE(x, y)+bf*UXE(x, y)+cf*UE(x, y)+UYYE(x, y)
            end do
        end do
        err2 = err
        !
        !     OBTAIN FOURTH ORDER APPROXIMATION
        !
        iorder = 4
        call SEPX4(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, d, &
            n, nbdcnd, dum, dum, get_coefficients_in_x_direction, grhs, usol, idmn, pertrb, ierror)
        !
        !     COMPUTE FOURTH ORDER DISCRETIZATION ERROR(RELATIVE)
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

        write( *, *) ''
        write( *, *) '    SEPEX4 TEST RUN *** '
        write( *, *) &
            '    Previous 64 bit floating point arithmetic result '
        write( *, *) '    IERROR = 0'
        write( *, *) '    Second Order Discretization Error = 1.5985E-4'
        write( *, *) '    Fourth Order Discretization Error = 1.8575E-6'

        write( *, *) '    The output from your computer is: '
        write( *, *) '    IERROR =', ierror
        write( *, *) '    Second Order Discretization Error =', err2
        write( *, *) '    Fourth Order Discretization Error =', err4

    contains
        !
        !*****************************************************************************************
        !
        pure function UE(s, t) result( return_value )
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            real, intent (in) :: s
            real, intent (in) :: t
            real (wp)          :: return_value

            return_value =(s*t)**3 + 1.0_wp

        end function UE
        !
        !*****************************************************************************************
        !
        pure function UXE(s, t) result( return_value )
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            real, intent (in) :: s
            real, intent (in) :: t
            real (wp)          :: return_value

            return_value = 3.0_wp *(s**2)*(t**3)

        end function UXE
        !
        !*****************************************************************************************
        !
        pure function UXXE(s, t) result( return_value )
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            real, intent (in) :: s
            real, intent (in) :: t
            real (wp)          :: return_value

            return_value = 6.0_wp * s *(t**3)

        end function UXXE
        !
        !*****************************************************************************************
        !
        pure function UYE(s, t) result( return_value )
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            real, intent (in) :: s
            real, intent (in) :: t
            real (wp) :: return_value

            return_value = 3.0_wp *(s**3) *(t**2)

        end function UYE
        !
        !*****************************************************************************************
        !
        pure function UYYE(s, t) result( return_value )
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            real, intent (in) :: s
            real, intent (in) :: t
            real (wp) :: return_value

            return_value = 6.0_wp *(s**3) * t

        end function UYYE
        !
        !*****************************************************************************************
        !
        subroutine get_coefficients_in_x_direction(x, af, bf, cf)
            !--------------------------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------------------------
            real (wp), intent (in) :: x
            real (wp), intent (out) :: af
            real (wp), intent (out) :: bf
            real (wp), intent (out) :: cf
            !-----------------------------------------------
            !
            !     SET COEFFICIENTS IN THE X-DIRECTION.
            !
            af =(x + 1.)**2
            bf = 2.0*(x + 1.)
            cf = -x

        end subroutine get_coefficients_in_x_direction


    end subroutine sepx4_unit_test


    subroutine SEPX4(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, &
        d, n, nbdcnd, bdc, bdd, cofx, grhs, usol, idmn, pertrb, &
        ierror)
        !
        !     file sepx4.f
        !
        !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        !     *                                                               *
        !     *                  copyright(c) 2005 by UCAR                   *
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
        !     *                Boulder, Colorado (80307)  U.S.A.             *
        !     *                                                               *
        !     *                   which is sponsored by                       *
        !     *                                                               *
        !     *              the National Science Foundation                  *
        !     *                                                               *
        !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        !
        !     SUBROUTINE SEPX4(IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA, C, D, N,
        !    +                  NBDCND, BDC, BDD, COFX, GRHS, USOL, IDMN, PERTRB,
        !    +                  IERROR)
        !
        !
        !
        ! DIMENSION OF           BDA(N+1), BDB(N+1), BDC(M+1), BDD(M+1),
        ! ARGUMENTS              USOL(IDMN, N+1),     GRHS(IDMN, N+1),
        !
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SEPX4 SOLVES FOR EITHER THE SECOND-ORDER
        !                        FINITE DIFFERENCE APPROXIMATION OR A
        !                        FOURTH-ORDER APPROXIMATION TO A SEPARABLE
        !                        ELLIPTIC EQUATION
        !
        !                          AF(X)*UXX+BF(X)*UX+CF(X)*U+UYY = G(X, Y)
        !
        !                        ON A RECTANGLE(X GREATER THAN OR EQUAL TO
        !                        A AND LESS THAN OR EQUAL TO B, Y GREATER THAN
        !                        OR EQUAL TO C AND LESS THAN OR EQUAL TO D).
        !                        ANY COMBINATION OF PERIODIC OR MIXED BOUNDARY
        !                        CONDITIONS IS ALLOWED.  IF BOUNDARY
        !                        CONDITIONS IN THE X DIRECTION ARE PERIODIC
        !                       (SEE MBDCND=0 BELOW) THEN THE COEFFICIENTS
        !                        MUST SATISFY
        !
        !                          AF(X)=C1, BF(X)=0, CF(X)=C2 FOR ALL X.
        !
        !                        HERE C1, C2 ARE CONSTANTS, C1.GT.0.
        !
        !                        THE POSSIBLE BOUNDARY CONDITIONS ARE:
        !                        IN THE X-DIRECTION:
        !                         (0) PERIODIC, U(X+B-A, Y)=U(X, Y) FOR
        !                              ALL Y, X
        !                         (1) U(A, Y), U(B, Y) ARE SPECIFIED FOR ALL Y
        !                         (2) U(A, Y), DU(B, Y)/DX+BETA*U(B, Y) ARE
        !                              SPECIFIED FOR ALL Y
        !                         (3) DU(A, Y)/DX+ALPHA*U(A, Y), DU(B, Y)/DX+
        !                              BETA*U(B, Y) ARE SPECIFIED FOR ALL Y
        !                         (4) DU(A, Y)/DX+ALPHA*U(A, Y), U(B, Y) ARE
        !                              SPECIFIED FOR ALL Y
        !
        !                        IN THE Y-DIRECTION:
        !                         (0) PERIODIC, U(X, Y+D-C)=U(X, Y) FOR ALL X, Y
        !                         (1) U(X, C), U(X, D) ARE SPECIFIED FOR ALL X
        !                         (2) U(X, C), DU(X, D)/DY ARE SPECIFIED FOR
        !                              ALL X
        !                         (3) DU(X, C)/DY, DU(X, D)/DY ARE SPECIFIED FOR
        !                              ALL X
        !                         (4) DU(X, C)/DY, U(X, D) ARE SPECIFIED FOR
        !                              ALL X
        !
        ! USAGE                  CALL SEPX4(IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB,
        !                                   BETA, C, D, N, NBDCND, BDC, BDD, COFX,
        !                                   GRHS, USOL, IDMN, W, PERTRB, IERROR)
        !
        ! ARGUMENTS
        ! ON INPUT               IORDER
        !                          = 2 IF A SECOND-ORDER APPROXIMATION IS
        !                              SOUGHT
        !                          = 4 IF A FOURTH-ORDER APPROXIMATION IS
        !                              SOUGHT
        !
        ! *** caution ***          GRHS SHOULD BE RESET IF SEPX4 WAS FIRST CALLED
        !                          WITH IORDER=2 AND WILL BE CALLED AGAIN WITH
        !                          IORDER=4.  VALUES IN GRHS ARE DESTROYED BY THE
        !                          IORDER=2 CALL.
        !
        !
        !                        A, B
        !                          THE RANGE OF THE X-INDEPENDENT VARIABLE,
        !                          I.E., X IS GREATER THAN OR EQUAL TO A
        !                          AND LESS THAN OR EQUAL TO B.  A MUST BE
        !                          LESS THAN B.
        !
        !                        M
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL(A, B) IS SUBDIVIDED.  HENCE,
        !                          THERE WILL BE M+1 GRID POINTS IN THE X-
        !                          DIRECTION GIVEN BY XI=A+(I-1)*DLX
        !                          FOR I=1, 2, ..., M+1 WHERE DLX=(B-A)/M IS
        !                          THE PANEL WIDTH.  M MUST BE LESS THAN
        !                          IDMN AND GREATER THAN 5.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITION
        !                          AT X=A AND X=B
        !                          = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.,
        !                              U(X+B-A, Y)=U(X, Y)  FOR ALL Y, X
        !                          = 1 IF THE SOLUTION IS SPECIFIED AT X=A
        !                              AND X=B, I.E., U(A, Y) AND U(B, Y) ARE
        !                              SPECIFIED FOR ALL Y
        !                          = 2 IF THE SOLUTION IS SPECIFIED AT X=A
        !                              AND THE BOUNDARY CONDITION IS MIXED AT
        !                              X=B, I.E., U(A, Y) AND
        !                              DU(B, Y)/DX+BETA*U(B, Y) ARE SPECIFIED
        !                              FOR ALL Y
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
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALUES OF
        !                          DU(A, Y)/DX+ ALPHA*U(A, Y) AT X=A, WHEN
        !                          MBDCND=3 OR 4.
        !                          BDA(J) = DU(A, YJ)/DX+ALPHA*U(A, YJ),
        !                          J=1, 2, ..., N+1
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDA IS
        !                          A DUMMY PARAMETER.
        !
        !                        ALPHA
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN CASE
        !                          OF A MIXED BOUNDARY CONDITION AT X=A
        !                         (SEE ARGUMENT BDA).  IF MBDCND IS NOT EQUAL
        !                          TO EITHER 3 OR 4, THEN ALPHA IS A DUMMY
        !                          PARAMETER.
        !
        !                        BDB
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH N+1 THAT
        !                          SPECIFIES THE VALUES OF
        !                          DU(B, Y)/DX+ BETA*U(B, Y) AT X=B.
        !                          WHEN MBDCND=2 OR 3
        !                          BDB(J) = DU(B, YJ)/DX+BETA*U(B, YJ),
        !                          J=1, 2, ..., N+1
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDB IS
        !                          A DUMMY PARAMETER.
        !
        !                        BETA
        !                          THE SCALAR MULTIPLYING THE SOLUTION IN
        !                          CASE OF A MIXED BOUNDARY CONDITION AT X=B
        !                         (SEE ARGUMENT BDB).  IF MBDCND IS NOT EQUAL
        !                          TO 2 OR 3, THEN BETA IS A DUMMY PARAMETER.
        !
        !                        C, D
        !                          THE RANGE OF THE Y-INDEPENDENT VARIABLE,
        !                          I.E., Y IS GREATER THAN OR EQUAL TO C AND
        !                          LESS THAN OR EQUAL TO D.  C MUST BE LESS
        !                          THAN D.
        !
        !                        N
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL(C, D) IS SUBDIVIDED.  HENCE,
        !                          THERE WILL BE N+1 GRID POINTS IN THE Y-
        !                          DIRECTION GIVEN BY YJ=C+(J-1)*DLY FOR
        !                          J=1, 2, ..., N+1 WHERE DLY=(D-C)/N IS THE
        !                          PANEL WIDTH.  IN ADDITION, N MUST BE
        !                          GREATER THAN 4.
        !
        !                        NBDCND
        !                          INDICATES THE TYPES OF BOUNDARY CONDITIONS
        !                          AT Y=C AND Y=D
        !                          = 0 IF THE SOLUTION IS PERIODIC IN Y,
        !                              I.E., U(X, Y+D-C)=U(X, Y) FOR ALL X, Y
        !                          = 1 IF THE SOLUTION IS SPECIFIED AT Y=C
        !                              AND Y = D, I.E., U(X, C)  AND U(X, D)
        !                              ARE SPECIFIED FOR ALL X
        !                          = 2 IF THE SOLUTION IS SPECIFIED AT Y=C
        !                              AND THE BOUNDARY CONDITION IS MIXED
        !                              AT Y=D, I.E., DU(X, C)/DY AND U(X, D)
        !                              ARE SPECIFIED FOR ALL X
        !                          = 3 IF THE BOUNDARY CONDITIONS ARE MIXED
        !                              AT Y=CAND Y=D I.E.,
        !                              DU(X, D)/DY AND DU(X, D)/DY ARE
        !                              SPECIFIED FOR ALL X
        !                          = 4 IF THE BOUNDARY CONDITION IS MIXED
        !                              AT Y=C AND THE SOLUTION IS SPECIFIED
        !                              AT Y=D, I.E. DU(X, C)/DY+GAMA*U(X, C)
        !                              AND U(X, D) ARE SPECIFIED FOR ALL X
        !
        !                        BDC
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIES THE VALUE DU(X, C)/DY AT Y=C.
        !
        !                          WHEN NBDCND=3 OR 4
        !                            BDC(I) = DU(XI, C)/DY I=1, 2, ..., M+1.
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDC IS
        !                          A DUMMY PARAMETER.
        !
        !                        BDD
        !                          A ONE-DIMENSIONAL ARRAY OF LENGTH M+1 THAT
        !                          SPECIFIED THE VALUE OF DU(X, D)/DY AT Y=D.
        !
        !                          WHEN NBDCND=2 OR 3
        !                            BDD(I)=DU(XI, D)/DY I=1, 2, ..., M+1.
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDD IS
        !                          A DUMMY PARAMETER.
        !
        !                        COFX
        !                          A USER-SUPPLIED SUBPROGRAM WITH PARAMETERS
        !                          X, AFUN, BFUN, CFUN WHICH RETURNS THE
        !                          VALUES OF THE X-DEPENDENT COEFFICIENTS
        !                          AF(X), BF(X), CF(X) IN THE ELLIPTIC
        !                          EQUATION AT X.  IF BOUNDARY CONDITIONS IN
        !                          THE X DIRECTION ARE PERIODIC THEN THE
        !                          COEFFICIENTS MUST SATISFY AF(X)=C1, BF(X)=0,
        !                          CF(X)=C2 FOR ALL X.  HERE C1.GT.0
        !                          AND C2 ARE CONSTANTS.
        !
        !                          NOTE THAT COFX MUST BE DECLARED EXTERNAL
        !                          IN THE CALLING ROUTINE.
        !
        !                        GRHS
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUES OF THE RIGHT-HAND SIDE OF THE
        !                          ELLIPTIC EQUATION, I.E., GRHS(I, J)=G(XI, YI),
        !                          FOR I=2, ..., M, J=2, ..., N.  AT THE
        !                          BOUNDARIES, GRHS IS DEFINED BY
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
        !                          WHERE * MEANS THESE QUANTITES ARE NOT USED.
        !                          GRHS SHOULD BE DIMENSIONED IDMN BY AT LEAST
        !                          N+1 IN THE CALLING ROUTINE.
        !
        ! *** caution              GRHS SHOULD BE RESET IF SEPX4 WAS FIRST CALLED
        !                          WITH IORDER=2 AND WILL BE CALLED AGAIN WITH
        !                          IORDER=4.  VALUES IN GRHS ARE DESTROYED BY THE
        !                          IORDER=2 CALL.
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
        !                          WHERE * MEANS THE QUANTITES ARE NOT USED
        !                          IN THE SOLUTION.
        !
        !                          IF IORDER=2, THE USER MAY EQUIVALENCE GRHS
        !                          AND USOL TO SAVE SPACE.  NOTE THAT IN THIS
        !                          CASE THE TABLES SPECIFYING THE BOUNDARIES
        !                          OF THE GRHS AND USOL ARRAYS DETERMINE THE
        !                          BOUNDARIES UNIQUELY EXCEPT AT THE CORNERS.
        !                          IF THE TABLES CALL FOR BOTH G(X, Y) AND
        !                          U(X, Y) AT A CORNER THEN THE SOLUTION MUST
        !                          BE CHOSEN.
        !                          FOR EXAMPLE, IF MBDCND=2 AND NBDCND=4,
        !                          THEN U(A, C), U(A, D), U(B, D) MUST BE CHOSEN
        !                          AT THE CORNERS IN ADDITION TO G(B, C).
        !
        !                          IF IORDER=4, THEN THE TWO ARRAYS, USOL AND
        !                          GRHS, MUST BE DISTINCT.
        !
        !                          USOL SHOULD BE DIMENSIONED IDMN BY AT LEAST
        !                          N+1 IN THE CALLING ROUTINE.
        !
        !                        IDMN
        !                          THE ROW(OR FIRST) DIMENSION OF THE ARRAYS
        !                          GRHS AND USOL AS IT APPEARS IN THE PROGRAM
        !                          CALLING SEPELI.  THIS PARAMETER IS USED
        !                          TO SPECIFY THE VARIABLE DIMENSION OF GRHS
        !                          AND USOL.  IDMN MUST BE AT LEAST 7 AND
        !                          GREATER THAN OR EQUAL TO M+1.
        !
        !
        ! ON OUTPUT              USOL
        !                          CONTAINS THE APPROXIMATE SOLUTION TO THE
        !                          ELLIPTIC EQUATION. USOL(I, J) IS THE
        !                          APPROXIMATION TO U(XI, YJ) FOR I=1, 2..., M+1
        !                          AND J=1, 2, ..., N+1.  THE APPROXIMATION HAS
        !                          ERROR O(DLX**2+DLY**2) IF CALLED WITH
        !                          IORDER=2 AND O(DLX**4+DLY**4) IF CALLED
        !                          WITH IORDER=4.
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
        !                          BOUNDARY CONDITIONS(I.E., ALPHA=BETA=0 IF
        !                          MBDCND=3) IS SPECIFIED AND IF CF(X)=0 FOR
        !                          ALL X THEN A SOLUTION TO THE DISCRETIZED
        !                          MATRIX EQUATION MAY NOT EXIST
        !                         (REFLECTING THE NON-UNIQUENESS OF SOLUTIONS
        !                          TO THE PDE).
        !                          PERTRB IS A CONSTANT CALCULATED AND
        !                          SUBTRACTED FROM THE RIGHT HAND SIDE OF THE
        !                          MATRIX EQUATION INSURING THE EXISTENCE OF A
        !                          SOLUTION.  SEPX4 COMPUTES THIS SOLUTION
        !                          WHICH IS A WEIGHTED MINIMAL LEAST SQUARES
        !                          SOLUTION TO THE ORIGINAL PROBLEM.  IF
        !                          SINGULARITY IS NOT DETECTED PERTRB=0.0 IS
        !                          RETURNED BY SEPX4.
        !
        !                        IERROR
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS OR FAILURE TO FIND A SOLUTION
        !
        !                          =  0 NO ERROR
        !                          =  1 IF A GREATER THAN B OR C GREATER
        !                               THAN D
        !                          =  2 IF MBDCND LESS THAN 0 OR MBDCND
        !                               GREATER THAN 4
        !                          =  3 IF NBDCND LESS THAN 0 OR NBDCND
        !                               GREATER THAN 4
        !                          =  4 IF ATTEMPT TO FIND A SOLUTION FAILS.
        !                              (THE LINEAR SYSTEM GENERATED IS NOT
        !                               DIAGONALLY DOMINANT.)
        !                          =  5 IF IDMN IS TOO SMALL(SEE DISCUSSION
        !                               OF IDMN)
        !                          =  6 IF M IS TOO SMALL OR TOO LARGE
        !                              (SEE DISCUSSION OF M)
        !                          =  7 IF N IS TOO SMALL(SEE DISCUSSION OF N)
        !                          =  8 IF IORDER IS NOT 2 OR 4
        !                          =  9 IF INTL IS NOT 0 OR 1
        !                          = 10 IF AFUN IS LESS THAN OR EQUAL TO ZERO
        !                               FOR SOME INTERIOR MESH POINT XI SOME
        !                               INTERIOR MESH POINT(XI, YJ)
        !                          = 12 IF MBDCND=0 AND AF(X)=CF(X)=CONSTANT
        !                               OR BF(X)=0 FOR ALL X IS NOT TRUE.
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails(for example if N, M are too large
        !                               for your computer)
        !
        ! SPECIAL CONDITIONS     NONE
        !
        ! I/O                    NONE
        !
        ! REQUIRED files         fish.f, comf.f, genbun.f, gnbnaux.f, sepaux.f
        !
        !
        ! PRECISION              SINGLE
        !
        !
        ! LANGUAGE               FORTRAN 90
        !
        ! HISTORY                SEPX4 WAS DEVELOPED AT NCAR BY JOHN C.
        !                        ADAMS OF THE SCIENTIFIC COMPUTING DIVISION
        !                        IN OCTOBER 1978.  THE BASIS OF THIS CODE IS
        !                        NCAR ROUTINE SEPELI.  BOTH PACKAGES WERE
        !                        RELEASED ON NCAR'S PUBLIC LIBRARIES IN
        !                        JANUARY 1980. SEPX4 was modified in June 2004
        !                        incorporating fortran 90 dynamical storage
        !                        allocation for work space requirements
        !
        ! PORTABILITY            FORTRAN 90
        !
        ! ALGORITHM              SEPX4 AUTOMATICALLY DISCRETIZES THE SEPARABLE
        !                        ELLIPTIC EQUATION WHICH IS THEN SOLVED BY A
        !                        GENERALIZED CYCLIC REDUCTION ALGORITHM IN THE
        !                        SUBROUTINE POIS.  THE FOURTH ORDER SOLUTION
        !                        IS OBTAINED USING THE TECHNIQUE OF DEFFERRED
        !                        CORRECTIONS REFERENCED BELOW.
        !
        ! TIMING                 WHEN POSSIBLE, SEPX4 SHOULD BE USED INSTEAD
        !                        OF PACKAGE SEPELI.  THE INCREASE IN SPEED
        !                        IS AT LEAST A FACTOR OF THREE.
        !
        ! REFERENCES             KELLER, H.B., NUMERICAL METHODS FOR TWO-POINT
        !                        BOUNDARY-VALUE PROBLEMS, BLAISDEL(1968),
        !                        WALTHAM, MASS.
        !
        !                        SWARZTRAUBER, P., AND R. SWEET(1975):
        !                        EFFICIENT FORTRAN SUBPROGRAMS FOR THE
        !                        SOLUTION OF ELLIPTIC PARTIAL DIFFERENTIAL
        !                        EQUATIONS.  NCAR TECHNICAL NOTE
        !                          NCAR-TN/IA-109, PP. 135-137.
        !***********************************************************************
        type(FishpackWorkspace) :: w
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
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
        real  :: pertrb
        real  :: bda(*)
        real  :: bdb(*)
        real  :: bdc(*)
        real  :: bdd(*)
        real  :: grhs(idmn, *)
        real  :: usol(idmn, *)
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer :: l, k, log2n, length, irwk, icwk, i1, i2, i3, i4, i5, i6 &
            , i7, i8, i9, i10, i11, i12, i13
        !-----------------------------------------------
        !   E x t e r n a l   F u n c t i o n s
        !-----------------------------------------------
        external cofx
        !-----------------------------------------------
        !
        !     CHECK INPUT PARAMETERS
        !
        call C4KPRM(iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, cofx, idmn, ierror)

        if (ierror /= 0) return
        !
        !     COMPUTE MINIMUM WORK SPACE AND CHECK WORK SPACE LENGTH INPUT
        !
        l = n + 1
        if (nbdcnd == 0) l = n
        k = m + 1
        l = n + 1
        !     ESTIMATE LOG BASE 2 OF N
        log2n = INT(log(real(n + 1))/log(2.0) + 0.5)
        !     set required work space estimate
        length = 4*(n + 1) +(10 + log2n)*(m + 1)
        irwk = length + 6*(k + l) + 1
        icwk = 0
        !     allocate work space
        call w%create( irwk, icwk, ierror )
        !     return if allocation failed
        if (ierror == 20) return
        ierror = 0
        !
        !     SET WORK SPACE INDICES
        !
        i1 = length + 1
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
        i13 = 1

        call S4ELIP(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, d, n, &
            nbdcnd, bdc, bdd, cofx, w%rew(i1), w%rew(i2), w%rew(i3), &
            w%rew(i4), w%rew(i5), w%rew(i6), w%rew(i7), w%rew(i8), &
            w%rew(i9), w%rew(i10), w%rew(i11), w%rew(i12), &
            grhs, usol, idmn, w%rew(i13), pertrb, ierror)

        !     release dynamically allocated work space
        call w%destroy()

    end subroutine SEPX4

    subroutine S4ELIP(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, &
        c, d, n, nbdcnd, bdc, bdd, cofx, an, bn, cn, dn, un, zn, am, bm &
        , cm, dm, um, zm, grhs, usol, idmn, w, pertrb, ierror)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: iorder
        integer , intent (in) :: m
        integer  :: mbdcnd
        integer , intent (in) :: n
        integer  :: nbdcnd
        integer  :: idmn
        integer , intent (in out) :: ierror
        real , intent (in) :: a
        real , intent (in) :: b
        real  :: alpha
        real  :: beta
        real , intent (in) :: c
        real , intent (in) :: d
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
        integer :: i, j, i1, mp, np, ieror
        real :: xi, ai, bi, ci, axi, bxi, cxi, dyj, eyj, fyj, ax1, cxm, &
            dy1, fyn, gama, xnu, prtrb
        logical :: singlr
        !-----------------------------------------------
        !   E x t e r n a l   F u n c t i o n s
        !-----------------------------------------------
        external COFX
        !-----------------------------------------------
        !
        !     S4ELIP SETS UP VECTORS AND ARRAYS FOR INPUT TO BLKTRI
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
        dly =(dit - cit)/real(n)
        !
        !     SET RIGHT HAND SIDE VALUES FROM GRHS IN USOL ON THE INTERIOR
        !     AND NON-SPECIFIED BOUNDARIES.
        !
        usol(2:m, 2:n) = dly**2*GRHS(2:m, 2:n)
        if (kswx/=2 .and. kswx/=3) then
            usol(1, 2:n) = dly**2*GRHS(1, 2:n)
        end if
        if (kswx/=2 .and. kswx/=5) then
            usol(k, 2:n) = dly**2*GRHS(k, 2:n)
        end if
        if (kswy/=2 .and. kswy/=3) then
            usol(2:m, 1) = dly**2*GRHS(2:m, 1)
        end if
        if (kswy/=2 .and. kswy/=5) then
            usol(2:m, l) = dly**2*GRHS(2:m, l)
        end if
        if (kswx/=2 .and. kswx/=3 .and. kswy/=2 .and. kswy/=3) usol(1, 1) &
            = dly**2*GRHS(1, 1)
        if (kswx/=2 .and. kswx/=5 .and. kswy/=2 .and. kswy/=3) usol(k, 1) &
            = dly**2*GRHS(k, 1)
        if (kswx/=2 .and. kswx/=3 .and. kswy/=2 .and. kswy/=5) usol(1, l) &
            = dly**2*GRHS(1, l)
        if (kswx/=2 .and. kswx/=5 .and. kswy/=2 .and. kswy/=5) usol(k, l) &
            = dly**2*GRHS(k, l)
        i1 = 1
        !
        !     SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
        !
        mp = 1
        if (kswx == 1) mp = 0
        np = nbdcnd
        !
        !     SET DLX, DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
        !     IN NINT, MINT
        !
        dlx =(bit - ait)/real(m)
        mit = k - 1
        if (kswx == 2) mit = k - 2
        if (kswx == 4) mit = k
        dly =(dit - cit)/real(n)
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
            call COFX(xi, ai, bi, ci)
            axi =(ai/dlx - 0.5*bi)/dlx
            bxi =(-2.*ai/dlx**2) + ci
            cxi =(ai/dlx + 0.5*bi)/dlx
            am(i) = dly**2*axi
            bm(i) = dly**2*bxi
            cm(i) = dly**2*cxi
        end do
        !
        !     SET Y DIRECTION
        !
        dyj = 1.0
        eyj = -2.0
        fyj = 1.0
        an(:nit) = dyj
        bn(:nit) = eyj
        cn(:nit) = fyj
        !
        !     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
        !
        ax1 = AM(1)
        cxm = CM(mit)
        select case(kswx)
            case(2)
                !
                !     DIRICHLET-DIRICHLET IN X DIRECTION
                !
                am(1) = 0.0
                cm(mit) = 0.0
            case(3)
                !
                !     DIRICHLET-MIXED IN X DIRECTION
                !
                am(1) = 0.0
                am(mit) = AM(mit) + cxm
                bm(mit) = BM(mit) - 2.*beta*dlx*cxm
                cm(mit) = 0.0
            case(4)
                !
                !     MIXED - MIXED IN X DIRECTION
                !
                am(1) = 0.0
                bm(1) = BM(1) + 2.*dlx*alpha*ax1
                cm(1) = CM(1) + ax1
                am(mit) = AM(mit) + cxm
                bm(mit) = BM(mit) - 2.*dlx*beta*cxm
                cm(mit) = 0.0
            case(5)
                !
                !     MIXED-DIRICHLET IN X DIRECTION
                !
                am(1) = 0.0
                bm(1) = BM(1) + 2.*alpha*dlx*ax1
                cm(1) = CM(1) + ax1
                cm(mit) = 0.0
        end select
        !
        !     ADJUST IN Y DIRECTION UNLESS PERIODIC
        !
        dy1 = AN(1)
        fyn = CN(nit)
        gama = 0.0
        xnu = 0.0
        select case(kswy)
            case(2)
                !
                !     DIRICHLET-DIRICHLET IN Y DIRECTION
                !
                an(1) = 0.0
                cn(nit) = 0.0
            case(3)
                !
                !     DIRICHLET-MIXED IN Y DIRECTION
                !
                an(1) = 0.0
                an(nit) = AN(nit) + fyn
                bn(nit) = BN(nit) - 2.*dly*xnu*fyn
                cn(nit) = 0.0
            case(4)
                !
                !     MIXED - MIXED DIRECTION IN Y DIRECTION
                !
                an(1) = 0.0
                bn(1) = BN(1) + 2.*dly*gama*dy1
                cn(1) = CN(1) + dy1
                an(nit) = AN(nit) + fyn
                bn(nit) = BN(nit) - 2.0*dly*xnu*fyn
                cn(nit) = 0.0
            case(5)
                !
                !     MIXED-DIRICHLET IN Y DIRECTION
                !
                an(1) = 0.0
                bn(1) = BN(1) + 2.*dly*gama*dy1
                cn(1) = CN(1) + dy1
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
        call C4KSNG(mbdcnd, nbdcnd, alpha, beta, cofx, singlr)
        !
        !     COMPUTE NON-ZERO EIGENVECTOR IN NULL SPACE OF TRANSPOSE
        !     IF SINGULAR
        !
        if (singlr) call SEPTRI(mit, am, bm, cm, dm, um, zm)
        if (singlr) call SEPTRI(nit, an, bn, cn, dn, un, zn)
        !
        !     ADJUST RIGHT HAND SIDE IF NECESSARY
        !
        if (singlr) call SEPORT(usol, idmn, zn, zm, pertrb)
        !
        !     COMPUTE SOLUTION
        !
        !     SAVE ADJUSTED RIGHT HAND SIDE IN GRHS
        grhs(is:ms, js:ns) = USOL(is:ms, js:ns)
        call GENBUN(np, nit, mp, mit, am, bm, cm, idmn, USOL(is, js), ieror)
        !
        !     CHECK IF ERROR DETECTED IN POIS
        !     THIS CAN ONLY CORRESPOND TO IERROR=12
        if (ieror /= 0) then
            !       SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS
            ierror = 12
            return
        end if
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
        if (singlr) call SEPMIN(usol, idmn, zn, zm, prtrb)
        !
        !     RETURN IF DEFERRED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
        !     NOT FLAGGED
        !
        if (iorder == 2) return
        !
        !     COMPUTE NEW RIGHT HAND SIDE FOR FOURTH ORDER SOLUTION
        !
        call D4FER(cofx, idmn, usol, grhs)
        if (singlr) call SEPORT(usol, idmn, zn, zm, pertrb)
        !
        !     COMPUTE SOLUTION
        !
        !     SAVE ADJUSTED RIGHT HAND SIDE IN GRHS
        grhs(is:ms, js:ns) = USOL(is:ms, js:ns)
        call GENBUN(np, nit, mp, mit, am, bm, cm, idmn, USOL(is, js), ieror)
        !     CHECK IF ERROR DETECTED IN POIS
        !     THIS CAN ONLY CORRESPOND TO IERROR=12
        if (ieror /= 0) then
            !       SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS
            ierror = 12
            return
        end if
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
        if (singlr) call SEPMIN(usol, idmn, zn, zm, prtrb)

    end subroutine S4ELIP

    subroutine C4KPRM(iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, cofx, &
        idmn, ierror)
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
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
        integer :: i
        real :: dlx, xi, ai, bi, ci
        !-----------------------------------------------
        !   E x t e r n a l   F u n c t i o n s
        !-----------------------------------------------
        !-----------------------------------------------
        !
        !     THIS PROGRAM CHECKS THE INPUT PARAMETERS FOR ERRORS
        !
        !
        !
        !     CHECK DEFINITION OF SOLUTION REGION
        !
        if (a>=b .or. c>=d) then
            ierror = 1
            return
        end if
        !
        !     CHECK BOUNDARY SWITCHES
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
        !     CHECK M
        !
        if (m>idmn - 1 .or. m<6) then
            ierror = 6
            return
        end if
        !
        !     CHECK N
        !
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
        !     CHECK THAT EQUATION IS ELLIPTIC
        !
        dlx =(b - a)/real(m)
        do i = 2, m
            xi = a + real(i - 1)*dlx
            call COFX(xi, ai, bi, ci)
            if (ai > 0.0) cycle
            ierror = 10
            return
        end do
        !
        !     NO ERROR FOUND
        !
        ierror = 0

    end subroutine C4KPRM

    subroutine C4KSNG(mbdcnd, nbdcnd, alpha, beta, cofx, singlr)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer , intent (in) :: mbdcnd
        integer , intent (in) :: nbdcnd
        real , intent (in) :: alpha
        real , intent (in) :: beta
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
        integer :: i
        real :: xi, ai, bi, ci
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
        if (mbdcnd/=0.and.mbdcnd/=3.or.nbdcnd/=0.and.nbdcnd/=3)return
        !
        !     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
        !
        if (mbdcnd == 3) then
            if (alpha/=0.0 .or. beta/=0.0) return
        end if
        !
        !     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
        !     ARE ZERO
        !
        do i = is, ms
            xi = ait + real(i - 1)*dlx
            call COFX(xi, ai, bi, ci)
            if (ci == 0.0) cycle
            return
        end do
        !
        !     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
        !
        singlr = .true.

    end subroutine C4KSNG

    subroutine D4FER(cofx, idmn, usol, grhs)

        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer                :: idmn
        real                   :: usol(idmn, *)
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
        integer :: i, j
        real :: xi, ai, bi, ci, uxxx, uxxxx, uyyy, uyyyy, tx, ty
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
        !    (VIA SECOND ORDER FINITE DIFFERENCING) THE TRUNCATION ERROR
        !     AND THE RESULT IS ADDED TO THE RIGHT HAND SIDE IN GRHS
        !     AND THEN TRANSFERRED TO USOL TO BE USED AS A NEW RIGHT
        !     HAND SIDE WHEN CALLING BLKTRI FOR A FOURTH ORDER SOLUTION.
        !
        !
        !
        !     COMPUTE TRUNCATION ERROR APPROXIMATION OVER THE ENTIRE MESH
        !
        do i = is, ms
            xi = ait + real(i - 1)*dlx
            call COFX(xi, ai, bi, ci)
            do j = js, ns
                !
                !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT(XI, YJ)
                !
                call SEPDX(usol, idmn, i, j, uxxx, uxxxx)
                call SEPDY(usol, idmn, i, j, uyyy, uyyyy)
                tx = ai*uxxxx/12.0 + bi*uxxx/6.0
                ty = uyyyy/12.0
                !
                !     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
                !
                if (kswx/=1 .and.(i==1 .or. i==k)) tx = ai/3.0*(uxxxx/4.0 &
                    + uxxx/dlx)
                if (kswy/=1.and.(j==1.or.j==l))ty=(uyyyy/4.0+uyyy/dly)/3.0
                grhs(i, j) = GRHS(i, j) + dly**2*(dlx**2*tx + dly**2*ty)
            end do
        end do
        !
        !     RESET THE RIGHT HAND SIDE IN USOL
        !
        usol(is:ms, js:ns) = GRHS(is:ms, js:ns)

    end subroutine D4FER

end module module_sepx4
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
