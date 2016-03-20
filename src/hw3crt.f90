module module_hw3crt

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_pois3d, only: &
        pois3dd

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hw3crt
    public :: test_hw3crt

contains

    subroutine test_hw3crt()
        !
        !     file thw3crt.f
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
        integer (ip) :: lbdcnd, mbdcnd, nbdcnd, l, m, n, ldimf, mdimf, lp1, i, &
            mp1, j, np1, k, ierror
        real (wp), dimension(11, 41, 16) :: f
        real (wp), dimension(11, 41) :: bdzf, bdxs, bdxf, bdys, bdyf, bdzs
        real (wp), dimension(11) :: x
        real (wp), dimension(41) :: y
        real (wp), dimension(16) :: z
        real (wp) :: elmbda, xs, xf, ys, pi, yf, zs, zf, dx, dy, dz, pertrb, discretization_error, t
        !-----------------------------------------------

        !
        !        FROM THE DESCRIPTION OF THE PROBLEM GIVEN ABOVE, WE DEFINE
        !     THE FOLLOWING QUANTITIES
        !
        elmbda = -3.
        xs = 0.
        xf = 1.
        lbdcnd = 1
        ys = 0.
        pi = acos( -1.0 )
        yf = 2.*pi
        mbdcnd = 0
        zs = 0.
        zf = pi/2.
        nbdcnd = 2
        l = 10
        m = 40
        n = 15
        !
        !     FROM THE DIMENSION STATEMENT ABOVE WE DEFINE
        !
        ldimf = 11
        mdimf = 41
        !
        !     WE DEFINE THE GRID POINTS FOR LATER USE.
        !
        lp1 = l + 1
        dx = (xf - xs)/real(l)
        do i = 1, lp1
            x(i) = xs + real(i - 1)*dx
        end do
        mp1 = m + 1
        dy = (yf - ys)/real(m)
        do j = 1, mp1
            y(j) = ys + real(j - 1)*dy
        end do
        np1 = n + 1
        dz = (zf - zs)/real(n)
        do k = 1, np1
            z(k) = zs + real(k - 1)*dz
        end do
        !
        !     WE DEFINE THE ARRAY OF DERIVATIVE BOUNDARY VALUES.
        !
        do i = 1, lp1
            do j = 1, mp1
                bdzf(i, j) = -X(i)**4*SIN(Y(j))
            end do
        end do
        !
        !     NOTE THAT FOR THIS EXAMPLE ALL OTHER BOUNDARY ARRAYS ARE
        !     DUMMY VARIABLES.
        !     WE DEFINE THE FUNCTION BOUNDARY VALUES IN THE F ARRAY.
        !
        do j = 1, mp1
            do k = 1, np1
                f(1, j, k) = 0.
                f(lp1, j, k) = SIN(Y(j))*COS(Z(k))
            end do
        end do
        do i = 1, lp1
            do j = 1, mp1
                f(i, j, 1) = X(i)**4*SIN(Y(j))
            end do
        end do
        !
        !     WE NOW DEFINE THE VALUES OF THE RIGHT SIDE OF THE HELMHOLTZ
        !     EQUATION.
        !
        do i = 2, l
            do j = 1, mp1
                do k = 2, np1
                    f(i, j, k) = 4.*X(i)**2*(3. - X(i)**2)*SIN(Y(j))*COS(Z(k))
                end do
            end do
        end do
        !
        !     CALL hw3crt TO GENERATE AND SOLVE THE FINITE DIFFERENCE EQUATION.
        !
        call hw3crt (xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, mbdcnd, &
            bdys, bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, ldimf, mdimf, &
            f, pertrb, ierror)
        !
        !     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION TO THE
        !     PROBLEM IS
        !
        !        U(X, Y, Z) = X**4*SIN(Y)*COS(Z)
        !
        discretization_error = 0.
        do i = 1, lp1
            do j = 1, mp1
                do k = 1, np1
                    t = abs(F(i, j, k)-X(i)**4*SIN(Y(j))*COS(Z(k)))
                    discretization_error = max(t, discretization_error)
                end do
            end do
        end do
        !     Print earlier output from platforms with 32 and 64 bit floating point
        !     arithemtic followed by the output from this computer
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     hw3crt *** TEST RUN *** '
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') '     ierror = 0,  discretization error = 9.6480E-3'
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)') &
            '     ierror =', ierror, ' discretization error = ', discretization_error

    end subroutine test_hw3crt



    subroutine hw3crt(xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, mbdcnd, &
        bdys, bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, ldimf, &
        mdimf, f, pertrb, ierror)
        !
        !     file hw3crt.f
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
        !     SUBROUTINE hw3crt (XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M, MBDCND, BDYS,
        !    +                   BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, ELMBDA, LDIMF,
        !    +                   MDIMF, F, PERTRB, ierror)
        !
        !
        ! DIMENSION OF           BDXS(MDIMF, N+1),    BDXF(MDIMF, N+1),
        ! ARGUMENTS              BDYS(LDIMF, N+1),    BDYF(LDIMF, N+1),
        !                        BDZS(LDIMF, M+1),    BDZF(LDIMF, M+1),
        !                        F(LDIMF, MDIMF, N+1)
        !
        ! LATEST REVISION        June 2004
        !
        ! PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
        !                        DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
        !                        EQUATION IN CARTESIAN COORDINATES.  THIS
        !                        EQUATION IS
        !
        !                          (D/DX)(DU/DX) + (D/DY)(DU/DY) +
        !                          (D/DZ)(DU/DZ) + LAMBDA*U = F(X, Y, Z) .
        !
        ! USAGE                  CALL hw3crt (XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M,
        !                                     MBDCND, BDYS, BDYF, ZS, ZF, N, NBDCND,
        !                                     BDZS, BDZF, ELMBDA, LDIMF, MDIMF, F,
        !                                     PERTRB, ierror)
        !
        ! ARGUMENTS
        !
        ! ON INPUT               XS, XF
        !
        !                          THE RANGE OF X, I.E. XS .LE. X .LE. XF .
        !                          XS MUST BE LESS THAN XF.
        !
        !                        L
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL (XS, XF) IS SUBDIVIDED.
        !                          HENCE, THERE WILL BE L+1 GRID POINTS
        !                          IN THE X-DIRECTION GIVEN BY
        !                          X(I) = XS+(I-1)DX FOR I=1, 2, ..., L+1,
        !                          WHERE DX = (XF-XS)/L IS THE PANEL WIDTH.
        !                          L MUST BE AT LEAST 5.
        !
        !                        LBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT X = XS AND X = XF.
        !
        !                          = 0  IF THE SOLUTION IS PERIODIC IN X,
        !                               I.E. U(L+I, J, K) = U(I, J, K).
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               X = XS AND X = XF.
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               X = XS AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO X IS
        !                               SPECIFIED AT X = XF.
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO X IS SPECIFIED AT
        !                               X = XS AND X = XF.
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO X IS SPECIFIED AT
        !                               X = XS AND THE SOLUTION IS SPECIFIED
        !                               AT X=XF.
        !
        !                        BDXS
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUES OF THE DERIVATIVE OF THE SOLUTION
        !                          WITH RESPECT TO X AT X = XS.
        !
        !                          WHEN LBDCND = 3 OR 4,
        !
        !                            BDXS(J, K) = (D/DX)U(XS, Y(J), Z(K)),
        !                            J=1, 2, ..., M+1,      K=1, 2, ..., N+1.
        !
        !                          WHEN LBDCND HAS ANY OTHER VALUE, BDXS
        !                          IS A DUMMY VARIABLE. BDXS MUST BE
        !                          DIMENSIONED AT LEAST (M+1)*(N+1).
        !
        !                        BDXF
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE
        !                          VALUES OF THE DERIVATIVE OF THE SOLUTION
        !                          WITH RESPECT TO X AT X = XF.
        !
        !                          WHEN LBDCND = 2 OR 3,
        !
        !                            BDXF(J, K) = (D/DX)U(XF, Y(J), Z(K)),
        !                            J=1, 2, ..., M+1,      K=1, 2, ..., N+1.
        !
        !                          WHEN LBDCND HAS ANY OTHER VALUE, BDXF IS
        !                          A DUMMY VARIABLE.  BDXF MUST BE
        !                          DIMENSIONED AT LEAST (M+1)*(N+1).
        !
        !                        YS, YF
        !                          THE RANGE OF Y, I.E. YS .LE. Y .LE. YF.
        !                          YS MUST BE LESS THAN YF.
        !
        !                        M
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL (YS, YF) IS SUBDIVIDED.
        !                          HENCE, THERE WILL BE M+1 GRID POINTS IN
        !                          THE Y-DIRECTION GIVEN BY Y(J) = YS+(J-1)DY
        !                          FOR J=1, 2, ..., M+1,
        !                          WHERE DY = (YF-YS)/M IS THE PANEL WIDTH.
        !                          M MUST BE AT LEAST 5.
        !
        !                        MBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT Y = YS AND Y = YF.
        !
        !                          = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.
        !                               U(I, M+J, K) = U(I, J, K).
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               Y = YS AND Y = YF.
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               Y = YS AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO Y IS
        !                               SPECIFIED AT Y = YF.
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO Y IS SPECIFIED AT
        !                               Y = YS AND Y = YF.
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO Y IS SPECIFIED AT
        !                               AT Y = YS AND THE SOLUTION IS
        !                               SPECIFIED AT Y=YF.
        !
        !                        BDYS
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
        !                          THE VALUES OF THE DERIVATIVE OF THE
        !                          SOLUTION WITH RESPECT TO Y AT Y = YS.
        !
        !                          WHEN MBDCND = 3 OR 4,
        !
        !                            BDYS(I, K) = (D/DY)U(X(I), YS, Z(K)),
        !                            I=1, 2, ..., L+1,      K=1, 2, ..., N+1.
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDYS
        !                          IS A DUMMY VARIABLE. BDYS MUST BE
        !                          DIMENSIONED AT LEAST (L+1)*(N+1).
        !
        !                        BDYF
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
        !                          THE VALUES OF THE DERIVATIVE OF THE
        !                          SOLUTION WITH RESPECT TO Y AT Y = YF.
        !
        !                          WHEN MBDCND = 2 OR 3,
        !
        !                            BDYF(I, K) = (D/DY)U(X(I), YF, Z(K)),
        !                            I=1, 2, ..., L+1,      K=1, 2, ..., N+1.
        !
        !                          WHEN MBDCND HAS ANY OTHER VALUE, BDYF
        !                          IS A DUMMY VARIABLE. BDYF MUST BE
        !                          DIMENSIONED AT LEAST (L+1)*(N+1).
        !
        !                        ZS, ZF
        !                          THE RANGE OF Z, I.E. ZS .LE. Z .LE. ZF.
        !                          ZS MUST BE LESS THAN ZF.
        !
        !                        N
        !                          THE NUMBER OF PANELS INTO WHICH THE
        !                          INTERVAL (ZS, ZF) IS SUBDIVIDED.
        !                          HENCE, THERE WILL BE N+1 GRID POINTS
        !                          IN THE Z-DIRECTION GIVEN BY
        !                          Z(K) = ZS+(K-1)DZ FOR K=1, 2, ..., N+1,
        !                          WHERE DZ = (ZF-ZS)/N IS THE PANEL WIDTH.
        !                          N MUST BE AT LEAST 5.
        !
        !                        NBDCND
        !                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
        !                          AT Z = ZS AND Z = ZF.
        !
        !                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
        !                               U(I, J, N+K) = U(I, J, K).
        !                          = 1  IF THE SOLUTION IS SPECIFIED AT
        !                               Z = ZS AND Z = ZF.
        !                          = 2  IF THE SOLUTION IS SPECIFIED AT
        !                               Z = ZS AND THE DERIVATIVE OF THE
        !                               SOLUTION WITH RESPECT TO Z IS
        !                               SPECIFIED AT Z = ZF.
        !                          = 3  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO Z IS SPECIFIED AT
        !                               Z = ZS AND Z = ZF.
        !                          = 4  IF THE DERIVATIVE OF THE SOLUTION
        !                               WITH RESPECT TO Z IS SPECIFIED AT
        !                               Z = ZS AND THE SOLUTION IS SPECIFIED
        !                               AT Z=ZF.
        !
        !                        BDZS
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
        !                          THE VALUES OF THE DERIVATIVE OF THE
        !                          SOLUTION WITH RESPECT TO Z AT Z = ZS.
        !
        !                          WHEN NBDCND = 3 OR 4,
        !
        !                            BDZS(I, J) = (D/DZ)U(X(I), Y(J), ZS),
        !                            I=1, 2, ..., L+1,      J=1, 2, ..., M+1.
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDZS
        !                          IS A DUMMY VARIABLE. BDZS MUST BE
        !                          DIMENSIONED AT LEAST (L+1)*(M+1).
        !
        !                        BDZF
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
        !                          THE VALUES OF THE DERIVATIVE OF THE
        !                          SOLUTION WITH RESPECT TO Z AT Z = ZF.
        !
        !                          WHEN NBDCND = 2 OR 3,
        !
        !                            BDZF(I, J) = (D/DZ)U(X(I), Y(J), ZF),
        !                            I=1, 2, ..., L+1,      J=1, 2, ..., M+1.
        !
        !                          WHEN NBDCND HAS ANY OTHER VALUE, BDZF
        !                          IS A DUMMY VARIABLE. BDZF MUST BE
        !                          DIMENSIONED AT LEAST (L+1)*(M+1).
        !
        !                        ELMBDA
        !                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
        !                          EQUATION. IF LAMBDA .GT. 0, A SOLUTION
        !                          MAY NOT EXIST.  HOWEVER, hw3crt WILL
        !                          ATTEMPT TO FIND A SOLUTION.
        !
        !                        LDIMF
        !                          THE ROW (OR FIRST) DIMENSION OF THE
        !                          ARRAYS F, BDYS, BDYF, BDZS, AND BDZF AS IT
        !                          APPEARS IN THE PROGRAM CALLING hw3crt.
        !                          THIS PARAMETER IS USED TO SPECIFY THE
        !                          VARIABLE DIMENSION OF THESE ARRAYS.
        !                          LDIMF MUST BE AT LEAST L+1.
        !
        !                        MDIMF
        !                          THE COLUMN (OR SECOND) DIMENSION OF THE
        !                          ARRAY F AND THE ROW (OR FIRST) DIMENSION
        !                          OF THE ARRAYS BDXS AND BDXF AS IT APPEARS
        !                          IN THE PROGRAM CALLING hw3crt.  THIS
        !                          PARAMETER IS USED TO SPECIFY THE VARIABLE
        !                          DIMENSION OF THESE ARRAYS.
        !                          MDIMF MUST BE AT LEAST M+1.
        !
        !                        F
        !                          A THREE-DIMENSIONAL ARRAY OF DIMENSION AT
        !                          AT LEAST (L+1)*(M+1)*(N+1), SPECIFYING THE
        !                          VALUES OF THE RIGHT SIDE OF THE HELMHOLZ
        !                          EQUATION AND BOUNDARY VALUES (IF ANY).
        !
        !                          ON THE INTERIOR, F IS DEFINED AS FOLLOWS:
        !                          FOR I=2, 3, ..., L,  J=2, 3, ..., M,
        !                          AND K=2, 3, ..., N
        !                          F(I, J, K) = F(X(I), Y(J), Z(K)).
        !
        !                          ON THE BOUNDARIES, F IS DEFINED AS FOLLOWS:
        !                          FOR J=1, 2, ..., M+1,  K=1, 2, ..., N+1,
        !                          AND I=1, 2, ..., L+1
        !
        !                          LBDCND      F(1, J, K)         F(L+1, J, K)
        !                          ------   ---------------   ---------------
        !
        !                            0      F(XS, Y(J), Z(K))   F(XS, Y(J), Z(K))
        !                            1      U(XS, Y(J), Z(K))   U(XF, Y(J), Z(K))
        !                            2      U(XS, Y(J), Z(K))   F(XF, Y(J), Z(K))
        !                            3      F(XS, Y(J), Z(K))   F(XF, Y(J), Z(K))
        !                            4      F(XS, Y(J), Z(K))   U(XF, Y(J), Z(K))
        !
        !                          MBDCND      F(I, 1, K)         F(I, M+1, K)
        !                          ------   ---------------   ---------------
        !
        !                            0      F(X(I), YS, Z(K))   F(X(I), YS, Z(K))
        !                            1      U(X(I), YS, Z(K))   U(X(I), YF, Z(K))
        !                            2      U(X(I), YS, Z(K))   F(X(I), YF, Z(K))
        !                            3      F(X(I), YS, Z(K))   F(X(I), YF, Z(K))
        !                            4      F(X(I), YS, Z(K))   U(X(I), YF, Z(K))
        !
        !                          NBDCND      F(I, J, 1)         F(I, J, N+1)
        !                          ------   ---------------   ---------------
        !
        !                            0      F(X(I), Y(J), ZS)   F(X(I), Y(J), ZS)
        !                            1      U(X(I), Y(J), ZS)   U(X(I), Y(J), ZF)
        !                            2      U(X(I), Y(J), ZS)   F(X(I), Y(J), ZF)
        !                            3      F(X(I), Y(J), ZS)   F(X(I), Y(J), ZF)
        !                            4      F(X(I), Y(J), ZS)   U(X(I), Y(J), ZF)
        !
        !                          NOTE:
        !                          IF THE TABLE CALLS FOR BOTH THE SOLUTION
        !                          U AND THE RIGHT SIDE F ON A BOUNDARY,
        !                          THEN THE SOLUTION MUST BE SPECIFIED.
        !
        !
        ! ON OUTPUT              F
        !                          CONTAINS THE SOLUTION U(I, J, K) OF THE
        !                          FINITE DIFFERENCE APPROXIMATION FOR THE
        !                          GRID POINT (X(I), Y(J), Z(K)) FOR
        !                          I=1, 2, ..., L+1, J=1, 2, ..., M+1,
        !                          AND K=1, 2, ..., N+1.
        !
        !                        PERTRB
        !                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
        !                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
        !                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
        !                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
        !                          CALCULATED AND SUBTRACTED FROM F, WHICH
        !                          ENSURES THAT A SOLUTION EXISTS.  PWSCRT
        !                          THEN COMPUTES THIS SOLUTION, WHICH IS A
        !                          LEAST SQUARES SOLUTION TO THE ORIGINAL
        !                          APPROXIMATION.  THIS SOLUTION IS NOT
        !                          UNIQUE AND IS UNNORMALIZED.  THE VALUE OF
        !                          PERTRB SHOULD BE SMALL COMPARED TO THE
        !                          THE RIGHT SIDE F.  OTHERWISE, A SOLUTION
        !                          IS OBTAINED TO AN ESSENTIALLY DIFFERENT
        !                          PROBLEM.  THIS COMPARISON SHOULD ALWAYS
        !                          BE MADE TO INSURE THAT A MEANINGFUL
        !                          SOLUTION HAS BEEN OBTAINED.
        !
        !                        ierror
        !                          AN ERROR FLAG THAT INDICATES INVALID INPUT
        !                          PARAMETERS.  EXCEPT FOR NUMBERS 0 AND 12,
        !                          A SOLUTION IS NOT ATTEMPTED.
        !
        !                          =  0  NO ERROR
        !                          =  1  XS .GE. XF
        !                          =  2  L .LT. 5
        !                          =  3  LBDCND .LT. 0 .OR. LBDCND .GT. 4
        !                          =  4  YS .GE. YF
        !                          =  5  M .LT. 5
        !                          =  6  MBDCND .LT. 0 .OR. MBDCND .GT. 4
        !                          =  7  ZS .GE. ZF
        !                          =  8  N .LT. 5
        !                          =  9  NBDCND .LT. 0 .OR. NBDCND .GT. 4
        !                          = 10  LDIMF .LT. L+1
        !                          = 11  MDIMF .LT. M+1
        !                          = 12  LAMBDA .GT. 0
        !                          = 20 If the dynamic allocation of real and
        !                               complex work space required for solution
        !                               fails (for example if N, M are too large
        !                               for your computer)
        !
        !                          SINCE THIS IS THE ONLY MEANS OF INDICATING
        !                          A POSSIBLY INCORRECT CALL TO hw3crt, THE
        !                          USER SHOULD TEST ierror AFTER THE CALL.
        !
        ! SPECIAL CONDITIONS     NONE
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED Files         fish.f, pois3d.f, fftpack.f, comf.f
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
        ! ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE DIFFERENCE
        !                        EQUATIONS, INCORPORATES BOUNDARY DATA, AND
        !                        ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS AND
        !                        THEN CALLS POIS3D TO SOLVE THE SYSTEM.
        !
        ! TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
        !                        IS ROUGHLY PROPORTIONAL TO
        !                          L*M*N*(LOG2(L)+LOG2(M)+5),
        !                        BUT ALSO DEPENDS ON INPUT PARAMETERS LBDCND
        !                        AND MBDCND.
        !
        ! ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN
        !                        A LOSS OF NO MORE THAN FOUR SIGNIFICANT
        !                        DIGITS FOR L, M AND N AS LARGE AS 32.
        !                        MORE DETAILED INFORMATION ABOUT ACCURACY
        !                        CAN BE FOUND IN THE DOCUMENTATION FOR
        !                        ROUTINE POIS3D WHICH IS THE ROUTINE THAT
        !                        ACTUALLY SOLVES THE FINITE DIFFERENCE
        !                        EQUATIONS.
        !
        ! REFERENCES             NONE
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)      :: l
        integer (ip), intent (in)      :: lbdcnd
        integer (ip), intent (in)      :: m
        integer (ip), intent (in)      :: mbdcnd
        integer (ip), intent (in)      :: n
        integer (ip), intent (in)      :: nbdcnd
        integer (ip), intent (in)      :: ldimf
        integer (ip), intent (in)      :: mdimf
        integer (ip), intent (out)     :: ierror
        real (wp),    intent (in)      :: xs
        real (wp),    intent (in)      :: xf
        real (wp),    intent (in)      :: ys
        real (wp),    intent (in)      :: yf
        real (wp),    intent (in)      :: zs
        real (wp),    intent (in)      :: zf
        real (wp),    intent (in)      :: elmbda
        real (wp),    intent (out)     :: pertrb
        real (wp), contiguous, intent (in)      :: bdxs(:,:)
        real (wp), contiguous, intent (in)      :: bdxf(:,:)
        real (wp), contiguous, intent (in)      :: bdys(:,:)
        real (wp), contiguous, intent (in)      :: bdyf(:,:)
        real (wp), contiguous, intent (in)      :: bdzs(:,:)
        real (wp), contiguous, intent (in)      :: bdzf(:,:)
        real (wp), contiguous, intent (in out)  :: f(:,:,:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        type (FishpackWorkspace) :: workspace
        !-----------------------------------------------

        ! initialize error flag
        ierror = 0

        ! Check if input values are valid
        if (xf <= xs) ierror = 1
        if (l < 5) ierror = 2
        if (lbdcnd<0 .or. lbdcnd>4) ierror = 3
        if (yf <= ys) ierror = 4
        if (m < 5) ierror = 5
        if (mbdcnd<0 .or. mbdcnd>4) ierror = 6
        if (zf <= zs) ierror = 7
        if (n < 5) ierror = 8
        if (nbdcnd<0 .or. nbdcnd>4) ierror = 9
        if (ldimf < l + 1) ierror = 10
        if (mdimf < m + 1) ierror = 11
        if (ierror /= 0) return

        ! Allocate real workspace array
        associate( &
            irwk => 30+l+m+5*n+max(l, m, n)+7*(int((l+1)/2)+int((m+1)/2)), & ! Estimate required workspace length (generous estimate)
            icwk => 0 &
            )
            call workspace%create( irwk, icwk, ierror )
        end associate

        ! Check if allocation was succcessful
        if (ierror == 20) return

        ! solve system
        associate( rew => workspace%rew )
            call hw3crtt(xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, mbdcnd, bdys, &
                bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, ldimf, &
                mdimf, f, pertrb, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hw3crt

    subroutine hw3crtT(xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, &
        mbdcnd, bdys, bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, &
        ldimf, mdimf, f, pertrb, ierror, w)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)      :: l
        integer (ip), intent (in)      :: lbdcnd
        integer (ip), intent (in)      :: m
        integer (ip), intent (in)      :: mbdcnd
        integer (ip), intent (in)      :: n
        integer (ip), intent (in)      :: nbdcnd
        integer (ip), intent (in)      :: ldimf
        integer (ip), intent (in)      :: mdimf
        integer (ip), intent (out)     :: ierror
        real (wp),    intent (in)      :: xs
        real (wp),    intent (in)      :: xf
        real (wp),    intent (in)      :: ys
        real (wp),    intent (in)      :: yf
        real (wp),    intent (in)      :: zs
        real (wp),    intent (in)      :: zf
        real (wp),    intent (in)      :: elmbda
        real (wp),    intent (out)     :: pertrb
        real (wp),    intent (in)      :: bdxs(mdimf, *)
        real (wp),    intent (in)      :: bdxf(mdimf, *)
        real (wp),    intent (in)      :: bdys(ldimf, *)
        real (wp),    intent (in)      :: bdyf(ldimf, *)
        real (wp),    intent (in)      :: bdzs(ldimf, *)
        real (wp),    intent (in)      :: bdzf(ldimf, *)
        real (wp),    intent (in out)  :: f(ldimf, mdimf, *)
        real (wp),    intent (in out)  :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: mstart, mstop, mp1, mp, munk, np, np1, nstart, nstop, &
            nunk, lp1, lp, lstart, lstop, j, k, lunk, i, iwb, iwc, iww, &
            mstpm1, lstpm1, nstpm1, nperod, ir
        real::dy, twbydy, c2, dz, twbydz, c3, dx, c1, twbydx, xlp, ylp, zlp, s1, s2, s
        !-----------------------------------------------

        dy = (yf - ys)/real(m)
        twbydy = 2./dy
        c2 = 1./dy**2
        mstart = 1
        mstop = m
        mp1 = m + 1
        mp = mbdcnd + 1
        go to (104, 101, 101, 102, 102) mp
101 continue
    mstart = 2
102 continue
    go to (104, 104, 103, 103, 104) mp
103 continue
    mstop = mp1
104 continue
    munk = mstop - mstart + 1
    dz = (zf - zs)/real(n)
    twbydz = 2./dz
    np = nbdcnd + 1
    c3 = 1./dz**2
    np1 = n + 1
    nstart = 1
    nstop = n
    go to (108, 105, 105, 106, 106) np
105 continue
    nstart = 2
106 continue
    go to (108, 108, 107, 107, 108) np
107 continue
    nstop = np1
108 continue
    nunk = nstop - nstart + 1
    lp1 = l + 1
    dx = (xf - xs)/real(l)
    c1 = 1./dx**2
    twbydx = 2./dx
    lp = lbdcnd + 1
    lstart = 1
    lstop = l
    !
    !     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
    !
    go to (122, 109, 109, 112, 112) lp
109 continue
    lstart = 2
    f(2, mstart:mstop, nstart:nstop) = F(2, mstart:mstop, nstart:nstop) - &
        c1*F(1, mstart:mstop, nstart:nstop)
    go to 115
112 continue
    f(1, mstart:mstop, nstart:nstop) = F(1, mstart:mstop, nstart:nstop) + &
        twbydx*BDXS(mstart:mstop, nstart:nstop)
115 continue
    go to (122, 116, 119, 119, 116) lp
116 continue
    f(l, mstart:mstop, nstart:nstop) = F(l, mstart:mstop, nstart:nstop) - &
        c1*F(lp1, mstart:mstop, nstart:nstop)
    go to 122
119 continue
    lstop = lp1
    f(lp1, mstart:mstop, nstart:nstop) = F(lp1, mstart:mstop, nstart:nstop &
        ) - twbydx*BDXF(mstart:mstop, nstart:nstop)
122 continue
    lunk = lstop - lstart + 1
    !
    !     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
    !
    go to (136, 123, 123, 126, 126) mp
123 continue
    f(lstart:lstop, 2, nstart:nstop) = F(lstart:lstop, 2, nstart:nstop) - &
        c2*F(lstart:lstop, 1, nstart:nstop)
    go to 129
126 continue
    f(lstart:lstop, 1, nstart:nstop) = F(lstart:lstop, 1, nstart:nstop) + &
        twbydy*BDYS(lstart:lstop, nstart:nstop)
129 continue
    go to (136, 130, 133, 133, 130) mp
130 continue
    f(lstart:lstop, m, nstart:nstop) = F(lstart:lstop, m, nstart:nstop) - &
        c2*F(lstart:lstop, mp1, nstart:nstop)
    go to 136
133 continue
    f(lstart:lstop, mp1, nstart:nstop) = F(lstart:lstop, mp1, nstart:nstop &
        ) - twbydy*BDYF(lstart:lstop, nstart:nstop)
136 continue
    go to (150, 137, 137, 140, 140) np
137 continue
    f(lstart:lstop, mstart:mstop, 2) = F(lstart:lstop, mstart:mstop, 2) - &
        c3*F(lstart:lstop, mstart:mstop, 1)
    go to 143
140 continue
    f(lstart:lstop, mstart:mstop, 1) = F(lstart:lstop, mstart:mstop, 1) + &
        twbydz*BDZS(lstart:lstop, mstart:mstop)
143 continue
    go to (150, 144, 147, 147, 144) np
144 continue
    f(lstart:lstop, mstart:mstop, n) = F(lstart:lstop, mstart:mstop, n) - &
        c3*F(lstart:lstop, mstart:mstop, np1)
    go to 150
147 continue
    f(lstart:lstop, mstart:mstop, np1) = F(lstart:lstop, mstart:mstop, np1 &
        ) - twbydz*BDZF(lstart:lstop, mstart:mstop)
!
!     DEFINE A, B, C COEFFICIENTS IN W-ARRAY.
!
150 continue
    iwb = nunk + 1
    iwc = iwb + nunk
    iww = iwc + nunk
    w(:nunk) = c3
    w(iwc:nunk-1+iwc) = c3
    w(iwb:nunk-1+iwb) = (-2.*c3) + elmbda
    go to (155, 155, 153, 152, 152) np
152 continue
    w(iwc) = 2.*c3
153 continue
    go to (155, 155, 154, 154, 155) np
154 continue
    w(iwb-1) = 2.*c3
155 continue
    pertrb = 0.
    !
    !     FOR SINGULAR PROBLEMS ADJUST DATA TO INSURE A SOLUTION WILL EXIST.
    !
    go to (156, 172, 172, 156, 172) lp
156 continue
    go to (157, 172, 172, 157, 172) mp
157 continue
    go to (158, 172, 172, 158, 172) np
158 continue
    if (elmbda >= 0.) then
        if (elmbda /= 0.) then
            ierror = 12
        else
            mstpm1 = mstop - 1
            lstpm1 = lstop - 1
            nstpm1 = nstop - 1
            xlp = (2 + lp)/3
            ylp = (2 + mp)/3
            zlp = (2 + np)/3
            s1 = 0.
            do k = 2, nstpm1
                do j = 2, mstpm1
                    s1 = s1 + SUM(F(2:lstpm1, j, k))
                    s1 = s1 + (F(1, j, k)+F(lstop, j, k))/xlp
                end do
                s2 = SUM(F(2:lstpm1, 1, k)+F(2:lstpm1, mstop, k))
                s2 = (s2 + (F(1, 1, k)+F(1, mstop, k)+F(lstop, 1, k)+F(lstop, &
                    mstop, k))/xlp)/ylp
                s1 = s1 + s2
            end do
            s = (F(1, 1, 1)+F(lstop, 1, 1)+F(1, 1, nstop)+F(lstop, 1, nstop)+F(1 &
                , mstop, 1)+F(lstop, mstop, 1)+F(1, mstop, nstop)+F(lstop, mstop &
                , nstop))/(xlp*ylp)
            do j = 2, mstpm1
                s = s + SUM(F(2:lstpm1, j, 1)+F(2:lstpm1, j, nstop))
            end do
            s2 = 0.
            s2 = SUM(F(2:lstpm1, 1, 1)+F(2:lstpm1, 1, nstop)+F(2:lstpm1, &
                mstop, 1)+F(2:lstpm1, mstop, nstop))
            s = s2/ylp + s
            s2 = 0.
            s2 = SUM(F(1, 2:mstpm1, 1)+F(1, 2:mstpm1, nstop)+F(lstop, 2: &
                mstpm1, 1)+F(lstop, 2:mstpm1, nstop))
            s = s2/xlp + s
            pertrb = (s/zlp + s1)/((real(lunk + 1) - xlp)*(real(munk &
                + 1) - ylp)*(real(nunk + 1) - zlp))
            f(:lunk, :munk, :nunk) = F(:lunk, :munk, :nunk) - pertrb
        end if
    end if
172 continue
    nperod = 0
    if (nbdcnd /= 0) then
        nperod = 1
        w(1) = 0.
        w(iww-1) = 0.
    end if

    call POIS3DD(lbdcnd, lunk, c1, mbdcnd, munk, c2, nperod, nunk, w, &
        W(iwb), W(iwc), ldimf, mdimf, F(lstart, mstart, nstart), ir, W(iww))
    !
    !     FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS.
    !
    if (lp == 1) then
        if (mp == 1) then
            f(1, mp1, nstart:nstop) = F(1, 1, nstart:nstop)
            mstop = mp1
        end if
        if (np == 1) then
            f(1, mstart:mstop, np1) = F(1, mstart:mstop, 1)
            nstop = np1
        end if
        f(lp1, mstart:mstop, nstart:nstop) = F(1, mstart:mstop, nstart: nstop)
    end if
    if (mp == 1) then
        if (np == 1) then
            f(lstart:lstop, 1, np1) = F(lstart:lstop, 1, 1)
            nstop = np1
        end if
        f(lstart:lstop, mp1, nstart:nstop) = F(lstart:lstop, 1, nstart:nstop)
    end if
    if (np == 1) then
        f(lstart:lstop, mstart:mstop, np1) = F(lstart:lstop, mstart:mstop,1)
    end if

end subroutine hw3crtT

end module module_hw3crt
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
