module module_blktri

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_comf, only: &
        epmach, &
        psgf, &
        ppspf, &
        ppsgf

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: blktri
    public :: blktrii
    public :: blktri_unit_test

contains
    !
    !*****************************************************************************************
    !
    subroutine blktri_unit_test()
        !     file tblktri.f
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
        !     program to illustrate the use of subroutine blktri to
        !     solve the equation
        !
        !     .5/s*(d/ds)(.5/s*du/ds)+.5/t*(d/dt)(.5/t*du/dt)
        !                                                          (1)
        !                   = 15/4*s*t*(s**4+t**4)
        !
        !     on the rectangle 0 .lt. s .lt. 1 and 0 .lt. t .lt. 1
        !     with the boundary conditions
        !
        !     u(0, t) = 0
        !                            0 .le. t .le. 1
        !     u(1, t) = t**5
        !
        !     and
        !
        !     u(s, 0) = 0
        !                            0 .le. s .le. 1
        !     u(s, 1) = s**5
        !
        !     the exact solution of this problem is u(s, t) = (s*t)**5
        !
        !     define the integers m = 50 and n = 63. then define the
        !     grid increments deltas = 1/(m+1) and deltat = 1/(n+1).
        !
        !     the grid is then given by s(i) = i*deltas for i = 1, ..., m
        !     and t(j) = j*deltat for j = 1, ..., n.
        !
        !     the approximate solution is given as the solution to
        !     the following finite difference approximation of equation (1).
        !
        !     .5/(s(i)*deltas)*((u(i+1, j)-u(i, j))/(2*s(i+.5)*deltas)
        !                     -(u(i, j)-u(i-1, j))/(2*s(i-.5)*deltas))
        !     +.5/(t(i)*deltat)*((u(i, j+1)-u(i, j))/(2*t(i+.5)*deltat) (2)
        !                     -(u(i, j)-u(i, j-1))/(2*t(i-.5)*deltat))
        !               = 15/4*s(i)*t(j)*(s(i)**4+t(j)**4)
        !
        !             where s(i+.5) = .5*(s(i+1)+s(i))
        !                   s(i-.5) = .5*(s(i)+s(i-1))
        !                   t(i+.5) = .5*(t(i+1)+t(i))
        !                   t(i-.5) = .5*(t(i)+t(i-1))
        !
        !     the approach is to write equation (2) in the form
        !
        !     am(i)*u(i-1, j)+bm(i, j)*u(i, j)+cm(i)*u(i+1, j)
        !       +an(j)*u(i, j-1)+bn(j)*u(i, j)+cn(j)*u(i, j+1)      (3)
        !           = y(i, j)
        !
        !     and then call subroutine blktri to determine u(i, j)
        !
        !-----------------------------------------------
        ! Dictionary: Local variables
        !-----------------------------------------------
        integer (ip) :: iflg, np, n, mp, m, idimy, i, j, ierror
        real (wp), dimension(75, 105) :: y
        real (wp), dimension(75) :: am, bm, cm
        real (wp), dimension(105) :: an, bn, cn
        real (wp), dimension(75) :: s
        real (wp), dimension(105) :: t
        real (wp) :: ds, dt, discretization_error
        type (FishpackWorkspace)  :: workspace
        !-----------------------------------------------

        iflg = 0
        np = 1
        n = 63
        mp = 1
        m = 50
        idimy = 75
        !
        !     Generate and store grid points for the purpose of computing the
        !     coefficients and the array y.
        !
        ds = 1.0_wp/(m + 1)

        do i = 1, m
            s(i) = real(i, kind = wp) * ds
        end do

        dt = 1.0_wp/(n + 1)

        do j = 1, n
            t(j) = real(j, kind = wp) * dt
        end do
        !
        !     Compute the coefficients am, bm, cm corresponding to the s direction
        !
        associate( &
            half_ds => ds/2, &
            two_ds => 2.0_wp * ds &
            )
            do i = 1, m
                associate( &
                    temp1 => 1.0_wp/(s(i)*two_ds), &
                    temp2 => 1.0_wp/((s(i)-half_ds)*two_ds), &
                    temp3 => 1.0_wp/((s(i)+half_ds)*two_ds) &
                    )
                    am(i) = temp1*temp2
                    cm(i) = temp1*temp3
                    bm(i) = -(am(i)+cm(i))
                end associate
            end do
        end associate

        !
        !     COMPUTE THE COEFFICIENTS AN, BN, CN CORRESPONDING TO THE T DIRECTION
        !
        associate( &
            half_dt => dt/2, &
            two_dt => 2.0_wp * dt &
            )

            do j = 1, n
                associate( &
                    temp1 => 1.0_wp / ( t(j) * two_dt ), &
                    temp2 => 1.0_wp / ( (t(j)-half_dt) * two_dt ), &
                    temp3 => 1.0_wp / ( (t(j)+half_dt) * two_dt ) &
                    )
                    an(j) = temp1*temp2
                    cn(j) = temp1*temp3
                    bn(j) = -(an(j)+cn(j))
                end associate
            end do
        end associate

        !
        !     Compute right side of equation
        !
        do j = 1, n
            y(:m, j) = 3.75_wp * s(:m) * t(j) * ( s(:m)**4 + t(j)**4 )
        end do
        !
        !     The nonzero boundary conditions enter the linear system via
        !     the right side y(i, j). if the equations (3) given above
        !     are evaluated at i=m and j=1, ..., n then the term cm(m)*u(m+1, j)
        !     is known from the boundary condition to be cm(m)*t(j)**5.
        !     therefore this term can be included in the right side y(m, j).
        !     the same analysis applies at j=n and i=1, .., m. note that the
        !     corner at j=n, i=m includes contributions from both boundaries.
        !
        y(m, :n) = y(m, :n) - cm(m)*t(:n)**5
        y(:m, n) = y(:m, n) - cn(n)*s(:m)**5
        !
        !     Determine the approximate solution u(i, j)
        !
        call blktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, idimy, y, ierror, workspace)

        iflg = iflg + 1

        do while(iflg - 1 <= 0)
            call blktri (iflg, np, n, an, bn, cn, mp, m, am, bm, cm, idimy &
                , y, ierror, workspace )
            iflg = iflg + 1
        end do

        discretization_error = 0.0_wp

        do j = 1, n
            do i = 1, m
                associate( local_error => abs(y(i, j)-(s(i)*t(j))**5) )
                    discretization_error = max(local_error, discretization_error)
                end associate
            end do
        end do

        !     Print earlier output from platform with 64 bit floating point
        !     arithemtic followed by the output from this computer
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     BLKTRI UNIT TEST ******* '
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') '     IERROR = 0,  Discretization Error = 1.6478E-05'
        write( stdout, '(A)') ''
        write( stdout, '(A)') '    The output from your computer is: '
        write( stdout, *) '    IERROR =', IERROR, ' Discretization Error = ', &
            discretization_error

        !! release dynamically allocated work space
        call workspace%destroy()

    end subroutine blktri_unit_test
    !
    !*****************************************************************************************
    !
    subroutine blktri( iflg, np, n, an, bn, cn, mp, m, am, bm, cm, &
        idimy, y, ierror, w )
        !
        !     file blktri.f
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
        !     SUBROUTINE BLKTRI (IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, IDIMY, Y,
        !    +                   IERROR, W)
        !
        !
        !
        ! DIMENSION OF           AN(N), BN(N), CN(N), AM(M), BM(M), CM(M), Y(IDIMY, N),
        ! ARGUMENTS
        !
        ! LATEST REVISION        JUNE 2004
        !
        ! USAGE                  CALL BLKTRI (IFLG, NP, N, AN, BN, CN, MP, M, AM, BM,
        !                                     CM, IDIMY, Y, IERROR, W)
        !
        ! PURPOSE                BLKTRI SOLVES A SYSTEM OF LINEAR EQUATIONS
        !                        OF THE FORM
        !
        !                        AN(J)*X(I, J-1) + AM(I)*X(I-1, J) +
        !                        (BN(J)+BM(I))*X(I, J) + CN(J)*X(I, J+1) +
        !                        CM(I)*X(I+1, J) = Y(I, J)
        !
        !                        FOR I = 1, 2, ..., M  AND  J = 1, 2, ..., N.
        !
        !                        I+1 AND I-1 ARE EVALUATED MODULO M AND
        !                        J+1 AND J-1 MODULO N, I.E.,
        !
        !                        X(I, 0) = X(I, N),  X(I, N+1) = X(I, 1),
        !                        X(0, J) = X(M, J),  X(M+1, J) = X(1, J).
        !
        !                        THESE EQUATIONS USUALLY RESULT FROM THE
        !                        DISCRETIZATION OF SEPARABLE ELLIPTIC
        !                        EQUATIONS.  BOUNDARY CONDITIONS MAY BE
        !                        DIRICHLET, NEUMANN, OR PERIODIC.
        !
        ! ARGUMENTS
        !
        ! ON INPUT               IFLG
        !
        !                          = 0  INITIALIZATION ONLY.
        !                               CERTAIN QUANTITIES THAT DEPEND ON NP,
        !                               N, AN, BN, AND CN ARE COMPUTED AND
        !                               STORED IN DERIVED data type w (see
        !                               description of w below)
        !
        !                          = 1  THE QUANTITIES THAT WERE COMPUTED
        !                               IN THE INITIALIZATION ARE USED
        !                               TO OBTAIN THE SOLUTION X(I, J).
        !
        !                               NOTE:
        !                               A CALL WITH IFLG=0 TAKES
        !                               APPROXIMATELY ONE HALF THE TIME
        !                               AS A CALL WITH IFLG = 1.
        !                               HOWEVER, THE INITIALIZATION DOES
        !                               NOT HAVE TO BE REPEATED UNLESS NP,
        !                               N, AN, BN, OR CN CHANGE.
        !
        !                        NP
        !                          = 0  IF AN(1) AND CN(N) ARE NOT ZERO,
        !                               WHICH CORRESPONDS TO PERIODIC
        !                               BOUNARY CONDITIONS.
        !
        !                          = 1  IF AN(1) AND CN(N) ARE ZERO.
        !
        !                        N
        !                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
        !                          N MUST BE GREATER THAN 4.
        !                          THE OPERATION COUNT IS PROPORTIONAL TO
        !                          MNLOG2(N), HENCE N SHOULD BE SELECTED
        !                          LESS THAN OR EQUAL TO M.
        !
        !                        AN, BN, CN
        !                          ONE-DIMENSIONAL ARRAYS OF LENGTH N
        !                          THAT SPECIFY THE COEFFICIENTS IN THE
        !                          LINEAR EQUATIONS GIVEN ABOVE.
        !
        !                        MP
        !                          = 0  IF AM(1) AND CM(M) ARE NOT ZERO,
        !                               WHICH CORRESPONDS TO PERIODIC
        !                               BOUNDARY CONDITIONS.
        !
        !                          = 1  IF AM(1) = CM(M) = 0  .
        !
        !                        M
        !                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
        !                           M MUST BE GREATER THAN 4.
        !
        !                        AM, BM, CM
        !                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
        !                          SPECIFY THE COEFFICIENTS IN THE LINEAR
        !                          EQUATIONS GIVEN ABOVE.
        !
        !                        IDIMY
        !                          THE ROW (OR FIRST) DIMENSION OF THE
        !                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS
        !                          IN THE PROGRAM CALLING BLKTRI.
        !                          THIS PARAMETER IS USED TO SPECIFY THE
        !                          VARIABLE DIMENSION OF Y.
        !                          IDIMY MUST BE AT LEAST M.
        !
        !                        Y
        !                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
        !                          THE VALUES OF THE RIGHT SIDE OF THE LINEAR
        !                          SYSTEM OF EQUATIONS GIVEN ABOVE.
        !                          Y MUST BE DIMENSIONED AT LEAST M*N.
        !
        !                        W
        !                          A fortran 90 derived TYPE (FishpackWorkspace) variable
        !                          that must be declared by the user.  The first
        !                          two declarative statements in the user program
        !                          calling BLKTRI must be:
        !
        !                               use type_FishpackWorkspace
        !                               TYPE (FishpackWorkspace) :: W
        !
        !                          The first statement makes the fishpack module
        !                          defined in the file "fish.f" available to the
        !                          user program calling BLKTRI.  The second statement
        !                          declares a derived type variable (defined in
        !                          the module "fish.f") which is used internally
        !                          in BLKTRI to dynamically allocate real and complex
        !                          work space used in solution.  An error flag
        !                          (IERROR = 20) is set if the required work space
        !                          allocation fails (for example if N, M are too large)
        !                          Real and complex values are set in the components
        !                          of W on a initial (IFLG=0) call to BLKTRI.  These
        !                          must be preserved on non-initial calls (IFLG=1)
        !                          to BLKTRI.  This eliminates redundant calculations
        !                          and saves compute time.
        !               ****       IMPORTANT!  The user program calling BLKTRI should
        !                          include the statement:
        !
        !                               CALL FISHFIN(W)
        !
        !                          after the final approximation is generated by
        !                          BLKTRI.  The will deallocate the real and complex
        !                          work space of W.  Failure to include this statement
        !                          could result in serious memory leakage.
        !
        !
        ! ARGUMENTS
        !
        ! ON OUTPUT              Y
        !                          CONTAINS THE SOLUTION X.
        !
        !                        IERROR
        !                          AN ERROR FLAG THAT INDICATES INVALID
        !                          INPUT PARAMETERS.  EXCEPT FOR NUMBER ZER0,
        !                          A SOLUTION IS NOT ATTEMPTED.
        !
        !                        = 0  NO ERROR.
        !                        = 1  M IS LESS THAN 5
        !                        = 2  N IS LESS THAN 5
        !                        = 3  IDIMY IS LESS THAN M.
        !                        = 4  BLKTRI FAILED WHILE COMPUTING RESULTS
        !                             THAT DEPEND ON THE COEFFICIENT ARRAYS
        !                             AN, BN, CN.  CHECK THESE ARRAYS.
        !                        = 5  AN(J)*CN(J-1) IS LESS THAN 0 FOR SOME J.
        !
        !                             POSSIBLE REASONS FOR THIS CONDITION ARE
        !                             1. THE ARRAYS AN AND CN ARE NOT CORRECT
        !                             2. TOO LARGE A GRID SPACING WAS USED
        !                                IN THE DISCRETIZATION OF THE ELLIPTIC
        !                                EQUATION.
        !                             3. THE LINEAR EQUATIONS RESULTED FROM A
        !                                PARTIAL DIFFERENTIAL EQUATION WHICH
        !                                WAS NOT ELLIPTIC.
        !
        !                        = 20 If the dynamic allocation of real and
        !                             complex work space in the derived type
        !                             (FishpackWorkspace) variable W fails (e.g.,
        !                             if N, M are too large for the platform used)
        !
        !
        !                        W
        !                             The derived type (FishpackWorkspace) variable W
        !                             contains real and complex values that must not
        !                             be destroyed if BLKTRI is called again with
        !                             IFLG=1.
        !
        !
        ! SPECIAL CONDITIONS     THE ALGORITHM MAY FAIL IF abs(BM(I)+BN(J))
        !                        IS LESS THAN abs(AM(I))+abs(AN(J))+
        !                        abs(CM(I))+abs(CN(J))
        !                        FOR SOME I AND J. THE ALGORITHM WILL ALSO
        !                        FAIL IF AN(J)*CN(J-1) IS LESS THAN ZERO FOR
        !                        SOME J.
        !                        SEE THE DESCRIPTION OF THE OUTPUT PARAMETER
        !                        IERROR.
        !
        ! I/O                    NONE
        !
        ! PRECISION              SINGLE
        !
        ! REQUIRED FILES         fish.f, comf.f
        !
        ! LANGUAGE               FORTRAN 90
        !
        ! HISTORY                WRITTEN BY PAUL SWARZTRAUBER AT NCAR IN THE
        !                        EARLY 1970'S.  REWRITTEN AND RELEASED IN
        !                        LIBRARIES IN JANUARY 1980. Revised in June
        !                        2004 using Fortan 90 dynamically allocated work
        !                        space and derived data types to eliminate mixed
        !                        mode conflicts in the earlier versions.
        !
        ! ALGORITHM              GENERALIZED CYCLIC REDUCTION
        !
        ! PORTABILITY            FORTRAN 90.  APPROXIMATE MACHINE ACCURACY
        !                        IS COMPUTED IN FUNCTION EPMACH.
        !
        ! REFERENCES             SWARZTRAUBER, P. AND R. SWEET, 'EFFICIENT
        !                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF
        !                        ELLIPTIC EQUATIONS'
        !                        NCAR TN/IA-109, JULY, 1975, 138 PP.
        !
        !                        SWARZTRAUBER P. N., A DIRECT METHOD FOR
        !                        THE DISCRETE SOLUTION OF SEPARABLE
        !                        ELLIPTIC EQUATIONS, S.I.A.M.
        !                        J. NUMER. ANAL., 11(1974) PP. 1136-1150.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip) :: iflg
        integer (ip) :: np
        integer (ip) :: n
        integer (ip) :: mp
        integer (ip) :: m
        integer (ip) :: idimy
        integer (ip) :: ierror
        real (wp) :: an(*)
        real (wp) :: bn(*)
        real (wp) :: cn(*)
        real (wp) :: am(*)
        real (wp) :: bm(*)
        real (wp) :: cm(*)
        real (wp) :: y(idimy, *)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)             :: irwk, icwk
        type (FishpackWorkspace) :: w
        !--------------------------------------------------------------------------------

        if (m < 5) then
            ierror = 1
            return
        end if
        if (n < 3) then
            ierror = 2
            return
        end if
        if (idimy < m) then
            ierror = 3
            return
        end if
        if (iflg == 0) then

            !! compute and allocate real and complex required work space
            call w%get_block_tridiagonal_workpace_dimensions (n, m, irwk, icwk)

            !! Create workspace
            call w%create( irwk, icwk, ierror )
            if (ierror == 20) then
                return
            end if
        end if

        call blktrii( iflg, np, n, an, bn, cn, mp, m, am, bm, cm, idimy, y, &
            ierror, w%rew, w%cxw )

    end subroutine blktri
    !
    !*****************************************************************************************
    !
    subroutine blktrii(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, &
        idimy, y, ierror, w, wc)
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer (ip), intent (in) :: iflg
        integer (ip), intent (in) :: np
        integer (ip), intent (in) :: n
        integer (ip), intent (in) :: mp
        integer (ip) :: m
        integer (ip) :: idimy
        integer (ip) :: ierror
        real (wp) :: an(*)
        real (wp) :: bn(*)
        real (wp) :: cn(*)
        real (wp) :: am(*)
        real (wp) :: bm(*)
        real (wp) :: cm(*)
        real (wp) :: y(idimy, *)
        real (wp) :: w(*)
        complex (wp) :: wc(*)
        !-----------------------------------------------
        !   c o m m o n   b l o c k s
        !-----------------------------------------------
        !...  /cblkt/
        common /cblkt/ npp, k, eps, cnv, nm, ncmplx, ik
        integer   npp, k, nm, ncmplx, ik
        real   eps, cnv
        !-----------------------------------------------
        !   l o c a l   v a r i a b l e s
        !-----------------------------------------------
        integer (ip) :: iw1, iw2, iw3, iww, iwu, iwd, nl, nh, iwah, iwbh

        save iw1, iw2, iw3, iww, iwu, iwd, nl
        !-----------------------------------------------

        !
        ! test m and n for the proper form
        !
        nm = n
        !     check again for solvers which call blktrii directly
        if (m < 5) then
            ierror = 1
            return
        end if
        if (nm < 3) then
            ierror = 2
            return
        end if
        if (idimy < m) then
            ierror = 3
            return
        end if

        if (iflg == 0) then
            nh = n
            npp = np
            if (npp /= 0) then
                nh = nh + 1
            end if
            ik = 2
            k = 1
            ik = ik + ik
            k = k + 1
            do while(nh - ik > 0)
                ik = ik + ik
                k = k + 1
            end do
            nl = ik
            ik = ik + ik
            nl = nl - 1
            iwah = (k - 2)*ik + k + 5
            if (npp == 0) then
                iwbh = iwah + nm + nm
                iw1 = iwbh
                nm = n - 1
            else
                iw1 = iwah
                iwbh = iw1 + nm
            end if
            !     set pointers in real, complex work space
            iw2 = iw1 + m
            iw3 = iw2 + m
            iwd = iw3 + m
            iww = iwd + m
            iwu = iww + m
            call compb (nl, ierror, an, bn, cn, w, wc, w(iwah), w(iwbh))
            return
        end if
        ! *** important to reset nm for np = 0
        if (npp == 0) nm = n - 1

        if (mp /= 0) then
            call blktr1 (nl, an, bn, cn, m, am, bm, cm, idimy, y, w, wc, w( &
                iw1), w(iw2), w(iw3), w(iwd), w(iww), w(iwu), wc(iw1), wc( &
                iw2), wc(iw3), prod, cprod)
        else
            call blktr1 (nl, an, bn, cn, m, am, bm, cm, idimy, y, w, wc, w( &
                iw1), w(iw2), w(iw3), w(iwd), w(iww), w(iwu), wc(iw1), wc( &
                iw2), wc(iw3), prodp, cprodp)
        end if

    end subroutine blktrii
    !
    !*****************************************************************************************
    !
    subroutine BLKTR1(N, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, B, BC, &
        W1, W2, W3, WD, WW, WU, CW1, CW2, CW3, PRDCT, CPRDCT)
        !-----------------------------------------------
        !   D u m m y   A r g u m e n t s
        !-----------------------------------------------
        integer (ip) :: N
        integer (ip) :: M
        integer (ip), intent (in) :: IDIMY
        real (wp) :: AN(*)
        real (wp) :: BN(*)
        real (wp) :: CN(*)
        real (wp) :: AM(*)
        real (wp) :: BM(*)
        real (wp) :: CM(*)
        real (wp) :: Y(IDIMY, 1)
        real (wp) :: B(*)
        real (wp) :: W1(*)
        real (wp) :: W2(*)
        real (wp) :: W3(*)
        real (wp) :: WD(*)
        real (wp) :: WW(*)
        real (wp) :: WU(*)
        complex (wp) :: BC(*)
        complex (wp) :: CW1(*)
        complex (wp) :: CW2(*)
        complex (wp) :: CW3(*)
        !-----------------------------------------------
        !   C o m m o n   B l o c k s
        !-----------------------------------------------
        !...  /CBLKT/
        common /CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
        integer   NPP, K, NM, NCMPLX, IK
        real   EPS, CNV
        !-----------------------------------------------
        !   L o c a l   V a r i a b l e s
        !-----------------------------------------------
        integer (ip) :: KDO, L, IR, I2, I1, I3, I4, IRM1, IM2, NM2, IM3, NM3, &
            IM1, NM1, I0, IF, I, IPI1, IPI2, IPI3, IDXC, NC, IDXA, NA, IP2 &
            , NP2, IP1, NP1, IP3, NP3, J, IZ, NZ, IZR, LL, IFD, IP, NP, &
            IMI1, IMI2
        real (wp) :: DUM
        !-----------------------------------------------
        !
        ! BLKTR1 SOLVES THE LINEAR SYSTEM
        !
        ! B  CONTAINS THE ROOTS OF ALL THE B POLYNOMIALS
        ! W1, W2, W3, WD, WW, WU  ARE ALL WORKING ARRAYS
        ! PRDCT  IS EITHER PRODP OR PROD DEPENDING ON WHETHER THE BOUNDARY
        ! CONDITIONS IN THE M DIRECTION ARE PERIODIC OR NOT
        ! CPRDCT IS EITHER CPRODP OR CPROD WHICH ARE THE COMPLEX VERSIONS
        ! OF PRODP AND PROD. THESE ARE CALLED IN THE EVENT THAT SOME
        ! OF THE ROOTS OF THE B SUB P POLYNOMIAL ARE COMPLEX
        !
        !
        !
        ! BEGIN REDUCTION PHASE
        !
        KDO = K - 1
        do L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I3 = I2 + I1
            I4 = I2 + I2
            IRM1 = IR - 1
            call INDXB (I2, IR, IM2, NM2)
            call INDXB (I1, IRM1, IM3, NM3)
            call INDXB (I3, IRM1, IM1, NM1)
            I0 = 0
            call PRDCT (NM2, B(IM2), NM3, B(IM3), NM1, B(IM1), I0, DUM, Y(1 &
                , I2), W3, M, AM, BM, CM, WD, WW, WU)
            IF = 2**K
            do I = I4, IF, I4
                if (I - NM > 0) cycle
                IPI1 = I + I1
                IPI2 = I + I2
                IPI3 = I + I3
                call INDXC (I, IR, IDXC, NC)
                if (I - IF >= 0) cycle
                call INDXA (I, IR, IDXA, NA)
                call INDXB (I - I1, IRM1, IM1, NM1)
                call INDXB (IPI2, IR, IP2, NP2)
                call INDXB (IPI1, IRM1, IP1, NP1)
                call INDXB (IPI3, IRM1, IP3, NP3)
                call PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W3, &
                    W1, M, AM, BM, CM, WD, WW, WU)
                if (IPI2 - NM > 0) then
                    W3(:M) = 0.
                    W2(:M) = 0.
                else
                    call PRDCT (NP2, B(IP2), NP1, B(IP1), NP3, B(IP3), 0, DUM &
                        , Y(1, IPI2), W3, M, AM, BM, CM, WD, WW, WU)
                    call PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W3 &
                        , W2, M, AM, BM, CM, WD, WW, WU)
                end if
                Y(:M, I) = W1(:M) + W2(:M) + Y(:M, I)
            end do
        end do
        if (NPP == 0) then
            IF = 2**K
            I = IF/2
            I1 = I/2
            call INDXB (I - I1, K - 2, IM1, NM1)
            call INDXB (I + I1, K - 2, IP1, NP1)
            call INDXB (I, K - 1, IZ, NZ)
            call PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, Y(1, I) &
                , W1, M, AM, BM, CM, WD, WW, WU)
            IZR = I
            W2(:M) = W1(:M)
            do LL = 2, K
                L = K - LL + 1
                IR = L - 1
                I2 = 2**IR
                I1 = I2/2
                I = I2
                call INDXC (I, IR, IDXC, NC)
                call INDXB (I, IR, IZ, NZ)
                call INDXB (I - I1, IR - 1, IM1, NM1)
                call INDXB (I + I1, IR - 1, IP1, NP1)
                call PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W1, &
                    W1, M, AM, BM, CM, WD, WW, WU)
                W1(:M) = Y(:M, I) + W1(:M)
                call PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W1 &
                    , W1, M, AM, BM, CM, WD, WW, WU)
            end do
            L118: do LL = 2, K
                L = K - LL + 1
                IR = L - 1
                I2 = 2**IR
                I1 = I2/2
                I4 = I2 + I2
                IFD = IF - I2
                do I = I2, IFD, I4
                    if (I - I2 - IZR /= 0) cycle
                    if (I - NM > 0) cycle  L118
                    call INDXA (I, IR, IDXA, NA)
                    call INDXB (I, IR, IZ, NZ)
                    call INDXB (I - I1, IR - 1, IM1, NM1)
                    call INDXB (I + I1, IR - 1, IP1, NP1)
                    call PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W2 &
                        , W2, M, AM, BM, CM, WD, WW, WU)
                    W2(:M) = Y(:M, I) + W2(:M)
                    call PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, &
                        W2, W2, M, AM, BM, CM, WD, WW, WU)
                    IZR = I
                    if (I - NM == 0) exit  L118
                end do
            end do L118
119     continue
        Y(:M, NM+1) = Y(:M, NM+1) - CN(NM+1)*W1(:M) - AN(NM+1)*W2(:M)
        call INDXB (IF/2, K - 1, IM1, NM1)
        call INDXB (IF, K - 1, IP, NP)
        if (NCMPLX /= 0) then
            call CPRDCT (NM + 1, BC(IP), NM1, B(IM1), 0, DUM, 0, DUM, Y( &
                1, NM+1), Y(1, NM+1), M, AM, BM, CM, CW1, CW2, CW3)
        else
            call PRDCT (NM + 1, B(IP), NM1, B(IM1), 0, DUM, 0, DUM, Y(1, &
                NM+1), Y(1, NM+1), M, AM, BM, CM, WD, WW, WU)
        end if
        W1(:M) = AN(1)*Y(:M, NM+1)
        W2(:M) = CN(NM)*Y(:M, NM+1)
        Y(:M, 1) = Y(:M, 1) - W1(:M)
        Y(:M, NM) = Y(:M, NM) - W2(:M)
        do L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I4 = I2 + I2
            I1 = I2/2
            I = I4
            call INDXA (I, IR, IDXA, NA)
            call INDXB (I - I2, IR, IM2, NM2)
            call INDXB (I - I2 - I1, IR - 1, IM3, NM3)
            call INDXB (I - I1, IR - 1, IM1, NM1)
            call PRDCT (NM2, B(IM2), NM3, B(IM3), NM1, B(IM1), 0, DUM, &
                W1, W1, M, AM, BM, CM, WD, WW, WU)
            call PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W1, &
                W1, M, AM, BM, CM, WD, WW, WU)
            Y(:M, I) = Y(:M, I) - W1(:M)
        end do
        !
        IZR = NM
        L131: do L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I3 = I2 + I1
            I4 = I2 + I2
            IRM1 = IR - 1
            do I = I4, IF, I4
                IPI1 = I + I1
                IPI2 = I + I2
                IPI3 = I + I3
                if (IPI2 - IZR /= 0) then
                    if (I - IZR /= 0) cycle
                    cycle  L131
                end if
                call INDXC (I, IR, IDXC, NC)
                call INDXB (IPI2, IR, IP2, NP2)
                call INDXB (IPI1, IRM1, IP1, NP1)
                call INDXB (IPI3, IRM1, IP3, NP3)
                call PRDCT (NP2, B(IP2), NP1, B(IP1), NP3, B(IP3), 0, DUM &
                    , W2, W2, M, AM, BM, CM, WD, WW, WU)
                call PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W2 &
                    , W2, M, AM, BM, CM, WD, WW, WU)
                Y(:M, I) = Y(:M, I) - W2(:M)
                IZR = I
                cycle  L131
            end do
        end do L131
    end if
    !
    ! BEGIN BACK SUBSTITUTION PHASE
    !
    do LL = 1, K
        L = K - LL + 1
        IR = L - 1
        IRM1 = IR - 1
        I2 = 2**IR
        I1 = I2/2
        I4 = I2 + I2
        IFD = IF - I2
        do I = I2, IFD, I4
            if (I - NM > 0) cycle
            IMI1 = I - I1
            IMI2 = I - I2
            IPI1 = I + I1
            IPI2 = I + I2
            call INDXA (I, IR, IDXA, NA)
            call INDXC (I, IR, IDXC, NC)
            call INDXB (I, IR, IZ, NZ)
            call INDXB (IMI1, IRM1, IM1, NM1)
            call INDXB (IPI1, IRM1, IP1, NP1)
            if (I - I2 <= 0) then
                W1(:M) = 0.
            else
                call PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), Y( &
                    1, IMI2), W1, M, AM, BM, CM, WD, WW, WU)
            end if
            if (IPI2 - NM > 0) then
                W2(:M) = 0.
            else
                call PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), Y( &
                    1, IPI2), W2, M, AM, BM, CM, WD, WW, WU)
            end if
            W1(:M) = Y(:M, I) + W1(:M) + W2(:M)
            call PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W1 &
                , Y(1, I), M, AM, BM, CM, WD, WW, WU)
        end do
    end do

end subroutine BLKTR1
    !
    !*****************************************************************************************
    !
real function BSRH (XLL, XRR, IZ, C, A, BH, F, SGN)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip) :: IZ
    real (wp), intent (in) :: XLL
    real (wp), intent (in) :: XRR
    real (wp) :: F
    real (wp), intent (in) :: SGN
    real (wp) :: C(*)
    real (wp) :: A(*)
    real (wp) :: BH(*)
    !-----------------------------------------------
    !   C o m m o n   B l o c k s
    !-----------------------------------------------
    !...  /CBLKT/
    common /CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
    integer   NPP, K, NM, NCMPLX, IK
    real   EPS, CNV
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    real (wp) :: R1, XL, XR, DX, X
    !-----------------------------------------------
    !   E x t e r n a l   F u n c t i o n s
    !-----------------------------------------------
    !-----------------------------------------------
    XL = XLL
    XR = XRR
    DX = 0.5*abs(XR - XL)
    X = 0.5*(XL + XR)
    R1 = SGN*F(X, IZ, C, A, BH)
    if (R1 >= 0.) then
        if (R1 == 0.) go to 105
        XR = X
    else
        XL = X
    end if
    DX = 0.5*DX
    do while(DX - CNV > 0.)
        X = 0.5*(XL + XR)
        R1 = SGN*F(X, IZ, C, A, BH)
        if (R1 >= 0.) then
            if (R1 == 0.) go to 105
            XR = X
        else
            XL = X
        end if
        DX = 0.5*DX
    end do
105 continue
    BSRH = 0.5*(XL + XR)

end function BSRH
    !
    !*****************************************************************************************
    !
subroutine COMPB(N, IERROR, AN, BN, CN, B, BC, AH, BH)
    real (wp) :: epmach
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip) :: N
    integer (ip) :: IERROR
    real (wp) :: AN(*)
    real (wp), intent (in) :: BN(*)
    real (wp) :: CN(*)
    real (wp) :: B(*)
    real (wp) :: AH(*)
    real (wp) :: BH(*)
    complex (wp) :: BC(*)
    !-----------------------------------------------
    !   C o m m o n   B l o c k s
    !-----------------------------------------------
    !...  /CBLKT/
    common /CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
    integer   NPP, K, NM, NCMPLX, IK
    real   EPS, CNV
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer (ip) :: J, IF, KDO, L, IR, I2, I4, IPL, IFD, I, IB, NB, JS, JF
    integer (ip) ::  LS, LH, NMP, L1, L2, J2, J1, N2M2
    real (wp) :: DUM, BNORM, ARG, D1, D2, D3
    !-----------------------------------------------
    !
    !     COMPB COMPUTES THE ROOTS OF THE B POLYNOMIALS USING SUBROUTINE
    !     TEVLS WHICH IS A MODIFICATION THE EISPACK PROGRAM TQLRAT.
    !     IERROR IS SET TO 4 IF EITHER TEVLS FAILS OR IF A(J+1)*C(J) IS
    !     LESS THAN ZERO FOR SOME J.  AH, BH ARE TEMPORARY WORK ARRAYS.
    !
    EPS = epsilon(DUM)
    BNORM = abs(BN(1))
    do J = 2, NM
        BNORM = max(BNORM, abs(BN(J)))
        ARG = AN(J)*CN(J-1)
        if (ARG < 0.) go to 119
        B(J) = SIGN(SQRT(ARG), AN(J))
    end do
    CNV = EPS*BNORM
    IF = 2**K
    KDO = K - 1
    L108: do L = 1, KDO
        IR = L - 1
        I2 = 2**IR
        I4 = I2 + I2
        IPL = I4 - 1
        IFD = IF - I4
        do I = I4, IFD, I4
            call INDXB (I, L, IB, NB)
            if (NB <= 0) cycle  L108
            JS = I - IPL
            JF = JS + NB - 1
            LS = 0
            BH(:JF-JS+1) = BN(JS:JF)
            AH(:JF-JS+1) = B(JS:JF)
            call TEVLS (NB, BH, AH, IERROR)
            if (IERROR /= 0) go to 118
            LH = IB - 1
            if (NB > 0) then
                B(LH+1:NB+LH) = -BH(:NB)
                LH = NB + LH
            end if
        end do
    end do L108
    B(:NM) = -BN(:NM)
    if (NPP == 0) then
        NMP = NM + 1
        NB = NM + NMP
        do J = 1, NB
            L1 = mod(J - 1, NMP) + 1
            L2 = mod(J + NM - 1, NMP) + 1
            ARG = AN(L1)*CN(L2)
            if (ARG < 0.) go to 119
            BH(J) = SIGN(SQRT(ARG), (-AN(L1)))
            AH(J) = -BN(L1)
        end do
        call TEVLS (NB, AH, BH, IERROR)
        if (IERROR /= 0) go to 118
        call INDXB (IF, K - 1, J2, LH)
        call INDXB (IF/2, K - 1, J1, LH)
        J2 = J2 + 1
        LH = J2
        N2M2 = J2 + NM + NM - 2
114 continue
    D1 = abs(B(J1)-B(J2-1))
    D2 = abs(B(J1)-B(J2))
    D3 = abs(B(J1)-B(J2+1))
    if (D2>=D1 .or. D2>=D3) then
        B(LH) = B(J2)
        J2 = J2 + 1
        LH = LH + 1
        if (J2 - N2M2 <= 0) go to 114
    else
        J2 = J2 + 1
        J1 = J1 + 1
        if (J2 - N2M2 <= 0) go to 114
    end if
    B(LH) = B(N2M2+1)
    call INDXB (IF, K - 1, J1, J2)
    J2 = J1 + NMP + NMP
    call PPADD (NM + 1, IERROR, AN, CN, BC(J1), B(J1), B(J2))
end if
return 
118 continue
    IERROR = 4
    return
119 continue
    IERROR = 5

end subroutine COMPB
    !
    !*****************************************************************************************
    !
subroutine CPROD(ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, YY, M, A, B, C, D, W, Y)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip), intent (in) :: ND
    integer (ip), intent (in) :: NM1
    integer (ip), intent (in) :: NM2
    integer (ip), intent (in) :: NA
    integer (ip), intent (in) :: M
    real (wp), intent (in) :: BM1(*)
    real (wp), intent (in) :: BM2(*)
    real (wp), intent (in) :: AA(*)
    real (wp), intent (in) :: X(*)
    real (wp), intent (out) :: YY(*)
    real (wp), intent (in) :: A(*)
    real (wp), intent (in) :: B(*)
    real (wp), intent (in) :: C(*)
    complex , intent (in) :: BD(*)
    complex , intent (in out) :: D(*)
    complex , intent (in out) :: W(*)
    complex , intent (in out) :: Y(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer (ip) :: J, MM, ID, M1, M2, IA, IFLG, K
    real (wp) :: RT
    complex :: CRT, DEN, Y1, Y2
    !-----------------------------------------------
    !
    ! PROD APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
    ! STORES THE RESULT IN YY           (COMPLEX CASE)
    ! AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
    ! ND, NM1, NM2 ARE THE LENGTHS OF THE ARRAYS BD, BM1, BM2 RESPECTIVELY
    ! BD, BM1, BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
    ! NA IS THE LENGTH OF THE ARRAY AA
    ! X, YY THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS YY
    ! A, B, C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
    ! M  IS THE ORDER OF THE MATRIX
    ! D, W, Y ARE WORKING ARRAYS
    ! ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
    !
    do J = 1, M
        Y(J) = CMPLX(X(J), 0.)
    end do
    MM = M - 1
    ID = ND
    M1 = NM1
    M2 = NM2
    IA = NA
102 continue
    IFLG = 0
    if (ID > 0) then
        CRT = BD(ID)
        ID = ID - 1
        !
        ! BEGIN SOLUTION TO SYSTEM
        !
        D(M) = A(M)/(B(M)-CRT)
        W(M) = Y(M)/(B(M)-CRT)
        do J = 2, MM
            K = M - J
            DEN = B(K+1) - CRT - C(K+1)*D(K+2)
            D(K+1) = A(K+1)/DEN
            W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
        end do
        DEN = B(1) - CRT - C(1)*D(2)
        if (abs(DEN) /= 0.) then
            Y(1) = (Y(1)-C(1)*W(2))/DEN
        else
            Y(1) = (1., 0.)
        end if
        do J = 2, M
            Y(J) = W(J) - D(J)*Y(J-1)
        end do
    end if
    if (M1 <= 0) then
        if (M2 <= 0) go to 121
        RT = BM2(M2)
        M2 = M2 - 1
    else
        if (M2 <= 0) then
            RT = BM1(M1)
            M1 = M1 - 1
        else
            if (abs(BM1(M1)) - abs(BM2(M2)) > 0.) then
                RT = BM1(M1)
                M1 = M1 - 1
            else
                RT = BM2(M2)
                M2 = M2 - 1
            end if
        end if
    end if
    Y1 = (B(1)-RT)*Y(1) + C(1)*Y(2)
    if (MM - 2 >= 0) then
        do J = 2, MM
            Y2 = A(J)*Y(J-1) + (B(J)-RT)*Y(J) + C(J)*Y(J+1)
            Y(J-1) = Y1
            Y1 = Y2
        end do
    end if
    Y(M) = A(M)*Y(M-1) + (B(M)-RT)*Y(M)
    Y(M-1) = Y1
    IFLG = 1
    go to 102
121 continue
    if (IA > 0) then
        RT = AA(IA)
        IA = IA - 1
        IFLG = 1
        !
        ! SCALAR MULTIPLICATION
        !
        Y(:M) = RT*Y(:M)
    end if
    if (IFLG > 0) go to 102
    do J = 1, M
        YY(J) = REAL(Y(J))
    end do
    return
end subroutine CPROD
    !
    !*****************************************************************************************
    !
subroutine CPRODP(ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, YY, M, A &
    , B, C, D, U, Y)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip), intent (in) :: ND
    integer (ip), intent (in) :: NM1
    integer (ip), intent (in) :: NM2
    integer (ip), intent (in) :: NA
    integer (ip), intent (in) :: M
    real (wp), intent (in) :: BM1(*)
    real (wp), intent (in) :: BM2(*)
    real (wp), intent (in) :: AA(*)
    real (wp), intent (in) :: X(*)
    real (wp), intent (out) :: YY(*)
    real (wp), intent (in) :: A(*)
    real (wp), intent (in) :: B(*)
    real (wp), intent (in) :: C(*)
    complex , intent (in) :: BD(*)
    complex , intent (in out) :: D(*)
    complex , intent (in out) :: U(*)
    complex , intent (in out) :: Y(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer (ip) :: J, MM, MM2, ID, M1, M2, IA, IFLG, K
    real (wp) :: RT
    complex :: V, DEN, BH, YM, AM, Y1, Y2, YH, CRT
    !-----------------------------------------------
    !
    ! PRODP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
    ! STORES THE RESULT IN YY       PERIODIC BOUNDARY CONDITIONS
    ! AND  COMPLEX  CASE
    !
    ! BD, BM1, BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
    ! ND, NM1, NM2 ARE THE LENGTHS OF THE ARRAYS BD, BM1, BM2 RESPECTIVELY
    ! AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
    ! NA IS THE LENGTH OF THE ARRAY AA
    ! X, YY THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS YY
    ! A, B, C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
    ! M  IS THE ORDER OF THE MATRIX
    ! D, U, Y ARE WORKING ARRAYS
    ! ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
    !
    do J = 1, M
        Y(J) = CMPLX(X(J), 0.)
    end do
    MM = M - 1
    MM2 = M - 2
    ID = ND
    M1 = NM1
    M2 = NM2
    IA = NA
102 continue
    IFLG = 0
    if (ID > 0) then
        CRT = BD(ID)
        ID = ID - 1
        IFLG = 1
        !
        ! BEGIN SOLUTION TO SYSTEM
        !
        BH = B(M) - CRT
        YM = Y(M)
        DEN = B(1) - CRT
        D(1) = C(1)/DEN
        U(1) = A(1)/DEN
        Y(1) = Y(1)/DEN
        V = CMPLX(C(M), 0.)
        if (MM2 - 2 >= 0) then
            do J = 2, MM2
                DEN = B(J) - CRT - A(J)*D(J-1)
                D(J) = C(J)/DEN
                U(J) = -A(J)*U(J-1)/DEN
                Y(J) = (Y(J)-A(J)*Y(J-1))/DEN
                BH = BH - V*U(J-1)
                YM = YM - V*Y(J-1)
                V = -V*D(J-1)
            end do
        end if
        DEN = B(M-1) - CRT - A(M-1)*D(M-2)
        D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
        Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/DEN
        AM = A(M) - V*D(M-2)
        BH = BH - V*U(M-2)
        YM = YM - V*Y(M-2)
        DEN = BH - AM*D(M-1)
        if (abs(DEN) /= 0.) then
            Y(M) = (YM - AM*Y(M-1))/DEN
        else
            Y(M) = (1., 0.)
        end if
        Y(M-1) = Y(M-1) - D(M-1)*Y(M)
        do J = 2, MM
            K = M - J
            Y(K) = Y(K) - D(K)*Y(K+1) - U(K)*Y(M)
        end do
    end if
    if (M1 <= 0) then
        if (M2 <= 0) go to 123
        RT = BM2(M2)
        M2 = M2 - 1
    else
        if (M2 <= 0) then
            RT = BM1(M1)
            M1 = M1 - 1
        else
            if (abs(BM1(M1)) - abs(BM2(M2)) > 0.) then
                RT = BM1(M1)
                M1 = M1 - 1
            else
                RT = BM2(M2)
                M2 = M2 - 1
            !
            ! MATRIX MULTIPLICATION
            !
            end if
        end if
    end if
    YH = Y(1)
    Y1 = (B(1)-RT)*Y(1) + C(1)*Y(2) + A(1)*Y(M)
    if (MM - 2 >= 0) then
        do J = 2, MM
            Y2 = A(J)*Y(J-1) + (B(J)-RT)*Y(J) + C(J)*Y(J+1)
            Y(J-1) = Y1
            Y1 = Y2
        end do
    end if
    Y(M) = A(M)*Y(M-1) + (B(M)-RT)*Y(M) + C(M)*YH
    Y(M-1) = Y1
    IFLG = 1
    go to 102
123 continue
    if (IA > 0) then
        RT = AA(IA)
        IA = IA - 1
        IFLG = 1
        !
        ! SCALAR MULTIPLICATION
        !
        Y(:M) = RT*Y(:M)
    end if
    if (IFLG > 0) go to 102
    do J = 1, M
        YY(J) = REAL(Y(J))
    end do

end subroutine CPRODP
    !
    !*****************************************************************************************
    !
subroutine INDXA(I, IR, IDXA, NA)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip), intent (in) :: I
    integer (ip), intent (in) :: IR
    integer (ip), intent (out) :: IDXA
    integer (ip), intent (out) :: NA
    !-----------------------------------------------
    !   C o m m o n   B l o c k s
    !-----------------------------------------------
    !...  /CBLKT/
    common /CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
    integer   NPP, K, NM, NCMPLX, IK
    real   EPS, CNV
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    !-----------------------------------------------
    NA = 2**IR
    IDXA = I - NA + 1
    if (I - NM > 0) then
        NA = 0
    end if
    return
end subroutine INDXA
    !
    !*****************************************************************************************
    !
subroutine INDXB(I, IR, IDX, IDP)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip), intent (in) :: I
    integer (ip), intent (in) :: IR
    integer (ip), intent (out) :: IDX
    integer (ip), intent (out) :: IDP
    !-----------------------------------------------
    !   C o m m o n   B l o c k s
    !-----------------------------------------------
    !...  /CBLKT/
    common /CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
    integer   NPP, K, NM, NCMPLX, IK
    real   EPS, CNV
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer (ip) :: IZH, ID, IPL
    !-----------------------------------------------
    !
    ! B(IDX) IS THE LOCATION OF THE FIRST ROOT OF THE B(I, IR) POLYNOMIAL
    !
    IDP = 0
    if (IR >= 0) then
        if (IR <= 0) then
            if (I - NM > 0) go to 107
            IDX = I
            IDP = 1
            return
        end if
        IZH = 2**IR
        ID = I - IZH - IZH
        IDX = ID + ID + (IR - 1)*IK + IR + (IK - I)/IZH + 4
        IPL = IZH - 1
        IDP = IZH + IZH - 1
        if (I - IPL - NM > 0) then
            IDP = 0
            return
        end if
        if (I + IPL - NM > 0) then
            IDP = NM + IPL - I + 1
        end if
    end if
107 continue

end subroutine INDXB
    !
    !*****************************************************************************************
    !
subroutine INDXC(I, IR, IDXC, NC)

    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip), intent (in) :: I
    integer (ip), intent (in) :: IR
    integer (ip), intent (out) :: IDXC
    integer (ip), intent (out) :: NC
    !-----------------------------------------------
    !   C o m m o n   B l o c k s
    !-----------------------------------------------
    !...  /CBLKT/
    common /CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
    integer   NPP, K, NM, NCMPLX, IK
    real   EPS, CNV
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    !-----------------------------------------------
    NC = 2**IR
    IDXC = I
    if (IDXC + NC - 1 - NM > 0) then
        NC = 0
    end if
    return
end subroutine INDXC
    !
    !*****************************************************************************************
    !
subroutine PPADD(N, IERROR, A, C, CBP, BP, BH)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip), intent (in) :: N
    integer (ip), intent (out) :: IERROR
    real (wp) :: A(*)
    real (wp) :: C(*)
    real (wp), intent (in out) :: BP(*)
    real (wp) :: BH(*)
    complex , intent (in out) :: CBP(*)
    !-----------------------------------------------
    !   C o m m o n   B l o c k s
    !-----------------------------------------------
    !...  /CBLKT/
    common /CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
    integer   NPP, K, NM, NCMPLX, IK
    real   EPS, CNV
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer::IZ, IZM, IZM2, J, NT, MODIZ, IS, IF, IG, IT, ICV, I3, I2, NHALF
    real (wp) :: R4, R5, R6, SCNV, XL, DB, SGN, XR, XM, PSG
    complex :: CF, CX, FSG, HSG, DD, F, FP, FPP, CDIS, R1, R2, R3
    !-----------------------------------------------
    !
    !     PPADD COMPUTES THE EIGENVALUES OF THE PERIODIC TRIDIAGONAL MATRIX
    !     WITH COEFFICIENTS AN, BN, CN
    !
    ! N IS THE ORDER OF THE BH AND BP POLYNOMIALS
    !     ON OUTPUT BP CONTIANS THE EIGENVALUES
    ! CBP IS THE SAME AS BP EXCEPT TYPE COMPLEX
    ! BH IS USED TO TEMPORARILY STORE THE ROOTS OF THE B HAT POLYNOMIAL
    ! WHICH ENTERS THROUGH BP
    !
    SCNV = SQRT(CNV)
    IZ = N
    IZM = IZ - 1
    IZM2 = IZ - 2
    if (BP(N) - BP(1) <= 0.) then
        if (BP(N) - BP(1) == 0.) go to 142
        BH(:N) = BP(N:1:(-1))
    else
        BH(:N) = BP(:N)
    end if
    NCMPLX = 0
    MODIZ = mod(IZ, 2)
    IS = 1
    if (MODIZ /= 0) then
        if (A(1) < 0.) go to 110
        if (A(1) == 0.) go to 142
    end if
    XL = BH(1)
    DB = BH(3) - BH(1)
    XL = XL - DB
    R4 = PSGF(XL, IZ, C, A, BH)
    do while(R4 <= 0.)
        XL = XL - DB
        R4 = PSGF(XL, IZ, C, A, BH)
    end do
    SGN = -1.
    CBP(1) = CMPLX(BSRH(XL, BH(1), IZ, C, A, BH, PSGF, SGN), 0.)
    BP(1) = REAL(CBP(1))
    IS = 2
110 continue
    IF = IZ - 1
    if (MODIZ /= 0) then
        if (A(1) > 0.) go to 115
        if (A(1) == 0.) go to 142
    end if
    XR = BH(IZ)
    DB = BH(IZ) - BH(IZ-2)
    XR = XR + DB
    R5 = PSGF(XR, IZ, C, A, BH)
    do while(R5 < 0.)
        XR = XR + DB
        R5 = PSGF(XR, IZ, C, A, BH)
    end do
    SGN = 1.
    CBP(IZ) = CMPLX(BSRH(BH(IZ), XR, IZ, C, A, BH, PSGF, SGN), 0.)
    IF = IZ - 2
115 continue
    do IG = IS, IF, 2
        XL = BH(IG)
        XR = BH(IG+1)
        SGN = -1.
        XM = BSRH(XL, XR, IZ, C, A, BH, PPSPF, SGN)
        PSG = PSGF(XM, IZ, C, A, BH)
        if (abs(PSG) - EPS <= 0.) go to 118
        R6 = PSG*PPSGF(XM, IZ, C, A, BH)
        if (R6 > 0.) go to 119
        if (R6 == 0.) go to 118
        SGN = 1.
        CBP(IG) = CMPLX(BSRH(BH(IG), XM, IZ, C, A, BH, PSGF, SGN), 0.)
        !        bp(ig) = real(cbp(ig))
        SGN = -1.
        CBP(IG+1) = CMPLX(BSRH(XM, BH(IG+1), IZ, C, A, BH, PSGF, SGN), 0.)
        !        bp(ig) = real(cbp(ig))
        !        bp(ig+1) = real(cbp(ig+1))
        cycle
    !
    !     CASE OF A MULTIPLE ZERO
    !
118 continue
    CBP(IG) = CMPLX(XM, 0.)
    CBP(IG+1) = CMPLX(XM, 0.)
    !        bp(ig) = real(cbp(ig))
    !        bp(ig+1) = real(cbp(ig+1))
    cycle
!
!     CASE OF A COMPLEX ZERO
!
119 continue
    IT = 0
    ICV = 0
    CX = CMPLX(XM, 0.)
120 continue
    FSG = (1., 0.)
    HSG = (1., 0.)
    FP = (0., 0.)
    FPP = (0., 0.)
    do J = 1, IZ
        DD = 1./(CX - BH(J))
        FSG = FSG*A(J)*DD
        HSG = HSG*C(J)*DD
        FP = FP + DD
        FPP = FPP - DD*DD
    end do
    if (MODIZ == 0) then
        F = (1., 0.) - FSG - HSG
    else
        F = (1., 0.) + FSG + HSG
    end if
    I3 = 0
    if (abs(FP) > 0.) then
        I3 = 1
        R3 = -F/FP
    end if
    I2 = 0
    if (abs(FPP) > 0.) then
        I2 = 1
        CDIS = CSQRT(FP**2 - 2.*F*FPP)
        R1 = CDIS - FP
        R2 = (-FP) - CDIS
        if (abs(R1) - abs(R2) > 0.) then
            R1 = R1/FPP
        else
            R1 = R2/FPP
        end if
        R2 = 2.*F/FPP/R1
        if (abs(R2) < abs(R1)) R1 = R2
        if (I3 <= 0) go to 133
        if (abs(R3) < abs(R1)) R1 = R3
        go to 133
    end if
    R1 = R3
133 continue
    CX = CX + R1
    IT = IT + 1
    if (IT > 50) go to 142
    if (abs(R1) > SCNV) go to 120
    if (ICV > 0) go to 135
    ICV = 1
    go to 120
135 continue
    CBP(IG) = CX
    CBP(IG+1) = CONJG(CX)
end do
if (abs(CBP(N)) - abs(CBP(1)) <= 0.) then
    if (abs(CBP(N)) - abs(CBP(1)) == 0.) go to 142
    NHALF = N/2
    do J = 1, NHALF
        NT = N - J
        CX = CBP(J)
        CBP(J) = CBP(NT+1)
        CBP(NT+1) = CX
    end do
end if
NCMPLX = 1
do J = 2, IZ
    if (AIMAG(CBP(J)) /= 0.) go to 143
end do
NCMPLX = 0
BP(1) = REAL(CBP(1))
do J = 2, IZ
    BP(J) = REAL(CBP(J))
end do
go to 143
142 continue
    IERROR = 4
143 continue

end subroutine PPADD
    !
    !*****************************************************************************************
    !
subroutine PROD(ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, Y, M, A, B, C, D, W, U)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip), intent (in) :: ND
    integer (ip), intent (in) :: NM1
    integer (ip), intent (in) :: NM2
    integer (ip), intent (in) :: NA
    integer (ip), intent (in) :: M
    real (wp), intent (in) :: BD(*)
    real (wp), intent (in) :: BM1(*)
    real (wp), intent (in) :: BM2(*)
    real (wp), intent (in) :: AA(*)
    real (wp), intent (in) :: X(*)
    real (wp), intent (in out) :: Y(*)
    real (wp), intent (in) :: A(*)
    real (wp), intent (in) :: B(*)
    real (wp), intent (in) :: C(*)
    real (wp), intent (in out) :: D(*)
    real (wp), intent (in out) :: W(*)
    real (wp) :: U(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer (ip) :: J, MM, ID, IBR, M1, M2, IA, K
    real (wp) :: RT, DEN
    !-----------------------------------------------
    !
    ! PROD APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
    ! STORES THE RESULT IN Y
    ! BD, BM1, BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
    ! ND, NM1, NM2 ARE THE LENGTHS OF THE ARRAYS BD, BM1, BM2 RESPECTIVELY
    ! AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
    ! NA IS THE LENGTH OF THE ARRAY AA
    ! X, Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
    ! A, B, C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
    ! M  IS THE ORDER OF THE MATRIX
    ! D, W, U ARE WORKING ARRAYS
    ! IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
    !
    W(:M) = X(:M)
    Y(:M) = W(:M)
    MM = M - 1
    ID = ND
    IBR = 0
    M1 = NM1
    M2 = NM2
    IA = NA
102 continue
    if (IA > 0) then
        RT = AA(IA)
        if (ND == 0) RT = -RT
        IA = IA - 1
        !
        ! SCALAR MULTIPLICATION
        !
        Y(:M) = RT*W(:M)
    end if
    if (ID <= 0) go to 125
    RT = BD(ID)
    ID = ID - 1
    if (ID == 0) IBR = 1
    !
    ! BEGIN SOLUTION TO SYSTEM
    !
    D(M) = A(M)/(B(M)-RT)
    W(M) = Y(M)/(B(M)-RT)
    do J = 2, MM
        K = M - J
        DEN = B(K+1) - RT - C(K+1)*D(K+2)
        D(K+1) = A(K+1)/DEN
        W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
    end do
    DEN = B(1) - RT - C(1)*D(2)
    W(1) = 1.
    if (DEN /= 0.) then
        W(1) = (Y(1)-C(1)*W(2))/DEN
    end if
    do J = 2, M
        W(J) = W(J) - D(J)*W(J-1)
    end do
    if (NA > 0) go to 102
    go to 113
111 continue
    Y(:M) = W(:M)
    IBR = 1
    go to 102
113 continue
    if (M1 <= 0) then
        if (M2 <= 0) go to 111
    else
        if (M2 > 0) then
            if (abs(BM1(M1)) - abs(BM2(M2)) <= 0.) go to 120
        end if
        if (IBR <= 0) then
            if (abs(BM1(M1)-BD(ID)) - abs(BM1(M1)-RT) < 0.) go to 111
        end if
        RT = RT - BM1(M1)
        M1 = M1 - 1
        go to 123
    end if
120 continue
    if (IBR <= 0) then
        if (abs(BM2(M2)-BD(ID)) - abs(BM2(M2)-RT) < 0.) go to 111
    end if
    RT = RT - BM2(M2)
    M2 = M2 - 1
123 continue
    Y(:M) = Y(:M) + RT*W(:M)
    go to 102
125 continue

end subroutine PROD
    !
    !*****************************************************************************************
    !
subroutine PRODP(ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, Y, M, A, B, C, D, U, W)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip), intent (in) :: ND
    integer (ip), intent (in) :: NM1
    integer (ip), intent (in) :: NM2
    integer (ip), intent (in) :: NA
    integer (ip), intent (in) :: M
    real (wp), intent (in) :: BD(*)
    real (wp), intent (in) :: BM1(*)
    real (wp), intent (in) :: BM2(*)
    real (wp), intent (in) :: AA(*)
    real (wp), intent (in) :: X(*)
    real (wp), intent (in out) :: Y(*)
    real (wp), intent (in) :: A(*)
    real (wp), intent (in) :: B(*)
    real (wp), intent (in) :: C(*)
    real (wp), intent (in out) :: D(*)
    real (wp), intent (in out) :: U(*)
    real (wp), intent (in out) :: W(*)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer (ip) :: J, MM, MM2, ID, IBR, M1, M2, IA, K
    real (wp) :: RT, BH, YM, DEN, V, AM
    !-----------------------------------------------
    !
    ! PRODP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
    ! STORES THE RESULT IN Y        PERIODIC BOUNDARY CONDITIONS
    !
    ! BD, BM1, BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
    ! ND, NM1, NM2 ARE THE LENGTHS OF THE ARRAYS BD, BM1, BM2 RESPECTIVELY
    ! AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
    ! NA IS THE LENGTH OF THE ARRAY AA
    ! X, Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
    ! A, B, C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
    ! M  IS THE ORDER OF THE MATRIX
    ! D, U, W ARE WORKING ARRAYS
    ! IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
    !

    Y(:M) = X(:M)
    W(:M) = Y(:M)
    MM = M - 1
    MM2 = M - 2
    ID = ND
    IBR = 0
    M1 = NM1
    M2 = NM2
    IA = NA
102 continue
    if (IA > 0) then
        RT = AA(IA)
        if (ND == 0) RT = -RT
        IA = IA - 1
        Y(:M) = RT*W(:M)
    end if
    if (ID <= 0) go to 128
    RT = BD(ID)
    ID = ID - 1
    if (ID == 0) IBR = 1
    !
    ! BEGIN SOLUTION TO SYSTEM
    !
    BH = B(M) - RT
    YM = Y(M)
    DEN = B(1) - RT
    D(1) = C(1)/DEN
    U(1) = A(1)/DEN
    W(1) = Y(1)/DEN
    V = C(M)
    if (MM2 - 2 >= 0) then
        do J = 2, MM2
            DEN = B(J) - RT - A(J)*D(J-1)
            D(J) = C(J)/DEN
            U(J) = -A(J)*U(J-1)/DEN
            W(J) = (Y(J)-A(J)*W(J-1))/DEN
            BH = BH - V*U(J-1)
            YM = YM - V*W(J-1)
            V = -V*D(J-1)
        end do
    end if
    DEN = B(M-1) - RT - A(M-1)*D(M-2)
    D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
    W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/DEN
    AM = A(M) - V*D(M-2)
    BH = BH - V*U(M-2)
    YM = YM - V*W(M-2)
    DEN = BH - AM*D(M-1)
    if (DEN /= 0.) then
        W(M) = (YM - AM*W(M-1))/DEN
    else
        W(M) = 1.
    end if
    W(M-1) = W(M-1) - D(M-1)*W(M)
    do J = 2, MM
        K = M - J
        W(K) = W(K) - D(K)*W(K+1) - U(K)*W(M)
    end do
    if (NA > 0) go to 102
    go to 116
114 continue
    Y(:M) = W(:M)
    IBR = 1
    go to 102
116 continue
    if (M1 <= 0) then
        if (M2 <= 0) go to 114
    else
        if (M2 > 0) then
            if (abs(BM1(M1)) - abs(BM2(M2)) <= 0.) go to 123
        end if
        if (IBR <= 0) then
            if (abs(BM1(M1)-BD(ID)) - abs(BM1(M1)-RT) < 0.) go to 114
        end if
        RT = RT - BM1(M1)
        M1 = M1 - 1
        go to 126
    end if
123 continue
    if (IBR <= 0) then
        if (abs(BM2(M2)-BD(ID)) - abs(BM2(M2)-RT) < 0.) go to 114
    end if
    RT = RT - BM2(M2)
    M2 = M2 - 1
126 continue
    Y(:M) = Y(:M) + RT*W(:M)
    go to 102
128 continue
    return
end subroutine PRODP
    !
    !*****************************************************************************************
    !
subroutine TEVLS(N, D, E2, IERR)
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer (ip), intent (in) :: N
    integer (ip), intent (out) :: IERR
    real (wp), intent (in out) :: D(N)
    real (wp), intent (in out) :: E2(N)
    !-----------------------------------------------
    !   C o m m o n   B l o c k s
    !-----------------------------------------------
    !...  /CBLKT/
    common /CBLKT/ NPP, K, MACHEP, CNV, NM, NCMPLX, IK
    integer   NPP, K, NM, NCMPLX, IK
    real   MACHEP, CNV
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer (ip) :: I, J, L, M, II, L1, MML, NHALF, NTOP
    real (wp) :: B, C, F, G, H, P, R, S, DHOLD
    !-----------------------------------------------
    !
    !
    !     REAL SQRT, ABS, SIGN
    !
    !
    !     THIS SUBROUTINE IS A MODIFICATION OF THE EISPACK SUBROUTINE TQLRAT
    !     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
    !
    !     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
    !     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
    !
    !     ON INPUT-
    !
    !        N IS THE ORDER OF THE MATRIX, 
    !
    !        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX, 
    !
    !        E2 CONTAINS THE                SUBDIAGONAL ELEMENTS OF THE
    !          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
    !
    !      ON OUTPUT-
    !
    !        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
    !          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
    !          ORDERED FOR INDICES 1, 2, ...IERR-1, BUT MAY NOT BE
    !          THE SMALLEST EIGENVALUES, 
    !
    !        E2 HAS BEEN DESTROYED, 
    !
    !        IERR IS SET TO
    !          ZERO       FOR NORMAL RETURN, 
    !          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
    !                     DETERMINED AFTER 30 ITERATIONS.
    !
    !     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, 
    !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
    !
    !
    !     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
    !                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
    !
    !                **********
    !
    IERR = 0
    if (N /= 1) then
        !
        E2(:N-1) = E2(2:N)*E2(2:N)
        !
        F = 0.0
        B = 0.0
        E2(N) = 0.0
        !
        do L = 1, N
            J = 0
            H = MACHEP*(abs(D(L))+SQRT(E2(L)))
            if (B <= H) then
                B = H
                C = B*B
            end if
            !
            !     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
            !
            do M = L, N
                if (E2(M) > C) cycle
                exit
            !
            !     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
            !                THROUGH THE BOTTOM OF THE LOOP **********
            !
            end do
            !
            if (M /= L) then
105         continue
            if (J == 30) go to 114
            J = J + 1
            !
            !     ********** FORM SHIFT **********
            !
            L1 = L + 1
            S = SQRT(E2(L))
            G = D(L)
            P = (D(L1)-G)/(2.0*S)
            R = SQRT(P*P + 1.0)
            D(L) = S/(P + SIGN(R, P))
            H = G - D(L)
            !
            D(L1:N) = D(L1:N) - H
            !
            F = F + H
            !
            !     ********** RATIONAL QL TRANSFORMATION **********
            !
            G = D(M)
            if (G == 0.0) G = B
            H = G
            S = 0.0
            MML = M - L
            !
            !     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
            !
            do II = 1, MML
                I = M - II
                P = G*H
                R = P + E2(I)
                E2(I+1) = S*R
                S = E2(I)/R
                D(I+1) = H + S*(H + D(I))
                G = D(I) - E2(I)/G
                if (G == 0.0) G = B
                H = G*P/R
            end do
            !
            E2(L) = S*G
            D(L) = H
            !
            !     ********** GUARD AGAINST UNDERFLOWED H **********
            !
            if (H == 0.0) go to 108
            if (abs(E2(L)) <= abs(C/H)) go to 108
            E2(L) = H*E2(L)
            if (E2(L) /= 0.0) go to 105
        end if
108 continue
    P = D(L) + F
    !
    !     ********** ORDER EIGENVALUES **********
    !
    if (L /= 1) then
        !
        !     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
        !
        do II = 2, L
            I = L + 2 - II
            if (P >= D(I-1)) go to 111
            D(I) = D(I-1)
        end do
    end if
    !
    I = 1
111 continue
    D(I) = P
end do
!
if (abs(D(N)) >= abs(D(1))) go to 115
NHALF = N/2
do I = 1, NHALF
    NTOP = N - I
    DHOLD = D(I)
    D(I) = D(NTOP+1)
    D(NTOP+1) = DHOLD
end do
go to 115
!
!     ********** SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS **********
!
114 continue
    IERR = L
end if
115 continue
    return
!
!     ********** LAST CARD OF TQLRAT **********
!
end subroutine TEVLS
    !
    !*****************************************************************************************
    !
end module module_blktri
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
