module module_blktri

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_comf, only: &
        psgf, &
        ppspf, &
        ppsgf, &
        comf_interface

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: blktri
    public :: blktrii
    public :: test_blktri

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    integer (ip), save :: npp, k, nm, ncmplx, ik
    real (wp),    save :: eps, cnv
    !---------------------------------------------------------------------------------

contains


    subroutine test_blktri()
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
            s(i) = real(i, kind=wp) * ds
        end do

        dt = 1.0_wp/(n + 1)

        do j = 1, n
            t(j) = real(j, kind=wp) * dt
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
            call blktri (iflg, np, n, an, bn, cn, mp, m, am, bm, cm, idimy, &
                y, ierror, workspace )
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
        write( stdout, '(A)') '     blktri *** TEST RUN *** '
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') '     ierror = 0,  discretization error = 1.6478E-05'
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,I3,A,1pe15.6)') '     ierror =', ierror, ' discretization error = ', &
            discretization_error

        ! Release memory
        call workspace%destroy()

    end subroutine test_blktri


    subroutine blktri( iflg, np, n, an, bn, cn, mp, m, am, bm, cm, &
        idimy, y, ierror, workspace )
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
        !     SUBROUTINE blktri (IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, IDIMY, Y,
        !    +                   ierror, W)
        !
        !
        !
        ! DIMENSION OF           AN(N), BN(N), CN(N), AM(M), BM(M), CM(M), Y(IDIMY, N),
        ! ARGUMENTS
        !
        ! LATEST REVISION        JUNE 2004
        !
        ! USAGE                  CALL blktri (IFLG, NP, N, AN, BN, CN, MP, M, AM, BM,
        !                                     CM, IDIMY, Y, ierror, W)
        !
        ! PURPOSE                blktri SOLVES A SYSTEM OF LINEAR EQUATIONS
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
        !                          IN THE PROGRAM CALLING blktri.
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
        !                          calling blktri must be:
        !
        !                               use type_FishpackWorkspace
        !                               TYPE (FishpackWorkspace) :: W
        !
        !                          The first statement makes the fishpack module
        !                          defined in the file "fish.f" available to the
        !                          user program calling blktri.  The second statement
        !                          declares a derived type variable (defined in
        !                          the module "fish.f") which is used internally
        !                          in blktri to dynamically allocate real and complex
        !                          work space used in solution.  An error flag
        !                          (ierror = 20) is set if the required work space
        !                          allocation fails (for example if N, M are too large)
        !                          Real and complex values are set in the components
        !                          of W on a initial (IFLG=0) call to blktri.  These
        !                          must be preserved on non-initial calls (IFLG=1)
        !                          to blktri.  This eliminates redundant calculations
        !                          and saves compute time.
        !               ****       IMPORTANT!  The user program calling blktri should
        !                          include the statement:
        !
        !                               CALL FISHFIN(W)
        !
        !                          after the final approximation is generated by
        !                          blktri.  The will deallocate the real and complex
        !                          work space of W.  Failure to include this statement
        !                          could result in serious memory leakage.
        !
        !
        ! ARGUMENTS
        !
        ! ON OUTPUT              Y
        !                          CONTAINS THE SOLUTION X.
        !
        !                        ierror
        !                          AN ERROR FLAG THAT INDICATES INVALID
        !                          INPUT PARAMETERS.  EXCEPT FOR NUMBER ZER0,
        !                          A SOLUTION IS NOT ATTEMPTED.
        !
        !                        = 0  NO ERROR.
        !                        = 1  M IS LESS THAN 5
        !                        = 2  N IS LESS THAN 5
        !                        = 3  IDIMY IS LESS THAN M.
        !                        = 4  blktri FAILED WHILE COMPUTING RESULTS
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
        !                             be destroyed if blktri is called again with
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
        !                        ierror.
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
        integer (ip),          intent (in)     :: iflg
        integer (ip),          intent (in)     :: np
        integer (ip),          intent (in)     :: n
        integer (ip),          intent (in)     :: mp
        integer (ip),          intent (in)     :: m
        integer (ip),          intent (in)     :: idimy
        integer (ip),          intent (out)    :: ierror
        real (wp), contiguous, intent (in out) :: an(:)
        real (wp), contiguous, intent (in out) :: bn(:)
        real (wp), contiguous, intent (in out) :: cn(:)
        real (wp), contiguous, intent (in out) :: am(:)
        real (wp), contiguous, intent (in out) :: bm(:)
        real (wp), contiguous, intent (in out) :: cm(:)
        real (wp), contiguous, intent (in out) :: y(:,:)
        class (FishpackWorkspace)              :: workspace
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)  :: irwk, icwk
        !--------------------------------------------------------------------------------

        ! Initialize error flag
        ierror = 0

        ! Check if input values are valid - case 1
        if (m < 5) then
            ierror = 1
            return
        end if

        ! Check if input values are valid - case 2
        if (n < 3) then
            ierror = 2
            return
        end if

        ! Check if input values are valid - case 3
        if (idimy < m) then
            ierror = 3
            return
        end if


        if (iflg == 0) then

            ! compute and allocate real and complex required work space
            call workspace%get_block_tridiagonal_workpace_dimensions (n, m, irwk, icwk)

            ! Allocate memory for workspace
            call workspace%create( irwk, icwk, ierror )

        end if

        ! Solve system
        associate( &
            rew => workspace%real_workspace, &
            cxw => workspace%complex_workspace &
            )
            call blktrii( iflg, np, n, an, bn, cn, mp, m, am, bm, cm, idimy, y, &
                ierror, rew, cxw )
        end associate

    end subroutine blktri


    subroutine blktrii( iflg, np, n, an, bn, cn, mp, m, am, bm, cm, &
        idimy, y, ierror, w, wc)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: iflg
        integer (ip), intent (in)     :: np
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: mp
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: idimy
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in out) :: an(*)
        real (wp),    intent (in out) :: bn(*)
        real (wp),    intent (in out) :: cn(*)
        real (wp),    intent (in out) :: am(*)
        real (wp),    intent (in out) :: bm(*)
        real (wp),    intent (in out) :: cm(*)
        real (wp),    intent (in out) :: y(idimy,*)
        real (wp),    intent (in out) :: w(*)
        complex (wp), intent (in out) :: wc(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip)       ::  nh, iwah, iwbh
        integer (ip), save :: iw1, iw2, iw3, iww, iwu, iwd, nl
        !----------------------------------------------

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
            call blktr1(nl, an, bn, cn, m, am, bm, cm, idimy, y, w, wc, &
                w(iw1), w(iw2), w(iw3), w(iwd), w(iww), w(iwu), wc(iw1), &
                wc(iw2), wc(iw3), prod, cprod)
        else
            call blktr1(nl, an, bn, cn, m, am, bm, cm, idimy, y, w, wc, &
                w(iw1), w(iw2), w(iw3), w(iwd), w(iww), w(iwu), wc(iw1), &
                wc(iw2), wc(iw3), prodp, cprodp)
        end if

    end subroutine blktrii


    subroutine blktr1(n, an, bn, cn, m, am, bm, cm, idimy, y, b, bc, &
        w1, w2, w3, wd, ww, wu, cw1, cw2, cw3, prdct, cprdct)
        !
        ! Purpose:
        !
        ! blktr1 solves the linear system
        !
        ! b  contains the roots of all the b polynomials
        ! w1, w2, w3, wd, ww, wu  are all working arrays
        ! prdct  is either prodp or prod depending on whether the boundary
        ! conditions in the m direction are periodic or not
        ! cprdct is either cprodp or cprod which are the complex versions
        ! of prodp and prod. these are called in the event that some
        ! of the roots of the b sub p polynomial are complex
        !
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: idimy
        real (wp),    intent (in)     :: an(*)
        real (wp),    intent (in)     :: bn(*)
        real (wp),    intent (in)     :: cn(*)
        real (wp),    intent (in)     :: am(*)
        real (wp),    intent (in)     :: bm(*)
        real (wp),    intent (in)     :: cm(*)
        real (wp),    intent (in out) :: y(idimy,1)
        real (wp),    intent (in)     :: b(*)
        real (wp),    intent (in out) :: w1(*)
        real (wp),    intent (in out) :: w2(*)
        real (wp),    intent (in out) :: w3(*)
        real (wp),    intent (in)     :: wd(*)
        real (wp),    intent (in)     :: ww(*)
        real (wp),    intent (in)     :: wu(*)
        complex (wp), intent (in)     :: bc(*)
        complex (wp), intent (in)     :: cw1(*)
        complex (wp), intent (in)     :: cw2(*)
        complex (wp), intent (in)     :: cw3(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: kdo, l, ir, i2, i1, i3, i4, irm1, im2, nm2, im3, nm3
        integer (ip) :: im1, nm1, i0, if_rename, i, ipi1, ipi2, ipi3, idxc, nc, idxa, na, ip2
        integer (ip) :: np2, ip1, np1, ip3, np3, iz, nz, izr, ll, ifd, ip_rename, np
        integer (ip) :: imi1, imi2
        real (wp) :: dum
        !-----------------------------------------------

        !
        ! begin reduction phase
        !
        kdo = k - 1
        do l = 1, kdo
            ir = l - 1
            i2 = 2**ir
            i1 = i2/2
            i3 = i2 + i1
            i4 = i2 + i2
            irm1 = ir - 1
            call indxb(i2, ir, im2, nm2)
            call indxb(i1, irm1, im3, nm3)
            call indxb(i3, irm1, im1, nm1)
            i0 = 0
            call prdct(nm2, b(im2), nm3, b(im3), nm1, b(im1), i0, dum, &
                y(1, i2), w3, m, am, bm, cm, wd, ww, wu)
            if_rename = 2**k
            do i = i4, if_rename, i4
                if (i - nm > 0) cycle
                ipi1 = i + i1
                ipi2 = i + i2
                ipi3 = i + i3
                call indxc (i, ir, idxc, nc)
                if (i - if_rename >= 0) cycle
                call indxa (i, ir, idxa, na)
                call indxb(i - i1, irm1, im1, nm1)
                call indxb(ipi2, ir, ip2, np2)
                call indxb(ipi1, irm1, ip1, np1)
                call indxb(ipi3, irm1, ip3, np3)
                call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), w3, &
                    w1, m, am, bm, cm, wd, ww, wu)
                if (ipi2 - nm > 0) then
                    w3(:m) = 0.
                    w2(:m) = 0.
                else
                    call prdct(np2, b(ip2), np1, b(ip1), np3, b(ip3), 0, dum &
                        , y(1, ipi2), w3, m, am, bm, cm, wd, ww, wu)
                    call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), w3 &
                        , w2, m, am, bm, cm, wd, ww, wu)
                end if
                y(:m, i) = w1(:m) + w2(:m) + y(:m, i)
            end do
        end do
        if (npp == 0) then
            if_rename = 2**k
            i = if_rename/2
            i1 = i/2
            call indxb(i - i1, k - 2, im1, nm1)
            call indxb(i + i1, k - 2, ip1, np1)
            call indxb(i, k - 1, iz, nz)
            call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dum, y(1, i) &
                , w1, m, am, bm, cm, wd, ww, wu)
            izr = i
            w2(:m) = w1(:m)
            do ll = 2, k
                l = k - ll + 1
                ir = l - 1
                i2 = 2**ir
                i1 = i2/2
                i = i2
                call indxc (i, ir, idxc, nc)
                call indxb(i, ir, iz, nz)
                call indxb(i - i1, ir - 1, im1, nm1)
                call indxb(i + i1, ir - 1, ip1, np1)
                call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), w1, &
                    w1, m, am, bm, cm, wd, ww, wu)
                w1(:m) = y(:m, i) + w1(:m)
                call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dum, w1 &
                    , w1, m, am, bm, cm, wd, ww, wu)
            end do
            l118: do ll = 2, k
                l = k - ll + 1
                ir = l - 1
                i2 = 2**ir
                i1 = i2/2
                i4 = i2 + i2
                ifd = if_rename - i2
                do i = i2, ifd, i4
                    if (i - i2 - izr /= 0) cycle
                    if (i - nm > 0) cycle  l118
                    call indxa (i, ir, idxa, na)
                    call indxb(i, ir, iz, nz)
                    call indxb(i - i1, ir - 1, im1, nm1)
                    call indxb(i + i1, ir - 1, ip1, np1)
                    call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), w2 &
                        , w2, m, am, bm, cm, wd, ww, wu)
                    w2(:m) = y(:m, i) + w2(:m)
                    call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dum, &
                        w2, w2, m, am, bm, cm, wd, ww, wu)
                    izr = i
                    if (i - nm == 0) exit  l118
                end do
            end do l118
119     continue
        y(:m, nm+1) = y(:m, nm+1) - cn(nm+1)*w1(:m) - an(nm+1)*w2(:m)
        call indxb(if_rename/2, k - 1, im1, nm1)
        call indxb(if_rename, k - 1, ip_rename, np)
        if (ncmplx /= 0) then
            call cprdct(nm + 1, bc(ip_rename), nm1, b(im1), 0, dum, 0, dum, y( &
                1, nm+1), y(1, nm+1), m, am, bm, cm, cw1, cw2, cw3)
        else
            call prdct(nm + 1, b(ip_rename), nm1, b(im1), 0, dum, 0, dum, y(1, &
                nm+1), y(1, nm+1), m, am, bm, cm, wd, ww, wu)
        end if
        w1(:m) = an(1)*y(:m, nm+1)
        w2(:m) = cn(nm)*y(:m, nm+1)
        y(:m, 1) = y(:m, 1) - w1(:m)
        y(:m, nm) = y(:m, nm) - w2(:m)
        do l = 1, kdo
            ir = l - 1
            i2 = 2**ir
            i4 = i2 + i2
            i1 = i2/2
            i = i4
            call indxa (i, ir, idxa, na)
            call indxb(i - i2, ir, im2, nm2)
            call indxb(i - i2 - i1, ir - 1, im3, nm3)
            call indxb(i - i1, ir - 1, im1, nm1)
            call prdct(nm2, b(im2), nm3, b(im3), nm1, b(im1), 0, dum, &
                w1, w1, m, am, bm, cm, wd, ww, wu)
            call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), w1, &
                w1, m, am, bm, cm, wd, ww, wu)
            y(:m, i) = y(:m, i) - w1(:m)
        end do
        !
        izr = nm
        l131: do l = 1, kdo
            ir = l - 1
            i2 = 2**ir
            i1 = i2/2
            i3 = i2 + i1
            i4 = i2 + i2
            irm1 = ir - 1
            do i = i4, if_rename, i4
                ipi1 = i + i1
                ipi2 = i + i2
                ipi3 = i + i3
                if (ipi2 - izr /= 0) then
                    if (i - izr /= 0) cycle
                    cycle  l131
                end if
                call indxc (i, ir, idxc, nc)
                call indxb(ipi2, ir, ip2, np2)
                call indxb(ipi1, irm1, ip1, np1)
                call indxb(ipi3, irm1, ip3, np3)
                call prdct(np2, b(ip2), np1, b(ip1), np3, b(ip3), 0, dum &
                    , w2, w2, m, am, bm, cm, wd, ww, wu)
                call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), w2 &
                    , w2, m, am, bm, cm, wd, ww, wu)
                y(:m, i) = y(:m, i) - w2(:m)
                izr = i
                cycle  l131
            end do
        end do l131
    end if
    !
    ! begin back substitution phase
    !
    do ll = 1, k
        l = k - ll + 1
        ir = l - 1
        irm1 = ir - 1
        i2 = 2**ir
        i1 = i2/2
        i4 = i2 + i2
        ifd = if_rename - i2
        do i = i2, ifd, i4
            if (i - nm > 0) cycle
            imi1 = i - i1
            imi2 = i - i2
            ipi1 = i + i1
            ipi2 = i + i2
            call indxa (i, ir, idxa, na)
            call indxc (i, ir, idxc, nc)
            call indxb(i, ir, iz, nz)
            call indxb(imi1, irm1, im1, nm1)
            call indxb(ipi1, irm1, ip1, np1)
            if (i - i2 <= 0) then
                w1(:m) = 0.
            else
                call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), y( &
                    1, imi2), w1, m, am, bm, cm, wd, ww, wu)
            end if
            if (ipi2 - nm > 0) then
                w2(:m) = 0.
            else
                call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), y( &
                    1, ipi2), w2, m, am, bm, cm, wd, ww, wu)
            end if
            w1(:m) = y(:m, i) + w1(:m) + w2(:m)
            call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dum, w1 &
                , y(1, i), m, am, bm, cm, wd, ww, wu)
        end do
    end do

end subroutine blktr1


function bsrh(xll, xrr, iz, c, a, bh, f, sgn) result(return_value)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    real (wp), intent (in)     :: xll
    real (wp), intent (in)     :: xrr
    integer (ip), intent (in)  :: iz
    real (wp), intent (in)     :: c(*)
    real (wp), intent (in)     :: a(*)
    real (wp), intent (in)     :: bh(*)
    procedure (comf_interface) :: f
    real (wp), intent (in)     :: sgn
    real (wp)                  :: return_value
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    real (wp) :: r1, xl, xr, dx, x
    !-----------------------------------------------

    xl = xll
    xr = xrr
    dx = 0.5_wp * abs(xr - xl)
    x = 0.5_wp * (xl + xr)
    r1 = sgn * f(x, iz, c, a, bh)
    if (r1 >= 0.0_wp) then
        if (r1 == 0.) go to 105
        xr = x
    else
        xl = x
    end if
    dx = 0.5_wp * dx
    do while(dx - cnv > 0.0_wp)
        x = 0.5_wp * (xl + xr)
        r1 = sgn * f(x, iz, c, a, bh)
        if (r1 >= 0.0_wp) then
            if (r1 == 0.0_wp) go to 105
            xr = x
        else
            xl = x
        end if
        dx = 0.5_wp * dx
    end do
105 continue

    return_value = 0.5_wp * (xl + xr)

end function bsrh


subroutine compb(n, ierror, an, bn, cn, b, bc, ah, bh)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: n
    integer (ip), intent (out)    :: ierror
    real (wp),    intent (in)     :: an(*)
    real (wp),    intent (in)     :: bn(*)
    real (wp),    intent (in)     :: cn(*)
    real (wp),    intent (in out) :: b(*)
    real (wp),    intent (in out) :: ah(*)
    real (wp),    intent (in out) :: bh(*)
    complex (wp), intent (in out) :: bc(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: j, if, kdo, l, ir, i2, i4, ipl, ifd, i, ib, nb, js, jf
    integer (ip) ::  ls, lh, nmp, l1, l2, j2, j1, n2m2
    real (wp)    :: dum, bnorm, arg, d1, d2, d3
    !-----------------------------------------------
    !
    !     compb computes the roots of the b polynomials using subroutine
    !     tevls which is a modification the eispack program tqlrat.
    !     ierror is set to 4 if either tevls fails or if a(j+1)*c(j) is
    !     less than zero for some j.  ah, bh are temporary work arrays.
    !
    eps = epsilon(dum)
    bnorm = abs(bn(1))
    do j = 2, nm
        bnorm = max(bnorm, abs(bn(j)))
        arg = an(j)*cn(j-1)
        if (arg < 0.) go to 119
        b(j) = sign(sqrt(arg), an(j))
    end do
    cnv = eps*bnorm
    if = 2**k
    kdo = k - 1
    l108: do l = 1, kdo
        ir = l - 1
        i2 = 2**ir
        i4 = i2 + i2
        ipl = i4 - 1
        ifd = if - i4
        do i = i4, ifd, i4
            call indxb(i, l, ib, nb)
            if (nb <= 0) cycle  l108
            js = i - ipl
            jf = js + nb - 1
            ls = 0
            bh(:jf-js+1) = bn(js:jf)
            ah(:jf-js+1) = b(js:jf)
            call tevls (nb, bh, ah, ierror)
            if (ierror /= 0) go to 118
            lh = ib - 1
            if (nb > 0) then
                b(lh+1:nb+lh) = -bh(:nb)
                lh = nb + lh
            end if
        end do
    end do l108
    b(:nm) = -bn(:nm)
    if (npp == 0) then
        nmp = nm + 1
        nb = nm + nmp
        do j = 1, nb
            l1 = mod(j - 1, nmp) + 1
            l2 = mod(j + nm - 1, nmp) + 1
            arg = an(l1)*cn(l2)
            if (arg < 0.) go to 119
            bh(j) = sign(sqrt(arg), (-an(l1)))
            ah(j) = -bn(l1)
        end do
        call tevls (nb, ah, bh, ierror)
        if (ierror /= 0) go to 118
        call indxb(if, k - 1, j2, lh)
        call indxb(if/2, k - 1, j1, lh)
        j2 = j2 + 1
        lh = j2
        n2m2 = j2 + nm + nm - 2
114 continue
    d1 = abs(b(j1)-b(j2-1))
    d2 = abs(b(j1)-b(j2))
    d3 = abs(b(j1)-b(j2+1))
    if (d2>=d1 .or. d2>=d3) then
        b(lh) = b(j2)
        j2 = j2 + 1
        lh = lh + 1
        if (j2 - n2m2 <= 0) go to 114
    else
        j2 = j2 + 1
        j1 = j1 + 1
        if (j2 - n2m2 <= 0) go to 114
    end if
    b(lh) = b(n2m2+1)
    call indxb(if, k - 1, j1, j2)
    j2 = j1 + nmp + nmp
    call ppadd(nm + 1, ierror, an, cn, bc(j1), b(j1), b(j2))
end if
return 
118 continue
    ierror = 4
    return
119 continue
    ierror = 5

end subroutine compb


subroutine cprod(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, yy, m, a, b, c, d, w, y)
    !
    ! Purpose:
    !
    ! cprod applies a sequence of matrix operations to the vector x and
    ! stores the result in yy (complex case)
    !
    ! aa   array containing scalar multipliers of the vector x
    !
    ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
    !
    ! bd, bm1, bm2 are arrays containing roots of certian b polynomials
    !
    ! na is the length of the array aa
    !
    ! x, yy the matrix operations are applied to x and the result is yy
    !
    ! a, b, c  are arrays which contain the tridiagonal matrix
    !
    ! m  is the order of the matrix
    !
    ! d, w, y are working arrays
    !
    ! isgn  determines whether or not a change in sign is made
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: nd
    integer (ip), intent (in)     :: nm1
    integer (ip), intent (in)     :: nm2
    integer (ip), intent (in)     :: na
    integer (ip), intent (in)     :: m
    real (wp),    intent (in)     :: bm1(*)
    real (wp),    intent (in)     :: bm2(*)
    real (wp),    intent (in)     :: aa(*)
    real (wp),    intent (in)     :: x(*)
    real (wp),    intent (out)    :: yy(*)
    real (wp),    intent (in)     :: a(*)
    real (wp),    intent (in)     :: b(*)
    real (wp),    intent (in)     :: c(*)
    complex (wp), intent (in)     :: bd(*)
    complex (wp), intent (in out) :: d(*)
    complex (wp), intent (in out) :: w(*)
    complex (wp), intent (in out) :: y(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: j, mm, id, m1, m2, ia, iflg, k
    real (wp)    :: rt
    complex (wp) :: crt, den, y1, y2
    !-----------------------------------------------

    do j = 1, m
        y(j) = cmplx(x(j), 0.)
    end do
    mm = m - 1
    id = nd
    m1 = nm1
    m2 = nm2
    ia = na
102 continue
    iflg = 0
    if (id > 0) then
        crt = bd(id)
        id = id - 1
        !
        ! begin solution to system
        !
        d(m) = a(m)/(b(m)-crt)
        w(m) = y(m)/(b(m)-crt)
        do j = 2, mm
            k = m - j
            den = b(k+1) - crt - c(k+1)*d(k+2)
            d(k+1) = a(k+1)/den
            w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
        end do
        den = b(1) - crt - c(1)*d(2)
        if (abs(den) /= 0.) then
            y(1) = (y(1)-c(1)*w(2))/den
        else
            y(1) = (1., 0.)
        end if
        do j = 2, m
            y(j) = w(j) - d(j)*y(j-1)
        end do
    end if
    if (m1 <= 0) then
        if (m2 <= 0) go to 121
        rt = bm2(m2)
        m2 = m2 - 1
    else
        if (m2 <= 0) then
            rt = bm1(m1)
            m1 = m1 - 1
        else
            if (abs(bm1(m1)) - abs(bm2(m2)) > 0.) then
                rt = bm1(m1)
                m1 = m1 - 1
            else
                rt = bm2(m2)
                m2 = m2 - 1
            end if
        end if
    end if
    y1 = (b(1)-rt)*y(1) + c(1)*y(2)
    if (mm - 2 >= 0) then
        do j = 2, mm
            y2 = a(j)*y(j-1) + (b(j)-rt)*y(j) + c(j)*y(j+1)
            y(j-1) = y1
            y1 = y2
        end do
    end if
    y(m) = a(m)*y(m-1) + (b(m)-rt)*y(m)
    y(m-1) = y1
    iflg = 1
    go to 102
121 continue
    if (ia > 0) then
        rt = aa(ia)
        ia = ia - 1
        iflg = 1
        !
        ! scalar multiplication
        !
        y(:m) = rt*y(:m)
    end if
    if (iflg > 0) go to 102
    do j = 1, m
        yy(j) = real(y(j))
    end do

end subroutine cprod


subroutine cprodp(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, yy, m, a, &
    b, c, d, u, y)
        !
    ! Purpose:
    !
    ! cprodp applies a sequence of matrix operations to the vector x and
    ! stores the result in yy       periodic boundary conditions
    ! and  complex  case
    !
    ! bd, bm1, bm2 are arrays containing roots of certian b polynomials
    ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
    ! aa   array containing scalar multipliers of the vector x
    ! na is the length of the array aa
    ! x, yy the matrix operations are applied to x and the result is yy
    ! a, b, c  are arrays which contain the tridiagonal matrix
    ! m  is the order of the matrix
    ! d, u, y are working arrays
    ! isgn  determines whether or not a change in sign is made
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: nd
    integer (ip), intent (in)     :: nm1
    integer (ip), intent (in)     :: nm2
    integer (ip), intent (in)     :: na
    integer (ip), intent (in)     :: m
    real (wp),    intent (in)     :: bm1(*)
    real (wp),    intent (in)     :: bm2(*)
    real (wp),    intent (in)     :: aa(*)
    real (wp),    intent (in)     :: x(*)
    real (wp),    intent (out)    :: yy(*)
    real (wp),    intent (in)     :: a(*)
    real (wp),    intent (in)     :: b(*)
    real (wp),    intent (in)     :: c(*)
    complex (wp), intent (in)     :: bd(*)
    complex (wp), intent (in out) :: d(*)
    complex (wp), intent (in out) :: u(*)
    complex (wp), intent (in out) :: y(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: j, mm, mm2, id, m1, m2, ia, iflg, k
    real (wp)    :: rt
    complex (wp) :: v, den, bh, ym, am, y1, y2, yh, crt
    !-----------------------------------------------

    do j = 1, m
        y(j) = cmplx(x(j), 0.0_wp, kind=wp)
    end do

    mm = m - 1
    mm2 = m - 2
    id = nd
    m1 = nm1
    m2 = nm2
    ia = na
102 continue
    iflg = 0
    if (id > 0) then
        crt = bd(id)
        id = id - 1
        iflg = 1
        !
        ! begin solution to system
        !
        bh = b(m) - crt
        ym = y(m)
        den = b(1) - crt
        d(1) = c(1)/den
        u(1) = a(1)/den
        y(1) = y(1)/den
        v = cmplx(c(m), 0.0_wp, kind=wp)
        if (mm2 - 2 >= 0) then
            do j = 2, mm2
                den = b(j) - crt - a(j)*d(j-1)
                d(j) = c(j)/den
                u(j) = -a(j)*u(j-1)/den
                y(j) = (y(j)-a(j)*y(j-1))/den
                bh = bh - v*u(j-1)
                ym = ym - v*y(j-1)
                v = -v*d(j-1)
            end do
        end if
        den = b(m-1) - crt - a(m-1)*d(m-2)
        d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
        y(m-1) = (y(m-1)-a(m-1)*y(m-2))/den
        am = a(m) - v*d(m-2)
        bh = bh - v*u(m-2)
        ym = ym - v*y(m-2)
        den = bh - am*d(m-1)
        if (abs(den) /= 0.) then
            y(m) = (ym - am*y(m-1))/den
        else
            y(m) = (1., 0.)
        end if
        y(m-1) = y(m-1) - d(m-1)*y(m)
        do j = 2, mm
            k = m - j
            y(k) = y(k) - d(k)*y(k+1) - u(k)*y(m)
        end do
    end if
    if (m1 <= 0) then
        if (m2 <= 0) go to 123
        rt = bm2(m2)
        m2 = m2 - 1
    else
        if (m2 <= 0) then
            rt = bm1(m1)
            m1 = m1 - 1
        else
            if (abs(bm1(m1)) - abs(bm2(m2)) > 0.) then
                rt = bm1(m1)
                m1 = m1 - 1
            else
                rt = bm2(m2)
                m2 = m2 - 1
            !
            ! matrix multiplication
            !
            end if
        end if
    end if
    yh = y(1)
    y1 = (b(1)-rt)*y(1) + c(1)*y(2) + a(1)*y(m)
    if (mm - 2 >= 0) then
        do j = 2, mm
            y2 = a(j)*y(j-1) + (b(j)-rt)*y(j) + c(j)*y(j+1)
            y(j-1) = y1
            y1 = y2
        end do
    end if
    y(m) = a(m)*y(m-1) + (b(m)-rt)*y(m) + c(m)*yh
    y(m-1) = y1
    iflg = 1
    go to 102
123 continue
    if (ia > 0) then
        rt = aa(ia)
        ia = ia - 1
        iflg = 1
        !
        ! scalar multiplication
        !
        y(:m) = rt*y(:m)
    end if
    if (iflg > 0) go to 102
    do j = 1, m
        yy(j) = real(y(j))
    end do

end subroutine cprodp

pure subroutine indxa(i, ir, idxa, na)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)  :: i
    integer (ip), intent (in)  :: ir
    integer (ip), intent (out) :: idxa
    integer (ip), intent (out) :: na
    !-----------------------------------------------

    na = 2**ir
    idxa = i - na + 1
    if (i - nm > 0) then
        na = 0
    end if

end subroutine indxa

pure subroutine indxb(i, ir, idx, idp)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)  :: i
    integer (ip), intent (in)  :: ir
    integer (ip), intent (out) :: idx
    integer (ip), intent (out) :: idp
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: izh, id, ipl
    !-----------------------------------------------
    !
    ! b(idx) is the location of the first root of the b(i, ir) polynomial
    !
    idp = 0
    if (ir >= 0) then
        if (ir <= 0) then
            if (i - nm > 0) go to 107
            idx = i
            idp = 1
            return
        end if
        izh = 2**ir
        id = i - izh - izh
        idx = id + id + (ir - 1)*ik + ir + (ik - i)/izh + 4
        ipl = izh - 1
        idp = izh + izh - 1
        if (i - ipl - nm > 0) then
            idp = 0
            return
        end if
        if (i + ipl - nm > 0) then
            idp = nm + ipl - i + 1
        end if
    end if
107 continue

end subroutine indxb


pure subroutine indxc(i, ir, idxc, nc)

    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)  :: i
    integer (ip), intent (in)  :: ir
    integer (ip), intent (out) :: idxc
    integer (ip), intent (out) :: nc
    !-----------------------------------------------

    nc = 2**ir
    idxc = i
    if (idxc + nc - 1 - nm > 0) then
        nc = 0
    end if

end subroutine indxc


subroutine ppadd(n, ierror, a, c, cbp, bp, bh)
    !
    ! Purpose
    !
    ! ppadd computes the eigenvalues of the periodic tridiagonal matrix
    ! with coefficients an, bn, cn
    !
    ! n is the order of the bh and bp polynomials
    ! on output bp contains the eigenvalues
    ! cbp is the same as bp except type complex
    ! bh is used to temporarily store the roots of the b hat polynomial
    ! which enters through bp
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: n
    integer (ip), intent (out)    :: ierror
    real (wp),    intent (in)     :: a(*)
    real (wp),    intent (in)     :: c(*)
    real (wp),    intent (in out) :: bp(*)
    real (wp),    intent (in out) :: bh(*)
    complex (wp), intent (in out) :: cbp(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: iz, izm, izm2, j, nt, modiz, is, if, ig, it, icv, i3, i2, nhalf
    real (wp)    :: r4, r5, r6, scnv, xl, db, sgn, xr, xm, psg
    complex (wp) :: cx, fsg, hsg, dd, f, fp, fpp, cdis, r1, r2, r3
    !-----------------------------------------------

    scnv = sqrt(cnv)
    iz = n
    izm = iz - 1
    izm2 = iz - 2
    if (bp(n) - bp(1) <= 0.) then
        if (bp(n) - bp(1) == 0.) go to 142
        bh(:n) = bp(n:1:(-1))
    else
        bh(:n) = bp(:n)
    end if
    ncmplx = 0
    modiz = mod(iz, 2)
    is = 1
    if (modiz /= 0) then
        if (a(1) < 0.) go to 110
        if (a(1) == 0.) go to 142
    end if
    xl = bh(1)
    db = bh(3) - bh(1)
    xl = xl - db
    r4 = psgf(xl, iz, c, a, bh)
    do while(r4 <= 0.)
        xl = xl - db
        r4 = psgf(xl, iz, c, a, bh)
    end do
    sgn = -1.
    cbp(1) = cmplx(bsrh(xl, bh(1), iz, c, a, bh, psgf, sgn), 0.0_wp, kind=wp)
    bp(1) = real(cbp(1))
    is = 2
110 continue
    if = iz - 1
    if (modiz /= 0) then
        if (a(1) > 0.0_wp) go to 115
        if (a(1) == 0.0_wp) go to 142
    end if
    xr = bh(iz)
    db = bh(iz) - bh(iz-2)
    xr = xr + db
    r5 = psgf(xr, iz, c, a, bh)
    do while(r5 < 0.0_wp)
        xr = xr + db
        r5 = psgf(xr, iz, c, a, bh)
    end do
    sgn = 1.
    cbp(iz) = cmplx(bsrh(bh(iz), xr, iz, c, a, bh, psgf, sgn), 0.0_wp, kind=wp)
    if = iz - 2
115 continue
    do ig = is, if, 2
        xl = bh(ig)
        xr = bh(ig+1)
        sgn = -1.
        xm = bsrh(xl, xr, iz, c, a, bh, ppspf, sgn)
        psg = psgf(xm, iz, c, a, bh)
        if (abs(psg) - eps <= 0.) go to 118
        r6 = psg*ppsgf(xm, iz, c, a, bh)
        if (r6 > 0.) go to 119
        if (r6 == 0.) go to 118
        sgn = 1.
        cbp(ig) = cmplx(bsrh(bh(ig), xm, iz, c, a, bh, psgf, sgn), 0.0_wp, kind=wp)
        !        bp(ig) = real(cbp(ig))
        sgn = -1.
        cbp(ig+1) = cmplx(bsrh(xm, bh(ig+1), iz, c, a, bh, psgf, sgn), 0.0_wp, kind=wp)
        !        bp(ig) = real(cbp(ig))
        !        bp(ig+1) = real(cbp(ig+1))
        cycle
    !
    !     case of a multiple zero
    !
118 continue
    cbp(ig) = cmplx(xm, 0.0_wp, kind=wp)
    cbp(ig+1) = cmplx(xm, 0.0_wp, kind=wp)
    !        bp(ig) = real(cbp(ig))
    !        bp(ig+1) = real(cbp(ig+1))
    cycle
!
!     case of a complex zero
!
119 continue
    it = 0
    icv = 0
    cx = cmplx(xm, 0.0_wp, kind=wp)
120 continue
    fsg = (1.0_wp, 0.0_wp)
    hsg = (1.0_wp, 0.0_wp)
    fp = (0.0_wp, 0.0_wp)
    fpp = (0.0_wp, 0.0_wp)
    do j = 1, iz
        dd = 1./(cx - bh(j))
        fsg = fsg*a(j)*dd
        hsg = hsg*c(j)*dd
        fp = fp + dd
        fpp = fpp - dd*dd
    end do
    if (modiz == 0) then
        f = (1.0_wp, 0.0_wp) - fsg - hsg
    else
        f = (1.0_wp, 0.0_wp) + fsg + hsg
    end if
    i3 = 0
    if (abs(fp) > 0.0_wp) then
        i3 = 1
        r3 = -f/fp
    end if
    i2 = 0
    if (abs(fpp) > 0.0_wp) then
        i2 = 1
        cdis = csqrt(fp**2 - 2.*f*fpp)
        r1 = cdis - fp
        r2 = (-fp) - cdis
        if (abs(r1) - abs(r2) > 0.) then
            r1 = r1/fpp
        else
            r1 = r2/fpp
        end if
        r2 = 2.*f/fpp/r1
        if (abs(r2) < abs(r1)) r1 = r2
        if (i3 <= 0) go to 133
        if (abs(r3) < abs(r1)) r1 = r3
        go to 133
    end if
    r1 = r3
133 continue
    cx = cx + r1
    it = it + 1
    if (it > 50) go to 142
    if (abs(r1) > scnv) go to 120
    if (icv > 0) go to 135
    icv = 1
    go to 120
135 continue
    cbp(ig) = cx
    cbp(ig+1) = conjg(cx)
end do
if (abs(cbp(n)) - abs(cbp(1)) <= 0.) then
    if (abs(cbp(n)) - abs(cbp(1)) == 0.) go to 142
    nhalf = n/2
    do j = 1, nhalf
        nt = n - j
        cx = cbp(j)
        cbp(j) = cbp(nt+1)
        cbp(nt+1) = cx
    end do
end if
ncmplx = 1
do j = 2, iz
    if (aimag(cbp(j)) /= 0.) go to 143
end do
ncmplx = 0
bp(1) = real(cbp(1))
do j = 2, iz
    bp(j) = real(cbp(j))
end do
go to 143
142 continue
    ierror = 4
143 continue

end subroutine ppadd


subroutine prod(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, w, u)
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: nd
    integer (ip), intent (in)     :: nm1
    integer (ip), intent (in)     :: nm2
    integer (ip), intent (in)     :: na
    integer (ip), intent (in)     :: m
    real (wp),    intent (in)     :: bd(*)
    real (wp),    intent (in)     :: bm1(*)
    real (wp),    intent (in)     :: bm2(*)
    real (wp),    intent (in)     :: aa(*)
    real (wp),    intent (in)     :: x(*)
    real (wp),    intent (in out) :: y(*)
    real (wp),    intent (in)     :: a(*)
    real (wp),    intent (in)     :: b(*)
    real (wp),    intent (in)     :: c(*)
    real (wp),    intent (in out) :: d(*)
    real (wp),    intent (in out) :: w(*)
    real (wp),    intent (in out) :: u(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: j, mm, id, ibr, m1, m2, ia, k
    real (wp)    :: rt, den
    !-----------------------------------------------
    !
    ! prod applies a sequence of matrix operations to the vector x and
    ! stores the result in y
    ! bd, bm1, bm2 are arrays containing roots of certian b polynomials
    ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
    ! aa   array containing scalar multipliers of the vector x
    ! na is the length of the array aa
    ! x, y  the matrix operations are applied to x and the result is y
    ! a, b, c  are arrays which contain the tridiagonal matrix
    ! m  is the order of the matrix
    ! d, w, u are working arrays
    ! is  determines whether or not a change in sign is made
    !
    w(:m) = x(:m)
    y(:m) = w(:m)
    mm = m - 1
    id = nd
    ibr = 0
    m1 = nm1
    m2 = nm2
    ia = na
102 continue
    if (ia > 0) then
        rt = aa(ia)
        if (nd == 0) rt = -rt
        ia = ia - 1
        !
        ! scalar multiplication
        !
        y(:m) = rt*w(:m)
    end if
    if (id <= 0) go to 125
    rt = bd(id)
    id = id - 1
    if (id == 0) ibr = 1
    !
    ! begin solution to system
    !
    d(m) = a(m)/(b(m)-rt)
    w(m) = y(m)/(b(m)-rt)
    do j = 2, mm
        k = m - j
        den = b(k+1) - rt - c(k+1)*d(k+2)
        d(k+1) = a(k+1)/den
        w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
    end do
    den = b(1) - rt - c(1)*d(2)
    w(1) = 1.
    if (den /= 0.) then
        w(1) = (y(1)-c(1)*w(2))/den
    end if
    do j = 2, m
        w(j) = w(j) - d(j)*w(j-1)
    end do
    if (na > 0) go to 102
    go to 113
111 continue
    y(:m) = w(:m)
    ibr = 1
    go to 102
113 continue
    if (m1 <= 0) then
        if (m2 <= 0) go to 111
    else
        if (m2 > 0) then
            if (abs(bm1(m1)) - abs(bm2(m2)) <= 0.) go to 120
        end if
        if (ibr <= 0) then
            if (abs(bm1(m1)-bd(id)) - abs(bm1(m1)-rt) < 0.) go to 111
        end if
        rt = rt - bm1(m1)
        m1 = m1 - 1
        go to 123
    end if
120 continue
    if (ibr <= 0) then
        if (abs(bm2(m2)-bd(id)) - abs(bm2(m2)-rt) < 0.) go to 111
    end if
    rt = rt - bm2(m2)
    m2 = m2 - 1
123 continue
    y(:m) = y(:m) + rt*w(:m)
    go to 102
125 continue

end subroutine prod


subroutine prodp( nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, u, w)
    !
    ! purpose:
    !
    ! prodp applies a sequence of matrix operations to the vector x and
    ! stores the result in y periodic boundary conditions
    !
    ! bd, bm1, bm2 are arrays containing roots of certain b polynomials
    !
    ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
    !
    ! aa   array containing scalar multipliers of the vector x
    !
    ! na is the length of the array aa
    !
    ! x, y  the matrix operations are applied to x and the result is y
    !
    ! a, b, c  are arrays which contain the tridiagonal matrix
    !
    ! m  is the order of the matrix
    !
    ! d, u, w are working arrays
    !
    ! is  determines whether or not a change in sign is made
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: nd
    integer (ip), intent (in)     :: nm1
    integer (ip), intent (in)     :: nm2
    integer (ip), intent (in)     :: na
    integer (ip), intent (in)     :: m
    real (wp),    intent (in)     :: bd(*)
    real (wp),    intent (in)     :: bm1(*)
    real (wp),    intent (in)     :: bm2(*)
    real (wp),    intent (in)     :: aa(*)
    real (wp),    intent (in)     :: x(*)
    real (wp),    intent (in out) :: y(*)
    real (wp),    intent (in)     :: a(*)
    real (wp),    intent (in)     :: b(*)
    real (wp),    intent (in)     :: c(*)
    real (wp),    intent (in out) :: d(*)
    real (wp),    intent (in out) :: u(*)
    real (wp),    intent (in out) :: w(*)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: j, mm, mm2, id, ibr, m1, m2, ia, k
    real (wp)    :: rt, bh, ym, den, v, am
    !-----------------------------------------------

    y(:m) = x(:m)
    w(:m) = y(:m)
    mm = m - 1
    mm2 = m - 2
    id = nd
    ibr = 0
    m1 = nm1
    m2 = nm2
    ia = na
102 continue
    if (ia > 0) then
        rt = aa(ia)
        if (nd == 0) rt = -rt
        ia = ia - 1
        y(:m) = rt*w(:m)
    end if
    if (id <= 0) go to 128
    rt = bd(id)
    id = id - 1
    if (id == 0) ibr = 1
    !
    ! begin solution to system
    !
    bh = b(m) - rt
    ym = y(m)
    den = b(1) - rt
    d(1) = c(1)/den
    u(1) = a(1)/den
    w(1) = y(1)/den
    v = c(m)
    if (mm2 - 2 >= 0) then
        do j = 2, mm2
            den = b(j) - rt - a(j)*d(j-1)
            d(j) = c(j)/den
            u(j) = -a(j)*u(j-1)/den
            w(j) = (y(j)-a(j)*w(j-1))/den
            bh = bh - v*u(j-1)
            ym = ym - v*w(j-1)
            v = -v*d(j-1)
        end do
    end if
    den = b(m-1) - rt - a(m-1)*d(m-2)
    d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
    w(m-1) = (y(m-1)-a(m-1)*w(m-2))/den
    am = a(m) - v*d(m-2)
    bh = bh - v*u(m-2)
    ym = ym - v*w(m-2)
    den = bh - am*d(m-1)
    if (den /= 0.) then
        w(m) = (ym - am*w(m-1))/den
    else
        w(m) = 1.
    end if
    w(m-1) = w(m-1) - d(m-1)*w(m)
    do j = 2, mm
        k = m - j
        w(k) = w(k) - d(k)*w(k+1) - u(k)*w(m)
    end do
    if (na > 0) go to 102
    go to 116
114 continue
    y(:m) = w(:m)
    ibr = 1
    go to 102
116 continue
    if (m1 <= 0) then
        if (m2 <= 0) go to 114
    else
        if (m2 > 0) then
            if (abs(bm1(m1)) - abs(bm2(m2)) <= 0.) go to 123
        end if
        if (ibr <= 0) then
            if (abs(bm1(m1)-bd(id)) - abs(bm1(m1)-rt) < 0.) go to 114
        end if
        rt = rt - bm1(m1)
        m1 = m1 - 1
        go to 126
    end if
123 continue
    if (ibr <= 0) then
        if (abs(bm2(m2)-bd(id)) - abs(bm2(m2)-rt) < 0.) go to 114
    end if
    rt = rt - bm2(m2)
    m2 = m2 - 1
126 continue
    y(:m) = y(:m) + rt*w(:m)
    go to 102
128 continue

end subroutine prodp


subroutine tevls(n, d, e2, ierr)
    !
    !
    !     real sqrt, abs, sign
    !
    !
    !     this subroutine is a modification of the eispack subroutine tqlrat
    !     algorithm 464, comm. acm 16, 689(1973) by reinsch.
    !
    !     this subroutine finds the eigenvalues of a symmetric
    !     tridiagonal matrix by the rational ql method.
    !
    !     on input-
    !
    !        n is the order of the matrix,
    !
    !        d contains the diagonal elements of the input matrix,
    !
    !        e2 contains the                subdiagonal elements of the
    !          input matrix in its last n-1 positions.  e2(1) is arbitrary.
    !
    !      on output-
    !
    !        d contains the eigenvalues in ascending order.  if an
    !          error exit is made, the eigenvalues are correct and
    !          ordered for indices 1, 2, ...ierr-1, but may not be
    !          the smallest eigenvalues,
    !
    !        e2 has been destroyed,
    !
    !        ierr is set to
    !          zero       for normal return,
    !          j          if the j-th eigenvalue has not been
    !                     determined after 30 iterations.
    !
    !     questions and comments should be directed to b. s. garbow,
    !     applied mathematics division, argonne national laboratory
    !
    !
    !     ********** eps is a machine dependent parameter specifying
    !                the relative precision of floating point arithmetic.
    !
    !                **********
    !
    !-----------------------------------------------
    ! Dictionary: calling arguments
    !-----------------------------------------------
    integer (ip), intent (in)     :: n
    integer (ip), intent (out)    :: ierr
    real (wp),    intent (in out) :: d(n)
    real (wp),    intent (in out) :: e2(n)
    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    integer (ip) :: i, j, l, m, ii, l1, mml, nhalf, ntop
    real (wp) :: b, c, f, g, h, p, r, s, dhold
    !-----------------------------------------------

    ierr = 0
    if (n /= 1) then
        !
        e2(:n-1) = e2(2:n)*e2(2:n)
        !
        f = 0.0
        b = 0.0
        e2(n) = 0.0
        !
        do l = 1, n
            j = 0
            h = eps*(abs(d(l))+sqrt(e2(l)))
            if (b <= h) then
                b = h
                c = b*b
            end if
            !
            !     ********** look for small squared sub-diagonal element **********
            !
            do m = l, n
                if (e2(m) > c) cycle
                exit
            !
            !     ********** e2(n) is always zero, so there is no exit
            !                through the bottom of the loop **********
            !
            end do
            !
            if (m /= l) then
105         continue
            if (j == 30) go to 114
            j = j + 1
            !
            !     ********** form shift **********
            !
            l1 = l + 1
            s = sqrt(e2(l))
            g = d(l)
            p = (d(l1)-g)/(2.0*s)
            r = sqrt(p*p + 1.0)
            d(l) = s/(p + sign(r, p))
            h = g - d(l)
            !
            d(l1:n) = d(l1:n) - h
            !
            f = f + h
            !
            !     ********** rational ql transformation **********
            !
            g = d(m)
            if (g == 0.0) g = b
            h = g
            s = 0.0
            mml = m - l
            !
            !     ********** for i=m-1 step -1 until l do -- **********
            !
            do ii = 1, mml
                i = m - ii
                p = g*h
                r = p + e2(i)
                e2(i+1) = s*r
                s = e2(i)/r
                d(i+1) = h + s*(h + d(i))
                g = d(i) - e2(i)/g
                if (g == 0.0) g = b
                h = g*p/r
            end do
            !
            e2(l) = s*g
            d(l) = h
            !
            !     ********** guard against underflowed h **********
            !
            if (h == 0.0) go to 108
            if (abs(e2(l)) <= abs(c/h)) go to 108
            e2(l) = h*e2(l)
            if (e2(l) /= 0.0) go to 105
        end if
108 continue
    p = d(l) + f
    !
    !     ********** order eigenvalues **********
    !
    if (l /= 1) then
        !
        !     ********** for i=l step -1 until 2 do -- **********
        !
        do ii = 2, l
            i = l + 2 - ii
            if (p >= d(i-1)) go to 111
            d(i) = d(i-1)
        end do
    end if
    !
    i = 1
111 continue
    d(i) = p
end do
!
if (abs(d(n)) >= abs(d(1))) go to 115
nhalf = n/2
do i = 1, nhalf
    ntop = n - i
    dhold = d(i)
    d(i) = d(ntop+1)
    d(ntop+1) = dhold
end do
go to 115
!
!     ********** set error -- no convergence to an
!                eigenvalue after 30 iterations **********
!
114 continue
    ierr = l
end if
115 continue
    return
!
!     ********** last card of tqlrat **********
!
end subroutine tevls


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
