!
!     file blktri.f90
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  Version 1.1                    *
!     *                                                               *
!     *                      A Package of Fortran                     *
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
! SUBROUTINE blktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm,
!            idimy, y, ierror, workspace)
!
!
!
! DIMENSION OF           an(n), bn(n), cn(n), am(m), bm(m), cm(m), y(idimy, n)
! ARGUMENTS
!
! LATEST REVISION        April 2016
!
! USAGE                  call blktri(iflg, np, n, an, bn, cn, mp, m, &
!                                   am, bm, cm, idimy, y, ierror, workspace)
!
! PURPOSE                blktri solves a system of linear equations
!                        of the form
!
!                        an(j)*x(i, j-1) + am(i)*x(i-1, j) +
!                        (bn(j)+bm(i))*x(i, j) + cn(j)*x(i, j+1) +
!                        cm(i)*x(i+1, j) = y(i, j)
!
!                        for i = 1, 2, ..., m  and  j = 1, 2, ..., n.
!
!                        i+1 and i-1 are evaluated modulo m and
!                        j+1 and j-1 modulo n, i.e.,
!
!                        x(i, 0) = x(i, n),  x(i, n+1) = x(i, 1),
!                        x(0, j) = x(m, j),  x(m+1, j) = x(1, j).
!
!                        These equations usually result from the
!                        discretization of separable elliptic
!                        equations. Boundary conditions may be
!                        dirichlet, neumann, or periodic.
!
! ARGUMENTS
!
! ON INPUT               iflg
!
!                          = 0  Unitialization only.
!                               certain quantities that depend on np,
!                               n, an, bn, and cn are computed and
!                               stored in derived data type w (see
!                               description of w below)
!
!                          = 1  The quantities that were computed
!                               in the initialization are used
!                               to obtain the solution x(i, j).
!
!                               note:
!                               A call with iflg=0 takes
!                               approximately one half the time
!                               as a call with iflg = 1.
!                               however, the initialization does
!                               not have to be repeated unless np,
!                               n, an, bn, or cn change.
!
!                        np
!                          = 0  If an(1) and cn(n) are not zero,
!                               which corresponds to periodic
!                               bounary conditions.
!
!                          = 1  If an(1) and cn(n) are zero.
!
!                        n
!                          The number of unknowns in the j-direction.
!                          n must be greater than 4.
!                          The operation count is proportional to
!                          m*n*log2(n), hence n should be selected
!                          less than or equal to m.
!
!                        an, bn, cn
!                          One-dimensional arrays of length n
!                          that specify the coefficients in the
!                          linear equations given above.
!
!                        mp
!                          = 0  If am(1) and cm(m) are not zero,
!                               which corresponds to periodic
!                               boundary conditions.
!
!                          = 1  If am(1) = cm(m) = 0  .
!
!                        m
!                          The number of unknowns in the i-direction.
!                           m must be greater than 4.
!
!                        am, bm, cm
!                          One-dimensional arrays of length m that
!                          specify the coefficients in the linear
!                          equations given above.
!
!                        idimy
!                          The row (or first) dimension of the
!                          two-dimensional array y as it appears
!                          in the program calling blktri.
!                          This parameter is used to specify the
!                          variable dimension of y.
!                          idimy must be at least m.
!
!                        y
!                          A two-dimensional array that specifies
!                          the values of the right side of the linear
!                          system of equations given above.
!                          y must be dimensioned at least m*n.
!
!                        workspace
!                          An object of class(FishpackWorkspace)
!                          that must be declared by the user.  The first
!                          two declarative statements in the user program
!                          calling blktri must be:
!
!                               use type_FishpackWorkspace
!                               type(Fishpackworkspace) :: workspace
!
!                          The first statement makes the fishpack module
!                          defined in the file "type_FishpackWorkspace.f90"
!                          available to the user program calling blktri.
!                          The second statement declares a derived type variable
!                          (defined in the module "type_FishpackWorkspace.f90")
!                          which is used internally in blktri to dynamically
!                          allocate real and complex workspace used in solution.
!                          An error flag (ierror = 20) is set if the required
!                          workspace allocation fails (for example if n, m
!                          are too large). Real and complex values are set in
!                          the components of workspace on a initial (iflg=0)
!                          call to blktri.  These must be preserved on
!                          non-initial calls (iflg=1) to blktri.
!                          This eliminates redundant calculations
!                          and saves compute time.
!
!               ****       IMPORTANT!  The user program calling blktri should
!                          include the statement:
!
!                              call workspace%destroy()
!
!                          after the final approximation is generated by
!                          blktri. This will deallocate the real and complex
!                          array components of workspace. Failure to include this
!                          statement could result in serious memory leakage.
!
!
! ARGUMENTS
!
! ON OUTPUT              y
!                          Contains the solution x.
!
!                        ierror
!                          An error flag that indicates invalid
!                          input parameters.  except for number zer0,
!                          a solution is not attempted.
!
!                        = 0  no error.
!                        = 1  m < than 5
!                        = 2  n < than 5
!                        = 3  idimy < m.
!                        = 4  blktri failed while computing results
!                             that depend on the coefficient arrays
!                             an, bn, cn. Check these arrays.
!                        = 5  an(j)*cn(j-1) is less than 0 for some j.
!
!                             Possible reasons for this condition are
!                             1. The arrays an and cn are not correct
!                             2. Too large a grid spacing was used
!                                in the discretization of the elliptic
!                                equation.
!                             3. The linear equations resulted from a
!                                partial differential equation which
!                                was not elliptic.
!
!                        = 20 If the dynamic allocation of real and
!                             complex work space in the derived type
!                             (FishpackWorkspace) variable W fails (e.g.,
!                             if N, M are too large for the platform used)
!
!
!                        workspace
!                             The derived type(FishpackWorkspace) variable
!                             contains real and complex array components that
!                             must not be destroyed if blktri is called again with
!                             iflg=1.
!
!
! SPECIAL CONDITIONS     The algorithm may fail if abs(bm(i)+bn(j))
!                        is less than abs(am(i))+abs(an(j))+
!                        abs(cm(i))+abs(cn(j))
!                        for some i and j. the algorithm will also
!                        fail if an(j)*cn(j-1) is less than zero for
!                        some j.
!                        see the description of the output parameter
!                        ierror.
!
! I/O                    None
!
! PRECISION              64-bit precision float and 32-bit precision integer
!
! REQUIRED FILES         type_FishpackWorkspace.f90, type_ComfAux.f90
!
! STANDARD               Fortran 2008
!
! HISTORY                * Written by Paul Swarztrauber at NCAR in the
!                          early 1970's.
!                        * Rewritten and released in libraries in January 1980.
!                        * Revised in June 2004 using Fortan 90 dynamically
!                          allocated workspace and derived data types to
!                          eliminate mixed mode conflicts in the earlier versions.
!                        * Revised in April 2016 to implement features of
!                          Fortran 2008
!
! ALGORITHM              Generalized cyclic reduction
!
! PORTABILITY            Approximate machine accuracy is obtained
!                        using EPS which is set by the intrinsic epsilon function
!
! REFERENCES             Swarztrauber, P. and R. Sweet, 'Efficient
!                        fortran subprograms for the solution of
!                        elliptic equations'
!                        NCAR TN/IA-109, July, 1975, 138 pp.
!
!                        Swarztrauber P. N., A direct method for
!                        the discrete solution of separable
!                        elliptic equations, SIAM
!                        J. Numer. Anal., 11(1974) pp. 1136-1150.
!
module module_blktri

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        MACHINE_EPSILON ! Machine epsilon

    use type_FishpackWorkspace, only: &
        Fish => FishpackWorkspace

    use type_ComfAux, only: &
        ComfAux, &
        psgf, &
        ppspf, &
        comf_interface

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: blktri
    public :: blktrii
    public :: BlktriAux

    type, public, extends(ComfAux) :: BlktriAux
        !-------------------------------------------------
        ! Type components
        !-------------------------------------------------
        integer(ip) :: indices(6), nl
        integer(ip) :: npp, k, nm, ncmplx, ik
        real(wp)    :: cnv
        real(wp)    :: MACHINE_EPSILON = MACHINE_EPSILON
        !-------------------------------------------------
    contains
        !-------------------------------------------------
        ! Type-bound procedures
        !-------------------------------------------------
        procedure, public  :: blktrii
        procedure, private :: blktri_lower_routine
        procedure, private :: bsrh
        procedure, private :: compb
        procedure, private :: ppadd
        procedure, private :: tevls
        procedure, private :: indxa
        procedure, private :: indxb
        procedure, private :: indxc
        !-------------------------------------------------
    end type BlktriAux

    !---------------------------------------------------------------------------------
    ! Parameters confined to the module
    !---------------------------------------------------------------------------------
    real(wp),    parameter :: ZERO = 0.0_wp
    real(wp),    parameter :: ONE = 1.0_wp
    real(wp),    parameter :: TWO = 2.0_wp
    !---------------------------------------------------------------------------------

contains

    subroutine blktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, &
        idimy, y, ierror, workspace)
        !--------------------------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------------------------
        integer(ip), intent(in)    :: iflg
        integer(ip), intent(in)    :: np
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: mp
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: idimy
        integer(ip), intent(out)   :: ierror
        real(wp),    intent(inout) :: an(:)
        real(wp),    intent(inout) :: bn(:)
        real(wp),    intent(inout) :: cn(:)
        real(wp),    intent(inout) :: am(:)
        real(wp),    intent(inout) :: bm(:)
        real(wp),    intent(inout) :: cm(:)
        real(wp),    intent(inout) :: y(:,:)
        class(Fish), intent(inout) :: workspace
        !--------------------------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------------------------
        integer(ip)           :: irwk, icwk
        type(BlktriAux), save :: self
        !--------------------------------------------------------------------------------

        ! Check input arguments
        call check_input_arguments(n, m, idimy, ierror)

        ! Check error flag
        if (ierror /= 0) return

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv, &
            EPS => self%MACHINE_EPSILON &
            )

            if (iflg == 0) then

                ! compute and allocate real and complex required work space
                call workspace%compute_blktri_workspace_lengths(n, m, irwk, icwk)

                ! Allocate memory for workspace
                call workspace%create(irwk, icwk)

            end if

            ! Solve system
            associate( &
                rew => workspace%real_workspace, &
                cxw => workspace%complex_workspace &
                )
                call self%blktrii(iflg, np, n, an, bn, cn, &
                    mp, m, am, bm, cm, idimy, y, ierror, rew, cxw)
            end associate

        end associate common_variables

    end subroutine blktri

    subroutine check_input_arguments(n, m, idimy, ierror)
        !--------------------------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------------------------
        integer(ip), intent(in)  :: n, m, idimy
        integer(ip), intent(out) :: ierror
        !--------------------------------------------------------------------------------

        if (m < 5) then
            ierror = 1
            return
        else if (n < 3) then
            ierror = 2
            return
        else if (idimy < m) then
            ierror = 3
            return
        else
            ierror = 0
        end if

    end subroutine check_input_arguments

    subroutine blktrii(self, iflg, np, n, an, bn, cn, mp, m, am, bm, cm, &
        idimy, y, ierror, w, wc)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        class(BlktriAux), intent(inout) :: self
        integer(ip),      intent(in)    :: iflg
        integer(ip),      intent(in)    :: np
        integer(ip),      intent(in)    :: n
        integer(ip),      intent(in)    :: mp
        integer(ip),      intent(in)    :: m
        integer(ip),      intent(in)    :: idimy
        integer(ip),      intent(out)   :: ierror
        real(wp),         intent(inout) :: an(:)
        real(wp),         intent(inout) :: bn(:)
        real(wp),         intent(inout) :: cn(:)
        real(wp),         intent(inout) :: am(:)
        real(wp),         intent(inout) :: bm(:)
        real(wp),         intent(inout) :: cm(:)
        real(wp),         intent(inout) :: y(:,:)
        real(wp),         intent(inout) :: w(:)
        complex(wp),      intent(inout) :: wc(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: nh, iwah, iwbh
        !----------------------------------------------

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv, &
            EPS => self%MACHINE_EPSILON &
            )

            ! test m and n for the proper form
            nm = n

            ! Check again for solvers which call blktrii directly
            call check_input_arguments(nm, m, idimy, ierror)

            ! Check error flag
            if (ierror /= 0) return

            associate( &
                iw1 => self%indices(1), &
                iw2 => self%indices(2), &
                iw3 => self%indices(3), &
                iww => self%indices(4), &
                iwu => self%indices(5), &
                iwd => self%indices(6), &
                nl => self%nl &
                )

                select case (iflg)
                    case (0)
                        nh = n
                        npp = np
                        if (npp /= 0)  nh = nh + 1
                        ik = 4
                        k = 2

                        do
                            if (nh <= ik) exit
                            ik = 2*ik
                            k = k + 1
                        end do

                        nl = ik
                        ik = 2*ik
                        nl = nl - 1
                        iwah = (k - 2)*ik + k + 5

                        if (npp == 0) then
                            iwbh = iwah + 2*nm
                            iw1 = iwbh
                            nm = n - 1
                        else
                            iw1 = iwah
                            iwbh = iw1 + nm
                        end if
                        !
                        !==> Set workspace indices
                        !
                        iw2 = iw1 + m
                        iw3 = iw2 + m
                        iwd = iw3 + m
                        iww = iwd + m
                        iwu = iww + m

                        call self%compb(nl, ierror, an, bn, cn, w, wc, w(iwah:), w(iwbh:))

                    case default

                        ! *** Important to reset nm for np = 0
                        if (npp == 0) nm = n - 1

                        select case (mp)
                            case (0)
                                call self%blktri_lower_routine(nl, an, bn, cn, m, am, bm, cm, idimy, y, w, wc, &
                                    w(iw1:), w(iw2:), w(iw3:), w(iwd:), w(iww:), w(iwu:), wc(iw1:), &
                                    wc(iw2:), wc(iw3:), prodp, cprodp)
                            case default
                                call self%blktri_lower_routine(nl, an, bn, cn, m, am, bm, cm, idimy, y, w, wc, &
                                    w(iw1:), w(iw2:), w(iw3:), w(iwd:), w(iww:), w(iwu:), wc(iw1:), &
                                    wc(iw2:), wc(iw3:), prod, cprod)
                        end select
                end select

            end associate

        end associate common_variables

    end subroutine blktrii

    subroutine blktri_lower_routine(self, n, an, bn, cn, m, am, bm, cm, idimy, y, b, bc, &
        w1, w2, w3, wd, ww, wu, cw1, cw2, cw3, prdct, cprdct)
        !
        ! Purpose:
        !
        ! blktri_lower_routine solves the linear system
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
        ! Dummy arguments
        !-----------------------------------------------
        class(BlktriAux), intent(inout) :: self
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: idimy
        real(wp),    intent(in)     :: an(:)
        real(wp),    intent(in)     :: bn(:)
        real(wp),    intent(in)     :: cn(:)
        real(wp),    intent(in)     :: am(:)
        real(wp),    intent(in)     :: bm(:)
        real(wp),    intent(in)     :: cm(:)
        real(wp),    intent(inout)  :: y(:,:)
        real(wp),    intent(in)     :: b(*)
        real(wp),    intent(inout)  :: w1(*)
        real(wp),    intent(inout)  :: w2(*)
        real(wp),    intent(inout)  :: w3(*)
        real(wp),    intent(in)     :: wd(*)
        real(wp),    intent(in)     :: ww(*)
        real(wp),    intent(in)     :: wu(*)
        complex(wp), intent(in)     :: bc(*)
        complex(wp), intent(in)     :: cw1(*)
        complex(wp), intent(in)     :: cw2(*)
        complex(wp), intent(in)     :: cw3(*)
        external :: prdct, cprdct
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: kdo, l, ir, i2, i1, i3, i4, irm1, im2, nm2, im3, nm3
        integer(ip) :: im1, nm1, i0, iif, i, ipi1, ipi2, ipi3, idxc, nc, idxa, na, ip2
        integer(ip) :: np2, ip1, np1, ip3, np3, iz, nz, izr, ll, ifd, iip, np
        integer(ip) :: imi1, imi2
        real(wp)    :: dum
        !-----------------------------------------------

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv, &
            EPS => self%MACHINE_EPSILON &
            )
            !
            !==> begin reduction phase
            !
            kdo = k - 1
            do l = 1, kdo
                ir = l - 1
                i2 = 2**ir
                i1 = i2/2
                i3 = i2 + i1
                i4 = i2 + i2
                irm1 = ir - 1
                call self%indxb(i2, ir, im2, nm2)
                call self%indxb(i1, irm1, im3, nm3)
                call self%indxb(i3, irm1, im1, nm1)

                i0 = 0

                call prdct(nm2, b(im2), nm3, b(im3), nm1, b(im1), i0, dum, &
                    y(1, i2), w3, m, am, bm, cm, wd, ww, wu)

                iif = 2**k

                do i = i4, iif, i4

                    if (i > nm) cycle

                    ipi1 = i + i1
                    ipi2 = i + i2
                    ipi3 = i + i3

                    call self%indxc(i, ir, idxc, nc)

                    if (i >= iif) cycle

                    call self%indxa(i, ir, idxa, na)
                    call self%indxb(i - i1, irm1, im1, nm1)
                    call self%indxb(ipi2, ir, ip2, np2)
                    call self%indxb(ipi1, irm1, ip1, np1)
                    call self%indxb(ipi3, irm1, ip3, np3)
                    call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), w3, &
                        w1, m, am, bm, cm, wd, ww, wu)

                    if (ipi2 > nm) then
                        w3(:m) = ZERO
                        w2(:m) = ZERO
                    else
                        call prdct(np2, b(ip2), np1, b(ip1), np3, b(ip3), 0, dum, &
                            y(1, ipi2), w3, m, am, bm, cm, wd, ww, wu)
                        call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), w3, &
                            w2, m, am, bm, cm, wd, ww, wu)
                    end if
                    y(:m, i) = w1(:m) + w2(:m) + y(:m, i)
                end do
            end do

            if (npp == 0) then
                iif = 2**k
                i = iif/2
                i1 = i/2

                call self%indxb(i - i1, k - 2, im1, nm1)
                call self%indxb(i + i1, k - 2, ip1, np1)
                call self%indxb(i, k - 1, iz, nz)
                call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dum, &
                    y(1, i), w1, m, am, bm, cm, wd, ww, wu)

                izr = i
                w2(:m) = w1(:m)

                do ll = 2, k
                    l = k - ll + 1
                    ir = l - 1
                    i2 = 2**ir
                    i1 = i2/2
                    i = i2

                    call self%indxc(i, ir, idxc, nc)
                    call self%indxb(i, ir, iz, nz)
                    call self%indxb(i - i1, ir - 1, im1, nm1)
                    call self%indxb(i + i1, ir - 1, ip1, np1)
                    call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), w1, &
                        w1, m, am, bm, cm, wd, ww, wu)

                    w1(:m) = y(:m, i) + w1(:m)

                    call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, &
                        dum, w1, w1, m, am, bm, cm, wd, ww, wu)
                end do

                outer_loop: do ll = 2, k
                    l = k - ll + 1
                    ir = l - 1
                    i2 = 2**ir
                    i1 = i2/2
                    i4 = i2 + i2
                    ifd = iif - i2
                    inner_loop: do i = i2, ifd, i4

                        if (i - i2 - izr /= 0) cycle inner_loop

                        if (i > nm) cycle outer_loop

                        call self%indxa(i, ir, idxa, na)
                        call self%indxb(i, ir, iz, nz)
                        call self%indxb(i - i1, ir - 1, im1, nm1)
                        call self%indxb(i + i1, ir - 1, ip1, np1)
                        call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), &
                            w2, w2, m, am, bm, cm, wd, ww, wu)

                        w2(:m) = y(:m, i) + w2(:m)

                        call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dum, &
                            w2, w2, m, am, bm, cm, wd, ww, wu)

                        izr = i

                        if (i == nm) exit outer_loop

                    end do inner_loop
                end do outer_loop

                y(:m, nm+1) = y(:m, nm+1) - cn(nm+1)*w1(:m) - an(nm+1)*w2(:m)

                call self%indxb(iif/2, k - 1, im1, nm1)
                call self%indxb(iif, k - 1, iip, np)

                if (ncmplx /= 0) then
                    call cprdct(nm + 1, bc(iip), nm1, b(im1), 0, dum, 0, dum, &
                        y(1, nm+1), y(1, nm+1), m, am, bm, cm, cw1, cw2, cw3)
                else
                    call prdct(nm + 1, b(iip), nm1, b(im1), 0, dum, 0, dum, &
                        y(1, nm+1), y(1, nm+1), m, am, bm, cm, wd, ww, wu)
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

                    call self%indxa(i, ir, idxa, na)
                    call self%indxb(i - i2, ir, im2, nm2)
                    call self%indxb(i - i2 - i1, ir - 1, im3, nm3)
                    call self%indxb(i - i1, ir - 1, im1, nm1)
                    call prdct(nm2, b(im2), nm3, b(im3), nm1, b(im1), 0, dum, &
                        w1, w1, m, am, bm, cm, wd, ww, wu)
                    call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), w1, &
                        w1, m, am, bm, cm, wd, ww, wu)

                    y(:m, i) = y(:m, i) - w1(:m)
                end do

                izr = nm

                loop_131: do l = 1, kdo
                    ir = l - 1
                    i2 = 2**ir
                    i1 = i2/2
                    i3 = i2 + i1
                    i4 = i2 + i2
                    irm1 = ir - 1
                    loop_132: do i = i4, iif, i4
                        ipi1 = i + i1
                        ipi2 = i + i2
                        ipi3 = i + i3

                        if (ipi2 /= izr) then

                            if (i /= izr) cycle loop_132

                            cycle loop_131

                        end if

                        call self%indxc(i, ir, idxc, nc)
                        call self%indxb(ipi2, ir, ip2, np2)
                        call self%indxb(ipi1, irm1, ip1, np1)
                        call self%indxb(ipi3, irm1, ip3, np3)

                        call prdct(np2, b(ip2), np1, b(ip1), np3, b(ip3), 0, &
                            dum, w2, w2, m, am, bm, cm, wd, ww, wu)

                        call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), &
                            w2, w2, m, am, bm, cm, wd, ww, wu)

                        y(:m, i) = y(:m, i) - w2(:m)
                        izr = i

                        cycle loop_131
                    end do loop_132
                end do loop_131
            end if
            !
            !==> begin back substitution phase
            !
            do ll = 1, k
                l = k - ll + 1
                ir = l - 1
                irm1 = ir - 1
                i2 = 2**ir
                i1 = i2/2
                i4 = i2 + i2
                ifd = iif - i2
                inner_back_sub: do i = i2, ifd, i4

                    if (i > nm) cycle inner_back_sub

                    imi1 = i - i1
                    imi2 = i - i2
                    ipi1 = i + i1
                    ipi2 = i + i2

                    call self%indxa(i, ir, idxa, na)
                    call self%indxc(i, ir, idxc, nc)
                    call self%indxb(i, ir, iz, nz)
                    call self%indxb(imi1, irm1, im1, nm1)
                    call self%indxb(ipi1, irm1, ip1, np1)

                    if (i <= i2) then
                        w1(:m) = ZERO
                    else
                        call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), &
                            y(1,imi2), w1, m, am, bm, cm, wd, ww, wu)
                    end if

                    if (ipi2 > nm) then
                        w2(:m) = ZERO
                    else
                        call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), &
                            y(1,ipi2), w2, m, am, bm, cm, wd, ww, wu)
                    end if

                    w1(:m) = y(:m, i) + w1(:m) + w2(:m)

                    call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dum, &
                        w1, y(1, i), m, am, bm, cm, wd, ww, wu)
                end do inner_back_sub
            end do

        end associate common_variables

    end subroutine blktri_lower_routine

    function bsrh(self, xll, xrr, iz, c, a, bh, sgn, f) result (return_value)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        class(BlktriAux), intent(inout) :: self
        real(wp),         intent(in)    :: xll
        real(wp),         intent(in)    :: xrr
        integer(ip),      intent(in)    :: iz
        real(wp),         intent(in)    :: c(:)
        real(wp),         intent(in)    :: a(:)
        real(wp),         intent(in)    :: bh(:)
        real(wp),         intent(in)    :: sgn
        procedure(comf_interface)      :: f
        real(wp)                        :: return_value
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        real(wp) :: r1, xl, xr, dx, x
        !-----------------------------------------------

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv, &
            EPS => self%MACHINE_EPSILON &
            )

            xl = xll
            xr = xrr
            dx = abs(xr - xl)/2
            x = (xl + xr)/2
            r1 = sgn * f(x, iz, c, a, bh)

            if (r1 >= ZERO) then
                if (r1 == ZERO) then
                    return_value = (xl + xr)/2
                    return
                end if
                xr = x
            else
                xl = x
            end if

            dx = dx/2

            do
                if (dx - cnv > ZERO) exit
                x = (xl + xr)/2
                r1 = sgn * f(x, iz, c, a, bh)
                if (r1 >= ZERO) then
                    if (r1 == ZERO) then
                        return_value = (xl + xr)/2
                        return
                    end if
                    xr = x
                else
                    xl = x
                end if
                dx = dx/2
            end do

            return_value = (xl + xr)/2

        end associate common_variables

    end function bsrh

    subroutine compb(self, n, ierror, an, bn, cn, b, bc, ah, bh)
        !
        ! Purpose:
        !
        !     compb computes the roots of the b polynomials using subroutine
        !     tevls which is a modification the eispack program tqlrat.
        !     ierror is set to 4 if either tevls fails or if a(j+1)*c(j) is
        !     less than zero for some j.  ah, bh are temporary work arrays.
        !
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        class(BlktriAux), intent(inout) :: self
        integer(ip),      intent(in)    :: n
        integer(ip),      intent(out)   :: ierror
        real(wp),         intent(in)    :: an(:)
        real(wp),         intent(in)    :: bn(:)
        real(wp),         intent(in)    :: cn(:)
        real(wp),         intent(inout) :: b(:)
        real(wp),         intent(inout) :: ah(:)
        real(wp),         intent(inout) :: bh(:)
        complex(wp),      intent(inout) :: bc(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip)  :: j, iif, kdo, l, ir, i2, i4, ipl, ifd, i, ib, nb, js, jf
        integer(ip)  :: ls, lh, nmp, l1, l2, j2, j1, n2m2
        real(wp)     :: bnorm, arg, d1, d2, d3
        !-----------------------------------------------

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv, &
            EPS => self%MACHINE_EPSILON &
            )

            bnorm = abs(bn(1))

            do j = 2, nm
                bnorm = max(bnorm, abs(bn(j)))
                arg = an(j)*cn(j-1)

                if (arg < ZERO) then
                    ierror = 5
                    return
                end if

                b(j) = sign(sqrt(arg), an(j))
            end do

            cnv = EPS*bnorm
            iif = 2**k
            kdo = k - 1

            outer_loop: do l = 1, kdo
                ir = l - 1
                i2 = 2**ir
                i4 = i2 + i2
                ipl = i4 - 1
                ifd = iif - i4
                do i = i4, ifd, i4
                    call self%indxb(i, l, ib, nb)

                    if (nb <= 0) cycle outer_loop

                    js = i - ipl
                    jf = js + nb - 1
                    ls = 0
                    bh(:jf-js+1) = bn(js:jf)
                    ah(:jf-js+1) = b(js:jf)

                    associate( order => nb )
                        call self%tevls(bh(1:order), ah(1:order), ierror)
                    end associate

                    if (ierror /= 0) then
                        ierror = 4
                        return
                    end if

                    lh = ib - 1
                    if (nb > 0) then
                        b(lh+1:nb+lh) = -bh(:nb)
                        lh = nb + lh
                    end if
                end do
            end do outer_loop

            b(:nm) = -bn(:nm)

            if (npp == 0) then
                nmp = nm + 1
                nb = nm + nmp
                do j = 1, nb
                    l1 = mod(j - 1, nmp) + 1
                    l2 = mod(j + nm - 1, nmp) + 1
                    arg = an(l1)*cn(l2)

                    if (arg < ZERO) then
                        ierror = 5
                        return
                    end if

                    bh(j) = sign(sqrt(arg), (-an(l1)))
                    ah(j) = -bn(l1)
                end do

                associate( order => nb )
                    call self%tevls(ah(1:order), bh(1:order), ierror)
                end associate

                if (ierror /= 0) then
                    ierror = 4
                    return
                end if

                call self%indxb(iif, k - 1, j2, lh)
                call self%indxb(iif/2, k - 1, j1, lh)

                j2 = j2 + 1
                lh = j2
                n2m2 = j2 + 2*nm - 2

                do
                    if (j2 <= n2m2) exit

                    d1 = abs(b(j1)-b(j2-1))
                    d2 = abs(b(j1)-b(j2))
                    d3 = abs(b(j1)-b(j2+1))

                    if (d1 <= d2 .or. d3 <= d2) then
                        b(lh) = b(j2)
                        j2 = j2 + 1
                        lh = lh + 1
                    else
                        j2 = j2 + 1
                        j1 = j1 + 1
                    end if
                end do

                b(lh) = b(n2m2+1)

                call self%indxb(iif, k - 1, j1, j2)

                j2 = j1 + nmp + nmp

                associate( order => nm + 1 )
                    call self%ppadd(ierror, an(1:order), cn(1:order), bc(j1:order), b(j1:order), b(j2:order))
                end associate

            end if

        end associate common_variables

    end subroutine compb

    pure subroutine cprod(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, yy, &
        m, a, b, c, d, w, y)
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
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in)  :: nd
        integer(ip), intent(in)  :: nm1
        integer(ip), intent(in)  :: nm2
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: bm1(nm1)
        real(wp),    intent(in)  :: bm2(nm2)
        real(wp),    intent(in)  :: aa(na)
        real(wp),    intent(in)  :: x(m)
        real(wp),    intent(out) :: yy(m)
        real(wp),    intent(in)  :: a(m)
        real(wp),    intent(in)  :: b(m)
        real(wp),    intent(in)  :: c(m)
        complex(wp), intent(in)  :: bd(nd)
        complex(wp), intent(out) :: d(m)
        complex(wp), intent(out) :: w(m)
        complex(wp), intent(out) :: y(m)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: j, mm, id, m1, m2, ia, iflg, k
        real(wp)    :: rt
        complex(wp) :: crt, den, y1, y2
        !-----------------------------------------------

        y = cmplx(x, ZERO, kind=wp)

        mm = m - 1
        id = nd
        m1 = nm1
        m2 = nm2
        ia = na

        main_loop: do

            iflg = 0

            if (id > 0) then
                crt = bd(id)
                id = id - 1
                !
                !==> begin solution to system
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

                if (abs(den) /= ZERO) then
                    y(1) = (y(1)-c(1)*w(2))/den
                else
                    y(1) = cmplx(ONE, ZERO, kind=wp)
                end if

                y(2:m) = w(2:m) - d(2:m)*y(1:m-1)

            end if

            if (m1 > 0 .and. m2 > 0) then
                if (m1 <= 0) then
                    rt = bm2(m2)
                    m2 = m2 - 1
                else
                    if (m2 <= 0) then
                        rt = bm1(m1)
                        m1 = m1 - 1
                    else
                        if (abs(bm1(m1)) - abs(bm2(m2)) > ZERO) then
                            rt = bm1(m1)
                            m1 = m1 - 1
                        else
                            rt = bm2(m2)
                            m2 = m2 - 1
                        end if
                    end if
                end if
                y1 = (b(1)-rt)*y(1) + c(1)*y(2)
                if (mm >= 2) then
                    do j = 2, mm
                        y2 = a(j)*y(j-1) + (b(j)-rt)*y(j) + c(j)*y(j+1)
                        y(j-1) = y1
                        y1 = y2
                    end do
                end if
                y(m) = a(m)*y(m-1) + (b(m)-rt)*y(m)
                y(m-1) = y1
                iflg = 1
                cycle main_loop
            end if

            if (ia > 0) then
                rt = aa(ia)
                ia = ia - 1
                iflg = 1
                !
                !==> scalar multiplication
                !
                y = rt*y
            end if

            if (iflg <= 0) exit main_loop

        end do main_loop

        yy = real(y, kind=wp)

    end subroutine cprod

    pure subroutine cprodp(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, yy, m, a, &
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
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in)  :: nd
        integer(ip), intent(in)  :: nm1
        integer(ip), intent(in)  :: nm2
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: bm1(nm1)
        real(wp),    intent(in)  :: bm2(nm2)
        real(wp),    intent(in)  :: aa(na)
        real(wp),    intent(in)  :: x(m)
        real(wp),    intent(out) :: yy(m)
        real(wp),    intent(in)  :: a(m)
        real(wp),    intent(in)  :: b(m)
        real(wp),    intent(in)  :: c(m)
        complex(wp), intent(in)  :: bd(nd)
        complex(wp), intent(out) :: d(m)
        complex(wp), intent(out) :: u(m)
        complex(wp), intent(out) :: y(m)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: j, mm, mm2, id, m1, m2, ia, iflg, k
        real(wp)    :: rt
        complex(wp) :: v, den, bh, ym, am, y1, y2, yh, crt
        !-----------------------------------------------

        y = cmplx(x, ZERO, kind=wp)

        mm = m - 1
        mm2 = m - 2
        id = nd
        m1 = nm1
        m2 = nm2
        ia = na

        main_loop: do

            iflg = 0

            if (id > 0) then
                crt = bd(id)
                id = id - 1
                iflg = 1
                !
                !==> begin solution to system
                !
                bh = b(m) - crt
                ym = y(m)
                den = b(1) - crt
                d(1) = c(1)/den
                u(1) = a(1)/den
                y(1) = y(1)/den
                v = cmplx(c(m), ZERO, kind=wp)

                if (mm2 >= 2) then
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

                if (abs(den) /= ZERO) then
                    y(m) = (ym - am*y(m-1))/den
                else
                    y(m) = cmplx(ONE, ZERO, kind=wp)
                end if

                y(m-1) = y(m-1) - d(m-1)*y(m)

                do j = 2, mm
                    k = m - j
                    y(k) = y(k) - d(k)*y(k+1) - u(k)*y(m)
                end do
            end if

            if (m1 > 0 .and. m2 > 0) then
                if (m1 <= 0) then
                    rt = bm2(m2)
                    m2 = m2 - 1
                else
                    if (m2 <= 0) then
                        rt = bm1(m1)
                        m1 = m1 - 1
                    else
                        if (abs(bm1(m1)) - abs(bm2(m2)) > ZERO) then
                            rt = bm1(m1)
                            m1 = m1 - 1
                        else
                            rt = bm2(m2)
                            m2 = m2 - 1
                        !
                        !==> matrix multiplication
                        !
                        end if
                    end if
                end if
                yh = y(1)
                y1 = (b(1)-rt)*y(1) + c(1)*y(2) + a(1)*y(m)
                if (mm >= 2) then
                    do j = 2, mm
                        y2 = a(j)*y(j-1) + (b(j)-rt)*y(j) + c(j)*y(j+1)
                        y(j-1) = y1
                        y1 = y2
                    end do
                end if
                y(m) = a(m)*y(m-1) + (b(m)-rt)*y(m) + c(m)*yh
                y(m-1) = y1
                iflg = 1
                cycle main_loop
            end if

            if (ia > 0) then
                rt = aa(ia)
                ia = ia - 1
                iflg = 1
                !
                !==> scalar multiplication
                !
                y = rt*y
            end if

            if (iflg <= 0) exit main_loop

        end do main_loop

        yy = real(y, kind=wp)


    end subroutine cprodp

    pure subroutine indxa(self, i, ir, idxa, na)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        class(BlktriAux), intent(inout) :: self
        integer(ip),      intent(in)    :: i
        integer(ip),      intent(in)    :: ir
        integer(ip),      intent(out)   :: idxa
        integer(ip),      intent(out)   :: na
        !-----------------------------------------------

        associate( nm => self%nm )
            na = 2**ir
            idxa = i - na + 1
            if (i > nm) na = 0
        end associate

    end subroutine indxa

    pure subroutine indxb(self, i, ir, idx, idp)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        class(BlktriAux), intent(inout) :: self
        integer(ip),      intent(in)    :: i
        integer(ip),      intent(in)    :: ir
        integer(ip),      intent(out)   :: idx
        integer(ip),      intent(out)   :: idp
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: izh, id, ipl
        !-----------------------------------------------

        associate( &
            nm => self%nm, &
            ik => self%ik &
            )
            !
            ! b(idx) is the location of the first root of the b(i, ir) polynomial
            !
            idp = 0

            select case (ir)
                case (0)
                    if (i > nm) then
                        return
                    else
                        idx = i
                        idp = 1
                    end if
                case default
                    if (ir > 0) then
                        izh = 2**ir
                        id = i - 2*izh
                        idx = 2*id + (ir - 1)*ik + ir + (ik - i)/izh + 4
                        ipl = izh - 1
                        idp = 2*izh - 1
                        if (i - ipl - nm > 0) then
                            idp = 0
                            return
                        end if
                        if (i + ipl - nm > 0) idp = nm + ipl - i + 1
                    end if
            end select
        end associate

    end subroutine indxb

    pure subroutine indxc(self, i, ir, idxc, nc)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        class(BlktriAux), intent(inout) :: self
        integer(ip),      intent(in)    :: i
        integer(ip),      intent(in)    :: ir
        integer(ip),      intent(out)   :: idxc
        integer(ip),      intent(out)   :: nc
        !-----------------------------------------------

        associate( nm => self%nm )
            nc = 2**ir
            idxc = i
            if (idxc + nc - 1 - nm > 0) nc = 0
        end associate

    end subroutine indxc

    subroutine ppadd(self, ierror, a, c, cbp, bp, bh)
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
        ! Dummy arguments
        !-----------------------------------------------
        class(BlktriAux), intent(inout) :: self
        integer(ip),      intent(out)   :: ierror
        real(wp),         intent(in)    :: a(:)
        real(wp),         intent(in)    :: c(:)
        real(wp),         intent(inout) :: bp(:)
        real(wp),         intent(inout) :: bh(:)
        complex(wp),      intent(inout) :: cbp(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip)   :: iz, izm, izm2, j, nt, modiz, is
        integer(ip)   :: iif, ig, it, icv, i3, i2, nhalf
        real(wp)      :: r4, r5, r6, scnv, xl, db, sgn, xr, xm, psg
        real(wp)      :: temp
        complex(wp)   :: cx, fsg, hsg, dd, f, fp, fpp, cdis, r1, r2, r3
        !-----------------------------------------------

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv, &
            EPS => self%MACHINE_EPSILON &
            )

            associate( n => size(a) )

                scnv = sqrt(cnv)
                iz = n
                izm = iz - 1
                izm2 = iz - 2

                if (bp(n) <= bp(1)) then
                    if (bp(n) == bp(1)) then
                        ierror = 4
                        return
                    else
                        bh(:n) = bp(n:1:(-1))
                    end if
                else
                    bh(:n) = bp(:n)
                end if

                ncmplx = 0
                modiz = mod(iz, 2)
                is = 1

                if (modiz /= 0) then
                    if (a(1) >= ZERO) then
                        if (a(1) == ZERO) then
                            ierror = 4
                            return
                        end if
                    end if

                    xl = bh(1)
                    db = bh(3) - bh(1)
                    xl = xl - db
                    r4 = self%psgf(xl, iz, c, a, bh)

                    do
                        if (ZERO < r4) exit
                        xl = xl - db
                        r4 = self%psgf(xl, iz, c, a, bh)
                    end do

                    sgn = -ONE

                    temp = self%bsrh(xl, bh(1), iz, c, a, bh, sgn, psgf)

                    cbp(1) = cmplx(temp, ZERO, kind=wp)

                    bp(1) = real(cbp(1), kind=wp)

                    is = 2

                end if

                iif = iz - 1

                if (modiz /= 0) then

                    if (a(1) == ZERO) then
                        ierror = 4
                        return
                    end if

                    xr = bh(iz)
                    db = bh(iz) - bh(iz-2)
                    xr = xr + db
                    r5 = self%psgf(xr, iz, c, a, bh)

                    do
                        if (ZERO <= r5) exit
                        xr = xr + db
                        r5 = self%psgf(xr, iz, c, a, bh)
                    end do

                    sgn = ONE
                    temp = self%bsrh(bh(iz), xr, iz, c, a, bh, sgn, psgf)

                    cbp(iz) = cmplx(temp, ZERO, kind=wp)
                    iif = iz - 2

                end if

                main_loop: do ig = is, iif, 2

                    xl = bh(ig)
                    xr = bh(ig+1)
                    sgn = -ONE
                    xm = self%bsrh(xl, xr, iz, c, a, bh, sgn, ppspf)
                    psg = self%psgf(xm, iz, c, a, bh)

                    block_construct: block

                        if (abs(psg) > EPS) then

                            r6 = psg*self%psgf(xm, iz, c, a, bh)

                            if (r6 > ZERO) exit block_construct

                            if (r6 /= ZERO) then

                                sgn = ONE
                                temp = self%bsrh(bh(ig), xm, iz, c, a, bh, sgn, psgf)
                                cbp(ig) = cmplx(temp, ZERO, kind=wp)
                                sgn = -ONE
                                temp = self%bsrh(xm, bh(ig+1), iz, c, a, bh, sgn, psgf)
                                cbp(ig+1) = cmplx(temp, ZERO, kind=wp)

                                cycle main_loop
                            !
                            !==> Case of a multiple zero
                            !
                            end if
                        end if

                        cbp(ig) = cmplx(xm, ZERO, kind=wp)
                        cbp(ig+1) = cmplx(xm, ZERO, kind=wp)

                        cycle main_loop
                    !
                    !==> case of a complex zero
                    !
                    end block block_construct

                    it = 0
                    icv = 0
                    cx = cmplx(xm, ZERO, kind=wp)

                    loop_120: do

                        fsg = cmplx(ONE, ZERO, kind=wp)
                        hsg = cmplx(ONE, ZERO, kind=wp)
                        fp = ZERO
                        fpp = ZERO

                        do j = 1, iz
                            dd = ONE/(cx - bh(j))
                            fsg = fsg*a(j)*dd
                            hsg = hsg*c(j)*dd
                            fp = fp + dd
                            fpp = fpp - dd**2
                        end do

                        select case (modiz)
                            case (0)
                                f = cmplx(ONE, ZERO, kind=wp) - fsg - hsg
                            case default
                                f = cmplx(ONE, ZERO, kind=wp) + fsg + hsg
                        end select

                        i3 = 0

                        if (abs(fp) > ZERO) then
                            i3 = 1
                            r3 = -f/fp
                        end if

                        i2 = 0

                        if (abs(fpp) > ZERO) then

                            i2 = 1
                            cdis = sqrt((fp**2) - TWO * f * fpp)
                            r1 = cdis - fp
                            r2 = (-fp) - cdis

                            if (abs(r1) - abs(r2) > ZERO) then
                                r1 = r1/fpp
                            else
                                r1 = r2/fpp
                            end if

                            r2 = ((TWO*f)/fpp)/r1

                            if (abs(r2) < abs(r1)) r1 = r2

                            if (i3 > 0 .and. abs(r3) < abs(r1)) r1 = r3

                        else
                            r1 = r3
                        end if

                        cx = cx + r1
                        it = it + 1

                        if (it > 50) then
                            ierror = 4
                            return
                        end if

                        if (abs(r1) > scnv) cycle loop_120

                        if (icv <= 0) then
                            icv = 1
                            cycle loop_120
                        end if

                        exit loop_120
                    end do loop_120

                    cbp(ig) = cx
                    cbp(ig+1) = conjg(cx)

                end do main_loop

                if (abs(cbp(n)) - abs(cbp(1)) <= ZERO) then
                    if (abs(cbp(n)) - abs(cbp(1)) == ZERO) then
                        ierror = 4
                        return
                    end if

                    nhalf = n/2
                    do j = 1, nhalf
                        nt = n - j
                        cx = cbp(j)
                        cbp(j) = cbp(nt+1)
                        cbp(nt+1) = cx
                    end do
                end if

                ncmplx = 1

                if (any(aimag(cbp(2:iz)) /= ZERO)) return

                ncmplx = 0
                bp(1) = real(cbp(1), kind=wp)
                bp(2:iz) = real(cbp(2:iz))

            end associate

        end associate common_variables

    end subroutine ppadd


    pure subroutine prod(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, w, u)
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
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in)  :: nd
        integer(ip), intent(in)  :: nm1
        integer(ip), intent(in)  :: nm2
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: bd(nd)
        real(wp),    intent(in)  :: bm1(nm1)
        real(wp),    intent(in)  :: bm2(nm2)
        real(wp),    intent(in)  :: aa(na)
        real(wp),    intent(in)  :: x(m)
        real(wp),    intent(out) :: y(m)
        real(wp),    intent(in)  :: a(m)
        real(wp),    intent(in)  :: b(m)
        real(wp),    intent(in)  :: c(m)
        real(wp),    intent(out) :: d(m)
        real(wp),    intent(out) :: w(m)
        real(wp),    intent(out) :: u(m)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: j, mm, id, ibr, m1, m2, ia, k
        real(wp)    :: rt, den
        !-----------------------------------------------

        w = x
        y = w
        mm = m - 1
        id = nd
        ibr = 0
        m1 = nm1
        m2 = nm2
        ia = na

        main_loop: do

            if (ia > 0) then

                if (nd == 0) then
                    rt = -aa(ia)
                else
                    rt = aa(ia)
                end if

                ia = ia - 1
                !
                ! scalar multiplication
                !
                y = rt*w
            end if

            if (id <= 0) return

            rt = bd(id)
            id = id - 1

            if (id == 0) ibr = 1

            !
            !==> begin solution to system
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
            w(1) = ONE

            if (den /= ZERO) w(1) = (y(1)-c(1)*w(2))/den

            do j = 2, m
                w(j) = w(j) - d(j)*w(j-1)
            end do

            if (na > 0) cycle main_loop

            if (m1 <= 0) then
                if (m2 <= 0) then
                    y = w
                    ibr = 1
                    cycle main_loop
                end if
            else
                if (.not.(m2 > 0 .and. abs(bm1(m1)) <= abs(bm2(m2)))) then
                    if (ibr <= 0) then
                        if (abs(bm1(m1)-bd(id)) < abs(bm1(m1)-rt)) then
                            y(:m) = w(:m)
                            ibr = 1
                            cycle main_loop
                        end if
                    end if
                    rt = rt - bm1(m1)
                    m1 = m1 - 1
                    y = y + rt*w
                    cycle main_loop
                end if
            end if

            if (ibr <= 0 .and. abs(bm2(m2)-bd(id)) < abs(bm2(m2)-rt)) then
                y = w
                ibr = 1
                cycle main_loop
            end if

            rt = rt - bm2(m2)
            m2 = m2 - 1
            y = y + rt*w

        end do main_loop

    end subroutine prod

    pure subroutine prodp(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, u, w)
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
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in)  :: nd
        integer(ip), intent(in)  :: nm1
        integer(ip), intent(in)  :: nm2
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: bd(nd)
        real(wp),    intent(in)  :: bm1(nm1)
        real(wp),    intent(in)  :: bm2(nm2)
        real(wp),    intent(in)  :: aa(na)
        real(wp),    intent(in)  :: x(m)
        real(wp),    intent(out) :: y(m)
        real(wp),    intent(in)  :: a(m)
        real(wp),    intent(in)  :: b(m)
        real(wp),    intent(in)  :: c(m)
        real(wp),    intent(out) :: d(m)
        real(wp),    intent(out) :: u(m)
        real(wp),    intent(out) :: w(m)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: j, mm, mm2, id, ibr, m1, m2, ia, k
        real(wp)    :: rt, bh, ym, den, v, am
        !-----------------------------------------------

        y = x
        w = y
        mm = m - 1
        mm2 = m - 2
        id = nd
        ibr = 0
        m1 = nm1
        m2 = nm2
        ia = na

        main_loop: do

            if (ia > 0) then
                rt = aa(ia)

                if (nd == 0) rt = -rt

                ia = ia - 1
                y = rt*w
            end if

            if (id <= 0) return

            rt = bd(id)
            id = id - 1

            if (id == 0) ibr = 1

            !
            !==> begin solution to system
            !
            bh = b(m) - rt
            ym = y(m)
            den = b(1) - rt
            d(1) = c(1)/den
            u(1) = a(1)/den
            w(1) = y(1)/den
            v = c(m)

            if (mm2 >= 2) then
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

            if (den /= ZERO) then
                w(m) = (ym - am*w(m-1))/den
            else
                w(m) = ONE
            end if

            w(m-1) = w(m-1) - d(m-1)*w(m)

            do j = 2, mm
                k = m - j
                w(k) = w(k) - d(k)*w(k+1) - u(k)*w(m)
            end do

            if (na <= 0) then
                if (m1 <= 0 .and. m2 <= 0) then
                    y = w
                    ibr = 1
                    cycle main_loop
                else

                    if (m2 > 0 .and. abs(bm1(m1)) - abs(bm2(m2)) <= ZERO) then
                        if (ibr <= 0 .and. abs(bm2(m2)-bd(id)) - abs(bm2(m2)-rt) < ZERO) then
                            y = w
                            ibr = 1
                            cycle main_loop
                        else
                            rt = rt - bm2(m2)
                            m2 = m2 - 1
                            y = y + rt*w
                        end if
                    end if

                    if (ibr <= 0 .and. abs(bm1(m1)-bd(id)) - abs(bm1(m1)-rt) < ZERO) then
                        y = w
                        ibr = 1
                        cycle main_loop
                    else
                        rt = rt - bm1(m1)
                        m1 = m1 - 1
                        y = y + rt*w
                        cycle main_loop
                    end if
                end if

                if (ibr <= 0 .and. abs(bm2(m2)-bd(id)) - abs(bm2(m2)-rt) < ZERO) then
                    y = w
                    ibr = 1
                else
                    rt = rt - bm2(m2)
                    m2 = m2 - 1
                    y = y + rt*w
                end if
            end if
        end do main_loop

    end subroutine prodp

    subroutine tevls(self, diagonal, subdiagonal, error_flag)
        !
        ! Purpose:
        !
        !     This subroutine is a modification of the eispack subroutine tqlrat
        !     algorithm 464, comm. acm 16, 689(1973) by reinsch.
        !
        !     this subroutine finds the eigenvalues of a symmetric
        !     tridiagonal matrix by the rational ql method.
        !
        !     on input-
        !
        !        d contains the diagonal elements of the input matrix,
        !
        !        e2 contains the subdiagonal elements of the
        !          input matrix in its last n-1 positions.  e2(1) is arbitrary.
        !
        !      on output-
        !
        !        d contains the eigenvalues in ascending order.  if an
        !          error exit is made, the eigenvalues are correct and
        !          ordered for indices 1, 2, ...error_flag-1, but may not be
        !          the smallest eigenvalues,
        !
        !        e2 has been destroyed,
        !
        !        error_flag is set to
        !          zero       for normal return,
        !          j          if the j-th eigenvalue has not been
        !                     determined after 30 iterations.
        !
        !     questions and comments should be directed to b. s. garbow,
        !     applied mathematics division, argonne national laboratory
        !
        !
        !     eps is a machine dependent parameter specifying
        !     the relative precision of floating point arithmetic.
        !
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        class(BlktriAux), intent(inout) :: self
        integer(ip),      intent(out)   :: error_flag
        real(wp),         intent(inout) :: diagonal(:)
        real(wp),         intent(inout) :: subdiagonal(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: i, j, l, m, ii, l1, mml, nhalf, ntop
        real(wp)    :: b, c, f, g, h, p, r, s, dhold
        !-----------------------------------------------

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv, &
            EPS => self%MACHINE_EPSILON &
            )

            error_flag = 0

            associate( &
                n => size(diagonal), &
                d => diagonal, &
                e2 => subdiagonal &
                )

                if (n == 1) return

                e2(:n-1) = e2(2:n)**2
                f = ZERO
                b = ZERO
                e2(n) = ZERO

                main_loop: do l = 1, n
                    j = 0
                    h = EPS*(abs(d(l))+sqrt(e2(l)))

                    if (b <= h) then
                        b = h
                        c = b**2
                    end if
                    !
                    !==>  look for small squared sub-diagonal element
                    !
                    do m = l, n
                        if (e2(m) > c) then
                            cycle
                        end if
                        exit
                    !
                    !==> e2(n) is always zero, so there is no exit
                    !    through the bottom of the loop
                    !
                    end do

                    if_construct: if (m /= l) then

                        loop_105: do

                            if (j == 30) then
                                !
                                !==> set error: no convergence to an
                                !    eigenvalue after 30 iterations
                                !
                                error_flag = l
                                return
                            end if

                            j = j + 1
                            !
                            !==> form shift
                            !
                            l1 = l + 1
                            s = sqrt(e2(l))
                            g = d(l)
                            p = (d(l1)-g)/(2.0*s)
                            r = sqrt(p**2 + ONE)
                            d(l) = s/(p + sign(r, p))
                            h = g - d(l)
                            d(l1:n) = d(l1:n) - h
                            f = f + h
                            !
                            !==> rational ql transformation
                            !
                            if (g == ZERO) then
                                g = b
                            else
                                g = d(m)
                            end if

                            h = g
                            s = ZERO
                            mml = m - l
                            !
                            !==> for i = m-1 step -1 until l do
                            !
                            do ii = 1, mml
                                i = m - ii
                                p = g*h
                                r = p + e2(i)
                                e2(i+1) = s*r
                                s = e2(i)/r
                                d(i+1) = h + s*(h + d(i))
                                g = d(i) - e2(i)/g

                                if (g == ZERO) g = b

                                h = g*p/r
                            end do
                            !
                            e2(l) = s*g
                            d(l) = h
                            !
                            !==> guard against underflowed h
                            !
                            if (h == ZERO) exit if_construct

                            if (abs(e2(l)) <= abs(c/h)) exit if_construct

                            e2(l) = h*e2(l)

                            if (e2(l) == ZERO) exit loop_105

                        end do loop_105
                    end if if_construct

                    p = d(l) + f
                    !
                    !==> order eigenvalues
                    !
                    if (l /= 1) then
                        !
                        !==> for i=l step -1 until 2 do
                        !
                        do ii = 2, l
                            i = l + 2 - ii
                            if (p >= d(i-1)) then
                                d(i) = p
                                cycle main_loop
                            else
                                d(i) = d(i-1)
                            end if
                        end do
                    end if

                    i = 1
                    d(i) = p

                end do main_loop

                if (abs(d(n)) < abs(d(1))) then
                    nhalf = n/2
                    do i = 1, nhalf
                        ntop = n - i
                        dhold = d(i)
                        d(i) = d(ntop+1)
                        d(ntop+1) = dhold
                    end do
                end if

            end associate
        !
        !==> last card of tqlrat
        !
        end associate common_variables

    end subroutine tevls

end module module_blktri
!
! REVISION HISTORY
!
! September 1973    Version 1
! April     1976    Version 2
! January   1978    Version 3
! December  1979    Version 3.1
! February  1985    Documentation upgrade
! November  1988    Version 3.2, FORTRAN 77 changes
! June      2004    Version 5.0, Fortran 90 changes
! April     2016    Fortran 2008 changes
!
