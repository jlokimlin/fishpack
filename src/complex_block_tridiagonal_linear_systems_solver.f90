!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                           Fishpack                            *
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
module complex_block_tridiagonal_linear_systems_solver

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        MACHINE_EPSILON

    use type_GeneralizedCyclicReductionUtility, only: &
        psgf, &
        ppspf, &
        ppsgf, &
        comf_interface, &
        GeneralizedCyclicReductionUtility

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: cblktri

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp

    type, private, extends(GeneralizedCyclicReductionUtility) :: ComplexGeneralizedCyclicReductionUtility
    contains
        ! Public type-bound procedures
        procedure, public  :: cblktri_lower_routine
        ! Private type-bound procedures
        procedure, private :: cblktri_bsrh
        procedure, private :: cblktri_compute_roots_of_b_polynomials
        procedure, private :: cblktri_compute_eigenvalues
        procedure, private :: cblktri_tevls
        procedure, private :: cblktri_compute_index_a_coeff
        procedure, private :: cblktri_compute_index_b_coeff
        procedure, private :: cblktri_compute_index_c_coeff
    end type ComplexGeneralizedCyclicReductionUtility

contains

    ! PURPOSE                cblktri solves a system of linear equations
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
    !                        these equations usually result from the
    !                        discretization of separable elliptic
    !                        equations.  boundary conditions may be
    !                        dirichlet, neumann, or periodic.
    !
    !                        cblktri is a complex version of package
    !                        blktri on ulib.
    !
    ! USAGE                  call cblktri(iflg, np, n, an, bn, cn, mp, m, am, bm,
    !                                     cm, idimy, y, ierror, w)
    !
    !
    ! DIMENSION OF           an(n), bn(n), cn(n), am(m), bm(m), cm(m), y(idimy, n)
    ! ARGUMENTS
    !
    ! ARGUMENTS
    !
    ! ON INPUT               iflg
    !
    !                          = 0  initialization only.
    !                               certain quantities that depend on np,
    !                               n, an, bn, and cn are computed and
    !                               stored in the derived data type w
    !
    !                          = 1  the quantities that were computed
    !                               in the initialization are used
    !                               to obtain the solution x(i, j).
    !
    !                               note:
    !                               a call with iflg=0 takes
    !                               approximately one half the time
    !                               as a call with iflg = 1.
    !                               however, the initialization does
    !                               not have to be repeated unless np,
    !                               n, an, bn, or cn change.
    !
    !                        np
    !                          = 0  if an(1) and cn(n) are not zero,
    !                               which corresponds to periodic
    !                               bounary conditions.
    !
    !                          = 1  if an(1) and cn(n) are zero.
    !
    !                        n
    !                          the number of unknowns in the j-direction.
    !                          n must be greater than 4.
    !                          the operation count is proportional to
    !                          mnlog2(n), hence n should be selected
    !                          less than or equal to m.
    !
    !                        an, bn, cn
    !                          one-dimensional arrays of length n
    !                          that specify the coefficients in the
    !                          linear equations given above.
    !
    !                        mp
    !                          = 0  if am(1) and cm(m) are not zero,
    !                               which corresponds to periodic
    !                               boundary conditions.
    !
    !                          = 1  if am(1) = cm(m) = 0  .
    !
    !                        m
    !                          the number of unknowns in the i-direction.
    !                           m must be greater than 4.
    !
    !                        am, bm, cm
    !                          complex one-dimensional arrays of length m
    !                          that specify the coefficients in the linear
    !                          equations given above.
    !
    !                        idimy
    !                          the row (or first) dimension of the
    !                          two-dimensional array y as it appears
    !                          in the program calling cblktri.
    !                          this parameter is used to specify the
    !                          variable dimension of y.
    !                          idimy must be at least m.
    !
    !                        y
    !                          a complex two-dimensional array that
    !                          specifies the values of the right side of
    !                          the linear system of equations given above.
    !                          y must be dimensioned y(idimy, n) with
    !                          idimy >= m.
    !
    !                        w
    !                          A FishpackWorkspace derived data type variable
    !                          which is used internally in cblktri to
    !                          dynamically allocate real and complex workspace
    !                          arrays used in solution. An error flag
    !                          (ierror = 20) is set if the required workspace
    !                          allocation fails (for example if n, m are too large)
    !                          real and complex values are set in the components
    !                          of w on a initial (iflg=0) call to cblktri. These
    !                          must be preserved on non-initial calls (iflg=1)
    !                          to cblktri. This eliminates redundant calculations
    !                          and saves compute time.
    !               ****       IMPORTANT!  the user program calling cblktri should
    !                          include the statement:
    !
    !                               call w%destroy()
    !
    !                          after the final approximation is generated by
    !                          cblktri. The will deallocate the real and complex
    !                          workspace arrays of w. Tailure to include this statement
    !                          could result in serious memory leakage.
    !
    !
    ! ARGUMENTS
    !
    ! ON OUTPUT              y
    !                          contains the solution x.
    !
    !                        ierror
    !                          an error flag that indicates invalid
    !                          input parameters.  except for number zer0,
    !                          a solution is not attempted.
    !
    !                        = 0  no error.
    !                        = 1  m < 5
    !                        = 2  n < 5
    !                        = 3  idimy < m.
    !                        = 4  cblktri failed while computing results
    !                             that depend on the coefficient arrays
    !                             an, bn, cn.  check these arrays.
    !                        = 5  an(j)*cn(j-1) is less than 0 for some j.
    !
    !                             possible reasons for this condition are
    !                             1. the arrays an and cn are not correct
    !                             2. too large a grid spacing was used
    !                                in the discretization of the elliptic
    !                                equation.
    !                             3. the linear equations resulted from a
    !                                partial differential equation which
    !                                was not elliptic.
    !
    !                          = 20 if the dynamic allocation of real and
    !                               complex workspace in the derived type
    !                               (FishpackWorkspace) variable w fails (e.g.,
    !                               if n, m are too large for the platform used)
    !
    !
    !
    ! SPECIAL CONDITIONS     The algorithm may fail if
    !
    !                        abs(bm(i)+bn(j)) < abs(am(i))+abs(an(j)) +
    !                        abs(cm(i))+abs(cn(j))
    !
    !                        for some i and j. the algorithm will also
    !                        fail if an(j)*cn(j-1) < 0 zero for some j.
    !                        see the description of the output parameter
    !                        ierror.
    !
    !
    ! HISTORY               * Written by Paul Swarztrauber at NCAR in
    !                         the early 1970's.  Rewritten an released
    !                         on NCAR's public software libraries in
    !                         January, 1980.
    !                       * Revised in June 2004 by John Adams using
    !                         Fortran 90 dynamically allocated workspace
    !                         and derived data types to eliminate mixed
    !                         mode conflicts in the earlier versions.
    !
    ! ALGORITHM              Generalized cyclic reduction
    !                        (see reference below)
    !
    ! PORTABILITY            The approximate machine accuracy is computed
    !                        by calling the intrinsic function epsilon
    !
    ! REFERENCES             Swarztrauber, P. and R. Sweet, 'Efficient
    !                        FORTRAN subprograms for the solution of
    !                        elliptic equations'
    !                        NCAR TN/IA-109, July, 1975, 138 pp.
    !
    !                        Swarztrauber P. N., A direct method for
    !                        the discrete solution of separable
    !                        elliptic equations, S.I.A.M.
    !                        J. Numer. Anal., 11(1974) pp. 1136-1150.
    !
    subroutine cblktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, &
        idimy, y, ierror, w)

        ! Dummy arguments
        integer(ip), intent(in)    :: iflg
        integer(ip), intent(in)    :: np
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: mp
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: idimy
        integer(ip), intent(out)   :: ierror
        real(wp),    intent(in)    :: an(:)
        real(wp),    intent(in)    :: bn(:)
        real(wp),    intent(in)    :: cn(:)
        complex(wp), intent(in)    :: am(:)
        complex(wp), intent(in)    :: bm(:)
        complex(wp), intent(in)    :: cm(:)
        complex(wp), intent(inout) :: y(:,:)
        class(FishpackWorkspace), intent(inout) :: w

        ! Local variables
        type(ComplexGeneralizedCyclicReductionUtility) :: util
        integer(ip)  :: m2, nh, nl, iwah, iw1, iwbh
        integer(ip)  :: iw2, iw3, iwd, iww, iwu
        integer(ip)  :: irwk, icwk

        common_variables: associate( &
            npp => util%npp, &
            k => util%k, &
            nm => util%nm, &
            ncmplx=> util%ncmplx, &
            ik => util%ik, &
            cnv => util%cnv &
            )

            ! Test m and n for the proper form
            nm = n
            m2 = 2*m

            ! Check input arguments
            if (m < 5) then
                ierror = 1
                return
            else if (nm < 3) then
                ierror = 2
                return
            else if (idimy < m) then
                ierror = 3
                return
            else
                ierror = 0
            end if

            ! Compute workspace indices
            nh = n
            npp = np

            if (npp /= 0) nh = nh + 1

            ik = 4
            k = 3

            do
                if (nh <= ik) exit
                ik = 2*ik
                k = k + 1
            end do

            nl = ik
            ik = 2*ik
            nl = nl - 1
            iwah = (k - 2)*ik + k + 6

            select case (npp)
                case(0)
                    iwbh = iwah + 2*nm
                    iw1 = iwbh
                    nm = nm - 1
                case default
                    iw1 = iwah
                    iwbh = iw1 + nm
            end select

            iw2 = iw1 + m
            iw3 = iw2 + m
            iwd = iw3 + m
            iww = iwd + m
            iwu = iww + m

            select case (iflg)
                case (0)
                    ! Initialize solver

                    ! Set required workspace sizes
                    irwk = iw1 + 2*n
                    icwk = iw1 + 6*m

                    ! Allocate memory
                    call w%create(irwk, icwk)

                    ! Check if allocation was successful
                    if (ierror == 20) return

                    associate( &
                        rew => w%real_workspace, &
                        cxw => w%complex_workspace &
                        )
                        ! Compute roots of b polynomials
                        call util%cblktri_compute_roots_of_b_polynomials(ierror, an, bn, cn, rew, cxw, rew(iwah:), rew(iwbh:))
                    end associate
                case default
                    ! Solve system
                    associate( &
                        rew => w%real_workspace, &
                        cxw => w%complex_workspace &
                        )
                        select case (mp)
                            case (0)
                                call util%cblktri_lower_routine(nl, an, cn, m, am, bm, cm, &
                                    idimy, y, rew, cxw, &
                                    cxw(iw1), cxw(iw2), cxw(iw3), cxw(iwd), cxw(iww), &
                                    cxw(iwu), procp, cprocp)
                            case default
                                call util%cblktri_lower_routine(nl, an, cn, m, am, bm, cm, &
                                    idimy, y, rew, cxw, &
                                    cxw(iw1), cxw(iw2), cxw(iw3), cxw(iwd), cxw(iww), &
                                    cxw(iwu), proc, cproc)
                        end select
                    end associate
            end select
        end associate common_variables

    end subroutine cblktri

    ! Purpose:
    !
    ! Solves the linear system
    !
    ! Remarks:
    !
    ! b  contains the roots of all the b polynomials
    ! w1, w2, w3, wd, ww, wu  are all working arrays
    ! prdct is either procp or proc depending on whether the boundary
    ! conditions in the m direction are periodic or not
    ! cprdct is either cprocp or cproc which are called if some of the zeros
    ! of the b polynomials are complex.
    !
    subroutine cblktri_lower_routine(self, n, an, cn, m, am, bm, cm, idimy, y, b, bc, &
        w1, w2, w3, wd, ww, wu, prdct, cprdct)

        ! Dummy arguments
        class(ComplexGeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(in)    :: n
        integer(ip),      intent(in)    :: m
        integer(ip),      intent(in)    :: idimy
        real(wp),         intent(in)    :: an(:)
        real(wp),         intent(in)    :: cn(:)
        complex(wp),      intent(in)    :: am(:)
        complex(wp),      intent(in)    :: bm(:)
        complex(wp),      intent(in)    :: cm(:)
        real(wp),         intent(out)   :: b(*)
        complex(wp),      intent(inout) :: y(idimy,n)
        complex(wp),      intent(out)   :: bc(*)
        complex(wp),      intent(out)   :: w1(*)
        complex(wp),      intent(out)   :: w2(*)
        complex(wp),      intent(out)   :: w3(*)
        complex(wp),      intent(out)   :: wd(*)
        complex(wp),      intent(out)   :: ww(*)
        complex(wp),      intent(out)   :: wu(*)
        external                        :: prdct, cprdct

        ! Local variables
        integer(ip) :: kdo, l, ir, i2, i1, i3, i4, irm1, im2, nm2, im3, nm3
        integer(ip) :: im1, nm1, iif, i, ipi1, ipi2, ipi3, idxc, nc, idxa, na, ip2, np2
        integer(ip) :: ip1, np1, ip3, np3, iz, nz, izr, ll, ifd
        integer(ip) :: iip, np, imi1, imi2
        real(wp)    :: dummy_variable(1)

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv &
            )

            ! begin reduction phase
            kdo = k - 1
            do l = 1, kdo
                ir = l - 1
                i2 = 2**ir
                i1 = i2/2
                i3 = i2 + i1
                i4 = i2 + i2
                irm1 = ir - 1
                call self%cblktri_compute_index_b_coeff(i2, ir, im2, nm2)
                call self%cblktri_compute_index_b_coeff(i1, irm1, im3, nm3)
                call self%cblktri_compute_index_b_coeff(i3, irm1, im1, nm1)
                call prdct(nm2, b(im2), nm3, b(im3), nm1, b(im1), 0, dummy_variable, &
                    y(1,i2), w3, m, am, bm, cm, wd, ww, wu)
                iif = 2**k
                do i = i4, iif, i4
                    if (i > nm) cycle
                    ipi1 = i + i1
                    ipi2 = i + i2
                    ipi3 = i + i3
                    call self%cblktri_compute_index_c_coeff(i, ir, idxc, nc)
                    if (iif <= i) cycle
                    call self%cblktri_compute_index_a_coeff(i, ir, idxa, na)
                    call self%cblktri_compute_index_b_coeff(i - i1, irm1, im1, nm1)
                    call self%cblktri_compute_index_b_coeff(ipi2, ir, ip2, np2)
                    call self%cblktri_compute_index_b_coeff(ipi1, irm1, ip1, np1)
                    call self%cblktri_compute_index_b_coeff(ipi3, irm1, ip3, np3)
                    call prdct(nm1, b(im1), 0, dummy_variable, 0, dummy_variable, na, an(idxa), w3, &
                        w1, m, am, bm, cm, wd, ww, wu)
                    if (ipi2 > nm) then
                        w3(:m) = ZERO
                        w2(:m) = ZERO
                    else
                        call prdct(np2, b(ip2), np1, b(ip1), np3, b(ip3), 0, dummy_variable, &
                            y(1, ipi2), w3, m, am, bm, cm, wd, ww, wu)
                        call prdct(np1, b(ip1), 0, dummy_variable, 0, dummy_variable, nc, cn(idxc), w3, &
                            w2, m, am, bm, cm, wd, ww, wu)
                    end if
                    y(:m, i) = w1(:m) + w2(:m) + y(:m, i)
                end do
            end do

            if (npp == 0) then
                iif = 2**k
                i = iif/2
                i1 = i/2
                call self%cblktri_compute_index_b_coeff(i - i1, k - 2, im1, nm1)
                call self%cblktri_compute_index_b_coeff(i + i1, k - 2, ip1, np1)
                call self%cblktri_compute_index_b_coeff(i, k - 1, iz, nz)
                call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dummy_variable, &
                    y(1, i), w1, m, am, bm, cm, wd, ww, wu)

                izr = i
                w2(:m) = w1(:m)

                do ll = 2, k
                    l = k - ll + 1
                    ir = l - 1
                    i2 = 2**ir
                    i1 = i2/2
                    i = i2
                    call self%cblktri_compute_index_c_coeff(i, ir, idxc, nc)
                    call self%cblktri_compute_index_b_coeff(i, ir, iz, nz)
                    call self%cblktri_compute_index_b_coeff(i - i1, ir - 1, im1, nm1)
                    call self%cblktri_compute_index_b_coeff(i + i1, ir - 1, ip1, np1)
                    call prdct(np1, b(ip1), 0, dummy_variable, 0, dummy_variable, nc, cn(idxc), w1, &
                        w1, m, am, bm, cm, wd, ww, wu)
                    w1(:m) = y(:m, i) + w1(:m)
                    call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dummy_variable, w1 &
                        , w1, m, am, bm, cm, wd, ww, wu)
                end do

                loop_118: do ll = 2, k
                    l = k - ll + 1
                    ir = l - 1
                    i2 = 2**ir
                    i1 = i2/2
                    i4 = i2 + i2
                    ifd = iif - i2
                    do i = i2, ifd, i4

                        if (i - i2 /= izr) cycle

                        if (i > nm) cycle loop_118

                        call self%cblktri_compute_index_a_coeff(i, ir, idxa, na)
                        call self%cblktri_compute_index_b_coeff(i, ir, iz, nz)
                        call self%cblktri_compute_index_b_coeff(i - i1, ir - 1, im1, nm1)
                        call self%cblktri_compute_index_b_coeff(i + i1, ir - 1, ip1, np1)
                        call prdct(nm1, b(im1), 0, dummy_variable, 0, dummy_variable, na, an(idxa), w2 &
                            , w2, m, am, bm, cm, wd, ww, wu)
                        w2(:m) = y(:m, i) + w2(:m)
                        call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dummy_variable, &
                            w2, w2, m, am, bm, cm, wd, ww, wu)
                        izr = i
                        if (i == nm) exit  loop_118
                    end do
                end do loop_118

                y(:m, nm+1) = y(:m, nm+1) - cn(nm+1)*w1(:m) - an(nm+1)*w2(:m)

                call self%cblktri_compute_index_b_coeff(iif/2, k - 1, im1, nm1)
                call self%cblktri_compute_index_b_coeff(iif, k - 1, iip, np)

                select case (ncmplx)
                    case (0)
                        call prdct(nm + 1, b(iip), nm1, b(im1), 0, dummy_variable, 0, dummy_variable, &
                            y(1,nm+1), y(1, nm+1), m, am, bm, cm, wd, ww, wu)
                    case default
                        call cprdct(nm + 1, bc(iip), nm1, b(im1), 0, dummy_variable, 0, dummy_variable, &
                            y(1,nm+1), y(1, nm+1), m, am, bm, cm, w1, w3, ww)
                end select

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
                    call self%cblktri_compute_index_a_coeff(i, ir, idxa, na)
                    call self%cblktri_compute_index_b_coeff(i - i2, ir, im2, nm2)
                    call self%cblktri_compute_index_b_coeff(i - i2 - i1, ir - 1, im3, nm3)
                    call self%cblktri_compute_index_b_coeff(i - i1, ir - 1, im1, nm1)
                    call prdct(nm2, b(im2), nm3, b(im3), nm1, b(im1), 0, dummy_variable, &
                        w1, w1, m, am, bm, cm, wd, ww, wu)
                    call prdct(nm1, b(im1), 0, dummy_variable, 0, dummy_variable, na, an(idxa), w1, &
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
                    do i = i4, iif, i4
                        ipi1 = i + i1
                        ipi2 = i + i2
                        ipi3 = i + i3

                        if (ipi2 /= izr) then
                            if (i /= izr) cycle
                            cycle  loop_131
                        end if

                        call self%cblktri_compute_index_c_coeff(i, ir, idxc, nc)
                        call self%cblktri_compute_index_b_coeff(ipi2, ir, ip2, np2)
                        call self%cblktri_compute_index_b_coeff(ipi1, irm1, ip1, np1)
                        call self%cblktri_compute_index_b_coeff(ipi3, irm1, ip3, np3)
                        call prdct(np2, b(ip2), np1, b(ip1), np3, b(ip3), 0, dummy_variable, &
                            w2, w2, m, am, bm, cm, wd, ww, wu)
                        call prdct(np1, b(ip1), 0, dummy_variable, 0, dummy_variable, nc, cn(idxc), w2, &
                            w2, m, am, bm, cm, wd, ww, wu)
                        y(:m, i) = y(:m, i) - w2(:m)
                        izr = i
                        cycle  loop_131
                    end do
                end do loop_131
            end if

            ! begin back substitution phase
            do ll = 1, k
                l = k - ll + 1
                ir = l - 1
                irm1 = ir - 1
                i2 = 2**ir
                i1 = i2/2
                i4 = i2 + i2
                ifd = iif - i2
                do i = i2, ifd, i4
                    if (i > nm) cycle
                    imi1 = i - i1
                    imi2 = i - i2
                    ipi1 = i + i1
                    ipi2 = i + i2
                    call self%cblktri_compute_index_a_coeff(i, ir, idxa, na)
                    call self%cblktri_compute_index_c_coeff(i, ir, idxc, nc)
                    call self%cblktri_compute_index_b_coeff(i, ir, iz, nz)
                    call self%cblktri_compute_index_b_coeff(imi1, irm1, im1, nm1)
                    call self%cblktri_compute_index_b_coeff(ipi1, irm1, ip1, np1)

                    if (i <= i2) then
                        w1(:m) = ZERO
                    else
                        call prdct(nm1, b(im1), 0, dummy_variable, 0, dummy_variable, na, an(idxa), &
                            y(1,imi2), w1, m, am, bm, cm, wd, ww, wu)
                    end if

                    if (ipi2 > nm) then
                        w2(:m) = ZERO
                    else
                        call prdct(np1, b(ip1), 0, dummy_variable, 0, dummy_variable, nc, cn(idxc), y( &
                            1, ipi2), w2, m, am, bm, cm, wd, ww, wu)
                    end if
                    w1(:m) = y(:m, i) + w1(:m) + w2(:m)
                    call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dummy_variable, w1, &
                        y(1, i), m, am, bm, cm, wd, ww, wu)
                end do
            end do

        end associate common_variables

    end subroutine cblktri_lower_routine

    function cblktri_bsrh(self, xll, xrr, iz, c, a, bh, f, sgn) &
        result(return_value)

        ! Dummy arguments
        class(ComplexGeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip)                     :: iz
        real(wp),         intent(in)    :: xll
        real(wp),         intent(in)    :: xrr
        real(wp)                        :: c(*)
        real(wp)                        :: a(*)
        real(wp)                        :: bh(*)
        procedure(comf_interface)       :: f
        real(wp), intent(in)            :: sgn
        real(wp)                        :: return_value

        ! Local variables
        real(wp) :: r1, xl, xr, dx, x

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv &
            )

            xl = xll
            xr = xrr
            dx = HALF*abs(xr - xl)
            x = HALF*(xl + xr)
            r1 = sgn*f(x, iz, c, a, bh)

            if (r1 >= ZERO) then
                if (r1 == ZERO) then
                    return_value = HALF*(xl + xr)
                    return
                end if
                xr = x
            else
                xl = x
            end if

            dx = HALF * dx

            do
                if (dx <= cnv) exit
                x = HALF*(xl + xr)
                r1 = sgn*f(x, iz, c, a, bh)
                if (r1 >= ZERO) then
                    if (r1 == ZERO) then
                        return_value = HALF*(xl + xr)
                        return
                    end if
                    xr = x
                else
                    xl = x
                end if
                dx = HALF*dx
            end do

            return_value = HALF*(xl + xr)

        end associate common_variables

    end function cblktri_bsrh

    ! Purpose:
    !
    ! Computes the roots of the b polynomials using subroutine
    ! cblktri_tevls which is a modification the eispack program tqlrat.
    ! ierror is set to 4 if either cblktri_tevls fails or if a(j+1)*c(j) is
    ! less than zero for some j.  ah, bh are temporary work arrays.
    !
    subroutine cblktri_compute_roots_of_b_polynomials(self, ierror, an, bn, cn, b, bc, ah, bh)

        ! Dummy arguments
        class(ComplexGeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),                  intent(out)   :: ierror
        real(wp),                     intent(in)    :: an(:)
        real(wp),                     intent(in)    :: bn(:)
        real(wp),                     intent(in)    :: cn(:)
        real(wp), target, contiguous, intent(inout) :: b(:)
        real(wp),                     intent(inout) :: ah(:)
        real(wp),                     intent(inout) :: bh(:)
        complex(wp),                  intent(inout) :: bc(:)

        ! Local variables
        integer(ip) :: j, iif, kdo, l, ir, i2, i4
        integer(ip) :: ipl, ifd, i, ib, nb, js, jf
        integer(ip) :: ls, lh, nmp, l1, l2, j2, j1, n2m2
        real(wp)    :: bnorm, arg, d1, d2, d3

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv &
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

            cnv = MACHINE_EPSILON*bnorm
            iif = 2**k
            kdo = k - 1

            outer_loop: do l = 1, kdo

                ir = l - 1
                i2 = 2**ir
                i4 = i2 + i2
                ipl = i4 - 1
                ifd = iif - i4

                do i = i4, ifd, i4

                    call self%cblktri_compute_index_b_coeff(i, l, ib, nb)

                    if (nb <= 0) cycle outer_loop

                    js = i - ipl
                    jf = js + nb - 1
                    ls = 0
                    bh(:jf-js+1) = bn(js:jf)
                    ah(:jf-js+1) = b(js:jf)

                    call self%cblktri_tevls(bh(1:nb), ah(1:nb), ierror)

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

                call self%cblktri_tevls(ah(1:nb), bh(1:nb), ierror)

                if (ierror /= 0) then
                    ierror = 4
                    return
                end if

                call self%cblktri_compute_index_b_coeff(iif, k - 1, j2, lh)
                call self%cblktri_compute_index_b_coeff(iif/2, k - 1, j1, lh)

                j2 = j2 + 1
                lh = j2
                n2m2 = j2 + 2*nm - 2

                iteration: do

                    d1 = abs(b(j1)-b(j2-1))
                    d2 = abs(b(j1)-b(j2))
                    d3 = abs(b(j1)-b(j2+1))

                    if (d1 <= d2 .or. d3 <= d2) then
                        b(lh) = b(j2)
                        j2 = j2 + 1
                        lh = lh + 1
                        if (j2 <= n2m2) cycle iteration
                    else
                        j2 = j2 + 1
                        j1 = j1 + 1
                        if (j2 <= n2m2) cycle iteration
                    end if
                    exit iteration
                end do iteration

                b(lh) = b(n2m2+1)

                call self%cblktri_compute_index_b_coeff(iif, k - 1, j1, j2)

                j2 = j1 + nmp + nmp

                ! Compute eigenvalues of the periodic tridiagonal
                block
                    complex(wp) :: cbp_arg(1)
                    real(wp)    :: bp_arg(1)
                    real(wp), contiguous, pointer :: bh_arg(:) => null()

                    ! Associate arguments
                    cbp_arg = cmplx(b(j1), kind=wp)
                    bp_arg = real(bc(j1), kind=wp)
                    bh_arg(1:) => b(j2:)

                    ! Call solver
                    call self%cblktri_compute_eigenvalues(nm + 1, ierror, an, cn, cbp_arg, bp_arg, bh_arg)

                    ! Terminate association
                    nullify(bh_arg)
                end block
            end if

        end associate common_variables

    end subroutine cblktri_compute_roots_of_b_polynomials

    ! Purpose:
    !
    ! Applies a sequence of matrix operations to the vector x and
    ! stores the result in y
    ! aa   array containing scalar multipliers of the vector x
    ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
    ! bd, bm1, bm2 are arrays containing roots of certian b polynomials
    ! na is the length of the array aa
    ! x, y the matrix operations are applied to x and the result is y
    ! a, b, c  are arrays which contain the tridiagonal matrix
    ! m  is the order of the matrix
    ! d, w are work arrays
    ! isgn  determines whether or not a change in sign is made
    !
    pure subroutine cproc(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, w, yy)

        ! Dummy arguments
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

        ! Local variables
        integer(ip) :: j, mm, id, m1, m2, ia, iflg, k
        real(wp)    :: rt
        complex(wp) :: crt, den, y1, y2

        y(:m) = x(:m)
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

                ! begin solution to system
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

                do j = 2, m
                    y(j) = w(j) - d(j)*y(j-1)
                end do

            end if

            if (.not.(m1 <= 0 .and. m2 <= 0)) then
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

                ! scalar multiplication
                y(:m) = rt*y(:m)
            end if

            if (iflg <= 0) exit main_loop

        end do main_loop

    end subroutine cproc

    ! Purpose:
    !
    ! cprocp applies a sequence of matrix operations to the vector x and
    ! stores the result in y
    !
    ! bd, bm1, bm2 are arrays containing roots of certian b polynomials
    ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
    ! aa   array containing scalar multipliers of the vector x
    ! na is the length of the array aa
    ! x, y the matrix operations are applied to x and the result is y
    ! a, b, c  are arrays which contain the tridiagonal matrix
    ! m  is the order of the matrix
    ! d, u are work arrays
    ! isgn  determines whether or not a change in sign is made
    !
    pure subroutine cprocp(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, u, yy)

        ! Dummy arguments
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

        ! Local variables
        integer(ip) :: j, mm, mm2, id, m1, m2, ia, iflg, k
        real(wp)    :: rt
        complex(wp) :: v, den, bh, ym, am, y1, y2, yh, crt

        y(:m) = x(:m)
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

                !  begin solution to system
                bh = b(m) - crt
                ym = y(m)
                den = b(1) - crt
                d(1) = c(1)/den
                u(1) = a(1)/den
                y(1) = y(1)/den
                v = c(m)
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

            if (.not.(m1 <= 0 .and. m2 <= 0)) then
                if (m1 <= 0) then
                    rt = bm2(m2)
                    m2 = m2 - 1
                else
                    if (m2 <= 0) then
                        rt = bm1(m1)
                        m1 = m1 - 1
                    else
                        if (abs(bm1(m1)) > abs(bm2(m2))) then
                            rt = bm1(m1)
                            m1 = m1 - 1
                        else
                            rt = bm2(m2)
                            m2 = m2 - 1

                        ! matrix multiplication
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

                ! scalar multiplication
                y(:m) = rt*y(:m)
            end if

            if (iflg <= 0) exit main_loop

        end do main_loop

    end subroutine cprocp

    subroutine cblktri_compute_index_a_coeff(self, i, ir, idxa, na)

        ! Dummy arguments
        class(ComplexGeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(in)    :: i
        integer(ip),      intent(in)    :: ir
        integer(ip),      intent(out)   :: idxa
        integer(ip),      intent(out)   :: na

        associate( nm => self%nm )
            na = 2**ir
            idxa = i - na + 1
            if (i > nm) na = 0
        end associate

    end subroutine cblktri_compute_index_a_coeff

    subroutine cblktri_compute_index_b_coeff(self, i, ir, idx, idp)

        ! Dummy arguments
        class(ComplexGeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip), intent(in) :: i
        integer(ip), intent(in) :: ir
        integer(ip), intent(out) :: idx
        integer(ip), intent(out) :: idp

        ! Local variables
        integer(ip) :: izh, id, ipl

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv &
            )

            ! b(idx) is the location of the first root of the b(i, ir) polynomial
            idp = 0
            if (ir >= 0) then
                if (ir <= 0) then

                    if (i > nm) return

                    idx = i
                    idp = 1
                    return
                end if
                izh = 2**ir
                id = i - 2*izh
                idx = 2*id + (ir - 1)*ik + ir + (ik - i)/izh + 4
                ipl = izh - 1
                idp = 2*izh - 1

                if (i - ipl > nm) then
                    idp = 0
                    return
                end if

                if (i + ipl > nm) idp = nm + ipl - i + 1

            end if

        end associate common_variables

    end subroutine cblktri_compute_index_b_coeff

    subroutine cblktri_compute_index_c_coeff(self, i, ir, idxc, nc)

        ! Dummy arguments
        class(ComplexGeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(in)    :: i
        integer(ip),      intent(in)    :: ir
        integer(ip),      intent(out)   :: idxc
        integer(ip),      intent(out)   :: nc

        associate( nm => self%nm )
            nc = 2**ir
            idxc = i
            if (idxc + nc - 1 > nm) nc = 0
        end associate

    end subroutine cblktri_compute_index_c_coeff

    ! Purpose:
    !
    ! Computes the eigenvalues of the periodic tridiagonal
    ! matrix with coefficients an, bn, cn
    !
    ! n is the order of the bh and bp polynomials
    ! on output bp contains the eigenvalues
    ! cbp is the same as bp except type complex
    ! bh is used to temporarily store the roots of the b hat polynomial
    ! which enters through bp
    subroutine cblktri_compute_eigenvalues(self, n, ierror, a, c, cbp, bp, bh)

        ! Dummy arguments
        class(ComplexGeneralizedCyclicReductionUtility), intent(inout) :: self
        integer(ip),      intent(in)    :: n
        integer(ip),      intent(out)   :: ierror
        real(wp),         intent(in)    :: a(:)
        real(wp),         intent(in)    :: c(:)
        complex(wp)       :: cbp(:)
        real(wp)          :: bp(:)
        real(wp)          :: bh(:)

        ! Local variables
        integer(ip) :: iz, izm, izm2, j, nt, modiz
        integer(ip) :: iis, iif, ig, it, icv, i3, i2, nhalf
        real(wp)    :: r4, r5, r6, scnv, xl, db, sgn, xr, xm, psg
        real(wp)    :: temp
        complex(wp) :: cx, fsg, hsg, dd, f, fp, fpp, cdis, r1, r2, r3

        common_variables: associate( &
            npp => self%npp, &
            k => self%k, &
            nm => self%nm, &
            ncmplx=> self%ncmplx, &
            ik => self%ik, &
            cnv => self%cnv &
            )

            scnv = sqrt(cnv)
            iz = n
            izm = iz - 1
            izm2 = iz - 2

            main_block: block

                if (bp(n) <= bp(1)) then
                    if (bp(n) == bp(1)) exit main_block
                    bh(:n) = bp(n:1:(-1))
                else
                    bh(:n) = bp(:n)
                end if

                ncmplx = 0
                modiz = mod(iz, 2)
                iis = 1

                block_110: block
                    if (modiz /= 0) then
                        if (a(1) < ZERO) exit block_110
                        if (a(1) == ZERO) exit main_block
                    end if
                    xl = bh(1)
                    db = bh(3) - bh(1)
                    xl = xl - db
                    r4 = psgf(xl, iz, c, a, bh)
                    do while (r4 <= ZERO)
                        xl = xl - db
                        r4 = psgf(xl, iz, c, a, bh)
                    end do
                    sgn = -ONE
                    temp = self%cblktri_bsrh(xl, bh(1), iz, c, a, bh, psgf, sgn)
                    cbp(1) = cmplx(temp, ZERO, kind=wp)
                    iis = 2
                end block block_110

                iif = iz - 1

                block_115: block
                    if (modiz /= 0) then
                        if (a(1) > ZERO) exit block_115
                        if (a(1) == ZERO) exit main_block
                    end if
                    xr = bh(iz)
                    db = bh(iz) - bh(iz-2)
                    xr = xr + db
                    r5 = psgf(xr, iz, c, a, bh)
                    do while (r5 < ZERO)
                        xr = xr + db
                        r5 = psgf(xr, iz, c, a, bh)
                    end do
                    sgn = ONE
                    temp = self%cblktri_bsrh(bh(iz), xr, iz, c, a, bh, psgf, sgn)
                    cbp(iz) = cmplx(temp, ZERO, kind=wp)
                    iif = iz - 2
                end block block_115

                main_loop: do ig = iis, iif, 2
                    xl = bh(ig)
                    xr = bh(ig+1)
                    sgn = -1.
                    xm = self%cblktri_bsrh(xl, xr, iz, c, a, bh, ppspf, sgn)
                    psg = psgf(xm, iz, c, a, bh)

                    if_block: block
                        if (abs(psg) > MACHINE_EPSILON) then
                            r6 = psg*ppsgf(xm, iz, c, a, bh)
                            if (r6 > ZERO) exit if_block
                            if (r6 /= ZERO) then
                                sgn = ONE
                                cbp(ig) = cmplx(self%cblktri_bsrh(bh(ig), xm, iz, c, a, bh, psgf, sgn), ZERO, kind=wp)
                                sgn = -ONE
                                cbp(ig+1) = cmplx(self%cblktri_bsrh(xm, bh(ig+1), iz, c, a, bh, psgf, sgn), ZERO, kind=wp)
                                cycle main_loop

                            !     case of a multiple zero
                            end if
                        end if
                        cbp(ig) = cmplx(xm, ZERO, kind=wp)
                        cbp(ig+1) = cmplx(xm, ZERO, kind=wp)
                        cycle main_loop

                    !     case of a complex zero
                    end block if_block

                    it = 0
                    icv = 0
                    cx = cmplx(xm, ZERO, kind=wp)

                    loop_120: do
                        fsg = cmplx(ONE, ZERO, kind=wp)
                        hsg = cmplx(ONE, ZERO, kind=wp)
                        fp = ZERO
                        fpp = ZERO

                        do j = 1, iz
                            dd = ONE /(cx - bh(j))
                            fsg = fsg*a(j)*dd
                            hsg = hsg*c(j)*dd
                            fp = fp + dd
                            fpp = fpp - dd*dd
                        end do

                        if (modiz == 0) then
                            f = cmplx(ONE, ZERO, kind=wp) - fsg - hsg
                        else
                            f = cmplx(ONE, ZERO, kind=wp) + fsg + hsg
                        end if

                        i3 = 0

                        if (abs(fp) > ZERO) then
                            i3 = 1
                            r3 = -f/fp
                        end if

                        i2 = 0

                        if (abs(fpp) > ZERO) then
                            i2 = 1
                            cdis = sqrt(fp**2 - TWO*f*fpp)
                            r1 = cdis - fp
                            r2 = (-fp) - cdis
                            if (abs(r1) - abs(r2) > ZERO) then
                                r1 = r1/fpp
                            else
                                r1 = r2/fpp
                            end if
                            r2 = TWO*f/fpp/r1
                            if (abs(r2) < abs(r1)) r1 = r2
                            if (i3 > 0) then
                                if (abs(r3) < abs(r1)) r1 = r3
                            end if
                        else
                            r1 = r3
                        end if

                        cx = cx + r1
                        it = it + 1
                        if (it > 50) exit main_block
                        if (abs(r1) > scnv) cycle loop_120
                        if (icv > 0) exit loop_120
                        icv = 1
                    end do loop_120

                    cbp(ig) = cx
                    cbp(ig+1) = conjg(cx)
                end do main_loop

                if (abs(cbp(n)) - abs(cbp(1)) <= ZERO) then
                    if (abs(cbp(n)) - abs(cbp(1)) == ZERO) exit main_block
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
                    if (aimag(cbp(j)) /= ZERO) return
                end do

                ncmplx = 0

                do j = 2, iz
                    bp(j) = real(cbp(j), kind=wp)
                end do

                return
            end block main_block
            !
            ! Procedure failed
            !
            ierror = 4

        end associate common_variables

    end subroutine cblktri_compute_eigenvalues

    ! Purpose:
    !
    ! proc applies a sequence of matrix operations to the vector x and
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
    pure subroutine proc(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, w, u)

        ! Dummy arguments
        integer(ip), intent(in)  :: nd
        integer(ip), intent(in)  :: nm1
        integer(ip), intent(in)  :: nm2
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: bd(nd)
        real(wp),    intent(in)  :: bm1(nm1)
        real(wp),    intent(in)  :: bm2(nm2)
        real(wp),    intent(in)  :: aa(na)
        complex(wp), intent(in)  :: x(m)
        complex(wp), intent(out) :: y(m)
        complex(wp), intent(in)  :: a(m)
        complex(wp), intent(in)  :: b(m)
        complex(wp), intent(in)  :: c(m)
        complex(wp), intent(out) :: d(m)
        complex(wp), intent(out) :: w(m)
        complex(wp), intent(out) :: u(m)

        ! Local variables
        integer(ip) :: j, mm, id, ibr, m1, m2, ia, k
        real(wp)    :: rt
        complex(wp) :: den

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

                ! scalar multiplication
                y = rt*w
            end if

            if (id <= 0) return

            rt = bd(id)
            id = id - 1

            if (id == 0) ibr = 1

            ! begin solution to system
            d(m) = a(m)/(b(m)-rt)
            w(m) = y(m)/(b(m)-rt)

            do j = 2, mm
                k = m - j
                den = b(k+1) - rt - c(k+1)*d(k+2)
                d(k+1) = a(k+1)/den
                w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
            end do

            den = b(1) - rt - c(1)*d(2)
            w(1) = cmplx(ONE, ZERO, kind=wp)

            if (abs(den) /= ZERO) then
                w(1) = (y(1)-c(1)*w(2))/den
            end if

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
                    if (ibr <= 0 .and. abs(bm1(m1)-bd(id)) < abs(bm1(m1)-rt)) then
                        y = w
                        ibr = 1
                        cycle main_loop
                    end if
                end if
                rt = rt - bm1(m1)
                m1 = m1 - 1
                y = y + rt*w
                cycle main_loop
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

    end subroutine proc

    ! Purpose:
    !
    ! procp applies a sequence of matrix operations to the vector x and
    ! stores the result in y  periodic boundary conditions
    !
    ! bd, bm1, bm2 are arrays containing roots of certian b polynomials
    ! nd, nm1, nm2 are the lengths of the arrays bd, bm1, bm2 respectively
    ! aa   array containing scalar multipliers of the vector x
    ! na is the length of the array aa
    ! x, y  the matrix operations are applied to x and the result is y
    ! a, b, c  are arrays which contain the tridiagonal matrix
    ! m  is the order of the matrix
    ! d, u, w are working arrays
    ! is  determines whether or not a change in sign is made
    !
    pure subroutine procp(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, u, w)

        ! Dummy arguments
        integer(ip), intent(in)  :: nd
        integer(ip), intent(in)  :: nm1
        integer(ip), intent(in)  :: nm2
        integer(ip), intent(in)  :: na
        integer(ip), intent(in)  :: m
        real(wp),    intent(in)  :: bd(nd)
        real(wp),    intent(in)  :: bm1(nm1)
        real(wp),    intent(in)  :: bm2(nm2)
        real(wp),    intent(in)  :: aa(na)
        complex(wp), intent(in)  :: x(m)
        complex(wp), intent(out) :: y(m)
        complex(wp), intent(in)  :: a(m)
        complex(wp), intent(in)  :: b(m)
        complex(wp), intent(in)  :: c(m)
        complex(wp), intent(out) :: d(m)
        complex(wp), intent(out) :: u(m)
        complex(wp), intent(out) :: w(m)

        ! Local variables
        integer(ip) :: j, mm, mm2, id, ibr, m1, m2, ia, k
        real(wp)    :: rt
        complex(wp) :: den, ym, v, bh, am

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
                if (nd == 0) then
                    rt = -aa(ia)
                else
                    rt = aa(ia)
                end if
                ia = ia - 1
                y = rt*w
            end if

            if (id <= 0) return

            rt = bd(id)
            id = id - 1

            if (id == 0) ibr = 1

            ! begin solution to system
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

            if (abs(den) /= ZERO) then
                w(m) = (ym - am*w(m-1))/den
            else
                w(m) = cmplx(ONE, ZERO, kind=wp)
            end if

            w(m-1) = w(m-1) - d(m-1)*w(m)

            do j = 2, mm
                k = m - j
                w(k) = w(k) - d(k)*w(k+1) - u(k)*w(m)
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
                    if (ibr <= 0 .and. abs(bm1(m1)-bd(id)) < abs(bm1(m1)-rt)) then
                        y = w
                        ibr = 1
                        cycle main_loop
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

    end subroutine procp

    ! Purpose:
    !
    !     Finds the eigenvalues of a symmetric
    !     tridiagonal matrix by the rational ql method.
    !     This subroutine is a modification of the eispack subroutine tqlrat
    !     algorithm 464, comm. acm 16, 689(1973) by reinsch.
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
    !     eps is a machine dependent parameter specifying
    !     the relative precision of floating point arithmetic.
    !
    !
    subroutine cblktri_tevls(self, diagonal, subdiagonal, error_flag)

        ! Dummy arguments
        class(ComplexGeneralizedCyclicReductionUtility), intent(inout) :: self
        real(wp),         intent(inout) :: diagonal(:)
        real(wp),         intent(inout) :: subdiagonal(:)
        integer(ip),      intent(out)   :: error_flag

        ! Local variables
        integer(ip) :: i, j, l, m, ii, l1, mml, nhalf, ntop
        real(wp)    :: b, c, f, g, h, p, r, s, dhold

        associate( &
            n => size(diagonal), &
            d => diagonal, &
            e2 => subdiagonal &
            )

            error_flag = 0
            if (n /= 1) then

                e2(:n-1) = e2(2:n)*e2(2:n)
                f = ZERO
                b = ZERO
                e2(n) = ZERO

                main_loop: do l = 1, n
                    j = 0
                    h = MACHINE_EPSILON*(abs(d(l))+sqrt(e2(l)))

                    if (b <= h) then
                        b = h
                        c = b*b
                    end if

                    ! look for small squared sub-diagonal element
                    do m = l, n
                        if (e2(m) > c) cycle
                        exit

                    ! 2(n) is always zero, so there is no exit
                    ! through the bottom of the loop
                    end do

                    if_block: block
                        if (m /= l) then
                            loop_105: do
                                if (j == 30) then

                                    ! set error  no convergence to an
                                    ! eigenvalue after 30 iterations
                                    error_flag = l
                                    return
                                end if

                                j = j + 1

                                ! form shift
                                l1 = l + 1
                                s = sqrt(e2(l))
                                g = d(l)
                                p = (d(l1)-g)/(TWO*s)
                                r = sqrt(p**2 + ONE)
                                d(l) = s/(p + sign(r, p))
                                h = g - d(l)
                                d(l1:n) = d(l1:n) - h
                                f = f + h

                                ! rational ql transformation
                                g = d(m)

                                if (g == ZERO) g = b

                                h = g
                                s = ZERO
                                mml = m - l

                                ! for i=m-1 step -1 until l do
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

                                e2(l) = s*g
                                d(l) = h

                                !  guard against underflowed h
                                if (h == ZERO .or. abs(e2(l)) <= abs(c/h)) exit if_block

                                e2(l) = h*e2(l)

                                if (e2(l) == ZERO) exit loop_105
                            end do loop_105
                        end if
                    end block if_block

                    p = d(l) + f

                    ! order eigenvalues
                    if (l /= 1) then

                        ! for i=l step -1 until 2 do
                        do ii = 2, l
                            i = l + 2 - ii
                            if (p >= d(i-1)) then
                                d(i) = p
                                cycle main_loop
                            end if
                            d(i) = d(i-1)
                        end do
                    end if
                    i = 1
                    d(i) = p
                end do main_loop

                if (abs(d(n)) >= abs(d(1))) return

                nhalf = n/2

                do i = 1, nhalf
                    ntop = n - i
                    dhold = d(i)
                    d(i) = d(ntop+1)
                    d(ntop+1) = dhold
                end do
            end if
        end associate

    end subroutine cblktri_tevls

end module complex_block_tridiagonal_linear_systems_solver
