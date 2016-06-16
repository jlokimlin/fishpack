!
!     file cblktri.f90
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
!     subroutine cblktr(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, idimy, y,
!                       ierror)
!
!
! DIMENSION OF           an(n), bn(n), cn(n), am(m), bm(m), cm(m), y(idimy, n)
! ARGUMENTS
!
! LATEST REVISION        May 2016
!
! PURPOSE                cblktr solves a system of linear equations
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
! USAGE                  call cblktr(iflg, np, n, an, bn, cn, mp, m, am, bm,
!                                     cm, idimy, y, ierror, w)
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
!                          in the program calling cblktr.
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
!                          a fortran 90 derived type (FishpackWorkspace) variable
!                          that must be declared by the user.  the first
!                          two declarative statements in the user program
!                          calling cblktri must be:
!
!                               use type_FishpackWorkspace
!                               type (FishpackWorkspace) :: w
!
!                          the first statement makes the fishpack module
!                          defined in the file "type_FishpackWorkspace.f90" available to the
!                          user program calling cblktri.  the second statement
!                          declares a derived type variable (defined in
!                          the module "type_FishpackWorkspace.f90") which is used internally
!                          in cblktri to dynamically allocate real and complex
!                          work space used in solution.  an error flag
!                          (ierror = 20) is set if the required work space
!                          allocation fails (for example if n, m are too large)
!                          real and complex values are set in the components
!                          of w on a initial (iflg=0) call to cblktri.  these
!                          must be preserved on non-initial calls (iflg=1)
!                          to cblktri.  this eliminates redundant calculations
!                          and saves compute time.
!               ****       important!  the user program calling cblktri should
!                          include the statement:
!
!                               call w%destroy()
!
!                          after the final approximation is generated by
!                          cblktri.  the will deallocate the real and complex
!                          work space of w.  failure to include this statement
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
!                        = 4  cblktr failed while computing results
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
!                               complex work space in the derived type
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
! I/O                    None
!
! PRECISION              64-bit double precision
!
! REQUIRED LIBRARY       comf.f90, type_FishpackWorkspace.f90
! FILES
!
! STANDARD               Fortran 2008
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
module module_cblktri

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use module_comf, only: &
        psgf, &
        ppspf, &
        ppsgf

    use type_FishpackWorkspace, only: &
        Fish => FishpackWorkspace

    ! Explicit typing only
    implicit None

    ! Everything is private unless stated otherwise
    private
    public :: cblktri

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    integer (ip), save :: npp, k, nm, ncmplx, ik
    real (wp),    save :: cnv
    !---------------------------------------------------------------------------------


contains


    subroutine cblktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, &
        idimy, y, ierror, w)
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
        real (wp),    intent (in)     :: an(:)
        real (wp),    intent (in)     :: bn(:)
        real (wp),    intent (in)     :: cn(:)
        complex (wp), intent (in)     :: am(:)
        complex (wp), intent (in)     :: bm(:)
        complex (wp), intent (in)     :: cm(:)
        complex (wp), intent (in out) :: y(:,:)
        class (Fish), intent (in out) :: w
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: m2, nh, nl, iwah, iw1, iwbh
        integer (ip) :: iw2, iw3, iwd, iww, iwu, irwk, icwk
        !-----------------------------------------------
        !
        ! test m and n for the proper form
        !
        nm = n
        m2 = 2*m

        !
        !==> Check validity of input arguments
        !
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

        !
        !==> Compute workspace indices
        !
        nh = n
        npp = np

        if (npp /= 0) nh = nh + 1

        ik = 4
        k = 3

        do while (nh > ik)
            ik = 2*ik
            k = k + 1
        end do

        nl = ik
        ik = 2*ik
        nl = nl - 1
        iwah = (k - 2)*ik + k + 6

        if (npp /= 0) then
            iw1 = iwah
            iwbh = iw1 + nm
        else
            iwbh = iwah + 2*nm
            iw1 = iwbh
            nm = nm - 1
        end if

        iw2 = iw1 + m
        iw3 = iw2 + m
        iwd = iw3 + m
        iww = iwd + m
        iwu = iww + m

        select case (iflg)
            !
            !==> Initialize solver
            !
            case (0)
                ! Set required workspace sizes
                irwk = iw1 + 2*n
                icwk = iw1 + 6*m
                !
                !==> Allocate memory
                !
                call w%create(irwk, icwk, ierror)

                ! Check if allocation was successful
                if (ierror == 20) return

                associate( &
                    rew => w%real_workspace, &
                    cxw => w%complex_workspace &
                    )
                    !
                    !==> Computes roots of b polynomials
                    !
                    call ccompb(nl, ierror, an, bn, cn, rew,cxw, rew(iwah), rew(iwbh))

                end associate
            case default
                !
                !==> Solve system
                !
                associate( &
                    rew => w%real_workspace, &
                    cxw => w%complex_workspace &
                    )

                    select case (mp)
                        case (0)
                            call cblkt1(nl, an, bn, cn, m, am, bm, cm, &
                                idimy, y, rew, cxw, &
                                cxw(iw1), cxw(iw2), cxw(iw3), cxw(iwd), cxw(iww), &
                                cxw(iwu), procp, cprocp)
                        case default
                            call cblkt1(nl, an, bn, cn, m, am, bm, cm, &
                                idimy, y, rew, cxw, &
                                cxw(iw1), cxw(iw2), cxw(iw3), cxw(iwd), cxw(iww), &
                                cxw(iwu), proc, cproc)
                    end select
                end associate
        end select

    contains

        subroutine cblkt1(n, an, bn, cn, m, am, bm, cm, idimy, y, b, bc, &
            w1, w2, w3, wd, ww, wu, prdct, cprdct)
            !
            ! Purpose:
            !
            ! cblkt1 solves the linear system
            !
            ! b  contains the roots of all the b polynomials
            ! w1, w2, w3, wd, ww, wu  are all working arrays
            ! prdct is either procp or proc depending on whether the boundary
            ! conditions in the m direction are periodic or not
            ! cprdct is either cprocp or cproc which are called if some of the zeros
            ! of the b polynomials are complex.
            !
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in)     :: n
            integer (ip), intent (in)     :: m
            integer (ip), intent (in)     :: idimy
            real (wp),    intent (in)     :: an(:)
            real (wp),    intent (in)     :: bn(:)
            real (wp),    intent (in)     :: cn(:)
            complex (wp), intent (in)     :: am(:)
            complex (wp), intent (in)     :: bm(:)
            complex (wp), intent (in)     :: cm(:)
            real (wp),    intent (out)    :: b(*)
            complex (wp), intent (in out) :: y(idimy, *)
            complex (wp), intent (out)    :: bc(*)
            complex (wp), intent (out)    :: w1(*)
            complex (wp), intent (out)    :: w2(*)
            complex (wp), intent (out)    :: w3(*)
            complex (wp), intent (out)    :: wd(*)
            complex (wp), intent (out)    :: ww(*)
            complex (wp), intent (out)    :: wu(*)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip) :: kdo, l, ir, i2, i1, i3, i4, irm1, im2, nm2, im3, nm3
            integer (ip) :: im1, nm1, iif, i, ipi1, ipi2, ipi3, idxc, nc, idxa, na, ip2, np2
            integer (ip) :: ip1, np1, ip3, np3, iz, nz, izr, ll, ifd
            integer (ip) :: iip, np, imi1, imi2
            real (wp)    :: dum
            !-----------------------------------------------

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
                call cindxb(i2, ir, im2, nm2)
                call cindxb(i1, irm1, im3, nm3)
                call cindxb(i3, irm1, im1, nm1)
                call prdct(nm2, b(im2), nm3, b(im3), nm1, b(im1), 0, dum, &
                    y(1,i2), w3, m, am, bm, cm, wd, ww, wu)
                iif = 2**k
                do i = i4, iif, i4
                    if (i > nm) cycle
                    ipi1 = i + i1
                    ipi2 = i + i2
                    ipi3 = i + i3
                    call cindxc (i, ir, idxc, nc)
                    if (i >= iif) cycle
                    call cindxa(i, ir, idxa, na)
                    call cindxb(i - i1, irm1, im1, nm1)
                    call cindxb(ipi2, ir, ip2, np2)
                    call cindxb(ipi1, irm1, ip1, np1)
                    call cindxb(ipi3, irm1, ip3, np3)
                    call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), w3, &
                        w1, m, am, bm, cm, wd, ww, wu)
                    if (ipi2 > nm) then
                        w3(:m) = 0.0_wp
                        w2(:m) = 0.0_wp
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
                call cindxb(i - i1, k - 2, im1, nm1)
                call cindxb(i + i1, k - 2, ip1, np1)
                call cindxb(i, k - 1, iz, nz)
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
                    call cindxc(i, ir, idxc, nc)
                    call cindxb(i, ir, iz, nz)
                    call cindxb(i - i1, ir - 1, im1, nm1)
                    call cindxb(i + i1, ir - 1, ip1, np1)
                    call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), w1, &
                        w1, m, am, bm, cm, wd, ww, wu)
                    w1(:m) = y(:m, i) + w1(:m)
                    call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dum, w1 &
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

                        call cindxa(i, ir, idxa, na)
                        call cindxb(i, ir, iz, nz)
                        call cindxb(i - i1, ir - 1, im1, nm1)
                        call cindxb(i + i1, ir - 1, ip1, np1)
                        call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), w2 &
                            , w2, m, am, bm, cm, wd, ww, wu)
                        w2(:m) = y(:m, i) + w2(:m)
                        call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dum, &
                            w2, w2, m, am, bm, cm, wd, ww, wu)
                        izr = i
                        if (i == nm) exit  loop_118
                    end do
                end do loop_118

                y(:m, nm+1) = y(:m, nm+1) - cn(nm+1)*w1(:m) - an(nm+1)*w2(:m)

                call cindxb(iif/2, k - 1, im1, nm1)
                call cindxb(iif, k - 1, iip, np)

                select case (ncmplx)
                    case (0)
                        call prdct(nm + 1, b(iip), nm1, b(im1), 0, dum, 0, dum, &
                            y(1,nm+1), y(1, nm+1), m, am, bm, cm, wd, ww, wu)
                    case default
                        call cprdct(nm + 1, bc(iip), nm1, b(im1), 0, dum, 0, dum, &
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
                    call cindxa(i, ir, idxa, na)
                    call cindxb(i - i2, ir, im2, nm2)
                    call cindxb(i - i2 - i1, ir - 1, im3, nm3)
                    call cindxb(i - i1, ir - 1, im1, nm1)
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
                    do i = i4, iif, i4
                        ipi1 = i + i1
                        ipi2 = i + i2
                        ipi3 = i + i3

                        if (ipi2 /= izr) then
                            if (i /= izr) cycle
                            cycle  loop_131
                        end if

                        call cindxc (i, ir, idxc, nc)
                        call cindxb(ipi2, ir, ip2, np2)
                        call cindxb(ipi1, irm1, ip1, np1)
                        call cindxb(ipi3, irm1, ip3, np3)
                        call prdct(np2, b(ip2), np1, b(ip1), np3, b(ip3), 0, dum, &
                            w2, w2, m, am, bm, cm, wd, ww, wu)
                        call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), w2, &
                            w2, m, am, bm, cm, wd, ww, wu)
                        y(:m, i) = y(:m, i) - w2(:m)
                        izr = i
                        cycle  loop_131
                    end do
                end do loop_131
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
                ifd = iif - i2
                do i = i2, ifd, i4
                    if (i > nm) cycle
                    imi1 = i - i1
                    imi2 = i - i2
                    ipi1 = i + i1
                    ipi2 = i + i2
                    call cindxa(i, ir, idxa, na)
                    call cindxc(i, ir, idxc, nc)
                    call cindxb(i, ir, iz, nz)
                    call cindxb(imi1, irm1, im1, nm1)
                    call cindxb(ipi1, irm1, ip1, np1)

                    if (i <= i2) then
                        w1(:m) = 0.0_wp
                    else
                        call prdct(nm1, b(im1), 0, dum, 0, dum, na, an(idxa), &
                            y(1,imi2), w1, m, am, bm, cm, wd, ww, wu)
                    end if

                    if (ipi2 > nm) then
                        w2(:m) = 0.0_wp
                    else
                        call prdct(np1, b(ip1), 0, dum, 0, dum, nc, cn(idxc), y( &
                            1, ipi2), w2, m, am, bm, cm, wd, ww, wu)
                    end if
                    w1(:m) = y(:m, i) + w1(:m) + w2(:m)
                    call prdct(nz, b(iz), nm1, b(im1), np1, b(ip1), 0, dum, w1, &
                        y(1, i), m, am, bm, cm, wd, ww, wu)
                end do
            end do

        end subroutine cblkt1



        function cbsrh(xll, xrr, iz, c, a, bh, f, sgn) result(return_value)
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip)           :: iz
            real (wp), intent (in) :: xll
            real (wp), intent (in) :: xrr
            real (wp)              :: f
            real (wp), intent (in) :: sgn
            real (wp)              :: c(*)
            real (wp)              :: a(*)
            real (wp)              :: bh(*)
            real (wp)              :: return_value
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            real (wp) :: r1, xl, xr, dx, x
            !-----------------------------------------------

            xl = xll
            xr = xrr
            dx = 0.5_wp*abs(xr - xl)
            x = 0.5_wp*(xl + xr)
            r1 = sgn*f(x, iz, c, a, bh)

            if (r1 >= 0.0_wp) then
                if (r1 == 0.0_wp) then
                    return_value = 0.5_wp*(xl + xr)
                    return
                end if
                xr = x
            else
                xl = x
            end if

            dx = 0.5_wp * dx

            do while (dx > cnv)

                x = 0.5_wp*(xl + xr)
                r1 = sgn*f(x, iz, c, a, bh)

                if (r1 >= 0.0_wp) then
                    if (r1 == 0.0_wp) then
                        return_value = 0.5_wp*(xl + xr)
                        return
                    end if
                    xr = x
                else
                    xl = x
                end if
                dx = 0.5_wp*dx
            end do

            return_value = 0.5_wp*(xl + xr)

        end function cbsrh


        subroutine ccompb(n, ierror, an, bn, cn, b, bc, ah, bh)
            !
            ! Purpose:
            !
            ! ccompb computes the roots of the b polynomials using subroutine
            ! ctevls which is a modification the eispack program tqlrat.
            ! ierror is set to 4 if either ctevls fails or if a(j+1)*c(j) is
            ! less than zero for some j.  ah, bh are temporary work arrays.
            !
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in)  :: n
            integer (ip), intent (out) :: ierror
            real (wp)                  :: an(*)
            real (wp),    intent (in)  :: bn(*)
            real (wp) :: cn(*)
            real (wp) :: b(*)
            real (wp) :: ah(*)
            real (wp) :: bh(*)
            complex (wp) :: bc(*)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip) :: j, if_rename, kdo, l, ir, i2, i4
            integer (ip) ::  ipl, ifd, i, ib, nb, js, jf
            integer (ip) ::  ls, lh, nmp, l1, l2, j2, j1, n2m2
            real (wp)    :: bnorm, arg, d1, d2, d3
            real (wp), parameter :: eps = epsilon(1.0_wp)
            !-----------------------------------------------

            bnorm = abs(bn(1))

            do j = 2, nm
                bnorm = max(bnorm, abs(bn(j)))
                arg = an(j)*cn(j-1)
                if (arg < 0.0_wp) then
                    ierror = 5
                    return
                end if
                b(j) = sign(sqrt(arg), an(j))
            end do

            cnv = eps*bnorm
            if_rename = 2**k
            kdo = k - 1

            outer_loop: do l = 1, kdo

                ir = l - 1
                i2 = 2**ir
                i4 = i2 + i2
                ipl = i4 - 1
                ifd = if_rename - i4

                do i = i4, ifd, i4

                    call cindxb(i, l, ib, nb)

                    if (nb <= 0) cycle outer_loop

                    js = i - ipl
                    jf = js + nb - 1
                    ls = 0
                    bh(:jf-js+1) = bn(js:jf)
                    ah(:jf-js+1) = b(js:jf)

                    call ctevls(nb, bh, ah, ierror)

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

                    if (arg < 0.0_wp) then
                        ierror = 5
                        return
                    end if

                    bh(j) = sign(sqrt(arg), (-an(l1)))
                    ah(j) = -bn(l1)
                end do

                call ctevls (nb, ah, bh, ierror)

                if (ierror /= 0) then
                    ierror = 4
                    return
                end if

                call cindxb(if_rename, k - 1, j2, lh)
                call cindxb(if_rename/2, k - 1, j1, lh)

                j2 = j2 + 1
                lh = j2
                n2m2 = j2 + 2*nm - 2

                iteration: do

                    d1 = abs(b(j1)-b(j2-1))
                    d2 = abs(b(j1)-b(j2))
                    d3 = abs(b(j1)-b(j2+1))

                    if (d2 >= d1 .or. d2 >= d3) then
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

                call cindxb(if_rename, k - 1, j1, j2)

                j2 = j1 + nmp + nmp

                call cppadd(nm + 1, ierror, an, cn, cmplx(b(j1:j1)), real(bc(j1:j1)), b(j2))

            end if


        end subroutine ccompb


        subroutine cproc(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, w, yy)
            !
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
            complex (wp), intent (in)     :: bd(*)
            complex (wp), intent (in)     :: x(*)
            complex (wp), intent (in out) :: y(*)
            complex (wp), intent (in)     :: a(*)
            complex (wp), intent (in)     :: b(*)
            complex (wp), intent (in)     :: c(*)
            complex (wp), intent (in out) :: d(*)
            complex (wp), intent (in out) :: w(*)
            complex (wp)                  :: yy(*)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip) :: j, mm, id, m1, m2, ia, iflg, k
            real (wp)    :: rt
            complex (wp) :: crt, den, y1, y2
            !-----------------------------------------------

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

                    if (abs(den) /= 0.0_wp) then
                        y(1) = (y(1)-c(1)*w(2))/den
                    else
                        y(1) = cmplx(1.0_wp, 0.0_wp, kind=wp)
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
                            if (abs(bm1(m1)) - abs(bm2(m2)) > 0.0_wp) then
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
                    ! scalar multiplication
                    !
                    y(:m) = rt*y(:m)
                end if

                if (iflg <= 0) exit main_loop

            end do main_loop

        end subroutine cproc



        subroutine cprocp(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, u, yy)
            !
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
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in) :: nd
            integer (ip), intent (in) :: nm1
            integer (ip), intent (in) :: nm2
            integer (ip), intent (in) :: na
            integer (ip), intent (in) :: m
            real (wp), intent (in) :: bm1(*)
            real (wp), intent (in) :: bm2(*)
            real (wp), intent (in) :: aa(*)
            real (wp) :: yy(*)
            complex (wp), intent (in) :: bd(*)
            complex (wp), intent (in) :: x(*)
            complex (wp), intent (in out) :: y(*)
            complex (wp), intent (in) :: a(*)
            complex (wp), intent (in) :: b(*)
            complex (wp), intent (in) :: c(*)
            complex (wp), intent (in out) :: d(*)
            complex (wp), intent (in out) :: u(*)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip) :: j, mm, mm2, id, m1, m2, ia, iflg, k
            real (wp) :: rt
            complex (wp) :: v, den, bh, ym, am, y1, y2, yh, crt
            !-----------------------------------------------

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
                    !
                    !==>  begin solution to system
                    !
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
                    if (abs(den) /= 0.0_wp) then
                        y(m) = (ym - am*y(m-1))/den
                    else
                        y(m) = cmplx(1.0_wp, 0.0_wp, kind=wp)
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
                            !
                            ! matrix multiplication
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
                    ! scalar multiplication
                    !
                    y(:m) = rt*y(:m)
                end if

                if (iflg <= 0) exit main_loop

            end do main_loop

        end subroutine cprocp

        subroutine cindxa(i, ir, idxa, na)

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
            if (i > nm) na = 0


        end subroutine cindxa


        subroutine cindxb(i, ir, idx, idp)

            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in) :: i
            integer (ip), intent (in) :: ir
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

        end subroutine cindxb


        subroutine cindxc(i, ir, idxc, nc)
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
            if (idxc + nc - 1 > nm) nc = 0

        end subroutine cindxc


        subroutine cppadd(n, ierror, a, c, cbp, bp, bh)
            !
            ! Purpose:
            !
            !     cppadd computes the eigenvalues of the periodic tridiagonal
            !     matrix with coefficients an, bn, cn
            !
            !     n is the order of the bh and bp polynomials
            !     on output bp contains the eigenvalues
            !     cbp is the same as bp except type complex
            !     bh is used to temporarily store the roots of the b hat polynomial
            !       which enters through bp
            !
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in)  :: n
            integer (ip), intent (out) :: ierror
            real (wp)                  :: a(*)
            real (wp)                  :: c(*)
            real (wp)                  :: bp(*)
            real (wp)                  :: bh(*)
            complex (wp)               :: cbp(*)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip)         :: iz, izm, izm2, j, nt, modiz
            integer (ip)         :: iis, iif, ig, it, icv, i3, i2, nhalf
            real (wp)            :: r4, r5, r6, scnv, xl, db, sgn, xr, xm, psg
            real (wp)            :: temp
            complex (wp)         :: cx, fsg, hsg, dd, f, fp, fpp, cdis, r1, r2, r3
            real (wp), parameter :: eps = epsilon(1.0_wp)
            !-----------------------------------------------

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
                        if (a(1) < 0.0_wp) exit block_110
                        if (a(1) == 0.0_wp) exit main_block
                    end if
                    xl = bh(1)
                    db = bh(3) - bh(1)
                    xl = xl - db
                    r4 = psgf(xl, iz, c, a, bh)
                    do while (r4 <= 0.0_wp)
                        xl = xl - db
                        r4 = psgf(xl, iz, c, a, bh)
                    end do
                    sgn = -1.0_wp
                    temp = cbsrh(xl, bh(1), iz, c, a, bh, psgf, sgn)
                    cbp(1) = cmplx(temp, 0.0_wp, kind=wp)
                    iis = 2
                end block block_110

                iif = iz - 1

                block_115: block
                    if (modiz /= 0) then
                        if (a(1) > 0.0_wp) exit block_115
                        if (a(1) == 0.0_wp) exit main_block
                    end if
                    xr = bh(iz)
                    db = bh(iz) - bh(iz-2)
                    xr = xr + db
                    r5 = psgf(xr, iz, c, a, bh)
                    do while (r5 < 0.0_wp)
                        xr = xr + db
                        r5 = psgf(xr, iz, c, a, bh)
                    end do
                    sgn = 1.0_wp
                    temp = cbsrh(bh(iz), xr, iz, c, a, bh, psgf, sgn)
                    cbp(iz) = cmplx(temp, 0.0_wp, kind=wp)
                    iif = iz - 2
                end block block_115

                main_loop: do ig = iis, iif, 2
                    xl = bh(ig)
                    xr = bh(ig+1)
                    sgn = -1.
                    xm = cbsrh(xl, xr, iz, c, a, bh, ppspf, sgn)
                    psg = psgf(xm, iz, c, a, bh)

                    if_block: block
                        if (abs(psg) > eps) then
                            r6 = psg*ppsgf(xm, iz, c, a, bh)
                            if (r6 > 0.0_wp) exit if_block
                            if (r6 /= 0.0_wp) then
                                sgn = 1.0_wp
                                cbp(ig) = cmplx(cbsrh(bh(ig), xm, iz, c, a, bh, psgf, sgn), 0.)
                                sgn = -1.0_wp
                                cbp(ig+1) = cmplx(cbsrh(xm, bh(ig+1), iz, c, a, bh, psgf, sgn), 0.)
                                cycle main_loop
                            !
                            !     case of a multiple zero
                            !
                            end if
                        end if
                        cbp(ig) = cmplx(xm, 0.0_wp, kind=wp)
                        cbp(ig+1) = cmplx(xm, 0.0_wp, kind=wp)
                        cycle main_loop
                    !
                    !     case of a complex zero
                    !
                    end block if_block

                    it = 0
                    icv = 0
                    cx = cmplx(xm, 0.0_wp, kind=wp)

                    loop_120: do
                        fsg = cmplx(1.0_wp, 0.0_wp, kind=wp)
                        hsg = cmplx(1.0_wp, 0.0_wp, kind=wp)
                        fp = 0.0_wp
                        fpp = 0.0_wp

                        do j = 1, iz
                            dd = 1./(cx - bh(j))
                            fsg = fsg*a(j)*dd
                            hsg = hsg*c(j)*dd
                            fp = fp + dd
                            fpp = fpp - dd*dd
                        end do

                        if (modiz == 0) then
                            f = cmplx(1.0_wp, 0.0_wp, kind=wp) - fsg - hsg
                        else
                            f = cmplx(1.0_wp, 0.0_wp, kind=wp) + fsg + hsg
                        end if

                        i3 = 0

                        if (abs(fp) > 0.0_wp) then
                            i3 = 1
                            r3 = -f/fp
                        end if

                        i2 = 0

                        if (abs(fpp) > 0.0_wp) then
                            i2 = 1
                            cdis = csqrt(fp**2 - 2.0_wp*f*fpp)
                            r1 = cdis - fp
                            r2 = (-fp) - cdis
                            if (abs(r1) - abs(r2) > 0.0_wp) then
                                r1 = r1/fpp
                            else
                                r1 = r2/fpp
                            end if
                            r2 = 2.0_wp*f/fpp/r1
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

                if (abs(cbp(n)) - abs(cbp(1)) <= 0.0_wp) then
                    if (abs(cbp(n)) - abs(cbp(1)) == 0.0_wp) exit main_block
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
                    if (aimag(cbp(j)) /= 0.0_wp) return
                end do

                ncmplx = 0

                do j = 2, iz
                    bp(j) = real(cbp(j), kind=wp)
                end do

                return
            end block main_block
            !
            !==> Procedure failed
            !
            ierror = 4

        end subroutine cppadd


        subroutine proc(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, w, u)
            !
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
            complex (wp), intent (in)     :: x(*)
            complex (wp), intent (in out) :: y(*)
            complex (wp), intent (in)     :: a(*)
            complex (wp), intent (in)     :: b(*)
            complex (wp), intent (in)     :: c(*)
            complex (wp), intent (in out) :: d(*)
            complex (wp), intent (in out) :: w(*)
            complex (wp)                  :: u(*)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip) :: j, mm, id, ibr, m1, m2, ia, k
            real (wp)    :: rt
            complex (wp) :: den
            !-----------------------------------------------

            w(:m) = x(:m)
            y(:m) = w(:m)
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
                    y(:m) = rt*w(:m)
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
                w(1) = cmplx(1.0_wp, 0.0_wp, kind=wp)

                if (abs(den) /= 0.) then
                    w(1) = (y(1)-c(1)*w(2))/den
                end if

                do j = 2, m
                    w(j) = w(j) - d(j)*w(j-1)
                end do

                if (na > 0) cycle main_loop

                if (m1 <= 0) then
                    if (m2 <= 0) then
                        y(:m) = w(:m)
                        ibr = 1
                        cycle main_loop
                    end if
                else
                    if (.not.(m2 > 0 .and. abs(bm1(m1)) <= abs(bm2(m2)))) then
                        if (ibr <= 0 .and. abs(bm1(m1)-bd(id)) < abs(bm1(m1)-rt)) then
                            y(:m) = w(:m)
                            ibr = 1
                            cycle main_loop
                        end if
                    end if
                    rt = rt - bm1(m1)
                    m1 = m1 - 1
                    y(:m) = y(:m) + rt*w(:m)
                    cycle main_loop
                end if

                if (ibr <= 0 .and. abs(bm2(m2)-bd(id)) < abs(bm2(m2)-rt)) then
                    y(:m) = w(:m)
                    ibr = 1
                    cycle main_loop
                end if

                rt = rt - bm2(m2)
                m2 = m2 - 1
                y(:m) = y(:m) + rt*w(:m)

            end do main_loop

        end subroutine proc



        subroutine procp(nd, bd, nm1, bm1, nm2, bm2, na, aa, x, y, m, a, b, c, d, u, w)
            !
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
            complex (wp), intent (in)     :: x(*)
            complex (wp), intent (in out) :: y(*)
            complex (wp), intent (in)     :: a(*)
            complex (wp), intent (in)     :: b(*)
            complex (wp), intent (in)     :: c(*)
            complex (wp), intent (in out) :: d(*)
            complex (wp), intent (in out) :: u(*)
            complex (wp), intent (in out) :: w(*)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip) :: j, mm, mm2, id, ibr, m1, m2, ia, k
            real (wp)    :: rt
            complex (wp) :: den, ym, v, bh, am
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

            main_loop: do

                if (ia > 0) then
                    if (nd == 0) then
                        rt = -aa(ia)
                    else
                        rt = aa(ia)
                    end if
                    ia = ia - 1
                    y(:m) = rt*w(:m)
                end if

                if (id <= 0) return

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

                if (abs(den) /= 0.0_wp) then
                    w(m) = (ym - am*w(m-1))/den
                else
                    w(m) = cmplx(1.0_wp, 0.0_wp, kind=wp)
                end if

                w(m-1) = w(m-1) - d(m-1)*w(m)

                do j = 2, mm
                    k = m - j
                    w(k) = w(k) - d(k)*w(k+1) - u(k)*w(m)
                end do

                if (na > 0) cycle main_loop

                if (m1 <= 0) then
                    if (m2 <= 0) then
                        y(:m) = w(:m)
                        ibr = 1
                        cycle main_loop
                    end if
                else
                    if (.not.(m2 > 0 .and. abs(bm1(m1)) <= abs(bm2(m2)))) then
                        if (ibr <= 0 .and. abs(bm1(m1)-bd(id)) < abs(bm1(m1)-rt)) then
                            y(:m) = w(:m)
                            ibr = 1
                            cycle main_loop
                        end if
                        rt = rt - bm1(m1)
                        m1 = m1 - 1
                        y(:m) = y(:m) + rt*w(:m)
                        cycle main_loop
                    end if
                end if

                if (ibr <= 0 .and. abs(bm2(m2)-bd(id)) < abs(bm2(m2)-rt)) then
                    y(:m) = w(:m)
                    ibr = 1
                    cycle main_loop
                end if

                rt = rt - bm2(m2)
                m2 = m2 - 1
                y(:m) = y(:m) + rt*w(:m)

            end do main_loop

        end subroutine procp


        subroutine ctevls(n, d, e2, ierr)
            !
            ! Purpose:
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
            !     eps is a machine dependent parameter specifying
            !     the relative precision of floating point arithmetic.
            !
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
            real (wp), parameter :: eps = epsilon(1.0_wp)
            !-----------------------------------------------

            ierr = 0
            if (n /= 1) then

                e2(:n-1) = e2(2:n)*e2(2:n)
                f = 0.0_wp
                b = 0.0_wp
                e2(n) = 0.0_wp

                main_loop: do l = 1, n
                    j = 0
                    h = eps*(abs(d(l))+sqrt(e2(l)))

                    if (b <= h) then
                        b = h
                        c = b*b
                    end if
                    !
                    !==> look for small squared sub-diagonal element
                    !
                    do m = l, n
                        if (e2(m) > c) cycle
                        exit
                    !
                    !==> 2(n) is always zero, so there is no exit
                    !    through the bottom of the loop
                    !
                    end do

                    if_block: block
                        if (m /= l) then
                            loop_105: do
                                if (j == 30) then
                                    !
                                    !==> set error -- no convergence to an
                                    !    eigenvalue after 30 iterations
                                    !
                                    ierr = l
                                    return
                                end if

                                j = j + 1
                                !
                                !==> form shift
                                !
                                l1 = l + 1
                                s = sqrt(e2(l))
                                g = d(l)
                                p = (d(l1)-g)/(2.0_wp*s)
                                r = sqrt(p**2 + 1.0_wp)
                                d(l) = s/(p + sign(r, p))
                                h = g - d(l)
                                d(l1:n) = d(l1:n) - h
                                f = f + h
                                !
                                !==> rational ql transformation
                                !
                                g = d(m)

                                if (g == 0.0_wp) g = b

                                h = g
                                s = 0.0_wp
                                mml = m - l
                                !
                                !==> for i=m-1 step -1 until l do --
                                !
                                do ii = 1, mml
                                    i = m - ii
                                    p = g*h
                                    r = p + e2(i)
                                    e2(i+1) = s*r
                                    s = e2(i)/r
                                    d(i+1) = h + s*(h + d(i))
                                    g = d(i) - e2(i)/g
                                    if (g == 0.0_wp) g = b
                                    h = g*p/r
                                end do

                                e2(l) = s*g
                                d(l) = h
                                !
                                !==>  guard against underflowed h
                                !
                                if (h == 0.0_wp .or. abs(e2(l)) <= abs(c/h)) exit if_block

                                e2(l) = h*e2(l)

                                if (e2(l) == 0.0_wp) exit loop_105
                            end do loop_105
                        end if
                    end block if_block

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

        end subroutine ctevls

    end subroutine cblktri

end module module_cblktri
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
! May       2016    Fortran 2008 changes
!
