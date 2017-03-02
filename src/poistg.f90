!
!     file poistg.f90
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                         Fishpack                              *
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
! SUBROUTINE poistg(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
!
!
! DIMENSION OF           a(m),  b(m),  c(m),  y(idimy, n)
! ARGUMENTS
!
! LATEST REVISION        April 2016
!
! PURPOSE                Solves the linear system of equations
!                        for unknown x values, where i=1, 2, ..., m
!                        and j=1, 2, ..., n
!
!                        a(i)*x(i-1, j) + b(i)*x(i, j) + c(i)*x(i+1, j)
!                        + x(i, j-1) - TWO *x(i, j) + x(i, j+1)
!                        = y(i, j)
!
!                        The indices i+1 and i-1 are evaluated modulo m,
!                        i.e. x(0, j) = x(m, j) and x(m+1, j) = x(1, j), and
!                        x(i, 0) may be equal to x(i, 1) or -x(i, 1), and
!                        x(i, n+1) may be equal to x(i, n) or -x(i, n),
!                        depending on an input parameter.
!
! USAGE                  call poistg(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
!
! ARGUMENTS
! ON INPUT
!
!                        nperod
!                          Indicates values which x(i, 0) and x(i, n+1)
!                          are assumed to have.
!                          = 1 if x(i, 0) = -x(i, 1) and x(i, n+1) = -x(i, n
!                          = 2 if x(i, 0) = -x(i, 1) and x(i, n+1) =  x(i, n
!                          = 3 if x(i, 0) =  x(i, 1) and x(i, n+1) =  x(i, n
!                          = 4 if x(i, 0) =  x(i, 1) and x(i, n+1) = -x(i, n
!
!                        n
!                          the number of unknowns in the j-direction.
!                          n must be greater than 2.
!
!                        mperod
!                          = 0 if a(1) and c(m) are not zero
!                          = 1 if a(1) = c(m) = 0
!
!                        m
!                          the number of unknowns in the i-direction.
!                          m must be greater than 2.
!
!                        a, b, c
!                          One-dimensional arrays of length m that
!                          specify the coefficients in the linear
!                          equations given above. If mperod = 0 the
!                          array elements must not depend on index i,
!                          but must be constant. Specifically, the
!                          subroutine checks the following condition
!                            a(i) = c(1)
!                            b(i) = b(1)
!                            c(i) = c(1)
!                          for i = 1, 2, ..., m.
!
!                        idimy
!                          The row (or first) dimension of the two-
!                          dimensional array y as it appears in the
!                          program calling poistg. This parameter is
!                          used to specify the variable dimension of y.
!                          idimy must be at least m.
!
!                        y
!                          A two-dimensional array that specifies the
!                          values of the right side of the linear system
!                          of equations given above.
!                          y must be dimensioned at least m x n.
!
! ON OUTPUT
!
!                        Y
!                          Contains the solution x.
!
!                        ierror
!                          An error flag that indicates invalid input
!                          parameters.  except for number zero, a
!                          solution is not attempted.
!                          = 0  no error
!                          = 1  if m <= 2
!                          = 2  if n <= 2
!                          = 3  idimy < m
!                          = 4  if nperod < 1 or nperod > 4
!                          = 5  if mperod < 0 or mperod > 1
!                          = 6  if mperod = 0 and a(i) /= c(1)
!                               or b(i) /= b(1) or c(i) /= c(1)
!                               for some i = 1, 2, ..., m.
!                          = 7  if mperod == 1 and (a(1)/= 0 or c(m)/= 0)
!                          = 20 If the dynamic allocation of real and
!                               complex workspace required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
!
!                          Since this is the only means of indicating a
!                          possibly incorrect call to pois3d, the user
!                          should test ierror after the call.
!
!
!
! I/O                    None
!
! PRECISION              64-bit precision float and 32-bit precision integer
!
! REQUIRED FILES         type_FishpackWorkspace.f90, type_CyclicReductionUtility.f9090
!
! STANDARD               Fortran 2008
!
! HISTORY                Written by Roland Sweet at NCAR in the late
!                        1970's.  Released on NCAR'S public software
!                        libraries in January, 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated workspace.
!                        * Revised in April 2016 by Jon Lo Kim Lin to
!                        incorporate features of modern Fortran (2008+)
!
!
! ALGORITHM              This subroutine is an implementation of the
!                        algorithm presented in the reference below.
!
! TIMING                 For large m and n, the execution time is
!                        roughly proportional to m*n*log2(n).
!
! ACCURACY               To measure the accuracy of the algorithm a
!                        uniform random number generator was used to
!                        create a solution array x for the system given
!                        in the 'purpose' section above, with
!
!                          a(i) = c(i) = -HALF * b(i) = 1,    i=1, 2, ..., m
!                        and, when mperod = 1
!                          a(1) = c(m) = 0
!                          b(1) = b(m) =-1.
!
!                        The solution x was substituted into the given
!                        system and, using double precision, a right sid
!                        y was computed.  using this array y subroutine
!                        poistg was called to produce an approximate
!                        solution z.  then the relative error, defined a
!
!                         e = max(abs(z(i, j)-x(i, j)))/max(abs(x(i, j)))
!
!                        where the two maxima are taken over i=1, 2, ..., m
!                        and j=1, 2, ..., n, was computed.  values of e are
!                        given in the table below for some typical value
!                        of m and n.
!
!                        m (=n)    mperod    nperod      e
!                        ------    ------    ------    ------
!
!                          31        0-1       1-4     9.e-13
!                          31        1         1       4.e-13
!                          31        1         3       3.e-13
!                          32        0-1       1-4     3.e-12
!                          32        1         1       3.e-13
!                          32        1         3       1.e-13
!                          33        0-1       1-4     1.e-12
!                          33        1         1       4.e-13
!                          33        1         3       1.e-13
!                          63        0-1       1-4     3.e-12
!                          63        1         1       1.e-12
!                          63        1         3       2.e-13
!                          64        0-1       1-4     4.e-12
!                          64        1         1       1.e-12
!                          64        1         3       6.e-13
!                          65        0-1       1-4     2.e-13
!                          65        1         1       1.e-11
!                          65        1         3       4.e-13
!
! REFERENCES             Schumann, U. and R. Sweet, "A direct method
!                        for the solution of Poisson's equation with
!                        Neumann boundary conditions on a staggered
!                        grid of arbitrary size, " J. Comp. Phys.
!                        20 (1976), pp. 171-182.
!
module module_poistg

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_CyclicReductionUtility, only: &
        CyclicReductionUtility

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: poistg
    public :: poistgg

    !---------------------------------------------------------------
    ! Parameters confined to the module
    !---------------------------------------------------------------
    real(wp),    parameter :: ZERO = 0.0_wp
    real(wp),    parameter :: HALF = 0.5_wp
    real(wp),    parameter :: ONE = 1.0_wp
    real(wp),    parameter :: TWO = 2.0_wp
    integer(ip), parameter :: IIWK = 11_ip ! Size of workspace indices
    !---------------------------------------------------------------

contains

    subroutine poistg(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in)     :: nperod
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: mperod
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: idimy
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: a(:)
        real(wp),    intent(in)     :: b(:)
        real(wp),    intent(in)     :: c(:)
        real(wp),    intent(inout)  :: y(:,:)
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: irwk, icwk
        type(FishpackWorkspace)  :: workspace
        !--------------------------------------------------------------

        ! Compute workspace dimensions
        irwk = m * (9 + int(log(real(n, kind=wp))/log(TWO),kind=ip)) + 4 * n
        icwk = 0

        ! Allocate memory
        call workspace%create(irwk, icwk)

        ! Solve system
        associate( rew => workspace%real_workspace )
            call poistgg(nperod, n, mperod, m, a, b, c, idimy, y, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine poistg

    subroutine poistgg(nperod, n, mperod, m, a, b, c, idimy, y, ierror, w)
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip), intent(in)     :: nperod
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: mperod
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: idimy
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: a(m)
        real(wp),    intent(in)     :: b(m)
        real(wp),    intent(in)     :: c(m)
        real(wp),    intent(inout)  :: y(idimy,n)
        real(wp),    intent(inout)  :: w(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip)     :: workspace_indices(IIWK)
        integer(ip)     :: i, k, j, np, mp
        integer(ip)     :: nby2, mskip, ipstor, irev, mh, mhm1, m_odd
        real(wp)        :: temp
        !-----------------------------------------------

        ! Check validity of calling arguments
        call check_input_arguments(nperod, n, mperod, m, idimy, ierror, a, b, c)

        ! Check error flag
        if (ierror /= 0) return

        ! Compute workspace indices
        workspace_indices = get_workspace_indices(n,m)

        associate( &
            iwba => workspace_indices(1), &
            iwbb => workspace_indices(2), &
            iwbc => workspace_indices(3), &
            iwb2 => workspace_indices(4), &
            iwb3 => workspace_indices(5), &
            iww1 => workspace_indices(6), &
            iww2 => workspace_indices(7), &
            iww3 => workspace_indices(8), &
            iwd => workspace_indices(9), &
            iwtcos => workspace_indices(10), &
            iwp => workspace_indices(11) &
            )

            do i = 1, m
                k = iwba + i - 1
                w(k) = -a(i)
                k = iwbc + i - 1
                w(k) = -c(i)
                k = iwbb + i - 1
                w(k) = TWO - b(i)
                y(i, :n) = -y(i, :n)
            end do

            np = nperod
            mp = mperod + 1

            if (mp == 1) then
                mh = (m + 1)/2
                mhm1 = mh - 1

                if (mh*2 == m) then
                    m_odd = 2
                else
                    m_odd = 1
                end if

                do j = 1, n
                    do i = 1, mhm1
                        w(i) = y(mh-i, j) - y(i+mh, j)
                        w(i+mh) = y(mh-i, j) + y(i+mh, j)
                    end do
                    w(mh) = TWO * y(mh, j)
                    select case (m_odd)
                        case (1)
                            y(:m, j) = w(:m)
                        case (2)
                            w(m) = TWO * y(m, j)
                            y(:m, j) = w(:m)
                    end select
                end do

                k = iwbc + mhm1 - 1
                i = iwba + mhm1
                w(k) = ZERO
                w(i) = ZERO
                w(k+1) = TWO * w(k+1)

                select case (m_odd)
                    case (2)
                        w(iwbb-1) = w(k+1)
                    case default
                        k = iwbb + mhm1 - 1
                        w(k) = w(k) - w(i-1)
                        w(iwbc-1) = w(iwbc-1) + w(iwbb-1)
                end select
            end if

            loop_107: do
                loop_108: do
                    if (nperod /= 4) then

                        call solve_poisson_staggered_grid(np, n, m, w(iwba:), w(iwbb:), w(iwbc:), idimy, y, &
                            w, w(iwb2:), w(iwb3:), w(iww1:), w(iww2:), w(iww3:), &
                            w(iwd:), w(iwtcos:), w(iwp:))

                        ipstor = int(w(iww1), kind=ip)
                        irev = 2
                        select case (mp)
                            case (1)
                                exit loop_107
                            case (2)
                                w(1) = real(ipstor + iwp - 1, kind=wp)
                                return
                        end select
                        cycle loop_107
                    end if

                    irev = 1
                    nby2 = n/2
                    np = 2

                    do j = 1, nby2
                        mskip = n + 1 - j
                        do i = 1, m
                            temp = y(i, j)
                            y(i, j) = y(i, mskip)
                            y(i, mskip) = temp
                        end do
                    end do

                    select case (irev)
                        case (1)
                            cycle loop_108
                        case (2)
                            select case (mp)
                                case (1)
                                    exit loop_107
                                case (2)
                                    w(1) = real(ipstor + iwp - 1, kind=wp)
                                    return
                            end select
                            cycle loop_107
                    end select

                    exit loop_108
                end do loop_108
                exit loop_107
            end do loop_107

            do j = 1, n
                w(mh-1:mh-mhm1:(-1)) = HALF * (y(mh+1:mhm1+mh, j)+y(:mhm1, j))
                w(mh+1:mhm1+mh) = HALF * (y(mh+1:mhm1+mh, j)-y(:mhm1, j))
                w(mh) = HALF * y(mh, j)
                select case (m_odd)
                    case (1)
                        y(:m, j) = w(:m)
                    case (2)
                        w(m) = HALF * y(m, j)
                        y(:m, j) = w(:m)
                end select
            end do

            w(1) = real(ipstor + iwp - 1, kind=wp)

        end associate

    end subroutine poistgg

    pure subroutine check_input_arguments(nperod, n, mperod, m, idimy, ierror, a, b, c)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in)  :: nperod
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: mperod
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: idimy
        integer(ip), intent(out) :: ierror
        real(wp),    intent(in)  :: a(m)
        real(wp),    intent(in)  :: b(m)
        real(wp),    intent(in)  :: c(m)
        !--------------------------------------------------------------

        if (3 > m) then
            ierror = 1
            return
        else if (3 > n) then
            ierror = 2
            return
        else if (idimy < m) then
            ierror = 3
            return
        else if (nperod < 1 .or. nperod > 4) then
            ierror = 4
            return
        else if (mperod < 0 .or. mperod > 1) then
            ierror = 5
            return
        else if (mperod /= 1) then
            if (any(a /= c(1)) .or. any(c /= c(1)) .or. any(b /= b(1))) then
                ierror = 6
                return
            end if
        else if (a(1) /= ZERO .or. c(m) /= ZERO) then
            ierror = 7
            return
        else
            ierror = 0
        end if

    end subroutine check_input_arguments

    function get_workspace_indices(n, m) result (return_value)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in) :: n
        integer(ip), intent(in) :: m
        integer(ip)             :: return_value(IIWK)
        !--------------------------------------------------------------
        integer(ip) :: j !! Counter
        !--------------------------------------------------------------

        associate( i => return_value)
            i(1) = m + 1
            do j = 2, IIWK - 1
                i(j) = i(j-1) + m
            end do
            i(IIWK) = i(IIWK - 1) + 4 * n
        end associate

    end function get_workspace_indices

    subroutine solve_poisson_staggered_grid(nperod, n, m, a, bb, c, idimq, q, b, b2, b3, w, &
        w2, w3, d, tcos, p)
        !
        ! Purpose:
        !
        ! Subroutine to solve poisson's equation on a staggered grid.
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in)    :: nperod
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: idimq
        real(wp),    intent(in)    :: a(m)
        real(wp),    intent(in)    :: bb(m)
        real(wp),    intent(in)    :: c(m)
        real(wp),    intent(inout) :: q(idimq, n)
        real(wp),    intent(inout) :: b(m)
        real(wp),    intent(inout) :: b2(m)
        real(wp),    intent(inout) :: b3(m)
        real(wp),    intent(inout) :: w(m)
        real(wp),    intent(inout) :: w2(m)
        real(wp),    intent(inout) :: w3(m)
        real(wp),    intent(inout) :: d(m)
        real(wp),    intent(inout) :: tcos(m)
        real(wp),    intent(inout) :: p(4*n)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip)     :: k(4)
        integer(ip)     :: np, mr, ipp, ipstor, i2r, jr, nr, nlast
        integer(ip)     :: kr, lr, nrod, jstart, jstop, i2rby2
        integer(ip)     :: j, ijump, jp1, jp2, jp3
        integer(ip)     :: jm1, jm2, jm3, i, nrodpr, ii, nlastp, jstep
        real(wp)        :: fnum, fnum2, fi, t
        type(CyclicReductionUtility) :: genbun_aux
        !-----------------------------------------------

        associate( &
            k1 => k(1), &
            k2 => k(2), &
            k3 => k(3), &
            k4 => k(4) &
            )

            np = nperod
            fnum = HALF * real(np/3, kind=wp)
            fnum2 = HALF * real(np/2, kind=wp)
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

                loop_101: do

                    jr = 2*i2r

                    if ((nr/2)*2 == nr) then
                        nrod = 0
                    else
                        nrod = 1
                    end if

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
                                call genbun_aux%generate_cosines(i2r, 1, fnum, HALF, tcos)
                                select case (i2r)
                                    case (1)
                                        b(:mr) = q(:mr, 1)
                                        q(:mr, 1) = q(:mr, 2)
                                    case default
                                        b(:mr) = q(:mr, 1) + HALF * (q(:mr, jp2)-q(:mr, jp1)-q(:mr, &
                                            jp3))
                                        q(:mr, 1) = q(:mr, jp2) + q(:mr, 1) - q(:mr, jp1)
                                end select
                            else
                                if (ijump == 1) then
                                    ijump = 2
                                    call genbun_aux%generate_cosines(i2r, 1, HALF, ZERO, tcos)
                                end if

                                select case (i2r)
                                    case (1)
                                        b(:mr) = TWO*q(:mr, j)
                                        q(:mr, j) = q(:mr, jm2) + q(:mr, jp2)
                                    case default
                                        do i = 1, mr
                                            fi = q(i, j)
                                            q(i, j)=q(i, j)-q(i, jm1)-q(i, jp1)+q(i, jm2)+q(i, jp2)
                                            b(i) = fi + q(i, j) - q(i, jm3) - q(i, jp3)
                                        end do
                                end select
                            end if

                            call genbun_aux%solve_tridiag(i2r, 0, mr, a, bb, c, b, tcos, d, w)

                            q(:mr, j) = q(:mr, j) + b(:mr)
                            !
                            ! End of reduction for regular unknowns.
                            !
                        end do
                        !
                        ! Begin special reduction for last unknown.
                        !
                        j = jstop + jr
                    end if

                    nlast = j
                    jm1 = j - i2rby2
                    jm2 = j - i2r
                    jm3 = jm2 - i2rby2

                    if (nrod /= 0) then
                        !
                        ! Odd number of unknowns
                        !
                        select case (i2r)
                            case (1)
                                b(:mr) = q(:mr, j)
                                q(:mr, j) = q(:mr, jm2)
                            case default
                                b(:mr)=q(:mr, j)+HALF * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))

                                select case (nrodpr)
                                    case (0)
                                        q(:mr, j) = q(:mr, jm2) + p(ipp+1:mr+ipp)
                                        ipp = ipp - mr
                                    case default
                                        q(:mr, j) = q(:mr, j) - q(:mr, jm1) + q(:mr, jm2)
                                end select

                                if (lr /= 0) call genbun_aux%generate_cosines(lr, 1, fnum2, HALF, tcos(kr+1:))

                        end select

                        call genbun_aux%generate_cosines(kr, 1, fnum2, HALF, tcos)
                        call genbun_aux%solve_tridiag(kr, lr, mr, a, bb, c, b, tcos, d, w)

                        q(:mr, j) = q(:mr, j) + b(:mr)
                        kr = kr + i2r
                    else
                        jp1 = j + i2rby2
                        jp2 = j + i2r

                        select case (i2r)
                            case (1)
                                b(:mr) = q(:mr, j)
                                tcos(1) = ZERO

                                call genbun_aux%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)

                                ipp = 0
                                ipstor = mr
                                p(:mr) = b(:mr)
                                b(:mr) = b(:mr) + q(:mr, n)
                                tcos(1) = -ONE + TWO*real(np/2, kind=wp)
                                tcos(2) = ZERO

                                call genbun_aux%solve_tridiag(1, 1, mr, a, bb, c, b, tcos, d, w)

                                q(:mr, j) = q(:mr, jm2) + p(:mr) + b(:mr)

                            case default

                                b(:mr) = q(:mr, j) + HALF * (q(:mr, jm2)-q(:mr, jm1)-q(:mr, jm3))

                                select case (nrodpr)
                                    case (0)
                                        b(:mr) = b(:mr) + p(ipp+1:mr+ipp)
                                    case default
                                        b(:mr) = b(:mr) + q(:mr, jp2) - q(:mr, jp1)
                                end select

                                call genbun_aux%generate_cosines(i2r, 1, HALF, ZERO, tcos)
                                call genbun_aux%solve_tridiag(i2r, 0, mr, a, bb, c, b, tcos, d, w)

                                ipp = ipp + mr
                                ipstor = max(ipstor, ipp + mr)
                                p(ipp+1:mr+ipp) = b(:mr) + HALF * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
                                b(:mr) = p(ipp+1:mr+ipp) + q(:mr, jp2)

                                if (lr /= 0) then
                                    call genbun_aux%generate_cosines(lr, 1, fnum2, HALF, tcos(i2r+1:))
                                    call genbun_aux%merge_cosines(tcos, 0, i2r, i2r, lr, kr)
                                else
                                    do i = 1, i2r
                                        ii = kr + i
                                        tcos(ii) = tcos(i)
                                    end do
                                end if

                                call genbun_aux%generate_cosines(kr, 1, fnum2, HALF, tcos)
                                call genbun_aux%solve_tridiag(kr, kr, mr, a, bb, c, b, tcos, d, w)

                                q(:mr, j) = q(:mr, jm2) + p(ipp+1:mr+ipp) + b(:mr)

                        end select
                        lr = kr
                        kr = kr + jr
                    end if

                    nr = (nlast - 1)/jr + 1

                    if (nr <= 3) exit loop_101

                    i2r = jr
                    nrodpr = nrod

                end do loop_101
            end if


            j = 1 + jr
            jm1 = j - i2r
            jp1 = j + i2r
            jm2 = nlast - i2r

            block_construct: block
                if_nr:  if (nr /= 2) then
                    if_lr:  if (lr == 0) then
                        if_n: if (n == 3) then
                            !
                            ! case n = 3.
                            !
                            select case (np)
                                case (1,3)
                                    b(:mr) = q(:mr, 2)
                                    b2(:mr) = q(:mr, 1) + q(:mr, 3)
                                    b3(:mr) = ZERO

                                    select case (np)
                                        case (1:2)
                                            tcos(1) = -TWO
                                            tcos(2) = ONE
                                            tcos(3) = -ONE
                                            k1 = 2
                                        case default
                                            tcos(1) = -ONE
                                            tcos(2) = ONE
                                            k1 = 1
                                    end select

                                    k2 = 1
                                    k3 = 0
                                    k4 = 0
                                case (2)
                                    b(:mr) = q(:mr, 2)
                                    b2(:mr) = q(:mr, 3)
                                    b3(:mr) = q(:mr, 1)

                                    call genbun_aux%generate_cosines(3, 1, HALF, ZERO, tcos)

                                    tcos(4) = -ONE
                                    tcos(5) = ONE
                                    tcos(6) = -ONE
                                    tcos(7) = ONE

                                    k1 = 3
                                    k2 = 2
                                    k3 = 1
                                    k4 = 1
                            end select


                            call genbun_aux%solve_tridiag3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)

                            b(:mr) = b(:mr) + b2(:mr) + b3(:mr)

                            if (np == 3) then
                                tcos(1) = TWO
                                call genbun_aux%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)
                            end if

                            q(:mr, 2) = b(:mr)
                            b(:mr) = q(:mr, 1) + b(:mr)
                            tcos(1) = -ONE + 4.0_wp*fnum

                            call genbun_aux%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)

                            q(:mr, 1) = b(:mr)
                            jr = 1
                            i2r = 0
                            exit block_construct
                        end if if_n
                        !
                        !     case n = 2**p+1
                        !
                        b(:mr)=q(:mr, j)+q(:mr, 1)-q(:mr, jm1)+q(:mr, nlast)-q(:mr, jm2)

                        select case (np)
                            case (1, 3)
                                b2(:mr) = q(:mr, 1) + q(:mr, nlast) &
                                    + q(:mr, j) - q(:mr, jm1) - q(:mr, jp1)
                                b3(:mr) = ZERO
                                k1 = nlast - 1
                                k2 = nlast + jr - 1

                                call genbun_aux%generate_cosines(jr - 1, 1, ZERO, ONE, tcos(nlast:))

                                tcos(k2) = TWO*real(np - 2, kind=wp)

                                call genbun_aux%generate_cosines(jr, 1, HALF - fnum, HALF, tcos(k2+1:))

                                k3 = (3 - np)/2

                                call genbun_aux%merge_cosines(tcos, k1, jr - k3, k2 - k3, jr + k3, 0)

                                k1 = k1 - 1 + k3

                                call genbun_aux%generate_cosines(jr, 1, fnum, HALF, tcos(k1+1:))

                                k2 = jr
                                k3 = 0
                                k4 = 0
                            case (2)
                                do i = 1, mr
                                    fi = (q(i, j)-q(i, jm1)-q(i, jp1))/2
                                    b2(i) = q(i, 1) + fi
                                    b3(i) = q(i, nlast) + fi
                                end do

                                k1 = nlast + jr - 1
                                k2 = k1 + jr - 1

                                call genbun_aux%generate_cosines(jr - 1, 1, ZERO, ONE, tcos(k1+1:))
                                call genbun_aux%generate_cosines(nlast, 1, HALF, ZERO, tcos(k2+1:))
                                call genbun_aux%merge_cosines(tcos, k1, jr - 1, k2, nlast, 0)

                                k3 = k1 + nlast - 1
                                k4 = k3 + jr

                                call genbun_aux%generate_cosines(jr, 1, HALF, HALF, tcos(k3+1:))
                                call genbun_aux%generate_cosines(jr, 1, ZERO, HALF, tcos(k4+1:))
                                call genbun_aux%merge_cosines(tcos, k3, jr, k4, jr, k1)

                                k2 = nlast - 1
                                k3 = jr
                                k4 = jr
                        end select

                        call genbun_aux%solve_tridiag3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)
                        b(:mr) = b(:mr) + b2(:mr) + b3(:mr)

                        if (np == 3) then
                            tcos(1) = TWO
                            call genbun_aux%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)
                        end if

                        q(:mr, j) = b(:mr) + HALF * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
                        b(:mr) = q(:mr, j) + q(:mr, 1)

                        call genbun_aux%generate_cosines(jr, 1, fnum, HALF, tcos)
                        call genbun_aux%solve_tridiag(jr, 0, mr, a, bb, c, b, tcos, d, w)

                        q(:mr, 1) = q(:mr, 1) - q(:mr, jm1) + b(:mr)

                        exit block_construct
                    end if if_lr
                    !
                    ! case of general n with nr = 3 .
                    !
                    b(:mr) = q(:mr, 1) - q(:mr, jm1) + q(:mr, j)

                    if (nrod == 0) then
                        b(:mr) = b(:mr) + p(ipp+1:mr+ipp)
                    else
                        b(:mr) = b(:mr) + q(:mr, nlast) - q(:mr, jm2)
                    end if

                    do i = 1, mr
                        t = HALF * (q(i, j)-q(i, jm1)-q(i, jp1))
                        q(i, j) = t
                        b2(i) = q(i, nlast) + t
                        b3(i) = q(i, 1) + t
                    end do

                    k1 = kr + 2*jr

                    call genbun_aux%generate_cosines(jr - 1, 1, ZERO, ONE, tcos(k1+1:))

                    k2 = k1 + jr
                    tcos(k2) = TWO*real(np - 2, kind=wp)
                    k4 = (np - 1)*(3 - np)
                    k3 = k2 + 1 - k4

                    !generate_cosines(n, ijump, fnum, fden, a)
                    associate( &
                        n_arg => kr+jr+k4, &
                        fnum_arg => real(k4, kind=wp)/2, &
                        fden_arg => ONE-real(k4, kind=wp) &
                        )
                        call genbun_aux%generate_cosines(n_arg, 1, fnum_arg , fden_arg, tcos(k3:))
                    end associate

                    k4 = 1 - np/3

                    call genbun_aux%merge_cosines(tcos, k1, jr - k4, k2 - k4, kr + jr + k4, 0)

                    if (np == 3) k1 = k1 - 1

                    k2 = kr + jr
                    k4 = k1 + k2

                    call genbun_aux%generate_cosines(kr, 1, fnum2, HALF, tcos(k4+1:))

                    k3 = k4 + kr

                    call genbun_aux%generate_cosines(jr, 1, fnum, HALF, tcos(k3+1:))
                    call genbun_aux%merge_cosines(tcos, k4, kr, k3, jr, k1)

                    k4 = k3 + jr

                    call genbun_aux%generate_cosines(lr, 1, fnum2, HALF, tcos(k4+1:))
                    call genbun_aux%merge_cosines(tcos, k3, jr, k4, lr, k1 + k2)
                    call genbun_aux%generate_cosines(kr, 1, fnum2, HALF, tcos(k3+1:))

                    k3 = kr
                    k4 = kr

                    call genbun_aux%solve_tridiag3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)

                    b(:mr) = b(:mr) + b2(:mr) + b3(:mr)

                    if (np == 3) then
                        tcos(1) = TWO
                        call genbun_aux%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)
                    end if

                    q(:mr, j) = q(:mr, j) + b(:mr)
                    b(:mr) = q(:mr, 1) + q(:mr, j)

                    call genbun_aux%generate_cosines(jr, 1, fnum, HALF, tcos)
                    call genbun_aux%solve_tridiag(jr, 0, mr, a, bb, c, b, tcos, d, w)

                    if (jr == 1) then
                        q(:mr, 1) = b(:mr)
                        exit block_construct
                    end if

                    q(:mr, 1) = q(:mr, 1) - q(:mr, jm1) + b(:mr)

                    exit block_construct
                end if if_nr

                b3(:mr) = ZERO
                b(:mr) = q(:mr, 1) + p(ipp+1:mr+ipp)
                q(:mr, 1) = q(:mr, 1) - q(:mr, jm1)
                b2(:mr) = q(:mr, 1) + q(:mr, nlast)
                k1 = kr + jr
                k2 = k1 + jr

                call genbun_aux%generate_cosines(jr - 1, 1, ZERO, ONE, tcos(k1+1:))

                select case (np)
                    case (1, 3)
                        tcos(k2) = TWO*real(np - 2, kind=wp)
                        call genbun_aux%generate_cosines(kr, 1, ZERO, ONE, tcos(k2+1:))
                    case (2)
                        call genbun_aux%generate_cosines(kr + 1, 1, HALF, ZERO, tcos(k2:))
                end select

                k4 = 1 - np/3

                call genbun_aux%merge_cosines(tcos, k1, jr - k4, k2 - k4, kr + k4, 0)

                if (np == 3) k1 = k1 - 1

                k2 = kr

                call genbun_aux%generate_cosines(kr, 1, fnum2, HALF, tcos(k1+1:))

                k4 = k1 + kr

                call genbun_aux%generate_cosines(lr, 1, fnum2, HALF, tcos(k4+1:))

                k3 = lr
                k4 = 0

                call genbun_aux%solve_tridiag3(mr, a, bb, c, k, b, b2, b3, tcos, d, w, w2, w3)

                b(:mr) = b(:mr) + b2(:mr)

                if (np == 3) then
                    tcos(1) = TWO
                    call genbun_aux%solve_tridiag(1, 0, mr, a, bb, c, b, tcos, d, w)
                end if

                q(:mr, 1) = q(:mr, 1) + b(:mr)

            end block block_construct

            loop_188: do

                j = nlast - jr
                b(:mr) = q(:mr, nlast) + q(:mr, j)
                jm2 = nlast - i2r

                if (jr == 1) then
                    q(:mr, nlast) = ZERO
                else
                    select case (nrod)
                        case (0)
                            q(:mr, nlast) = p(ipp+1:mr+ipp)
                            ipp = ipp - mr
                        case default
                            q(:mr, nlast) = q(:mr, nlast) - q(:mr, jm2)
                    end select
                end if

                call genbun_aux%generate_cosines(kr, 1, fnum2, HALF, tcos)
                call genbun_aux%generate_cosines(lr, 1, fnum2, HALF, tcos(kr+1:))
                call genbun_aux%solve_tridiag(kr, lr, mr, a, bb, c, b, tcos, d, w)

                q(:mr, nlast) = q(:mr, nlast) + b(:mr)
                nlastp = nlast

                loop_197: do

                    jstep = jr
                    jr = i2r
                    i2r = i2r/2

                    if (jr == 0) then
                        w(1) = real(ipstor, kind=wp)
                        return
                    end if

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

                    call genbun_aux%generate_cosines(jr, 1, HALF, ZERO, tcos)

                    do j = jstart, jstop, jstep
                        jm2 = j - jr
                        jp2 = j + jr

                        if (j == jr) then
                            b(:mr) = q(:mr, j) + q(:mr, jp2)
                        else
                            b(:mr) = q(:mr, j) + q(:mr, jm2) + q(:mr, jp2)
                        end if

                        if (jr == 1) then
                            q(:mr, j) = ZERO
                        else
                            jm1 = j - i2r
                            jp1 = j + i2r
                            q(:mr, j) = HALF * (q(:mr, j)-q(:mr, jm1)-q(:mr, jp1))
                        end if

                        call genbun_aux%solve_tridiag(jr, 0, mr, a, bb, c, b, tcos, d, w)

                        q(:mr, j) = q(:mr, j) + b(:mr)
                    end do

                    if (nlast + i2r <= n) then
                        nrod = 0
                    else
                        nrod = 1
                    end if

                    if (nlastp /= nlast) cycle loop_188

                end do loop_197
            end do loop_188

        end associate

    end subroutine solve_poisson_staggered_grid

end module module_poistg
!
! Revision history
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
