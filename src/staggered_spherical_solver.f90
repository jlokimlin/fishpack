!
!     file hstssp.f90
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
!     SUBROUTINE hstssp(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd,
!                       elmbda, f, idimf, pertrb, ierror)
!
!
! DIMENSION OF           bda(n), bdb(n), bdc(m), bdd(m), f(idimf, n)
! ARGUMENTS
!
! LATEST REVISION        May 2016
!
! PURPOSE                Solves the standard five-point finite
!                        difference approximation on a staggered grid
!                        to the helmholtz equation in spherical
!                        coordinates and on the surface of the unit
!                        sphere (radius of 1).  the equation is
!
!                          (1/sin(theta))(d/dtheta)(sin(theta)
!                          (du/dtheta)) + (1/sin(theta)**2)
!                          (d/dphi)(du/dphi) + lambda*u = f(theta, phi)
!
!                        where theta is colatitude and phi is
!                        longitude.
!
! USAGE                  call hstssp (a, b, m, mbdcnd, bda, bdb, c, d, n,
!                                     nbdcnd, bdc, bdd, elmbda, f, idimf,
!                                     pertrb, ierror)
!
!
! ARGUMENTS
! ON INPUT
!
!                        a, b
!                          the range of theta (colatitude),
!                          i.e. a <= theta <= b.
!                          a must be less than b and a must be
!                          non-negative.  a and b are in radians.
!                          a = 0 corresponds to the north pole and
!                          b = pi corresponds to the south pole.
!
!
!                            * * *  important  * * *
!
!                          if b is equal to pi, then b must be
!                          computed using the statement
!                            b = pi_mach(dum)
!
!                          this insures that b in the user"s program
!                          is equal to pi in this program which
!                          permits several tests of the input
!                          parameters that otherwise would not be
!                          possible.
!
!                            * * * * * * * * * * * *
!                        m
!                          the number of grid points in the interval
!                          (a, b).  the grid points in the theta
!                          direction are given by
!                          theta(i) = a + (i-HALF)dtheta
!                          for i=1, 2, ..., m where dtheta =(b-a)/m.
!                          m must be greater than 2.
!
!                        mbdcnd
!                          indicates the type of boundary conditions
!                          at theta = a and theta = b.
!
!                          = 1  if the solution is specified at
!                               theta = a and theta = b.
!                               (see note 3 below)
!
!                          = 2  if the solution is specified at
!                               theta = a and the derivative of the
!                               solution with respect to theta is
!                               specified at theta = b
!                               (see notes 2 and 3 below).
!
!                          = 3  if the derivative of the solution
!                             with respect to theta is specified
!                             at theta = a
!                             (see notes 1, 2 below) and theta = b.
!
!                          = 4  if the derivative of the solution
!                               with respect to theta is specified
!                               at theta = a
!                               (see notes 1 and 2 below) and the
!                               solution is specified at theta = b.
!
!                          = 5  if the solution is unspecified at
!                               theta = a = 0 and the solution is
!                               specified at theta = b.
!                               (see note 3 below)
!
!                          = 6  if the solution is unspecified at
!                               theta = a = 0 and the derivative
!                               of the solution with respect to theta
!                               is specified at theta = b
!                               (see note 2 below).
!
!                          = 7  if the solution is specified at
!                               theta = a and the solution is
!                               unspecified at theta = b = pi.
!                               (see note 3 below)
!
!                          = 8  if the derivative of the solution
!                               with respect to theta is specified at
!                               theta = a (see note 1 below)
!                               and the solution is unspecified at
!                               theta = b = pi.
!
!                          = 9  if the solution is unspecified at
!                               theta = a = 0 and theta = b = pi.
!
!                          note 1:
!                          if a = 0, do not use mbdcnd = 3, 4, or 8,
!                          but instead use mbdcnd = 5, 6, or 9.
!
!                          note 2:
!                          if b = pi, do not use mbdcnd = 2, 3, or 6,
!                          but instead use mbdcnd = 7, 8, or 9.
!
!                          note 3:
!                          when the solution is specified at
!                          theta = 0 and/or theta = pi and the other
!                          boundary conditions are combinations
!                          of unspecified, normal derivative, or
!                          periodicity a singular system results.
!                          the unique solution is determined by
!                          extrapolation to the specification of the
!                          solution at either theta = 0 or theta = pi.
!                          but in these cases the right side of the
!                          system  will be perturbed by the constant
!                          pertrb.
!
!                        bda
!                          a one-dimensional array of length n that
!                          specifies the boundary values (if any) of
!                          the solution at theta = a.
!
!                          when mbdcnd = 1, 2, or 7,
!                            bda(j) = u(a, phi(j)) ,      j=1, 2, ..., n.
!
!                          when mbdcnd = 3, 4, or 8,
!                            bda(j) = (d/dtheta)u(a, phi(j)) ,
!                            j=1, 2, ..., n.
!
!                          when mbdcnd has any other value,
!                          bda is a dummy variable.
!
!                        bdb
!                          a one-dimensional array of length n that
!                          specifies the boundary values of the
!                          solution at theta = b.
!
!                          when mbdcnd = 1, 4, or 5,
!                            bdb(j) = u(b, phi(j)) ,       j=1, 2, ..., n.
!
!                          when mbdcnd = 2, 3, or 6,
!                            bdb(j) = (d/dtheta)u(b, phi(j)) ,
!                            j=1, 2, ..., n.
!
!                          when mbdcnd has any other value, bdb is
!                          a dummy variable.
!
!                        c, d
!                          the range of phi (longitude),
!                          i.e. c <= phi <= d.
!                          c must be less than d.  if d-c = 2*pi,
!                          periodic boundary conditions are usually
!                          usually prescribed.
!
!                        n
!                          the number of unknowns in the interval
!                          (c, d).  the unknowns in the phi-direction
!                          are given by phi(j) = c + (j-HALF)dphi,
!                          j=1, 2, ..., n, where dphi = (d-c)/n.
!                          n must be greater than 2.
!
!                        nbdcnd
!                          indicates the type of boundary conditions
!                          at phi = c  and phi = d.
!
!                          = 0  if the solution is periodic in phi,
!                               i.e.  u(i, j) = u(i, n+j).
!
!                          = 1  if the solution is specified at
!                               phi = c and phi = d
!                               (see note below).
!
!                          = 2  if the solution is specified at
!                               phi = c and the derivative of the
!                               solution with respect to phi is
!                               specified at phi = d
!                               (see note below).
!
!                          = 3  if the derivative of the solution
!                               with respect to phi is specified
!                               at phi = c and phi = d.
!
!                          = 4  if the derivative of the solution
!                               with respect to phi is specified
!                               at phi = c and the solution is
!                               specified at phi = d
!                               (see note below).
!
!                          note:
!                          when nbdcnd = 1, 2, or 4, do not use
!                          mbdcnd = 5, 6, 7, 8, or 9
!                          (the former indicates that the solution
!                          is specified at a pole; the latter
!                          indicates the solution is unspecified).
!                          use instead mbdcnd = 1 or 2.
!
!                        bdc
!                          a one dimensional array of length m that
!                          specifies the boundary values of the
!                          solution at phi = c.
!
!                          when nbdcnd = 1 or 2,
!                            bdc(i) = u(theta(i), c) ,     i=1, 2, ..., m.
!
!                          when nbdcnd = 3 or 4,
!                            bdc(i) = (d/dphi)u(theta(i), c),
!                            i=1, 2, ..., m.
!
!                          when nbdcnd = 0, bdc is a dummy variable.
!
!                        bdd
!                          a one-dimensional array of length m that
!                          specifies the boundary values of the
!                          solution at phi = d.
!
!                          when nbdcnd = 1 or 4,
!                            bdd(i) = u(theta(i), d) ,     i=1, 2, ..., m.
!
!                          when nbdcnd = 2 or 3,
!                            bdd(i) = (d/dphi)u(theta(i), d) ,
!                            i=1, 2, ..., m.
!
!                          when nbdcnd = 0, bdd is a dummy variable.
!
!                        elmbda
!                          the constant lambda in the helmholtz
!                          equation.  if lambda is greater than 0,
!                          a solution may not exist.  however,
!                          hstssp will attempt to find a solution.
!
!                        f
!                          a two-dimensional array that specifies
!                          the values of the right side of the
!                          helmholtz equation.
!                          for i=1, 2, ..., m and j=1, 2, ..., n
!
!                            f(i, j) = f(theta(i), phi(j)) .
!
!                          f must be dimensioned at least m x n.
!
!                        idimf
!                          the row (or first) dimension of the array
!                          f as it appears in the program calling
!                          hstssp.  this parameter is used to specify
!                          the variable dimension of f.
!                          idimf must be at least m.
!
!
! ON OUTPUT              f
!                          contains the solution u(i, j) of the finite
!                          difference approximation for the grid point
!                          (theta(i), phi(j)) for
!                          i=1, 2, ..., m, j=1, 2, ..., n.
!
!                        pertrb
!                          if a combination of periodic, derivative,
!                          or unspecified boundary conditions is
!                          specified for a poisson equation
!                          (lambda = 0), a solution may not exist.
!                          pertrb is a constant, calculated and
!                          subtracted from f, which ensures that a
!                          solution exists.  hstssp then computes
!                          this solution, which is a least squares
!                          solution to the original approximation.
!                          this solution plus any constant is also
!                          a solution; hence, the solution is not
!                          unique.  the value of pertrb should be
!                          small compared to the right side f.
!                          otherwise, a solution is obtained to an
!                          essentially different problem.
!                          this comparison should always be made to
!                          insure that a meaningful solution has been
!                          obtained.
!
!                        ierror
!                          an error flag that indicates invalid input
!                          parameters. except to numbers 0 and 14,
!                          a solution is not attempted.
!
!                          =  0  no error
!
!                          =  1  a < 0 or b > pi
!
!                          =  2  a >= b
!
!                          =  3  mbdcnd < 1 or mbdcnd > 9
!
!                          =  4  c >= d
!
!                          =  5  n <= 2
!
!                          =  6  nbdcnd < 0 or nbdcnd > 4
!
!                          =  7  a > 0 and mbdcnd = 5, 6, or 9
!
!                          =  8  a = 0 and mbdcnd = 3, 4, or 8
!
!                          =  9  b < pi and mbdcnd >= 7
!
!                          = 10  b = pi and mbdcnd = 2, 3, or 6
!
!                          = 11  mbdcnd >= 5 and ndbcnd = 1, 2, or 4
!
!                          = 12  idimf < m
!
!                          = 13  m <= 2
!
!                          = 14  lambda > 0
!
!                          = 20 if the dynamic allocation of real and
!                               complex workspace required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
!
!                          since this is the only means of indicating
!                          a possibly incorrect call to hstssp, the
!                          user should test ierror after the call.
!
! I/O                    None
!
! PRECISION              64-bit double precision
!
! REQUIRED FILES         type_FishpackWorkspace.f90, genbun.f90, type_CyclicReductionUtility.f9090, poistg.f90
!
! HISTORY                * Written by Roland Sweet at NCAR in 1977.
!                          released on NCAR's public software libraries
!                          in January 1980.
!                        * Revised in June 2004 by John Adams using
!                          Fortran 90 dynamically allocated workspace.
!
! PORTABILITY            Fortran 2008
!
! ALGORITHM              This subroutine defines the finite-
!                        difference equations, incorporates boundary
!                        data, adjusts the right side when the system
!                        is singular and calls either poistg or genbun
!                        which solves the linear system of equations.
!
! TIMING                 For large m and n, the operation count
!                        is roughly proportional to m*n*log2(n).
!
! ACCURACY               The solution process employed results in
!                        a loss of no more than four significant
!                        digits for n and m as large as 64.
!                        more detailed information about accuracy
!                        can be found in the documentation for
!                        routine poistg which is the routine that
!                        actually solves the finite difference
!                        equations.
!
! REFERENCES             U. Schumann and R. Sweet, "A direct method
!                        for the solution of poisson's equation with
!                        neumann boundary conditions on a staggered
!                        grid of arbitrary size, " J. Comp. Phys.
!                        20(1976), pp. 171-182.
!
submodule(staggered_helmholtz_solvers) staggered_spherical_solver

contains

    module subroutine hstssp(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror)

        ! Dummy arguments

        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: mbdcnd
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: nbdcnd
        integer(ip), intent(in)     :: idimf
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: a
        real(wp),    intent(in)     :: b
        real(wp),    intent(in)     :: c
        real(wp),    intent(in)     :: d
        real(wp),    intent(in)     :: elmbda
        real(wp),    intent(out)    :: pertrb
        real(wp),    intent(in)     :: bda(:)
        real(wp),    intent(in)     :: bdb(:)
        real(wp),    intent(in)     :: bdc(:)
        real(wp),    intent(in)     :: bdd(:)
        real(wp),    intent(inout) :: f(:,:)

        ! Local variables

        type(FishpackWorkspace) :: workspace


        ! Check input arguments
        call hstssp_check_input_arguments(a, b, m, mbdcnd, c, d, n, nbdcnd, idimf, ierror)

        if (ierror /= 0) return

        ! Allocate memory
        call workspace%initialize_staggered_workspace(n, m)

        ! Solve system
        associate( rew => workspace%real_workspace )
            call hstssp_lower_routine(a, b, m, mbdcnd, bda, bdb, c, d, n, &
                nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hstssp

    pure subroutine hstssp_check_input_arguments(a, b, m, mbdcnd, c, d, n, nbdcnd, idimf, ierror)

        ! Dummy arguments

        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: mbdcnd
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: nbdcnd
        integer(ip), intent(in)     :: idimf
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: a
        real(wp),    intent(in)     :: b
        real(wp),    intent(in)     :: c
        real(wp),    intent(in)     :: d


        ! Check input arguments
        if (a < ZERO .or. b > PI) then
            ierror = 1
            return
        else if (a >= b) then
            ierror = 2
            return
        else if (mbdcnd <= 0 .or. mbdcnd > 9) then
            ierror = 3
            return
        else if (c >= d) then
            ierror = 4
            return
        else if (3 > n) then
            ierror = 5
            return
        else if (nbdcnd < 0 .or. nbdcnd >= 5) then
            ierror = 6
            return
        else if (a > ZERO) then
            select case (mbdcnd)
                case (5:6,9)
                    ierror = 7
                    return
            end select
        else if (a == ZERO) then
            select case (mbdcnd)
                case (3:4, 8)
                    ierror=8
                    return
            end select
        else if (b < PI .and. mbdcnd >= 7) then
            ierror = 9
        else if (b == PI) then
            select case (mbdcnd)
                case (2:3, 6)
                    ierror = 10
                    return
            end select
        else if (mbdcnd >= 5) then
            select case (nbdcnd)
                case (1:2, 4)
                    ierror = 11
                    return
            end select
        else if (idimf < m) then
            ierror = 12
            return
        else if (3 > m) then
            ierror = 13
            return
        else
            ierror = 0
        end if

    end subroutine hstssp_check_input_arguments

    subroutine hstssp_lower_routine(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, w)

        ! Dummy arguments
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: mbdcnd
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: nbdcnd
        integer(ip), intent(in)    :: idimf
        integer(ip), intent(out)   :: ierror
        real(wp),    intent(in)    :: a
        real(wp),    intent(in)    :: b
        real(wp),    intent(in)    :: c
        real(wp),    intent(in)    :: d
        real(wp),    intent(in)    :: elmbda
        real(wp),    intent(out)   :: pertrb
        real(wp),    intent(in)    :: bda(:)
        real(wp),    intent(in)    :: bdb(:)
        real(wp),    intent(in)    :: bdc(:)
        real(wp),    intent(in)    :: bdd(:)
        real(wp),    intent(inout) :: f(:,:)
        real(wp),    intent(out)   :: w(:)

        ! Local variables
        integer(ip) :: np, isw, jsw, mb, iwb, iwc, iwr, iws
        integer(ip) :: i, j, mm1, lp, local_error_flag
        real(wp)    :: dr, dr2, dth, dth2, a1, a2, a3
        type(CenteredCyclicReductionUtility)  :: centered_util
        type(StaggeredCyclicReductionUtility) :: staggered_util

        dr = (b - a)/m
        dr2 = dr**2
        dth = (d - c)/n
        dth2 = dth**2
        np = nbdcnd + 1
        isw = 1
        jsw = 1
        mb = mbdcnd

        if (elmbda == ZERO) then
            case_construct: select case (mbdcnd)
                case (1, 5, 7)
                    if (a /= ZERO .or. b /= PI) exit case_construct
                    mb = 9
                    jsw = 2
                case (2)
                    if (a /= ZERO) exit case_construct
                    mb = 6
                    jsw = 2
                case (4)
                    if (b /= PI) exit case_construct
                    mb = 8
                    jsw = 2
            end select case_construct
        end if

        iwb = m
        iwc = iwb + m
        iwr = iwc + m
        iws = iwr + m

        do i = 1, m
            j = iwr + i
            w(j) = sin(a + (real(i, kind=wp) - HALF)*dr)
            w(i) = sin(a + real(i - 1, kind=wp)*dr)/dr2
        end do

        mm1 = m - 1
        w(iwc+1:mm1+iwc) = w(2:mm1+1)
        w(iwb+1:mm1+iwb) = elmbda*w(iwr+1:mm1+iwr) - (w(:mm1)+w(2:mm1+1))
        w(iwr) = sin(b)/dr2
        w(iwc) = elmbda*w(iws) - (w(m)+w(iwr))

        do i = 1, m
            j = iwr + i
            a1 = w(j)
            f(i,:n) = a1*f(i,:n)
        end do

        !
        !     enter boundary data for theta-boundaries.
        !
        select case (mb)
            case (1:2, 7)
                a1 = TWO*w(1)
                w(iwb+1) = w(iwb+1) - w(1)
                f(1,:n) = f(1,:n) - a1*bda(:n)
            case (3:4, 8)
                a1 = dr*w(1)
                w(iwb+1) = w(iwb+1) + w(1)
                f(1,:n) = f(1,:n) + a1*bda(:n)
        end select

        select case (mb)
            case (1, 4:5)
                a1 = TWO*w(iwr)
                w(iwc) = w(iwc) - w(iwr)
                f(m,:n) = f(m,:n) - a1*bdb(:n)
            case (2:3, 6)
                a1 = dr*w(iwr)
                w(iwc) = w(iwc) + w(iwr)
                f(m,:n) = f(m,:n) - a1*bdb(:n)
        end select

        !
        ! Enter boundary data for phi-boundaries.
        !
        a1 = TWO/dth2

        if (np /= 1) then
            select case (np)
                case (2:3)
                    f(:m, 1) = f(:m, 1) - a1*bdc(:m)/w(iwr+1:m+iwr)
                case (4:5)
                    a1 = ONE/dth
                    f(:m, 1) = f(:m, 1) + a1*bdc(:m)/w(iwr+1:m+iwr)
            end select

            a1 = TWO/dth2

            select case (np)
                case (2, 5)
                    f(:m, n) = f(:m, n) - a1*bdd(:m)/w(iwr+1:m+iwr)
                case (3:4)
                    a1 = ONE/dth
                    f(:m, n) = f(:m, n) - a1*bdd(:m)/w(iwr+1:m+iwr)
            end select
        end if

        pertrb = ZERO

        if_construct: if (elmbda >= ZERO) then
            if (elmbda /= ZERO) then
                ierror = 14
            else
                select case (mb)
                    case (1:2, 4:5, 7)
                        exit if_construct
                end select

                select case (np)
                    case (2:3, 5)
                        exit if_construct
                end select

                isw = 2
                pertrb = pertrb + sum(f(:m,1:n))

                a1 = real(n, kind=wp)*(cos(a) - cos(b))/(TWO*sin(HALF*dr))
                pertrb = pertrb/a1
                do i = 1, m
                    j = iwr + i
                    a1 = pertrb*w(j)
                    f(i,:n) = f(i,:n) - a1
                end do
                a2 = sum(f(1,:n))
                a3 = sum(f(m,:n))
                a2 = a2/w(iwr+1)
                a3 = a3/w(iws)
            end if
        end if if_construct

        do i = 1, m
            j = iwr + i
            a1 = dth2*w(j)
            w(i) = a1*w(i)
            j = iwc + i
            w(j) = a1*w(j)
            j = iwb + i
            w(j) = a1*w(j)
            f(i,:n) = a1*f(i,:n)
        end do

        lp = nbdcnd
        w(1) = ZERO
        w(iwr) = ZERO
        !
        ! call poistg or genbun to solve the system of equations.
        !
        local_error_flag = 0

        associate( &
            iw0 => 1, &
            iw1 => iwb + 1, &
            iw2 => iwc + 1, &
            iw3 => iwr + 1 &
            )

            if (nbdcnd /= 0) then

                ! Invoke poistg_lower_routine solver
                call staggered_util%poistg_lower_routine(lp, n, iw0, m, w, w(iw1:), w(iw2:), idimf, f, local_error_flag, w(iw3:))

                ! Check error flag
                if (local_error_flag /= 0) then
                    error stop 'fishpack library: poistg_lower_routine call failed in hstssp_lower_routine'
                end if
            else

                ! Invoke genbun_lower_routine solver
                call centered_util%genbun_lower_routine(lp, n, iw0, m, w, w(iw1:), w(iw2:), idimf, f, local_error_flag, w(iw3:))

                ! Check error flag
                if (local_error_flag /= 0) then
                    error stop 'fishpack library: genbun_lower_routine call failed in hstssp_lower_routine'
                end if
            end if

        end associate


        if (isw == 2 .and. jsw == 2) then
            if (mb == 8) then
                a1 = sum(f(m,:n))
                a1 = (a1 - dr2*a3/16)/n
                if (nbdcnd == 3) a1 = a1 + (bdd(m)-bdc(m))/(d - c)
                a1 = bdb(1) - a1
            else
                a1 = sum(f(1,:n))
                a1 = (a1 - dr2*a2/16)/n
                if (nbdcnd == 3) a1 = a1 + (bdd(1)-bdc(1))/(d - c)
                a1 = bda(1) - a1
            end if
            f(:m,:n) = f(:m,:n) + a1
        end if

    end subroutine hstssp_lower_routine

end submodule staggered_spherical_solver
