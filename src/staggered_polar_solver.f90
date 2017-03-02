!
!     file hstplr.f90
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
!     SUBROUTINE hstplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd,
!                       elmbda, f, idimf, pertrb, ierror)
!
! DIMENSION OF           bda(n), bdb(n), bdc(m), bdd(m), f(idimf, n)
! ARGUMENTS
!
! LATEST REVISION        May 2016
!
! PURPOSE                Solves the standard five-point finite
!                        difference approximation on a staggered
!                        grid to the helmholtz equation in polar
!                        coordinates.  the equation is
!
!                           (1/r)(d/dr)(r(du/dr)) +
!                           (1/r**2)(d/dtheta)(du/dtheta) +
!                           lambda*u = f(r, theta)
!
! USAGE                  call hstplr(a, b, m, mbdcnd, bda, bdb, c, d, n,
!                                    nbdcnd, bdc, bdd, elmbda, f,
!                                    idimf, pertrb, ierror)
!
! ARGUMENTS
! ON INPUT               a, b
!
!                          the range of r, i.e. a <= r <= b.
!                          a must be less than b and a must be
!                          non-negative.
!
!                        m
!                          the number of grid points in the interval
!                          (a, b).  the grid points in the r-direction
!                          are given by r(i) = a + (i-0.5)dr for
!                          i=1, 2, ..., m where dr =(b-a)/m.
!                          m must be greater than 2.
!
!                        mbdcnd
!                          indicates the type of boundary conditions
!                          at r = a and r = b.
!
!                          = 1  if the solution is specified at r = a
!                               and r = b.
!
!                          = 2  if the solution is specified at r = a
!                               and the derivative of the solution
!                               with respect to r is specified at r = b.
!                               (see note 1 below)
!
!                          = 3  if the derivative of the solution
!                               with respect to r is specified at
!                               r = a (see note 2 below) and r = b.
!
!                          = 4  if the derivative of the solution
!                               with respect to r is specified at
!                               specified at r = a (see note 2 below)
!                               and the solution is specified at r = b.
!
!
!                          = 5  if the solution is unspecified at
!                               r = a = 0 and the solution is
!                               specified at r = b.
!
!                          = 6  if the solution is unspecified at
!                               r = a = 0 and the derivative of the
!                               solution with respect to r is specified
!                               at r = b.
!
!                          note 1:
!                          if a = 0, mbdcnd = 2, and nbdcnd = 0 or 3,
!                          the system of equations to be solved is
!                          singular.  the unique solution is
!                          is determined by extrapolation to the
!                          specification of u(0, theta(1)).
!                          but in this case the right side of the
!                          system will be perturbed by the constant
!                          pertrb.
!
!                          note 2:
!                          if a = 0, do not use mbdcnd = 3 or 4,
!                          but instead use mbdcnd = 1, 2, 5, or 6.
!
!                        bda
!                          a one-dimensional array of length n that
!                          specifies the boundary values (if any) of
!                          the solution at r = a.
!
!                          when mbdcnd = 1 or 2,
!                            bda(j) = u(a, theta(j)) ,     j=1, 2, ..., n.
!
!                          when mbdcnd = 3 or 4,
!                            bda(j) = (d/dr)u(a, theta(j)) ,
!                            j=1, 2, ..., n.
!
!                          when mbdcnd = 5 or 6, bda is a dummy
!                          variable.
!
!                        bdb
!                          a one-dimensional array of length n that
!                          specifies the boundary values of the
!                          solution at r = b.
!
!                          when mbdcnd = 1, 4, or 5,
!                            bdb(j) = u(b, theta(j)) ,     j=1, 2, ..., n.
!
!                          when mbdcnd = 2, 3, or 6,
!                            bdb(j) = (d/dr)u(b, theta(j)) ,
!                            j=1, 2, ..., n.
!
!                        c, d
!                          the range of theta, i.e. c <= theta <= d.
!                          c must be less than d.
!
!                        n
!                          the number of unknowns in the interval
!                          (c, d).  the unknowns in the theta-
!                          direction are given by theta(j) = c +
!                          (j-0.5)dt,   j=1, 2, ..., n, where
!                          dt = (d-c)/n.  n must be greater than 2.
!
!                        nbdcnd
!                          indicates the type of boundary conditions
!                          at theta = c  and theta = d.
!
!                          = 0  if the solution is periodic in theta,
!                               i.e. u(i, j) = u(i, n+j).
!
!                          = 1  if the solution is specified at
!                               theta = c and theta = d
!                               (see note below).
!
!                          = 2  if the solution is specified at
!                               theta = c and the derivative of the
!                               solution with respect to theta is
!                               specified at theta = d
!                               (see note below).
!
!                          = 3  if the derivative of the solution
!                               with respect to theta is specified
!                               at theta = c and theta = d.
!
!                          = 4  if the derivative of the solution
!                               with respect to theta is specified
!                               at theta = c and the solution is
!                               specified at theta = d
!                               (see note below).
!
!                          note:
!                          when nbdcnd = 1, 2, or 4, do not use
!                          mbdcnd = 5 or 6 (the former indicates that
!                          the solution is specified at r =  0; the
!                          latter indicates the solution is unspecified
!                          at r = 0).  use instead mbdcnd = 1 or 2.
!
!                        bdc
!                          a one dimensional array of length m that
!                          specifies the boundary values of the
!                          solution at theta = c.
!
!                          when nbdcnd = 1 or 2,
!                            bdc(i) = u(r(i), c) ,        i=1, 2, ..., m.
!
!                          when nbdcnd = 3 or 4,
!                            bdc(i) = (d/dtheta)u(r(i), c),
!                            i=1, 2, ..., m.
!
!                          when nbdcnd = 0, bdc is a dummy variable.
!
!                        bdd
!                          a one-dimensional array of length m that
!                          specifies the boundary values of the
!                          solution at theta = d.
!
!                          when nbdcnd = 1 or 4,
!                            bdd(i) = u(r(i), d) ,         i=1, 2, ..., m.
!
!                          when nbdcnd = 2 or 3,
!                            bdd(i) =(d/dtheta)u(r(i), d), i=1, 2, ..., m.
!
!                          when nbdcnd = 0, bdd is a dummy variable.
!
!                        elmbda
!                          the constant lambda in the helmholtz
!                          equation.  if lambda is greater than 0,
!                          a solution may not exist.  however, hstplr
!                          will attempt to find a solution.
!
!                        f
!                          a two-dimensional array that specifies the
!                          values of the right side of the helmholtz
!                          equation.
!
!                          for i=1, 2, ..., m and j=1, 2, ..., n
!                            f(i, j) = f(r(i), theta(j)) .
!
!                          f must be dimensioned at least m x n.
!
!                        idimf
!                          the row (or first) dimension of the array
!                          f as it appears in the program calling
!                          hstplr.  this parameter is used to specify
!                          the variable dimension of f.
!                          idimf must be at least m.
!
!
! ON OUTPUT
!
!                        f
!                          contains the solution u(i, j) of the finite
!                          difference approximation for the grid point
!                          (r(i), theta(j)) for i=1, 2, ..., m,
!                          j=1, 2, ..., n.
!
!                        pertrb
!                          if a combination of periodic, derivative,
!                          or unspecified boundary conditions is
!                          specified for a poisson equation
!                          (lambda = 0), a solution may not exist.
!                          pertrb is a constant calculated and
!                          subtracted from f, which ensures that a
!                          solution exists.  hstplr then computes this
!                          solution, which is a least squares solution
!                          to the original approximation.
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
!                          parameters. except to numbers 0 and 11,
!                          a solution is not attempted.
!
!                          =  0  no error
!
!                          =  1  a < 0
!
!                          =  2  a >= b
!
!                          =  3  mbdcnd < 1 or mbdcnd > 6
!
!                          =  4  c >= d
!
!                          =  5  n <= 2
!
!                          =  6  nbdcnd < 0 or nbdcnd > 4
!
!                          =  7  a = 0 and mbdcnd = 3 or 4
!
!                          =  8  a > 0 and mbdcnd >= 5
!
!                          =  9  mbdcnd >= 5 and nbdcnd /= 0 or 3
!
!                          = 10  idimf < m
!
!                          = 11  lambda > 0
!
!                          = 12  m <= 2
!
!                          = 20 if the dynamic allocation of real and
!                               complex workspace required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
!
!                          since this is the only means of indicating
!                          a possibly incorrect call to hstplr, the
!                          user should test ierror after the call.
!
!
! I/O                    None
!
! PRECISION              64-bit double precision
!
! REQUIRED FILES         type_FishpackWorkspace.f90, genbun.f90, type_CyclicReductionUtility.f9090, poistg.f90
!
! STANDARD               Fortran 2008
!
! HISTORY                Written by Roland Sweet at NCAR in 1977.
!                        released on NCAR's public software libraries
!                        IN January 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated workspace.
!
! PORTABILITY            FORTRAN 90
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
submodule(staggered_helmholtz_solvers) staggered_polar_solver

contains

    module subroutine hstplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
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
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        type(FishpackWorkspace)  :: workspace
        !-----------------------------------------------

        ! Check input arguments
        call hstplr_check_input_arguments(a, b, m, mbdcnd, c, d, n, &
            nbdcnd, idimf, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Allocate memory
        call workspace%initialize_staggered_workspace(n, m)

        ! Solve system
        associate( rew => workspace%real_workspace )
            call hstplr_lower_routine(a, b, m, mbdcnd, &
                bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hstplr

    pure subroutine hstplr_check_input_arguments(a, b, m, mbdcnd, c, d, n, &
        nbdcnd, idimf, ierror)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
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
        !-----------------------------------------------

        ! Check validity of calling arguments
        if (a < ZERO) then
            ierror = 1
            return
        else if (a >= b) then
            ierror = 2
            return
        else if (mbdcnd <= 0 .or. mbdcnd >= 7) then
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
        else if (a == ZERO .and. (mbdcnd == 3 .or. mbdcnd == 4)) then
            ierror = 7
            return
        else if (a > ZERO .and. mbdcnd >= 5) then
            ierror = 8
            return
        else if (mbdcnd >= 5 .and. nbdcnd /= 0 .and. nbdcnd /= 3) then
            ierror = 9
            return
        else if (idimf < m) then
            ierror = 10
            return
        else if (3 > m) then
            ierror = 12
            return
        else
            ierror = 0
        end if

    end subroutine hstplr_check_input_arguments

    subroutine hstplr_lower_routine(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, w)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
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
        real(wp),    intent(inout)  :: f(:,:)
        real(wp),    intent(inout)  :: w(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: np, isw, mb, iwb, iwc, iwr
        integer(ip) :: i, j, k, lp, local_error_flag
        real(wp)    :: dr, dr2, dth, dth2, a1, a2
        type(CenteredCyclicReductionUtility)  :: centered_util
        type(StaggeredCyclicReductionUtility) :: staggered_util

        dr = (b - a)/m
        dr2 = dr**2
        dth = (d - c)/n
        dth2 = dth**2
        np = nbdcnd + 1
        isw = 1
        mb = mbdcnd

        if (a == ZERO .and. mbdcnd == 2) mb = 6
        !
        ! define a, b, c coefficients in w-array.
        !
        iwb = m
        iwc = iwb + m
        iwr = iwc + m

        do i = 1, m
            j = iwr + i
            w(j) = a + (real(i, kind=wp) - HALF)*dr
            w(i) = (a + real(i - 1, kind=wp)*dr)/dr2
            k = iwc + i
            w(k) = (a + real(i, kind=wp)*dr)/dr2
            k = iwb + i
            w(k) = (elmbda - TWO/dr2)*w(j)
        end do

        do i = 1, m
            j = iwr + i
            f(i,:n) = w(j)*f(i,:n)
        end do
        !
        ! Enter boundary data for r-boundaries.
        !
        select case (mb)
            case (1:2)
                a1 = TWO*w(1)
                w(iwb+1) = w(iwb+1) - w(1)
                f(1,:n) = f(1,:n) - a1*bda(:n)
            case (3:4)
                a1 = dr*w(1)
                w(iwb+1) = w(iwb+1) + w(1)
                f(1,:n) = f(1,:n) + a1*bda(:n)
        end select

        select case (mb)
            case (1, 4:5)
                a1 = TWO *w(iwr)
                w(iwc) = w(iwc) - w(iwr)
                f(m,:n) = f(m,:n) - a1*bdb(:n)
            case (2:3, 6)
                a1 = dr*w(iwr)
                w(iwc) = w(iwc) + w(iwr)
                f(m,:n) = f(m,:n) - a1*bdb(:n)
        end select

        !
        ! Enter boundary data for theta-boundaries.
        !

        a1 = TWO/dth2
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
                a1 = ONE /dth
                f(:m, n) = f(:m, n) - a1*bdd(:m)/w(iwr+1:m+iwr)
        end select

        pertrb = ZERO
        if (elmbda >= ZERO) then
            if (elmbda /= ZERO) then
                ierror = 11
                return
            else
                select case (mb)
                    case (3, 6)
                        select case (np)
                            case (1, 4)
                                isw = 2
                                do j = 1, n
                                    pertrb = pertrb + sum(f(:m, j))
                                end do
                                pertrb = pertrb/(real(m*n, kind=wp)*HALF*(a + b))
                                do i = 1, m
                                    j = iwr + i
                                    a1 = pertrb*w(j)
                                    f(i,:n) = f(i,:n) - a1
                                end do
                                a2 = sum(f(1,:n))
                                a2 = a2/w(iwr+1)
                        end select
                end select
            end if
        end if

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
        ! To solve the system of equations.
        !
        local_error_flag = 0

        set_arguments: associate( &
            a_arg => w(1:m), &
            b_arg => w(iwb+1:iwb+1+m), &
            c_arg => w(iwc+1:iwc+1+m), &
            w_arg => w(iwr+1:iwr+1+m) &
            )
            if (lp /= 0) then
                call staggered_util%poistg_lower_routine(lp, n, 1, m, a_arg, b_arg, c_arg, idimf, f, local_error_flag, w_arg)
            else
                call centered_util%genbun_lower_routine(lp, n, 1, m, a_arg, b_arg, w_arg, idimf, f, local_error_flag, w_arg)
            end if
        end associate set_arguments

        if (.not.(a /= ZERO .or. mbdcnd /= 2 .or. isw /= 2)) then
            a1 = sum(f(1,:n))
            a1 = (a1 - dr2*a2/16)/n

            if (nbdcnd == 3) a1 = a1 + (bdd(1)-bdc(1))/(d - c)

            a1 = bda(1) - a1
            f(:m,:n) = f(:m,:n) + a1
        end if

    end subroutine hstplr_lower_routine

end submodule staggered_polar_solver
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
