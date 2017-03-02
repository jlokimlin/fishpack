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

submodule(centered_helmholtz_solvers) centered_cartesian_solver

contains

    !     SUBROUTINE hwscrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd,
    !                        elmbda, f, idimf, pertrb, ierror)
    !
    ! DIMENSION OF           bda(n),      bdb(n),   bdc(m), bdd(m),
    ! ARGUMENTS              f(idimf, n)
    !
    ! PURPOSE                Solves the standard five-point finite
    !                        difference approximation to the helmholtz
    !                        equation in cartesian coordinates.  this
    !                        equation is
    !
    !                          (d/dx)(du/dx) + (d/dy)(du/dy)
    !                          + lambda*u = f(x, y).
    !
    ! USAGE                  call hwscrt (a, b, m, mbdcnd, bda, bdb, c, d, n,
    !                                     nbdcnd, bdc, bdd, elmbda, f, idimf,
    !                                     pertrb, ierror)
    !
    ! ARGUMENTS
    ! ON INPUT               a, b
    !
    !                          the range of x, i.e., a <= x <= b.
    !                          a must be less than b.
    !
    !                        m
    !                          the number of panels into which the
    !                          interval (a, b) is subdivided.
    !                          hence, there will be m+1 grid points
    !                          in the x-direction given by
    !                          x(i) = a+(i-1)dx for i = 1, 2, ..., m+1,
    !                          where dx = (b-a)/m is the panel width.
    !                          m must be greater than 3.
    !
    !                        mbdcnd
    !                          indicates the type of boundary conditions
    !                          at x = a and x = b.
    !
    !                          = 0  if the solution is periodic in x,
    !                               i.e., u(i, j) = u(m+i, j).
    !                          = 1  if the solution is specified at
    !                               x = a and x = b.
    !                          = 2  if the solution is specified at
    !                               x = a and the derivative of the
    !                               solution with respect to x is
    !                               specified at x = b.
    !                          = 3  if the derivative of the solution
    !                               with respect to x is specified at
    !                               at x = a and x = b.
    !                          = 4  if the derivative of the solution
    !                               with respect to x is specified at
    !                               x = a and the solution is specified
    !                               at x = b.
    !
    !                        bda
    !                          a one-dimensional array of length n+1 that
    !                          specifies the values of the derivative
    !                          of the solution with respect to x at x = a.
    !
    !                          when mbdcnd = 3 or 4,
    !
    !                            bda(j) = (d/dx)u(a, y(j)), j = 1, 2, ..., n+1.
    !
    !                          when mbdcnd has any other value, bda is
    !                          a dummy variable.
    !
    !                        bdb
    !                          a one-dimensional array of length n+1
    !                          that specifies the values of the derivative
    !                          of the solution with respect to x at x = b.
    !
    !                          when mbdcnd = 2 or 3,
    !
    !                            bdb(j) = (d/dx)u(b, y(j)), j = 1, 2, ..., n+1
    !
    !                          when mbdcnd has any other value bdb is a
    !                          dummy variable.
    !
    !                        c, d
    !                          the range of y, i.e., c <= y <= d.
    !                          c must be less than d.
    !
    !                        n
    !                          the number of panels into which the
    !                          interval (c, d) is subdivided.  hence,
    !                          there will be n+1 grid points in the
    !                          y-direction given by y(j) = c+(j-1)dy
    !                          for j = 1, 2, ..., n+1, where
    !                          dy = (d-c)/n is the panel width.
    !                          n must be greater than 3.
    !
    !                        nbdcnd
    !                          indicates the type of boundary conditions at
    !                          y = c and y = d.
    !
    !                          = 0  if the solution is periodic in y,
    !                               i.e., u(i, j) = u(i, n+j).
    !                          = 1  if the solution is specified at
    !                               y = c and y = d.
    !                          = 2  if the solution is specified at
    !                               y = c and the derivative of the
    !                               solution with respect to y is
    !                               specified at y = d.
    !                          = 3  if the derivative of the solution
    !                               with respect to y is specified at
    !                               y = c and y = d.
    !                          = 4  if the derivative of the solution
    !                               with respect to y is specified at
    !                               y = c and the solution is specified
    !                               at y = d.
    !
    !                        bdc
    !                          a one-dimensional array of length m+1 that
    !                          specifies the values of the derivative
    !                          of the solution with respect to y at y = c.
    !
    !                          when nbdcnd = 3 or 4,
    !
    !                            bdc(i) = (d/dy)u(x(i), c), i = 1, 2, ..., m+1
    !
    !                          when nbdcnd has any other value, bdc is
    !                          a dummy variable.
    !
    !                        bdd
    !                          a one-dimensional array of length m+1 that
    !                          specifies the values of the derivative
    !                          of the solution with respect to y at y = d.
    !
    !                          when nbdcnd = 2 or 3,
    !
    !                            bdd(i) = (d/dy)u(x(i), d), i = 1, 2, ..., m+1
    !
    !                          when nbdcnd has any other value, bdd is
    !                          a dummy variable.
    !
    !                        elmbda
    !                          the constant lambda in the helmholtz
    !                          equation.  if lambda > 0, a solution
    !                          may not exist.  however, hwscrt will
    !                          attempt to find a solution.
    !
    !                        f
    !                          a two-dimensional array, of dimension at
    !                          least (m+1)*(n+1), specifying values of the
    !                          right side of the helmholtz  equation and
    !                          boundary values (if any).
    !
    !                          on the interior, f is defined as follows:
    !                          for i = 2, 3, ..., m and j = 2, 3, ..., n
    !                          f(i, j) = f(x(i), y(j)).
    !
    !                          on the boundaries, f is defined as follows:
    !                          for j=1, 2, ..., n+1,  i=1, 2, ..., m+1,
    !
    !                          mbdcnd     f(1, j)        f(m+1, j)
    !                          ------     ---------     --------
    !
    !                            0        f(a, y(j))     f(a, y(j))
    !                            1        u(a, y(j))     u(b, y(j))
    !                            2        u(a, y(j))     f(b, y(j))
    !                            3        f(a, y(j))     f(b, y(j))
    !                            4        f(a, y(j))     u(b, y(j))
    !
    !
    !                          nbdcnd     f(i, 1)        f(i, n+1)
    !                          ------     ---------     --------
    !
    !                            0        f(x(i), c)     f(x(i), c)
    !                            1        u(x(i), c)     u(x(i), d)
    !                            2        u(x(i), c)     f(x(i), d)
    !                            3        f(x(i), c)     f(x(i), d)
    !                            4        f(x(i), c)     u(x(i), d)
    !
    !                          note:
    !                          if the table calls for both the solution u
    !                          and the right side f at a corner then the
    !                          solution must be specified.
    !
    !                        idimf
    !                          the row (or first) dimension of the array
    !                          f as it appears in the program calling
    !                          hwscrt.  this parameter is used to specify
    !                          the variable dimension of f.  idimf must
    !                          be at least m+1  .
    !
    !
    ! ON OUTPUT              f
    !                          contains the solution u(i, j) of the finite
    !                          difference approximation for the grid point
    !                          (x(i), y(j)), i = 1, 2, ..., m+1,
    !                          j = 1, 2, ..., n+1  .
    !
    !                        pertrb
    !                          if a combination of periodic or derivative
    !                          boundary conditions is specified for a
    !                          poisson equation (lambda = 0), a solution
    !                          may not exist.  pertrb is a constant,
    !                          calculated and subtracted from f, which
    !                          ensures that a solution exists.  hwscrt
    !                          then computes this solution, which is a
    !                          least squares solution to the original
    !                          approximation.  this solution plus any
    !                          constant is also a solution.  hence, the
    !                          solution is not unique.  the value of
    !                          pertrb should be small compared to the
    !                          right side f.  otherwise, a solution is
    !                          obtained to an essentially different
    !                          problem. this comparison should always
    !                          be made to insure that a meaningful
    !                          solution has been obtained.
    !
    !                        ierror
    !                          an error flag that indicates invalid input
    !                          parameters.  except for numbers 0 and 6,
    !                          a solution is not attempted.
    !
    !                          = 0  no error
    !                          = 1  a >= b
    !                          = 2  mbdcnd < 0 or mbdcnd > 4
    !                          = 3  c >= d
    !                          = 4  n <= 3
    !                          = 5  nbdcnd < 0 or nbdcnd > 4
    !                          = 6  lambda > 0
    !                          = 7  idimf < m+1
    !                          = 8  m <= 3
    !                          = 20 If the dynamic allocation of real and
    !                               complex workspace required for solution
    !                               fails (for example if n, m are too large
    !                               for your computer)
    !
    !                          since this is the only means of indicating
    !                          a possibly incorrect call to hwscrt, the
    !                          user should test ierror after the call.
    !
    !
    ! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
    !                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
    !                        LIBRARIES IN January 1980.
    !                        Revised in June 2004 by John Adams using
    !                        Fortran 90 dynamically allocated workspace.
    !
    !
    ! ALGORITHM              The routine defines the finite difference
    !                        equations, incorporates boundary data, and
    !                        adjusts the right side of singular systems
    !                        and then calls genbun to solve the system.
    !
    ! TIMING                 For large  m and n, the operation count
    !                        is roughly proportional to
    !
    !                        m*n*log2(n)
    !
    !                        but also depends on input parameters nbdcnd
    !                        and mbdcnd.
    !
    ! ACCURACY               The solution process employed results in a loss
    !                        of no more than three significant digits for n
    !                        and m as large as 64.  more details about
    !                        accuracy can be found in the documentation for
    !                        subroutine genbun which is the routine that
    !                        solves the finite difference equations.
    !
    ! REFERENCES             Swarztrauber, P. and R. Sweet, "Efficient
    !                        FORTRAN subprograms for the solution of
    !                        elliptic equations"
    !                        NCAR TN/IA-109, July, 1975, 138 pp.
    !
    module subroutine hwscrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
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
        real(wp),    intent(inout)  :: f(:,:)

        ! Local variables
        type(FishpackWorkspace)  :: workspace

        ! Check input arguments
        call hwscrt_check_input_arguments(a, b, m, mbdcnd, c, d, n, nbdcnd, idimf, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Allocate memory
        call workspace%initialize_centered_workspace(n, m)

        ! Solve system
        associate( rew => workspace%real_workspace)
            call hwscrt_lower_routine(a, b, m, mbdcnd, bda, bdb, c, d, &
                n, nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hwscrt

    subroutine hwscrt_lower_routine(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, w)

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
        real(wp),    intent(inout)  :: f(idimf,*)
        real(wp),    intent(inout)  :: w(:)

        ! Local variables
        integer(ip) :: nperod, mperod, np, np1, mp, mp1, nstart, nstop, nskip
        integer(ip) :: nunk, mstart, mstop, mskip, munk, id2, id3, id4
        integer(ip) :: local_error_flag
        real(wp)    :: dx, twdelx, delxsq, dy
        real(wp)    :: twdely, delysq, s, two_s
        type(CenteredCyclicReductionUtility) :: util

        nperod = nbdcnd

        if (mbdcnd > 0) then
            mperod = 1
        else
            mperod = 0
        end if

        dx = (b - a)/m
        twdelx = TWO/dx
        delxsq = ONE/dx**2
        dy = (d - c)/n
        twdely = TWO/dy
        delysq = ONE/dy**2
        np = nbdcnd + 1
        np1 = n + 1
        mp = mbdcnd + 1
        mp1 = m + 1
        nstart = 1
        nstop = n
        nskip = 1

        select case (np)
            case (2)
                nstart = 2
            case (3)
                nstart = 2
                nstop = np1
                nskip = 2
            case (4)
                nstop = np1
                nskip = 2
        end select

        nunk = nstop - nstart + 1
        !
        ! Enter boundary data for x-boundaries.
        !
        mstart = 1
        mstop = m
        mskip = 1

        if (mp /= 1) then
            select case (mp)
                case (2)
                    mstart = 2
                    f(2, nstart:nstop) = f(2, nstart:nstop) - f(1, nstart:nstop)*delxsq
                case (3)
                    mstart = 2
                    mstop = mp1
                    mskip = 2
                    f(2, nstart:nstop) = f(2, nstart:nstop) - f(1, nstart:nstop)*delxsq
                case (4)
                    mstop = mp1
                    mskip = 2
                    f(1, nstart:nstop) = f(1, nstart:nstop) + bda(nstart:nstop)*twdelx
                case (5)
                    f(1, nstart:nstop) = f(1, nstart:nstop) + bda(nstart:nstop)*twdelx
            end select

            select case (mskip)
                case default
                    f(m, nstart:nstop) = f(m, nstart:nstop) - f(mp1, nstart:nstop)* &
                        delxsq
                case (2)
                    f(mp1, nstart:nstop) = f(mp1, nstart:nstop) - bdb(nstart:nstop)* &
                        twdelx
            end select
        end if

        munk = mstop - mstart + 1
        !
        ! Enter boundary data for y-boundaries.
        !
        if (np /= 1) then
            select case (np)
                case (2:3)
                    f(mstart:mstop, 2) = f(mstart:mstop, 2) - f(mstart:mstop, 1)*delysq
                case (4:5)
                    f(mstart:mstop, 1) = f(mstart:mstop, 1) + bdc(mstart:mstop)*twdely
            end select

            select case (nskip)
                case default
                    f(mstart:mstop, n) = f(mstart:mstop, n) - f(mstart:mstop, np1)* &
                        delysq
                case (2)
                    f(mstart:mstop, np1) = f(mstart:mstop, np1) - bdd(mstart:mstop)* &
                        twdely
            end select
        end if

        ! Multiply right side by deltay**2.
        delysq = dy**2
        f(mstart:mstop, nstart:nstop) = f(mstart:mstop, nstart:nstop)*delysq

        !     Define the a, b, c coefficients in w-array.
        id2 = munk
        id3 = id2 + munk
        id4 = id3 + munk
        s = delysq*delxsq
        two_s = TWO*s
        w(:munk) = s
        w(id2+1:munk+id2) = (-two_s) + elmbda*delysq
        w(id3+1:munk+id3) = s

        if (mp /= 1) then
            w(1) = ZERO
            w(id4) = ZERO
        end if

        select case (mp)
            case (3)
                w(id2) = two_s
            case (4)
                w(id2) = two_s
                w(id3+1) = two_s
            case (5)
                w(id3+1) = two_s
        end select

        pertrb = ZERO
        if (elmbda >= ZERO) then
            if (elmbda /= ZERO) then
                ierror = 6
            end if
        end if

        ! Set worspace indices
        associate( &
            iw1 => 1, &
            iw2 => id2 + 1, &
            iw3 => id3 + 1, &
            iw4 => id4 + 1 &
            )

            ! Solve system
            call util%genbun_lower_routine(nperod, nunk, mperod, munk, w(iw1:), w(iw2:), w(iw3:), &
                idimf, f(mstart, nstart), local_error_flag, w(iw4:))

            ! Check error flag
            if (local_error_flag /= 0) then
                error stop 'fishpack library: genbun_lower_routine call failed in hwscrt_lower_routine'
            end if
        end associate

        !  Fill in identical values when have periodic boundary conditions.
        if (nbdcnd == 0) f(mstart:mstop, np1) = f(mstart:mstop, 1)

        if (mbdcnd == 0) then
            f(mp1, nstart:nstop) = f(1, nstart:nstop)
            if (nbdcnd == 0) f(mp1, np1) = f(1, np1)
        end if

    end subroutine hwscrt_lower_routine

    pure subroutine hwscrt_check_input_arguments(a, b, m, mbdcnd, c, d, n, nbdcnd, idimf, ierror)

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
        if (ZERO <= (a - b)) then
            ierror = 1
        else if (mbdcnd < 0 .or. mbdcnd > 4) then
            ierror = 2
        else if (ZERO <= (c - d)) then
            ierror = 3
        else if (n <= 3) then
            ierror = 4
        else if (nbdcnd < 0 .or. nbdcnd > 4) then
            ierror = 5
        else if (idimf < m + 1) then
            ierror = 7
        else if (m <= 3) then
            ierror = 8
        else
            ierror = 0
        end if

    end subroutine hwscrt_check_input_arguments

end submodule centered_cartesian_solver
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
