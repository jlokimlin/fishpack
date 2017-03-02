!
!     file hstcyl.f90
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
! SUBROUTINE hstcyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd,
!            elmbda, f, idimf, pertrb, ierror)
!
! DIMENSION OF           bda(n), bdb(n), bdc(m), bdd(m), f(idimf, n)
! ARGUMENTS
!
! LATEST REVISION        April 2016
!
! PURPOSE                Solves the standard five-point finite
!                        difference approximation on a staggered
!                        grid to the modified helmholtz equation
!                        in cylindrical coordinates. this equation
!
!                          (1/r)(d/dr)(r(du/dr)) + (d/dz)(du/dz)
!
!                          + lambda*(1/r**2)*u = f(r, z)
!
!                        Is a two-dimensional modified helmholtz
!                        equation resulting from the fourier transform
!                        of a three-dimensional poisson equation.
!
! USAGE                  call hstcyl(a, b, m, mbdcnd, bda, bdb, c, d, n, &
!                        nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror)
!
! ARGUMENTS
! ON INPUT               a, b
!
!                          the range of r, i.e. a <= r <= b.
!                          a must be less than b and a must be
!                          be non-negative.
!
!                        m
!                          the number of grid points in the interval
!                          (a, b).  the grid points in the r-direction
!                          r-direction are given by
!                          r(i) = a + (i-0.5)dr for i=1, 2, ..., m
!                          where dr =(b-a)/m.
!                          m must be greater than 2.
!
!                        mbdcnd
!                          indicates the type of boundary conditions
!                          at r = a and r = b.
!
!                          = 1  if the solution is specified at r = a
!                               (see note below) and r = b.
!
!                          = 2  If the solution is specified at r = a
!                               (see note below) and the derivative
!                               of the solution with respect to r is
!                               specified at r = b.
!
!                          = 3  If the derivative of the solution
!                               with respect to r is specified at
!                               r = a (see note below) and r = b.
!
!                          = 4  If the derivative of the solution
!                               with respect to r is specified at
!                               r = a (see note below) and the
!                               solution is specified at r = b.
!
!                          = 5  If the solution is unspecified at
!                               r = a = 0 and the solution is
!                               specified at r = b.
!
!                          = 6  If the solution is unspecified at
!                               r = a = 0 and the derivative of the
!                               solution with respect to r is specified
!                               at r = b.
!
!                          note:
!                          If a = 0, do not use mbdcnd = 1, 2, 3, or 4,
!                          but instead use mbdcnd = 5 or 6.
!                          The resulting approximation gives the only
!                          meaningful boundary condition,
!                          i.e. du/dr = 0
!                          (see d. greenspan, 'Introductory numerical
!                          analysis of elliptic boundary value
!                          problems, ' Harper and Row, 1965, Chapter 5.)
!
!                        bda
!                          A one-dimensional array of length n that
!                          specifies the boundary values (if any)
!                          of the solution at r = a.
!
!                          When mbdcnd = 1 or 2,
!                            bda(j) = u(a, z(j)) ,       j=1, 2, ..., n.
!
!                          When mbdcnd = 3 or 4,
!                            bda(j) = (d/dr)u(a, z(j)) ,   j=1, 2, ..., n.
!
!                          When mbdcnd = 5 or 6, bda is a dummy
!                          variable.
!
!                        bdb
!                          A one-dimensional array of length n that
!                          specifies the boundary values of the
!                          solution at r = b.
!
!                          When mbdcnd = 1, 4, or 5,
!                            bdb(j) = u(b, z(j)) ,        j=1, 2, ..., n.
!
!                          When mbdcnd = 2, 3, or 6,
!                            bdb(j) = (d/dr)u(b, z(j)) ,   j=1, 2, ..., n.
!
!                        c, d
!                          The range of z, i.e. c <= z <= d.
!                          c must be less than d.
!
!                        n
!                          The number of unknowns in the interval
!                          (c, d).  the unknowns in the z-direction
!                          are given by z(j) = c + (j-0.5)dz,
!                          j=1, 2, ..., n, where dz = (d-c)/n.
!                          n must be greater than 2.
!
!                        nbdcnd
!                          Indicates the type of boundary conditions
!                          at z = c  and z = d.
!
!                          = 0  If the solution is periodic in z, i.e.
!                               u(i, j) = u(i, n+j).
!
!                          = 1  If the solution is specified at z = c
!                               and z = d.
!
!                          = 2  If the solution is specified at z = c
!                               and the derivative of the solution with
!                               respect to z is specified at z = d.
!
!                          = 3  If the derivative of the solution with
!                               respect to z is specified at z = c
!                               and z = d.
!
!                          = 4  If the derivative of the solution with
!                               respect to z is specified at z = c and
!                               the solution is specified at z = d.
!
!                        bdc
!                          A one dimensional array of length m that
!                          specifies the boundary values of the
!                          solution at z = c.
!
!                          when nbdcnd = 1 or 2,
!                            bdc(i) = u(r(i), c) ,        i=1, 2, ..., m.
!
!                          when nbdcnd = 3 or 4,
!                            bdc(i) = (d/dz)u(r(i), c),    i=1, 2, ..., m.
!
!                          when nbdcnd = 0, bdc is a dummy variable.
!
!                        bdd
!                          A one-dimensional array of length m that
!                          specifies the boundary values of the
!                          solution at z = d.
!
!                          when nbdcnd = 1 or 4,
!                            bdd(i) = u(r(i), d) ,       i=1, 2, ..., m.
!
!                          when nbdcnd = 2 or 3,
!                            bdd(i) = (d/dz)u(r(i), d) ,   i=1, 2, ..., m.
!
!                          when nbdcnd = 0, bdd is a dummy variable.
!
!                        elmbda
!                          The constant lambda in the modified
!                          helmholtz equation. If lambda is greater
!                          than 0, a solution may not exist.
!                          however, hstcyl will attempt to find a
!                          solution.  lambda must be zero when
!                          mbdcnd = 5 or 6.
!
!                        f
!                          A two-dimensional array that specifies
!                          the values of the right side of the
!                          modified helmholtz equation.
!                          for i=1, 2, ..., m   and j=1, 2, ..., n
!                            f(i, j) = f(r(i), z(j)) .
!                          f must be dimensioned at least m x n.
!
!                        idimf
!                          the row (or first) dimension of the array
!                          f as it appears in the program calling
!                          hstcyl.  this parameter is used to specify
!                          the variable dimension of f.  idimf must
!                          be at least m.
!
! ON OUTPUT
!
!                        f
!                          Contains the solution u(i, j) of the finite
!                          difference approximation for the grid point
!                          (r(i), z(j)) for  i=1, 2, ..., m, j=1, 2, ..., n.
!
!                        pertrb
!                          If a combination of periodic, derivative,
!                          or unspecified boundary conditions is
!                          specified for a poisson equation
!                          (lambda = 0), a solution may not exist.
!                          pertrb is a constant, calculated and
!                          subtracted from f, which ensures that a
!                          solution exists.  hstcyl then computes
!                          this solution, which is a least squares
!                          solution to the original approximation.
!                          this solution plus any constant is also
!                          a solution; hence, the solution is not
!                          unique. The value of pertrb should be
!                          small compared to the right side f.
!                          Otherwise, a solution is obtained to an
!                          essentially different problem.
!                          This comparison should always be made to
!                          insure that a meaningful solution has been
!                          obtained.
!
!                        ierror
!                          An error flag that indicates invalid input
!                          parameters. Except to numbers 0 and 11,
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
!                          =  7  a = 0 and mbdcnd = 1, 2, 3, or 4
!
!                          =  8  a > 0 and mbdcnd >= 5
!
!                          =  9  m <= 2
!
!                          = 10  idimf < m
!
!                          = 11  lambda > 0
!
!                          = 12  a=0, mbdcnd >= 5, elmbda /= 0
!
!                          Since this is the only means of indicating
!                          a possibly incorrect call to hstcyl, the
!                          user should test ierror after the call.
!
!                          = 20 If the dynamic allocation of real and
!                               complex workspace required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
! I/O                    None
!
! PRECISION              64-bit precision float and 32-bit precision integer
!
! REQUIRED LIBRARY       type_FishpackWrapper.f90
! FILES                  genbun.f90, type_CyclicReductionUtility.f9090, poistg.f90
!
! STANDARD               Fortran 2008
!
! HISTORY                * Written by Roland Sweet at NCAR in 1977.
!                          released on NCAR's public software libraries
!                          in January 1980.
!                        * Revised in June 2004 by John Adams using
!                          Fortran 90 dynamically allocated workspace.
!                        * Revised in April 2016 by Jon Lo Kim Lin
!                          using object-oriented features of Fortran 2008
!
! ALGORITHM              This subroutine defines the finite-difference
!                        equations, incorporates boundary data, adjusts
!                        the right side when the system is singular and
!                        calls either poistg or genbun which solves the
!                        linear system of equations.
!
! TIMING                 For large m and n, the operation count
!                        is roughly proportional to m*n*log2(n).
!
! ACCURACY               The solution process results in a loss
!                        of no more than four significant digits
!                        for n and m as large as 64.
!                        more detailed information about accuracy
!                        can be found in the documentation for
!                        subroutine poistg which is the routine that
!                        actually solves the finite difference
!                        equations.
!
! REFERENCES             U. Schumann and R. Sweet, "A direct method for
!                        the solution of poisson's equation with neumann
!                        boundary conditions on a staggered grid of
!                        arbitrary size, " J. Comp. Phys. 20(1976), 
!                        pp. 171-182.
!
submodule(staggered_helmholtz_solvers) staggered_cylindrical_solver

contains

    module subroutine hstcyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
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

        ! Allocate memory
        call workspace%initialize_staggered_workspace(n, m)

        ! Solve system
        associate( rew => workspace%real_workspace )
            call hstcyl_lower_routine(a, b, m, mbdcnd, bda, bdb, c, d, n, &
                nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hstcyl

    subroutine hstcyl_lower_routine(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
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
        real(wp),    intent(inout) :: f(:,:)
        real(wp),    intent(inout) :: w(:)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: np, iwb, iwc, iwr, i, j, k, lp, local_error_flag
        real(wp)    :: dr, dr2, dt, dt2, temp
        type(CenteredCyclicReductionUtility)  :: centered_util
        type(StaggeredCyclicReductionUtility) :: staggered_util

        !
        ! Check validity of calling arguments
        !
        call hstcyl_check_input_arguments(a, b, m, mbdcnd, c, d, n, nbdcnd, &
            elmbda, idimf, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Set radial mesh
        dr = (b - a)/m
        dr2 = dr**2

        ! Set polar mesh
        dt = (d - c)/n
        dt2 = dt**2

        np = nbdcnd + 1

        !
        ! Define a, b, c coefficients in w-array.
        !
        iwb = m
        iwc = iwb + m
        iwr = iwc + m
        do i = 1, m
            j = iwr + i
            w(j) = a + (real(i, kind=wp) - HALF)*dr
            w(i) = (a + real(i - 1, kind=wp)*dr)/(dr2*w(j))
            k = iwc + i
            w(k) = (a + real(i, kind=wp)*dr)/(dr2*w(j))
            k = iwb + i
            w(k) = elmbda/w(j)**2 - TWO/dr2
        end do
        !
        ! Enter boundary data for r-boundaries.
        !
        select case (mbdcnd)
            case (1:2)
                temp = TWO *w(1)
                w(iwb+1) = w(iwb+1) - w(1)
                f(1, :n) = f(1, :n) - temp*bda
            case (3:4)
                temp = dr*w(1)
                w(iwb+1) = w(iwb+1) + w(1)
                f(1, :n) = f(1, :n) + temp * bda
        end select

        select case (mbdcnd)
            case (1, 4:5)
                w(iwc) = w(iwc) - w(iwr)
                temp = TWO * w(iwr)
                f(m, :n) = f(m, :n) - temp * bdb
            case (2:3, 6)
                w(iwc) = w(iwc) + w(iwr)
                temp = dr * w(iwr)
                f(m, :n) = f(m, :n) - temp * bdb
        end select

        !
        ! Enter boundary data for theta-boundaries.
        !
        temp = TWO/dt2

        if (n /= 1) then

            select case (np)
                case (2:3)
                    f(:m, 1) = f(:m, 1) - temp*bdc
                case (4:5)
                    temp = ONE/dt
                    f(:m, 1) = f(:m, 1) + temp*bdc
            end select

            temp = TWO/dt2

            select case (np)
                case (2, 5)
                    f(:m, n) = f(:m, n) - temp*bdd
                case (3:4)
                    temp = ONE/dt
                    f(:m, n) = f(:m, n) - temp*bdd
            end select

        end if

        pertrb = ZERO

        if (elmbda >= ZERO ) then
            if (elmbda /= ZERO ) then
                ierror = 11
                return
            else
                select case (mbdcnd)
                    case (3, 6)
                        select case (np)
                            case (1, 4)
                                do i = 1, m
                                    temp = ZERO
                                    temp = sum(f(i, :n))
                                    j = iwr + i
                                    pertrb = pertrb + temp * w(j)
                                end do
                                pertrb = pertrb/(real(m * n, kind=wp) * HALF * (a + b))
                                f(:m, :n) = f(:m, :n) - pertrb
                        end select
                end select
            end if
        end if

        w(:m) = w(:m) * dt2
        w(iwc+1:m+iwc) = w(iwc+1:m+iwc) * dt2
        w(iwb+1:m+iwb) = w(iwb+1:m+iwb) * dt2
        f(:m, :n) = f(:m, :n)*dt2
        lp = nbdcnd
        w(1) = ZERO
        w(iwr) = ZERO

        associate( &
            iw1 => iwb + 1, &
            iw2 => iwc + 1, &
            iw3 => iwr + 1 &
            )
            !
            ! Solve the system of equations.
            !
            select case (nbdcnd)
                case (0)
                    !
                    ! Solve system with call to genbun_lower_routine
                    !
                    call centered_util%genbun_lower_routine(lp, n, 1, m, w, w(iw1:), w(iw2:), idimf, f, local_error_flag, w(iw3:))

                    ! Check error flag
                    if (local_error_flag /= 0) then
                        error stop 'fishpack library: genbun_lower_routine call failed in hstcyl_lower_routine'
                    end if

                case default
                    !
                    ! Solve system with call to poistg_lower_routine
                    !
                    call staggered_util%poistg_lower_routine(lp, n, 1, m, w, w(iw1:), w(iw2:), idimf, f, local_error_flag, w(iw3:))

                    ! Check error flag
                    if (local_error_flag /= 0) then
                        error stop 'fishpack library: poistg call failed in hstcyl_lower_routine'
                    end if
            end select

        end associate

    end subroutine hstcyl_lower_routine

    pure subroutine hstcyl_check_input_arguments(a, b, m, mbdcnd, c, d, n, nbdcnd, &
        elmbda, idimf, ierror)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: mbdcnd
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: nbdcnd
        integer(ip), intent(in)  :: idimf
        integer(ip), intent(out) :: ierror
        real(wp),    intent(in)  :: a
        real(wp),    intent(in)  :: b
        real(wp),    intent(in)  :: c
        real(wp),    intent(in)  :: d
        real(wp),    intent(in)  :: elmbda
        !--------------------------------------------------------------

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
        else if (a == ZERO .and. mbdcnd /= 5 .and. mbdcnd /= 6) then
            ierror = 7
            return
        else if (a > ZERO .and. mbdcnd >= 5) then
            ierror = 8
            return
        else if (3 > m) then
            ierror = 9
            return
        else if (idimf < m) then
            ierror = 10
            return
        else if (a == ZERO .and. mbdcnd >= 5 .and. elmbda /= ZERO) then
            ierror = 12
            return
        else
            ierror = 0
        end if

    end subroutine hstcyl_check_input_arguments

end submodule staggered_cylindrical_solver
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
