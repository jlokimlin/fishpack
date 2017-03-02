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
!     SUBROUTINE hwscyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd,
!                       elmbda, f, idimf, pertrb, ierror)
!
!
! DIMENSION OF           bda(n), bdb(n), bdc(m), bdd(m), f(idimf, n+1)
! ARGUMENTS
!
! PURPOSE                Solves a finite difference approximation
!                        to the helmholtz equation in cylindrical
!                        coordinates.  this modified helmholtz equation
!
!                          (1/r)(d/dr)(r(du/dr)) + (d/dz)(du/dz)
!
!                          + (lambda/r**2)u = f(r, z)
!
!                        results from the fourier transform of the
!                        three-dimensional poisson equation.
!
! USAGE                  call hwscyl(a, b, m, mbdcnd, bda, bdb, c, d, n,
!                                    nbdcnd, bdc, bdd, elmbda, f, idimf,
!                                    pertrb, ierror, w)
!
! ARGUMENTS
! ON INPUT               a, b
!                          the range of r, i.e., a <= r <= b.
!                          a must be less than b and a must be
!                          non-negative.
!
!                        m
!                          the number of panels into which the
!                          interval (a, b) is subdivided.  hence,
!                          there will be m+1 grid points in the
!                          r-direction given by r(i) = a+(i-1)dr,
!                          for i = 1, 2, ..., m+1, where dr = (b-a)/m
!                          is the panel width. m must be greater
!                          than 3.
!
!                        mbdcnd
!                          indicates the type of boundary conditions
!                          at r = a and r = b.
!
!                          = 1  if the solution is specified at
!                               r = a and r = b.
!                          = 2  if the solution is specified at
!                               r = a and the derivative of the
!                               solution with respect to r is
!                               specified at r = b.
!                          = 3  if the derivative of the solution
!                               with respect to r is specified at
!                               r = a (see note below) and r = b.
!                          = 4  if the derivative of the solution
!                               with respect to r is specified at
!                               r = a (see note below) and the
!                               solution is specified at r = b.
!                          = 5  if the solution is unspecified at
!                               r = a = 0 and the solution is
!                               specified at r = b.
!                          = 6  if the solution is unspecified at
!                               r = a = 0 and the derivative of the
!                               solution with respect to r is specified
!                               at r = b.
!
!                          if a = 0, do not use mbdcnd = 3 or 4,
!                          but instead use mbdcnd = 1, 2, 5, or 6  .
!
!                        bda
!                          a one-dimensional array of length n+1 that
!                          specifies the values of the derivative of
!                          the solution with respect to r at r = a.
!
!                          when mbdcnd = 3 or 4,
!                            bda(j) = (d/dr)u(a, z(j)), j = 1, 2, ..., n+1.
!
!                          when mbdcnd has any other value, bda is
!                          a dummy variable.
!
!                        bdb
!                          a one-dimensional array of length n+1 that
!                          specifies the values of the derivative
!                          of the solution with respect to r at r = b.
!
!                          when mbdcnd = 2, 3, or 6,
!                            bdb(j) = (d/dr)u(b, z(j)), j = 1, 2, ..., n+1.
!
!                          when mbdcnd has any other value, bdb is
!                          a dummy variable.
!
!                        c, d
!                          the range of z, i.e., c <= z <= d.
!                          c must be less than d.
!
!                        n
!                          the number of panels into which the
!                          interval (c, d) is subdivided.  hence,
!                          there will be n+1 grid points in the
!                          z-direction given by z(j) = c+(j-1)dz,
!                          for j = 1, 2, ..., n+1,
!                          where dz = (d-c)/n is the panel width.
!                          n must be greater than 3.
!
!                        nbdcnd
!                          indicates the type of boundary conditions
!                          at z = c and z = d.
!
!                          = 0  if the solution is periodic in z,
!                               i.e., u(i, 1) = u(i, n+1).
!                          = 1  if the solution is specified at
!                               z = c and z = d.
!                          = 2  if the solution is specified at
!                               z = c and the derivative of
!                               the solution with respect to z is
!                               specified at z = d.
!                          = 3  if the derivative of the solution
!                               with respect to z is
!                               specified at z = c and z = d.
!                          = 4  if the derivative of the solution
!                               with respect to z is specified at
!                               z = c and the solution is specified
!                               at z = d.
!
!                        bdc
!                          a one-dimensional array of length m+1 that
!                          specifies the values of the derivative
!                          of the solution with respect to z at z = c.
!
!                          when nbdcnd = 3 or 4,
!                            bdc(i) = (d/dz)u(r(i), c), i = 1, 2, ..., m+1.
!
!                          when nbdcnd has any other value, bdc is
!                          a dummy variable.
!
!                        bdd
!                          a one-dimensional array of length m+1 that
!                          specifies the values of the derivative of
!                          the solution with respect to z at z = d.
!
!                          when nbdcnd = 2 or 3,
!                            bdd(i) = (d/dz)u(r(i), d), i = 1, 2, ..., m+1
!
!                          when nbdcnd has any other value, bdd is
!                          a dummy variable.
!
!                        elmbda
!                          the constant lambda in the helmholtz
!                          equation.  if lambda > 0, a solution
!                          may not exist.  however, hwscyl will
!                          attempt to find a solution.  lambda must
!                          be zero when mbdcnd = 5 or 6  .
!
!                        f
!                          a two-dimensional array, of dimension at
!                          least (m+1)*(n+1), specifying values
!                          of the right side of the helmholtz
!                          equation and boundary data (if any).
!
!                          on the interior, f is defined as follows:
!                          for i = 2, 3, ..., m and j = 2, 3, ..., n
!                          f(i, j) = f(r(i), z(j)).
!
!                          on the boundaries f is defined as follows:
!                          for j = 1, 2, ..., n+1 and i = 1, 2, ..., m+1
!
!                          mbdcnd   f(1, j)            f(m+1, j)
!                          ------   ---------         ---------
!
!                            1      u(a, z(j))         u(b, z(j))
!                            2      u(a, z(j))         f(b, z(j))
!                            3      f(a, z(j))         f(b, z(j))
!                            4      f(a, z(j))         u(b, z(j))
!                            5      f(0, z(j))         u(b, z(j))
!                            6      f(0, z(j))         f(b, z(j))
!
!                          nbdcnd   f(i, 1)            f(i, n+1)
!                          ------   ---------         ---------
!
!                            0      f(r(i), c)         f(r(i), c)
!                            1      u(r(i), c)         u(r(i), d)
!                            2      u(r(i), c)         f(r(i), d)
!                            3      f(r(i), c)         f(r(i), d)
!                            4      f(r(i), c)         u(r(i), d)
!
!                          note:
!                          if the table calls for both the solution
!                          u and the right side f at a corner then
!                          the solution must be specified.
!
!                        idimf
!                          the row (or first) dimension of the array
!                          f as it appears in the program calling
!                          hwscyl.  this parameter is used to specify
!                          the variable dimension of f.  idimf must
!                          be at least m+1  .
!
!
! ON OUTPUT              f
!                          contains the solution u(i, j) of the finite
!                          difference approximation for the grid point
!                          (r(i), z(j)), i =1, 2, ..., m+1, j =1, 2, ..., n+1.
!
!                        pertrb
!                          if one specifies a combination of periodic,
!                          derivative, and unspecified boundary
!                          conditions for a poisson equation
!                          (lambda = 0), a solution may not exist.
!                          pertrb is a constant, calculated and
!                          subtracted from f, which ensures that a
!                          solution exists.  hwscyl then computes
!                          this solution, which is a least squares
!                          solution to the original approximation.
!                          this solution plus any constant is also
!                          a solution.  hence, the solution is not
!                          unique. the value of pertrb should be
!                          small compared to the right side f.
!                          otherwise, a solution is obtained to an
!                          essentially different problem.  this
!                          comparison should always be made to insure
!                          that a meaningful solution has been obtained.
!
!                        ierror
!                          an error flag which indicates invalid input
!                          parameters.  except for numbers 0 and 11,
!                          a solution is not attempted.
!
!                          =  0  no error.
!                          =  1  a < 0  .
!                          =  2  a >= b.
!                          =  3  mbdcnd < 1 or mbdcnd > 6  .
!                          =  4  c >= d.
!                          =  5  n <= 3
!                          =  6  nbdcnd < 0 or nbdcnd > 4  .
!                          =  7  a = 0, mbdcnd = 3 or 4  .
!                          =  8  a > 0, mbdcnd >= 5  .
!                          =  9  a = 0, lambda /= 0, mbdcnd >= 5  .
!                          = 10  idimf < m+1  .
!                          = 11  lambda > 0  .
!                          = 12  m <= 3
!                          = 20 if the dynamic allocation of real and
!                               complex workspace required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
!
!                          since this is the only means of indicating
!                          a possibly incorrect call to hwscyl, the
!                          user should test ierror after the call.
!
!
! HISTORY                Written by Roland Sweet at NCAR in the late
!                        1970's.  Released on NCAR's public software
!                        libraries in January 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated workspace.
!
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
!                          m*n*(log2(n)
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
submodule(centered_helmholtz_solvers) centered_cylindrical_solver

contains

    module subroutine hwscyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror)

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

        ! Local variables
        type(FishpackWorkspace)  :: workspace

        ! Check input arguments
        call hwscyl_check_input_arguments(a, b, m, mbdcnd, c, d, n, &
            nbdcnd, elmbda, idimf, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Allocate memory
        call workspace%initialize_staggered_workspace(n, m)

        ! Solve system
        associate( rew => workspace%real_workspace )
            call hwscyl_lower_routine(a, b, m, mbdcnd, bda, bdb, c, d, n, &
                nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hwscyl

    subroutine hwscyl_check_input_arguments(a, b, m, mbdcnd, c, d, n, &
        nbdcnd, elmbda, idimf, ierror)

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

        if (a < ZERO) then
            ierror = 1
        else if (a >= b) then
            ierror = 2
        else if (mbdcnd <= 0 .or. mbdcnd >= 7) then
            ierror = 3
        else if (c >= d) then
            ierror = 4
        else if (n <= 3) then
            ierror = 5
        else if (nbdcnd <= (-1) .or. nbdcnd >= 5) then
            ierror = 6
        else if (a == ZERO .and. (mbdcnd == 3 .or. mbdcnd == 4)) then
            ierror = 7
        else if (a > ZERO .and. mbdcnd >= 5) then
            ierror = 8
        else if (a == ZERO .and. elmbda /= ZERO .and. mbdcnd >= 5) then
            ierror = 9
        else if (idimf < m + 1) then
            ierror = 10
        else if (m <= 3) then
            ierror = 12
        else
            ierror = 0
        end if

    end subroutine hwscyl_check_input_arguments

    subroutine hwscyl_lower_routine(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
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
        real(wp),    intent(inout) :: f(:,:)
        real(wp),    intent(inout) :: w(:)

        ! Local variables
        integer(ip) :: mp1, np1, np, mstart, mstop, munk, nstart, nstop, nunk
        integer(ip) :: id2, id3, id4, id5, id6, istart
        integer(ip) :: ij, i, j, k, l, nsp1, nstm1
        integer(ip) :: local_error_flag, i1
        real(wp)    :: dr, half_dr, dr2, dth, dth2
        real(wp)    :: a1, r, a2, s, s1, s2
        type(CenteredCyclicReductionUtility) :: util

        mp1 = m + 1
        dr = (b - a)/m
        half_dr = dr/2
        dr2 = dr**2
        np1 = n + 1
        dth = (d - c)/n
        dth2 = dth**2
        np = nbdcnd + 1

        ! Define range of indices i and j for unknowns u(i, j).
        mstart = 2
        mstop = m
        select case (mbdcnd)
            case (2)
                mstop = mp1
            case (3, 6)
                mstart = 1
                mstop = mp1
            case (4:5)
                mstart = 1
        end select

        munk = mstop - mstart + 1
        nstart = 1
        nstop = n
        select case (np)
            case (2)
                nstart = 2
            case (3)
                nstart = 2
                nstop = np1
            case (4)
                nstop = np1
        end select

        nunk = nstop - nstart + 1

        ! Define a, b, c coefficients in w-array.
        id2 = munk
        id3 = id2 + munk
        id4 = id3 + munk
        id5 = id4 + munk
        id6 = id5 + munk
        istart = 1
        a1 = TWO/dr2
        ij = 0

        if (mbdcnd==3 .or. mbdcnd==4) ij = 1

        if (mbdcnd > 4) then
            w(1) = 0.
            w(id2+1) = -TWO*a1
            w(id3+1) = TWO*a1
            istart = 2
            ij = 1
        end if

        do i = istart, munk
            r = a + real(i - ij, kind=wp)*dr
            j = id5 + i
            w(j) = r
            j = id6 + i
            w(j) = ONE/r**2
            w(i) = (r - half_dr)/(r*dr2)
            j = id3 + i
            w(j) = (r + half_dr)/(r*dr2)
            k = id6 + i
            j = id2 + i
            w(j) = (-a1) + elmbda*w(k)
        end do

        select case (mbdcnd)
            case (2)
                w(id2) = a1
            case (3, 6)
                w(id2) = a1
                w(id3+1) = a1*real(istart, kind=wp)
            case (4)
                w(id3+1) = a1*real(istart, kind=wp)
        end select

        select case (mbdcnd)
            case (1:2)
                a1 = w(1)
                f(2, nstart:nstop) = f(2, nstart:nstop) - a1*f(1, nstart:nstop)
            case (3:4)
                a1 = TWO*dr*w(1)
                f(1, nstart:nstop) = f(1, nstart:nstop) + a1*bda(nstart:nstop)
        end select

        select case (mbdcnd)
            case (1, 4:5)
                a1 = w(id4)
                f(m, nstart:nstop) = f(m, nstart:nstop) - a1*f(mp1, nstart:nstop)
            case (2:3, 6)
                a1 = TWO*dr*w(id4)
                f(mp1, nstart:nstop) = f(mp1, nstart:nstop) - a1*bdb(nstart:nstop)
        end select

        ! Enter boundary data for z-boundaries.
        a1 = ONE/dth2
        l = id5 - mstart + 1

        if (np /= 1) then
            select case (np)
                case (2:3)
                    f(mstart:mstop, 2) = f(mstart:mstop, 2) - a1*f(mstart:mstop, 1)
                case (4:5)
                    a1 = 2./dth
                    f(mstart:mstop, 1) = f(mstart:mstop, 1) + a1*bdc(mstart:mstop)
            end select

            a1 = ONE/dth2
            select case (np)
                case (2, 5)
                    f(mstart:mstop, n) = f(mstart:mstop, n) - a1*f(mstart:mstop, np1)
                case (3:4)
                    a1 = 2./dth
                    f(mstart:mstop, np1) = f(mstart:mstop, np1) - a1*bdd(mstart:mstop)
            end select
        end if

        if_block: block
            pertrb = ZERO
            if (elmbda >= ZERO) then
                if (elmbda /= ZERO) then
                    ierror = 11
                    return
                else
                    w(id5+1) = HALF*(w(id5+2)-half_dr)

                    select case (mbdcnd)
                        case (1:2, 4:5)
                            exit if_block
                        case (6)
                            w(id5+1) = HALF*w(id5+1)
                    end select

                    select case (np)
                        case (1)
                            a2 = ONE
                        case (2:3, 5)
                            exit if_block
                        case (4)
                            a2 = TWO
                    end select

                    k = id5 + munk
                    w(k) = HALF*(w(k-1)+half_dr)
                    s = ZERO

                    do i = mstart, mstop
                        s1 = ZERO
                        nsp1 = nstart + 1
                        nstm1 = nstop - 1
                        s1 = sum(f(i, nsp1:nstm1))
                        k = i + l
                        s = s + (a2*s1 + f(i, nstart)+f(i, nstop))*w(k)
                    end do

                    s2 = real(m, kind=wp)*a + (0.75_wp + real((m - 1)*(m + 1), kind=wp))*half_dr

                    if (mbdcnd == 3) s2 = s2 + 0.25_wp*half_dr

                    s1 = (TWO + a2*real(nunk-2, kind=wp))*s2

                    pertrb = s/s1
                    f(mstart:mstop, nstart:nstop) = &
                        f(mstart:mstop, nstart:nstop) - pertrb
                end if
            end if
        end block if_block

        w(:mstop-mstart+1) = w(:mstop-mstart+1)*dth2
        w(id2+1:mstop-mstart+1+id2) = w(id2+1:mstop-mstart+1+id2)*dth2
        w(id3+1:mstop-mstart+1+id3) = w(id3+1:mstop-mstart+1+id3)*dth2
        f(mstart:mstop, nstart:nstop) = f(mstart:mstop, nstart:nstop)*dth2
        w(1) = ZERO
        w(id4) = ZERO

        ! Solve the system of equations.
        local_error_flag = 0
        i1 = 1
        call util%genbun_lower_routine(nbdcnd, nunk, i1, munk, w(1:), w(id2+1:), w(id3+1:), &
            idimf, f(mstart:, nstart:), local_error_flag, w(id4+1:))

        if (local_error_flag /= 0) then
            error stop 'fishpack library: genbun_lower_routine call failed in hwscyl_lower_routine'
        end if

        if (nbdcnd == 0) f(mstart:mstop, np1) = f(mstart:mstop, 1)

    end subroutine hwscyl_lower_routine

end submodule centered_cylindrical_solver
