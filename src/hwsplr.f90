!
!     file hwsplr.f90
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
!     SUBROUTINE hwsplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd,
!                       elmbda, f, idimf, pertrb, ierror)
!
!
! DIMENSION OF           bda(n), bdb(n), bdc(m), bdd(m), f(idimf, n+1)
! ARGUMENTS
!
! LATEST REVISION        May 2016
!
! PURPOSE                Solves a finite difference approximation to
!                        the helmholtz equation in polar coordinates.
!                        the equation is
!
!                            (1/r)(d/dr)(r(du/dr)) +
!                            (1/r**2)(d/dtheta)(du/dtheta) +
!                            lambda*u = f(r, theta).
!
! USAGE                  call hwsplr(a, b, m, mbdcnd, bda, bdb, c, d, n,
!                                     nbdcnd, bdc, bdd, elmbda, f, idimf,
!                                     pertrb, ierror, w)
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
!                          for i = 1, 2, ..., m+1,
!                          where dr = (b-a)/m is the panel width.
!                          m must be greater than 3.
!
!                        mbdcnd
!                          indicates the type of boundary condition
!                          at r = a and r = b.
!
!                          = 1  if the solution is specified at
!                               r = a and r = b.
!                          = 2  if the solution is specified at
!                               r = a and the derivative of
!                               the solution with respect to r is
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
!                          note:
!                          if a = 0, do not use mbdcnd = 3 or 4, but
!                          instead use mbdcnd = 1, 2, 5, or 6  .
!
!                        bda
!                          a one-dimensional array of length n+1 that
!                          specifies the values of the derivative of
!                          the solution with respect to r at r = a.
!
!                          when mbdcnd = 3 or 4,
!                            bda(j) = (d/dr)u(a, theta(j)),
!                            j = 1, 2, ..., n+1  .
!
!                          when mbdcnd has any other value, bda is
!                          a dummy variable.
!
!                        bdb
!                          a one-dimensional array of length n+1 that
!                          specifies the values of the derivative of
!                          the solution with respect to r at r = b.
!
!                          when mbdcnd = 2, 3, or 6,
!                            bdb(j) = (d/dr)u(b, theta(j)),
!                            j = 1, 2, ..., n+1  .
!
!                          when mbdcnd has any other value, bdb is
!                          a dummy variable.
!
!                        c, d
!                          the range of theta, i.e., c <=
!                          theta <= d.  c must be less than d.
!
!                        n
!                          the number of panels into which the
!                          interval (c, d) is subdivided.  hence,
!                          there will be n+1 grid points in the
!                          theta-direction given by
!                          theta(j) = c+(j-1)dtheta for
!                          j = 1, 2, ..., n+1, where
!                          dtheta = (d-c)/n is the panel width.
!                          n must be greater than 3.
!
!                        nbdcnd
!                          indicates the type of boundary conditions
!                          at theta = c and at theta = d.
!
!                          = 0  if the solution is periodic in theta,
!                               i.e., u(i, j) = u(i, n+j).
!                          = 1  if the solution is specified at
!                               theta = c and theta = d
!                               (see note below).
!                          = 2  if the solution is specified at
!                               theta = c and the derivative of the
!                               solution with respect to theta is
!                               specified at theta = d
!                               (see note below).
!                          = 4  if the derivative of the solution
!                               with respect to theta is specified
!                               at theta = c and the solution is
!                               specified at theta = d
!                               (see note below).
!
!                          note:
!                          when nbdcnd = 1, 2, or 4, do not use
!                          mbdcnd = 5 or 6
!                          (the former indicates that the solution
!                          is specified at r = 0, the latter indicates
!                          the solution is unspecified at r = 0).
!                          use instead mbdcnd = 1 or 2  .
!
!                        bdc
!                          a one-dimensional array of length m+1 that
!                          specifies the values of the derivative
!                          of the solution with respect to theta at
!                          theta = c.  when nbdcnd = 3 or 4,
!
!                            bdc(i) = (d/dtheta)u(r(i), c),
!                            i = 1, 2, ..., m+1  .
!
!                          when nbdcnd has any other value, bdc is
!                          a dummy variable.
!
!                        bdd
!                          a one-dimensional array of length m+1 that
!                          specifies the values of the derivative
!                          of the solution with respect to theta at
!                          theta = d.  when nbdcnd = 2 or 3,
!
!                            bdd(i) = (d/dtheta)u(r(i), d),
!                            i = 1, 2, ..., m+1  .
!
!                          when nbdcnd has any other value, bdd is
!                          a dummy variable.
!
!                        elmbda
!                          the constant lambda in the helmholtz
!                          equation.  if lambda < 0, a solution
!                          may not exist.  however, hwsplr will
!                          attempt to find a solution.
!
!                        f
!                          a two-dimensional array, of dimension at
!                          least (m+1)*(n+1), specifying values
!                          of the right side of the helmholtz
!                          equation and boundary data (if any).
!
!                          on the interior, f is defined as follows:
!                          for i = 2, 3, ..., m and j = 2, 3, ..., n
!                          f(i, j) = f(r(i), theta(j)).
!
!                          on the boundaries f is defined as follows:
!                          for j = 1, 2, ..., n+1 and i = 1, 2, ..., m+1
!
!                          mbdcnd   f(1, j)            f(m+1, j)
!                          ------   -------------     -------------
!
!                            1      u(a, theta(j))     u(b, theta(j))
!                            2      u(a, theta(j))     f(b, theta(j))
!                            3      f(a, theta(j))     f(b, theta(j))
!                            4      f(a, theta(j))     u(b, theta(j))
!                            5      f(0, 0)            u(b, theta(j))
!                            6      f(0, 0)            f(b, theta(j))
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
!                          then the solution must be specified.
!
!                        idimf
!                          the row (or first) dimension of the array
!                          f as it appears in the program calling
!                          hwsplr.  this parameter is used to specify
!                          the variable dimension of f.  idimf must
!                          be at least m+1.
!
! ON OUTPUT              f
!                          contains the solution u(i, j) of the finite
!                          difference approximation for the grid point
!                          (r(i), theta(j)),
!                          i = 1, 2, ..., m+1, j = 1, 2, ..., n+1  .
!
!                        pertrb
!                          if a combination of periodic, derivative,
!                          or unspecified boundary conditions is
!                          specified for a poisson equation
!                          (lambda = 0), a solution may not exist.
!                          pertrb is a constant, calculated and
!                          subtracted from f, which ensures that a
!                          solution exists.  hwsplr then computes
!                          this solution, which is a least squares
!                          solution to the original approximation.
!                          this solution plus any constant is also
!                          a solution.  hence, the solution is not
!                          unique.  pertrb should be small compared
!                          to the right side. otherwise, a solution
!                          is obtained to an essentially different
!                          problem.  this comparison should always
!                          be made to insure that a meaningful
!                          solution has been obtained.
!
!                        ierror
!                          An error flag that indicates invalid input
!                          parameters. Except for numbers 0 and 11,
!                          a solution is not attempted.
!
!                          =  0  no error.
!                          =  1  a < 0  .
!                          =  2  a >= b.
!                          =  3  mbdcnd < 1 or mbdcnd > 6  .
!                          =  4  c >= d.
!                          =  5  n <= 3
!                          =  6  nbdcnd < 0 or > 4  .
!                          =  7  a = 0, mbdcnd = 3 or 4  .
!                          =  8  a > 0, mbdcnd >= 5  .
!                          =  9  mbdcnd >= 5, nbdcnd /= 0
!                                and nbdcnd /= 3  .
!                          = 10  idimf < m+1  .
!                          = 11  lambda > 0  .
!                          = 12  m <= 3
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
!
!                          Since this is the only means of indicating
!                          a possibly incorrect call to hwsplr, the
!                          user should test ierror after the call.
!
! SPECIAL CONDITIONS     None
!
! I/O                    None
!
! PRECISION              64-bit double precision
!
! REQUIRED files         type_FishpackWorkspace.f90, genbun.f90, gnbnaux.f90, comf.f90
!
! STANDARD               Fortran 2008
!
! HISTORY                * Written by Roland Sweet at NCAR in the late
!                          1970's. Released on NCAR's public software
!                          libraries in January 1980.
!                        * Revised in June 2004 by John Adams using
!                          Fortran 90 dynamically allocated work space.
!
! ALGORITHM              The routine defines the finite difference
!                        equations, incorporates boundary data, and
!                        adjusts the right side of singular systems
!                        and then calls genbun to solve the system.
!
! TIMING                 For large m and n, the operation count
!                        is roughly proportional to
!
!                          m*n*log2(n)
!
!                        but also depends on input parameters nbdcnd
!                        and mbdcnd.
!
! ACCURACY               The solution process employed results in a loss
!                        of no more than three significant digits for n
!                        and m as large as 64. More details about
!                        accuracy can be found in the documentation for
!                        subroutine genbun which is the routine that
!                        solves the finite difference equations.
!
! REFERENCES             Swarztrauber, P. and R. Sweet, "Efficient
!                        FORTRAN subprograms for the solution of
!                        elliptic equations"
!                          NCAR TN/IA-109, July, 1975, 138 pp.
!
module module_hwsplr

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_FishpackWorkspace, only: &
        Fish => FishpackWorkspace

    use module_genbun, only: &
        genbunn

    ! Explicit typing only
    implicit None

    ! Everything is private unless stated otherwise
    private
    public :: hwsplr


contains


    subroutine hwsplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: mbdcnd
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: nbdcnd
        integer (ip), intent (in)     :: idimf
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in)     :: a
        real (wp),    intent (in)     :: b
        real (wp),    intent (in)     :: c
        real (wp),    intent (in)     :: d
        real (wp),    intent (in)     :: elmbda
        real (wp),    intent (out)    :: pertrb
        real (wp),    intent (in)     :: bda(:)
        real (wp),    intent (in)     :: bdb(:)
        real (wp),    intent (in)     :: bdc(:)
        real (wp),    intent (in)     :: bdd(:)
        real (wp),    intent (in out) :: f(:,:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        type (Fish)  :: workspace
        integer (ip) :: irwk, icwk
        !-----------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if (a < 0.0_wp) then
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
        else if (n <= 3) then
            ierror = 5
            return
        else if (nbdcnd <= -1 .or. nbdcnd >= 5) then
            ierror = 6
            return
        else if (a == 0. .and. (mbdcnd==3 .or. mbdcnd==4)) then
            ierror = 7
            return
        else if (a > 0.0_wp .and. mbdcnd >= 5) then
            ierror = 8
            return
        else if (mbdcnd >= 5 .and. nbdcnd /= 0 .and. nbdcnd /= 3) then
            ierror = 9
            return
        else if (idimf < m + 1) then
            ierror = 10
            return
        else if (m <= 3) then
            ierror = 12
            return
        else
            ierror = 0
        end if

        !
        !==> Compute workspace arrays sizes
        !
        irwk = 4*(n+1)+(m+1)*(13+int(log(real(n+1, kind=wp))/log(2.0_wp)))
        icwk = 0

        !
        !==> Allocate memory
        !
        call workspace%create(irwk, icwk, ierror)

        associate( rew => workspace%real_workspace )
            !
            !==> Solve system
            !
            call hwsplrr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror, rew)

        end associate

        !
        !==> Release memory
        !
        call workspace%destroy()

    end subroutine hwsplr



    subroutine hwsplrr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc &
        , bdd, elmbda, f, idimf, pertrb, ierror, w)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: mbdcnd
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: nbdcnd
        integer (ip), intent (in)     :: idimf
        integer (ip), intent (out)    :: ierror
        real (wp),    intent (in)     :: a
        real (wp),    intent (in)     :: b
        real (wp),    intent (in)     :: c
        real (wp),    intent (in)     :: d
        real (wp),    intent (in)     :: elmbda
        real (wp),    intent (out)    :: pertrb
        real (wp),    intent (in)     :: bda(:)
        real (wp),    intent (in)     :: bdb(:)
        real (wp),    intent (in)     :: bdc(:)
        real (wp),    intent (in)     :: bdd(:)
        real (wp),    intent (in out) :: f(idimf,*)
        real (wp),    intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: mp1, np1, np, mstart, mstop, munk, nstart, nstop, nunk
        integer (ip) :: id2, id3, id4, id5, id6, ij, i
        integer (ip) :: j, l, lp, k, i1, local_error_flag, iip
        real (wp)    :: dr, half_dr, dr2, dt, dt2
        real (wp)    :: a1, r, s2, a2, s, s1, ypole
        !-----------------------------------------------

        mp1 = m + 1
        dr = (b - a)/m
        half_dr = dr/2
        dr2 = dr**2
        np1 = n + 1
        dt = (d - c)/n
        dt2 = dt**2
        np = nbdcnd + 1
        !
        !==> Define range of indices i and j for unknowns u(i, j).
        !
        mstart = 2
        mstop = mp1

        select case (mbdcnd)
            case (1)
                mstop = m
            case (3)
                mstart = 1
            case (4)
                mstart = 1
                mstop = m
            case (5)
                mstop = m
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
        !
        !==> Define a, b, c coefficients in w-array.
        !
        id2 = munk
        id3 = id2 + munk
        id4 = id3 + munk
        id5 = id4 + munk
        id6 = id5 + munk
        a1 = 2.0_wp/dr2
        ij = 0

        if (mbdcnd==3 .or. mbdcnd==4) ij = 1

        do i = 1, munk
            r = a + real(i - ij, kind=wp)*dr
            j = id5 + i
            w(j) = r
            j = id6 + i
            w(j) = 1.0_wp/r**2
            w(i) = (r - half_dr)/(r*dr2)
            j = id3 + i
            w(j) = (r + half_dr)/(r*dr2)
            j = id2 + i
            w(j) = (-a1) + elmbda
        end do

        select case (mbdcnd)
            case (2, 6)
                w(id2) = a1
            case (3)
                w(id2) = a1
                w(id3+1) = a1
            case (4)
                w(id3+1) = a1
        end select

        select case (mbdcnd)
            case (1:2)
                a1 = w(1)
                f(2, nstart:nstop) = f(2, nstart:nstop) - a1*f(1, nstart:nstop)
            case (3:4)
                a1 = 2.0_wp * dr*w(1)
                f(1, nstart:nstop) = f(1, nstart:nstop) + a1*bda(nstart:nstop)
        end select

        select case (mbdcnd)
            case (1, 4:5)
                a1 = w(id4)
                f(m, nstart:nstop) = f(m, nstart:nstop) - a1*f(mp1, nstart:nstop)
            case (2:3, 6)
                a1 = 2.0_wp * dr*w(id4)
                f(mp1, nstart:nstop) = f(mp1, nstart:nstop) - a1*bdb(nstart:nstop)
        end select

        !
        !==> Enter boundary data for theta-boundaries.
        !
        a1 = 1.0_wp/dt2
        l = id5 - mstart + 1
        lp = id6 - mstart + 1

        if (np /= 1) then
            select case (np)
                case (2:3)
                    f(mstart:mstop, 2) = f(mstart:mstop, 2) - a1*w(mstart+lp:mstop+lp)*f &
                        (mstart:mstop, 1)
                case (4:5)
                    a1 = 2.0_wp/dt
                    f(mstart:mstop, 1) = f(mstart:mstop, 1) + a1*w(mstart+lp:mstop+lp)* &
                        bdc(mstart:mstop)
            end select

            a1 = 1.0_wp/dt2

            select case (np)
                case (2, 5)
                    f(mstart:mstop, n) = f(mstart:mstop, n) &
                        - a1*w(mstart+lp:mstop+lp) * f(mstart:mstop, np1)
                case (3:4)
                    a1 = 2.0_wp/dt
                    f(mstart:mstop, np1) = f(mstart:mstop, np1) &
                        - a1 * w(mstart+lp:mstop+lp) * bdd(mstart:mstop)
            end select
        end if

        if ((mbdcnd >= 5) .and. (nbdcnd == 3)) then
            f(1, 1) = f(1, 1) - (bdd(2)-bdc(2))* 4.0_wp/(real(n, kind=wp)*dt*dr2)
        end if
        !
        !     adjust right side of singular problems to insure existence of a
        !     solution.
        !
        pertrb = 0.0_wp

        if_block: block
            if (elmbda >= 0.0_wp) then
                if (elmbda /= 0.0_wp) then
                    ierror = 11
                    return
                else
                    if (nbdcnd==0 .or. nbdcnd==3) then
                        s2 = 0.0_wp
                        select case (mbdcnd)
                            case (1:2, 4:5)
                                exit if_block
                            case (3)
                                w(id5+1) = 0.5_wp * (w(id5+2)-half_dr)
                                s2 = 0.25_wp * dr
                        end select

                        if (nbdcnd == 0) then
                            a2 = 1.0_wp
                        else
                            a2 = 2.0_wp
                        end if

                        j = id5 + munk
                        w(j) = 0.5_wp * (w(j-1)+half_dr)
                        s = 0.0_wp

                        do i = mstart, mstop
                            s1 = 0.0_wp
                            ij = nstart + 1
                            k = nstop - 1
                            s1 = sum(f(i, ij:k))
                            j = i + l
                            s = s + (a2*s1 + f(i, nstart)+f(i, nstop))*w(j)
                        end do

                        s2=real(m, kind=wp)*a+dr*(real((m-1)*(m+1), kind=wp)*0.5_wp+0.25_wp)+s2
                        s1 = (2.0_wp + a2*real(nunk - 2, kind=wp))*s2

                        if (mbdcnd /= 3) then
                            s2 = (real(n, kind=wp)*a2*dr)/8
                            s = s + f(1, 1)*s2
                            s1 = s1 + s2
                        end if

                        pertrb = s/s1
                        f(mstart:mstop, nstart:nstop) = &
                            f(mstart:mstop, nstart:nstop) - pertrb
                    end if
                end if
            end if
        end block if_block


        do i = mstart, mstop
            k = i - mstart + 1
            j = i + lp
            a1 = dt2/w(j)
            w(k) = a1*w(k)
            j = id2 + k
            w(j) = a1*w(j)
            j = id3 + k
            w(j) = a1*w(j)
            f(i, nstart:nstop) = a1*f(i, nstart:nstop)
        end do

        w(1) = 0.0_wp
        w(id4) = 0.0_wp
        !
        !==> Solve the system of equations.
        !
        i1 = 1
        local_error_flag = 0

        call genbunn(nbdcnd, nunk, i1, munk, w(1), w(id2+1), w(id3+1), &
            idimf, f(mstart, nstart), local_error_flag, w(id4+1))

        if (local_error_flag /= 0) then
            error stop 'fishpack library: genbunn call failed in hwsplrr'
        end if

        select case (mbdcnd)
            case (1:4)
                if (nbdcnd == 0) then
                    f(mstart:mstop, np1) = f(mstart:mstop, 1)
                end if
                return
            case (5)
                j = id5 + munk
                w(j) = w(id2)/w(id3)

                do iip = 3, munk
                    i = munk - iip + 2
                    j = id5 + i
                    lp = id2 + i
                    k = id3 + i
                    w(j) = w(i)/(w(lp)-w(k)*w(j+1))
                end do

                w(id5+1) = -0.5_wp * dt2/(w(id2+1)-w(id3+1)*w(id5+2))

                do i = 2, munk
                    j = id5 + i
                    w(j) = -w(j)*w(j-1)
                end do

                s = 0.0_wp
                s = sum(f(2, nstart:nstop))
                a2 = nunk

                if (nbdcnd /= 0) then
                    s = s - 0.5_wp * (f(2, nstart)+f(2, nstop))
                    a2 = a2 - 1.0_wp
                end if

                ypole = (0.25*dr2*f(1, 1)-s/a2)/(w(id5+1)-1.+elmbda*dr2*0.25)

                do i = mstart, mstop
                    k = l + i
                    f(i, nstart:nstop) = f(i, nstart:nstop) + ypole*w(k)
                end do

                f(1, :np1) = ypole

                if (nbdcnd == 0) f(mstart:mstop, np1) = f(mstart:mstop, 1)

            case (6)
                !
                !==> Adjust the solution as necessary for the problems where a = 0.
                !
                if (elmbda == 0.0_wp) then

                    ypole = 0.0_wp
                    f(1, :np1) = ypole

                    if (nbdcnd == 0) f(mstart:mstop, np1) = f(mstart:mstop, 1)

                end if
        end select

    end subroutine hwsplrr


end module module_hwsplr
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
