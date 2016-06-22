!
!     file hstcrt.f90
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
!     SUBROUTINE hstcrt (a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd,
!                        elmbda, f, idimf, pertrb, ierror)
!
! DIMENSION OF           bda(n), bdb(n), bdc(m), bdd(m), f(idimf, n)
! ARGUMENTS
!
! LATEST REVISION        May 2016
!
! PURPOSE                 Solves the standard five-point finite
!                         difference approximation to the helmholtz
!                         equation
!                           (d/dx)(du/dx) + (d/dy)(du/dy) + lambda*u
!                           = f(x, y)
!                         on a staggered grid in cartesian coordinates.
!
! USAGE                   call hstcrt (a, b, m, mbdcnd, bda, bdb, c, d
!                                      n, nbdcnd, bdc, bdd, elmbda,
!                                      f, idimf, pertrb, ierror)
!
! ARGUMENTS
! ON INPUT
!
!                        a, b
!                          the range of x, i.e. a <= x <= b.
!                          a must be less than b.
!
!                        m
!                          the number of grid points in the
!                          interval (a, b).  the grid points
!                          in the x-direction are given by
!                          x(i) = a + (i-0.5)dx for i=1, 2, ..., m
!                          where dx =(b-a)/m.  m must be greater
!                          than 2.
!
!                        mbdcnd
!                          indicates the type of boundary conditions
!                          at x = a and x = b.
!
!                          = 0  if the solution is periodic in x,
!                               u(m+i, j) = u(i, j).
!
!                          = 1  if the solution is specified at
!                               x = a and x = b.
!
!                          = 2  if the solution is specified at
!                               x = a and the derivative
!                               of the solution with respect to x
!                               is specified at x = b.
!
!                          = 3  if the derivative of the solution
!                               with respect to x is specified
!                               at x = a  and x = b.
!
!                          = 4  if the derivative of the solution
!                               with respect to x is specified
!                               at x = a  and the solution is
!                               specified at x = b.
!
!                        bda
!                          a one-dimensional array of length n
!                          that specifies the boundary values
!                          (if any) of the solution at x = a.
!
!                          when mbdcnd = 1 or 2,
!                            bda(j) = u(a, y(j)) ,         j=1, 2, ..., n.
!
!                          when mbdcnd = 3 or 4,
!                            bda(j) = (d/dx)u(a, y(j)) ,   j=1, 2, ..., n.
!
!                        bdb
!                          a one-dimensional array of length n
!                          that specifies the boundary values
!                          of the solution at x = b.
!
!                          when mbdcnd = 1 or 4
!                            bdb(j) = u(b, y(j)) ,        j=1, 2, ..., n.
!
!                          when mbdcnd = 2 or 3
!                            bdb(j) = (d/dx)u(b, y(j)) ,  j=1, 2, ..., n.
!
!                        c, d
!                          the range of y, i.e. c <= y <= d.
!                          c must be less than d.
!
!
!                        n
!                          the number of unknowns in the interval
!                          (c, d).  the unknowns in the y-direction
!                          are given by y(j) = c + (j-0.5)dy,
!                          j=1, 2, ..., n, where dy = (d-c)/n.
!                          n must be greater than 2.
!
!                        nbdcnd
!                          indicates the type of boundary conditions
!                          at y = c   and y = d.
!
!
!                          = 0  if the solution is periodic in y, i.e.
!                               u(i, j) = u(i, n+j).
!
!                          = 1  if the solution is specified at y = c
!                               and y = d.
!
!                          = 2  if the solution is specified at y = c
!                               and the derivative of the solution
!                               with respect to y is specified at
!                               y = d.
!
!                          = 3  if the derivative of the solution
!                               with respect to y is specified at
!                               y = c and y = d.
!
!                          = 4  if the derivative of the solution
!                               with respect to y is specified at
!                               y = c and the solution is specified
!                               at y = d.
!
!                        bdc
!                          a one dimensional array of length m that
!                          specifies the boundary values of the
!                          solution at y = c.
!
!                          when nbdcnd = 1 or 2,
!                            bdc(i) = u(x(i), c) ,        i=1, 2, ..., m.
!
!                          when nbdcnd = 3 or 4,
!                            bdc(i) = (d/dy)u(x(i), c),   i=1, 2, ..., m.
!
!                          when nbdcnd = 0, bdc is a dummy variable.
!
!                        bdd
!                          a one-dimensional array of length m that
!                          specifies the boundary values of the
!                          solution at y = d.
!
!                          when nbdcnd = 1 or 4,
!                            bdd(i) = u(x(i), d) ,        i=1, 2, ..., m.
!
!                          when nbdcnd = 2 or 3,
!                            bdd(i) = (d/dy)u(x(i), d) ,  i=1, 2, ..., m.
!
!                          when nbdcnd = 0, bdd is a dummy variable.
!
!                        elmbda
!                          the constant lambda in the helmholtz
!                          equation.  if lambda is greater than 0,
!                          a solution may not exist. however,
!                          hstcrt will  attempt to find a solution.
!
!                        f
!                          a two-dimensional array that specifies
!                          the values of the right side of the
!                          helmholtz equation.  for i=1, 2, ..., m
!                          and j=1, 2, ..., n
!
!                            f(i, j) = f(x(i), y(j)) .
!
!                          f must be dimensioned at least m x n.
!
!                        idimf
!                          the row (or first) dimension of the array
!                          f as it appears in the program calling
!                          hstcrt.  this parameter is used to specify
!                          the variable dimension of f.
!                          idimf must be at least m.
!
!
! ON OUTPUT              f
!                          contains the solution u(i, j) of the finite
!                          difference approximation for the grid point
!                          (x(i), y(j)) for  i=1, 2, ..., m, j=1, 2, ..., n.
!
!                        pertrb
!                          If a combination of periodic or derivative
!                          boundary conditions is specified for a
!                          poisson equation (lambda = 0), a solution
!                          may not exist.  pertrb is a constant,
!                          calculated and subtracted from f, which
!                          ensures that a solution exists.  hstcrt
!                          then computes this solution, which is a
!                          least squares solution to the original
!                          approximation.  this solution plus any
!                          constant is also a solution; hence, the
!                          solution is not unique.  the value of
!                          pertrb should be small compared to the
!                          right side f.  otherwise, a solution is
!                          obtained to an essentially different problem.
!                          this comparison should always be made to
!                          insure that a meaningful solution has been
!                          obtained.
!
!                        ierror
!                          an error flag that indicates invalid input
!                          parameters.  except to numbers 0 and  6,
!                          a solution is not attempted.
!
!                          =  0  no error
!
!                          =  1  a >= b
!
!                          =  2  mbdcnd < 0 or mbdcnd > 4
!
!                          =  3  c >= d
!
!                          =  4  n <= 2
!
!                         =  5  nbdcnd < 0 or nbdcnd > 4
!
!                         =  6  lambda > 0
!
!                         =  7  idimf < m
!
!                         =  8  m <= 2
!
!                         Since this is the only means of indicating
!                         a possibly incorrect call to hstcrt, the
!                         user should test ierror after the call.
!
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
!
!
! I/O                    None
!
! PRECISION              64-bit double precision
!
! REQUIRED LIBRARY       type_FishpackWorkspace.f90, comf.f90, genbun.f90, gnbnaux.f90, poistg.f90
! FILES
!
! LANGUAGE               Fortran
!
! HISTORY                * Written by Roland Sweet at NCAR in 1977.
!                          released on NCAR's public software libraries
!                          in January 1980.
!                        * Revised in June 2004 by John Adams using
!                          Fortran 90 dynamically allocated work space.
!
! PORTABILITY            Fortran 2008
!
! ALGORITHM              This subroutine defines the finite-difference
!                        equations, incorporates boundary data, adjusts
!                        the right side when the system is singular
!                        and calls either poistg or genbun which solves
!                        the linear system of equations.
!
! TIMING                 For large m and n, the operation count
!                        is roughly proportional to m*n*log2(n).
!
! ACCURACY               The solution process employed results in a
!                        loss of no more than four significant digits
!                        for n and m as large as 64.  more detailed
!                        information about accuracy can be found in
!                        the documentation for package poistg which
!                        solves the finite difference equations.
!
! REFERENCES             U. Schumann and R. Sweet, "A direct method
!                        for the solution of Poisson's equation with
!                        boundary conditions on a staggered grid of
!                        arbitrary size, " J. Comp. Phys. 20(1976), 
!                        PP. 171-182.
!
module module_hstcrt

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_genbun, only: &
        genbunn

    use module_poistg, only: &
        poistgg

    ! Explicit typing only
    implicit None

    ! Everything is private unless stated otherwise
    private
    public :: hstcrt


contains


    subroutine hstcrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
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
        real (wp),             intent (out)    :: pertrb
        real (wp), contiguous, intent (in)     :: bda(:)
        real (wp), contiguous, intent (in)     :: bdb(:)
        real (wp), contiguous, intent (in)     :: bdc(:)
        real (wp), contiguous, intent (in)     :: bdd(:)
        real (wp), contiguous, intent (in out) :: f(:,:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        type (FishpackWorkspace) :: workspace
        integer (ip)             :: irwk, icwk
        !-----------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if ((a-b) >= 0.0_wp) then
            ierror = 1
            return
        else if (mbdcnd < 0 .or. mbdcnd > 4) then
            ierror = 2
            return
        else if ((c-d) >= 0.0_wp) then
            ierror = 3
            return
        else if(n <= 2) then
            ierror = 4
            return
        else if (nbdcnd < 0 .or. nbdcnd > 4) then
            ierror = 5
            return
        else if (idimf < m) then
            ierror = 7
            return
        else if (m <= 2) then
            ierror = 8
            return
        else
            ierror = 0
        end if

        ! compute required real work space
        call workspace%get_genbun_workspace_dimensions(n, m, irwk)

        ! Adjust workspace sizes for solver
        irwk = irwk + 3 * m
        icwk = 0

        !
        !==> Allocate memory
        !
        call workspace%create(irwk, icwk, ierror)

        associate( rew => workspace%real_workspace )

            !
            !==> Solver system
            !
            call hstcrtt(a, b, m, mbdcnd, bda, bdb, &
                c, d, n, nbdcnd, bdc, bdd, &
                elmbda, f, idimf, pertrb, ierror, rew)

        end associate

        !
        !==> Release memory
        !
        call workspace%destroy()

    end subroutine hstcrt


    subroutine hstcrtt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, w)
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
        real (wp),    intent (in out) :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: nperod, mperod, np, mp
        integer (ip) :: id2, id3, id4, local_error_flag
        real (wp)    :: dx, twdelx, delxsq, dy
        real (wp)    :: twdely, dy2, twdysq, s, st2
        !-----------------------------------------------

        nperod = nbdcnd

        if (mbdcnd > 0) then
            mperod = 1
        else
            mperod = 0
        end if

        dx = (b - a)/m
        twdelx = 1.0_wp/dx
        delxsq = 2.0_wp/dx**2
        dy = (d - c)/n
        twdely = 1.0_wp/dy
        dy2 = dy**2
        twdysq = 2.0_wp/dy2
        np = nbdcnd + 1
        mp = mbdcnd + 1
        !
        !     define the a, b, c coefficients in w-array.
        !
        id2 = m
        id3 = id2 + m
        id4 = id3 + m
        s = (dy/dx)**2
        st2 = 2.0_wp*s
        w(:m) = s
        w(id2+1:m+id2) = (-st2) + elmbda*dy2
        w(id3+1:m+id3) = s
        !
        !==> Set boundary data for x-boundaries.
        !
        if (mp /= 1) then
            select case (mp)
                case (2:3)
                    f(1, :n) = f(1, :n) - bda(:n)*delxsq
                    w(id2+1) = w(id2+1) - w(1)
                case (4:5)
                    f(1, :n) = f(1, :n) + bda(:n)*twdelx
                    w(id2+1) = w(id2+1) + w(1)
            end select

            select case (mp)
                case (2, 5)
                    f(m, :n) = f(m, :n) - bdb(:n)*delxsq
                    w(id3) = w(id3) - w(1)
                case (3:4)
                    f(m, :n) = f(m, :n) - bdb(:n)*twdelx
                    w(id3) = w(id3) + w(1)
            end select
        end if

        if (np /= 1) then
            select case (np)
                case (2:3)
                    f(:m, 1) = f(:m, 1) - bdc(:m)*twdysq
                case (4:5)
                    f(:m, 1) = f(:m, 1) + bdc(:m)*twdely
            end select

            select case (np)
                case (2, 5)
                    f(:m, n) = f(:m, n) - bdd(:m)*twdysq
                case (3:4)
                    f(:m, n) = f(:m, n) - bdd(:m)*twdely
            end select
        end if

        f(:m, :n) = f(:m, :n)*dy2

        if (mperod /= 0) then
            w(1) = 0.0_wp
            w(id4) = 0.0_wp
        end if

        pertrb = 0.0_wp

        if (elmbda >= 0.0_wp) then
            if (elmbda /= 0.0_wp) then
                ierror = 6
                return
            else
                select case (mp)
                    case (1, 4)
                        select case (np)
                            case (1, 4)
                                !
                                !==> For singular problems must adjust data to
                                !    insure that a solution will exist.
                                !
                                s = sum(f(:m, 1:n))

                                pertrb = s/(m*n)
                                f(:m, :n) = f(:m, :n) - pertrb
                                pertrb = pertrb/dy2
                        end select
                end select
            end if
        end if

        !
        !==> Solve system
        !
        local_error_flag = 0

        if (nperod /= 0) then
            call poistgg(nperod, n, mperod, m, w(1), w(id2+1), w(id3+1), &
                idimf, f, local_error_flag, w(id4+1))
        else
            call genbunn(nperod, n, mperod, m, w(1), w(id2+1), w(id3+1), &
                idimf, f, local_error_flag, w(id4+1))
        end if

    end subroutine hstcrtt


end module module_hstcrt
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
