module module_hw3crt

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use module_pois3d, only: &
        pois3dd

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hw3crt


contains


    subroutine hw3crt(xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, mbdcnd, &
        bdys, bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, ldimf, &
        mdimf, f, pertrb, ierror)
        !
        !     file hw3crt.f90
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
        !     *                 A Package of Fortran 77 and 90                *
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
        ! SUBROUTINE hw3crt (xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, mbdcnd, bdys,
        !                    bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, ldimf,
        !                    mdimf, f, pertrb, ierror)
        !
        !
        ! DIMENSION OF           bdxs(mdimf, n+1),    bdxf(mdimf, n+1),
        ! ARGUMENTS              bdys(ldimf, n+1),    bdyf(ldimf, n+1),
        !                        bdzs(ldimf, m+1),    bdzf(ldimf, m+1),
        !                        f(ldimf, mdimf, n+1)
        !
        ! LATEST REVISION        April 2016
        !
        ! PURPOSE                Solves the standard five-point finite
        !                        difference approximation to the helmholtz
        !                        equation in cartesian coordinates.  this
        !                        equation is
        !
        !                          (d/dx)(du/dx) + (d/dy)(du/dy) +
        !                          (d/dz)(du/dz) + lambda*u = f(x, y, z).
        !
        ! USAGE                 call hw3crt(xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m,
        !                            mbdcnd, bdys, bdyf, zs, zf, n, nbdcnd,
        !                            bdzs, bdzf, elmbda, ldimf, mdimf, f,
        !                            pertrb, ierror)
        !
        ! ARGUMENTS
        !
        ! ON INPUT               xs, xf
        !
        !                          the range of x, i.e. xs .le. x .le. xf .
        !                          xs must be less than xf.
        !
        !                        l
        !                          the number of panels into which the
        !                          interval (xs, xf) is subdivided.
        !                          hence, there will be l+1 grid points
        !                          in the x-direction given by
        !                          x(i) = xs+(i-1)dx for i=1, 2, ..., l+1,
        !                          where dx = (xf-xs)/l is the panel width.
        !                          l must be at least 5.
        !
        !                        lbdcnd
        !                          indicates the type of boundary conditions
        !                          at x = xs and x = xf.
        !
        !                          = 0  if the solution is periodic in x,
        !                               i.e. u(l+i, j, k) = u(i, j, k).
        !                          = 1  if the solution is specified at
        !                               x = xs and x = xf.
        !                          = 2  if the solution is specified at
        !                               x = xs and the derivative of the
        !                               solution with respect to x is
        !                               specified at x = xf.
        !                          = 3  if the derivative of the solution
        !                               with respect to x is specified at
        !                               x = xs and x = xf.
        !                          = 4  if the derivative of the solution
        !                               with respect to x is specified at
        !                               x = xs and the solution is specified
        !                               at x=xf.
        !
        !                        bdxs
        !                          a two-dimensional array that specifies the
        !                          values of the derivative of the solution
        !                          with respect to x at x = xs.
        !
        !                          when lbdcnd = 3 or 4,
        !
        !                            bdxs(j, k) = (d/dx)u(xs, y(j), z(k)),
        !                            j=1, 2, ..., m+1,      k=1, 2, ..., n+1.
        !
        !                          when lbdcnd has any other value, bdxs
        !                          is a dummy variable. bdxs must be
        !                          dimensioned at least (m+1)*(n+1).
        !
        !                        bdxf
        !                          a two-dimensional array that specifies the
        !                          values of the derivative of the solution
        !                          with respect to x at x = xf.
        !
        !                          when lbdcnd = 2 or 3,
        !
        !                            bdxf(j, k) = (d/dx)u(xf, y(j), z(k)),
        !                            j=1, 2, ..., m+1,      k=1, 2, ..., n+1.
        !
        !                          when lbdcnd has any other value, bdxf is
        !                          a dummy variable.  bdxf must be
        !                          dimensioned at least (m+1)*(n+1).
        !
        !                        ys, yf
        !                          the range of y, i.e. ys .le. y .le. yf.
        !                          ys must be less than yf.
        !
        !                        m
        !                          the number of panels into which the
        !                          interval (ys, yf) is subdivided.
        !                          hence, there will be m+1 grid points in
        !                          the y-direction given by y(j) = ys+(j-1)dy
        !                          for j=1, 2, ..., m+1,
        !                          where dy = (yf-ys)/m is the panel width.
        !                          m must be at least 5.
        !
        !                        mbdcnd
        !                          indicates the type of boundary conditions
        !                          at y = ys and y = yf.
        !
        !                          = 0  if the solution is periodic in y, i.e.
        !                               u(i, m+j, k) = u(i, j, k).
        !                          = 1  if the solution is specified at
        !                               y = ys and y = yf.
        !                          = 2  if the solution is specified at
        !                               y = ys and the derivative of the
        !                               solution with respect to y is
        !                               specified at y = yf.
        !                          = 3  if the derivative of the solution
        !                               with respect to y is specified at
        !                               y = ys and y = yf.
        !                          = 4  if the derivative of the solution
        !                               with respect to y is specified at
        !                               at y = ys and the solution is
        !                               specified at y=yf.
        !
        !                        bdys
        !                          A two-dimensional array that specifies
        !                          the values of the derivative of the
        !                          solution with respect to y at y = ys.
        !
        !                          when mbdcnd = 3 or 4,
        !
        !                            bdys(i, k) = (d/dy)u(x(i), ys, z(k)),
        !                            i=1, 2, ..., l+1,      k=1, 2, ..., n+1.
        !
        !                          when mbdcnd has any other value, bdys
        !                          is a dummy variable. bdys must be
        !                          dimensioned at least (l+1)*(n+1).
        !
        !                        bdyf
        !                          A two-dimensional array that specifies
        !                          the values of the derivative of the
        !                          solution with respect to y at y = yf.
        !
        !                          when mbdcnd = 2 or 3,
        !
        !                            bdyf(i, k) = (d/dy)u(x(i), yf, z(k)),
        !                            i=1, 2, ..., l+1,      k=1, 2, ..., n+1.
        !
        !                          when mbdcnd has any other value, bdyf
        !                          is a dummy variable. bdyf must be
        !                          dimensioned at least (l+1)*(n+1).
        !
        !                        zs, zf
        !                          The range of z, i.e. zs .le. z .le. zf.
        !                          zs must be less than zf.
        !
        !                        n
        !                          The number of panels into which the
        !                          interval (zs, zf) is subdivided.
        !                          hence, there will be n+1 grid points
        !                          in the z-direction given by
        !                          z(k) = zs+(k-1)dz for k=1, 2, ..., n+1,
        !                          where dz = (zf-zs)/n is the panel width.
        !                          n must be at least 5.
        !
        !                        nbdcnd
        !                          Indicates the type of boundary conditions
        !                          at z = zs and z = zf.
        !
        !                          = 0  if the solution is periodic in z, i.e.
        !                               u(i, j, n+k) = u(i, j, k).
        !                          = 1  if the solution is specified at
        !                               z = zs and z = zf.
        !                          = 2  if the solution is specified at
        !                               z = zs and the derivative of the
        !                               solution with respect to z is
        !                               specified at z = zf.
        !                          = 3  if the derivative of the solution
        !                               with respect to z is specified at
        !                               z = zs and z = zf.
        !                          = 4  if the derivative of the solution
        !                               with respect to z is specified at
        !                               z = zs and the solution is specified
        !                               at z=zf.
        !
        !                        bdzs
        !                          A two-dimensional array that specifies
        !                          the values of the derivative of the
        !                          solution with respect to z at z = zs.
        !
        !                          When nbdcnd = 3 or 4,
        !
        !                            bdzs(i, j) = (d/dz)u(x(i), y(j), zs),
        !                            i=1, 2, ..., l+1,      j=1, 2, ..., m+1.
        !
        !                          When nbdcnd has any other value, bdzs
        !                          is a dummy variable. bdzs must be
        !                          dimensioned at least (l+1)*(m+1).
        !
        !                        bdzf
        !                          A two-dimensional array that specifies
        !                          the values of the derivative of the
        !                          solution with respect to z at z = zf.
        !
        !                          when nbdcnd = 2 or 3,
        !
        !                            bdzf(i, j) = (d/dz)u(x(i), y(j), zf),
        !                            i=1, 2, ..., l+1,      j=1, 2, ..., m+1.
        !
        !                          when nbdcnd has any other value, bdzf
        !                          is a dummy variable. bdzf must be
        !                          dimensioned at least (l+1)*(m+1).
        !
        !                        elmbda
        !                          The constant lambda in the helmholtz
        !                          equation. if lambda > 0, a solution
        !                          may not exist.  however, hw3crt will
        !                          attempt to find a solution.
        !
        !                        ldimf
        !                          The row (or first) dimension of the
        !                          arrays f, bdys, bdyf, bdzs, and bdzf as it
        !                          appears in the program calling hw3crt.
        !                          this parameter is used to specify the
        !                          variable dimension of these arrays.
        !                          ldimf must be at least l+1.
        !
        !                        mdimf
        !                          the column (or second) dimension of the
        !                          array f and the row (or first) dimension
        !                          of the arrays bdxs and bdxf as it appears
        !                          in the program calling hw3crt.  this
        !                          parameter is used to specify the variable
        !                          dimension of these arrays.
        !                          mdimf must be at least m+1.
        !
        !                        f
        !                          a three-dimensional array of dimension at
        !                          at least (l+1)*(m+1)*(n+1), specifying the
        !                          values of the right side of the helmholz
        !                          equation and boundary values (if any).
        !
        !                          on the interior, f is defined as follows:
        !                          for i=2, 3, ..., l,  j=2, 3, ..., m,
        !                          and k=2, 3, ..., n
        !                          f(i, j, k) = f(x(i), y(j), z(k)).
        !
        !                          on the boundaries, f is defined as follows:
        !                          for j=1, 2, ..., m+1,  k=1, 2, ..., n+1,
        !                          and i=1, 2, ..., l+1
        !
        !                          lbdcnd      f(1, j, k)         f(l+1, j, k)
        !                          ------   ---------------   ---------------
        !
        !                            0      f(xs, y(j), z(k))   f(xs, y(j), z(k))
        !                            1      u(xs, y(j), z(k))   u(xf, y(j), z(k))
        !                            2      u(xs, y(j), z(k))   f(xf, y(j), z(k))
        !                            3      f(xs, y(j), z(k))   f(xf, y(j), z(k))
        !                            4      f(xs, y(j), z(k))   u(xf, y(j), z(k))
        !
        !                          mbdcnd      f(i, 1, k)         f(i, m+1, k)
        !                          ------   ---------------   ---------------
        !
        !                            0      f(x(i), ys, z(k))   f(x(i), ys, z(k))
        !                            1      u(x(i), ys, z(k))   u(x(i), yf, z(k))
        !                            2      u(x(i), ys, z(k))   f(x(i), yf, z(k))
        !                            3      f(x(i), ys, z(k))   f(x(i), yf, z(k))
        !                            4      f(x(i), ys, z(k))   u(x(i), yf, z(k))
        !
        !                          nbdcnd      f(i, j, 1)         f(i, j, n+1)
        !                          ------   ---------------   ---------------
        !
        !                            0      f(x(i), y(j), zs)   f(x(i), y(j), zs)
        !                            1      u(x(i), y(j), zs)   u(x(i), y(j), zf)
        !                            2      u(x(i), y(j), zs)   f(x(i), y(j), zf)
        !                            3      f(x(i), y(j), zs)   f(x(i), y(j), zf)
        !                            4      f(x(i), y(j), zs)   u(x(i), y(j), zf)
        !
        !                          Note:
        !                          If the table calls for both the solution
        !                          u and the right side f on a boundary,
        !                          then the solution must be specified.
        !
        !
        ! ON OUTPUT              f
        !                          Contains the solution u(i, j, k) of the
        !                          finite difference approximation for the
        !                          grid point (x(i), y(j), z(k)) for
        !                          i=1, 2, ..., l+1, j=1, 2, ..., m+1,
        !                          and k=1, 2, ..., n+1.
        !
        !                        pertrb
        !                          If a combination of periodic or derivative
        !                          boundary conditions is specified for a
        !                          poisson equation (lambda = 0), a solution
        !                          may not exist.  pertrb is a constant,
        !                          calculated and subtracted from f, which
        !                          ensures that a solution exists.  pwscrt
        !                          then computes this solution, which is a
        !                          least squares solution to the original
        !                          approximation. This solution is not
        !                          unique and is unnormalized. The value of
        !                          pertrb should be small compared to the
        !                          the right side f. Otherwise, a solution
        !                          is obtained to an essentially different
        !                          problem. This comparison should always
        !                          be made to insure that a meaningful
        !                          solution has been obtained.
        !
        !                        ierror
        !                          An error flag that indicates invalid input
        !                          parameters.  except for numbers 0 and 12,
        !                          a solution is not attempted.
        !
        !                          =  0  no error
        !                          =  1  xs >= xf
        !                          =  2  l < 5
        !                          =  3  lbdcnd < 0 .or. lbdcnd > 4
        !                          =  4  ys >= yf
        !                          =  5  m < 5
        !                          =  6  mbdcnd < 0 .or. mbdcnd > 4
        !                          =  7  zs >= zf
        !                          =  8  n < 5
        !                          =  9  nbdcnd < 0 .or. nbdcnd > 4
        !                          = 10  ldimf < l+1
        !                          = 11  mdimf < m+1
        !                          = 12  lambda > 0
        !                          = 20  If the dynamic allocation of real and
        !                                complex work space required for solution
        !                                fails (for example if n, m are too large
        !                                for your computer)
        !
        !                          Since this is the only means of indicating
        !                          a possibly incorrect call to hw3crt, the
        !                          user should test ierror after the call.
        !
        ! SPECIAL CONDITIONS     None
        !
        ! I/O                    None
        !
        ! PRECISION              64-bit precision real and 32-bit precision float
        !
        ! REQUIRED Files         type_FishpackWorkspace.f90,
        !                        pois3d.f90, fftpack.f90, comf.f90
        !
        ! STANDARD               Fortran 2008
        !
        ! HISTORY                * Written by Roland Sweet at NCAR in the late
        !                          1970's.
        !                        * Released on ncar's public software
        !                          libraries in January 1980.
        !                        * Revised in June 2004 by John Adams using
        !                          Fortran 90 dynamically allocated work space.
        !                        * Revised in April 2016 by Jon Lo Kim Lin to
        !                          incorporate features of Fortran 2008
        !
        ! ALGORITHM              This subroutine defines the finite difference
        !                        equations, incorporates boundary data, and
        !                        adjusts the right side of singular systems and
        !                        then calls pois3d to solve the system.
        !
        ! TIMING                 For large l, m and n, the operation count
        !                        is roughly proportional to
        !
        !                          l*m*n*(log2(l)+log2(m)+5),
        !
        !                        but also depends on input parameters lbdcnd
        !                        and mbdcnd.
        !
        ! ACCURACY               The solution process employed results in
        !                        a loss of no more than four significant
        !                        digits for l, m and n as large as 32.
        !                        more detailed information about accuracy
        !                        can be found in the documentation for
        !                        routine pois3d which is the routine that
        !                        actually solves the finite difference
        !                        equations.
        !
        ! REFERENCES             None
        !
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)      :: l
        integer (ip), intent (in)      :: lbdcnd
        integer (ip), intent (in)      :: m
        integer (ip), intent (in)      :: mbdcnd
        integer (ip), intent (in)      :: n
        integer (ip), intent (in)      :: nbdcnd
        integer (ip), intent (in)      :: ldimf
        integer (ip), intent (in)      :: mdimf
        integer (ip), intent (out)     :: ierror
        real (wp),    intent (in)      :: xs
        real (wp),    intent (in)      :: xf
        real (wp),    intent (in)      :: ys
        real (wp),    intent (in)      :: yf
        real (wp),    intent (in)      :: zs
        real (wp),    intent (in)      :: zf
        real (wp),    intent (in)      :: elmbda
        real (wp),    intent (out)     :: pertrb
        real (wp),    intent (in)      :: bdxs(:,:)
        real (wp),    intent (in)      :: bdxf(:,:)
        real (wp),    intent (in)      :: bdys(:,:)
        real (wp),    intent (in)      :: bdyf(:,:)
        real (wp),    intent (in)      :: bdzs(:,:)
        real (wp),    intent (in)      :: bdzf(:,:)
        real (wp),    intent (in out)  :: f(:,:,:)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        type (FishpackWorkspace) :: workspace
        integer (ip)             :: irwk, icwk
        !-----------------------------------------------

        !
        !==> Check validity of input arguments
        !
        call check_input_arguments(l, lbdcnd, m, mbdcnd, n, nbdcnd, &
            ldimf, mdimf, xs, xf, ys, yf, zs, zf, ierror)

        ! Check error flag
        if (ierror /= 0) then
            return
        end if

        !
        !==> Compute required workspace dimensions for hw3crtt
        !
        irwk = 30+l+m+5*n+max(l, m, n) + 7*((l+1)/2 + (m+1)/2)
        icwk = 0
        call workspace%create(irwk, icwk)

        !
        !==> Solve system
        !
        associate( rew => workspace%real_workspace )

            call hw3crtt(xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, mbdcnd, bdys, &
                bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, ldimf, &
                mdimf, f, pertrb, ierror, rew)

        end associate

        !
        !==> Release memory
        !
        call workspace%destroy()

    end subroutine hw3crt


    subroutine hw3crtt(xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, &
        mbdcnd, bdys, bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, &
        ldimf, mdimf, f, pertrb, ierror, w)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)      :: l
        integer (ip), intent (in)      :: lbdcnd
        integer (ip), intent (in)      :: m
        integer (ip), intent (in)      :: mbdcnd
        integer (ip), intent (in)      :: n
        integer (ip), intent (in)      :: nbdcnd
        integer (ip), intent (in)      :: ldimf
        integer (ip), intent (in)      :: mdimf
        integer (ip), intent (out)     :: ierror
        real (wp),    intent (in)      :: xs
        real (wp),    intent (in)      :: xf
        real (wp),    intent (in)      :: ys
        real (wp),    intent (in)      :: yf
        real (wp),    intent (in)      :: zs
        real (wp),    intent (in)      :: zf
        real (wp),    intent (in)      :: elmbda
        real (wp),    intent (out)     :: pertrb
        real (wp),    intent (in)      :: bdxs(mdimf, *)
        real (wp),    intent (in)      :: bdxf(mdimf, *)
        real (wp),    intent (in)      :: bdys(ldimf, *)
        real (wp),    intent (in)      :: bdyf(ldimf, *)
        real (wp),    intent (in)      :: bdzs(ldimf, *)
        real (wp),    intent (in)      :: bdzf(ldimf, *)
        real (wp),    intent (in out)  :: f(ldimf, mdimf, *)
        real (wp),    intent (in out)  :: w(*)
        !-----------------------------------------------
        ! Dictionary: local variables
        !-----------------------------------------------
        integer (ip) :: exit_counter
        integer (ip) :: mstart, mstop, mp1, mp, munk, np, np1
        integer (ip) :: nstart, nstop, nunk, lp1, lp, lstart
        integer (ip) :: lstop, j, k, lunk, i, iwb, iwc, iww
        integer (ip) :: mstpm1, lstpm1, nstpm1, nperod, ir
        real (wp)    :: dy, twbydy, c2, dz, twbydz, c3, dx
        real (wp)    :: c1, twbydx, xlp, ylp, zlp, s1, s2, s
        !-----------------------------------------------

        dy = (yf - ys)/m
        twbydy = 2.0_wp/dy
        c2 = 1.0_wp/dy**2
        mstart = 1
        mstop = m
        mp1 = m + 1
        mp = mbdcnd + 1

        !
        ! GCC 5.3 does not support Fortran 2008 the exit statement within
        ! a select case construct yet.
        !
        exit_select_case_1: do exit_counter = 1,1
            select case (mp)
                case (1)
                    exit exit_select_case_1
                case (2:3)
                    mstart = 2
                    select case (mp)
                        case (1:2, 5)
                            exit exit_select_case_1
                        case (3:4)
                            mstop = mp1
                    end select
                case (4:5)
                    select case (mp)
                        case (1:2, 5)
                            exit exit_select_case_1
                        case (3:4)
                            mstop = mp1
                    end select
            end select
        end do exit_select_case_1

        munk = mstop - mstart + 1
        dz = (zf - zs)/n
        twbydz = 2.0_wp/dz
        np = nbdcnd + 1
        c3 = 1.0_wp/dz**2
        np1 = n + 1
        nstart = 1
        nstop = n

        !
        ! GCC 5.3 does not support Fortran 2008 the exit statement within
        ! a select case construct yet.
        !
        exit_select_case_2: do exit_counter = 1,1
            select case (np)
                case (1)
                    exit exit_select_case_2
                case (2:3)
                    nstart = 2
                    select case (np)
                        case (2) ! case(3)
                            exit exit_select_case_2
                        case default
                            nstop = np1
                    end select
                case (4:5)
                    select case (np)
                        case (4) ! case(5)
                            exit exit_select_case_2
                        case default
                            nstop = np1
                    end select
            end select
        end do exit_select_case_2

        nunk = nstop - nstart + 1
        lp1 = l + 1
        dx = (xf - xs)/l
        c1 = 1.0_wp/dx**2
        twbydx = 2.0_wp/dx
        lp = lbdcnd + 1
        lstart = 1
        lstop = l
        !
        !==> Enter boundary data for x-boundaries.
        !
        !
        if (lp /= 1) then
            select case (lp)
                case (2:3)
                    lstart = 2
                    f(2, mstart:mstop, nstart:nstop) = &
                        f(2, mstart:mstop, nstart:nstop) &
                        - c1*f(1, mstart:mstop, nstart:nstop)
                case (4:5)
                    f(1, mstart:mstop, nstart:nstop) = &
                        f(1, mstart:mstop, nstart:nstop) &
                        + twbydx*bdxs(mstart:mstop, nstart:nstop)
            end select

            select case (lp)
                case (2, 5)
                    f(l, mstart:mstop, nstart:nstop) = &
                        f(l, mstart:mstop, nstart:nstop) &
                        - c1*f(lp1, mstart:mstop, nstart:nstop)
                case (3:4)
                    lstop = lp1
                    f(lp1, mstart:mstop, nstart:nstop) = &
                        f(lp1, mstart:mstop, nstart:nstop) &
                        - twbydx*bdxf(mstart:mstop, nstart:nstop)
            end select
        end if

        lunk = lstop - lstart + 1

        !
        !==> Enter boundary data for y-boundaries.
        !
        if (mp /= 1) then
            select case (mp)
                case (2:3)
                    f(lstart:lstop, 2, nstart:nstop) = &
                        f(lstart:lstop, 2, nstart:nstop)&
                        - c2*f(lstart:lstop, 1, nstart:nstop)
                case (4:5)
                    f(lstart:lstop, 1, nstart:nstop) = &
                        f(lstart:lstop, 1, nstart:nstop) &
                        + twbydy*bdys(lstart:lstop, nstart:nstop)
            end select

            select case (mp)
                case (2, 5)
                    f(lstart:lstop, m, nstart:nstop) = &
                        f(lstart:lstop, m, nstart:nstop) &
                        - c2*f(lstart:lstop, mp1, nstart:nstop)
                case (3:4)
                    f(lstart:lstop, mp1, nstart:nstop) = &
                        f(lstart:lstop, mp1, nstart:nstop) &
                        - twbydy*bdyf(lstart:lstop, nstart:nstop)
            end select
        end if


        if (np /= 1) then
            select case (np)
                case (2:3)
                    f(lstart:lstop, mstart:mstop, 2) = &
                        f(lstart:lstop, mstart:mstop, 2) &
                        - c3*f(lstart:lstop, mstart:mstop, 1)
                case (4:5)
                    f(lstart:lstop, mstart:mstop, 1) = &
                        f(lstart:lstop, mstart:mstop, 1) &
                        + twbydz*bdzs(lstart:lstop, mstart:mstop)
            end select


            select case (np)
                case (2, 5)
                    f(lstart:lstop, mstart:mstop, n) = &
                        f(lstart:lstop, mstart:mstop, n) &
                        - c3*f(lstart:lstop, mstart:mstop, np1)
                case (3:4)
                    f(lstart:lstop, mstart:mstop, np1) = &
                        f(lstart:lstop, mstart:mstop, np1) &
                        - twbydz*bdzf(lstart:lstop, mstart:mstop)
            end select
        end if

        !
        !==> Define a, b, c coefficients in w-array.
        !
        iwb = nunk + 1
        iwc = iwb + nunk
        iww = iwc + nunk
        w(:nunk) = c3
        w(iwc:nunk-1+iwc) = c3
        w(iwb:nunk-1+iwb) = (-2.*c3) + elmbda

        !
        ! GCC 5.3 does not support Fortran 2008 the exit statement within
        ! a select case construct yet.
        !
        exit_select_case_3: do exit_counter = 1,1
            select case (np)
                case (1:2)
                    exit exit_select_case_3
                case (3)
                    w(iwb-1) = 2.0_wp*c3
                case (4:5)
                    w(iwc) = 2.0_wp*c3
                    select case (np)
                        case (4)
                            w(iwb-1) = 2.0_wp*c3
                        case (5)
                            exit exit_select_case_3
                    end select
            end select
        end do exit_select_case_3

        pertrb = 0.0_wp
        !
        !==> For singular problems adjust data to insure a solution will exist.
        !
        !
        ! GCC 5.3 does not support Fortran 2008 the exit statement within
        ! a select case construct yet.
        !
        exit_select_case_4: do exit_counter = 1,1
            select case (lp)
                case (1, 4)
                    select case (mp)
                        case (1, 4)
                            select case (np)
                                case (1, 4)
                                    if (elmbda >= 0.0_wp) then
                                        if (elmbda /= 0.0_wp) then
                                            ierror = 12
                                            return
                                        else
                                            mstpm1 = mstop - 1
                                            lstpm1 = lstop - 1
                                            nstpm1 = nstop - 1
                                            xlp = (2 + lp)/3
                                            ylp = (2 + mp)/3
                                            zlp = (2 + np)/3
                                            s1 = 0.0_wp

                                            do k = 2, nstpm1
                                                do j = 2, mstpm1
                                                    s1 = s1 + sum(f(2:lstpm1, j, k))
                                                    s1 = s1 + (f(1, j, k)+f(lstop, j, k))/xlp
                                                end do
                                                s2 = sum(f(2:lstpm1, 1, k)+f(2:lstpm1, mstop, k))
                                                s2 = (s2 + (f(1, 1, k) + f(1, mstop, k) &
                                                    +f(lstop, 1, k) + f(lstop,mstop, k))/xlp)/ylp
                                                s1 = s1 + s2
                                            end do

                                            s = (f(1, 1, 1)+f(lstop, 1, 1) &
                                                + f(1, 1, nstop)+f(lstop, 1, nstop) &
                                                + f(1, mstop, 1)+f(lstop, mstop, 1) &
                                                + f(1, mstop, nstop)+f(lstop, mstop, nstop))/(xlp*ylp)

                                            do j = 2, mstpm1
                                                s = s + sum(f(2:lstpm1, j, 1)+f(2:lstpm1, j, nstop))
                                            end do

                                            s2 = 0.0_wp
                                            s2 = sum(f(2:lstpm1, 1, 1)+f(2:lstpm1, 1, nstop) &
                                                + f(2:lstpm1, mstop, 1)+f(2:lstpm1, mstop, nstop))
                                            s = s2/ylp + s
                                            s2 = 0.0_wp
                                            s2 = sum(f(1, 2:mstpm1, 1)+f(1, 2:mstpm1, nstop) &
                                                + f(lstop, 2:mstpm1, 1)+f(lstop, 2:mstpm1, nstop))
                                            s = s2/xlp + s
                                            pertrb = &
                                                (s/zlp + s1)/((real(lunk + 1) - xlp) &
                                                *(real(munk + 1) - ylp)*(real(nunk + 1) - zlp))
                                            f(:lunk, :munk, :nunk) = &
                                                f(:lunk, :munk, :nunk) - pertrb
                                        end if
                                    end if
                                case (2:3, 5)
                                    exit exit_select_case_4
                            end select
                        case (2:3, 5)
                            exit exit_select_case_4
                    end select
                case (2:3, 5)
                    exit exit_select_case_4
            end select
        end do exit_select_case_4

        if (nbdcnd /= 0) then
            nperod = 1
            w(1) = 0.0_wp
            w(iww-1) = 0.0_wp
        else
            nperod = 0
        end if

        call pois3dd(lbdcnd, lunk, c1, mbdcnd, munk, c2, nperod, nunk, w, &
            w(iwb), w(iwc), ldimf, mdimf, f(lstart, mstart, nstart), ir, w(iww))
        !
        !==> Fill in sides for periodic boundary conditions.
        !
        if (lp == 1) then
            if (mp == 1) then
                f(1, mp1, nstart:nstop) = f(1, 1, nstart:nstop)
                mstop = mp1
            end if
            if (np == 1) then
                f(1, mstart:mstop, np1) = f(1, mstart:mstop, 1)
                nstop = np1
            end if
            f(lp1, mstart:mstop, nstart:nstop) = &
                f(1, mstart:mstop, nstart: nstop)
        end if

        if (mp == 1) then
            if (np == 1) then
                f(lstart:lstop, 1, np1) = f(lstart:lstop, 1, 1)
                nstop = np1
            end if
            f(lstart:lstop, mp1, nstart:nstop) = &
                f(lstart:lstop, 1, nstart:nstop)
        end if

        if (np == 1) then
            f(lstart:lstop, mstart:mstop, np1) = &
                f(lstart:lstop, mstart:mstop,1)
        end if

    end subroutine hw3crtt



    pure subroutine check_input_arguments(l, lbdcnd, m, mbdcnd, n, nbdcnd, &
        ldimf, mdimf, xs, xf, ys, yf, zs, zf, ierror)
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------
        integer (ip), intent (in)      :: l
        integer (ip), intent (in)      :: lbdcnd
        integer (ip), intent (in)      :: m
        integer (ip), intent (in)      :: mbdcnd
        integer (ip), intent (in)      :: n
        integer (ip), intent (in)      :: nbdcnd
        integer (ip), intent (in)      :: ldimf
        integer (ip), intent (in)      :: mdimf
        real (wp),    intent (in)      :: xs
        real (wp),    intent (in)      :: xf
        real (wp),    intent (in)      :: ys
        real (wp),    intent (in)      :: yf
        real (wp),    intent (in)      :: zs
        real (wp),    intent (in)      :: zf
        integer (ip), intent (out)     :: ierror
        !-----------------------------------------------
        ! Dictionary: calling arguments
        !-----------------------------------------------

        ! initialize error flag
        ierror = 0

        ! Case 1
        if (xf <= xs) then
            ierror = 1
            return
        end if

        if (l < 5) then
            ierror = 2
            return
        end if

        if (lbdcnd < 0 .or. lbdcnd > 4) then
            ierror = 3
            return
        end if

        if (yf <= ys) then
            ierror = 4
            return
        end if

        if (m < 5) then
            ierror = 5
            return
        end if

        if (mbdcnd < 0 .or. mbdcnd > 4) then
            ierror = 6
            return
        end if

        if (zf <= zs) then
            ierror = 7
            return
        end if

        if (n < 5) then
            ierror = 8
            return
        end if

        if (nbdcnd < 0 .or. nbdcnd > 4) then
            ierror = 9
            return
        end if

        if (ldimf < l + 1) then
            ierror = 10
            return
        end if

        if (mdimf < m + 1) then
            ierror = 11
            return
        end if

    end subroutine check_input_arguments



end module module_hw3crt
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
