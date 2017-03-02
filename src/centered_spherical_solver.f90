!
!     file hwsssp.f90
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
!     SUBROUTINE hwssp(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, bdps,
!                      bdpf, elmbda, f, idimf, pertrb, ierror)
!
! DIMENSION OF           bdts(n+1),    bdtf(n+1), bdps(m+1), bdpf(m+1),
! ARGUMENTS              f(idimf, n+1)
!
! LATEST REVISION        May 2016
!
! PURPOSE                Solves a finite difference approximation to
!                        the helmholtz equation in spherical
!                        coordinates and on the surface of the unit
!                        sphere (radius of 1).  the equation is
!
!                          (1/sin(theta))(d/dtheta)(sin(theta)
!                          (du/dtheta)) + (1/sin(theta)**2)(d/dphi)
!                          (du/dphi)  + lambda*u = f(theta, phi)
!
!                        where theta is colatitude and phi is
!                        longitude.
!
! USAGE                  call hwsssp(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf,
!                                     n, nbdcnd, bdps, bdpf, elmbda, f,
!                                     idimf, pertrb, ierror, w)
!
! ARGUMENTS
! ON INPUT               ts, tf
!
!                          the range of theta (colatitude), i.e.,
!                          ts <= theta <= tf. ts must be less
!                          than tf.  ts and tf are in radians.
!                          a ts of zero corresponds to the north
!                          pole and a tf of pi corresponds to
!                          the south pole.
!
!                          * * * important * * *
!
!                          if tf is equal to pi then it must be
!                          computed using the statement
!                          tf = pi_mach(dum). this insures that tf
!                          in the user's program is equal to pi in
!                          this program which permits several tests
!                          of the input parameters that otherwise
!                          would not be possible.
!
!
!                        m
!                          the number of panels into which the
!                          interval (ts, tf) is subdivided.
!                          hence, there will be m+1 grid points in the
!                          theta-direction given by
!                          theta(i) = (i-1)dtheta+ts for
!                          i = 1, 2, ..., m+1, where
!                          dtheta = (tf-ts)/m is the panel width.
!                          m must be greater than 5
!
!                        mbdcnd
!                          indicates the type of boundary condition
!                          at theta = ts and theta = tf.
!
!                          = 1  if the solution is specified at
!                               theta = ts and theta = tf.
!                          = 2  if the solution is specified at
!                               theta = ts and the derivative of
!                               the solution with respect to theta is
!                               specified at theta = tf
!                               (see note 2 below).
!                          = 3  if the derivative of the solution
!                               with respect to theta is specified
!                               specified at theta = ts and
!                               theta = tf (see notes 1, 2 below).
!                          = 4  if the derivative of the solution
!                               with respect to theta is specified
!                               at theta = ts (see note 1 below)
!                               and the solution is specified at
!                               theta = tf.
!                          = 5  if the solution is unspecified at
!                               theta = ts = 0 and the solution
!                               is specified at theta = tf.
!                          = 6  if the solution is unspecified at
!                               theta = ts = 0 and the derivative
!                               of the solution with respect to theta
!                               is specified at theta = tf
!                               (see note 2 below).
!                          = 7  if the solution is specified at
!                               theta = ts and the solution is
!                               is unspecified at theta = tf = pi.
!                          = 8  if the derivative of the solution
!                               with respect to theta is specified
!                               at theta = ts (see note 1 below) and
!                               the solution is unspecified at
!                               theta = tf = pi.
!                          = 9  if the solution is unspecified at
!                               theta = ts = 0 and theta = tf = pi.
!
!                          notes:
!                          if ts = 0, do not use mbdcnd = 3, 4, or 8,
!                          but instead use mbdcnd = 5, 6, or 9  .
!
!                          if tf = pi, do not use mbdcnd = 2, 3, or 6,
!                          but instead use mbdcnd = 7, 8, or 9  .
!
!                        bdts
!                          a one-dimensional array of length n+1 that
!                          specifies the values of the derivative of
!                          the solution with respect to theta at
!                          theta = ts.  when mbdcnd = 3, 4, or 8,
!
!                          bdts(j) = (d/dtheta)u(ts, phi(j)),
!                          j = 1, 2, ..., n+1  .
!
!                          when mbdcnd has any other value, bdts is
!                          a dummy variable.
!
!                        bdtf
!                          a one-dimensional array of length n+1
!                          that specifies the values of the derivative
!                          of the solution with respect to theta at
!                          theta = tf.  when mbdcnd = 2, 3, or 6,
!
!                          bdtf(j) = (d/dtheta)u(tf, phi(j)),
!                          j = 1, 2, ..., n+1  .
!
!                          when mbdcnd has any other value, bdtf is
!                          a dummy variable.
!
!                        ps, pf
!                          the range of phi (longitude), i.e.,
!                          ps <= phi <= pf.  ps must be less
!                          than pf.  ps and pf are in radians.
!                          if ps = 0 and pf = 2*pi, periodic
!                          boundary conditions are usually prescribed.
!
!                          * * * IMPORTANT * * *
!
!                          if pf is equal to 2*pi then it must be
!                          computed using the statement
!                          pf = TWO *pi_mach(dum). this insures that
!                          pf in the users program is equal to
!                          2*pi in this program which permits tests
!                          of the input parameters that otherwise
!                          would not be possible.
!
!                        n
!                          the number of panels into which the
!                          interval (ps, pf) is subdivided.
!                          hence, there will be n+1 grid points
!                          in the phi-direction given by
!                          phi(j) = (j-1)dphi+ps  for
!                          j = 1, 2, ..., n+1, where
!                          dphi = (pf-ps)/n is the panel width.
!                          n must be greater than 4
!
!                        nbdcnd
!                          indicates the type of boundary condition
!                          at phi = ps and phi = pf.
!
!                          = 0  if the solution is periodic in phi,
!                               i.u., u(i, j) = u(i, n+j).
!                          = 1  if the solution is specified at
!                               phi = ps and phi = pf
!                               (see note below).
!                          = 2  if the solution is specified at
!                               phi = ps (see note below)
!                               and the derivative of the solution
!                               with respect to phi is specified
!                               at phi = pf.
!                          = 3  if the derivative of the solution
!                               with respect to phi is specified
!                               at phi = ps and phi = pf.
!                          = 4  if the derivative of the solution
!                               with respect to phi is specified
!                               at ps and the solution is specified
!                               at phi = pf
!
!                          note:
!                          nbdcnd = 1, 2, or 4 cannot be used with
!                          mbdcnd = 5, 6, 7, 8, or 9.  the former indicates
!                          that the solution is specified at a pole, the
!                          latter indicates that the solution is not
!                          specified.  use instead  mbdcnd = 1 or 2.
!
!                        bdps
!                          a one-dimensional array of length m+1 that
!                          specifies the values of the derivative
!                          of the solution with respect to phi at
!                          phi = ps.  when nbdcnd = 3 or 4,
!
!                            bdps(i) = (d/dphi)u(theta(i), ps),
!                            i = 1, 2, ..., m+1  .
!
!                          when nbdcnd has any other value, bdps is
!                          a dummy variable.
!
!                        bdpf
!                          a one-dimensional array of length m+1 that
!                          specifies the values of the derivative
!                          of the solution with respect to phi at
!                          phi = pf.  when nbdcnd = 2 or 3,
!
!                            bdpf(i) = (d/dphi)u(theta(i), pf),
!                            i = 1, 2, ..., m+1  .
!
!                          when nbdcnd has any other value, bdpf is
!                          a dummy variable.
!
!                        elmbda
!                          the constant lambda in the helmholtz
!                          equation.  if lambda > 0, a solution
!                          may not exist.  however, hwsssp will
!                          attempt to find a solution.
!
!                        f
!                          a two-dimensional array that specifies the
!                          value of the right side of the helmholtz
!                          equation and boundary values (if any).
!                          f must be dimensioned at least (m+1)*(n+1).
!
!                          on the interior, f is defined as follows:
!                          for i = 2, 3, ..., m and j = 2, 3, ..., n
!                          f(i, j) = f(theta(i), phi(j)).
!
!                          on the boundaries f is defined as follows:
!                          for j = 1, 2, ..., n+1 and i = 1, 2, ..., m+1
!
!                          mbdcnd   f(1, j)            f(m+1, j)
!                          ------   ------------      ------------
!
!                            1      u(ts, phi(j))      u(tf, phi(j))
!                            2      u(ts, phi(j))      f(tf, phi(j))
!                            3      f(ts, phi(j))      f(tf, phi(j))
!                            4      f(ts, phi(j))      u(tf, phi(j))
!                            5      f(0, ps)           u(tf, phi(j))
!                            6      f(0, ps)           f(tf, phi(j))
!                            7      u(ts, phi(j))      f(pi, ps)
!                            8      f(ts, phi(j))      f(pi, ps)
!                            9      f(0, ps)           f(pi, ps)
!
!                          nbdcnd   f(i, 1)            f(i, n+1)
!                          ------   --------------    --------------
!
!                            0      f(theta(i), ps)    f(theta(i), ps)
!                            1      u(theta(i), ps)    u(theta(i), pf)
!                            2      u(theta(i), ps)    f(theta(i), pf)
!                            3      f(theta(i), ps)    f(theta(i), pf)
!                            4      f(theta(i), ps)    u(theta(i), pf)
!
!                          note:
!                          if the table calls for both the solution u
!                          and the right side f at a corner then the
!                          solution must be specified.
!
!                        idimf
!                          the row (or first) dimension of the array
!                          f as it appears in the program calling
!                          hwsssp.  this parameter is used to specify
!                          the variable dimension of f. idimf must be
!                          at least m+1  .
!
!
! ON OUTPUT              f
!                          contains the solution u(i, j) of the finite
!                          difference approximation for the grid point
!                          (theta(i), phi(j)),  i = 1, 2, ..., m+1  and
!                          j = 1, 2, ..., n+1  .
!
!                        pertrb
!                          if one specifies a combination of periodic,
!                          derivative or unspecified boundary
!                          conditions for a poisson equation
!                          (lambda = 0), a solution may not exist.
!                          pertrb is a constant, calculated and
!                          subtracted from f, which ensures that a
!                          solution exists.  hwsssp then computes
!                          this solution, which is a least squares
!                          solution to the original approximation.
!                          this solution is not unique and is
!                          unnormalized. the value of pertrb should
!                          be small compared to the right side f.
!                          otherwise , a solution is obtained to an
!                          essentially different problem. this
!                          comparison should always be made to insure
!                          that a meaningful solution has been
!                          obtained
!
!                        ierror
!                          an error flag that indicates invalid input
!                          parameters.  except for numbers 0 and 8,
!                          a solution is not attempted.
!
!                          = 0  no error
!                          = 1  ts<0 or tf>pi
!                          = 2  ts>=tf
!                          = 3  mbdcnd<1 or mbdcnd>9
!                          = 4  ps<0 or ps>pi+pi
!                          = 5  ps>=pf
!                          = 6  n<5
!                          = 7  m<5
!                          = 8  nbdcnd<0 or nbdcnd>4
!                          = 9  elmbda>0
!                          = 10 idimf<m+1
!                          = 11 nbdcnd equals 1, 2 or 4 and mbdcnd>=5
!                          = 12 ts==0 and mbdcnd equals 3, 4 or 8
!                          = 13 tf==pi and mbdcnd equals 2, 3 or 6
!                          = 14 mbdcnd equals 5, 6 or 9 and ts/=0
!                          = 15 mbdcnd>=7 and tf/=pi
!                          = 20 if the dynamic allocation of real and
!                               complex workspace required for solution
!                               fails (for example if n, m are too large
!                               for your computer)
!
!
! SPECIAL CONDITIONS     None
!
! I/O                    None
!
! PRECISION              64-bit double precision
!
! REQUIRED files         type_FishpackWorkspace.f90, genbun.f90, type_CyclicReductionUtility.f9090
!
! HISTORY                * Written by Roland Sweet at NCAR in the late
!                          1970's.  released on NCAR's public software
!                          libraries in January 1980.
!                        * Revised in June 2004 by John Adams using
!                          Fortran 90 dynamically allocated workspace.
!
! PORTABILITY            Fortran 2008
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
! REFERENCES             P. N. Swarztrauber, "The direct solution of
!                        the discrete poisson equation on the surface of
!                        a sphere", S.I.A.M. J. Numer. Anal., 15(1974), 
!                        pp 212-215.
!
!                        Swarztrauber, P. and R. Sweet, "Efficient
!                        FORTRAN subprograms for the solution of
!                        elliptic equations", NCAR TN/IA-109, July, 
!                        1975, 138 pp.
!
submodule(centered_helmholtz_solvers) centered_spherical_solver


! Parameters confined to the submodule

integer(ip), parameter :: IIWK = 7 ! Size for workspace_indices


contains

    module subroutine hwsssp(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
        bdps, bdpf, elmbda, f, idimf, pertrb, ierror)

        ! Dummy arguments

        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: mbdcnd
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: nbdcnd
        integer(ip), intent(in)     :: idimf
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: ts
        real(wp),    intent(in)     :: tf
        real(wp),    intent(in)     :: ps
        real(wp),    intent(in)     :: pf
        real(wp),    intent(in)     :: elmbda
        real(wp),    intent(out)    :: pertrb
        real(wp),    intent(in)     :: bdts(:)
        real(wp),    intent(in)     :: bdtf(:)
        real(wp),    intent(in)     :: bdps(:)
        real(wp),    intent(in)     :: bdpf(:)
        real(wp),    intent(inout)  :: f(:,:)

        ! Local variables
        type(FishpackWorkspace) :: workspace

        ! Check input arguments
        call hwsssp_check_input_arguments(ts, tf, m, mbdcnd, ps, pf, n, &
            nbdcnd, elmbda, idimf, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Allocate memory
        call hwsssp_initialize_workspace(n, m, workspace)

        ! Solve system
        associate( &
            indx => workspace%workspace_indices, &
            rew => workspace%real_workspace &
            )
            associate( &
                am => rew(indx(1):), &
                bm => rew(indx(2):), &
                cm => rew(indx(3):), &
                sn => rew(indx(4):), &
                ss => rew(indx(5):), &
                sint => rew(indx(6):), &
                d => rew(indx(7):) &
                )
                call hwsssp_lower_routine(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
                    bdps, bdpf, elmbda, f, idimf, pertrb, am, bm, cm, sn, ss, &
                    sint, d, ierror)
            end associate
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine hwsssp

    subroutine hwsssp_check_input_arguments(ts, tf, m, mbdcnd, ps, pf, n, &
        nbdcnd, elmbda, idimf, ierror)

        ! Dummy arguments
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: mbdcnd
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: nbdcnd
        integer(ip), intent(in)     :: idimf
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: elmbda
        real(wp),    intent(in)     :: ts
        real(wp),    intent(in)     :: tf
        real(wp),    intent(in)     :: ps
        real(wp),    intent(in)     :: pf

        if (ts < ZERO .or. tf > PI) then
            ierror = 1
        else if (ts >= tf) then
            ierror = 2
        else if (mbdcnd < 1 .or. mbdcnd > 9) then
            ierror = 3
        else if (ps < ZERO .or. pf > TWO_PI) then
            ierror = 4
        else if (ps >= pf) then
            ierror = 5
        else if (n < 5) then
            ierror = 6
        else if (m < 5) then
            ierror = 7
        else if (nbdcnd < 0 .or. nbdcnd > 4) then
            ierror = 8
        else if (elmbda > ZERO) then
            ierror = 9
        else if (idimf < m + 1) then
            ierror = 10
        else if ((nbdcnd == 1 .or. nbdcnd == 2 .or. nbdcnd == 4) .and. mbdcnd >= 5) then
            ierror = 11
        else if (ts == ZERO .and. (mbdcnd==3.or. mbdcnd==4 .or. mbdcnd == 8)) then
            ierror = 12
        else if (tf == PI .and. (mbdcnd==2.or. mbdcnd == 3 .or. mbdcnd == 6)) then
            ierror = 13
        else if ((mbdcnd == 5.or.mbdcnd == 6.or.mbdcnd == 9) .and. ts /= ZERO) then
            ierror = 14
        else if (mbdcnd >= 7 .and. tf /= PI) then
            ierror = 15
        else
            ierror = 0
        end if

    end subroutine hwsssp_check_input_arguments

    subroutine hwsssp_initialize_workspace(n, m, workspace)

        ! Dummy arguments
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: m
        class(FishpackWorkspace), intent(out) :: workspace

        ! Local variables
        integer(ip) :: irwk, icwk

        ! Compute workspace lengths for hwsssp
        irwk = 4 * (n + 1)  &
            + (16 + int(log(real(n+1,kind=wp))/log(TWO), kind=ip)) &
            * (m + 1)
        icwk = 0

        ! Allocate memory
        call workspace%create(irwk, icwk, IIWK)

        ! Set workspace indices
        workspace%workspace_indices = hwsssp_get_workspace_indices(m)

    end subroutine hwsssp_initialize_workspace

    pure function hwsssp_get_workspace_indices(m) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: m
        integer(ip)             :: return_value(IIWK)

        associate( indx => return_value)
            indx(1) = 1
            indx(2) = m+2
            indx(3) = 2*m+3
            indx(4) = 3*m+4
            indx(5) = 4*m+5
            indx(6) = 5*m+6
            indx(7) = 6*m+7
        end associate

    end function hwsssp_get_workspace_indices

    subroutine hwsssp_lower_routine(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
        bdps, bdpf, elmbda, f, idimf, pertrb, am, bm, cm, sn, ss, &
        sint, d, error_flag)

        ! Dummy arguments
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: mbdcnd
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: nbdcnd
        integer(ip), intent(in)     :: idimf
        real(wp),    intent(in)     :: ts
        real(wp),    intent(in)     :: tf
        real(wp),    intent(in)     :: ps
        real(wp),    intent(in)     :: pf
        real(wp),    intent(in)     :: elmbda
        real(wp),    intent(out)    :: pertrb
        real(wp),    intent(in)     :: bdts(:)
        real(wp),    intent(in)     :: bdtf(:)
        real(wp),    intent(in)     :: bdps(:)
        real(wp),    intent(in)     :: bdpf(:)
        real(wp),    intent(inout)  :: f(:,:)
        real(wp),    intent(out)    :: am(:)
        real(wp),    intent(out)    :: bm(:)
        real(wp),    intent(out)    :: cm(:)
        real(wp),    intent(out)    :: sn(:)
        real(wp),    intent(out)    :: ss(:)
        real(wp),    intent(out)    :: sint(:)
        real(wp),    intent(out)    :: d(:)
        integer(ip), intent(out)    :: error_flag

        ! Local variables
        integer(ip) :: mp1, np1, i, inp, isp, mbr, its, itf, itsp, itfm, munk
        integer(ip) :: iid, ii, nbr, jps, jpf, jpsp, jpfm, nunk, ising
        real(wp)    :: fn, fm, dth, half_dth, two_dth, dphi, two_dphi
        real(wp)    :: dphi2, edp2, dth2, cp, wpp, fim1, theta, t1, at, ct, wts, wtf
        real(wp)    :: wps, wpf, fjj, cf, summation, sum1, hne, yhld, sum2, dfn, dnn, dsn
        real(wp)    :: cnp, hld, dfs, dss, dns, csp, rtn, rts, den
        type(CenteredCyclicReductionUtility) :: util

        mp1 = m + 1
        np1 = n + 1
        fn = n
        fm = m
        dth = (tf - ts)/fm
        half_dth = dth/2
        two_dth = TWO * dth
        dphi = (pf - ps)/fn
        two_dphi = TWO * dphi
        dphi2 = dphi**2
        edp2 = elmbda*dphi2
        dth2 = dth**2
        cp = 4.0_wp/(fn*dth2)
        wpp = fn*sin(half_dth)/4

        ! Initialize constants
        dfn = ZERO
        dnn = ZERO
        dsn = ZERO

        do i = 1, mp1
            fim1 = i - 1
            theta = fim1*dth + ts
            sint(i) = sin(theta)

            if (sint(i) == ZERO) cycle

            t1 = ONE/(dth2*sint(i))
            am(i) = t1*sin(theta - half_dth)
            cm(i) = t1*sin(theta + half_dth)
            bm(i) = (-am(i)) - cm(i) + elmbda
        end do
        inp = 0
        isp = 0
        !
        ! boundary condition at theta=ts
        !
        mbr = mbdcnd + 1
        select case (mbr)
            case (1)
                its = 1
            case (2:3, 8)
                at = am(2)
                its = 2
            case (4:5, 9)
                at = am(1)
                its = 1
                cm(1) = am(1) + cm(1)
            case (6:7, 10)
                at = am(2)
                inp = 1
                its = 2
        end select

        !
        ! boundary condition theta=tf
        !
        select case (mbr)
            case (1)
                itf = m
            case (2, 5:6)
                ct = cm(m)
                itf = m
            case (3:4, 7)
                ct = cm(m+1)
                am(m+1) = am(m+1) + cm(m+1)
                itf = m + 1
            case (8:10)
                itf = m
                isp = 1
                ct = cm(m)
        end select

        !
        ! compute homogeneous solution with solution at pole equal to one
        !
        itsp = its + 1
        itfm = itf - 1
        wts = sint(its+1)*am(its+1)/cm(its)
        wtf = sint(itf-1)*cm(itf-1)/am(itf)
        munk = itf - its + 1

        if (isp > 0) then
            d(its) = cm(its)/bm(its)
            do i = itsp, m
                d(i) = cm(i)/(bm(i)-am(i)*d(i-1))
            end do
            ss(m) = -d(m)
            iid = m - its
            do ii = 1, iid
                i = m - ii
                ss(i) = -d(i)*ss(i+1)
            end do
            ss(m+1) = ONE
        end if

        if (inp > 0) then
            sn(1) = ONE
            d(itf) = am(itf)/bm(itf)
            iid = itf - 2
            do ii = 1, iid
                i = itf - ii
                d(i) = am(i)/(bm(i)-cm(i)*d(i+1))
            end do
            sn(2) = -d(2)
            do i = 3, itf
                sn(i) = -d(i)*sn(i-1)
            end do
        end if
        !
        ! boundary conditions at phi=ps
        !
        nbr = nbdcnd + 1
        wps = ONE
        wpf = ONE
        select case (nbr)
            case default
                jps = 1
            case (2:3)
                jps = 2
            case (4:5)
                jps = 1
                wps = HALF
        end select
        !
        ! boundary condition at phi=pf
        !
        select case (nbr)
            case (1)
                jpf = n
            case (2, 5)
                jpf = n
            case (3:4)
                wpf = HALF
                jpf = n + 1
        end select

        jpsp = jps + 1
        jpfm = jpf - 1
        nunk = jpf - jps + 1
        fjj = jpfm - jpsp + 1
        !
        ! scale coefficients for subroutine genbun
        !
        do i = its, itf
            cf = dphi2*sint(i)*sint(i)
            am(i) = cf*am(i)
            bm(i) = cf*bm(i)
            cm(i) = cf*cm(i)
        end do

        am(its) = ZERO
        cm(itf) = ZERO
        ising = 0

        select case (mbr)
            case (1, 4, 7, 9:10)
                select case (nbr)
                    case (1, 4)
                        if (elmbda >= ZERO) then
                            ising = 1
                            summation = wts*wps + wts*wpf + wtf*wps + wtf*wpf

                            if (inp > 0) summation = summation + wpp

                            if (isp > 0) summation = summation + wpp

                            sum1 = sum(sint(itsp:itfm))
                            summation = summation + fjj*(sum1 + wts + wtf)
                            summation = summation + (wps + wpf)*sum1
                            hne = summation
                        end if
                end select
        end select

        select case (mbr)
            case (2:3, 8)
                f(2,jps:jpf) = f(2,jps:jpf) - at*f(1,jps:jpf)
            case (4:5, 9)
                f(1,jps:jpf) = f(1,jps:jpf) + two_dth*bdts(jps:jpf)*at
            case (6:7, 10)
                if (nbdcnd == 3) then
                    yhld = f(1,jps) - 4.0_wp/(fn*dphi*dth2)*(bdpf(2)-bdps(2))
                    f(1,:np1) = yhld
                end if
        end select

        select case (mbr)
            case (2, 5:6)
                f(m,jps:jpf) = f(m,jps:jpf) - ct*f(m+1,jps:jpf)
            case (3:4, 7)
                f(m+1,jps:jpf) = f(m+1,jps:jpf) - two_dth*bdtf(jps:jpf)*ct
            case (8:10)
                if (nbdcnd == 3) then
                    yhld = f(m+1,jps) - 4.0_wp/(fn*dphi*dth2)*(bdpf(m)-bdps(m))
                    f(m+1,:np1) = yhld
                end if
        end select

        select case (nbr)
            case (2:3)
                f(its:itf,2) = f(its:itf,2) &
                    - f(its:itf,1)/(dphi2*sint(its:itf)*sint(its:itf))
            case (4:5)
                f(its:itf,1) = f(its:itf,1) &
                    + two_dphi*bdps(its:itf)/(dphi2*sint(its:itf)*sint(its:itf))
        end select

        select case (nbr)
            case (2, 5)
                f(its:itf,n) = f(its:itf,n) &
                    - f(its:itf,n+1)/(dphi2*sint(its:itf)*sint(its:itf))
            case (3:4)
                f(its:itf,n+1) = f(its:itf,n+1) &
                    - two_dphi*bdpf(its:itf)/(dphi2*sint(its:itf)*sint(its:itf))
        end select

        pertrb = ZERO

        if (ising /= 0) then
            summation = &
                wts*wps*f(its,jps) + wts*wpf*f(its,jpf) &
                + wtf*wps*f(itf,jps) + wtf*wpf*f(itf,jpf)

            if (inp > 0) summation = summation + wpp*f(1,jps)
            if (isp > 0) summation = summation + wpp*f(m+1,jps)

            do i = itsp, itfm
                summation = summation + sint(i)*sum(f(i,jpsp:jpfm))
            end do

            sum1 = sum(f(its,jpsp:jpfm))
            sum2 = sum(f(itf,jpsp:jpfm))
            summation = summation + wts*sum1 + wtf*sum2
            sum1 = ZERO
            sum2 = ZERO
            sum1 = dot_product(sint(itsp:itfm),f(itsp:itfm,jps))
            sum2 = dot_product(sint(itsp:itfm),f(itsp:itfm,jpf))
            summation = summation + wps*sum1 + wpf*sum2
            pertrb = summation/hne
            f(:mp1,:np1) = f(:mp1,:np1) - pertrb
        end if

        ! scale right side for subroutine genbun_lower_routine
        do i = its, itf
            cf = dphi2*sint(i)*sint(i)
            f(i,jps:jpf) = cf*f(i,jps:jpf)
        end do

        call_arguments: associate( &
            a_arg => am(its:), &
            b_arg => bm(its:), &
            c_arg => cm(its:), &
            y_arg => f(its:,jps:) &
            )

            ! Invoke genbun_lower_routine solver
            call util%genbun_lower_routine(nbdcnd, nunk, 1, munk, a_arg, b_arg, c_arg, &
                idimf, y_arg, error_flag, d)
        end associate call_arguments

        ! Check error flag
        if (error_flag /= 0) then
            error stop 'fishpack library: genbun_lower_routine call failed in hwsssp_lower_routine'
        end if

        if (ising > 0 .and. inp > 0 .and. isp <= 0) then
            f(1,:np1) = ZERO
        else if (isp > 0) then
            f(m+1,:np1) = ZERO
        else if (inp > 0) then

            summation = wps*f(its,jps) + wpf*f(its,jpf) + sum(f(its,jpsp:jpfm))
            dfn = cp*summation
            dnn = cp*((wps + wpf + fjj)*(sn(2)-ONE)) + elmbda
            dsn = cp*(wps + wpf + fjj)*sn(m)

            if (isp > 0) then
                cnp = (f(1,1)-dfn)/dnn
                do i = its, itf
                    hld = cnp*sn(i)
                    f(i,jps:jpf) = f(i,jps:jpf) + hld
                end do
                f(1,:np1) = cnp
            end if
        else
            summation = wps*f(itf,jps) + wpf*f(itf,jpf) + sum(f(itf,jpsp:jpfm))
            dfs = cp*summation
            dss = cp*((wps + wpf + fjj)*(ss(m)-ONE)) + elmbda
            dns = cp*(wps + wpf + fjj)*ss(2)

            if (inp <= 0) then
                csp = (f(m+1,1)-dfs)/dss
                do i = its, itf
                    hld = csp*ss(i)
                    f(i,jps:jpf) = f(i,jps:jpf) + hld
                end do
                f(m+1,:np1) = csp
            else
                rtn = f(1,1) - dfn
                rts = f(m+1,1) - dfs
                if (ising > 0) then
                    csp = ZERO
                    cnp = rtn/dnn
                else
                    if (abs(dnn) - abs(dsn) > ZERO) then
                        den = dss - dns*dsn/dnn
                        rts = rts - rtn*dsn/dnn
                        csp = rts/den
                        cnp = (rtn - csp*dns)/dnn
                    else
                        den = dns - dss*dnn/dsn
                        rtn = rtn - rts*dnn/dsn
                        csp = rtn/den
                        cnp = (rts - dss*csp)/dsn
                    end if
                end if
                do i = its, itf
                    hld = cnp*sn(i) + csp*ss(i)
                    f(i,jps:jpf) = f(i,jps:jpf) + hld
                end do
                f(1,:np1) = cnp
                f(m+1,:np1) = csp
            end if
        end if

        if (nbdcnd == 0) f(:mp1,jpf+1) = f(:mp1,jps)

    end subroutine hwsssp_lower_routine

end submodule centered_spherical_solver
