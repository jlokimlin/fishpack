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

submodule(centered_helmholtz_solvers) centered_axisymmetric_spherical_solver

! Parameters confined to the submodule
integer(ip), parameter :: SIZE_OF_WORKSPACE_INDICES = 10

contains

    ! SUBROUTINE hwscsp(intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, nbdcnd, &
    !                       bdrs, bdrf, elmbda, f, idimf, pertrb, ierror, w)
    !
    !
    ! DIMENSION OF           bdts(n+1),  bdtf(n+1), bdrs(m+1), bdrf(m+1),
    ! ARGUMENTS              f(idimf, n+1)
    !
    !
    ! PURPOSE                Solves a finite difference approximation
    !                        to the modified helmholtz equation in
    !                        spherical coordinates assuming axisymmetry
    !                        (no dependence on longitude).  the equation
    !                        is
    !
    !                          (1/r**2)(d/dr)((r**2)(d/dr)u) +
    !
    !                          (1/(r**2)sin(theta))(d/dtheta)
    !
    !                          (sin(theta)(d/dtheta)u) +
    !
    !                          (lambda/(rsin(theta))**2)u = f(theta, r).
    !
    !                        this two dimensional modified helmholtz
    !                        equation results from the fourier transform
    !                        of the three dimensional poisson equation.
    !
    ! USAGE                  call hwscsp(intl, ts, tf, m, mbdcnd, bdts, bdtf,
    !                                    rs, rf, n, nbdcnd, bdrs, bdrf, elmbda,
    !                                    f, idimf, pertrb, ierror, w)
    !
    ! ARGUMENTS
    ! ON INPUT               intl
    !                          = 0  on initial entry to hwscsp or if any
    !                               of the arguments rs, rf, n, nbdcnd
    !                               are changed from a previous call.
    !                          = 1  if rs, rf, n, nbdcnd are all unchanged
    !                               from previous call to hwscsp.
    !
    !                          note:
    !                          a call with intl=0 takes approximately
    !                          1.5 times as much time as a call with
    !                          intl = 1  .  once a call with intl = 0
    !                          has been made then subsequent solutions
    !                          corresponding to different f, bdts, bdtf,
    !                          bdrs, bdrf can be obtained faster with
    !                          intl = 1 since initialization is not
    !                          repeated.
    !
    !                        ts, tf
    !                          the range of theta (colatitude), i.e.,
    !                          ts <= theta <= tf. ts must be less
    !                          than tf.  ts and tf are in radians. a ts of
    !                          zero corresponds to the north pole and a
    !                          tf of pi corresponds to the south pole.
    !
    !                          **** important ****
    !
    !                          if tf is equal to pi then it must be
    !                          computed using the statement
    !                          tf = pi_mach(dum). this insures that tf
    !                          in the user's program is equal to pi in
    !                          this program which permits several tests
    !                          of the  input parameters that otherwise
    !                          would not be possible.
    !
    !                        m
    !                          the number of panels into which the
    !                          interval (ts, tf) is subdivided.
    !                          hence, there will be m+1 grid points
    !                          in the theta-direction given by
    !                          theta(k) = (i-1)dtheta+ts for
    !                          i = 1, 2, ..., m+1, where dtheta = (tf-ts)/m
    !                          is the panel width.
    !
    !                        mbdcnd
    !                          indicates the type of boundary condition
    !                          at theta = ts and  theta = tf.
    !
    !                          = 1  if the solution is specified at
    !                               theta = ts and theta = tf.
    !                          = 2  if the solution is specified at
    !                               theta = ts and the derivative of the
    !                               solution with respect to theta is
    !                               specified at theta = tf
    !                               (see note 2 below).
    !                          = 3  if the derivative of the solution
    !                               with respect to theta is specified
    !                               at theta = ts and theta = tf
    !                               (see notes 1, 2 below).
    !                          = 4  if the derivative of the solution
    !                               with respect to theta is specified
    !                               at theta = ts (see note 1 below) and
    !                               solution is specified at theta = tf.
    !                          = 5  if the solution is unspecified at
    !                               theta = ts = 0 and the solution is
    !                                specified at theta = tf.
    !                          = 6  if the solution is unspecified at
    !                               theta = ts = 0 and the derivative
    !                               of the solution with respect to theta
    !                               is specified at theta = tf
    !                               (see note 2 below).
    !                          = 7  if the solution is specified at
    !                               theta = ts and the solution is
    !                                unspecified at theta = tf = pi.
    !                          = 8  if the derivative of the solution
    !                               with respect to theta is specified
    !                               at theta = ts (see note 1 below)
    !                               and the solution is unspecified at
    !                               theta = tf = pi.
    !                          = 9  if the solution is unspecified at
    !                               theta = ts = 0 and theta = tf = pi.
    !
    !                          note 1:
    !                          if ts = 0, do not use mbdcnd = 3, 4, or 8,
    !                          but instead use mbdcnd = 5, 6, or 9  .
    !
    !                          note 2:
    !                          if tf = pi, do not use mbdcnd = 2, 3, or 6,
    !                          but instead use mbdcnd = 7, 8, or 9  .
    !
    !                        bdts
    !                          a one-dimensional array of length n+1 that
    !                          specifies the values of the derivative of
    !                          the solution with respect to theta at
    !                          theta = ts.  when mbdcnd = 3, 4, or 8,
    !
    !                            bdts(j) = (d/dtheta)u(ts, r(j)),
    !                            j = 1, 2, ..., n+1  .
    !
    !                          when mbdcnd has any other value, bdts is
    !                          a dummy variable.
    !
    !                        bdtf
    !                          a one-dimensional array of length n+1 that
    !                          specifies the values of the derivative of
    !                          the solution with respect to theta at
    !                          theta = tf.  when mbdcnd = 2, 3, or 6,
    !
    !                          bdtf(j) = (d/dtheta)u(tf, r(j)),
    !                          j = 1, 2, ..., n+1  .
    !
    !                          when mbdcnd has any other value, bdtf is
    !                          a dummy variable.
    !
    !                        rs, rf
    !                          the range of r, i.e., rs <= r < rf.
    !                          rs must be less than rf.  rs must be
    !                          non-negative.
    !
    !                        n
    !                          the number of panels into which the
    !                          interval (rs, rf) is subdivided.
    !                          hence, there will be n+1 grid points in the
    !                          r-direction given by r(j) = (j-1)dr+rs
    !                          for j = 1, 2, ..., n+1, where dr = (rf-rs)/n
    !                          is the panel width.
    !                          n must be greater than 2
    !
    !                        nbdcnd
    !                          indicates the type of boundary condition
    !                          at r = rs and r = rf.
    !
    !                          = 1  if the solution is specified at
    !                               r = rs and r = rf.
    !                          = 2  if the solution is specified at
    !                               r = rs and the derivative
    !                               of the solution with respect to r
    !                               is specified at r = rf.
    !                          = 3  if the derivative of the solution
    !                               with respect to r is specified at
    !                               r = rs and r = rf.
    !                          = 4  if the derivative of the solution
    !                               with respect to r is specified at
    !                               rs and the solution is specified at
    !                               r = rf.
    !                          = 5  if the solution is unspecified at
    !                               r = rs = 0 (see note below)  and the
    !                               solution is specified at r = rf.
    !                          = 6  if the solution is unspecified at
    !                               r = rs = 0 (see note below) and the
    !                               derivative of the solution with
    !                               respect to r is specified at r = rf.
    !
    !                          note:
    !                          nbdcnd = 5 or 6 cannot be used with
    !                          mbdcnd = 1, 2, 4, 5, or 7.  the former
    !                          indicates that the solution is unspecified
    !                          at r = 0, the latter indicates that the
    !                          solution is specified).
    !                          use instead   nbdcnd = 1 or 2  .
    !
    !                        bdrs
    !                          a one-dimensional array of length m+1 that
    !                          specifies the values of the derivative of
    !                          the solution with respect to r at r = rs.
    !
    !                          when nbdcnd = 3 or 4,
    !                            bdrs(i) = (d/dr)u(theta(i), rs),
    !                            i = 1, 2, ..., m+1  .
    !
    !                          when nbdcnd has any other value, bdrs is
    !                          a dummy variable.
    !
    !                        bdrf
    !                          a one-dimensional array of length m+1
    !                          that specifies the values of the
    !                          derivative of the solution with respect
    !                          to r at r = rf.
    !
    !                          when nbdcnd = 2, 3, or 6,
    !                            bdrf(i) = (d/dr)u(theta(i), rf),
    !                            i = 1, 2, ..., m+1  .
    !
    !                          when nbdcnd has any other value, bdrf is
    !                          a dummy variable.
    !
    !                        elmbda
    !                          the constant lambda in the helmholtz
    !                          equation.  if lambda > 0, a solution
    !                          may not exist.  however, hwscsp will
    !                          attempt to find a solution.  if nbdcnd = 5
    !                          or 6 or  mbdcnd = 5, 6, 7, 8, or 9, elmbda
    !                          must be zero.
    !
    !                        f
    !                          a two-dimensional array, of dimension at
    !                          least (m+1)*(n+1), specifying values of the
    !                          right side of the helmholtz equation and
    !                          boundary values (if any).
    !
    !                          on the interior, f is defined as follows:
    !                          for i = 2, 3, ..., m and j = 2, 3, ..., n
    !                          f(i, j) = f(theta(i), r(j)).
    !
    !                          on the boundaries, f is defined as follows:
    !                          for j=1, 2, ..., n+1,  i=1, 2, ..., m+1,
    !
    !                          mbdcnd   f(1, j)            f(m+1, j)
    !                          ------   ----------        ----------
    !
    !                            1      u(ts, r(j))        u(tf, r(j))
    !                            2      u(ts, r(j))        f(tf, r(j))
    !                            3      f(ts, r(j))        f(tf, r(j))
    !                            4      f(ts, r(j))        u(tf, r(j))
    !                            5      f(0, r(j))         u(tf, r(j))
    !                            6      f(0, r(j))         f(tf, r(j))
    !                            7      u(ts, r(j))        f(pi, r(j))
    !                            8      f(ts, r(j))        f(pi, r(j))
    !                            9      f(0, r(j))         f(pi, r(j))
    !
    !                            nbdcnd   f(i, 1)            f(i, n+1)
    !                            ------   --------------    --------------
    !
    !                              1      u(theta(i), rs)    u(theta(i), rf)
    !                              2      u(theta(i), rs)    f(theta(i), rf)
    !                              3      f(theta(i), rs)    f(theta(i), rf)
    !                              4      f(theta(i), rs)    u(theta(i), rf)
    !                              5      f(ts, 0)           u(theta(i), rf)
    !                              6      f(ts, 0)           f(theta(i), rf)
    !
    !                          note:
    !                          if the table calls for both the solution
    !                          u and the right side f at a corner then
    !                          the solution must be specified.
    !
    !                        idimf
    !                          the row (or first) dimension of the array
    !                          f as it appears in the program calling
    !                          hwscsp.  this parameter is used to specify
    !                          the variable dimension of f.  idimf must
    !                          be at least m+1  .
    !
    !                        w
    !                          A FishpackWorkspace derived data type variable
    !                          that must be declared by the user.  the first
    !                          two declarative statements in the user program
    !                          calling sepeli must be:
    !
    !                               use type_fishpackworkspace
    !                               type(fishpackworkspace) :: w
    !
    !                          The first statement makes the fishpack module
    !                          defined in the file "type_FishpackWorkspace.f90" available to the
    !                          user program calling hwscsp. The second statement
    !                          declares a derived type variable (defined in
    !                          the module "type_FishpackWorkspace.f90") which is used internally
    !                          in hwscsp to dynamically allocate real and complex
    !                          workspace used in solution.  an error flag
    !                          (ierror = 20) is set if the required workspace
    !                          allocation fails (for example if n, m are too large)
    !                          real and complex values are set in the components
    !                          of w on a initial (intl=0) call to hwscsp.  these
    !                          must be preserved on non-initial calls (intl=1)
    !                          to hwscsp.  this eliminates redundant calculations
    !                          and saves compute time.
    !               ****       IMPORTANT!  The user program calling hwscsp should
    !                          include the statement:
    !
    !                               call w%destroy()
    !
    !                          after the final approximation is generated by
    !                          hwscsp.  The will deallocate the real and complex
    !                          workspace of W.  Failure to include this statement
    !                          could result in serious memory leakage.
    !
    !
    ! ON OUTPUT              f
    !                          contains the solution u(i, j) of the finite
    !                          difference approximation for the grid point
    !                          (theta(i), r(j)),  i = 1, 2, ..., m+1,
    !                                            j = 1, 2, ..., n+1  .
    !
    !                        pertrb
    !                          if a combination of periodic or derivative
    !                          boundary conditions is specified for a
    !                          poisson equation (lambda = 0), a solution
    !                          may not exist.  pertrb is a constant,
    !                          calculated and subtracted from f, which
    !                          ensures that a solution exists.  hwscsp
    !                          then computes this solution, which is a
    !                          least squares solution to the original
    !                          approximation. this solution is not unique
    !                          and is unnormalized. the value of pertrb
    !                          should be small compared to the right side
    !                          f. otherwise , a solution is obtained to
    !                          an essentially different problem. this
    !                          comparison should always be made to insure
    !                          that a meaningful solution has been obtained.
    !
    !                        ierror
    !                          an error flag that indicates invalid input
    !                          parameters.  except for numbers 0 and 10,
    !                          a solution is not attempted.
    !
    !                          = 1  ts<0. or tf>pi
    !                          = 2  ts>=tf
    !                          = 3  m<5
    !                          = 4  mbdcnd<1 or mbdcnd>9
    !                          = 5  rs<0
    !                          = 6  rs>=rf
    !                          = 7  n<5
    !                          = 8  nbdcnd<1 or nbdcnd>6
    !                          = 9  elmbda>0
    !                          = 10 idimf<m+1
    !                          = 11 elmbda/=0 and mbdcnd>=5
    !                          = 12 elmbda/=0 and nbdcnd equals 5 or 6
    !                          = 13 mbdcnd equals 5, 6 or 9 and ts/=0
    !                          = 14 mbdcnd>=7 and tf/=pi
    !                          = 15 ts.eq.0 and mbdcnd equals 3, 4 or 8
    !                          = 16 tf.eq.pi and mbdcnd equals 2, 3 or 6
    !                          = 17 nbdcnd>=5 and rs/=0
    !                          = 18 nbdcnd>=5 and mbdcnd equals 1, 2, 4, 5 or
    !                          = 20 if the dynamic allocation of real and
    !                               complex workspace in the derived type
    !                               (fishpackworkspace) variable w fails (e.g.,
    !                               if n, m are too large for the platform used)
    !
    !                          since this is the only means of indicating
    !                          a possliby incorrect call to hwscsp, the
    !                          user should test ierror after a call.
    !
    !                        w
    !                          the derived type(fishpackworkspace) variable w
    !                          contains real and complex values that must not
    !                          be destroyed if hwscsp is called again with
    !                          intl=1.
    !
    !
    ! HISTORY                Written by Roland Sweet AT NCAR in the late
    !                        1970'S.  Released on NCAR's public software
    !                        libraries in January 1980. Revised by John
    !                        Adams in June 2004 using Fortran 90 dynamically
    !                        allocated workspace and derived datat types
    !                        to eliminate mixed mode conflicts in the earlier
    !                        versions.
    !
    !
    ! ALGORITHM              The routine defines the finite difference
    !                        equations, incorporates boundary data, and
    !                        adjusts the right side of singular systems
    !                        and then calls blktri to solve the system.
    !
    ! REFERENCES             Swarztrauber, P. and R. Sweet, "Efficient
    !                        FORTRAN subprograms for the solution of
    !                        elliptic equations"
    !                          NCAR TN/IA-109, July, 1975, 138 pp.
    !
    module subroutine hwscsp(intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, &
        nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, ierror, workspace)

        ! Dummy arguments
        integer(ip), intent(in)     :: intl
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: mbdcnd
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: nbdcnd
        integer(ip), intent(in)     :: idimf
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: ts
        real(wp),    intent(in)     :: tf
        real(wp),    intent(in)     :: rs
        real(wp),    intent(in)     :: rf
        real(wp),    intent(in)     :: elmbda
        real(wp),    intent(out)    :: pertrb
        real(wp),    intent(in)     :: bdts(:)
        real(wp),    intent(in)     :: bdtf(:)
        real(wp),    intent(in)     :: bdrs(:)
        real(wp),    intent(in)     :: bdrf(:)
        real(wp),    intent(inout)  :: f(:,:)
        class(FishpackWorkspace), intent(inout)  :: workspace

        ! Check input arguments
        call hwscsp_check_input_arguments(ts, tf, m, mbdcnd, rs, rf, &
            n, nbdcnd, elmbda, idimf, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Set up workspace on initial call only
        if (intl == 0) call hwscsp_initialize_workspace(n, m, nbdcnd, workspace)

        ! Solve system
        associate( &
            indx => workspace%workspace_indices, &
            rew => workspace%real_workspace, &
            cxw => workspace%complex_workspace &
            )
            associate( &
                w => rew, &
                wc => cxw, &
                s => rew(indx(1):), &
                an => rew(indx(2):), &
                bn => rew(indx(3):), &
                cn => rew(indx(4):), &
                r => rew(indx(5):), &
                am => rew(indx(6):), &
                bm => rew(indx(7):), &
                cm => rew(indx(8):), &
                sint => rew(indx(9):), &
                bmh => rew(indx(10):) &
                )
                call hwscsp_lower_routine(intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, &
                    nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, w, wc, s, an, bn, &
                    cn, r, am, bm, cm, sint, bmh, ierror)
            end associate
        end associate

    end subroutine hwscsp

    subroutine hwscsp_initialize_workspace(n, m, nbdcnd, workspace)

        ! Dummy arguments
        integer(ip), intent(in)  :: nbdcnd
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: m
        class(FishpackWorkspace), intent(out) :: workspace

        ! Local variables
        integer(ip) :: irwk, icwk, indx(SIZE_OF_WORKSPACE_INDICES)

        ! Compute workspace indices
        indx = hwscsp_get_workspace_indices(n, m, nbdcnd)

        ! Compute required blktri workspace lengths
        call workspace%compute_blktri_workspace_lengths(n, m, irwk, icwk)

        ! Adjust workspace requirements for hwscsp
        irwk = indx(10) + m + 1
        icwk = icwk + 3 * m

        ! Allocate memory
        call workspace%create(irwk, icwk, SIZE_OF_WORKSPACE_INDICES)

        ! Set workspace indices
        workspace%workspace_indices = indx

    end subroutine hwscsp_initialize_workspace

    pure function hwscsp_get_workspace_indices(n, m, nbdcnd) result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: n, m, nbdcnd
        integer(ip)             :: return_value(SIZE_OF_WORKSPACE_INDICES)

        ! Local variables
        integer(ip)  :: j, nck, l, k

        nck = n

        select case (nbdcnd)
            case (1, 5)
                nck = nck - 1
            case (3)
                nck = nck + 1
        end select

        k = 1
        l = 4
        k = k + 1

        do
            if (nck <= l) exit
            l = 2*l
            k = k + 1
        end do

        l = 2*l

        associate( indx => return_value)

            indx(1) = (k - 2)*l + k + max(2*n, 6*m) + 13

            do j = 1, 5
                indx(j+1) = indx(j) + n + 1
            end do

            do j = 6, 9
                indx(j+1) = indx(j) + m + 1
            end do

        end associate

    end function hwscsp_get_workspace_indices

    subroutine hwscsp_check_input_arguments(ts, tf, m, mbdcnd, rs, rf, &
        n, nbdcnd, elmbda, idimf, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: m
        integer(ip), intent(in)  :: mbdcnd
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: nbdcnd
        integer(ip), intent(in)  :: idimf
        integer(ip), intent(out) :: ierror
        real(wp),    intent(in)  :: ts
        real(wp),    intent(in)  :: tf
        real(wp),    intent(in)  :: rs
        real(wp),    intent(in)  :: rf
        real(wp),    intent(in)  :: elmbda

        if (ts < ZERO .or. tf > PI) then
            ierror = 1
            return
        else if (ts >= tf) then
            ierror = 2
            return
        else if (m < 5) then
            ierror = 3
            return
        else if (mbdcnd < 1 .or. mbdcnd > 9) then
            ierror = 4
            return
        else if (rs < ZERO) then
            ierror = 5
            return
        else if (rs >= rf) then
            ierror = 6
            return
        else if (n < 5) then
            ierror = 7
            return
        else if (nbdcnd < 1 .or. nbdcnd > 6) then
            ierror = 8
            return
        else if (elmbda > ZERO) then
            ierror = 9
            return
        else if (idimf < m + 1) then
            ierror = 10
            return
        else if (elmbda /= ZERO .and. mbdcnd >= 5) then
            ierror = 11
            return
        else if (elmbda /= ZERO) then
            select case (nbdcnd)
                case (5, 6)
                    ierror = 12
                    return
            end select
        else if (ts /= ZERO) then
            select case (mbdcnd)
                case (5, 6, 9)
                    ierror = 13
                    return
            end select
        else if (mbdcnd >= 7 .and. tf /= PI) then
            ierror = 14
            return
        else if (ts == ZERO) then
            select case (mbdcnd)
                case (3:4, 8)
                    ierror = 15
                    return
            end select
        else if (tf == PI) then
            select case (mbdcnd)
                case (2, 3, 6)
                    ierror = 16
                    return
            end select
        else if (nbdcnd >= 5 .and. rs /= ZERO) then
            ierror = 17
            return
        else if (5 <= nbdcnd) then
            select case (mbdcnd)
                case (1:2, 5, 7)
                    ierror = 18
                    return
            end select
        else
            ierror = 0
        end if

    end subroutine hwscsp_check_input_arguments

    subroutine hwscsp_lower_routine(intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, &
        nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, w, wc, s, an, bn &
        , cn, r, am, bm, cm, sint, bmh, ierror)

        ! Dummy arguments

        integer(ip), intent(in)     :: intl
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: mbdcnd
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: nbdcnd
        integer(ip), intent(in)     :: idimf
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: ts
        real(wp),    intent(in)     :: tf
        real(wp),    intent(in)     :: rs
        real(wp),    intent(in)     :: rf
        real(wp),    intent(in)     :: elmbda
        real(wp),    intent(out)    :: pertrb
        real(wp),    intent(in)     :: bdts(:)
        real(wp),    intent(in)     :: bdtf(:)
        real(wp),    intent(in)     :: bdrs(:)
        real(wp),    intent(in)     :: bdrf(:)
        real(wp),    intent(inout) :: f(idimf,n + 1)
        real(wp),    intent(out) :: w(:)
        real(wp),    intent(out) :: s(:)
        real(wp),    intent(out) :: an(:)
        real(wp),    intent(out) :: bn(:)
        real(wp),    intent(out) :: cn(:)
        real(wp),    intent(out) :: r(:)
        real(wp),    intent(out) :: am(:)
        real(wp),    intent(out) :: bm(:)
        real(wp),    intent(out) :: cm(:)
        real(wp),    intent(out) :: sint(:)
        real(wp),    intent(out) :: bmh(:)
        complex(wp), intent(out) :: wc(:)

        ! Local variables
        integer(ip)         :: mp1, i, np1, j, mp, np
        integer(ip)         :: its, itf, itsp, itfm, ictr, jrs
        integer(ip)         :: l, jrf, jrsp, jrfm, munk, nunk, ising, iflg
        real(wp)            :: dth, tdt, hdth, sdts
        real(wp)            :: theta, t1, dr, hdr
        real(wp)            :: tdr, dr2, czr, at, ct, wts, wtf
        real(wp)            :: ar, wtnm, yps, cr, wrs, wrf
        real(wp)            :: wrz, summation, r2, hne, yhld
        real(wp)            :: rs2, rf2, rsq, xp, yph, xps
        real(wp), parameter :: FOUR = 4.0_wp
        real(wp), parameter :: SIX = 6.0_wp
        type(GeneralizedCyclicReductionUtility) :: util


        mp1 = m + 1
        dth = (tf - ts)/m
        tdt = dth + dth
        hdth = dth/2
        sdts = ONE/(dth**2)

        do i = 1, mp1
            theta = ts + real(i - 1, kind=wp)*dth
            sint(i) = sin(theta)
            if (sint(i) == ZERO) cycle
            t1 = sdts/sint(i)
            am(i) = t1*sin(theta - hdth)
            cm(i) = t1*sin(theta + hdth)
            bm(i) = -(am(i)+cm(i))
        end do

        np1 = n + 1
        dr = (rf - rs)/n
        hdr = dr/2
        tdr = dr + dr
        dr2 = dr**2
        czr = SIX*dth/(dr2*(cos(ts) - cos(tf)))

        do j = 1, np1
            r(j) = rs + real(j - 1, kind=wp)*dr
            an(j) = (r(j)-hdr)**2/dr2
            cn(j) = (r(j)+hdr)**2/dr2
            bn(j) = -(an(j)+cn(j))
        end do

        mp = 1
        np = 1

        !
        ! Boundary condition at phi=ps
        !
        select case (mbdcnd)
            case (1:2,7)
                at = am(2)
                its = 2
            case (3:4,8)
                at = am(1)
                its = 1
                cm(1) = cm(1) + am(1)
            case (5:6,9)
                its = 1
                bm(1) = -FOUR*sdts
                cm(1) = -bm(1)
        end select

        !
        ! Boundary condition at phi=pf
        !
        select case (mbdcnd)
            case (1,4:5)
                ct = cm(m)
                itf = m
            case (2:3,6)
                ct = cm(m+1)
                am(m+1) = am(m+1) + cm(m+1)
                itf = m + 1
            case (7:9)
                itf = m + 1
                am(m+1) = FOUR*sdts
                bm(m+1) = -am(m+1)
        end select

        wts = sint(its+1)*am(its+1)/cm(its)
        wtf = sint(itf-1)*cm(itf-1)/am(itf)
        itsp = its + 1
        itfm = itf - 1
        !
        ! Boundary condition at r=rs
        !
        ictr = 0
        select case (nbdcnd)
            case default
                ar = an(2)
                jrs = 2
            case (3:4)
                ar = an(1)
                jrs = 1
                cn(1) = cn(1) + an(1)
            case (5:6)
                jrs = 2
                ictr = 1
                s(n) = an(n)/bn(n)
                do j = 3, n
                    l = n - j + 2
                    s(l) = an(l)/(bn(l)-cn(l)*s(l+1))
                end do
                s(2) = -s(2)
                do j = 3, n
                    s(j) = -s(j)*s(j-1)
                end do
                wtnm = wts + wtf
                do i = itsp, itfm
                    wtnm = wtnm + sint(i)
                end do
                yps = czr*wtnm*(s(2)-ONE)
        end select

        !
        ! Boundary condition at r=rf
        !
        select case (nbdcnd)
            case (1,4:5)
                cr = cn(n)
                jrf = n
            case (2:3,6)
                cr = cn(n+1)
                an(n+1) = an(n+1) + cn(n+1)
                jrf = n + 1
        end select

        wrs = an(jrs+1)*r(jrs)**2/cn(jrs)
        wrf = cn(jrf-1)*r(jrf)**2/an(jrf)
        wrz = an(jrs)/czr
        jrsp = jrs + 1
        jrfm = jrf - 1
        munk = itf - its + 1
        nunk = jrf - jrs + 1
        bmh(its:itf) = bm(its:itf)
        ising = 0

        if  (nbdcnd == 3 .or. nbdcnd == 6 ) then
            select case (mbdcnd)
                case (3,6,8:9)
                    if (elmbda >= ZERO) then
                        ising = 1
                        summation = wts*wrs + wts*wrf + wtf*wrs + wtf*wrf
                        if (ictr /= 0) then
                            summation = summation + wrz
                        end if
                        do j = jrsp, jrfm
                            r2 = r(j)**2
                            do i = itsp, itfm
                                summation = summation + r2*sint(i)
                            end do
                        end do
                        do j = jrsp, jrfm
                            summation = summation + (wts + wtf)*r(j)**2
                        end do
                        do i = itsp, itfm
                            summation = summation + (wrs + wrf)*sint(i)
                        end do
                        hne = summation
                    end if
            end select
        end if

        select case (mbdcnd)
            case (1:4, 7:8)
                bm(its) = bmh(its) + elmbda/sint(its)**2
        end select

        select case (mbdcnd)
            case (1:6)
                bm(itf) = bmh(itf) + elmbda/sint(itf)**2
        end select

        bm(itsp:itfm) = bmh(itsp:itfm) + elmbda/sint(itsp:itfm)**2

        select case (mbdcnd)
            case (1:2, 7)
                f(2, jrs:jrf) = f(2, jrs:jrf) - at*f(1, jrs:jrf)/r(jrs:jrf)**2
            case (3:4, 8)
                f(1, jrs:jrf) = f(1, jrs:jrf) + tdt*bdts(jrs:jrf)*at/r(jrs:jrf)**2
        end select

        select case (mbdcnd)
            case (1, 4:5)
                f(m, jrs:jrf) = f(m, jrs:jrf) - ct*f(m+1, jrs:jrf)/r(jrs:jrf)**2
            case (2:3, 6)
                f(m+1, jrs:jrf)=f(m+1, jrs:jrf)-tdt*bdtf(jrs:jrf)*ct/r(jrs:jrf)**2
        end select

        case_block: block
            select case (nbdcnd)
                case default
                    if (mbdcnd /= 3) exit case_block
                    yhld = f(its, 1) - czr/tdt*(sin(tf)*bdtf(2)-sin(ts)*bdts(2))
                    f(:mp1, 1) = yhld
                case (1:2)
                    rs2 = (rs + dr)**2
                    f(its:itf, 2) = f(its:itf, 2) - ar*f(its:itf, 1)/rs2
                case (3:4)
                    f(its:itf, 1) = f(its:itf, 1) + tdr*bdrs(its:itf)*ar/rs**2
            end select
        end block case_block

        select case (nbdcnd)
            case (1, 4:5)
                rf2 = (rf - dr)**2
                f(its:itf, n) = f(its:itf, n) - cr*f(its:itf, n+1)/rf2
            case (2:3, 6)
                f(its:itf, n+1) = f(its:itf, n+1) - tdr*bdrf(its:itf)*cr/rf**2
        end select

        pertrb = ZERO

        if (ising /= 0) then
            summation = wts*wrs*f(its, jrs) + wts*wrf*f(its, jrf) + wtf*wrs*f(itf, &
                jrs) + wtf*wrf*f(itf, jrf)

            if (ictr /= 0) summation = summation + wrz*f(its, 1)

            do j = jrsp, jrfm
                r2 = r(j)**2
                do i = itsp, itfm
                    summation = summation + r2*sint(i)*f(i, j)
                end do
            end do

            summation = summation &
                + dot_product(r(jrsp:jrfm)**2, wts*f(its, jrsp:jrfm)+ &
                wtf*f(itf, jrsp:jrfm))

            summation = summation &
                + dot_product(sint(itsp:itfm), wrs*f(itsp:itfm, jrs)+ &
                wrf*f(itsp:itfm, jrf))

            pertrb = summation/hne

            f(:mp1, :np1) = f(:mp1, :np1) - pertrb

        end if

        do j = jrs, jrf
            rsq = r(j)**2
            f(its:itf, j) = rsq*f(its:itf, j)
        end do

        iflg = intl

        call util%blktrii(iflg, np, nunk, an(jrs:), bn(jrs:), cn(jrs:), mp, munk, &
            am(its:), bm(its:), cm(its:), idimf, f(its:, jrs:), ierror, w, wc)

        if (ierror /= 0) then
            error stop 'fishpack library: blktrii call failed in hwscsp_lower_routine'
        end if

        iflg = iflg + 1

        do while(iflg == 1)
            !
            ! Iterate solver
            !
            call util%blktrii(iflg, np, nunk, an(jrs:), bn(jrs:), cn(jrs:), mp, &
                munk, am(its:), bm(its:), cm(its:), idimf, f(its:, jrs:), ierror, w, wc)

            if (ierror /= 0) then
                error stop 'fishpack library: blktrii call failed in hwscsp_lower_routine'
            end if
            !
            ! Increment solver flag
            !
            iflg = iflg + 1
        end do

        if (nbdcnd == 0) f(:mp1, jrf+1) = f(:mp1, jrs)

        if (mbdcnd == 0) f(itf+1, :np1) = f(its, :np1)

        xp = ZERO

        if (ictr /= 0) then
            if (ising == 0) then
                summation = wts*f(its, 2) + wtf*f(itf, 2)
                summation = summation + dot_product(sint(itsp:itfm), f(itsp:itfm, 2))
                yph = czr*summation
                xp = (f(its, 1)-yph)/yps
                do j = jrs, jrf
                    xps = xp*s(j)
                    f(its:itf, j) = f(its:itf, j) + xps
                end do
            end if
            f(:mp1, 1) = xp
        end if

    end subroutine hwscsp_lower_routine

end submodule centered_axisymmetric_spherical_solver
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
!
