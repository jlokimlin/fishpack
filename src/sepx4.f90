!
!     file sepx4.f90
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright(c) 2005 by UCAR                   *
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
!     *                Boulder, Colorado (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     SUBROUTINE sepx4(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, d, n,
!    +                  nbdcnd, bdc, bdd, cofx, grhs, usol, idmn, pertrb,
!    +                  ierror)
!
!
!
! DIMENSION OF           bda(n+1), bdb(n+1), bdc(m+1), bdd(m+1),
! ARGUMENTS              usol(idmn, n+1),     grhs(idmn, n+1),
!
!
! LATEST REVISION        May 2016
!
! PURPOSE                sepx4 solves for either the second-order
!                        finite difference approximation or a
!                        fourth-order approximation to a separable
!                        elliptic equation
!
!                          af(x)*uxx+bf(x)*ux+cf(x)*u+uyy = g(x, y)
!
!                        on a rectangle(x greater than or equal to
!                        a and less than or equal to b, y greater than
!                        or equal to c and less than or equal to d).
!                        any combination of periodic or mixed boundary
!                        conditions is allowed.  if boundary
!                        conditions in the x direction are periodic
!                       (see mbdcnd=0 below) then the coefficients
!                        must satisfy
!
!                          af(x)=c1, bf(x)=0, cf(x)=c2 for all x.
!
!                        here c1, c2 are constants, c1>0.
!
!                        the possible boundary conditions are:
!                        in the x-direction:
!                         (0) periodic, u(x+b-a, y)=u(x, y) for
!                              all y, x
!                         (1) u(a, y), u(b, y) are specified for all y
!                         (2) u(a, y), du(b, y)/dx+beta*u(b, y) are
!                              specified for all y
!                         (3) du(a, y)/dx+alpha*u(a, y), du(b, y)/dx+
!                              beta*u(b, y) are specified for all y
!                         (4) du(a, y)/dx+alpha*u(a, y), u(b, y) are
!                              specified for all y
!
!                        in the y-direction:
!                         (0) periodic, u(x, y+d-c)=u(x, y) for all x, y
!                         (1) u(x, c), u(x, d) are specified for all x
!                         (2) u(x, c), du(x, d)/dy are specified for
!                              all x
!                         (3) du(x, c)/dy, du(x, d)/dy are specified for
!                              all x
!                         (4) du(x, c)/dy, u(x, d) are specified for
!                              all x
!
! USAGE                  call sepx4(iorder, a, b, m, mbdcnd, bda, alpha, bdb,
!                                   beta, c, d, n, nbdcnd, bdc, bdd, cofx,
!                                   grhs, usol, idmn, w, pertrb, ierror)
!
! ARGUMENTS
! ON INPUT               iorder
!                          = 2 if a second-order approximation is
!                              sought
!                          = 4 if a fourth-order approximation is
!                              sought
!
! *** CAUTION ***          grhs should be reset if sepx4 was first called
!                          with iorder=2 and will be called again with
!                          iorder=4.  values in grhs are destroyed by the
!                          iorder=2 call.
!
!
!                        a, b
!                          the range of the x-independent variable,
!                          i.e., x is greater than or equal to a
!                          and less than or equal to b.  a must be
!                          less than b.
!
!                        m
!                          the number of panels into which the
!                          interval(a, b) is subdivided.  hence,
!                          there will be m+1 grid points in the x-
!                          direction given by xi=a+(i-1)*dlx
!                          for i=1, 2, ..., m+1 where dlx=(b-a)/m is
!                          the panel width.  m must be less than
!                          idmn and greater than 5.
!
!                        mbdcnd
!                          indicates the type of boundary condition
!                          at x=a and x=b
!                          = 0 if the solution is periodic in x, i.e.,
!                              u(x+b-a, y)=u(x, y)  for all y, x
!                          = 1 if the solution is specified at x=a
!                              and x=b, i.e., u(a, y) and u(b, y) are
!                              specified for all y
!                          = 2 if the solution is specified at x=a
!                              and the boundary condition is mixed at
!                              x=b, i.e., u(a, y) and
!                              du(b, y)/dx+beta*u(b, y) are specified
!                              for all y
!                          = 3 if the boundary conditions at x=a and
!                              x=b are mixed, i.e.,
!                              du(a, y)/dx+alpha*u(a, y) and
!                              du(b, y)/dx+beta*u(b, y) are specified
!                              for all y
!                          = 4 if the boundary condition at x=a is
!                              mixed and the solution is specified
!                              at x=b, i.e., du(a, y)/dx+alpha*u(a, y)
!                              and u(b, y) are specified for all y
!
!                        bda
!                          a one-dimensional array of length n+1 that
!                          specifies the values of
!                          du(a, y)/dx+ alpha*u(a, y) at x=a, when
!                          mbdcnd=3 or 4.
!                          bda(j) = du(a, yj)/dx+alpha*u(a, yj),
!                          j=1, 2, ..., n+1
!                          when mbdcnd has any other value, bda is
!                          a dummy parameter.
!
!                        alpha
!                          the scalar multiplying the solution in case
!                          of a mixed boundary condition at x=a
!                         (see argument bda).  if mbdcnd is not equal
!                          to either 3 or 4, then alpha is a dummy
!                          parameter.
!
!                        bdb
!                          a one-dimensional array of length n+1 that
!                          specifies the values of
!                          du(b, y)/dx+ beta*u(b, y) at x=b.
!                          when mbdcnd=2 or 3
!                          bdb(j) = du(b, yj)/dx+beta*u(b, yj),
!                          j=1, 2, ..., n+1
!                          when mbdcnd has any other value, bdb is
!                          a dummy parameter.
!
!                        beta
!                          the scalar multiplying the solution in
!                          case of a mixed boundary condition at x=b
!                         (see argument bdb).  if mbdcnd is not equal
!                          to 2 or 3, then beta is a dummy parameter.
!
!                        c, d
!                          the range of the y-independent variable,
!                          i.e., y is greater than or equal to c and
!                          less than or equal to d.  c must be less
!                          than d.
!
!                        n
!                          the number of panels into which the
!                          interval(c, d) is subdivided.  hence,
!                          there will be n+1 grid points in the y-
!                          direction given by yj=c+(j-1)*dly for
!                          j=1, 2, ..., n+1 where dly=(d-c)/n is the
!                          panel width.  in addition, n must be
!                          greater than 4.
!
!                        nbdcnd
!                          indicates the types of boundary conditions
!                          at y=c and y=d
!                          = 0 if the solution is periodic in y,
!                              i.e., u(x, y+d-c)=u(x, y) for all x, y
!                          = 1 if the solution is specified at y=c
!                              and y = d, i.e., u(x, c)  and u(x, d)
!                              are specified for all x
!                          = 2 if the solution is specified at y=c
!                              and the boundary condition is mixed
!                              at y=d, i.e., du(x, c)/dy and u(x, d)
!                              are specified for all x
!                          = 3 if the boundary conditions are mixed
!                              at y=cand y=d i.e.,
!                              du(x, d)/dy and du(x, d)/dy are
!                              specified for all x
!                          = 4 if the boundary condition is mixed
!                              at y=c and the solution is specified
!                              at y=d, i.e. du(x, c)/dy+gama*u(x, c)
!                              and u(x, d) are specified for all x
!
!                        bdc
!                          a one-dimensional array of length m+1 that
!                          specifies the value du(x, c)/dy at y=c.
!
!                          when nbdcnd=3 or 4
!                            bdc(i) = du(xi, c)/dy i=1, 2, ..., m+1.
!
!                          when nbdcnd has any other value, bdc is
!                          a dummy parameter.
!
!                        bdd
!                          a one-dimensional array of length m+1 that
!                          specified the value of du(x, d)/dy at y=d.
!
!                          when nbdcnd=2 or 3
!                            bdd(i)=du(xi, d)/dy i=1, 2, ..., m+1.
!
!                          when nbdcnd has any other value, bdd is
!                          a dummy parameter.
!
!                        cofx
!                          a user-supplied subprogram with parameters
!                          x, afun, bfun, cfun which returns the
!                          values of the x-dependent coefficients
!                          af(x), bf(x), cf(x) in the elliptic
!                          equation at x.  if boundary conditions in
!                          the x direction are periodic then the
!                          coefficients must satisfy af(x)=c1, bf(x)=0,
!                          cf(x)=c2 for all x.  here c1>0
!                          and c2 are constants.
!
!                          note that cofx must be declared external
!                          in the calling routine.
!
!                        grhs
!                          a two-dimensional array that specifies the
!                          values of the right-hand side of the
!                          elliptic equation, i.e., grhs(i, j)=g(xi, yi),
!                          for i=2, ..., m, j=2, ..., n.  at the
!                          boundaries, grhs is defined by
!
!                          mbdcnd   grhs(1, j)   grhs(m+1, j)
!                          ------   ---------   -----------
!                            0      g(a, yj)     g(b, yj)
!                            1         *           *
!                            2         *        g(b, yj)  j=1, 2, ..., n+1
!                            3      g(a, yj)     g(b, yj)
!                            4      g(a, yj)        *
!
!                          nbdcnd   grhs(i, 1)   grhs(i, n+1)
!                          ------   ---------   -----------
!                            0      g(xi, c)     g(xi, d)
!                            1         *           *
!                            2         *        g(xi, d)  i=1, 2, ..., m+1
!                            3      g(xi, c)     g(xi, d)
!                            4      g(xi, c)        *
!
!                          where * means these quantites are not used.
!                          grhs should be dimensioned idmn by at least
!                          n+1 in the calling routine.
!
! *** CAUTION              grhs should be reset if sepx4 was first called
!                          with iorder=2 and will be called again with
!                          iorder=4.  values in grhs are destroyed by the
!                          iorder=2 call.
!
!                        usol
!                          a two-dimensional array that specifies the
!                          values of the solution along the boundaries.
!                          at the boundaries, usol is defined by
!
!                          mbdcnd   usol(1, j)   usol(m+1, j)
!                          ------   ---------   -----------
!                            0         *           *
!                            1      u(a, yj)     u(b, yj)
!                            2      u(a, yj)        *     j=1, 2, ..., n+1
!                            3         *           *
!                            4         *        u(b, yj)
!
!                          nbdcnd   usol(i, 1)   usol(i, n+1)
!                          ------   ---------   -----------
!                            0         *           *
!                            1      u(xi, c)     u(xi, d)
!                            2      u(xi, c)        *     i=1, 2, ..., m+1
!                            3         *           *
!                            4         *        u(xi, d)
!
!                          where * means the quantites are not used
!                          in the solution.
!
!                          if iorder=2, the user may equivalence grhs
!                          and usol to save space.  note that in this
!                          case the tables specifying the boundaries
!                          of the grhs and usol arrays determine the
!                          boundaries uniquely except at the corners.
!                          if the tables call for both g(x, y) and
!                          u(x, y) at a corner then the solution must
!                          be chosen.
!                          for example, if mbdcnd=2 and nbdcnd=4,
!                          then u(a, c), u(a, d), u(b, d) must be chosen
!                          at the corners in addition to g(b, c).
!
!                          if iorder=4, then the two arrays, usol and
!                          grhs, must be distinct.
!
!                          usol should be dimensioned idmn by at least
!                          n+1 in the calling routine.
!
!                        idmn
!                          the row(or first) dimension of the arrays
!                          grhs and usol as it appears in the program
!                          calling sepeli.  this parameter is used
!                          to specify the variable dimension of grhs
!                          and usol.  idmn must be at least 7 and
!                          greater than or equal to m+1.
!
!
! ON OUTPUT              usol
!                          contains the approximate solution to the
!                          elliptic equation. usol(i, j) is the
!                          approximation to u(xi, yj) for i=1, 2..., m+1
!                          and j=1, 2, ..., n+1.  the approximation has
!                          error o(dlx**2+dly**2) if called with
!                          iorder=2 and o(dlx**4+dly**4) if called
!                          with iorder=4.
!
!                        pertrb
!                          if a combination of periodic or derivative
!                          boundary conditions(i.e., alpha=beta=0 if
!                          mbdcnd=3) is specified and if cf(x)=0 for
!                          all x then a solution to the discretized
!                          matrix equation may not exist
!                         (reflecting the non-uniqueness of solutions
!                          to the pde).
!                          pertrb is a constant calculated and
!                          subtracted from the right hand side of the
!                          matrix equation insuring the existence of a
!                          solution.  sepx4 computes this solution
!                          which is a weighted minimal least squares
!                          solution to the original problem.  if
!                          singularity is not detected pertrb=ZERO is
!                          returned by sepx4.
!
!                        ierror
!                          an error flag that indicates invalid input
!                          parameters or failure to find a solution
!
!                          =  0 no error
!                          =  1 if a greater than b or c greater
!                               than d
!                          =  2 if mbdcnd less than 0 or mbdcnd
!                               greater than 4
!                          =  3 if nbdcnd less than 0 or nbdcnd
!                               greater than 4
!                          =  4 if attempt to find a solution fails.
!                              (the linear system generated is not
!                               diagonally dominant.)
!                          =  5 if idmn is too small(see discussion
!                               of idmn)
!                          =  6 if m is too small or too large
!                              (see discussion of m)
!                          =  7 if n is too small(see discussion of n)
!                          =  8 if iorder is not 2 or 4
!                          =  9 if intl is not 0 or 1
!                          = 10 if afun is less than or equal to zero
!                               for some interior mesh point xi some
!                               interior mesh point(xi, yj)
!                          = 12 if mbdcnd=0 and af(x)=cf(x)=constant
!                               or bf(x)=0 for all x is not true.
!                          = 20 if the dynamic allocation of real and
!                               complex workspace required for solution
!                               fails(for example if n, m are too large
!                               for your computer)
!
! SPECIAL CONDITIONS     None
!
! I/O                    None
!
! REQUIRED FILES         type_FishpackWorkspace.f90, genbun.f90, type_CyclicReductionUtility.f9090, type_SepAux.f90
!
!
! PRECISION              64-bit double precision
!
!
! STANDARD               Fortran 2008
!
! HISTORY                sepx4 was developed at NCAR by John C.
!                        Adams of the scientific computing division
!                        in October 1978.  The basis of this code is
!                        NCAR routine sepeli.  Both packages were
!                        released on NCAR's public libraries in
!                        January 1980. sepx4 was modified in June 2004
!                        incorporating Fortran 90 dynamical storage
!                        allocation for workspace requirements
!
!
! ALGORITHM              sepx4 automatically discretizes the separable
!                        elliptic equation which is then solved by a
!                        generalized cyclic reduction algorithm in the
!                        subroutine pois. The fourth order solution
!                        is obtained using the technique of defferred
!                        corrections referenced below.
!
! TIMING                 When possible, sepx4 should be used instead
!                        of package sepeli. The increase in speed
!                        is at least a factor of three.
!
! REFERENCES             Keller, H.B., Numerical methods for two-point
!                        boundary-value problems, BLAISDEL (1968),
!                        Waltham, Mass.
!
!                        Swarztrauber, P., and R. Sweet (1975):
!                        Efficient FORTRAN subprograms for the
!                        solution of elliptic partial differential
!                        equations.  NCAR Technical note
!                          NCAR-TN/IA-109, pp. 135-137.
!
module module_sepx4

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use centered_real_linear_systems_solver, only: &
        genbun

    use type_SepAux, only: &
        SepAux, &
        get_coefficients

    ! Explicit typing only!
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: sepx4

    type, private, extends(SepAux) :: Sepx4Aux
        !---------------------------------------------------------------
        ! Type components
        !---------------------------------------------------------------
        type(FishpackWorkspace), public :: workspace
        !---------------------------------------------------------------
    contains
        !---------------------------------------------------------------
        ! Type-bound procedures
        !---------------------------------------------------------------
        procedure, public  :: initialize_workspace
        procedure, public  :: s4elip
        procedure, private :: is_PDE_singular
        procedure, private :: defer
        !---------------------------------------------------------------
    end type Sepx4Aux

    !---------------------------------------------------------------
    ! Parameters confined to the module
    !---------------------------------------------------------------
    real(wp),    parameter :: ZERO = 0.0_wp
    real(wp),    parameter :: HALF = 0.5_wp
    real(wp),    parameter :: ONE = 1.0_wp
    real(wp),    parameter :: TWO = 2.0_wp
    integer(ip), parameter :: IIWK = 12 !! Size of workspace indices
    !---------------------------------------------------------------

contains

    subroutine sepx4(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, &
        d, n, nbdcnd, bdc, bdd, cofx, grhs, usol, idmn, pertrb, &
        ierror)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in)     :: iorder
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: mbdcnd
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: nbdcnd
        integer(ip), intent(in)     :: idmn
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: a
        real(wp),    intent(in)     :: b
        real(wp),    intent(in)     :: alpha
        real(wp),    intent(in)     :: beta
        real(wp),    intent(in)     :: c
        real(wp),    intent(in)     :: d
        real(wp),    intent(out)    :: pertrb
        real(wp),    intent(in)     :: bda(:)
        real(wp),    intent(in)     :: bdb(:)
        real(wp),    intent(in)     :: bdc(:)
        real(wp),    intent(in)     :: bdd(:)
        real(wp),    intent(inout)  :: grhs(:,:)
        real(wp),    intent(out)    :: usol(:,:)
        procedure(get_coefficients) :: cofx
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        type(Sepx4Aux) :: aux
        !--------------------------------------------------------------

        ! Check input parameters
        call check_input_parameters(iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, cofx, idmn, ierror)

        if (ierror /= 0) return

        ! Initialize workspace arrays and indices
        call aux%initialize_workspace(n, m, nbdcnd)

        ! Solve system
        associate( &
            i => aux%workspace%workspace_indices, &
            rew => aux%workspace%real_workspace &
            )
            associate( &
                an => rew(i(1):), &
                bn => rew(i(2):), &
                cn => rew(i(3):), &
                dn => rew(i(4):), &
                un => rew(i(5):), &
                zn => rew(i(6):), &
                am => rew(i(7):i(7)), &
                bm => rew(i(8):i(8)), &
                cm => rew(i(9):i(9)), &
                dm => rew(i(10):), &
                um => rew(i(11):), &
                zm => rew(i(12):) &
                )
                call aux%s4elip(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, &
                    c, d, n, nbdcnd, bdc, bdd, cofx, an, bn, cn, dn, un, zn, am, bm, &
                    cm, dm, um, zm, grhs, usol, idmn, pertrb, ierror)
            end associate
        end associate

        !
        ! Release memory
        !
        call aux%workspace%destroy()

    end subroutine sepx4

    subroutine initialize_workspace(self, n, m, nbdcnd)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(Sepx4Aux), intent(inout) :: self
        integer(ip), intent(in)        :: n, m, nbdcnd
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip)            :: l, k, length, irwk, icwk
        !--------------------------------------------------------------

        associate( w => self%workspace )
            ! Compute minimum workspace and check workspace length input
            select case (nbdcnd)
                case (0)
                    l = n
                    k = m + 1
                case default
                    l = n + 1
                    k = m + 1
            end select

            ! Compute required real and complex workspace sizes
            call compute_workspace_dimensions(n, l, k, length, irwk, icwk)

            ! Allocate memory for workspace arrays
            call w%create(irwk, icwk, IIWK)

            ! Set workspace indices
            w%workspace_indices = get_workspace_indices(length, l, k)
        end associate

    end subroutine initialize_workspace

    pure subroutine compute_workspace_dimensions(n, l, k, length, irwk, icwk)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: l
        integer(ip), intent(in)  :: k
        integer(ip), intent(out) :: length
        integer(ip), intent(out) :: irwk
        integer(ip), intent(out) :: icwk
        !--------------------------------------------------------------
        integer(ip) :: log2n
        !--------------------------------------------------------------

        log2n = int(log(real(n + 1, kind=wp))/log(TWO) + HALF, kind=ip)
        length = 4*(n + 1) +(10 + log2n) * k

        ! set real and complex workspace sizes
        irwk = length + 6 * (k + l) + 1
        icwk = 0

    end subroutine compute_workspace_dimensions

    pure function get_workspace_indices(length, l, k) result (return_value)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in) :: length
        integer(ip), intent(in) :: l
        integer(ip), intent(in) :: k
        integer(ip)             :: return_value(IIWK)
        !--------------------------------------------------------------
        integer(ip) :: j !! Counter
        !--------------------------------------------------------------

        associate( i => return_value)
            i(1) = length + 1

            do j = 1, 6
                i(j+1) = i(j) + l
            end do

            do j = 7, 11
                i(j+1) = i(j) + k
            end do
        end associate

    end function get_workspace_indices

    subroutine s4elip(self, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, &
        c, d, n, nbdcnd, bdc, bdd, cofx, an, bn, cn, dn, un, zn, am, bm, &
        cm, dm, um, zm, grhs, usol, idmn, pertrb, ierror)
        !
        ! Purpose:
        !
        !     s4elip sets up vectors and arrays for input to blktri
        !     and computes a second order solution in usol.  a return jump to
        !     sepeli occurrs if iorder=2.  if iorder=4 a fourth order
        !     solution is generated in usol.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(Sepx4Aux), intent(inout) :: self
        integer(ip),   intent(in)     :: iorder
        integer(ip),   intent(in)     :: m
        integer(ip),   intent(in)     :: mbdcnd
        integer(ip),   intent(in)     :: n
        integer(ip),   intent(in)     :: nbdcnd
        integer(ip),   intent(in)     :: idmn
        integer(ip),   intent(inout) :: ierror
        real(wp),      intent(in)     :: a
        real(wp),      intent(in)     :: b
        real(wp),      intent(in)     :: alpha
        real(wp),      intent(in)     :: beta
        real(wp),      intent(in)     :: c
        real(wp),      intent(in)     :: d
        real(wp),      intent(out)    :: pertrb
        real(wp),      intent(in)     :: bda(:)
        real(wp),      intent(in)     :: bdb(:)
        real(wp),      intent(in)     :: bdc(:)
        real(wp),      intent(in)     :: bdd(:)
        real(wp),      intent(inout) :: an(:)
        real(wp),      intent(inout) :: bn(:)
        real(wp),      intent(inout) :: cn(:)
        real(wp),      intent(inout) :: dn(:)
        real(wp),      intent(inout) :: un(:)
        real(wp),      intent(inout) :: zn(:)
        real(wp),      intent(inout) :: am(:)
        real(wp),      intent(inout) :: bm(:)
        real(wp),      intent(inout) :: cm(:)
        real(wp),      intent(inout) :: dm(:)
        real(wp),      intent(inout) :: um(:)
        real(wp),      intent(inout) :: zm(:)
        real(wp),      intent(inout) :: grhs(:,:)
        real(wp),      intent(inout) :: usol(:,:)
        procedure(get_coefficients)    :: cofx
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: i, i1, mp, np, local_error_flag
        real(wp)    :: xi, ai, bi, ci, axi, bxi, cxi
        real(wp)    :: dyj, eyj, fyj, ax1, cxm
        real(wp)    :: dy1, fyn, gama, xnu, prtrb
        logical      :: singular
        !-----------------------------------------------

        ! Associate various quantities
        associate( &
            kswx => self%kswx, &
            kswy => self%kswy, &
            k => self%k, &
            l=>self%l, &
            mit=>self%mit, &
            nit=> self%nit, &
            is=> self%is, &
            ms=> self%ms, &
            js=> self%js, &
            ns=> self%ns, &
            ait => self%ait, &
            bit => self%bit, &
            cit => self%cit, &
            dit => self%dit, &
            dlx => self%dlx, &
            dly => self%dly, &
            tdlx3 => self%tdlx3, &
            tdly3 => self%tdly3, &
            dlx4 => self%dlx4, &
            dly4 => self%dly4 &
            )

            !     set parameters internally
            !
            kswx = mbdcnd + 1
            kswy = nbdcnd + 1
            k = m + 1
            l = n + 1
            ait = a
            bit = b
            cit = c
            dit = d
            dly =(dit - cit)/n
            !
            ! set right hand side values from grhs in usol on the interior
            !     and non-specified boundaries.
            !
            usol(2:m, 2:n) = (dly**2) * grhs(2:m, 2:n)

            if (kswx /= 2 .and. kswx /= 3) then
                usol(1, 2:n) = (dly**2) * grhs(1, 2:n)
            end if

            if (kswx /= 2 .and. kswx /= 5) then
                usol(k, 2:n) = (dly**2) * grhs(k, 2:n)
            end if

            if (kswy /= 2 .and. kswy /= 3) then
                usol(2:m, 1) = (dly**2) * grhs(2:m, 1)
            end if

            if (kswy /= 2 .and. kswy /= 5) then
                usol(2:m, l) = (dly**2) * grhs(2:m, l)
            end if

            if (kswx /= 2 .and. kswx /= 3 .and. kswy /= 2 .and. kswy /= 3) then
                usol(1, 1) = (dly**2) * grhs(1, 1)
            end if

            if (kswx /= 2 .and. kswx /= 5 .and. kswy /= 2 .and. kswy /= 3) then
                usol(k, 1) = (dly**2) * grhs(k, 1)
            end if

            if (kswx /= 2 .and. kswx /= 3 .and. kswy /= 2 .and. kswy /= 5) then
                usol(1, l) = (dly**2) * grhs(1, l)
            end if

            if (kswx /= 2 .and. kswx /= 5 .and. kswy /= 2 .and. kswy /= 5) then
                usol(k, l) = (dly**2) * grhs(k, l)
            end if

            i1 = 1
            !
            ! set switches for periodic or non-periodic boundaries
            !
            if (kswx == 1) then
                mp = 0
            else
                mp = 1
            end if

            np = nbdcnd
            !
            !     set dlx, dly and size of block tri-diagonal system generated
            !     in nint, mint
            !
            dlx =(bit - ait)/m
            mit = k - 1

            if (kswx == 2) then
                mit = k - 2
            end if

            if (kswx == 4) then
                mit = k
            end if

            dly =(dit - cit)/n
            nit = l - 1

            if (kswy == 2) then
                nit = l - 2
            end if

            if (kswy == 4) then
                nit = l
            end if

            tdlx3 = TWO * (dlx**3)
            dlx4 = dlx**4
            tdly3 = TWO * (dly**3)
            dly4 = dly**4
            !
            ! set subscript limits for portion of array to input to blktri
            !
            if (kswx==2 .or. kswx==3) then
                is = 2
            else
                is = 1
            end if

            if (kswy==2 .or. kswy==3) then
                js = 2
            else
                js = 1
            end if

            ns = nit + js - 1
            ms = mit + is - 1
            !
            !     set x - direction
            !
            do i = 1, mit
                xi = ait + real(is + i - 2, kind=wp)*dlx
                call cofx(xi, ai, bi, ci)
                axi =(ai/dlx - HALF*bi)/dlx
                bxi =(-TWO * ai/dlx**2) + ci
                cxi =(ai/dlx + HALF*bi)/dlx
                am(i) = (dly**2) * axi
                bm(i) = (dly**2)*bxi
                cm(i) = (dly**2)*cxi
            end do
            !
            !     set y direction
            !
            dyj = ONE
            eyj = -TWO
            fyj = ONE
            an(:nit) = dyj
            bn(:nit) = eyj
            cn(:nit) = fyj
            !
            !     adjust edges in x direction unless periodic
            !
            ax1 = am(1)
            cxm = cm(mit)
            select case(kswx)
                case(2)
                    !
                    !     dirichlet-dirichlet in x direction
                    !
                    am(1) = ZERO
                    cm(mit) = ZERO
                case(3)
                    !
                    !     dirichlet-mixed in x direction
                    !
                    am(1) = ZERO
                    am(mit) = am(mit) + cxm
                    bm(mit) = bm(mit) - TWO * beta*dlx*cxm
                    cm(mit) = ZERO
                case(4)
                    !
                    !     mixed - mixed in x direction
                    !
                    am(1) = ZERO
                    bm(1) = bm(1) + TWO * dlx * alpha * ax1
                    cm(1) = cm(1) + ax1
                    am(mit) = am(mit) + cxm
                    bm(mit) = bm(mit) - TWO * dlx * beta * cxm
                    cm(mit) = ZERO
                case(5)
                    !
                    !     mixed-dirichlet in x direction
                    !
                    am(1) = ZERO
                    bm(1) = bm(1) + TWO * alpha * dlx * ax1
                    cm(1) = cm(1) + ax1
                    cm(mit) = ZERO
            end select
            !
            !     adjust in y direction unless periodic
            !
            dy1 = an(1)
            fyn = cn(nit)
            gama = ZERO
            xnu = ZERO
            select case(kswy)
                case(2)
                    !
                    !     dirichlet-dirichlet in y direction
                    !
                    an(1) = ZERO
                    cn(nit) = ZERO
                case(3)
                    !
                    !     dirichlet-mixed in y direction
                    !
                    an(1) = ZERO
                    an(nit) = an(nit) + fyn
                    bn(nit) = bn(nit) - TWO * dly * xnu * fyn
                    cn(nit) = ZERO
                case(4)
                    !
                    !     mixed - mixed direction in y direction
                    !
                    an(1) = ZERO
                    bn(1) = bn(1) + TWO * dly * gama * dy1
                    cn(1) = cn(1) + dy1
                    an(nit) = an(nit) + fyn
                    bn(nit) = bn(nit) - TWO * dly * xnu * fyn
                    cn(nit) = ZERO
                case(5)
                    !
                    !     mixed-dirichlet in y direction
                    !
                    an(1) = ZERO
                    bn(1) = bn(1) + TWO * dly * gama * dy1
                    cn(1) = cn(1) + dy1
                    cn(nit) = ZERO
            end select

            if (kswx /= 1) then
                !
                !     adjust usol along x edge
                !
                if (kswx==2 .or. kswx==3) then
                    if (kswx==2 .or. kswx==5) then
                        usol(is,js:ns) = usol(is,js:ns) - ax1*usol(1,js:ns)
                        usol(ms,js:ns) = usol(ms,js:ns) - cxm*usol(k,js:ns)
                    else
                        usol(is,js:ns) = usol(is,js:ns) - ax1*usol(1,js:ns)
                        usol(ms,js:ns) = usol(ms,js:ns) - TWO * dlx*cxm*bdb(js:ns)
                    end if
                else
                    if (kswx==2 .or. kswx==5) then
                        usol(is,js:ns) = usol(is,js:ns) + TWO * dlx*ax1*bda(js:ns)
                        usol(ms,js:ns) = usol(ms,js:ns) - cxm*usol(k,js:ns)
                    else
                        usol(is,js:ns) = usol(is,js:ns) + TWO * dlx*ax1*bda(js:ns)
                        usol(ms,js:ns) = usol(ms,js:ns) - TWO * dlx*cxm*bdb(js:ns)
                    end if
                end if
            end if
            if (kswy /= 1) then
                !
                !     adjust usol along y edge
                !
                if (kswy==2 .or. kswy==3) then
                    if (kswy==2 .or. kswy==5) then
                        usol(is:ms,js) = usol(is:ms,js) - dy1*usol(is:ms, 1)
                        usol(is:ms, ns) = usol(is:ms, ns) - fyn*usol(is:ms, l)
                    else
                        usol(is:ms,js) = usol(is:ms,js) - dy1*usol(is:ms, 1)
                        usol(is:ms, ns) = usol(is:ms, ns) - TWO * dly*fyn*bdd(is:ms)
                    end if
                else
                    if (kswy==2 .or. kswy==5) then
                        usol(is:ms,js) = usol(is:ms,js) + TWO * dly*dy1*bdc(is:ms)
                        usol(is:ms, ns) = usol(is:ms, ns) - fyn*usol(is:ms, l)
                    else
                        usol(is:ms,js) = usol(is:ms,js) + TWO * dly*dy1*bdc(is:ms)
                        usol(is:ms, ns) = usol(is:ms, ns) - TWO * dly*fyn*bdd(is:ms)
                    end if
                end if
            end if
            !
            ! save adjusted edges in grhs if iorder=4
            !
            if (iorder == 4) then
                grhs(is,js:ns) = usol(is,js:ns)
                grhs(ms,js:ns) = usol(ms,js:ns)
                grhs(is:ms,js) = usol(is:ms,js)
                grhs(is:ms, ns) = usol(is:ms, ns)
            end if

            pertrb = ZERO
            !
            ! check if operator is singular
            !
            call self%is_PDE_singular(mbdcnd, nbdcnd, alpha, beta, cofx, singular)
            !
            !     compute non-zero eigenvector in null space of transpose
            !     if singular
            !
            if (singular) then
                call self%septri(mit, am, bm, cm, dm, um, zm)
            end if

            if (singular) then
                call self%septri(nit, an, bn, cn, dn, un, zn)
            end if
            !
            !     adjust right hand side if necessary
            !
            if (singular) then
                call self%seport(usol, zn, zm, pertrb)
            end if

            !
            ! compute solution
            !
            !     save adjusted right hand side in grhs
            grhs(is:ms,js:ns) = usol(is:ms,js:ns)

            call genbun(np, nit, mp, mit, am, bm, cm, idmn, usol(is:,js:), local_error_flag)
            !
            !  Check if error detected in pois
            !     this can only correspond to ierror=12
            if (local_error_flag /= 0) then
                !       set error flag if improper coefficients input to pois
                ierror = 12
                return
            end if

            if (ierror /= 0) return
            !
            ! set periodic boundaries if necessary
            !
            if (kswx == 1) usol(k, :l) = usol(1, :l)

            if (kswy == 1) usol(:k, l) = usol(:k, 1)
            !
            !     minimize solution with respect to weighted least squares
            !     norm if operator is singular
            !
            if (singular) call self%sepmin(usol, zn, zm, prtrb)
            !
            !     return if deferred corrections and a fourth order solution are
            !     not flagged
            !
            if (iorder == 2) return
            !
            ! compute new right hand side for fourth order solution
            !
            call self%defer(cofx, idmn, usol, grhs)

            if (singular) call self%seport(usol, zn, zm, pertrb)
            !
            ! compute solution
            !
            !     save adjusted right hand side in grhs
            grhs(is:ms,js:ns) = usol(is:ms,js:ns)

            call genbun(np, nit, mp, mit, am, bm, cm, idmn, usol(is:,js:), local_error_flag)

            !
            ! check if error detected in pois
            !    this can only correspond to ierror=12
            !
            if (local_error_flag /= 0) then
                !       set error flag if improper coefficients input to pois
                ierror = 12
                return
            end if

            if (ierror /= 0) return
            !
            ! set periodic boundaries if necessary
            !
            if (kswx == 1) usol(k, :l) = usol(1, :l)
            if (kswy == 1) usol(:k, l) = usol(:k, 1)

            !
            !     minimize solution with respect to weighted least squares
            !     norm if operator is singular
            !
            if (singular) call self%sepmin(usol, zn, zm, prtrb)
        end associate

    end subroutine s4elip

    subroutine check_input_parameters(iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, cofx, &
        idmn, ierror)
        !
        ! Purpose:
        !
        ! This program checks the input parameters for errors
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in)    :: iorder
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: mbdcnd
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: nbdcnd
        integer(ip), intent(in)    :: idmn
        integer(ip), intent(out)   :: ierror
        real(wp),    intent(in)    :: a
        real(wp),    intent(in)    :: b
        real(wp),    intent(in)    :: c
        real(wp),    intent(in)    :: d
        procedure(get_coefficients) :: cofx
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: i
        real(wp)    :: xi, ai, bi, ci
        real(wp)    :: dlx
        !--------------------------------------------------------------

        if (a >= b .or. c >= d) then ! check definition of solution region
            ierror = 1
            return
        else if (mbdcnd < 0 .or. mbdcnd > 4) then ! check boundary switches
            ierror = 2
            return
        else if (nbdcnd < 0 .or. nbdcnd > 4) then
            ierror = 3
            return
        else if (idmn < 7) then ! check first dimension in calling routine
            ierror = 5
            return
        else if (m > idmn - 1 .or. m < 6) then ! check m
            ierror = 6
            return
        else if (n < 5) then ! check n
            ierror = 7
            return
        else if (iorder /= 2 .and. iorder /= 4) then ! Check iorder
            ierror = 8
            return
        end if
        !
        ! Check that equation is elliptic
        !
        dlx =(b - a)/m
        do i = 2, m
            xi = a + real(i - 1, kind=wp) * dlx
            call cofx(xi, ai, bi, ci)

            if (ai > ZERO) cycle

            ierror = 10
            return
        end do
        !
        ! no error found
        !
        ierror = 0


    end subroutine check_input_parameters

    subroutine is_PDE_singular(self, mbdcnd, nbdcnd, alpha, beta, cofx, singlr)
        !
        ! Purpose:
        !
        !     this subroutine checks if the pde sepeli
        !     must solve is a singular operator
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(Sepx4Aux), intent(inout) :: self
        integer(ip),   intent(in)     :: mbdcnd
        integer(ip),   intent(in)     :: nbdcnd
        real(wp),      intent(in)     :: alpha
        real(wp),      intent(in)     :: beta
        logical ,       intent(out)    :: singlr
        procedure(get_coefficients)    :: cofx
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: i
        real(wp)    :: xi, ai, bi, ci
        !--------------------------------------------------------------

        ! Associate various quantities
        associate( &
            kswx => self%kswx, &
            kswy => self%kswy, &
            k => self%k, &
            l=>self%l, &
            mit=>self%mit, &
            nit=> self%nit, &
            is=> self%is, &
            ms=> self%ms, &
            js=> self%js, &
            ns=> self%ns, &
            ait => self%ait, &
            bit => self%bit, &
            cit => self%cit, &
            dit => self%dit, &
            dlx => self%dlx, &
            dly => self%dly, &
            tdlx3 => self%tdlx3, &
            tdly3 => self%tdly3, &
            dlx4 => self%dlx4, &
            dly4 => self%dly4 &
            )

            singlr = .false.
            !
            !     check if the boundary conditions are
            !     entirely periodic and/or mixed
            !
            if (mbdcnd /=0 .and. mbdcnd /=3 .or. nbdcnd /=0 .and. nbdcnd /= 3) return
            !
            !     check that mixed conditions are pure neuman
            !
            if (mbdcnd == 3 .and. (alpha /= ZERO .or. beta /= ZERO)) return
            !
            !     check that non-derivative coefficient functions
            !     are zero
            !
            do i = is, ms
                xi = ait + real(i - 1, kind=wp)*dlx
                call cofx(xi, ai, bi, ci)
                if (ci == ZERO) cycle
                return
            end do
            !
            !     the operator must be singular if this point is reached
            !
            singlr = .true.

        end associate

    end subroutine is_PDE_singular

    subroutine defer(self, cofx, idmn, usol, grhs)
        !
        ! Purpose:
        !
        !     this subroutine first approximates the truncation error given by
        !     trun1(x, y)=dlx**2*tx+dly**2*ty where
        !     tx=afun(x)*uxxxx/12 + bfun(x)*uxxx/6 on the interior and
        !     at the boundaries if periodic(here uxxx, uxxxx are the third
        !     and fourth partial derivatives of u with respect to x).
        !     tx is of the form afun(x)/3 * (uxxxx/4+uxxx/dlx)
        !     at x=a or x=b if the boundary condition there is mixed.
        !     tx=ZERO along specified boundaries.  ty has symmetric form
        !     in y with x, afun(x), bfun(x) replaced by y, dfun(y), efun(y).
        !     the second order solution in usol is used to approximate
        !    (via second order finite differencing) the truncation error
        !     and the result is added to the right hand side in grhs
        !     and then transferred to usol to be used as a new right
        !     hand side when calling blktri for a fourth order solution.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(Sepx4Aux), intent(inout) :: self
        integer(ip),   intent(in)     :: idmn
        real(wp),      intent(inout) :: usol(:,:)
        real(wp),      intent(inout) :: grhs(:,:)
        procedure(get_coefficients)    :: cofx
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: i, j
        real(wp)    :: xi, ai, bi, ci
        real(wp)    :: uxxx, uxxxx, uyyy, uyyyy, tx, ty
        !--------------------------------------------------------------

        ! Associate various quantities
        associate( &
            kswx => self%kswx, &
            kswy => self%kswy, &
            k => self%k, &
            l=>self%l, &
            mit=>self%mit, &
            nit=> self%nit, &
            is=> self%is, &
            ms=> self%ms, &
            js=> self%js, &
            ns=> self%ns, &
            ait => self%ait, &
            bit => self%bit, &
            cit => self%cit, &
            dit => self%dit, &
            dlx => self%dlx, &
            dly => self%dly, &
            tdlx3 => self%tdlx3, &
            tdly3 => self%tdly3, &
            dlx4 => self%dlx4, &
            dly4 => self%dly4 &
            )

            !     compute truncation error approximation over the entire mesh
            !
            do i = is, ms
                xi = ait + real(i - 1, kind=wp)*dlx
                call cofx(xi, ai, bi, ci)
                do j = js, ns
                    !
                    !     compute partial derivative approximations at(xi, yj)
                    !
                    call self%sepdx(usol, i, j, uxxx, uxxxx)
                    call self%sepdy(usol, idmn, i, j, uyyy, uyyyy)
                    tx = ai*(uxxxx/12) + bi*(uxxx/6)
                    ty = uyyyy/12
                    !
                    !     reset form of truncation if at boundary which is non-periodic
                    !
                    if (kswx /= 1 .and. (i==1 .or. i==k)) then
                        tx = (ai/3) * ((uxxxx/4) + uxxx/dlx)
                    end if

                    if (kswy /= 1 .and. (j==1 .or. j==l)) then
                        ty = ((uyyyy/4)+uyyy/dly)/3
                    end if

                    grhs(i, j) = grhs(i, j) + (dly**2)*((dlx**2)*tx + (dly**2)*ty)
                end do
            end do
            !
            !     reset the right hand side in usol
            !
            usol(is:ms,js:ns) = grhs(is:ms,js:ns)

        end associate

    end subroutine defer

end module module_sepx4
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
