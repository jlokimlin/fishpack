!
!     file sepeli.f90
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
!     SUBROUTINE sepeli(intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c,
!                       d, n, nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, grhs,
!                       usol, idmn, workspace, pertrb, ierror)
!
! DIMENSION OF           bda(n+1), bdb(n+1), bdc(m+1), bdd(m+1),
! ARGUMENTS              usol(idmn, n+1), grhs(idmn, n+1),
!
! LATEST REVISION        April 2016
!
! PURPOSE                sepeli solves for either the second-order
!                        finite difference approximation or a
!                        fourth-order approximation to a separable
!                        elliptic equation
!
!                                 2    2
!                          af(x)*d u/dx + bf(x)*du/dx  + cf(x)*u +
!                                 2    2
!                          df(y)*d u/dy  + ef(y)*du/dy + ff(y)*u
!
!                          = g(x, y)
!
!                        on a rectangle (x greater than or equal to a
!                        and less than or equal to b; y greater than
!                        or equal to c and less than or equal to d).
!                        any combination of periodic or mixed boundary
!                        conditions is allowed.
!
!                        The possible boundary conditions are:
!                        in the x-direction:
!                        (0) periodic, u(x+b-a, y)=u(x, y) for all
!                            y, x (1) u(a, y), u(b, y) are specified for
!                            all y
!                        (2) u(a, y), du(b, y)/dx+beta*u(b, y) are
!                            specified for all y
!                        (3) du(a, y)/dx+alpha*u(a, y), du(b, y)/dx+
!                            beta*u(b, y) are specified for all y
!                        (4) du(a, y)/dx+alpha*u(a, y), u(b, y) are
!                            specified for all y
!
!                        in the y-direction:
!                        (0) periodic, u(x, y+d-c)=u(x, y) for all x, y
!                        (1) u(x, c), u(x, d) are specified for all x
!                        (2) u(x, c), du(x, d)/dy+xnu*u(x, d) are
!                            specified for all x
!                        (3) du(x, c)/dy+gama*u(x, c), du(x, d)/dy+
!                            xnu*u(x, d) are specified for all x
!                        (4) du(x, c)/dy+gama*u(x, c), u(x, d) are
!                            specified for all x
!
! USAGE                  call sepeli (intl, iorder, a, b, m, mbdcnd, bda,
!                                     alpha, bdb, beta, c, d, n, nbdcnd, bdc,
!                                     gama, bdd, xnu, cofx, cofy, grhs, usol,
!                                     idmn, w, pertrb, ierror)
!
! ARGUMENTS
! ON INPUT               intl
!                          = 0 on initial entry to sepeli or if any
!                              of the arguments c, d, n, nbdcnd, cofy
!                              are changed from a previous call
!                          = 1 if c, d, n, nbdcnd, cofy are unchanged
!                              from the previous call.
!
!                        iorder
!                          = 2 if a second-order approximation
!                              is sought
!                          = 4 if a fourth-order approximation
!                              is sought
!
!                        a, b
!                          the range of the x-independent variable,
!                          i.e., x is greater than or equal to a
!                          and less than or equal to b.  a must be
!                          less than b.
!
!                        m
!                          the number of panels into which the
!                          interval [a, b] is subdivided. hence,
!                          there will be m+1 grid points in the x-
!                          direction given by xi=a+(i-1)*dlx
!                          for i=1, 2, ..., m+1 where dlx=(b-a)/m is
!                          the panel width.  m must be less than
!                          idmn and greater than 5.
!
!                        mbdcnd
!                          indicates the type of boundary condition
!                          at x=a and x=b
!
!                          = 0 if the solution is periodic in x, i.e.,
!                              u(x+b-a, y)=u(x, y)  for all y, x
!                          = 1 if the solution is specified at x=a
!                              and x=b, i.e., u(a, y) and u(b, y) are
!                              specified for all y
!                          = 2 if the solution is specified at x=a and
!                              the boundary condition is mixed at x=b,
!                              i.e., u(a, y) and du(b, y)/dx+beta*u(b, y)
!                              are specified for all y
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
!                          a one-dimensional array of length n+1
!                          that specifies the values of
!                          du(a, y)/dx+ alpha*u(a, y) at x=a, when
!                          mbdcnd=3 or 4.
!                          bda(j) = du(a, yj)/dx+alpha*u(a, yj),
!                          j=1, 2, ..., n+1. when mbdcnd has any other
!                          other value, bda is a dummy parameter.
!
!                        alpha
!                          the scalar multiplying the solution in
!                          case of a mixed boundary condition at x=a
!                          (see argument bda).  if mbdcnd is not
!                          equal to 3 or 4 then alpha is a dummy
!                          parameter.
!
!                        bdb
!                          a one-dimensional array of length n+1
!                          that specifies the values of
!                          du(b, y)/dx+ beta*u(b, y) at x=b.
!                          when mbdcnd=2 or 3
!                          bdb(j) = du(b, yj)/dx+beta*u(b, yj),
!                          j=1, 2, ..., n+1. when mbdcnd has any other
!                          other value, bdb is a dummy parameter.
!
!                        beta
!                          the scalar multiplying the solution in
!                          case of a mixed boundary condition at
!                          x=b (see argument bdb).  if mbdcnd is
!                          not equal to 2 or 3 then beta is a dummy
!                          parameter.
!
!                        c, d
!                          the range of the y-independent variable,
!                          i.e., y is greater than or equal to c
!                          and less than or equal to d.  c must be
!                          less than d.
!
!                        n
!                          the number of panels into which the
!                          interval [c, d] is subdivided.
!                          hence, there will be n+1 grid points
!                          in the y-direction given by
!                          yj=c+(j-1)*dly for j=1, 2, ..., n+1 where
!                          dly=(d-c)/n is the panel width.
!                          in addition, n must be greater than 4.
!
!                        nbdcnd
!                          indicates the types of boundary conditions
!                          at y=c and y=d
!
!                          = 0 if the solution is periodic in y,
!                              i.e., u(x, y+d-c)=u(x, y)  for all x, y
!                          = 1 if the solution is specified at y=c
!                              and y = d, i.e., u(x, c) and u(x, d)
!                              are specified for all x
!                          = 2 if the solution is specified at y=c
!                              and the boundary condition is mixed
!                              at y=d, i.e., u(x, c) and
!                              du(x, d)/dy+xnu*u(x, d) are specified
!                              for all x
!                          = 3 if the boundary conditions are mixed
!                              at y=c and y=d, i.e.,
!                              du(x, d)/dy+gama*u(x, c) and
!                              du(x, d)/dy+xnu*u(x, d) are specified
!                              for all x
!                          = 4 if the boundary condition is mixed
!                              at y=c and the solution is specified
!                              at y=d, i.e. du(x, c)/dy+gama*u(x, c)
!                              and u(x, d) are specified for all x
!
!                        bdc
!                          a one-dimensional array of length m+1
!                          that specifies the value of
!                          du(x, c)/dy+gama*u(x, c) at y=c.
!                          when nbdcnd=3 or 4 bdc(i) = du(xi, c)/dy +
!                          gama*u(xi, c), i=1, 2, ..., m+1.
!                          when nbdcnd has any other value, bdc
!                          is a dummy parameter.
!
!                        gama
!                          the scalar multiplying the solution in
!                          case of a mixed boundary condition at
!                          y=c (see argument bdc).  if nbdcnd is
!                          not equal to 3 or 4 then gama is a dummy
!                          parameter.
!
!                        bdd
!                          a one-dimensional array of length m+1
!                          that specifies the value of
!                          du(x, d)/dy + xnu*u(x, d) at y=c.
!                          when nbdcnd=2 or 3 bdd(i) = du(xi, d)/dy +
!                          xnu*u(xi, d), i=1, 2, ..., m+1.
!                          when nbdcnd has any other value, bdd
!                          is a dummy parameter.
!
!                        xnu
!                          the scalar multiplying the solution in
!                          case of a mixed boundary condition at
!                          y=d (see argument bdd).  if nbdcnd is
!                          not equal to 2 or 3 then xnu is a
!                          dummy parameter.
!
!                        cofx
!                          a user-supplied subprogram with
!                          parameters x, afun, bfun, cfun which
!                          returns the values of the x-dependent
!                          coefficients af(x), bf(x), cf(x) in the
!                          elliptic equation at x.
!
!                        cofy
!                          a user-supplied subprogram with parameters
!                          y, dfun, efun, ffun which returns the
!                          values of the y-dependent coefficients
!                          df(y), ef(y), ff(y) in the elliptic
!                          equation at y.
!
!                          note:  cofx and cofy must be declared
!                          external in the calling routine.
!                          the values returned in afun and dfun
!                          must satisfy afun*dfun greater than 0
!                          for a less than x less than b, c less
!                          than y less than d (see ierror=10).
!                          the coefficients provided may lead to a
!                          matrix equation which is not diagonally
!                          dominant in which case solution may fail
!                          (see ierror=4).
!
!                        grhs
!                          a two-dimensional array that specifies the
!                          values of the right-hand side of the
!                          elliptic equation, i.e.,
!                          grhs(i, j)=g(xi, yi), for i=2, ..., m,
!                          j=2, ..., n.  at the boundaries, grhs is
!                          defined by
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
!                          where * means these quantities are not used.
!                          grhs should be dimensioned idmn by at least
!                          n+1 in the calling routine.
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
!                          where * means the quantities are not used
!                          in the solution.
!
!                          if iorder=2, the user may equivalence grhs
!                          and usol to save space.  note that in this
!                          case the tables specifying the boundaries
!                          of the grhs and usol arrays determine the
!                          boundaries uniquely except at the corners.
!                          if the tables call for both g(x, y) and
!                          u(x, y) at a corner then the solution must
!                          be chosen.  for example, if mbdcnd=2 and
!                          nbdcnd=4, then u(a, c), u(a, d), u(b, d) must
!                          be chosen at the corners in addition
!                          to g(b, c).
!
!                          if iorder=4, then the two arrays, usol and
!                          grhs, must be distinct.
!
!                          usol should be dimensioned idmn by at least
!                          n+1 in the calling routine.
!
!                        idmn
!                          the row (or first) dimension of the arrays
!                          grhs and usol as it appears in the program
!                          calling sepeli.  this parameter is used
!                          to specify the variable dimension of grhs
!                          and usol.  idmn must be at least 7 and
!                          greater than or equal to m+1.
!
!                        workspace
!                          An object of class(FishpackWorkspace) variable
!                          which is used internally in sepeli to dynamically
!                          allocate real and complex workspace arrays used
!                          in the solver. An error flag (ierror = 20) is
!                          set if the required workspace
!                          allocation fails (for example if n, m are too large)
!                          real and complex values are set in the components
!                          of workspace on a initial (intl=0) call to sepeli.
!                          These must be preserved on non-initial calls (intl=1)
!                          to sepeli. This eliminates redundant calculations
!                          and saves compute time.
!
!               ****       IMPORTANT!  The user program calling sepeli should
!                          include the statement:
!
!                               call workspace%destroy()
!
!                          after the final approximation is generated by
!                          sepeli. This will deallocate the real and complex
!                          workspace of workspace. Failure to include this statement
!                          could result in serious memory leakage.
!
! ON OUTPUT              usol
!                          Contains the approximate solution to the
!                          elliptic equation.
!                          usol(i, j) is the approximation to u(xi, yj)
!                          for i=1, 2..., m+1 and j=1, 2, ..., n+1.
!                          the approximation has error
!                          o(dlx**2+dly**2) if called with iorder=2
!                          and o(dlx**4+dly**4) if called with
!                          iorder=4.
!
!                        workwpace
!                          The derived type(FishpackWorkspace) variable
!                          contains real and complex values that must not
!                          be destroyed if sepeli is called again with
!                          intl=1.
!
!                        pertrb
!                          if a combination of periodic or derivative
!                          boundary conditions
!                          (i.e., alpha=beta=0 if mbdcnd=3;
!                          gama=xnu=0 if nbdcnd=3) is specified
!                          and if the coefficients of u(x, y) in the
!                          separable elliptic equation are zero
!                          (i.e., cf(x)=0 for x greater than or equal
!                          to a and less than or equal to b;
!                          ff(y)=0 for y greater than or equal to c
!                          and less than or equal to d) then a
!                          solution may not exist.  pertrb is a
!                          constant calculated and subtracted from
!                          the right-hand side of the matrix equations
!                          generated by sepeli which insures that a
!                          solution exists. sepeli then computes this
!                          solution which is a weighted minimal least
!                          squares solution to the original problem.
!
!                        ierror
!                          an error flag that indicates invalid input
!                          parameters or failure to find a solution
!                          = 0 no error
!                          = 1 if a greater than b or c greater than d
!                          = 2 if mbdcnd less than 0 or mbdcnd greater
!                              than 4
!                          = 3 if nbdcnd less than 0 or nbdcnd greater
!                              than 4
!                          = 4 if attempt to find a solution fails.
!                              (the linear system generated is not
!                              diagonally dominant.)
!                          = 5 if idmn is too small
!                              (see discussion of idmn)
!                          = 6 if m is too small or too large
!                              (see discussion of m)
!                          = 7 if n is too small (see discussion of n)
!                          = 8 if iorder is not 2 or 4
!                          = 9 if intl is not 0 or 1
!                          = 10 if afun*dfun less than or equal to 0
!                               for some interior mesh point (xi, yj)
!                          = 20 If the dynamic allocation of real and
!                               complex workspace in the derived type
!                               (FishpackWorkspace) variable W fails (e.g.,
!                               if N, M are too large for the platform used)
!
!                          Note (concerning ierror=4):  for the
!                          coefficients input through cofx, cofy,
!                          the discretization may lead to a block
!                          tridiagonal linear system which is not
!                          diagonally dominant (for example, this
!                          happens if cfun=0 and bfun/(TWO *dlx) greater
!                          than afun/dlx**2).  in this case solution
!                          may fail.  this cannot happen in the limit
!                          as dlx, dly approach zero.  hence, the
!                          condition may be remedied by taking larger
!                          values for m or n.
!
! SPECIAL CONDITIONS     See cofx, cofy argument descriptions above.
!
! I/O                    None
!
! PRECISION              Set by the instrinsic module ISO_Fortran_env to 64-bit double precision
!
! REQUIRED FILES         blktri.f90, type_SepAux.f90, type_FishpackWorkspace.f90
!
! STANDARD               Fortran 2008
!
! HISTORY                Developed at NCAR during 1975-76 by
!                        John c. Adams of the scientific computing
!                        division.  Released on NCAR's public software
!                        libraries in January 1980. Revised in June
!                        2004 using Fortan 90 dynamically allocated work
!                        space and derived data types to eliminate mixed
!                        mode conflicts in the earlier versions. All
!                        statement labels, arithmetic if statements and
!                        computed go to statements have been removed from
!                        the current version of sepeli.
!
! ALGORITHM              sepeli automatically discretizes the
!                        separable elliptic equation which is then
!                        solved by a generalized cyclic reduction
!                        algorithm in the subroutine, blktri. The
!                        fourth-order solution is obtained using
!                        'deferred corrections' which is described
!                        and referenced in sections, references and
!                        method.
!
! TIMING                 The operational count is proportional to
!                        m*n*log2(n).
!
! ACCURACY               The following accuracy results were obtained
!                        using 64 bit floating point arithmetic. note
!                        that the fourth-order accuracy is not realized
!                        until the mesh is sufficiently refined.
!
!                                     second-order  fourth-order
!                            m    n     error         error
!
!                             6    6    6.8e-1        1.2e0
!                            14   14    1.4e-1        1.8e-1
!                            30   30    3.2e-2        9.7e-3
!                            62   62    7.5e-3        3.0e-4
!                           126  126    1.8e-3        3.5e-6
!
!
! REFERENCES             Keller, H.B., Numerical methods for two-point
!                        boundary-value problems, Blaisdel (1968),
!                        Waltham, Mass.
!
!                        Swarztrauber, P., and R. Sweet (1975):
!                        Efficient FORTRAN subprograms for the
!                        solution of elliptic partial differential
!                        equations.  NCAR Technical note
!                        NCAR-TN/IA-109, PP. 135-137.
!
module module_sepeli

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_GeneralizedCyclicReductionUtility, only: &
        GeneralizedCyclicReductionUtility

    use type_SepAux, only: &
        SepAux, &
        get_coefficients

    ! Explicit typing only!
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: sepeli

    type, private, extends(SepAux) :: SepeliAux
        !---------------------------------------------------------------
        ! Type components
        !---------------------------------------------------------------
        type(GeneralizedCyclicReductionUtility), private :: blktri_aux
        !---------------------------------------------------------------
    contains
        !---------------------------------------------------------------
        ! Type-bound procedures
        !---------------------------------------------------------------
        procedure, public  :: spelip
        procedure, private :: is_PDE_singular
        procedure, private :: defer
        !---------------------------------------------------------------
    end type SepeliAux

    !---------------------------------------------------------------
    ! Parameters confined to the module
    !---------------------------------------------------------------
    real(wp),    parameter :: ZERO = 0.0_wp
    real(wp),    parameter :: HALF = 0.5_wp
    real(wp),    parameter :: TWO = 2.0_wp
    integer(ip), parameter :: IIWK = 12 !! Size of workspace indices
    !---------------------------------------------------------------

contains

    subroutine sepeli(intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, &
        beta, c, d, n, nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, grhs, &
        usol, idmn, workspace, pertrb, ierror)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in)     :: intl
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
        real(wp),    intent(in)     :: gama
        real(wp),    intent(in)     :: xnu
        real(wp),    intent(out)    :: pertrb
        real(wp),    intent(in)     :: bda(:)
        real(wp),    intent(in)     :: bdb(:)
        real(wp),    intent(in)     :: bdc(:)
        real(wp),    intent(in)     :: bdd(:)
        real(wp),    intent(inout)  :: grhs(:,:)
        real(wp),    intent(inout)  :: usol(:,:)
        class(FishpackWorkspace), intent(inout)  :: workspace
        procedure(get_coefficients) :: cofx
        procedure(get_coefficients) :: cofy
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        type(SepeliAux), save :: aux
        !--------------------------------------------------------------

        !
        ! Check input arguments
        !
        call check_input_arguments(intl, iorder, a, b, m, mbdcnd, c, d, n, &
            nbdcnd, cofx, cofy, idmn, ierror)

        ! Check error flag
        if (ierror /= 0) return

        !
        ! allocate workspace arrays on initial call only
        !
        if (intl == 0) call initialize_workspace(n, m, workspace)

        associate( &
            indx => workspace%workspace_indices, &
            rew => workspace%real_workspace, &
            cxw => workspace%complex_workspace &
            )

            !
            ! Compute 2nd or 4th order solution
            !
            call aux%spelip(intl, iorder, a, b, m, mbdcnd, &
                bda, alpha, bdb, beta, c, d, n, &
                nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, &
                rew(indx(1):), rew(indx(2):), rew(indx(3):), rew(indx(4):), &
                rew(indx(5):), rew(indx(6):), rew(indx(7):), rew(indx(8):), &
                rew(indx(9):), rew(indx(10):), rew(indx(11):), rew(indx(12):), &
                grhs, usol, idmn, rew, cxw, pertrb, ierror)

        end associate

    end subroutine sepeli

    subroutine initialize_workspace(n, m, workspace)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: m
        class(FishpackWorkspace), intent(out) :: workspace
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: irwk, icwk, indx(IIWK)
        !-----------------------------------------------

        ! Compute required blktri workspace lengths
        call workspace%compute_blktri_workspace_lengths(n, m, irwk, icwk)

        ! TODO **************
        ! Try to eliminate this local variable altogether
        ! Compute workspace indices
        indx = get_workspace_indices(irwk, n, m)

        ! Adjust workspace requirements for sepeli
        irwk = indx(12) + m + 1
        icwk = icwk + 3 * (m + 1)

        ! Allocate required memory
        call workspace%create(irwk, icwk, IIWK)

        ! Set workspace indices
        workspace%workspace_indices = indx

    end subroutine initialize_workspace

    pure function get_workspace_indices(irwk, n, m) result (return_value)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in) :: irwk
        integer(ip), intent(in) :: n
        integer(ip), intent(in) :: m
        integer(ip)             :: return_value(IIWK)
        !--------------------------------------------------------------
        integer(ip) :: j !! Counter
        !--------------------------------------------------------------

        associate( indx => return_value)

            indx(1) = irwk + 1

            do j = 1, 6
                indx(j+1) = indx(j) + n + 1
            end do

            do j = 7, 11
                indx(j+1) = indx(j) + m + 1
            end do

        end associate

    end function get_workspace_indices

    subroutine spelip(self, intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, &
        beta, c, d, n, nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, an, bn, &
        cn, dn, un, zn, am, bm, cm, dm, um, zm, grhs, usol, idmn, w, &
        wc, pertrb, ierror)
        !
        ! Purpose:
        !
        !     spelip sets up vectors and arrays for input to blktri
        !     and computes a second order solution in usol.  a return jump to
        !     sepeli occurrs if iorder=2.  if iorder=4 a fourth order
        !     solution is generated in usol.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(SepeliAux), intent(inout) :: self
        integer(ip),   intent(in)     :: intl
        integer(ip),   intent(in)     :: iorder
        integer(ip),   intent(in)     :: m
        integer(ip),   intent(in)     :: mbdcnd
        integer(ip),   intent(in)     :: n
        integer(ip),   intent(in)     :: nbdcnd
        integer(ip),   intent(in)     :: idmn
        integer(ip),   intent(out)    :: ierror
        real(wp),      intent(in)     :: a
        real(wp),      intent(in)     :: b
        real(wp),      intent(in)     :: alpha
        real(wp),      intent(in)     :: beta
        real(wp),      intent(in)     :: c
        real(wp),      intent(in)     :: d
        real(wp),      intent(in)     :: gama
        real(wp),      intent(in)     :: xnu
        real(wp),      intent(out)    :: pertrb
        real(wp),      intent(in)     :: bda(:)
        real(wp),      intent(in)     :: bdb(:)
        real(wp),      intent(in)     :: bdc(:)
        real(wp),      intent(in)     :: bdd(:)
        real(wp),      intent(out)    :: an(:)
        real(wp),      intent(out)    :: bn(:)
        real(wp),      intent(out)    :: cn(:)
        real(wp),      intent(out)    :: dn(:)
        real(wp),      intent(out)    :: un(:)
        real(wp),      intent(out)    :: zn(:)
        real(wp),      intent(out)    :: am(:)
        real(wp),      intent(out)    :: bm(:)
        real(wp),      intent(out)    :: cm(:)
        real(wp),      intent(out)    :: dm(:)
        real(wp),      intent(out)    :: um(:)
        real(wp),      intent(out)    :: zm(:)
        real(wp),      intent(inout)  :: grhs(:,:)
        real(wp),      intent(inout)  :: usol(:,:)
        real(wp),      intent(out)    :: w(:)
        complex(wp),   intent(out)    :: wc(:)
        procedure(get_coefficients)   :: cofx
        procedure(get_coefficients)   :: cofy
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: i, j, i1, mp, np
        real(wp)    :: xi, ai, bi, ci, axi, bxi, cxi
        real(wp)    :: yj, dj, ej, fj, dyj, eyj
        real(wp)    :: fyj, ax1, cxm, dy1, fyn, prtrb
        logical     :: singular
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
            !
            ! set right hand side values from grhs in usol on the interior
            !    and non-specified boundaries.
            !
            usol(2:m, 2:n) = grhs(2:m, 2:n)

            if (kswx /= 2 .and. kswx /= 3) then
                usol(1, 2:n) = grhs(1, 2:n)
            end if

            if (kswx /= 2 .and. kswx /= 5) then
                usol(k, 2:n) = grhs(k, 2:n)
            end if

            if (kswy /= 2 .and. kswy /= 3) then
                usol(2:m, 1) = grhs(2:m, 1)
            end if

            if (kswy /= 2 .and. kswy /= 5) then
                usol(2:m, l) = grhs(2:m, l)
            end if

            if (kswx/=2 .and. kswx/=3 .and. kswy/=2 .and. kswy/=3) then
                usol(1, 1)  = grhs(1, 1)
            end if

            if (kswx/=2 .and. kswx/=5 .and. kswy/=2 .and. kswy/=3) then
                usol(k, 1) = grhs(k, 1)
            end if

            if (kswx/=2 .and. kswx/=3 .and. kswy/=2 .and. kswy/=5) then
                usol(1, l) = grhs(1, l)
            end if

            if (kswx/=2 .and. kswx/=5 .and. kswy/=2 .and. kswy/=5) then
                usol(k, l) = grhs(k, l)
            end if

            i1 = 1
            !
            ! set switches for periodic or non-periodic boundaries
            !
            mp = 1
            np = 1

            if (kswx == 1) then
                mp = 0
            end if

            if (kswy == 1) then
                np = 0
            end if
            !
            !     set dlx, dly and size of block tri-diagonal system generated
            !     in nint, mint
            !
            dlx = (bit - ait)/m
            mit = k - 1

            select case (kswx)
                case (2)
                    mit = k - 2
                case (4)
                    mit = k
            end select

            dly = (dit - cit)/n
            nit = l - 1

            select case (kswy)
                case (2)
                    nit = l - 2
                case (4)
                    nit = l
            end select

            tdlx3 = TWO * (dlx**3)
            dlx4 = dlx**4
            tdly3 = TWO * (dly**3)
            dly4 = dly**4
            !
            ! set subscript limits for portion of array to input to blktri
            !
            is = 1
            js = 1

            if (kswx==2 .or. kswx==3) then
                is = 2
            end if

            if (kswy==2 .or. kswy==3) then
                js = 2
            end if

            ns = nit + js - 1
            ms = mit + is - 1
            !
            !     set x - direction
            !
            do i = 1, mit
                xi = ait + real(is + i - 2, kind=wp)*dlx
                call cofx(xi, ai, bi, ci)
                axi = (ai/dlx - HALF*bi)/dlx
                bxi = (-TWO*ai/dlx**2) + ci
                cxi = (ai/dlx + HALF*bi)/dlx
                am(i) = axi
                bm(i) = bxi
                cm(i) = cxi
            end do
            !
            !     set y direction
            !
            do j = 1, nit
                yj = cit + real(js + j - 2, kind=wp)*dly
                call cofy(yj, dj, ej, fj)
                dyj = (dj/dly - HALF*ej)/dly
                eyj = (-TWO*dj/dly**2) + fj
                fyj = (dj/dly + HALF*ej)/dly
                an(j) = dyj
                bn(j) = eyj
                cn(j) = fyj
            end do
            !
            !     adjust edges in x direction unless periodic
            !
            ax1 = am(1)
            cxm = cm(mit)
            select case (kswx)
                case (2)
                    !
                    !     dirichlet-dirichlet in x direction
                    !
                    am(1) = ZERO
                    cm(mit) = ZERO
                case (5)
                    !
                    !     mixed-dirichlet in x direction
                    !
                    am(1) = ZERO
                    bm(1) = bm(1) + TWO*alpha*dlx*ax1
                    cm(1) = cm(1) + ax1
                    cm(mit) = ZERO
                case (3)
                    !
                    !     dirichlet-mixed in x direction
                    !
                    am(1) = ZERO
                    am(mit) = am(mit) + cxm
                    bm(mit) = bm(mit) - TWO*beta*dlx*cxm
                    cm(mit) = ZERO
                !
                !     mixed - mixed in x direction
                !
                case (4)
                    am(1) = ZERO
                    bm(1) = bm(1) + TWO*dlx*alpha*ax1
                    cm(1) = cm(1) + ax1
                    am(mit) = am(mit) + cxm
                    bm(mit) = bm(mit) - TWO*dlx*beta*cxm
                    cm(mit) = ZERO
            end select
            !
            !     adjust in y direction unless periodic
            !
            dy1 = an(1)
            fyn = cn(nit)

            select case (kswy)
                case (2)
                    !
                    !     dirichlet-dirichlet in y direction
                    !
                    an(1) = ZERO
                    cn(nit) = ZERO
                case (5)
                    !
                    !     mixed-dirichlet in y direction
                    !
                    an(1) = ZERO
                    bn(1) = bn(1) + TWO*dly*gama*dy1
                    cn(1) = cn(1) + dy1
                    cn(nit) = ZERO
                case (3)
                    !
                    !     dirichlet-mixed in y direction
                    !
                    an(1) = ZERO
                    an(nit) = an(nit) + fyn
                    bn(nit) = bn(nit) - TWO*dly*xnu*fyn
                    cn(nit) = ZERO
                case (4)
                    !
                    !     mixed - mixed direction in y direction
                    !
                    an(1) = ZERO
                    bn(1) = bn(1) + TWO*dly*gama*dy1
                    cn(1) = cn(1) + dy1
                    an(nit) = an(nit) + fyn
                    bn(nit) = bn(nit) - TWO * dly*xnu*fyn
                    cn(nit) = ZERO
            end select

            if (kswx /= 1) then
                !
                !     adjust usol along x edge
                !
                if (kswx==2 .or. kswx==3) then
                    if (kswx==2 .or. kswx==5) then
                        usol(is, js:ns) = usol(is, js:ns) - ax1*usol(1, js:ns)
                        usol(ms, js:ns) = usol(ms, js:ns) - cxm*usol(k, js:ns)
                    else
                        usol(is, js:ns) = usol(is, js:ns) - ax1*usol(1, js:ns)
                        usol(ms, js:ns) = usol(ms, js:ns) - TWO * dlx*cxm*bdb(js:ns)
                    end if
                else
                    if (kswx==2 .or. kswx==5) then
                        usol(is, js:ns) = usol(is, js:ns) + TWO * dlx*ax1*bda(js:ns)
                        usol(ms, js:ns) = usol(ms, js:ns) - cxm*usol(k, js:ns)
                    else
                        usol(is, js:ns) = usol(is, js:ns) + TWO * dlx*ax1*bda(js:ns)
                        usol(ms, js:ns) = usol(ms, js:ns) - TWO * dlx*cxm*bdb(js:ns)
                    end if
                end if
            end if

            if (kswy /= 1) then
                !
                !     adjust usol along y edge
                !
                if (kswy==2 .or. kswy==3) then
                    if (kswy==2 .or. kswy==5) then
                        usol(is:ms, js) = usol(is:ms, js) - dy1*usol(is:ms, 1)
                        usol(is:ms, ns) = usol(is:ms, ns) - fyn*usol(is:ms, l)
                    else
                        usol(is:ms, js) = usol(is:ms, js) - dy1*usol(is:ms, 1)
                        usol(is:ms, ns) = usol(is:ms, ns) - TWO * dly*fyn*bdd(is:ms)
                    end if
                else
                    if (kswy==2 .or. kswy==5) then
                        usol(is:ms, js) = usol(is:ms, js) + TWO * dly*dy1*bdc(is:ms)
                        usol(is:ms, ns) = usol(is:ms, ns) - fyn*usol(is:ms, l)
                    else
                        usol(is:ms, js) = usol(is:ms, js) + TWO * dly*dy1*bdc(is:ms)
                        usol(is:ms, ns) = usol(is:ms, ns) - TWO * dly*fyn*bdd(is:ms)
                    end if
                end if
            end if
            !
            ! save adjusted edges in grhs if iorder=4
            !
            if (iorder == 4) then
                grhs(is, js:ns) = usol(is, js:ns)
                grhs(ms, js:ns) = usol(ms, js:ns)
                grhs(is:ms, js) = usol(is:ms, js)
                grhs(is:ms, ns) = usol(is:ms, ns)
            end if

            ! Initialize perturbation
            pertrb = ZERO
            !
            ! check if operator is singular
            !
            call self%is_PDE_singular(mbdcnd, nbdcnd, alpha, beta, &
                gama, xnu, cofx, cofy, singular)
            !
            ! compute non-zero eigenvector in null space of transpose
            !     if singular
            !
            if (singular) then
                call self%septri(mit, am, bm, cm, dm, um, zm)
            end if

            if (singular) then
                call self%septri(nit, an, bn, cn, dn, un, zn)
            end if
            !
            ! make initialization call to blktrii
            !
            if (intl == 0) then
                call self%blktri_aux%blktrii(intl, np, nit, an, bn, cn, mp, mit, am, bm, cm, &
                    idmn, usol(is:, js:), ierror, w, wc)

                ! Check error flag
                if (ierror /= 0) then
                    return
                end if
            end if
            !
            !     adjust right hand side if necessary
            !
            if (singular) then
                call self%seport(usol, zn, zm, pertrb)
            end if
            !
            !     compute solution
            !
            call self%blktri_aux%blktrii(i1, np, nit, an, bn, cn, mp, mit, am, bm, cm, idmn, &
                usol(is:, js:), ierror, w, wc)

            if (ierror /= 0) then
                return
            end if
            !
            !     set periodic boundaries if necessary
            !
            if (kswx == 1) then
                usol(k, :l) = usol(1, :l)
            end if

            if (kswy == 1) then
                usol(:k, l) = usol(:k, 1)
            end if
            !
            !     minimize solution with respect to weighted least squares
            !     norm if operator is singular
            !
            if (singular) then
                call self%sepmin(usol, zn, zm, prtrb)
            end if
            !
            !     return if deferred corrections and a fourth order solution are
            !     not flagged
            !
            if (iorder == 2) return
            !
            !     compute new right hand side for fourth order solution
            !
            call self%defer(cofx, cofy, idmn, usol, grhs)

            if (singular) call self%seport(usol, zn, zm, pertrb)
            !
            !     compute fourth order solution
            !
            call self%blktri_aux%blktrii(i1, np, nit, an, bn, cn, mp, mit, am, bm, cm, idmn, &
                usol(is:, js:), ierror, w, wc)

            if (ierror /= 0) return
            !
            !     set periodic boundaries if necessary
            !
            if (kswx == 1) then
                usol(k, :l) = usol(1, :l)
            end if
            if (kswy == 1) then
                usol(:k, l) = usol(:k, 1)
            end if
            !
            !     minimize solution with respect to weighted least squares
            !     norm if operator is singular
            !
            if (singular) call self%sepmin(usol, zn, zm, prtrb)

        end associate

    end subroutine spelip

    subroutine check_input_arguments(intl, iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, &
        cofx, cofy, idmn, ierror)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip), intent(in)    :: intl
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
        !--------------------------------------------------------------
        ! Dummy procedure arguments
        !--------------------------------------------------------------
        procedure(get_coefficients) :: cofx
        procedure(get_coefficients) :: cofy
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer(ip) :: i, j
        real(wp)    :: dlx, dly, xi, ai, bi, ci, yj, dj, ej, fj
        !-----------------------------------------------

        !     check definition of solution region
        !
        if (a>=b .or. c>=d) then
            ierror = 1
            return
        end if
        !
        !     check boundary condition arguments
        !
        if (mbdcnd<0 .or. mbdcnd>4) then
            ierror = 2
            return
        end if
        if (nbdcnd<0 .or. nbdcnd>4) then
            ierror = 3
            return
        end if
        !
        !     check first dimension in calling routine
        !
        if (idmn < 7) then
            ierror = 5
            return
        end if
        !
        !     check m, n
        !
        if (m>idmn - 1 .or. m<6) then
            ierror = 6
            return
        end if
        if (n < 5) then
            ierror = 7
            return
        end if
        !
        !     check iorder
        !
        if (iorder/=2 .and. iorder/=4) then
            ierror = 8
            return
        end if
        !
        !     check intl
        !
        if (intl/=0 .and. intl/=1) then
            ierror = 9
            return
        end if
        !
        !     check that equation is elliptic (only on initial call)
        !
        if (intl == 0) then
            dlx = (b - a)/m
            dly = (d - c)/n
            outer_loop: do i = 2, m
                xi = a + real(i - 1)*dlx

                call cofx(xi, ai, bi, ci)

                inner_loop: do j = 2, n
                    yj = c + real(j - 1)*dly
                    call cofy(yj, dj, ej, fj)

                    if (ai*dj > ZERO) cycle inner_loop

                    ierror = 10
                    return
                end do inner_loop
            end do outer_loop
        end if
        !
        !     no error found
        !
        ierror = 0

    end subroutine check_input_arguments

    subroutine is_PDE_singular(self, mbdcnd, nbdcnd, alpha, beta, gama, xnu, cofx, cofy, singlr)
        !
        ! Purpose:
        !
        ! Checks if the PDE that sepeli must solve is a singular operator
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(SepeliAux), intent(inout) :: self
        integer(ip),   intent(in)     :: mbdcnd
        integer(ip),   intent(in)     :: nbdcnd
        real(wp),      intent(in)     :: alpha
        real(wp),      intent(in)     :: beta
        real(wp),      intent(in)     :: gama
        real(wp),      intent(in)     :: xnu
        logical ,       intent(out)    :: singlr
        procedure(get_coefficients)    :: cofx
        procedure(get_coefficients)    :: cofy
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: i, j
        real(wp)    :: xi, ai, bi, ci, yj, dj, ej, fj
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

            ! Initialize flag
            singlr = .false.
            !
            !     check if the boundary conditions are
            !     entirely periodic and/or mixed
            !
            if(mbdcnd/=0.and.mbdcnd/=3.or.nbdcnd/=0.and.nbdcnd/=3) then
                return
            end if
            !
            !     check that mixed conditions are pure neuman
            !
            if (mbdcnd == 3) then
                if (alpha/=ZERO .or. beta/=ZERO) return
            end if

            if (nbdcnd == 3) then
                if (gama/=ZERO .or. xnu/=ZERO) return
            end if
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

            do j = js, ns
                yj = cit + real(j - 1, kind=wp)*dly
                call cofy(yj, dj, ej, fj)
                if (fj == ZERO) cycle
                return
            end do
            !
            !     the operator must be singular if this point is reached
            !
            singlr = .true.

        end associate

    end subroutine is_PDE_singular

    subroutine defer(self, cofx, cofy, idmn, usol, grhs)
        !
        ! Purpose:
        !
        !     this subroutine first approximates the truncation error given by
        !     trun1(x, y)=dlx**2*tx+dly**2*ty where
        !     tx=afun(x)*uxxxx/12+bfun(x)*uxxx/6 on the interior and
        !     at the boundaries if periodic(here uxxx, uxxxx are the third
        !     and fourth partial derivatives of u with respect to x).
        !     tx is of the form afun(x)/3 * (uxxxx/4 +uxxx/dlx)
        !     at x=a or x=b if the boundary condition there is mixed.
        !     tx=0.0 along specified boundaries.  ty has symmetric form
        !     in y with x, afun(x), bfun(x) replaced by y, dfun(y), efun(y).
        !     the second order solution in usol is used to approximate
        !     (via second order finite differencing) the truncation error
        !     and the result is added to the right hand side in grhs
        !     and then transferred to usol to be used as a new right
        !     hand side when calling blktri for a fourth order solution.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(SepeliAux), intent(inout) :: self
        integer(ip),   intent(in)     :: idmn
        real(wp),      intent(inout) :: usol(:,:)
        real(wp),      intent(inout) :: grhs(:,:)
        procedure(get_coefficients)    :: cofx
        procedure(get_coefficients)    :: cofy
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: j, i
        real(wp)    :: yj, dj, ej, fj, xi, ai, bi, ci
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

            !
            !     compute truncation error approximation over the entire mesh
            !
            do j = js, ns
                yj = cit + real(j - 1, kind=wp)*dly
                call cofy(yj, dj, ej, fj)
                do i = is, ms
                    xi = ait + real(i - 1, kind=wp)*dlx
                    call cofx(xi, ai, bi, ci)
                    !
                    !     compute partial derivative approximations at (xi, yj)
                    !
                    call self%sepdx(usol, i, j, uxxx, uxxxx)
                    call self%sepdy(usol, idmn, i, j, uyyy, uyyyy)
                    tx = ai*uxxxx/12 + bi*uxxx/6
                    ty = dj*uyyyy/12 + ej*uyyy/6
                    !
                    ! reset form of truncation if at boundary which is non-periodic
                    !
                    if (kswx/=1 .and. (i==1 .or. i==k)) then
                        tx = (ai/3) * (uxxxx/4 + uxxx/dlx)
                    end if

                    if (kswy/=1 .and. (j==1 .or. j==l)) then
                        ty = (dj/3) * (uyyyy/4 + uyyy/dly)
                    end if

                    grhs(i, j) = grhs(i, j) + (dlx**2)*tx + (dly**2)*ty
                end do
            end do
            !
            ! reset the right hand side in usol
            !
            usol(is:ms, js:ns) = grhs(is:ms, js:ns)

        end associate

    end subroutine defer

end module module_sepeli
!
! REVISION HISTORY
!
! September 1973    Version 1
! April     1976    Version 2
! January   1978    Version 3
! December  1979    Version 3.1
! February  1985    Documentation upgrade
! November  1988    Version 3.2, FORTRAN 77 changes
! June      2004    Version 5.0, fortran 90 changes
! May       2016    Fortran 2008 changes
!
