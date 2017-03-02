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
submodule(three_dimensional_solvers) general_linear_systems_solver_3d

contains

    ! SUBROUTINE pois3d(lperod, l, c1, mperod, m, c2, nperod, n, &
    !            a, b, c, ldimf, mdimf, f, ierror)
    !
    !
    ! DIMENSION OF           a(n), b(n), c(n), f(ldimf, mdimf, n)
    ! ARGUMENTS
    !
    ! PURPOSE                Solves the linear system of equations
    !                        for unknown x values, where i=1, 2, ..., l,
    !                        j=1, 2, ..., m, and k=1, 2, ..., n
    !
    !                        c1*(x(i-1, j, k) -TWO *x(i, j, k) +x(i+1, j, k)) +
    !                        c2*(x(i, j-1, k) -TWO *x(i, j, k) +x(i, j+1, k)) +
    !                        a(k)*x(i, j, k-1) +b(k)*x(i, j, k)+ c(k)*x(i, j, k+1)
    !                        = f(i, j, k)
    !
    !                        the indices k-1 and k+1 are evaluated modulo n,
    !                        i.e. x(i, j, 0)=x(i, j, n) and x(i, j, n+1)=x(i, j, 1).
    !                        the unknowns
    !                        x(0, j, k), x(l+1, j, k), x(i, 0, k), and x(i, m+1, k)
    !                        are assumed to take on certain prescribed
    !                        values described below.
    !
    ! USAGE                  call pois3d (lperod, l, c1, mperod, m, c2, nperod,
    !                        n, a, b, c, ldimf, mdimf, f, ierror)
    !
    ! ARGUMENTS
    !
    ! ON INPUT
    !                        lperod
    !                          indicates the values that x(0, j, k) and
    !                          x(l+1, j, k) are assumed to have.
    !                          = 0  x(0, j, k)=x(l, j, k), x(l+1, j, k)=x(1, j, k)
    !                          = 1  x(0, j, k) = 0,      x(l+1, j, k) = 0
    !                          = 2  x(0, j, k)=0,        x(l+1, j, k)=x(l-1, j, k)
    !                          = 3  x(0, j, k)=x(2, j, k), x(l+1, j, k)=x(l-1, j, k)
    !                          = 4  x(0, j, k)=x(2, j, k), x(l+1, j, k) = 0.
    !
    !                        l
    !                          the number of unknowns in the i-direction.
    !                          l must be at least 3.
    !
    !                        c1
    !                          real constant in the above linear system
    !                          of equations to be solved.
    !
    !                        mperod
    !                          indicates the values that x(i, 0, k) and
    !                          x(i, m+1, k) are assumed to have.
    !                          = 0  x(i, 0, k)=x(i, m, k), x(i, m+1, k)=x(i, 1, k)
    !                          = 1  x(i, 0, k)=0,        x(i, m+1, k)=0
    !                          = 2  x(i, 0, k)=0,        x(i, m+1, k)=x(i, m-1, k)
    !                          = 3  x(i, 0, k)=x(i, 2, k)  x(i, m+1, k)=x(i, m-1, k)
    !                          = 4  x(i, 0, k)=x(i, 2, k)  x(i, m+1, k)=0
    !
    !                        m
    !                          the number of unknowns in the j-direction.
    !                          m must be at least 3.
    !
    !                        c2
    !                          real constant in the above linear system
    !                          of equations to be solved.
    !
    !                        nperod
    !                          = 0  if a(1) and c(n) are not zero.
    !                          = 1  if a(1) = c(n) = 0.
    !
    !                        n
    !                          the number of unknowns in the k-direction.
    !                          n must be at least 3.
    !
    !                        a, b, c
    !                          one-dimensional arrays of length n that
    !                          specify the coefficients in the linear
    !                          equations given above.
    !
    !                          if nperod = 0 the array elements must not
    !                          depend upon index k, but must be constant.
    !                          specifically, the subroutine checks the
    !                          following condition
    !                            a(k) = c(1)
    !                            c(k) = c(1)
    !                            b(k) = b(1)
    !                          for k=1, 2, ..., n.
    !
    !                        ldimf
    !                          The row (or first) dimension of the three-
    !                          dimensional array f as it appears in the
    !                          program calling pois3d.  this parameter is
    !                          used to specify the variable dimension
    !                          of f.  ldimf must be at least l.
    !
    !                        mdimf
    !                          The column (or second) dimension of the three
    !                          dimensional array f as it appears in the
    !                          program calling pois3d.  this parameter is
    !                          used to specify the variable dimension
    !                          of f.  mdimf must be at least m.
    !
    !                        f
    !                          A three-dimensional array that specifies the
    !                          values of the right side of the linear system
    !                          of equations given above.  f must be
    !                          dimensioned at least l x m x n.
    !
    ! ON OUTPUT
    !
    !                        f
    !                          Contains the solution x.
    !
    !                        ierror
    !                          An error flag that indicates invalid input
    !                          parameters.  except for number zero, a
    !                          solution is not attempted.
    !                          = 0  no error
    !                          = 1  if lperod < 0 or > 4
    !                          = 2  if l < 3
    !                          = 3  if mperod < 0 or > 4
    !                          = 4  if m < 3
    !                          = 5  if nperod < 0 or > 1
    !                          = 6  if n < 3
    !                          = 7  if ldimf < l
    !                          = 8  if mdimf < m
    !                          = 9  if a(k) /= c(1) or c(k) /= c(1)
    !                               or b(i) /=b(1) for some k=1, 2, ..., n.
    !                          = 10 if nperod = 1 and a(1) /= 0
    !                               or c(n) /= 0
    !                          = 20 If the dynamic allocation of real and
    !                               complex workspace required for solution
    !                               fails (for example if n, m are too large
    !                               for your computer)
    !
    !                          Since this is the only means of indicating a
    !                          possibly incorrect call to pois3d, the user
    !                          should test ierror after the call.
    !
    ! HISTORY                * Written by Roland Sweet at NCAR in the late
    !                          1970's.  released on NCAR's public software
    !                          libraries in January, 1980.
    !                        * Revised in June 2004 by John Adams using
    !                          Fortran 90 dynamically allocated workspace.
    !
    ! ALGORITHM              This subroutine solves three-dimensional block
    !                        tridiagonal linear systems arising from finite
    !                        difference approximations to three-dimensional
    !                        poisson equations using the FFT package
    !                        fftpack written by Paul Swarztrauber.
    !
    ! TIMING                 For large l, m and n, the operation count
    !                        is roughly proportional to
    !                          l*m*n*(log2(l)+log2(m)+5)
    !                        but also depends on input parameters lperod
    !                        and mperod.
    !
    ! ACCURACY               To measure the accuracy of the algorithm a
    !                        uniform random number generator was used to
    !                        create a solution array x for the system given
    !                        in the 'purpose' section with
    !                          a(k) = c(k) = -HALF *b(k) = 1,  k=1, 2, ..., n
    !                        and, when nperod = 1
    !                          a(1) = c(n) = 0
    !                          a(n) = c(1) = 2.
    !
    !                        The solution x was substituted into the given
    !                        system and, using double precision, a right
    !                        side y was computed.  using this array y
    !                        subroutine pois3d was called to produce an
    !                        approximate solution z.  Relative error
    !
    !                         = max(abs(z(i, j, k)-x(i, j, k)))/max(abs(x(i, j, k
    !
    !                        was computed, where the two maxima are taken
    !                        over i=1, 2, ..., l, j=1, 2, ..., m and k=1, 2, ..., n.
    !                        values of e are given in the table below for
    !                        some typical values of l, m and n.
    !
    !                        l(=m=n)   lperod    mperod       e
    !                        ------    ------    ------     ------
    !
    !                          16        0         0        1.e-13
    !                          15        1         1        4.e-13
    !                          17        3         3        2.e-13
    !                          32        0         0        2.e-13
    !                          31        1         1        2.e-12
    !                          33        3         3        7.e-13
    !
    module subroutine pois3d(lperod, l, c1, mperod, m, c2, nperod, n, a, b, c, &
        ldimf, mdimf, f, ierror)

        ! Dummy arguments
        integer(ip), intent(in)    :: lperod
        integer(ip), intent(in)    :: l
        integer(ip), intent(in)    :: mperod
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: nperod
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: ldimf
        integer(ip), intent(in)    :: mdimf
        integer(ip), intent(out)   :: ierror
        real(wp),    intent(in)    :: c1
        real(wp),    intent(in)    :: c2
        real(wp),    intent(inout) :: a(:)
        real(wp),    intent(inout) :: b(:)
        real(wp),    intent(inout) :: c(:)
        real(wp),    intent(inout) :: f(:,:,:)

        ! Local variables
        type(FishpackWorkspace) :: workspace
        type(SolverUtility3D)   :: util3d

        ! Check input arguments
        call util3d%check_input_arguments(lperod, l, mperod, m, nperod, n, &
            a, b, c, ldimf, mdimf, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Allocate memory
        call initialize_workspace(n, m, l, workspace)

        ! Solve system
        associate( &
            rew => workspace%real_workspace, &
            indx => workspace%workspace_indices &
            )
            call util3d%pois3dd(lperod, l, c1, mperod, m, c2, nperod, n, &
                a, b, c, ldimf, mdimf, f, ierror, rew, indx)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine pois3d

    subroutine initialize_workspace(n, m, l, workspace)

        ! Dummy arguments
        integer(ip),              intent(in)  :: n
        integer(ip),              intent(in)  :: m
        integer(ip),              intent(in)  :: l
        class(FishpackWorkspace), intent(out) :: workspace

        ! Local variables
        integer(ip)     :: irwk, icwk
        type(SolverUtility3D) :: aux

        ! Adjust workspace for pois3d
        irwk = 30 + l + m + (2 * n) + max(l, m, n) + 7 * ( (l + 1)/2 + (m + 1)/2)
        icwk = 0

        ! Allocate memory
        call workspace%create(irwk, icwk, aux%IIWK)

        ! Set workspace indices
        workspace%workspace_indices = aux%get_workspace_indices(l, m, n)

    end subroutine initialize_workspace

end submodule general_linear_systems_solver_3d
