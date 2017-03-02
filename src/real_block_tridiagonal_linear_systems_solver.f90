!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                           Fishpack                            *
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
module real_block_tridiagonal_linear_systems_solver

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_GeneralizedCyclicReductionUtility, only: &
        GeneralizedCyclicReductionUtility

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: blktri

contains

    !
    ! SUBROUTINE blktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm,
    !            idimy, y, ierror, workspace)
    !
    !
    !
    ! DIMENSION OF           an(n), bn(n), cn(n), am(m), bm(m), cm(m), y(idimy, n)
    ! ARGUMENTS
    !
    ! USAGE                  call blktri(iflg, np, n, an, bn, cn, mp, m, &
    !                                   am, bm, cm, idimy, y, ierror, workspace)
    !
    ! PURPOSE                blktri solves a system of linear equations
    !                        of the form
    !
    !                        an(j)*x(i, j-1) + am(i)*x(i-1, j) +
    !                        (bn(j)+bm(i))*x(i, j) + cn(j)*x(i, j+1) +
    !                        cm(i)*x(i+1, j) = y(i, j)
    !
    !                        for i = 1, 2, ..., m  and  j = 1, 2, ..., n.
    !
    !                        i+1 and i-1 are evaluated modulo m and
    !                        j+1 and j-1 modulo n, i.e.,
    !
    !                        x(i, 0) = x(i, n),  x(i, n+1) = x(i, 1),
    !                        x(0, j) = x(m, j),  x(m+1, j) = x(1, j).
    !
    !                        These equations usually result from the
    !                        discretization of separable elliptic
    !                        equations. Boundary conditions may be
    !                        dirichlet, neumann, or periodic.
    !
    ! ARGUMENTS
    !
    ! ON INPUT               iflg
    !
    !                          = 0  Unitialization only.
    !                               certain quantities that depend on np,
    !                               n, an, bn, and cn are computed and
    !                               stored in derived data type w (see
    !                               description of w below)
    !
    !                          = 1  The quantities that were computed
    !                               in the initialization are used
    !                               to obtain the solution x(i, j).
    !
    !                               note:
    !                               A call with iflg=0 takes
    !                               approximately one half the time
    !                               as a call with iflg = 1.
    !                               however, the initialization does
    !                               not have to be repeated unless np,
    !                               n, an, bn, or cn change.
    !
    !                        np
    !                          = 0  If an(1) and cn(n) are not zero,
    !                               which corresponds to periodic
    !                               bounary conditions.
    !
    !                          = 1  If an(1) and cn(n) are zero.
    !
    !                        n
    !                          The number of unknowns in the j-direction.
    !                          n must be greater than 4.
    !                          The operation count is proportional to
    !                          m*n*log2(n), hence n should be selected
    !                          less than or equal to m.
    !
    !                        an, bn, cn
    !                          One-dimensional arrays of length n
    !                          that specify the coefficients in the
    !                          linear equations given above.
    !
    !                        mp
    !                          = 0  If am(1) and cm(m) are not zero,
    !                               which corresponds to periodic
    !                               boundary conditions.
    !
    !                          = 1  If am(1) = cm(m) = 0  .
    !
    !                        m
    !                          The number of unknowns in the i-direction.
    !                           m must be greater than 4.
    !
    !                        am, bm, cm
    !                          One-dimensional arrays of length m that
    !                          specify the coefficients in the linear
    !                          equations given above.
    !
    !                        idimy
    !                          The row (or first) dimension of the
    !                          two-dimensional array y as it appears
    !                          in the program calling blktri.
    !                          This parameter is used to specify the
    !                          variable dimension of y.
    !                          idimy must be at least m.
    !
    !                        y
    !                          A two-dimensional array that specifies
    !                          the values of the right side of the linear
    !                          system of equations given above.
    !                          y must be dimensioned at least m*n.
    !
    !                        workspace
    !                          A FishpackWorkspace derived data type
    !                          that must be declared by the user which is used
    !                          internally in blktri to dynamically allocate real
    !                          and complex workspace arrays used in calls to lower routines.
    !                          An error flag (ierror = 20) is set if the required
    !                          workspace allocation fails (for example if n, m
    !                          are too large). Real and complex values are set in
    !                          the components of workspace on a initial (iflg=0)
    !                          call to blktri.  These must be preserved on
    !                          non-initial calls (iflg=1) to blktri.
    !                          This eliminates redundant calculations
    !                          and saves compute time.
    !
    !               ****       IMPORTANT! The user program calling blktri should
    !                          include the statement:
    !
    !                              call workspace%destroy()
    !
    !                          after the final approximation is generated by
    !                          blktri. This will deallocate the real and complex
    !                          array components of workspace. Failure to include this
    !                          statement could result in serious memory leakage.
    !
    !
    ! ARGUMENTS
    !
    ! ON OUTPUT              y
    !                          Contains the solution x.
    !
    !                        ierror
    !                          An error flag that indicates invalid
    !                          input parameters.  except for number zer0,
    !                          a solution is not attempted.
    !
    !                        = 0  no error.
    !                        = 1  m < than 5
    !                        = 2  n < than 5
    !                        = 3  idimy < m.
    !                        = 4  blktri failed while computing results
    !                             that depend on the coefficient arrays
    !                             an, bn, cn. Check these arrays.
    !                        = 5  an(j)*cn(j-1) is less than 0 for some j.
    !
    !                             Possible reasons for this condition are
    !                             1. The arrays an and cn are not correct
    !                             2. Too large a grid spacing was used
    !                                in the discretization of the elliptic
    !                                equation.
    !                             3. The linear equations resulted from a
    !                                partial differential equation which
    !                                was not elliptic.
    !
    !                        = 20 If the dynamic allocation of real and
    !                             complex workspace in the derived type
    !                             (FishpackWorkspace) variable W fails (e.g.,
    !                             if N, M are too large for the platform used)
    !
    !
    !                        workspace
    !                             The derived type(FishpackWorkspace) variable
    !                             contains real and complex array components that
    !                             must not be destroyed if blktri is called again with
    !                             iflg=1.
    !
    !
    ! SPECIAL CONDITIONS     The algorithm may fail if abs(bm(i)+bn(j))
    !                        is less than abs(am(i))+abs(an(j))+
    !                        abs(cm(i))+abs(cn(j))
    !                        for some i and j. the algorithm will also
    !                        fail if an(j)*cn(j-1) is less than zero for
    !                        some j.
    !                        see the description of the output parameter
    !                        ierror.
    !
    !
    ! HISTORY                * Written by Paul Swarztrauber at NCAR in the
    !                          early 1970's.
    !                        * Rewritten and released in libraries in January 1980.
    !                        * Revised in June 2004 using Fortan 90 dynamically
    !                          allocated workspace and derived data types to
    !                          eliminate mixed mode conflicts in the earlier versions.
    !                        * Revised in April 2016 to implement features of
    !                          Fortran 2008+
    !
    ! ALGORITHM              Generalized cyclic reduction
    !
    ! PORTABILITY            Approximate machine accuracy is obtained
    !                        using EPS which is set by the intrinsic epsilon function
    !
    ! REFERENCES             Swarztrauber, P. and R. Sweet, 'Efficient
    !                        fortran subprograms for the solution of
    !                        elliptic equations'
    !                        NCAR TN/IA-109, July, 1975, 138 pp.
    !
    !                        Swarztrauber P. N., A direct method for
    !                        the discrete solution of separable
    !                        elliptic equations, SIAM
    !                        J. Numer. Anal., 11(1974) pp. 1136-1150.
    !
    subroutine blktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, &
        idimy, y, ierror, workspace)

        ! Dummy arguments
        integer(ip), intent(in)    :: iflg
        integer(ip), intent(in)    :: np
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: mp
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: idimy
        integer(ip), intent(out)   :: ierror
        real(wp),    intent(in)    :: an(:)
        real(wp),    intent(in)    :: bn(:)
        real(wp),    intent(in)    :: cn(:)
        real(wp),    intent(in)    :: am(:)
        real(wp),    intent(in)    :: bm(:)
        real(wp),    intent(in)    :: cm(:)
        real(wp),    intent(inout) :: y(:,:)
        class(FishpackWorkspace), intent(inout) :: workspace

        ! Local variables
        type(GeneralizedCyclicReductionUtility), save :: blktri_aux

        ! Check input arguments
        call blktri_aux%check_input_arguments(n, m, idimy, ierror)

        ! Check error flag
        if (ierror /= 0) return

        ! Allocate memory on first call only
        if (iflg == 0) call initialize_blktri_workspace(n, m, workspace)

        ! Solve system
        associate( &
            rew => workspace%real_workspace, &
            cxw => workspace%complex_workspace &
            )
            call blktri_aux%blktrii(iflg, np, n, an, bn, cn, &
                mp, m, am, bm, cm, idimy, y, ierror, rew, cxw)
        end associate

    end subroutine blktri

    subroutine initialize_blktri_workspace(n, m, workspace)

        ! Dummy arguments
        integer(ip),              intent(in)  :: n
        integer(ip),              intent(in)  :: m
        class(FishpackWorkspace), intent(out) :: workspace

        ! Local variables
        integer(ip) :: irwk, icwk

        ! compute real and complex required workspace sizes
        call workspace%compute_blktri_workspace_lengths(n, m, irwk, icwk)

        ! Allocate memory
        call workspace%create(irwk, icwk)

    end subroutine initialize_blktri_workspace

end module real_block_tridiagonal_linear_systems_solver
