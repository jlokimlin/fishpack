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
module staggered_real_linear_systems_solver

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_StaggeredCyclicReductionUtility, only: &
        StaggeredCyclicReductionUtility

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: poistg

    ! Parameters confined to the module
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: LOG_TWO = log(TWO)

contains

    ! SUBROUTINE poistg(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
    !
    !
    ! DIMENSION OF           a(m),  b(m),  c(m),  y(idimy, n)
    ! ARGUMENTS
    !
    ! PURPOSE                Solves the linear system of equations
    !                        for unknown x values, where i=1, 2, ..., m
    !                        and j=1, 2, ..., n
    !
    !                        a(i)*x(i-1, j) + b(i)*x(i, j) + c(i)*x(i+1, j)
    !                        + x(i, j-1) - TWO *x(i, j) + x(i, j+1)
    !                        = y(i, j)
    !
    !                        The indices i+1 and i-1 are evaluated modulo m,
    !                        i.e. x(0, j) = x(m, j) and x(m+1, j) = x(1, j), and
    !                        x(i, 0) may be equal to x(i, 1) or -x(i, 1), and
    !                        x(i, n+1) may be equal to x(i, n) or -x(i, n),
    !                        depending on an input parameter.
    !
    ! USAGE                  call poistg(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
    !
    ! ARGUMENTS
    ! ON INPUT
    !
    !                        nperod
    !                          Indicates values which x(i, 0) and x(i, n+1)
    !                          are assumed to have.
    !                          = 1 if x(i, 0) = -x(i, 1) and x(i, n+1) = -x(i, n
    !                          = 2 if x(i, 0) = -x(i, 1) and x(i, n+1) =  x(i, n
    !                          = 3 if x(i, 0) =  x(i, 1) and x(i, n+1) =  x(i, n
    !                          = 4 if x(i, 0) =  x(i, 1) and x(i, n+1) = -x(i, n
    !
    !                        n
    !                          the number of unknowns in the j-direction.
    !                          n must be greater than 2.
    !
    !                        mperod
    !                          = 0 if a(1) and c(m) are not zero
    !                          = 1 if a(1) = c(m) = 0
    !
    !                        m
    !                          the number of unknowns in the i-direction.
    !                          m must be greater than 2.
    !
    !                        a, b, c
    !                          One-dimensional arrays of length m that
    !                          specify the coefficients in the linear
    !                          equations given above. If mperod = 0 the
    !                          array elements must not depend on index i,
    !                          but must be constant. Specifically, the
    !                          subroutine checks the following condition
    !                            a(i) = c(1)
    !                            b(i) = b(1)
    !                            c(i) = c(1)
    !                          for i = 1, 2, ..., m.
    !
    !                        idimy
    !                          The row (or first) dimension of the two-
    !                          dimensional array y as it appears in the
    !                          program calling poistg. This parameter is
    !                          used to specify the variable dimension of y.
    !                          idimy must be at least m.
    !
    !                        y
    !                          A two-dimensional array that specifies the
    !                          values of the right side of the linear system
    !                          of equations given above.
    !                          y must be dimensioned at least m x n.
    !
    ! ON OUTPUT
    !
    !                        Y
    !                          Contains the solution x.
    !
    !                        ierror
    !                          An error flag that indicates invalid input
    !                          parameters.  except for number zero, a
    !                          solution is not attempted.
    !                          = 0  no error
    !                          = 1  if m <= 2
    !                          = 2  if n <= 2
    !                          = 3  idimy < m
    !                          = 4  if nperod < 1 or nperod > 4
    !                          = 5  if mperod < 0 or mperod > 1
    !                          = 6  if mperod = 0 and a(i) /= c(1)
    !                               or b(i) /= b(1) or c(i) /= c(1)
    !                               for some i = 1, 2, ..., m.
    !                          = 7  if mperod == 1 and (a(1)/= 0 or c(m)/= 0)
    !                          = 20 If the dynamic allocation of real and
    !                               complex workspace required for solution
    !                               fails (for example if n, m are too large
    !                               for your computer)
    !
    !                          Since this is the only means of indicating a
    !                          possibly incorrect call to pois3d, the user
    !                          should test ierror after the call.
    !
    !
    ! HISTORY                Written by Roland Sweet at NCAR in the late
    !                        1970's.  Released on NCAR'S public software
    !                        libraries in January, 1980.
    !                        Revised in June 2004 by John Adams using
    !                        Fortran 90 dynamically allocated workspace.
    !
    !
    ! ALGORITHM              This subroutine is an implementation of the
    !                        algorithm presented in the reference below.
    !
    ! TIMING                 For large m and n, the execution time is
    !                        roughly proportional to m*n*log2(n).
    !
    ! ACCURACY               To measure the accuracy of the algorithm a
    !                        uniform random number generator was used to
    !                        create a solution array x for the system given
    !                        in the 'purpose' section above, with
    !
    !                          a(i) = c(i) = -HALF * b(i) = 1,    i=1, 2, ..., m
    !                        and, when mperod = 1
    !                          a(1) = c(m) = 0
    !                          b(1) = b(m) =-1.
    !
    !                        The solution x was substituted into the given
    !                        system and, using double precision, a right sid
    !                        y was computed.  using this array y subroutine
    !                        poistg was called to produce an approximate
    !                        solution z.  then the relative error, defined a
    !
    !                         e = max(abs(z(i, j)-x(i, j)))/max(abs(x(i, j)))
    !
    !                        where the two maxima are taken over i=1, 2, ..., m
    !                        and j=1, 2, ..., n, was computed.  values of e are
    !                        given in the table below for some typical value
    !                        of m and n.
    !
    !                        m (=n)    mperod    nperod      e
    !                        ------    ------    ------    ------
    !
    !                          31        0-1       1-4     9.e-13
    !                          31        1         1       4.e-13
    !                          31        1         3       3.e-13
    !                          32        0-1       1-4     3.e-12
    !                          32        1         1       3.e-13
    !                          32        1         3       1.e-13
    !                          33        0-1       1-4     1.e-12
    !                          33        1         1       4.e-13
    !                          33        1         3       1.e-13
    !                          63        0-1       1-4     3.e-12
    !                          63        1         1       1.e-12
    !                          63        1         3       2.e-13
    !                          64        0-1       1-4     4.e-12
    !                          64        1         1       1.e-12
    !                          64        1         3       6.e-13
    !                          65        0-1       1-4     2.e-13
    !                          65        1         1       1.e-11
    !                          65        1         3       4.e-13
    !
    ! REFERENCES             Schumann, U. and R. Sweet, "A direct method
    !                        for the solution of Poisson's equation with
    !                        Neumann boundary conditions on a staggered
    !                        grid of arbitrary size, " J. Comp. Phys.
    !                        20 (1976), pp. 171-182.
    !
    subroutine poistg(nperod, n, mperod, m, a, b, c, idimy, y, ierror)

        ! Dummy arguments
        integer(ip), intent(in)     :: nperod
        integer(ip), intent(in)     :: n
        integer(ip), intent(in)     :: mperod
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: idimy
        integer(ip), intent(out)    :: ierror
        real(wp),    intent(in)     :: a(:)
        real(wp),    intent(in)     :: b(:)
        real(wp),    intent(in)     :: c(:)
        real(wp),    intent(inout)  :: y(:,:)

        ! Local variables
        integer(ip)                           :: irwk, icwk
        type(FishpackWorkspace)               :: workspace
        type(StaggeredCyclicReductionUtility) :: util

        ! Compute workspace dimensions
        irwk = m * (9 + int(log(real(n, kind=wp))/LOG_TWO,kind=ip)) + 4 * n
        icwk = 0

        ! Allocate memory
        call workspace%create(irwk, icwk)

        ! Solve system
        associate( rew => workspace%real_workspace )
            call util%poistg_lower_routine(nperod, n, mperod, m, a, b, c, idimy, y, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine poistg

end module staggered_real_linear_systems_solver
