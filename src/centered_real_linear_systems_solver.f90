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

module centered_real_linear_systems_solver

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_CenteredCyclicReductionUtility, only: &
        CenteredCyclicReductionUtility

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: genbun

contains

    ! SUBROUTINE genbun(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
    !
    !
    ! DIMENSION OF           a(m), b(m), c(m), y(idimy, n)
    ! ARGUMENTS
    !
    ! PURPOSE                The name of this package is a mnemonic for the
    !                        generalized buneman algorithm.
    !
    !                        It solves the real linear system of equations
    !
    !                        a(i)*x(i-1, j) + b(i)*x(i, j) + c(i)*x(i+1, j)
    !                        + x(i, j-1) - 2.0*x(i, j) + x(i, j+1) = y(i, j)
    !
    !                        for i = 1, 2, ..., m  and  j = 1, 2, ..., n.
    !
    !                        indices i+1 and i-1 are evaluated modulo m,
    !                        i.e., x(0, j) = x(m, j) and x(m+1, j) = x(1, j),
    !                        and x(i, 0) may equal 0, x(i, 2), or x(i, n),
    !                        and x(i, n+1) may equal 0, x(i, n-1), or x(i, 1)
    !                        depending on an input parameter.
    !
    ! USAGE                  call genbun(nperod, n, mperod, m, a, b, c, idimy, y, ierror)
    !
    ! ARGUMENTS
    !
    ! ON INPUT               nperod
    !
    !                          Indicates the values that x(i, 0) and
    !                          x(i, n+1) are assumed to have.
    !
    !                          = 0  if x(i, 0) = x(i, n) and x(i, n+1) =
    !                               x(i, 1).
    !                          = 1  if x(i, 0) = x(i, n+1) = 0  .
    !                          = 2  if x(i, 0) = 0 and x(i, n+1) = x(i, n-1).
    !                          = 3  if x(i, 0) = x(i, 2) and x(i, n+1) =
    !                               x(i, n-1).
    !                          = 4  if x(i, 0) = x(i, 2) and x(i, n+1) = 0.
    !
    !                        n
    !                          The number of unknowns in the j-direction.
    !                          n must be greater than 2.
    !
    !                        mperod
    !                          = 0 if a(1) and c(m) are not zero
    !                          = 1 if a(1) = c(m) = 0
    !
    !                        m
    !                          The number of unknowns in the i-direction.
    !                          n must be greater than 2.
    !
    !                        a, b, c
    !                          One-dimensional arrays of length m that
    !                          specify the coefficients in the linear
    !                          equations given above.  if mperod = 0
    !                          the array elements must not depend upon
    !                          the index i, but must be constant.
    !                          specifically, the subroutine checks the
    !                          following condition .
    !
    !                            a(i) = c(1)
    !                            c(i) = c(1)
    !                            b(i) = b(1)
    !
    !                          for i=1, 2, ..., m.
    !
    !                        idimy
    !                          The row (or first) dimension of the
    !                          two-dimensional array y as it appears
    !                          in the program calling genbun.
    !                          this parameter is used to specify the
    !                          variable dimension of y.
    !                          idimy must be at least m.
    !
    !                        y
    !                          A two-dimensional complex array that
    !                          specifies the values of the right side
    !                          of the linear system of equations given
    !                          above.
    !                          y must be dimensioned at least m*n.
    !
    !
    !  ON OUTPUT             y
    !
    !                          Contains the solution x.
    !
    !                        ierror
    !                          An error flag which indicates invalid
    !                          input parameters  except for number
    !                          zero, a solution is not attempted.
    !
    !                          = 0  no error
    !                          = 1  m <= 2
    !                          = 2  n <= 2
    !                          = 3  idimy <= m
    !                          = 4  nperod <= 0 or nperod > 4
    !                          = 5  mperod < 0 or mperod > 1
    !                          = 6  a(i) /= c(1) or c(i) /= c(1) or
    !                               b(i) /= b(1) for
    !                               some i=1, 2, ..., m.
    !                          = 7  a(1) /= 0 or c(m) /= 0 and
    !                                 mperod = 1
    !                          = 20 If the dynamic allocation of real and
    !                               complex workspace required for solution
    !                               fails (for example if n, m are too large
    !                               for your computer)
    !
    ! HISTORY                Written in 1979 by Roland Sweet of NCAR'S
    !                        scientific computing division. made available
    !                        on NCAR'S public libraries in January, 1980.
    !
    !                        Revised in June 2004 by John Adams using
    !                        Fortran 90 dynamically allocated workspace.
    !
    ! ALGORITHM              The linear system is solved by a cyclic
    !                        reduction algorithm described in the
    !                        reference.
    !
    ! REFERENCES             Sweet, R., "A cyclic reduction algorithm for
    !                        solving block tridiagonal systems of arbitrary
    !                        dimensions, " SIAM J. On Numer. Anal., 14 (1977)
    !                        PP. 706-720.
    !
    ! ACCURACY               This test was performed on a platform with
    !                        64 bit floating point arithmetic.
    !                        a uniform random number generator was used
    !                        to create a solution array x for the system
    !                        given in the 'purpose' description above
    !                        with
    !                          a(i) = c(i) = -HALF *b(i) = 1, i=1, 2, ..., m
    !
    !                        and, when mperod = 1
    !
    !                          a(1) = c(m) = 0
    !                          a(m) = c(1) = 2.
    !
    !                        The solution x was substituted into the
    !                        given system  and, using double precision
    !                        a right side y was computed.
    !                        using this array y, subroutine genbun
    !                        was called to produce approximate
    !                        solution z.  then relative error
    !                          e = max(abs(z(i, j)-x(i, j)))/
    !                              max(abs(x(i, j)))
    !                        was computed, where the two maxima are taken
    !                        over i=1, 2, ..., m and j=1, ..., n.
    !
    !                        The value of e is given in the table
    !                        below for some typical values of m and n.
    !
    !                   m (=n)    mperod    nperod        e
    !                   ------    ------    ------      ------
    !
    !                     31        0         0         6.e-14
    !                     31        1         1         4.e-13
    !                     31        1         3         3.e-13
    !                     32        0         0         9.e-14
    !                     32        1         1         3.e-13
    !                     32        1         3         1.e-13
    !                     33        0         0         9.e-14
    !                     33        1         1         4.e-13
    !                     33        1         3         1.e-13
    !                     63        0         0         1.e-13
    !                     63        1         1         1.e-12
    !                     63        1         3         2.e-13
    !                     64        0         0         1.e-13
    !                     64        1         1         1.e-12
    !                     64        1         3         6.e-13
    !                     65        0         0         2.e-13
    !                     65        1         1         1.e-12
    !                     65        1         3         4.e-13
    !
    subroutine genbun(nperod, n, mperod, m, a, b, c, idimy, y, ierror)

        ! Dummy arguments
        integer(ip), intent(in)    :: nperod
        integer(ip), intent(in)    :: n
        integer(ip), intent(in)    :: mperod
        integer(ip), intent(in)    :: m
        integer(ip), intent(in)    :: idimy
        integer(ip), intent(out)   :: ierror
        real(wp),    intent(in)    :: a(:)
        real(wp),    intent(in)    :: b(:)
        real(wp),    intent(in)    :: c(:)
        real(wp),    intent(inout) :: y(idimy,n)

        ! Local variables
        type(FishpackWorkspace)              :: workspace
        type(CenteredCyclicReductionUtility) :: util

        ! Allocate memory
        workspace = genbun_get_workspace(n, m)

        ! Solve system
        associate( rew => workspace%real_workspace )
            call util%genbun_lower_routine(nperod, n, mperod, m, a, b, c, idimy, y, ierror, rew)
        end associate

        ! Release memory
        call workspace%destroy()

    end subroutine genbun

    function genbun_get_workspace(n, m) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: n, m
        type(FishpackWorkspace)  :: return_value

        ! Local variables
        integer(ip)  :: irwk, icwk

        ! Get workspace dimensions for genbun
        call return_value%compute_genbun_workspace_lengths(n, m, irwk)

        ! No need to allocate complex arrays
        icwk = 0

        ! Allocate memory
        call return_value%create(irwk, icwk)

    end function genbun_get_workspace

end module centered_real_linear_systems_solver
