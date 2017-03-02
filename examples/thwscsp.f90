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
program test_hwscsp

    use fishpack

    ! Explicit typing only
    implicit none

    ! Dictionary
    type(FishpackWorkspace)  :: workspace
    integer(ip), parameter   :: M = 36, N = 32
    integer(ip), parameter   :: MP1 = M + 1, NP1 = N + 1
    integer(ip), parameter   :: IDIMF = M + 12
    integer(ip)              :: intl, mbdcnd, nbdcnd, i, j, ierror
    real(wp)                 :: f(IDIMF, NP1), theta(IDIMF)
    real(wp), dimension(NP1) :: bdtf, r
    real(wp), dimension(1)   :: bdts, bdrs, bdrf
    real(wp)                 :: ts, tf, rs, rf, elmbda
    real(wp)                 :: dtheta, dr, pertrb, dphi
    real(wp)                 :: ZERO = 0.0_wp, ONE = 1.0_wp, TWO = 2.0_wp

    ! Initialization flag
    intl = 0

    ! Set domain
    ts = ZERO
    tf = HALF_PI
    rs = ZERO
    rf = ONE

    ! Set boundary conditions
    mbdcnd = 6
    nbdcnd = 5

    ! Set helmholtz constant
    elmbda = ZERO

    ! Set mesh sizes
    dtheta = tf/M
    dr = ONE /N

    ! Generate and store grid points for the purpose of computing the
    ! boundary data and the right side of the equation.
    do i = 1, MP1
        theta(i) = real(i - 1, kind=wp) * dtheta
    end do

    do j = 1, NP1
        r(j) = real(j - 1, kind=wp) * dr
    end do

    ! Generate normal derivative data at equator
    ! In our example, bdts, bdrs, and bdrf are 1-dimensional dummy variables
    bdtf = ZERO

    ! Compute boundary data on the surface of the sphere
    do i = 1, MP1
        f(i, NP1) = cos(theta(i))**4
    end do

    ! Compute right hand side of equation
    block
        real(wp)            :: cost2
        real(wp), parameter :: TWELVE = 12.0_wp

        do i = 1, MP1
            cost2 = cos(theta(i))**2
            f(i, :N) = TWELVE * cost2 * (r(:N)**2)
        end do
    end block

    ! Solve 2D axisymmetric Helmholtz equation on centered grid
    call hwscsp(intl, ts, tf, M, mbdcnd, bdts, bdtf, rs, rf, N, &
        nbdcnd, bdrs, bdrf, elmbda, f, IDIMF, pertrb, ierror, workspace)

    ! Compute discretization error
    block
        real(wp), parameter :: KNOWN_ERROR = 0.799841637481730e-3_wp
        real(wp) :: discretization_error, exact_solution(MP1, N)

        do j = 1, N
            do i = 1, MP1
                exact_solution(i,j) = (cos(theta(i))**4) * (r(j)**4)
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:MP1,:N)))

        call check_output('hwscsp example 1', ierror, KNOWN_ERROR, discretization_error)
    end block

    ! The following program illustrates the use of hwscsp to solve
    ! a three dimensional problem which has longitudinal dependence
    mbdcnd = 2
    nbdcnd = 1
    dphi = PI/72
    elmbda = -TWO *(ONE - cos(dphi))/dphi**2

    ! Compute boundary data on the surface of the sphere
    do i = 1, MP1
        f(i, NP1) = sin(theta(i))
    end do

    ! Compute right side of the equation
    f(:MP1, :N) = ZERO

    ! Solve system
    call hwscsp(intl, ts, tf, M, mbdcnd, bdts, bdtf, rs, rf, N, &
        nbdcnd, bdrs, bdrf, elmbda, f, IDIMF, pertrb, ierror, workspace)

    ! Compute discretization error (fourier coefficients)
    block
        real(wp), parameter :: KNOWN_ERROR = 0.586824289033339e-4_wp
        real(wp) :: exact_solution(MP1,NP1), discretization_error

        discretization_error = ZERO
        do j = 1, NP1
            do i = 1, MP1
                exact_solution(i,j) = r(j) * sin(theta(i))
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:MP1,:NP1)))

        call check_output('hwscsp example 2', ierror, KNOWN_ERROR, discretization_error)
    end block

    ! Release memory
    call workspace%destroy()

end program test_hwscsp
