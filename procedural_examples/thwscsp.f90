!
!     file thwscsp.f90
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
program thwscsp

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        FishpackWorkspace, ip, wp, HALF_PI, hwscsp, PI

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    type(FishpackWorkspace)  :: workspace
    integer(ip), parameter   :: M = 36, N = 32
    integer(ip), parameter   :: MP1 = M + 1, NP1 = N + 1
    integer(ip), parameter   :: IDIMF = M + 12
    integer(ip)              :: intl, mbdcnd, nbdcnd, i, j, ierror
    real(wp)                 :: f(IDIMF, NP1), theta(IDIMF)
    real(wp), dimension(NP1) :: bdtf, bdts, bdrs, bdrf, r
    real(wp)                 :: ts, tf, rs, rf, elmbda
    real(wp)                 :: dtheta, dr, pertrb
    real(wp)                 :: discretization_error, dphi
    real(wp)                 :: ZERO = 0.0_wp, ONE = 1.0_wp, TWO = 2.0_wp
    !-----------------------------------------------

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
        real(wp) :: cost4, exact_solution, local_error

        discretization_error = ZERO
        do j = 1, N
            do i = 1, MP1
                cost4 = cos(theta(i))**4
                exact_solution = cost4 * (r(j)**4)
                local_error = abs(f(i, j) - exact_solution)
                discretization_error = max(local_error, discretization_error)
            end do
        end do
    end block

    ! Print earlier output from platforms with 64 bit floating point
    ! arithmetic followed by the output from this computer
    write( stdout, '(/a)') '     hwscsp *** TEST RUN, EXAMPLE 1 *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 7.9984e-4 '
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        '     discretization error =', discretization_error
    !
    ! The following program illustrates the use of hwscsp to solve
    ! a three dimensional problem which has longitudinal dependence
    !
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
        real(wp) :: exact_solution, local_error

        discretization_error = ZERO
        do j = 1, NP1
            do i = 1, MP1
                exact_solution = r(j) * sin(theta(i))
                local_error = abs(f(i, j) - exact_solution)
                discretization_error = max(local_error, discretization_error)
            end do
        end do
    end block

    ! Print earlier output from platforms with 64-bit floating point
    ! arithmetic followed by the output from this computer
    write( stdout, '(/a)') '     hwscsp *** TEST RUN, EXAMPLE 2 *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0, discretization error = 5.8682e-5 '
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        '     discretization error =', discretization_error

    ! Release memory
    call workspace%destroy()

end program thwscsp
