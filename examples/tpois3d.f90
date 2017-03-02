!
!     file tpois3d.f90
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
program test_pois3d

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer(ip), parameter :: L = 30, M = 30, N = 10
    integer(ip), parameter :: LDIMF = L + 2, MDIMF = M + 3
    integer(ip)            :: lperod, mperod, nperod, i, j, k, ierror
    real(wp)               :: x(L), y(M), f(LDIMF, MDIMF, N)
    real(wp), dimension(M) :: a, b, c, z
    real(wp)               :: dx, c1, dy, c2, dz, dz2
    real(wp), parameter    :: ZERO = 0.0_wp, ONE = 1.0_wp, TWO = 2.0_wp
    !-----------------------------------------------

    ! Set boundary conditions
    lperod = 0
    mperod = 0
    nperod = 1

    ! Set mesh sizes
    dx = TWO_PI/L
    dy = TWO_PI/M
    dz = ONE/N
    dz2 = ONE/dz**2

    ! Define constants
    c1 = ONE/dx**2
    c2 = ONE/dy**2

    ! Generate grid points for later use
    do i = 1, L
        x(i) = -PI + real(i - 1, kind=wp)*dx
    end do

    do j = 1, M
        y(j) = -PI + real(j - 1, kind=wp)*dy
    end do

    ! Generate coefficients
    block
        real(wp) :: t, t2

        a(1) = ZERO
        b(1) = -TWO *dz2
        c(1) = -b(1)
        z(1) = ZERO
        do k = 2, N
            z(k) = real(k - 1, kind=wp)*dz
            t = ONE + z(k)
            t2 = t**2
            a(k) = t2 * dz2 + t/dz
            b(k) = -TWO * t2 * dz2
            c(k) = t2 * dz2 - t/dz
        end do
    end block

    ! Generate right hand side of equation
    block
        real(wp), parameter :: EIGHT = 8.0_wp
        real(wp), parameter :: TEN = 10.0_wp
        real(wp), parameter :: SIXTEEN = 16.0_wp

        do i = 1, L
            do j = 1, M
                do k = 2, N
                    f(i, j, k) = TWO *sin(x(i))*sin(y(j))*(ONE + z(k))**4
                end do
            end do
        end do
        do i = 1, L
            do j = 1, L
                f(i,j,1) = (TEN + (EIGHT/dz)) * sin(x(i)) * sin(y(j))
                f(i,j,N) = f(i,j,N) - c(N) * SIXTEEN * sin(x(i)) * sin(y(j))
            end do
        end do
        c(N) = ZERO
    end block

    ! Solve 3D real linear system on centered grid
    call pois3d(lperod, L, c1, mperod, M, c2, nperod, N, a, b, c, &
        LDIMF, MDIMF, f, ierror)

    ! Compute discretization error. The exact solution is
    !
    ! u(x, y, z) = sin(x)*sin(y)*(1+z)**4
    !
    block
        real(wp), parameter :: KNOWN_ERROR = 0.293277049861086e-001_wp
        real(wp) :: discretization_error
        real(wp) :: exact_solution(L,M,N)

        do k = 1, N
            do j = 1, M
                do i = 1, L
                    exact_solution(i,j,k) = sin(x(i)) * sin(y(j)) * (ONE + z(k))**4
                end do
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:L,:M,:N)))

        call check_output('pois3d', ierror, KNOWN_ERROR, discretization_error)
    end block

end program test_pois3d
