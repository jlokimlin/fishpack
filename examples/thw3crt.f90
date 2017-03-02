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
program test_hw3crt

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer(ip), parameter   :: L = 10, M = 40, N = 15
    integer(ip), parameter   :: LP1 = L + 1, MP1 = M + 1, NP1 = N + 1
    integer(ip), parameter   :: LDIMF = LP1, MDIMF = MP1
    integer(ip)              :: lbdcnd, mbdcnd, nbdcnd, i, j, k, ierror
    real(wp)                 :: f(LDIMF, MDIMF, NP1), pertrb
    real(wp)                 :: bdzf(LP1,MP1), x(LP1), y(MP1), z(NP1)
    real(wp), dimension(1,1) :: bdxs, bdxf, bdys, bdyf, bdzs
    real(wp)                 :: elmbda, xs, xf, ys, yf, zs, zf
    real(wp)                 :: dx, dy, dz
    real(wp), parameter      :: ZERO = 0.0_wp, ONE = 1.0_wp
    real(wp), parameter      :: THREE = 3.0_wp, FOUR = 4.0_wp
    !-----------------------------------------------

    ! Set domain
    xs = ZERO
    xf = ONE
    ys = ZERO
    yf = TWO_PI
    zs = ZERO
    zf = HALF_PI

    ! Set boundary conditions
    lbdcnd = 1
    mbdcnd = 0
    nbdcnd = 2

    ! Set helmholtz constant
    elmbda = -THREE

    ! Set mesh sizes
    dx = (xf - xs)/L
    dy = (yf - ys)/M
    dz = (zf - zs)/N

    ! We define the grid points for later use.
    do i = 1, LP1
        x(i) = xs + real(i - 1, kind=wp) * dx
    end do

    do j = 1, MP1
        y(j) = ys + real(j - 1, kind=wp) * dy
    end do

    do k = 1, NP1
        z(k) = zs + real(k - 1, kind=wp) * dz
    end do

    ! We define the array of derivative boundary values.
    ! In our example all other boundary arrays are
    ! 1-by-1 dimensional dummy variables.
    do j = 1, MP1
        do i = 1, LP1
            bdzf(i, j) = (-x(i)**4) * sin(y(j))
        end do
    end do

    ! We define the function boundary values in the f array.
    do k = 1, NP1
        do j = 1, MP1
            f(1, j, k) = ZERO
            f(LP1, j, k) = sin(y(j))*cos(z(k))
        end do
    end do

    do j = 1, MP1
        do i = 1, LP1
            f(i, j, 1) = (x(i)**4)*sin(y(j))
        end do
    end do

    ! We now define the values of the right side of the helmholtz equation.
    do i = 2, L
        do j = 1, MP1
            do k = 2, NP1
                f(i, j, k) = FOUR*x(i)**2*(THREE - x(i)**2)*sin(y(j))*cos(z(k))
            end do
        end do
    end do

    ! Solve 3D Helmholtz in cartesian coordinates on centered grid
    call hw3crt(xs, xf, L, lbdcnd, bdxs, bdxf, ys, yf, M, mbdcnd, &
        bdys, bdyf, zs, zf, N, nbdcnd, bdzs, bdzf, elmbda, LDIMF, MDIMF, &
        f, pertrb, ierror)

    ! Compute discretization error. The exact solution is
    !
    ! u(x,y,z) = (x**4) * sin(y) * cos(z)
    !
    block
        real(wp), parameter :: KNOWN_ERROR = 0.964801952068239e-002_wp
        real(wp) :: discretization_error
        real(wp) :: exact_solution(LP1, MP1, NP1)

        do k = 1, NP1
            do j = 1, MP1
                do i = 1, LP1
                    exact_solution(i,j,k) = (x(i)**4) * sin(y(j)) * cos(z(k))
                end do
            end do
        end do

        discretization_error = maxval(abs(exact_solution - f(:LP1,:MP1, :NP1)))

        call check_output('hw3crt', ierror, KNOWN_ERROR, discretization_error)
    end block

end program test_hw3crt
