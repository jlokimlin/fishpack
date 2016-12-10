!
!     file thwscrt.f90
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
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
program thwscrt

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        ip, wp, PI, HALF_PI, hwscrt

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer(ip), parameter   :: M = 40, N = 80
    integer(ip), parameter   :: MP1 = M + 1, NP1 = N + 1
    integer(ip), parameter   :: IDIMF = M + 5
    integer(ip)              :: mbdcnd, nbdcnd, i, j, ierror
    real(wp)                 :: f(IDIMF, NP1), x(MP1)
    real(wp), dimension(NP1) :: bdb, bda, bdc, bdd, y
    real(wp)                 :: a, b, c, d, elmbda, pertrb, discretization_error
    real(wp), parameter      :: PI2 = PI**2
    real(wp), parameter      :: ZERO = 0.0_wp, ONE = 1.0_wp
    real(wp), parameter      :: TWO = 2.0_wp, THREE = 3.0_wp, FOUR = 4.0_wp
    !-----------------------------------------------

    ! Set domain
    a = ZERO
    b = TWO
    c = -ONE
    d = THREE

    ! Set boundary conditions
    mbdcnd = 2
    nbdcnd = 0

    ! Set helmholtz constant
    elmbda = -FOUR

    ! Generate and store grid points for the purpose of computing
    ! boundary data and the right side of the helmholtz equation.
    do i = 1, MP1
        x(i) = TWO * real(i - 1, kind=wp)/M
    end do

    do j = 1, NP1
        y(j) = -ONE + FOUR * real(j - 1, kind=wp)/N
    end do

    ! Generate boundary data. bda, bdc, and bdd are dummy variables.
    do j = 1, NP1
        bdb(j) = FOUR * cos((y(j) + ONE)*HALF_PI)
    end do

    ! Generate right side of equation.
    f = ZERO
    do j = 1, NP1
        do i = 2, MP1
            f(i, j) = (TWO - (FOUR + PI2/4)*(x(i)**2)) * cos((y(j) + ONE)*HALF_PI)
        end do
    end do

    call hwscrt(a, b, M, mbdcnd, bda, bdb, c, d, N, nbdcnd, bdc, bdd, &
        elmbda, f, IDIMF, pertrb, ierror)

    ! Compute discretization error. The exact solution is
    !
    ! u(x, y) = (x**2) * cos((y+1)*(pi/2))
    !
    block
        real(wp) :: local_error, exact_solution

        discretization_error = ZERO
        do i = 1, MP1
            do j = 1, NP1
                exact_solution = (x(i)**2) * cos((y(j)+ONE) * HALF_PI)
                local_error = abs(f(i, j) - exact_solution)
                discretization_error = max(local_error, discretization_error)
            end do
        end do
    end block

    ! Print earlier output from platforms with 64-bit floating point
    ! arithmetic followed by the output from this computer
    write( stdout, '(/a)') '     hwscrt *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 5.36508e-4'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        ' discretization error = ', discretization_error

end program thwscrt
