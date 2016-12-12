!
!     file thwsssp.f90
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
!     program to illustrate the use of hwsssp
!
program thwsssp

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        ip, wp, HALF_PI, TWO_PI, hwsssp

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer(ip), parameter     :: M = 18, N = 72
    integer(ip), parameter     :: MP1 = M + 1, NP1 = N + 1
    integer(ip), parameter     :: IDIMF = MP1
    integer(ip)                :: mbdcnd, nbdcnd, i, j, ierror
    real(wp)                   :: f(IDIMF, NP1), sint(IDIMF)
    real(wp), dimension(NP1)   :: bdtf, bdts, bdps, bdpf, sinp
    real(wp)                   :: ts, tf, ps, pf, elmbda
    real(wp)                   :: dtheta, dphi, pertrb
    real(wp)                   :: ZERO = 0.0_wp
    !-----------------------------------------------

    ! Set domain
    ts = ZERO
    tf = HALF_PI
    ps = ZERO
    pf = TWO_PI

    ! Set boundary conditions
    nbdcnd = 0
    mbdcnd = 6

    ! Set helmholtz constant
    elmbda = ZERO

    ! Set mesh sizes
    dtheta = tf/M
    dphi = TWO_PI/N

    ! Generate sines for use in subsequent computations
    do i = 1, MP1
        sint(i) = sin(real(i - 1, kind=wp) * dtheta)
    end do

    do j = 1, NP1
        sinp(j) = sin(real(j - 1, kind=wp) * dphi)
    end do

    ! Compute right hand side of equation and store in f
    block
        real(wp), parameter :: TWO = 2.0_wp
        real(wp), parameter :: SIX = 6.0_wp

        do j = 1, NP1
            f(:MP1, j) = TWO - SIX * (sint(:MP1)*sinp(j))**2
        end do
    end block

    !
    !     store derivative data at the equator
    !
    bdtf(:NP1) = ZERO

    ! Solve 2D Helmholtz in spherical coordinates on a centered grid
    call hwsssp(ts, tf, M, mbdcnd, bdts, bdtf, ps, pf, N, nbdcnd, &
        bdps, bdpf, elmbda, f, IDIMF, pertrb, ierror)

    ! Compute discretization error.
    ! Since problem is singular, the solution must be normalized.
    block
        real(wp) :: discretization_error, exact_solution(MP1, NP1)

        do j = 1, NP1
            do i = 1, MP1
                exact_solution(i,j) = (sint(i)*sinp(j))**2 - f(1, 1)
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution - f(:MP1,:NP1)))

        ! Print earlier output from platforms with 64-bit floating point
        ! arithmetic followed by the output from this computer
        write( stdout, '(/a)') '     hwsssp *** TEST RUN *** '
        write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(a)') '     ierror = 0,  discretization error = 3.38107e-3'
        write( stdout, '(a)') '     The output from your computer is: '
        write( stdout, '(a,i3,a,1pe15.6/)') &
            '      ierror =', ierror, &
            ' discretization error = ', discretization_error
    end block

end program thwsssp
