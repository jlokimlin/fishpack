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
!
program test_cblktri

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack

    ! Explicit typing only
    implicit none

    ! Dictionary
    type(FishpackWorkspace)       :: workspace
    integer(ip), parameter        :: N = 63, M = 50
    integer(ip), parameter        :: IDIMY = 75, NT = 105
    integer(ip)                   :: iflg, np, mp, i, j, ierror
    real(wp), dimension(NT)       :: an, bn, cn, t
    real(wp), dimension(IDIMY)    :: s
    real(wp)                      :: ds, dt
    real(wp), parameter           :: ZERO = 0.0_wp, ONE = 1.0_wp, TWO = 2.0_wp
    complex(wp)                   :: y(IDIMY,NT)
    complex(wp), dimension(IDIMY) :: am, bm, cm
    complex(wp), parameter        :: IMAGINARY_UNIT = cmplx(ZERO, ONE, kind=wp)

    ! Set boundary conditions
    np = 1
    mp = 1

    ! Set mesh sizes
    ds = ONE/(M + 1)
    dt = ONE/(N + 1)

    ! Generate and store grid points for the purpose of computing the
    ! coefficients and the array y.
    do i = 1, M
        s(i) = real(i, kind=wp) * ds
    end do

    do j = 1, N
        t(j) = real(j, kind=wp) * dt
    end do

    ! Compute the coefficients am, bm, cm corresponding to the s direction
    block
        real(wp) :: half_ds, two_ds
        real(wp) :: temp1, temp2, temp3

        half_ds = ds/2
        two_ds = TWO * ds

        do i = 1, M
            temp1 = ONE /(s(i)*two_ds)
            temp2 = ONE /((s(i)-half_ds) * two_ds)
            temp3 = ONE /((s(i)+half_ds) * two_ds)
            am(i) = cmplx(temp1*temp2, ZERO, kind=wp)
            cm(i) = cmplx(temp1*temp3, ZERO, kind=wp)
            bm(i) = (-(am(i)+cm(i))) - IMAGINARY_UNIT
        end do
    end block

    ! Compute the coefficients an, bn, cn corresponding to the t direction
    block
        real(wp) :: half_dt, two_dt
        real(wp) :: temp1, temp2, temp3

        half_dt = dt/2
        two_dt = TWO * dt

        do j = 1, N
            temp1 = ONE/(t(j) * two_dt)
            temp2 = ONE/((t(j) - half_dt) * two_dt)
            temp3 = ONE/((t(j) + half_dt) * two_dt)
            an(j) = temp1 * temp2
            cn(j) = temp1 * temp3
            bn(j) = -(an(j) + cn(j))
        end do
    end block

    ! Compute right side of equation
    do j = 1, N
        y(:M, j) = 3.75_wp * s(:M) * t(j) * (s(:M)**4 + t(j)**4) &
            - IMAGINARY_UNIT * (s(:M)*t(j))**5
    end do

    ! The nonzero boundary conditions enter the linear system via
    ! the right side y(i, j). if the equations (3) given above
    ! are evaluated at i=m and j=1, ..., n then the term cm(m)*u(m+1, j)
    ! is known from the boundary condition to be cm(m)*t(j)**5.
    ! therefore this term can be included in the right side y(m, j).
    ! the same analysis applies at j=n and i=1, .., m. note that the
    ! corner at j=n, i=m includes contributions from both boundaries.
    y(M,:N) = y(M,:N) - cm(M) * (t(:N)**5)
    y(:M,N) = y(:M,N) - cn(N) * (s(:M)**5)

    ! Initialize workspace and lower solver routines
    iflg = 0
    call cblktri(iflg, np, N, an, bn, cn, mp, M, am, bm, cm, IDIMY, y, ierror, workspace)

    ! Solve complex block tridiagonal linear system
    iflg = iflg + 1
    do while(iflg <= 1)
        call cblktri(iflg, np, N, an, bn, cn, mp, M, am, bm, cm, IDIMY, y, ierror, workspace)
        iflg = iflg + 1
    end do

    ! Compute discretization error. The exact solution is
    ! u(s,t) = (st)**5
    block
        real(wp), parameter :: KNOWN_ERROR = 0.164571992877414e-004_wp
        real(wp) :: discretization_error
        real(wp) :: exact_solution(M, N)

        do j = 1, N
            do i = 1, M
                exact_solution(i,j) = (s(i)*t(j))**5
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution-y(:M,:N)))

        call check_output('cblktri', ierror, KNOWN_ERROR, discretization_error)
    end block

    ! Release memory
    call workspace%destroy()

end program test_cblktri
