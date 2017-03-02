!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *                                                               *
! *                  copyright (c) 2005 by UCAR                   *
! *                                                               *
! *       University Corporation for Atmospheric Research         *
! *                                                               *
! *                      all rights reserved                      *
! *                                                               *
! *                         Fishpack                              *
! *                                                               *
! *                      A Package of Fortran                     *
! *                                                               *
! *                Subroutines and Example Programs               *
! *                                                               *
! *               for Modeling Geophysical Processes              *
! *                                                               *
! *                             by                                *
! *                                                               *
! *        John Adams, Paul Swarztrauber and Roland Sweet         *
! *                                                               *
! *                             of                                *
! *                                                               *
! *         the National Center for Atmospheric Research          *
! *                                                               *
! *                Boulder, Colorado  (80307)  U.S.A.             *
! *                                                               *
! *                   which is sponsored by                       *
! *                                                               *
! *              the National Science Foundation                  *
! *                                                               *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! Purpose:
!
! To illustrate the use of subroutine blktri to
! solve the equation
!
! 0.5/s*(d/ds)(0.5/s*du/ds)+0.5/t*(d/dt)(0.5/t*du/dt)
!                                                      (1)
!               = 15/4*s*t*(s**4+t**4)
!
! on the rectangle 0 < s < 1 and 0 < t < 1
! with the boundary conditions
!
! u(0, t) = 0
!                        0 <= t <= 1
! u(1, t) = t**5
!
! and
!
! u(s, 0) = 0
!                        0 <= s <= 1
! u(s, 1) = s**5
!
! the exact solution of this problem is u(s, t) = (s*t)**5
!
! define the integers m = 50 and n = 63. then define the
! grid increments deltas = 1/(m+1) and deltat = 1/(n+1).
!
! the grid is then given by s(i) = i*deltas for i = 1, ..., m
! and t(j) = j*deltat for j = 1, ..., n.
!
! the approximate solution is given as the solution to
! the following finite difference approximation of equation (1).
!
! .5/(s(i)*deltas)*((u(i+1, j)-u(i, j))/(2*s(i+.5)*deltas)
!                 -(u(i, j)-u(i-1, j))/(2*s(i-.5)*deltas))
! +.5/(t(i)*deltat)*((u(i, j+1)-u(i, j))/(2*t(i+.5)*deltat) (2)
!                 -(u(i, j)-u(i, j-1))/(2*t(i-.5)*deltat))
!           = 15/4*s(i)*t(j)*(s(i)**4+t(j)**4)
!
!         where s(i+.5) = .5*(s(i+1)+s(i))
!               s(i-.5) = .5*(s(i)+s(i-1))
!               t(i+.5) = .5*(t(i+1)+t(i))
!               t(i-.5) = .5*(t(i)+t(i-1))
!
! the approach is to write equation (2) in the form
!
! am(i)*u(i-1, j)+bm(i, j)*u(i, j)+cm(i)*u(i+1, j)
!   +an(j)*u(i, j-1)+bn(j)*u(i, j)+cn(j)*u(i, j+1)      (3)
!       = y(i, j)
!
! and then call subroutine blktri to determine u(i, j)
!
program test_blktri

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack

    ! Explicit typing only
    implicit none

    ! Dictionary
    type(FishpackWorkspace)    :: workspace
    integer(ip), parameter     :: M = 50, N = 63
    integer(ip), parameter     :: IDIMY = 75, NT = 105
    integer(ip)                :: iflg, np, mp, i, j, ierror
    real(wp)                   :: y(IDIMY, NT)
    real(wp), dimension(IDIMY) :: am, bm, cm, s
    real(wp), dimension(NT)    :: an, bn, cn, t
    real(wp)                   :: ds, dt
    real(wp), parameter        :: ONE = 1.0_wp, TWO = 2.0_wp

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
            temp1 = ONE/(s(i) * two_ds)
            temp2 = ONE/((s(i) - half_ds) * two_ds)
            temp3 = ONE/((s(i) + half_ds) * two_ds)
            am(i) = temp1 * temp2
            cm(i) = temp1 * temp3
            bm(i) = -(am(i) + cm(i))
        end do
    end block

    ! Compute the coefficients an, bn, cn corresponding to the t direction
    block
        real(wp) :: half_dt, two_dt
        real(wp) :: temp1, temp2, temp3

        half_dt = dt/2
        two_dt = TWO * dt

        do j = 1, N
            temp1 = ONE / (t(j) * two_dt)
            temp2 = ONE / ((t(j)-half_dt) * two_dt)
            temp3 = ONE / ((t(j)+half_dt) * two_dt)
            an(j) = temp1 * temp2
            cn(j) = temp1 * temp3
            bn(j) = -(an(j) + cn(j))
        end do
    end block

    ! Compute right hand side of equation
    do j = 1, n
        y(:M, j) = 3.75_wp * s(:M) * t(j) * (s(:M)**4 + t(j)**4)
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
    call blktri(iflg, np, N, an, bn, cn, mp, M, am, bm, cm, IDIMY, y, ierror, workspace)

    iflg = iflg + 1
    do while (iflg <= 1)

        ! Solve real block tridiagonal linear system
        call blktri(iflg, np, N, an, bn, cn, mp, M, am, bm, cm, IDIMY, &
            y, ierror, workspace)

        iflg = iflg + 1 ! Increment flag
    end do

    ! Compute discretization error. The exact solution is
    !
    ! u(s,t) = (st)**5
    block
        real(wp), parameter :: KNOWN_ERROR = 0.164778571363905e-4_wp
        real(wp) :: discretization_error
        real(wp) :: exact_solution(M, N)

        do j = 1, N
            do i = 1, M
                exact_solution(i,j) = (s(i)*t(j))**5
            end do
        end do

        ! Set discretization error
        discretization_error = maxval(abs(exact_solution-y(:M,:N)))

        call check_output('blktri', ierror, KNOWN_ERROR, discretization_error)
    end block

    ! Release memory
    call workspace%destroy()

end program test_blktri
