!     file tcblktri.f90
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
!
program tcblktri

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        FishpackSolver, &
        FishpackWorkspace

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    type(FishpackSolver)       :: solver
    type(FishpackWorkspace)    :: workspace
    integer(ip), parameter     :: NS = 75, NT = 105, N = 63, M = 50
    integer(ip)                :: iflg, np, mp, i, j, ierror
    real(wp), dimension(NT)    :: an, bn, cn, t
    real(wp), dimension(NS)    :: s
    real(wp)                   :: ds, dt, half_ds, two_ds
    real(wp)                   :: half_dt, two_dt, discretization_error, z
    real(wp), parameter        :: ZERO = 0.0_wp, ONE = 1.0_wp, TWO = 2.0_wp
    complex(wp)                :: y(NS,NT)
    complex(wp), dimension(NS) :: am, bm, cm
    complex(wp), parameter     :: IMAGINARY_UNIT = cmplx(ZERO, ONE, kind=wp)
    !-----------------------------------------------

    ! Initialize flag
    iflg = 0

    ! Set boundary conditions
    np = 1
    mp = 1
    !
    !     generate and store grid points for the purpose of computing the
    !     coefficients and the array y.
    !
    ds = ONE/(M + 1)
    do i = 1, M
        s(i) = real(i, kind=wp)*ds
    end do

    dt = ONE/(N + 1)
    do j = 1, N
        t(j) = real(j, kind=wp)*dt
    end do
    !
    !     compute the coefficients am, bm, cm corresponding to the s direction
    !
    half_ds = ds/2
    two_ds = TWO * ds
    do i = 1, M
        associate( &
            temp1 => ONE /(s(i)*two_ds),&
            temp2 => ONE /((s(i)-half_ds)*two_ds), &
            temp3 => ONE /((s(i)+half_ds)*two_ds) &
            )
            am(i) = cmplx(temp1*temp2, ZERO, kind=wp)
            cm(i) = cmplx(temp1*temp3, ZERO, kind=wp)
            bm(i) = (-(am(i)+cm(i))) - IMAGINARY_UNIT
        end associate
    end do
    !
    !     compute the coefficients an, bn, cn corresponding to the t direction
    !
    half_dt = dt/2
    two_dt = TWO * dt
    do j = 1, N
        associate( &
            temp1 => ONE/(t(j)*two_dt), &
            temp2 => ONE/((t(j)-half_dt)*two_dt), &
            temp3 => ONE/((t(j)+half_dt)*two_dt) &
            )
            an(j) = temp1*temp2
            cn(j) = temp1*temp3
            bn(j) = -(an(j)+cn(j))
        end associate
    end do
    !
    !     compute right side of equation
    !
    do j = 1, N
        y(:M, j) = 3.75_wp * s(:M) * t(j) * (s(:M)**4 + t(j)**4) &
            - IMAGINARY_UNIT * (s(:M)*t(j))**5
    end do
    !
    !     the nonzero boundary conditions enter the linear system via
    !     the right side y(i, j). if the equations (3) given above
    !     are evaluated at i=m and j=1, ..., n then the term cm(m)*u(m+1, j)
    !     is known from the boundary condition to be cm(m)*t(j)**5.
    !     therefore this term can be included in the right side y(m, j).
    !     the same analysis applies at j=n and i=1, .., m. note that the
    !     corner at j=n, i=m includes contributions from both boundaries.
    !
    y(M,:N) = y(M,:N) - cm(M)*t(:N)**5
    y(:M,N) = y(:M,N) - cn(N)*s(:M)**5

    call solver%cblktri(iflg, np, N, an, bn, cn, mp, M, am, bm, cm, NS, y, ierror, workspace)

    iflg = iflg + 1
    do while(iflg <= 1)
        ! Solve system
        call solver%cblktri(iflg, np, N, an, bn, cn, mp, M, am, bm, cm, NS, y, ierror, workspace)

        ! Increment iteration flag
        iflg = iflg + 1
    end do

    ! Compute error
    discretization_error = ZERO
    do j = 1, N
        do i = 1, M
            z = abs(y(i, j)-(s(i)*t(j))**5)
            discretization_error = max(z, discretization_error)
        end do
    end do

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     cblktri *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 1.6457e-05'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        ' discretization error = ', discretization_error

    !
    !==> Release memory
    !
    call workspace%destroy()

end program tcblktri
