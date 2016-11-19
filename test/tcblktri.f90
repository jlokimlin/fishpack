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
    type(FishpackSolver)               :: solver
    type(FishpackWorkspace)            :: workspace
    integer(ip), parameter             :: IDIMY = 75
    integer(ip)                        :: iflg, np, n, mp, m, i, j, ierror
    real(wp), dimension(105)           :: an, bn, cn
    real(wp), dimension(IDIMY)         :: s
    real(wp), dimension(105)           :: t
    real(wp)                           :: ds, dt, half_ds, two_ds
    real(wp)                           :: temp1, temp2, temp3, half_dt, two_dt, discretization_error, z
    complex(wp), dimension(IDIMY, 105) :: y
    complex(wp), dimension(IDIMY)      :: am, bm, cm
    !-----------------------------------------------
    !
    iflg = 0
    np = 1
    n = 63
    mp = 1
    m = 50
    !
    !     generate and store grid points for the purpose of computing the
    !     coefficients and the array y.
    !
    ds = 1.0_wp/(m + 1)

    do i = 1, m
        s(i) = real(i, kind=wp)*ds
    end do

    dt = 1.0_wp/(n + 1)

    do j = 1, n
        t(j) = real(j, kind=wp)*dt
    end do
    !
    !     compute the coefficients am, bm, cm corresponding to the s direction
    !
    half_ds = ds/2
    two_ds = ds + ds
    do i = 1, m
        temp1 = 1.0_wp /(s(i)*two_ds)
        temp2 = 1.0_wp /((s(i)-half_ds)*two_ds)
        temp3 = 1.0_wp /((s(i)+half_ds)*two_ds)
        am(i) = cmplx(temp1*temp2, 0.0_wp, kind=wp)
        cm(i) = cmplx(temp1*temp3, 0.0_wp, kind=wp)
        bm(i) = (-(am(i)+cm(i))) - cmplx(0.0_wp, 1.0_wp, kind=wp)
    end do
    !
    !     compute the coefficients an, bn, cn corresponding to the t direction
    !
    half_dt = dt/2
    two_dt = 2.0_wp * dt
    do j = 1, n
        temp1 = 1.0_wp/(t(j)*two_dt)
        temp2 = 1.0_wp/((t(j)-half_dt)*two_dt)
        temp3 = 1.0_wp/((t(j)+half_dt)*two_dt)
        an(j) = temp1*temp2
        cn(j) = temp1*temp3
        bn(j) = -(an(j)+cn(j))
    end do
    !
    !     compute right side of equation
    !
    do j = 1, n
        y(:m, j) = 3.75_wp*s(:m)*t(j)*(s(:m)**4+t(j)**4) &
            - cmplx(0.0_wp, 1.0_wp, kind=wp)*(s(:m)*t(j))**5
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
    y(m, :n) = y(m, :n) - cm(m)*t(:n)**5
    y(:m, n) = y(:m, n) - cn(n)*s(:m)**5

    call solver%cblktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, IDIMY, y, ierror, workspace)

    iflg = iflg + 1

    do while(iflg <= 1)

        ! Solve system
        call solver%cblktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, IDIMY, y, ierror, workspace)

        ! Increment iteration flag
        iflg = iflg + 1

    end do

    discretization_error = 0.0_wp
    do j = 1, n
        do i = 1, m
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
