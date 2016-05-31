!     file tblktri.f
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
!     program to illustrate the use of subroutine blktri to
!     solve the equation
!
!     .5/s*(d/ds)(.5/s*du/ds)+.5/t*(d/dt)(.5/t*du/dt)
!                                                          (1)
!                   = 15/4*s*t*(s**4+t**4)
!
!     on the rectangle 0 < s < 1 and 0 < t < 1
!     with the boundary conditions
!
!     u(0, t) = 0
!                            0 <= t <= 1
!     u(1, t) = t**5
!
!     and
!
!     u(s, 0) = 0
!                            0 <= s <= 1
!     u(s, 1) = s**5
!
!     the exact solution of this problem is u(s, t) = (s*t)**5
!
!     define the integers m = 50 and n = 63. then define the
!     grid increments deltas = 1/(m+1) and deltat = 1/(n+1).
!
!     the grid is then given by s(i) = i*deltas for i = 1, ..., m
!     and t(j) = j*deltat for j = 1, ..., n.
!
!     the approximate solution is given as the solution to
!     the following finite difference approximation of equation (1).
!
!     .5/(s(i)*deltas)*((u(i+1, j)-u(i, j))/(2*s(i+.5)*deltas)
!                     -(u(i, j)-u(i-1, j))/(2*s(i-.5)*deltas))
!     +.5/(t(i)*deltat)*((u(i, j+1)-u(i, j))/(2*t(i+.5)*deltat) (2)
!                     -(u(i, j)-u(i, j-1))/(2*t(i-.5)*deltat))
!               = 15/4*s(i)*t(j)*(s(i)**4+t(j)**4)
!
!             where s(i+.5) = .5*(s(i+1)+s(i))
!                   s(i-.5) = .5*(s(i)+s(i-1))
!                   t(i+.5) = .5*(t(i+1)+t(i))
!                   t(i-.5) = .5*(t(i)+t(i-1))
!
!     the approach is to write equation (2) in the form
!
!     am(i)*u(i-1, j)+bm(i, j)*u(i, j)+cm(i)*u(i+1, j)
!       +an(j)*u(i, j-1)+bn(j)*u(i, j)+cn(j)*u(i, j+1)      (3)
!           = y(i, j)
!
!     and then call subroutine blktri to determine u(i, j)
!
program tblktri

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        FishpackSolver, &
        FishpackWorkspace

    ! Explicit typing only
    implicit None

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    type (FishpackSolver)     :: solver
    type (FishpackWorkspace)  :: workspace
    integer (ip), parameter   :: IDIMY = 75
    integer (ip)              :: iflg, np, n, mp, m, i, j, ierror
    real (wp), dimension(IDIMY, 105) :: y
    real (wp), dimension(IDIMY) :: am, bm, cm
    real (wp), dimension(105) :: an, bn, cn
    real (wp), dimension(IDIMY) :: s
    real (wp), dimension(105) :: t
    real (wp) :: ds, dt, discretization_error
    !-----------------------------------------------

    iflg = 0
    np = 1
    n = 63
    mp = 1
    m = 50
    !
    !     Generate and store grid points for the purpose of computing the
    !     coefficients and the array y.
    !
    ds = 1.0_wp/(m + 1)

    do i = 1, m
        s(i) = real(i, kind=wp) * ds
    end do

    dt = 1.0_wp/(n + 1)

    do j = 1, n
        t(j) = real(j, kind=wp) * dt
    end do
    !
    !     Compute the coefficients am, bm, cm corresponding to the s direction
    !
    associate( &
        half_ds => ds/2, &
        two_ds => 2.0_wp * ds &
        )
        do i = 1, m
            associate( &
                temp1 => 1.0_wp/(s(i)*two_ds), &
                temp2 => 1.0_wp/((s(i)-half_ds)*two_ds), &
                temp3 => 1.0_wp/((s(i)+half_ds)*two_ds) &
                )
                am(i) = temp1*temp2
                cm(i) = temp1*temp3
                bm(i) = -(am(i)+cm(i))
            end associate
        end do
    end associate

    !
    !     compute the coefficients an, bn, cn corresponding to the t direction
    !
    associate( &
        half_dt => dt/2, &
        two_dt => 2.0_wp * dt &
        )

        do j = 1, n
            associate( &
                temp1 => 1.0_wp / ( t(j) * two_dt ), &
                temp2 => 1.0_wp / ( (t(j)-half_dt) * two_dt ), &
                temp3 => 1.0_wp / ( (t(j)+half_dt) * two_dt ) &
                )
                an(j) = temp1*temp2
                cn(j) = temp1*temp3
                bn(j) = -(an(j)+cn(j))
            end associate
        end do
    end associate

    !
    !     Compute right side of equation
    !
    do j = 1, n
        y(:m, j) = 3.75_wp * s(:m) * t(j) * ( s(:m)**4 + t(j)**4 )
    end do
    !
    !     The nonzero boundary conditions enter the linear system via
    !     the right side y(i, j). if the equations (3) given above
    !     are evaluated at i=m and j=1, ..., n then the term cm(m)*u(m+1, j)
    !     is known from the boundary condition to be cm(m)*t(j)**5.
    !     therefore this term can be included in the right side y(m, j).
    !     the same analysis applies at j=n and i=1, .., m. note that the
    !     corner at j=n, i=m includes contributions from both boundaries.
    !
    y(m, :n) = y(m, :n) - cm(m)*t(:n)**5
    y(:m, n) = y(:m, n) - cn(n)*s(:m)**5
    !
    !==> Determine the approximate solution u(i, j)
    !
    call solver%blktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, IDIMY, y, ierror, workspace)

    iflg = iflg + 1

    do while(iflg - 1 <= 0)

        ! Solver system
        call solver%blktri(iflg, np, n, an, bn, cn, mp, m, am, bm, cm, IDIMY, &
            y, ierror, workspace )

        ! Increment flag
        iflg = iflg + 1

    end do

    discretization_error = 0.0_wp

    do j = 1, n
        do i = 1, m
            associate( local_error => abs(y(i, j)-(s(i)*t(j))**5) )

                discretization_error = max(local_error, discretization_error)

            end associate
        end do
    end do

    !     Print earlier output from platform with 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     blktri *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  discretization error = 1.6478e-05'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') '     ierror =', ierror, ' discretization error = ', &
        discretization_error

    !
    !==> Release memory
    !
    call workspace%destroy()

end program tblktri
