!
!     file thw3crt.f
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
program thw3crt

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        FishpackSolver

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    type(FishpackSolver)    :: solver
    integer(ip) :: lbdcnd, mbdcnd, nbdcnd, l, m, n
    integer(ip) :: ldimf, mdimf, lp1, i
    integer(ip) :: mp1, j, np1, k, ierror
    real(wp), dimension(11, 41, 16) :: f
    real(wp), dimension(11, 41) :: bdzf, bdxs, bdxf, bdys, bdyf, bdzs
    real(wp), dimension(11) :: x
    real(wp), dimension(41) :: y
    real(wp), dimension(16) :: z
    real(wp) :: elmbda, xs, xf, ys, pi, yf, zs, zf, dx, dy, dz
    real(wp) :: pertrb, discretization_error, t
    !-----------------------------------------------

    !
    !        from the description of the problem given above, we define
    !     the following quantities
    !
    elmbda = -3.0_wp
    xs = 0.0_wp
    xf = 1.0_wp
    lbdcnd = 1
    ys = 0.0_wp
    pi = acos(-1.0_wp)
    yf = 2.0_wp*pi
    mbdcnd = 0
    zs = 0.0_wp
    zf = acos(0.0_wp)
    nbdcnd = 2
    l = 10
    m = 40
    n = 15
    !
    !     from the dimension statement above we define
    !
    ldimf = 11
    mdimf = 41
    !
    !     we define the grid points for later use.
    !
    lp1 = l + 1
    dx = (xf - xs)/l
    do i = 1, lp1
        x(i) = xs + real(i - 1, kind=wp)*dx
    end do

    mp1 = m + 1
    dy = (yf - ys)/m
    do j = 1, mp1
        y(j) = ys + real(j - 1, kind=wp)*dy
    end do

    np1 = n + 1
    dz = (zf - zs)/n
    do k = 1, np1
        z(k) = zs + real(k - 1, kind=wp)*dz
    end do
    !
    !     we define the array of derivative boundary values.
    !
    do i = 1, lp1
        do j = 1, mp1
            bdzf(i, j) = (-x(i)**4)*sin(y(j))
        end do
    end do
    !
    !     note that for this example all other boundary arrays are
    !     dummy variables.
    !     we define the function boundary values in the f array.
    !
    do j = 1, mp1
        do k = 1, np1
            f(1, j, k) = 0.0_wp
            f(lp1, j, k) = sin(y(j))*cos(z(k))
        end do
    end do

    do i = 1, lp1
        do j = 1, mp1
            f(i, j, 1) = (x(i)**4)*sin(y(j))
        end do
    end do
    !
    !     we now define the values of the right side of the helmholtz
    !     equation.
    !
    do i = 2, l
        do j = 1, mp1
            do k = 2, np1
                f(i, j, k) = 4.0_wp*x(i)**2*(3.0_wp - x(i)**2)*sin(y(j))*cos(z(k))
            end do
        end do
    end do
    !
    !     call hw3crt to generate and solve the finite difference equation.
    !
    call solver%hw3crt(xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m, mbdcnd, &
        bdys, bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda, ldimf, mdimf, &
        f, pertrb, ierror)
    !
    !     compute discretization error.  the exact solution to the
    !     problem is
    !
    !        u(x, y, z) = x**4*sin(y)*cos(z)
    !
    discretization_error = 0.0_wp

    do i = 1, lp1
        do j = 1, mp1
            do k = 1, np1
                t = abs(f(i, j, k)-x(i)**4*sin(y(j))*cos(z(k)))
                discretization_error = max(t, discretization_error)
            end do
        end do
    end do

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     hw3crt *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 9.6480e-3'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        ' discretization error = ', discretization_error

end program thw3crt
