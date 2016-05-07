!
!     file tpois3d.f
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
!     *                 A Package of Fortran 77 and 90                *
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
program tpois3d

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        pois3d

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer :: ldimf, mdimf, lperod, l, mperod, m, nperod, n, i, j, k, ierror
    real (wp), dimension(32, 33, 10) :: f
    real (wp), dimension(10) :: a, b, c
    real (wp), dimension(30) :: x, y
    real (wp), dimension(10) :: z
    real (wp) :: pi, dx, c1, dy, c2, dz, dzsq, t, discretization_error
    !-----------------------------------------------
    !
    !     from the dimension statement we get that ldimf = 32, mdimf = 33,
    !
    ldimf = 32
    mdimf = 33
    pi = acos( -1.0 )
    lperod = 0
    l = 30
    dx = 2.*pi/real(l)
    c1 = 1./dx**2
    mperod = 0
    m = 30
    dy = 2.*pi/real(m)
    c2 = 1./dy**2
    nperod = 1
    n = 10
    dz = 1./real(n)
    dzsq = 1./dz**2
    !
    !     generate grid points for later use.
    !
    do i = 1, l
        x(i) = (-pi) + real(i - 1)*dx
    end do
    do j = 1, m
        y(j) = (-pi) + real(j - 1)*dy
    end do
    !
    !     generate coefficients
    !
    a(1) = 0.
    b(1) = -2.*dzsq
    c(1) = -b(1)
    z(1) = 0.
    do k = 2, n
        z(k) = real(k - 1)*dz
        t = 1. + z(k)
        a(k) = t**2*dzsq + t/dz
        b(k) = -2.*t**2*dzsq
        c(k) = t**2*dzsq - t/dz
    end do
    !
    !     generate right side of equation
    !
    do i = 1, l
        do j = 1, m
            do k = 2, n
                f(i, j, k) = 2.*sin(x(i))*sin(y(j))*(1. + z(k))**4
            end do
        end do
    end do
    do i = 1, l
        do j = 1, l
            f(i, j, 1) = (10. + 8./dz)*sin(x(i))*sin(y(j))
            f(i, j, n) = f(i, j, n) - c(n)*16.*sin(x(i))*sin(y(j))
        end do
    end do
    c(n) = 0.
    !
    !     call pois3d to solve equations.
    !
    call pois3d(lperod, l, c1, mperod, m, c2, nperod, n, a, b, c, &
        ldimf, mdimf, f, ierror)
    !
    !     compute discretization error.  the exact solution is
    !
    !              u(x, y, z) = sin(x)*sin(y)*(1+z)**4
    !
    discretization_error = 0.
    do i = 1, l
        do j = 1, m
            do k = 1, n
                t = abs(f(i, j, k)-sin(x(i))*sin(y(j))*(1.+z(k))**4)
                discretization_error = max(t, discretization_error)
            end do
        end do
    end do
    !     Print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     pois3d *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  discretization error = 2.93277e-2'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') '     ierror =', ierror, ' discretization error = ', &
        discretization_error

end program tpois3d
