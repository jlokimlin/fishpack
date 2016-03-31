!
!     file thwscyl.f
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  version 1.1                    *
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
program thwscyl

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use modern_fishpack_library, only: &
        hwscyl

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer (ip) :: idimf, m, mbdcnd, n, nbdcnd, mp1, np1, i, j, ierror
    real , dimension(75, 105) :: f
    real , dimension(101) :: bda, bdb
    real , dimension(51) :: bdc, bdd, r
    real , dimension(101) :: z
    real :: a, b, c, d, elmbda, pertrb, x, discretization_error
    !-----------------------------------------------
    !
    !     from dimension statement we get value of idimf.
    !
    idimf = 75
    a = 0.
    b = 1.
    m = 50
    mbdcnd = 6
    c = 0.
    d = 1.
    n = 100
    nbdcnd = 3
    elmbda = 0.
    !
    !     auxiliary quantities.
    !
    mp1 = m + 1
    np1 = n + 1
    !
    !     generate and store grid points for the purpose of computing
    !     boundary data and the right side of the poisson equation.
    !
    do i = 1, mp1
        r(i) = real(i - 1)/50.
    end do
    do j = 1, np1
        z(j) = real(j - 1)/100.
    end do
    !
    !     generate boundary data.
    !
    bdb(:np1) = 4.*z(:np1)**4
    !
    !     generate boundary data.
    !
    bdc(:mp1) = 0.
    bdd(:mp1) = 4.*r(:mp1)**4
    !
    !     bda is a dummy variable.
    !
    !
    !     generate right side of equation.
    !
    do i = 1, mp1
        f(i, :np1) = 4.*r(i)**2*z(:np1)**2*(4.*z(:np1)**2+3.*r(i)**2)
    end do

    ! Solve system
    call hwscyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error by minimizing over all a the function
    !     norm(f(i, j) - a*1 - u(r(i), z(j))).  the exact solution is
    !                u(r, z) = (r*z)**4 + arbitrary constant.
    !
    x = 0.
    do i = 1, mp1
        x = x + sum(f(i, :np1)-(r(i)*z(:np1))**4)
    end do
    x = x/real(np1*mp1)
    f(:mp1, :np1) = f(:mp1, :np1) - x
    discretization_error = 0.
    do i = 1, mp1
        do j = 1, np1
            x = abs(f(i, j)-(r(i)*z(j))**4)
            discretization_error = max(x, discretization_error)
        end do
    end do


    !     Print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     hwscyl *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  pertrb = 2.2674e-4'
    write( stdout, '(A)') '     discretization error = 3.7367e-4'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') &
        '     ierror =', ierror, ' pertrb = ', pertrb
    write( stdout, '(A,1pe15.6)') '     discretization error = ', discretization_error

end program thwscyl
