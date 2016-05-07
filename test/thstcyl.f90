!     file thstcyl.f
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
!
!     program to illustrate the use of hstcyl to solve the equation
!
!    (1/r)(d/dr)(r*du/dr) + (d/dz)(du/dz) = (2*r*z)**2*(4*z**2 + 3*r**2)
!
!     on the rectangle 0 .lt. r .lt. 1 , 0 .lt. z .lt. 1 with the
!     boundary conditions
!
!     (du/dr)(1, z) = 4*z**2  for  0 .le. z .le. 1
!
!     and
!
!     (du/dz)(r, 0) = 0 and (du/dz)(r, 1) = 4*r**2  for  0 .le. r .le. 1 .
!
!     the solution to this problem is not unique.  it is a
!     one-parameter family of solutions given by
!
!            u(r, z) = (r*z)**4 + arbitrary constant .
!
!     the r-interval will contain 50 unknowns and the z-interval will
!     contain 52 unknowns.
!
program thstcyl

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        FishpackSolver

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    type (FishpackSolver)        :: solver
    integer (ip)                 :: idimf, m, mbdcnd, n, nbdcnd, i, j, ierror
    real (wp), dimension(51, 52) :: f
    real (wp), dimension(52)     :: bda, bdb
    real (wp), dimension(50)     :: bdc, bdd, r
    real (wp), dimension(52)     :: z
    real (wp)                    :: a, b, c, d, elmbda, pertrb, x, discretization_error
    !-----------------------------------------------

    !     from dimension statement we get value of idimf.
    !
    idimf = 51
    a = 0.0_wp
    b = 1.0_wp
    m = 50
    mbdcnd = 6
    c = 0.0_wp
    d = 1.0_wp
    n = 52
    nbdcnd = 3
    elmbda = 0.0_wp
    !
    !     generate and store grid points for the purpose of computing
    !     boundary data and the right side of the poisson equation.
    !
    do i = 1, m
        r(i) = (real(i, kind=wp) - 0.5_wp)/50
    end do

    do j = 1, n
        z(j) = (real(j, kind=wp) - 0.5_wp)/52
    end do
    !
    !     generate boundary data.
    !
    bdb(:n) = 4.0_wp * ( z(:n)**4 )
    !
    !     generate boundary data.
    !
    bdc(:m) = 0.0_wp
    bdd(:m) = 4.0_wp * ( r(:m)**4 )
    !
    !     bda is a dummy variable.
    !
    !     generate right side of equation.
    !
    do i = 1, m
        f(i, :n) = 4.0_wp*(r(i)**2) * (z(:n)**2) * (4.0_wp * (z(:n)**2) &
            + 3.0_wp * (r(i)**2) )
    end do

    ! Solve system
    call solver%hstcyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error by minimizing over all a the function
    !     norm(f(i, j) - a*1 - u(r(i), z(j))).  the exact solution is
    !                u(r, z) = (r*z)**4 + arbitrary constant.
    !
    x = 0.0_wp
    do i = 1, m
        x = x + sum(f(i, :n)-(r(i)*z(:n))**4)
    end do

    x = x/(m*n)
    f(:m, :n) = f(:m, :n) - x

    discretization_error = 0.0_wp
    do i = 1, m
        do j = 1, n
            x = abs(f(i, j)-(r(i)*z(j))**4)
            discretization_error = max(x, discretization_error)
        end do
    end do

    !     Print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     hstcyl *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  pertrb = -4.4311e-4'
    write( stdout, '(a)') '     discretization error = 7.5280e-5'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') '     ierror =', ierror, ' pertrb = ', pertrb
    write( stdout, '(A,1pe15.6)') '     discretization error = ', discretization_error

end program thstcyl
