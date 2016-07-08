!     file thstcyl.f90
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
!     program to illustrate the use of hstcyl to solve the equation
!
!    (1/r)(d/dr)(r*du/dr) + (d/dz)(du/dz) = (2*r*z)**2*(4*z**2 + 3*r**2)
!
!     on the rectangle 0 < r < 1 , 0 < z < 1 with the
!     boundary conditions
!
!     (du/dr)(1, z) = 4*z**2  for  0 <= z <= 1
!
!     and
!
!     (du/dz)(r, 0) = 0 and (du/dz)(r, 1) = 4*r**2  for  0 <= r <= 1 .
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
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        FishpackSolver, &
        FishpackGrid

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    type (FishpackSolver)   :: solver
    type (FishpackGrid)     :: grid
    integer (ip), parameter :: m = 50
    integer (ip), parameter :: n = 52
    integer (ip), parameter :: idimf = 50 + 1
    integer (ip)            :: mbdcnd, nbdcnd, i, j, ierror
    real (wp)               :: f(idimf,n)
    real (wp), allocatable  :: r(:), z(:)
    real (wp), allocatable  :: bda(:), bdb(:), bdc(:), bdd(:)
    real (wp)               :: a, b, c, d, elmbda
    real (wp)               :: pertrb, x, discretization_error
    !-----------------------------------------------

    ! Set conditions in r
    a = 0.0_wp
    b = 1.0_wp
    mbdcnd = 6

    ! Set conditions in z
    c = 0.0_wp
    d = 1.0_wp
    nbdcnd = 3

    ! Set helmholtz constant
    elmbda = 0.0_wp
    !
    !     generate and store grid points for the purpose of computing
    !     boundary data and the right side of the poisson equation.
    !
    r = grid%get_staggered_grid(start=a, stop=b, num=m)
    z = grid%get_staggered_grid(start=c, stop=d, num=n)

    !
    !==> Generate boundary data. bda is a dummy variable
    !
    allocate( bda, bdb, mold=z )
    allocate( bdc, bdd, mold=r )

    bdb = 4.0_wp * (z**4)
    bdc = 0.0_wp
    bdd = 4.0_wp * (r**4)
    !
    !     generate right side of equation.
    !
    do i = 1, m
        f(i, :n) = 4.0_wp*(r(i)**2) * (z**2) * (4.0_wp * (z**2) + 3.0_wp * (r(i)**2) )
    end do

    !
    !==> Solve system
    !
    call solver%hstcyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error by minimizing over all a the function
    !     norm(f(i, j) - a*1 - u(r(i), z(j))).  the exact solution is
    !                u(r, z) = (r*z)**4 + arbitrary constant.
    !
    x = 0.0_wp
    do i = 1, m
        x = x + sum(f(i,:)-(r(i)*z)**4)
    end do

    x = x/(m*n)
    f(:m,:) = f(:m,:) - x

    discretization_error = 0.0_wp
    do i = 1, m
        do j = 1, n
            x = abs(f(i, j)-(r(i)*z(j))**4)
            discretization_error = max(x, discretization_error)
        end do
    end do

    !
    !==> Print earlier output from platforms with 64 bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     hstcyl *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  pertrb = -4.4311e-4'
    write( stdout, '(a)') '     discretization error = 7.5280e-5'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6)') '     ierror =', ierror, ' pertrb = ', pertrb
    write( stdout, '(a,1pe15.6/)') '     discretization error = ', discretization_error

    !
    !==> Release memory
    !
    deallocate( r, z, bda, bdb, bdc, bdd )

end program thstcyl
