!
!     file thwscyl.f90
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
program thwscyl

    use, intrinsic :: iso_fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        FishpackSolver, &
        FishpackGrid

    ! Explicit typing only
    implicit None

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    type (FishpackSolver)        :: solver
    type (FishpackGrid)          :: grid
    integer (ip), parameter      :: m = 50
    integer (ip), parameter      :: n = 100
    integer (ip), parameter      :: idimf = m + 25
    integer (ip)                 :: mbdcnd, nbdcnd, mp1, np1, i, j, ierror
    real (wp)                    :: f(idimf, n + 5)
    real (wp), allocatable       :: bda(:), bdb(:), bdc(:), bdd(:)
    real (wp), allocatable       :: z(:), r(:)
    real (wp)                    :: a, b, c, d
    real (wp)                    :: elmbda, pertrb, x, discretization_error
    !-----------------------------------------------

    ! Set domain and boundary conditions in r
    a = 0.0_wp
    b = 1.0_wp
    mbdcnd = 6

    ! Set domain and boundary conditions in z
    c = 0.0_wp
    d = 1.0_wp
    nbdcnd = 3

    ! Set helmholtz constant
    elmbda = 0.0_wp

    ! Set auxiliary quantities
    mp1 = m + 1
    np1 = n + 1

    !
    !==> Generate and store grid points for the purpose of computing
    !    boundary data and the right side of the poisson equation.
    !
    r = grid%get_centered_grid(start=a, stop=b, num=50)
    z = grid%get_centered_grid(start=c, stop=d, num=100)

    !
    !==> Generate boundary data in z.
    !    bda is a dummy variable
    !
    allocate( bda, bdb, mold=z )
    bdb = 4.0_wp * (z**4)
    !
    !==>  Generate boundary data in r.
    !
    allocate( bdc, bdd, mold=r )
    bdc = 0.0_wp
    bdd = 4.0_wp * (r**4)

    !
    !==> Generate right side of equation.
    !
    do i = 1, mp1
        f(i,:np1) = 4.0_wp * (r(i)**2) * (z(:np1)**2) * (4.0_wp * (z(:np1)**2) + 3.0_wp * (r(i)**2))
    end do

    !
    !==> Solve system
    !
    call solver%hwscyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !==> Compute discretization error by minimizing over all a the function
    !    norm(f(i, j) - a*1 - u(r(i), z(j))).  the exact solution is
    !                u(r, z) = (r*z)**4 + arbitrary constant.
    !
    x = 0.0_wp
    do i = 1, mp1
        x = x + sum(f(i,:np1)-(r(i)*z(:np1))**4)
    end do

    x = x/(np1*mp1)
    f(:mp1,:np1) = f(:mp1,:np1) - x
    discretization_error = 0.0_wp
    do i = 1, mp1
        do j = 1, np1
            x = abs(f(i, j)-(r(i)*z(j))**4)
            discretization_error = max(x, discretization_error)
        end do
    end do

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     hwscyl *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  pertrb = 2.2674e-4'
    write( stdout, '(a)') '     discretization error = 3.7367e-4'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6)') &
        '     ierror =', ierror, ' pertrb = ', pertrb
    write( stdout, '(a,1pe15.6/)') '     discretization error = ', discretization_error

    !
    !==> Release memory
    !
    deallocate( r, z, bda, bdb, bdc, bdd )

end program thwscyl
