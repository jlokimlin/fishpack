!     file thstcrt.f90
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
program thstcrt

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
    ! Dictionary
    !-----------------------------------------------
    type (FishpackSolver)   :: solver
    type (FishpackGrid)     :: grid
    integer (ip), parameter :: m = 48
    integer (ip), parameter :: n = 53
    integer (ip), parameter :: idimf = m + 2
    integer (ip)            :: mbdcnd, nbdcnd, i, j, ierror
    real (wp)               :: f(idimf, n)
    real (wp), allocatable  :: bda(:), bdb(:), bdc(:), bdd(:)
    real (wp), allocatable  :: x(:), y(:)
    real (wp)               :: a, b, c, d
    real (wp)               :: elmbda, local_error, pertrb, discretization_error
    real (wp), parameter    :: PI = acos(-1.0_wp)
    real (wp), parameter    :: PI2 = PI**2
    !-----------------------------------------------

    ! Set conditions in x
    a = 1.0_wp
    b = 3.0_wp
    mbdcnd = 2

    ! Set conditions in y
    c = -1.0_wp
    d = 1.0_wp
    nbdcnd = 0

    ! Set helmholtz constant
    elmbda = -2.0_wp

    !
    !==> Generate and store grid points for computation of boundary data
    !    and the right side of the helmholtz equation.
    !
    x = grid%get_staggered_grid(start=a, stop=b, num=m)
    y = grid%get_staggered_grid(start=c, stop=d, num=n)

    !
    !==> Generate boundary data.
    !
    allocate( bda, bdb, mold=y )

    bda = 0.0_wp
    bdb = -PI*cos(PI*y)

    ! bdc and bdd are dummy arguments in this example.
    allocate( bdc, bdd, mold=x )


    associate( CONST => -2.0_wp*(PI2 + 1.0_wp) )
        do i = 1, m
            !
            !==> Generate right side of equation.
            !
            f(i,:) = CONST*sin(PI*x(i))*cos(PI*y(:))

        end do
    end associate


    !
    !==> Solve system
    !
    call solver%hstcrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, &
        bdc, bdd, elmbda, f, idimf, pertrb, ierror)

    discretization_error = 0.0_wp
    do i = 1, m
        do j = 1, n
            !
            !==> Compute discretization error. The exact solution is
            !
            !    u(x, y) = sin(pi*x)*cos(pi*y) .
            !
            local_error = abs(f(i, j)-sin(PI*x(i))*cos(PI*y(j)))
            discretization_error = max(local_error, discretization_error)
        end do
    end do

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     hstcrt *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 1.2600e-3'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        ' discretization error = ', discretization_error

    !
    !==> Release memory
    !
    deallocate( x, y, bda, bdb, bdc, bdd )

end program thstcrt
