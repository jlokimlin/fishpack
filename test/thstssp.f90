!
!     file thstssp.f90
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
program thstssp

    use, intrinsic :: iso_fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, & ! machine precision PI
        TWO_PI, &
        FishpackGrid, &
        FishpackSolver

    ! Explicit typing only
    implicit None

    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    type (FishpackSolver)   :: solver
    type (FishpackGrid)     :: grid
    integer (ip), parameter :: m = 18
    integer (ip), parameter :: n = 72
    integer (ip), parameter :: idimf = 18
    integer (ip)            :: mbdcnd, nbdcnd, i, j, ierror
    real (wp)               :: f(m,n)
    real (wp), allocatable  :: bda(:), bdb(:), bdc(:), bdd(:)
    real (wp), allocatable  :: sint(:), sinp(:)
    real (wp)               :: a, b, c, d, elmbda
    real (wp)               :: pertrb, discretization_error, local_error
    !-----------------------------------------------

    a = 0.0_wp
    b = PI/2
    mbdcnd = 6

    c = 0.0_wp
    d = TWO_PI
    nbdcnd = 0

    ! Set helmholtz constant
    elmbda = 0.0_wp

    !
    !==> generate sines for use in subsequent computations
    !
    sint = grid%get_staggered_grid(start=a, stop=b, num=m)
    sint = sin(sint)

    sinp = grid%get_staggered_grid(start=c, stop=d, num=n)
    sinp = sin(sinp)
    !
    !     compute right side of equation and store in f
    !
    do j = 1, n
        f(:, j) = 2.0_wp - 6.0_wp*(sint*sinp(j))**2
    end do
    !
    !     store derivative data at the equator
    !
    allocate(bda, bdb, bdc, bdd, mold=sinp)
    bdb = 0.0_wp
    !
    !     bda, bdc, and bdd are dummy variables.
    !
    call solver%hstssp(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error. since problem is singular, the
    !     solution must be normalized.
    !
    discretization_error = 0.0_wp
    do j = 1, n
        do i = 1, m
            local_error = abs(f(i, j)-(sint(i)*sinp(j))**2-f(1, 1))
            discretization_error = max(local_error, discretization_error)
        end do
    end do

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     hstssp *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  pertrb = 6.35830e-4'
    write( stdout, '(a)') '     discretization error = 3.37523e-3'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6)') '     ierror =', ierror, ' pertrb = ', pertrb
    write( stdout, '(a,1pe15.6/)') '     discretization error = ', discretization_error

    !
    !==> Release memory
    !
    deallocate( sint, sinp )
    deallocate( bda, bdb, bdc, bdd )

end program thstssp
