!
!     file thstcsp.f90
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
program thstcsp

    use, intrinsic :: iso_fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        FishpackSolver, &
        FishpackGrid, &
        FishpackWorkspace

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    type (FishpackSolver)    :: solver
    type (FishpackGrid)      :: grid
    type (FishpackWorkspace) :: workspace
    integer (ip), parameter  :: m = 45
    integer (ip), parameter  :: n = 15
    integer (ip), parameter  :: idimf = m + 2
    integer (ip)             :: mbdcnd, i, nbdcnd, j, intl, ierror
    real (wp)                :: f(idimf, n + 1)
    real (wp), allocatable   :: theta(:), cost(:), r(:)
    real (wp), allocatable   :: bda(:), bdb(:), bdc(:), bdd(:)
    real (wp), parameter     :: PI = acos(-1.0_wp)
    real (wp)                :: a, b, c, d, elmbda
    real (wp)                :: pertrb, discretization_error, local_error
    !-----------------------------------------------

    ! Set conditions in theta
    a = 0.0_wp
    b = PI
    mbdcnd = 9

    ! Set conditions in r
    c = 0.0_wp
    d = 1.0_wp
    nbdcnd = 5

    !
    !==> Define grid points theta(i) and cos(theta(i))
    !
    theta = grid%get_staggered_grid(start=a, stop=b, num=m)
    ! allocate( cost, source=cos(theta ) ! *** bug in GCC 6.1
    allocate( cost, mold=theta )
    cost = cos(theta)
    !
    !==> Define grid points r(j)
    !
    r = grid%get_staggered_grid(start=c, stop=d, num=n)

    !
    !==> Define boundary array bdd.
    !    bda, bdb, and bdc are dummy
    !    variables in this example.
    !
    allocate( bda, bdb, bdc, bdd, mold=theta )

    bdd(:m) = cost(:m)**4
    elmbda = 0.0_wp
    !
    !     define right side f
    !
    do i = 1, m
        f(i,:n) = 12.0_wp*(r(:n)*cost(i))**2
    end do

    ! Initialize
    intl = 0

    call solver%hstcsp(intl, a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, workspace)

    discretization_error = 0.0_wp
    do i = 1, m
        do j = 1, n
            !
            !==> Compute discretization error.
            !    The exact solution is
            !
            !     u(theta, r) = (r*cos(theta))**4
            !
            local_error = abs(f(i,j)-(r(j)*cost(i))**4)
            discretization_error = max(local_error, discretization_error)
        end do
    end do

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer.
    !    In this example (contrast with blktri and sepeli) the extra precision
    !    does not reduce the discretization error
    !
    write( stdout, '(/a)') '     hstcsp *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0,  discretization error = 5.5843e-3'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3,a,1pe15.6/)') &
        '     ierror =', ierror, &
        ' discretization error = ', discretization_error

    !
    !==> Release memory
    !
    call workspace%destroy()
    deallocate( theta, cost, r, bda, bdb, bdc, bdd )

end program thstcsp
