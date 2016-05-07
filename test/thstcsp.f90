!
!     file thstcsp.f
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
program thstcsp

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        FishpackWorkspace, &
        hstcsp

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary: local variables
    !-----------------------------------------------
    type (FishpackWorkspace) :: workspace
    integer (ip) :: idimf, m, mbdcnd, i, n, nbdcnd, j, intl, ierror
    real (wp), dimension(47, 16) :: f
    real (wp), dimension(45) :: bda, bdb, bdc, bdd, theta
    real (wp), dimension(15) :: r
    real (wp), dimension(45) :: cost
    real (wp) :: a, b, dt, c, d, dr, elmbda, pertrb, discretization_error, z
    !-----------------------------------------------
    !
    !     note that from dimension statement we get that idimf = 47
    !
    idimf = 47
    a = 0.0_wp
    b = acos(-1.0_wp)
    !
    !
    m = 45
    mbdcnd = 9
    dt = (b - a)/real(m)
    !
    !     define grid points theta(i) and cos(theta(i))
    !
    do i = 1, m
        theta(i) = a + (real(i) - 0.5)*dt
        cost(i) = cos(theta(i))
    end do
    c = 0.
    d = 1.
    n = 15
    nbdcnd = 5
    dr = (d - c)/real(n)
    !
    !     define grid points r(j)
    !
    do j = 1, n
        r(j) = c + (real(j) - 0.5)*dr
    end do
    !
    !     define boundary array bdd.  bda, bdb, and bdc are dummy
    !     variables in this example.
    !
    bdd(:m) = cost(:m)**4
    elmbda = 0.
    !
    !     define right side f
    !
    do i = 1, m
        f(i, :n) = 12.*(r(:n)*cost(i))**2
    end do

    ! Initialize
    intl = 0

    call hstcsp(intl, a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
        bdd, elmbda, f, idimf, pertrb, ierror, workspace)
    !
    !     compute discretization error.  the exact solution is
    !
    !     u(theta, r) = (r*cos(theta))**4
    !
    discretization_error = 0.
    do i = 1, m
        do j = 1, n
            z = abs(f(i, j)-(r(j)*cost(i))**4)
            discretization_error = max(z, discretization_error)
        end do
    end do


    !     Print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    !     in this example (contrast with blktri and sepeli) the extra precision
    !     does not reduce the discretization error
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     hstcsp *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  discretization error = 5.5843e-3'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') &
        '     ierror =', ierror, ' discretization error = ', discretization_error

    !     release memory allocated by hstcsp (intl=0 call)
    call workspace%destroy()

end program thstcsp
