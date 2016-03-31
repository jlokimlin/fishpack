!
!     file thwscrt.f
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
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
program thwscrt

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use modern_fishpack_library, only: &
        hwscrt

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    integer :: idimf, m, mbdcnd, n, nbdcnd, mp1, np1, i, j, ierror
    real (wp), dimension(45, 82) :: f
    real (wp), dimension(81) :: bdb, bda, bdc, bdd, y
    real (wp), dimension(41) :: x
    real::a, b, c, d, elmbda, pi, dum, piby2, pisq, pertrb, discretization_error, z
    !-----------------------------------------------
    !
    !     from dimension statement we get value of idimf.
    !
    idimf = 45
    a = 0.
    b = 2.
    m = 40
    mbdcnd = 2
    c = -1.
    d = 3.
    n = 80
    nbdcnd = 0
    elmbda = -4.
    !
    !     auxiliary quantities.
    !
    pi = acos( -1.0 )
    piby2 = pi/2.
    pisq = pi**2
    mp1 = m + 1
    np1 = n + 1
    !
    !     generate and store grid points for the purpose of computing
    !     boundary data and the right side of the helmholtz equation.
    !
    do i = 1, mp1
        x(i) = real(i - 1)/20.
    end do
    do j = 1, np1
        y(j) = (-1.) + real(j - 1)/20.
    end do
    !
    !     generate boundary data.
    !
    do j = 1, np1
        bdb(j) = 4.*cos((y(j)+1.)*piby2)
    end do
    !
    !     bda, bdc, and bdd are dummy variables.
    !
    f(1, :np1) = 0.
    !
    !     generate right side of equation.
    !
    do i = 2, mp1
        do j = 1, np1
            f(i, j) = (2. - (4. + pisq/4.)*x(i)**2)*cos((y(j)+1.)*piby2)
        end do
    end do

    call hwscrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, bdd, &
        elmbda, f, idimf, pertrb, ierror)
    !
    !     compute discretization error.  the exact solution is
    !                u(x, y) = x**2*cos((y+1)*piby2)
    !
    discretization_error = 0.0_wp
    do i = 1, mp1
        do j = 1, np1
            z = abs(f(i, j)-(x(i)**2) * cos((y(j)+1.0_wp) * piby2))
            discretization_error = max(z, discretization_error)
        end do
    end do
    !     print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     hwscrt *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0,  discretization error = 5.36508e-4'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3,A,1pe15.6)') &
        '     ierror =', ierror, ' discretization error = ', discretization_error

end program thwscrt
