!
!     file tsepx4.f90
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright(c) 2005 by UCAR                   *
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
!     *                Boulder, Colorado (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
program tsepx4

    use, intrinsic :: iso_fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        FishpackSolver

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    type (FishpackSolver)   :: solver
    integer (ip), parameter :: m = 32
    integer (ip), parameter :: n = 32
    integer (ip)            :: nx, ny, i, j, mbdcnd, nbdcnd, idmn, iorder, ierror
    real (wp), dimension(33, 33) :: usol, grhs
    real (wp), dimension(33) :: bda, bdb
    real (wp) :: a, b, c, d, dlx, dly, x, af, bf, cf, y
    real (wp) :: alpha, beta, dum(1), pertrb, err, err2, err4
    !--------------------------------------------------------------------------

    !
    !     define arithmetic functions giving exact solution
    !
    !
    !     set limits on region
    !
    a = 0.0_wp
    b = 1.0_wp
    c = 0.0_wp
    d = 1.0_wp
    !
    !     set grid size
    !
    dlx =(b - a)/m
    dly =(d - c)/n
    nx = m + 1
    ny = n + 1
    do i = 1, nx
        x = a + real(i - 1, kind=wp) * dlx
        !
        !     set specified boundary conditions at y=c, d
        !
        usol(i, 1) = ue(x, c)
        usol(i, ny) = ue(x, d)
        call get_coefficients_in_x_direction(x, af, bf, cf)
        do j = 1, ny
            y = c + real(j - 1, kind=wp)*dly
            !
            !     set right hand side
            !
            grhs(i, j)=af*uxxe(x, y)+bf*uxe(x, y)+cf*ue(x, y)+uyye(x, y)
        end do
    end do
    !
    !     set mixed boundary conditions at x=a, b
    !
    alpha = 1.0_wp
    beta = 1.0_wp
    do j = 1, ny
        y = c + real(j - 1, kind=wp) * dly
        bda(j) = uxe(a, y) + alpha * ue(a, y)
        bdb(j) = uxe(b, y) + beta * ue(b, y)
    end do
    !
    !     set boundary swithces
    !
    mbdcnd = 3
    nbdcnd = 1
    !
    !     set first dimension of usol, grhs and work space length
    !
    idmn = 33
    !
    !     obtain second order approximation
    !
    iorder = 2

    call solver%sepx4(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, d, &
        n, nbdcnd, dum, dum, get_coefficients_in_x_direction, grhs, usol, idmn, pertrb, ierror)
    !
    !     compute second order discretization error(relative)
    !     also reset specified boundaries and right hand side.
    !
    err = 0.0_wp
    do i = 1, nx
        x = a + real(i - 1, kind=wp) * dlx
        usol(i, 1) = ue(x, c)
        usol(i, ny) = ue(x, d)
        call get_coefficients_in_x_direction(x, af, bf, cf)
        do j = 1, ny
            y = c + real(j - 1, kind=wp)*dly
            err = max(err, abs((usol(i, j)-ue(x, y))/ue(x, y)))
            !
            !     reset right hand side in grhs for fourth order approximation call
            !
            grhs(i, j)=af*uxxe(x, y)+bf*uxe(x, y)+cf*ue(x, y)+uyye(x, y)
        end do
    end do
    err2 = err
    !
    !     obtain fourth order approximation
    !
    iorder = 4
    call solver%sepx4(iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, c, d, &
        n, nbdcnd, dum, dum, get_coefficients_in_x_direction, &
        grhs, usol, idmn, pertrb, ierror)
    !
    !     compute fourth order discretization error(relative)
    !
    err = 0.0_wp
    do j = 1, ny
        y = c + real(j - 1, kind=wp)*dly
        do i = 1, nx
            x = a + real(i - 1, kind=wp)*dlx
            err = max(err, abs((usol(i, j)-ue(x, y))/ue(x, y)))
        end do
    end do
    err4 = err

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     sepx4 *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0'
    write( stdout, '(a)') '     Second Order discretization error = 1.5985e-4'
    write( stdout, '(a)') '     Fourth Order discretization error = 1.8575e-6'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3)') '     ierror =', ierror
    write( stdout, '(a,1pe15.6)') '     Second Order discretization error =', err2
    write( stdout, '(a,1pe15.6/)') '     Fourth Order discretization error =', err4


contains


    pure function ue(s, t) result (return_value)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real (wp), intent (in) :: s
        real (wp), intent (in) :: t
        real (wp)          :: return_value
        !--------------------------------------------------------------

        return_value =(s*t)**3 + 1.0_wp

    end function ue


    pure function uxe(s, t) result (return_value)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real (wp), intent (in) :: s
        real (wp), intent (in) :: t
        real (wp)          :: return_value
        !--------------------------------------------------------------

        return_value = 3.0_wp *(s**2)*(t**3)

    end function uxe


    pure function uxxe(s, t) result (return_value)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real (wp), intent (in) :: s
        real (wp), intent (in) :: t
        real (wp)          :: return_value
        !--------------------------------------------------------------

        return_value = 6.0_wp * s *(t**3)

    end function uxxe


    pure function uyye(s, t) result (return_value)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real (wp), intent (in) :: s
        real (wp), intent (in) :: t
        real (wp)              :: return_value
        !--------------------------------------------------------------

        return_value = 6.0_wp *(s**3) * t

    end function uyye


    pure subroutine get_coefficients_in_x_direction(x, af, bf, cf)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real (wp), intent (in)  :: x
        real (wp), intent (out) :: af
        real (wp), intent (out) :: bf
        real (wp), intent (out) :: cf
        !--------------------------------------------------------------

        !     set coefficients in the x-direction.
        !
        af =(x + 1.0_wp)**2
        bf = 2.0_wp * (x + 1.0_wp)
        cf = -x

    end subroutine get_coefficients_in_x_direction


end program tsepx4
