!
!     file tsepeli.f90
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
program tsepeli

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use fishpack_library, only: &
        wp, &
        ip, &
        FishpackSolver, &
        FishpackWorkspace

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    type(FishpackSolver)    :: solver
    type(FishpackWorkspace) :: workspace
    integer(ip), parameter  :: m = 32
    integer(ip), parameter  :: n = 32
    integer(ip), parameter  :: idmn = m + 1
    integer(ip)             :: nx, ny, i, j, mbdcnd, nbdcnd, intl, iorder, ierror
    real(wp), allocatable   :: usol(:,:), grhs(:,:), bda(:), bdb(:)
    real(wp)                :: a, b, c, d, dlx, dly
    real(wp)                :: x, af, bf, cf, y, df, ef, ff, alpha
    real(wp)                :: beta, dummy_variable(1), pertrb
    real(wp)                :: local_error, order2_error, order4_error
    real(wp), parameter     :: ZERO = 0.0_wp
    real(wp), parameter     :: ONE = 1.0_wp, TWO = 2.0_wp
    !------------------------------------------------------------------

    !
    !==> Allocate memory
    !
    allocate( usol(idmn,idmn), grhs(idmn,idmn), bda(idmn), bdb(idmn) )

    !     define arithmetic functions giving exact solution
    !
    !
    !     set limits on region
    !
    a = ZERO
    b = ONE
    c = ZERO
    d = ONE
    !
    !     set grid size
    !
    dlx = (b - a)/m
    dly = (d - c)/n
    nx = m + 1
    ny = n + 1
    do i = 1, nx
        x = a + real(i - 1, kind=wp)*dlx
        !
        !     set specified boundary conditions at y=c, d
        !
        usol(i, 1) = ue(x, c)
        usol(i, ny) = ue(x, d)
        call get_coefficients_in_x_direction (x, af, bf, cf)
        do j = 1, ny
            y = c + real(j - 1, kind=wp)*dly
            call get_coefficients_in_y_direction (y, df, ef, ff)
            !
            !     set right hand side
            !
            grhs(i, j) = af*uxxe(x, y) + bf*uxe(x, y) + cf*ue(x, y) + df* &
                uyye(x, y) + ef*uye(x, y) + ff*ue(x, y)
        end do
    end do
    !
    !     set mixed boundary conditions at x=a, b
    !
    alpha = ONE
    beta = ONE
    do j = 1, ny
        y = c + real(j - 1, kind=wp)*dly
        bda(j) = uxe(a, y) + alpha*ue(a, y)
        bdb(j) = uxe(b, y) + beta*ue(b, y)
    end do
    !
    !     set boundary swithces
    !
    mbdcnd = 3
    nbdcnd = 1
    !

    !     set for initialization of sepeli
    intl = 0
    !
    !     obtain second order approximation
    !
    iorder = 2

    call solver%sepeli(intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta, &
        c, d, n, nbdcnd, dummy_variable(1:1), dummy_variable(1), dummy_variable(1:1), dummy_variable(1), &
        get_coefficients_in_x_direction, get_coefficients_in_y_direction, grhs, usol, &
        idmn, workspace, pertrb, ierror)

    local_error = ZERO

    do i = 1, nx
        x = a + real(i - 1)*dlx
        do j = 1, ny
            y = c + real(j - 1)*dly
            local_error = max(local_error, abs((usol(i, j)-ue(x, y))/ue(x, y)))
        end do
    end do

    order2_error = local_error
    !
    !     obtain fourth order approximation
    !
    iorder = 4
    !
    !     non-initial call
    !
    intl = 1
    call solver%sepeli(intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta &
        , c, d, n, nbdcnd, dummy_variable(1:1), dummy_variable(1), dummy_variable(1:1), dummy_variable(1), &
        get_coefficients_in_x_direction, get_coefficients_in_y_direction, grhs, usol, &
        idmn, workspace, pertrb, ierror)
    !
    !     compute discretization error
    !
    local_error = ZERO
    do j = 1, ny
        y = c + real(j - 1, kind=wp)*dly
        do i = 1, nx
            x = a + real(i - 1, kind=wp)*dlx
            local_error = max(local_error, abs((usol(i, j)-ue(x, y))/ue(x, y)))
        end do
    end do
    order4_error = local_error

    !
    !==> Print earlier output from platforms with 64-bit floating point
    !    arithmetic followed by the output from this computer
    !
    write( stdout, '(/a)') '     sepeli *** TEST RUN *** '
    write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(a)') '     ierror = 0'
    write( stdout, '(a)') '     Second Order discretization error = 9.7891e-5'
    write( stdout, '(a)') '     Fourth Order discretization error = 1.4735e-6'
    write( stdout, '(a)') '     The output from your computer is: '
    write( stdout, '(a,i3)')  '     ierror =', ierror
    write( stdout, '(a,1pe15.6)') '     Second Order discretization error =', order2_error
    write( stdout, '(a,1pe15.6/)')  '     Fourth Order discretization error =', order4_error

    !
    !==> Release memory
    !
    call workspace%destroy()
    deallocate( usol, grhs, bda, bdb )


contains


    pure function ue (s, t) result (return_value)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real(wp), intent(in) :: s
        real(wp), intent(in) :: t
        real(wp)              :: return_value
        !--------------------------------------------------------------

        return_value = (s * t)**3 + ONE

    end function ue


    pure function uxe (s, t) result( return_value )
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real(wp), intent(in) :: s
        real(wp), intent(in) :: t
        real(wp)              :: return_value
        !--------------------------------------------------------------

        return_value = 3.0_wp * (s**2) * (t**3)

    end function uxe


    pure function uxxe(s, t) result( return_value )
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real(wp), intent(in) :: s
        real(wp), intent(in) :: t
        real(wp)              :: return_value
        !--------------------------------------------------------------

        return_value = 6.0_wp * s * (t**3)

    end function uxxe


    pure function uye(s, t) result( return_value )
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real(wp), intent(in) :: s
        real(wp), intent(in) :: t
        real(wp)              :: return_value

        return_value = 3.0_wp * (s**3) * (t**2)

    end function uye


    pure function uyye (s, t) result( return_value )
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real(wp), intent(in) :: s
        real(wp), intent(in) :: t
        real(wp)              :: return_value

        return_value = 6.0_wp * (s**3) * t

    end function uyye


    pure subroutine get_coefficients_in_x_direction(x, af, bf, cf)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real(wp), intent(in)  :: x
        real(wp), intent(out) :: af
        real(wp), intent(out) :: bf
        real(wp), intent(out) :: cf
        !-----------------------------------------------
        !
        !     set coefficients in the x-direction.
        !
        af = (x + ONE)**2
        bf = TWO * (x + ONE)
        cf = -x

    end subroutine get_coefficients_in_x_direction


    pure subroutine get_coefficients_in_y_direction(y, df, ef, ff)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        real(wp), intent(in)  :: y
        real(wp), intent(out) :: df
        real(wp), intent(out) :: ef
        real(wp), intent(out) :: ff
        !-----------------------------------------------
        !
        !     set coefficients in y direction
        !
        df = exp(y)
        ef = ZERO
        ff = -y

    end subroutine get_coefficients_in_y_direction


end program tsepeli
