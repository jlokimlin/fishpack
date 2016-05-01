!
!     file tsepeli.f
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
program tsepeli

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT

    use modern_fishpack_library, only: &
        FishpackWorkspace, &
        sepeli

    ! Explicit typing only
    implicit none

    !-----------------------------------------------
    ! Dictionary
    !-----------------------------------------------
    type (FishpackWorkspace) :: workspace
    integer (ip)             :: m, n, nx, ny, i, j, mbdcnd, nbdcnd, idmn, intl, iorder, ierror
    real (wp), allocatable   :: usol(:,:), grhs(:,:)
    real (wp), allocatable   :: bda(:), bdb(:)
    real (wp)                :: a, b, c, d, dlx, dly
    real (wp)                :: x, af, bf, cf, y, df, ef, ff, alpha
    real (wp)                :: beta, dum(1), pertrb, err, err2, err4
    !-----------------------------------------------

    ! allocate memory
    allocate( usol(33,33), grhs(33,33), bda(33), bdb(33) )

    !     define arithmetic functions giving exact solution
    !
    !
    !     set limits on region
    !
    a = 0.0
    b = 1.0
    c = 0.0
    d = 1.0
    !
    !     set grid size
    !
    m = 32
    n = 32
    dlx = (b - a)/real(m)
    dly = (d - c)/real(n)
    nx = m + 1
    ny = n + 1
    do i = 1, nx
        x = a + real(i - 1)*dlx
        !
        !     set specified boundary conditions at y=c, d
        !
        usol(i, 1) = ue(x, c)
        usol(i, ny) = ue(x, d)
        call get_coefficients_in_x_direction (x, af, bf, cf)
        do j = 1, ny
            y = c + real(j - 1)*dly
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
    alpha = 1.0
    beta = 1.0
    do j = 1, ny
        y = c + real(j - 1)*dly
        bda(j) = uxe(a, y) + alpha*ue(a, y)
        bdb(j) = uxe(b, y) + beta*ue(b, y)
    end do
    !
    !     set boundary swithces
    !
    mbdcnd = 3
    nbdcnd = 1
    !
    !     set first dimension of usol, grhs
    !
    idmn = 33
    !     set for initialization of sepeli
    intl = 0
    !
    !     obtain second order approximation
    !
    iorder = 2
    call sepeli(intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta &
        , c, d, n, nbdcnd, dum(1:1), dum(1), dum(1:1), dum(1), &
        get_coefficients_in_x_direction, get_coefficients_in_y_direction, grhs, usol, &
        idmn, workspace, pertrb, ierror)
    err = 0.0
    do i = 1, nx
        x = a + real(i - 1)*dlx
        do j = 1, ny
            y = c + real(j - 1)*dly
            err = max(err, abs((usol(i, j)-ue(x, y))/ue(x, y)))
        end do
    end do
    err2 = err
    !
    !     obtain fourth order approximation
    !
    iorder = 4
    !
    !     non-initial call
    !
    intl = 1
    call sepeli (intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta &
        , c, d, n, nbdcnd, dum(1:1), dum(1), dum(1:1), dum(1), &
        get_coefficients_in_x_direction, get_coefficients_in_y_direction, grhs, usol, &
        idmn, workspace, pertrb, ierror)
    !
    !     compute discretization error
    !
    err = 0.0
    do j = 1, ny
        y = c + real(j - 1)*dly
        do i = 1, nx
            x = a + real(i - 1)*dlx
            err = max(err, abs((usol(i, j)-ue(x, y))/ue(x, y)))
        end do
    end do
    err4 = err
    !     Print earlier output from platforms with 32 and 64 bit floating point
    !     arithemtic followed by the output from this computer
    write( stdout, '(A)') ''
    write( stdout, '(A)') '     sepeli *** TEST RUN *** '
    write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
    write( stdout, '(A)') '     ierror = 0'
    write( stdout, '(A)') '     Second Order discretization error = 9.7891e-5'
    write( stdout, '(A)') '     Fourth Order discretization error = 1.4735e-6'
    write( stdout, '(A)') '     The output from your computer is: '
    write( stdout, '(A,I3)')  '     ierror =', ierror
    write( stdout, '(A,1pe15.6)') '     Second Order discretization error =', err2
    write( stdout, '(A,1pe15.6)')  '     Fourth Order discretization error =', err4

    ! release dynamically allocated real and complex work space
    call workspace%destroy()

    ! Release memory
    deallocate( usol, grhs, bda, bdb )


contains


    pure function ue (s, t) result( return_value )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in) :: s
        real (wp), intent (in) :: t
        real (wp)              :: return_value
        !--------------------------------------------------------------------------------

        return_value = (s * t)**3 + 1.0_wp

    end function ue


    pure function uxe (s, t) result( return_value )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in) :: s
        real (wp), intent (in) :: t
        real (wp)              :: return_value
        !--------------------------------------------------------------------------------

        return_value = 3.0_wp * s**2*t**3

    end function uxe


    pure function uxxe(s, t) result( return_value )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in) :: s
        real (wp), intent (in) :: t
        real (wp)              :: return_value
        !--------------------------------------------------------------------------------

        return_value = 6.0_wp * s * (t**3)

    end function uxxe


    pure function uye(s, t) result( return_value )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in) :: s
        real (wp), intent (in) :: t
        real (wp)              :: return_value

        return_value = 3.0_wp * (s**3) * (t**2)

    end function uye


    pure function uyye (s, t) result( return_value )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in) :: s
        real (wp), intent (in) :: t
        real (wp)              :: return_value

        return_value = 6.0_wp * (s**3) * t

    end function uyye


    pure subroutine get_coefficients_in_x_direction(x, af, bf, cf)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in)  :: x
        real (wp), intent (out) :: af
        real (wp), intent (out) :: bf
        real (wp), intent (out) :: cf
        !-----------------------------------------------
        !
        !     set coefficients in the x-direction.
        !
        af = (x + 1.0_wp)**2
        bf = 2.0_wp * (x + 1.0_wp)
        cf = -x

    end subroutine get_coefficients_in_x_direction


    pure subroutine get_coefficients_in_y_direction(y, df, ef, ff)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in)  :: y
        real (wp), intent (out) :: df
        real (wp), intent (out) :: ef
        real (wp), intent (out) :: ff
        !-----------------------------------------------
        !
        !     set coefficients in y direction
        !
        df = exp(y)
        ef = 0.0
        ff = -y

    end subroutine get_coefficients_in_y_direction


end program tsepeli
