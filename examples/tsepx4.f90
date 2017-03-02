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
!     *                         Fishpack                              *
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
program test_sepx4

    use fishpack

    ! Explicit typing only
    implicit none

    ! Dictionary
    integer(ip), parameter      :: M = 32, N = M
    integer(ip), parameter      :: NX = M + 1, NY = N + 1
    integer(ip)                 :: i, j, mbdcnd, nbdcnd, idmn, iorder
    integer(ip)                 :: ierror, ierror2, ierror4
    real(wp), dimension(NX, NY) :: usol, grhs
    real(wp), dimension(NY)     :: bda, bdb
    real(wp)                    :: a, b, c, d, dx, dy, x, af, bf, cf, y
    real(wp)                    :: alpha, beta, dummy_variable(1)
    real(wp)                    :: pertrb, discretization_error
    real(wp)                    :: second_order_error, fourth_order_error
    real(wp), parameter         :: ZERO = 0.0_wp, ONE = 1.0_wp

    ! Set domain
    a = ZERO
    b = ONE
    c = ZERO
    d = ONE

    ! Set mesh size
    dx = (b - a)/M
    dy = (d - c)/N

    do i = 1, NX
        x = a + real(i - 1, kind=wp) * dx

        ! Set specified boundary conditions at y=c, d
        usol(i, 1) = ue(x, c)
        usol(i, NY) = ue(x, d)

        call get_coefficients_in_x_direction(x, af, bf, cf)

        do j = 1, NY
            y = c + real(j - 1, kind=wp)*dy

            ! Set right hand side
            grhs(i, j)=af*uxxe(x, y)+bf*uxe(x, y)+cf*ue(x, y)+uyye(x, y)
        end do
    end do

    ! Set mixed boundary conditions at x=a, b
    alpha = ONE
    beta = ONE
    do j = 1, NY
        y = c + real(j - 1, kind=wp) * dy
        bda(j) = uxe(a, y) + alpha * ue(a, y)
        bdb(j) = uxe(b, y) + beta * ue(b, y)
    end do

    ! Set boundary switches
    mbdcnd = 3
    nbdcnd = 1

    ! Set first dimension of usol, grhs and workspace length
    idmn = size(usol, dim=1)

    ! Obtain 2nd-order approximation
    iorder = 2

    call sepx4(iorder, a, b, M, mbdcnd, bda, alpha, bdb, beta, c, d, &
        N, nbdcnd, dummy_variable, dummy_variable, &
        get_coefficients_in_x_direction, grhs, usol, idmn, pertrb, ierror)

    ! Compute (relative) 2nd-order discretization error
    ! also reset specified boundaries and right hand side.
    discretization_error = ZERO
    do i = 1, NX
        x = a + real(i - 1, kind=wp) * dx
        usol(i, 1) = ue(x, c)
        usol(i, NY) = ue(x, d)
        call get_coefficients_in_x_direction(x, af, bf, cf)
        do j = 1, NY
            y = c + real(j - 1, kind=wp)*dy
            discretization_error = max(discretization_error, abs((usol(i, j)-ue(x, y))/ue(x, y)))

            ! Reset right hand side in grhs for 4th-order approximation call
            grhs(i, j)=af*uxxe(x, y)+bf*uxe(x, y)+cf*ue(x, y)+uyye(x, y)
        end do
    end do

    ! Set 2nd-order results
    second_order_error = discretization_error
    ierror2 = ierror

    ! Obtain 4th-order approximation
    iorder = 4
    call sepx4(iorder, a, b, M, mbdcnd, bda, alpha, bdb, beta, c, d, &
        N, nbdcnd, dummy_variable, dummy_variable, get_coefficients_in_x_direction, &
        grhs, usol, idmn, pertrb, ierror)

    ! Compute (relative) 4th-order discretization error
    discretization_error = ZERO
    do j = 1, NY
        y = c + real(j - 1, kind=wp)*dy
        do i = 1, NX
            x = a + real(i - 1, kind=wp)*dx
            discretization_error = max(discretization_error, abs((usol(i, j)-ue(x, y))/ue(x, y)))
        end do
    end do

    ! Set 4th-order results
    fourth_order_error = discretization_error
    ierror4 = ierror

    block
        real(wp), parameter :: KNOWN_2ND_ORDER_ERROR = 0.159850337215328e-3_wp
        real(wp), parameter :: KNOWN_4TH_ORDER_ERROR = 0.185748791237117e-5_wp

        call check_output('sepx4', &
            ierror2, KNOWN_2ND_ORDER_ERROR, second_order_error, &
            ierror4, KNOWN_4TH_ORDER_ERROR, fourth_order_error)

    end block
contains

    pure function ue(s, t) &
        result (return_value)

        ! Dummy arguments
        real(wp), intent(in) :: s
        real(wp), intent(in) :: t
        real(wp)             :: return_value

        return_value =(s*t)**3 + ONE

    end function ue

    pure function uxe(s, t) &
        result (return_value)

        ! Dummy arguments
        real(wp), intent(in) :: s
        real(wp), intent(in) :: t
        real(wp)          :: return_value

        ! Local variables
        real(wp), parameter :: THREE = 3.0_wp

        return_value = THREE *(s**2)*(t**3)

    end function uxe

    pure function uxxe(s, t) &
        result (return_value)

        ! Dummy arguments
        real(wp), intent(in) :: s
        real(wp), intent(in) :: t
        real(wp)             :: return_value

        ! Local variables
        real(wp), parameter :: SIX = 6.0_wp

        return_value = SIX * s *(t**3)

    end function uxxe

    pure function uyye(s, t) &
        result (return_value)

        ! Dummy arguments
        real(wp), intent(in) :: s
        real(wp), intent(in) :: t
        real(wp)             :: return_value

        ! Local variables
        real(wp), parameter :: SIX = 6.0_wp

        return_value = SIX *(s**3) * t

    end function uyye

    pure subroutine get_coefficients_in_x_direction(x, af, bf, cf)

        ! Dummy arguments
        real(wp), intent(in)  :: x
        real(wp), intent(out) :: af
        real(wp), intent(out) :: bf
        real(wp), intent(out) :: cf

        ! Local variables
        real(wp), parameter :: TWO = 2.0_wp

        !     set coefficients in the x-direction.
        af =(x + ONE)**2
        bf = TWO * (x + ONE)
        cf = -x

    end subroutine get_coefficients_in_x_direction

end program test_sepx4
