module fishpack

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT, &
        stderr => ERROR_UNIT

    use fishpack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, HALF_PI, TWO_PI

    use real_block_tridiagonal_linear_systems_solver, &
        only: blktri

    use complex_block_tridiagonal_linear_systems_solver, only: &
        cblktri

    use complex_linear_systems_solver, only: &
        cmgnbn

    use centered_real_linear_systems_solver, only: &
        genbun

    use three_dimensional_solvers, only: &
        pois3d, & ! general_linear_systems_solver_3d
        hw3crt ! centered_cartesian_helmholtz_solver_3d

    use staggered_real_linear_systems_solver, only: &
        poistg

    use module_sepeli, only: &
        sepeli

    use module_sepx4, only: &
        sepx4

    use staggered_helmholtz_solvers, only: &
        hstcrt, & ! Staggered cartesian solver
        hstplr, & ! Staggered polar solver
        hstcyl, & ! Staggered cylindrical solver
        hstssp, & ! Staggered spherical solver
        hstcsp ! Staggered axisymmetric spherical solver

    use centered_helmholtz_solvers, only: &
        hwscrt, & ! Centered cartesian solver
        hwsplr, & ! Centered polar solver
        hwscyl, & ! Centered cylindrical solver
        hwsssp, & ! Centered spherical solver
        hwscsp ! Centered axisymmetric spherical solver

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_PeriodicFastFourierTransform, only: &
        PeriodicFastFourierTransform

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: wp, ip
    public :: PI, HALF_PI, TWO_PI
    public :: FishpackWorkspace
    public :: PeriodicFastFourierTransform
    public :: blktri
    public :: cblktri
    public :: cmgnbn
    public :: genbun
    public :: hstcrt
    public :: hstcsp
    public :: hstcyl
    public :: hstplr
    public :: hstssp
    public :: hw3crt
    public :: hwscrt
    public :: hwscsp
    public :: hwscyl
    public :: hwsplr
    public :: hwsssp
    public :: pois3d
    public :: poistg
    public :: sepeli
    public :: sepx4
    public :: check_output, assert_equal

    interface check_output
        module procedure common_check_output
        module procedure check_output_with_pertrb
        module procedure check_output_with_2nd_and_4th_order_errors
    end interface

contains

    subroutine common_check_output(routine, ierror, input, output)

        ! Dummy arguments
        character(len=*), intent(in) :: routine
        integer(ip),      intent(in) :: ierror
        real(wp),         intent(in) :: input
        real(wp),         intent(in) :: output

        ! Print earlier output from platforms with 64-bit floating point
        ! arithmetic followed by the output from this computer
        if (assert_equal(input, output)) then
            write(stdout, '(a)') 'PASS'
        else
            write(stderr, '(a)')                '======================================================================'
            write(stderr, '(/a/)')              '     ***FAIL: '//routine//' test program***'
            write(stderr, '(a)')                '     Previous 64-bit floating point arithmetic result '
            write(stderr, '(a,e23.15e3)')       '     ierror = 0,  discretization error = ', input
            write(stderr, '(/a)')               '     The output from your computer is: '
            write(stderr, '(a,i3,a,e23.15e3/)') '     ierror = ', ierror, ' discretization error = ', output
            write(stderr, '(a)')                '======================================================================'
        end if

    end subroutine common_check_output

    subroutine check_output_with_pertrb(routine, ierror, &
        input, output, pertrb_in, pertrb_out)

        ! Dummy arguments
        character(len=*), intent(in) :: routine
        integer(ip),      intent(in) :: ierror
        real(wp),         intent(in) :: input, pertrb_in
        real(wp),         intent(in) :: output, pertrb_out

        ! Print earlier output from platforms with 64-bit floating point
        ! arithmetic followed by the output from this computer
        if (assert_equal(input, output)) then
            write(stdout, '(a)') 'PASS'
        else
            write(stderr, '(a)')               '======================================================================'
            write(stderr, '(/a/)')             '   ***FAIL: '//routine//' test program***'
            write(stderr, '(a)')               '   Previous 64-bit floating point arithmetic result '
            write(stderr, '(a,e23.15e3)')      '   ierror = 0, pertrb = ', pertrb_in
            write(stderr, '(a,e23.15e3)')      '   discretization error = ', input
            write(stdout, '(/a)')              '   The output from your computer is: '
            write(stdout, '(a,i3,a,e23.15e3)') '   ierror = ', ierror, ', pertrb = ', pertrb_out
            write(stdout, '(a,e23.15e3/)')     '   discretization error = ', output
            write(stderr, '(a)')               '======================================================================'
        end if

    end subroutine check_output_with_pertrb

    subroutine check_output_with_2nd_and_4th_order_errors(routine, &
        ierror2, input2, output2, ierror4, input4, output4)

        ! Dummy arguments
        character(len=*), intent(in) :: routine
        integer(ip),      intent(in) :: ierror2, ierror4
        real(wp),         intent(in) :: input2, output2, input4, output4

        ! Print earlier output from platforms with 64-bit floating point
        ! arithmetic followed by the output from this computer
        if (assert_equal(input2, output2) .and. assert_equal(input4, output4)) then
            write(stdout, '(a)') 'PASS'
        else
            write(stderr, '(a)')            '======================================================================'
            write(stderr, '(/a/)')          '   ***FAIL: '//routine//' test program***'
            write(stderr, '(a)')            '   Previous 64-bit floating point arithmetic result '
            write(stdout, '(a)')            '   2nd-order ierror = 0'
            write(stdout, '(a,e23.15e3)')   '   2nd-order discretization error =', input2
            write(stdout, '(a)')            '   4th-order ierror = 0'
            write(stdout, '(a,e23.15e3/)')  '   4th-order discretization error =', input4
            write(stdout, '(/a)')           '   The output from your computer is: '
            write(stdout, '(a,i3)')         '   2nd-order ierror =', ierror2
            write(stdout, '(a,e23.15e3)')   '   2nd-order discretization error =', output2
            write(stdout, '(a,i3)')         '   4th-order ierror =', ierror4
            write(stdout, '(a,e23.15e3/)')  '   4th-order discretization error =', output4
            write(stderr, '(a)')            '======================================================================'
        end if

    end subroutine check_output_with_2nd_and_4th_order_errors

    pure function assert_equal(first, second) &
        result(return_value)

        ! Dummy arguments
        real(wp), intent(in) :: first
        real(wp), intent(in) :: second
        logical              :: return_value

        ! Local variables
        real(wp), parameter :: ZERO = 0.0_wp

        return_value = (abs(first - second) > ZERO)

    end function assert_equal

end module fishpack
