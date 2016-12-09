module helmholtz_solvers

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    use type_FishpackWorkspace, only: &
        Fish => FishpackWorkspace

    use module_genbun, only: &
        genbunn

    use module_poistg, only: &
        poistgg

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hwsplr ! Centered polar solver
    !public :: hstcrt ! Staggered cartesian solver
    !public :: hwscrt ! Centered cartesian solver

    ! Declare module subroutine interfaces
    interface
        module subroutine hwsplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
            bdd, elmbda, f, idimf, pertrb, ierror)
            !-----------------------------------------------
            ! Dummy arguments
            !-----------------------------------------------
            integer(ip), intent(in)     :: m
            integer(ip), intent(in)     :: mbdcnd
            integer(ip), intent(in)     :: n
            integer(ip), intent(in)     :: nbdcnd
            integer(ip), intent(in)     :: idimf
            integer(ip), intent(out)    :: ierror
            real(wp),    intent(in)     :: a
            real(wp),    intent(in)     :: b
            real(wp),    intent(in)     :: c
            real(wp),    intent(in)     :: d
            real(wp),    intent(in)     :: elmbda
            real(wp),    intent(out)    :: pertrb
            real(wp),    intent(in)     :: bda(:)
            real(wp),    intent(in)     :: bdb(:)
            real(wp),    intent(in)     :: bdc(:)
            real(wp),    intent(in)     :: bdd(:)
            real(wp),    intent(inout) :: f(:,:)
            !-----------------------------------------------
        end subroutine hwsplr
    end interface

    !---------------------------------------------------------------
    ! Parameters confined to the module
    !---------------------------------------------------------------
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    !---------------------------------------------------------------

contains

end module helmholtz_solvers
