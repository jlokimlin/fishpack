module staggered_helmholtz_solvers

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        PI

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_CenteredCyclicReductionUtility, only: &
        CenteredCyclicReductionUtility

    use type_StaggeredCyclicReductionUtility, only: &
        StaggeredCyclicReductionUtility

    use type_GeneralizedCyclicReductionUtility, only:&
        GeneralizedCyclicReductionUtility

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hstcrt ! Staggered cartesian solver
    public :: hstplr ! Staggered polar solver
    public :: hstcyl ! Staggered cylindrical solver
    public :: hstssp ! Staggered spherical solver
    public :: hstcsp ! Staggered axisymmetric spherical solver

    ! Declare module subroutine interfaces
    interface
        module subroutine hstcrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
            bdd, elmbda, f, idimf, pertrb, ierror)
            !-----------------------------------------------
            ! Dummy arguments
            !-----------------------------------------------
            integer(ip), intent(in)    :: m
            integer(ip), intent(in)    :: mbdcnd
            integer(ip), intent(in)    :: n
            integer(ip), intent(in)    :: nbdcnd
            integer(ip), intent(in)    :: idimf
            integer(ip), intent(out)   :: ierror
            real(wp),    intent(in)    :: a
            real(wp),    intent(in)    :: b
            real(wp),    intent(in)    :: c
            real(wp),    intent(in)    :: d
            real(wp),    intent(in)    :: elmbda
            real(wp),    intent(out)   :: pertrb
            real(wp),    intent(in)    :: bda(:)
            real(wp),    intent(in)    :: bdb(:)
            real(wp),    intent(in)    :: bdc(:)
            real(wp),    intent(in)    :: bdd(:)
            real(wp),    intent(inout) :: f(:,:)
            !-----------------------------------------------
        end subroutine hstcrt

        module subroutine hstplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
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
        end subroutine hstplr

        module subroutine hstcyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
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
        end subroutine hstcyl

        module subroutine hstssp(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
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
        end subroutine hstssp

        module subroutine hstcsp(intl, a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, &
            bdc, bdd, elmbda, f, idimf, pertrb, ierror, workspace)
            !-----------------------------------------------
            ! Dummy arguments
            !-----------------------------------------------
            integer(ip), intent(inout) :: intl
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
            class(FishpackWorkspace), intent(inout) :: workspace
            !-----------------------------------------------
        end subroutine hstcsp
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

end module staggered_helmholtz_solvers
