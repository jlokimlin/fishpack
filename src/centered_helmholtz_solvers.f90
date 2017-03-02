module centered_helmholtz_solvers

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip, & ! Integer precision
        PI, &
        TWO_PI

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_CenteredCyclicReductionUtility, only: &
        CenteredCyclicReductionUtility

    use type_GeneralizedCyclicReductionUtility, only:&
        GeneralizedCyclicReductionUtility

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: hwscrt ! Centered cartesian solver
    public :: hwsplr ! Centered polar solver
    public :: hwscyl ! Centered cylindrical solver
    public :: hwsssp ! Centered spherical solver
    public :: hwscsp ! Centered axisymmetric spherical solver

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp

    ! Declare module subroutine interfaces
    interface
        module subroutine hwscrt(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
            bdd, elmbda, f, idimf, pertrb, ierror)

            ! Dummy arguments
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
            real(wp),    intent(inout)  :: f(:,:)
        end subroutine hwscrt

        module subroutine hwsplr(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
            bdd, elmbda, f, idimf, pertrb, ierror)

            ! Dummy arguments
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
        end subroutine hwsplr

        module subroutine hwscyl(a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd, bdc, &
            bdd, elmbda, f, idimf, pertrb, ierror)

            ! Dummy arguments
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
        end subroutine hwscyl

        module subroutine hwsssp(ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n, nbdcnd, &
            bdps, bdpf, elmbda, f, idimf, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)     :: m
            integer(ip), intent(in)     :: mbdcnd
            integer(ip), intent(in)     :: n
            integer(ip), intent(in)     :: nbdcnd
            integer(ip), intent(in)     :: idimf
            integer(ip), intent(out)    :: ierror
            real(wp),    intent(in)     :: ts
            real(wp),    intent(in)     :: tf
            real(wp),    intent(in)     :: ps
            real(wp),    intent(in)     :: pf
            real(wp),    intent(in)     :: elmbda
            real(wp),    intent(out)    :: pertrb
            real(wp),    intent(in)     :: bdts(:)
            real(wp),    intent(in)     :: bdtf(:)
            real(wp),    intent(in)     :: bdps(:)
            real(wp),    intent(in)     :: bdpf(:)
            real(wp),    intent(inout)  :: f(:,:)

        end subroutine hwsssp

        module subroutine hwscsp(intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n, &
            nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, ierror, workspace)

            ! Dummy arguments
            integer(ip), intent(in)     :: intl
            integer(ip), intent(in)     :: m
            integer(ip), intent(in)     :: mbdcnd
            integer(ip), intent(in)     :: n
            integer(ip), intent(in)     :: nbdcnd
            integer(ip), intent(in)     :: idimf
            integer(ip), intent(out)    :: ierror
            real(wp),    intent(in)     :: ts
            real(wp),    intent(in)     :: tf
            real(wp),    intent(in)     :: rs
            real(wp),    intent(in)     :: rf
            real(wp),    intent(in)     :: elmbda
            real(wp),    intent(out)    :: pertrb
            real(wp),    intent(in)     :: bdts(:)
            real(wp),    intent(in)     :: bdtf(:)
            real(wp),    intent(in)     :: bdrs(:)
            real(wp),    intent(in)     :: bdrf(:)
            real(wp),    intent(inout)  :: f(:,:)
            class(FishpackWorkspace), intent(inout)  :: workspace
        end subroutine hwscsp
    end interface

end module centered_helmholtz_solvers
