module fishpack_library

    use fishpack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, & ! machine precision pi
        HALF_PI, &
        TWO_PI

    use real_block_tridiagonal_linear_systems_solver, &
        only: blktri

    use module_cblktri, only: cblktri

    use module_cmgnbn, only: cmgnbn

    use module_genbun, only: genbun

    use module_hw3crt, only: hw3crt

    use module_pois3d, only: pois3d

    use module_poistg, only: poistg

    use module_sepeli, only: sepeli

    use module_sepx4, only: sepx4

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

    use type_FishpackGrid, only: &
        FishpackGrid

    use type_FishpackWorkspace, only: &
        FishpackWorkspace

    use type_FishpackSolver, only: &
        FishpackSolver

    use type_HelmholtzData, only: &
        HelmholtzData

    use type_HelmholtzSolver, only: &
        HelmholtzSolver

    use type_PoissonSolver, only: &
        PoissonSolver

    use type_TridiagonalData, only: &
        TridiagonalData

    use type_TridiagonalSolver, only: &
        TridiagonalSolver

    use type_Grid, only: &
        Grid

    use type_CenteredGrid, only: &
        CenteredGrid

    use type_StaggeredGrid, only: &
        StaggeredGrid

    use type_MacGrid, only: &
        MacGrid

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: wp, ip
    public :: PI, HALF_PI, TWO_PI
    public :: FishpackGrid
    public :: FishpackWorkspace
    public :: FishpackSolver
    public :: HelmholtzData
    public :: HelmholtzSolver
    public :: PoissonSolver
    public :: TridiagonalData
    public :: TridiagonalSolver
    public :: Grid
    public :: CenteredGrid
    public :: StaggeredGrid
    public :: MacGrid
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

end module fishpack_library
