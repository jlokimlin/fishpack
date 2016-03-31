module modern_fishpack_library

    use module_blktri, only: &
        blktri

    use module_cblktri, only: &
        cblktri

    use module_cmgnbn, only: &
        cmgnbn

    use module_genbun, only: &
        genbun

    use module_hstcrt, only: &
        hstcrt

    use module_hstcsp, only: &
        hstcsp

    use module_hstcyl, only: &
        hstcyl

    use module_hstplr, only: &
        hstplr

    use module_hstssp, only: &
        hstssp

    use module_hw3crt, only: &
        hw3crt

    use module_hwscrt, only: &
        hwscrt

    use module_hwscsp, only: &
        hwscsp

    use module_hwscyl, only: &
        hwscyl

    use module_hwsplr, only: &
        hwsplr

    use module_hwsssp, only: &
        hwsssp

    use module_pois3d, only: &
        pois3d

    use module_poistg, only: &
        poistg

    use module_sepeli, only: &
        sepeli

    use module_sepx4, only: &
        sepx4

    use type_FishpackWorkspace, only: &
    FishpackWorkspace

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
    public :: FishpackWorkspace
    public :: HelmholtzData
    public :: HelmholtzSolver
    public :: PoissonSolver
    public :: TridiagonalData
    public :: TridiagonalSolver
    public :: Grid
    public :: CenteredGrid
    public :: StaggeredGrid
    public :: MacGrid

end module modern_fishpack_library
