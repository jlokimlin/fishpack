module modern_fishpack_library

    use module_blktri, only: &
        blktri, &
        test_blktri

    use module_cblktri, only: &
        cblktri, &
        test_cblktri

    use module_cmgnbn, only: &
        cmgnbn, &
        test_cmgnbn

    use module_genbun, only: &
        genbun, &
        test_genbun

    use module_hstcrt, only: &
        hstcrt, &
        test_hstcrt

    use module_hstcsp, only: &
        hstcsp, &
        test_hstcsp

    use module_hstcyl, only: &
        hstcyl, &
        test_hstcyl

    use module_hstplr, only: &
        hstplr, &
        test_hstplr

    use module_hstssp, only: &
        hstssp, &
        test_hstssp

    use module_hw3crt, only: &
        hw3crt, &
        test_hw3crt

    use module_hwscrt, only: &
        hwscrt, &
        test_hwscrt

    use module_hwscsp, only: &
        hwscsp, &
        test_hwscsp

    use module_hwscyl, only: &
        hwscyl, &
        test_hwscyl

    use module_hwsplr, only: &
        hwsplr, &
        test_hwsplr

    use module_hwsssp, only: &
        hwsssp, &
        test_hwsssp

    use module_pois3d, only: &
        pois3d, &
        test_pois3d

    use module_poistg, only: &
        poistg, &
        test_poistg

    use module_sepeli, only: &
        sepeli, &
        test_sepeli

    use module_sepx4, only: &
        sepx4, &
        test_sepx4

    use type_HelmholtzData, only: &
        HelmholtzData

    use type_HelmholtzSolver, only: &
        HelmholtzSolver

    use type_TridiagonalData, only: &
        TridiagonalData

    use type_TridiagonalSolver, only: &
        TridiagonalSolver

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
    public :: test_blktri
    public :: cblktri
    public :: test_cblktri
    public :: cmgnbn
    public :: test_cmgnbn
    public :: genbun
    public :: test_genbun
    public :: hstcrt
    public :: test_hstcrt
    public :: hstcsp
    public :: test_hstcsp
    public :: hstcyl
    public :: test_hstcyl
    public :: hstplr
    public :: test_hstplr
    public :: hstssp
    public :: test_hstssp
    public :: hw3crt
    public :: test_hw3crt
    public :: hwscrt
    public :: test_hwscrt
    public :: hwscsp
    public :: test_hwscsp
    public :: hwscyl
    public :: test_hwscyl
    public :: hwsplr
    public :: test_hwsplr
    public :: hwsssp
    public :: test_hwsssp
    public :: pois3d
    public :: test_pois3d
    public :: poistg
    public :: test_poistg
    public :: sepeli
    public :: test_sepeli
    public :: sepx4
    public :: test_sepx4
    public :: HelmholtzData
    public :: HelmholtzSolver
    public :: TridiagonalData
    public :: TridiagonalSolver
    public :: CenteredGrid
    public :: StaggeredGrid
    public :: MacGrid

end module modern_fishpack_library
