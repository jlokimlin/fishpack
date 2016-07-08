module fishpack_library

    use fishpack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, & ! machine precision pi
        TWO_PI

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
    public :: PI, TWO_PI
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

end module fishpack_library
