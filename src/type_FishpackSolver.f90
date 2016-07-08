!
! Purpose:
!
! This derived data type encloses all of fishpack's procedural solvers inside a single class
!
module type_FishpackSolver

    use module_blktri, only: blktri

    use module_cblktri, only: cblktri

    use module_cmgnbn, only: cmgnbn

    use module_genbun, only: genbun

    use module_hstcrt, only: hstcrt

    use module_hstcsp, only: hstcsp

    use module_hstcyl, only: hstcyl

    use module_hstplr, only: hstplr

    use module_hstssp, only: hstssp

    use module_hw3crt, only: hw3crt

    use module_hwscrt, only: hwscrt

    use module_hwscsp, only: hwscsp

    use module_hwscyl, only: hwscyl

    use module_hwsplr, only: hwsplr

    use module_hwsssp, only: hwsssp

    use module_pois3d, only: pois3d

    use module_poistg, only: poistg

    use module_sepeli, only: sepeli

    use module_sepx4, only: sepx4

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: FishpackSolver

    ! Declare derived data type
    type, public :: FishpackSolver
    contains
        !---------------------------------------------------------------
        ! Type-bound procedures
        !---------------------------------------------------------------
        procedure, nopass, public :: blktri
        procedure, nopass, public :: cblktri
        procedure, nopass, public :: cmgnbn
        procedure, nopass, public :: genbun
        procedure, nopass, public :: hstcrt
        procedure, nopass, public :: hstcsp
        procedure, nopass, public :: hstcyl
        procedure, nopass, public :: hstplr
        procedure, nopass, public :: hstssp
        procedure, nopass, public :: hw3crt
        procedure, nopass, public :: hwscrt
        procedure, nopass, public :: hwscsp
        procedure, nopass, public :: hwscyl
        procedure, nopass, public :: hwsplr
        procedure, nopass, public :: hwsssp
        procedure, nopass, public :: pois3d
        procedure, nopass, public :: poistg
        procedure, nopass, public :: sepeli
        procedure, nopass, public :: sepx4
        !---------------------------------------------------------------
    end type FishpackSolver


end module type_FishpackSolver



