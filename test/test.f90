program test

    use modern_fishpack_library, only: &
        test_blktri, &
        test_cblktri, &
        test_cmgnbn, &
        test_genbun, &
        test_hstcrt, &
        test_hstcsp, &
        test_hstcyl, &
        test_hstplr, &
        test_hstssp, &
        test_hw3crt, &
        test_hwscrt, &
        test_hwscsp, &
        test_hwscyl, &
        test_hwsplr, &
        test_hwsssp, &
        test_pois3d, &
        test_poistg, &
        test_sepeli, &
        test_sepx4, &
        HelmholtzSolver, &
        TridiagonalSolver

    ! Explicit typing only
    implicit none

    call test_procedural_library()
    !call test_object_oriented_library()


contains


    subroutine test_object_oriented_library()
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        type (HelmholtzSolver)   :: helmholtz_solver
        type (TridiagonalSolver) :: tridiagonal_solver
        !--------------------------------------------------------------------------------

        call helmholtz_solver%unit_test()
        call tridiagonal_solver%unit_test()

    end subroutine test_object_oriented_library


    subroutine test_procedural_library()

        call test_blktri()
        call test_cblktri()
        call test_cmgnbn()
        call test_genbun()
        call test_hstcrt()
        call test_hstcsp()
        call test_hstcyl()
        call test_hstplr()
        call test_hstssp()
        call test_hw3crt()
        call test_hwscrt()
        call test_hwscsp()
        call test_hwscyl()
        call test_hwsplr()
        call test_hwsssp()
        call test_pois3d()
        call test_poistg()
        call test_sepeli()
        call test_sepx4()

    end subroutine test_procedural_library

end program test
