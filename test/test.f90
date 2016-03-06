program test

    use module_blktri, only: &
        blktri_unit_test

!    use module_cblktri, only: &
!        cblktri_unit_test

    use module_genbun, only: &
        genbun_unit_test

    use module_hstcrt, only: &
        hstcrt_unit_test

    use module_hstcsp, only: &
        hstcsp_unit_test

    use module_hstcyl, only: &
        hstcyl_unit_test

    use module_hstplr, only: &
        hstplr_unit_test

    use module_hstssp, only: &
        hstssp_unit_test

    use module_hwscrt, only: &
        hwscrt_unit_test

    use module_poistg, only: &
        poistg_unit_test

    call blktri_unit_test()
    !call cblktri_unit_test()
    call genbun_unit_test()
    call hstcrt_unit_test()
    call hstcsp_unit_test()
    call hstcyl_unit_test()
    call hstplr_unit_test()
    call hstssp_unit_test()
    call hwscrt_unit_test()
    call poistg_unit_test()

end program test
