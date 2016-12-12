module type_FFTpack

    use fishpack_precision, only: &
        wp ! Working precision

    use fftpack_routines, only: &
        rffti, & ! initialize  rfftf and rfftb
        rfftf, & ! forward transform of a real periodic sequence
        rfftb, & ! backward transform of a real coefficient array
        ezffti, & !initialize ezfftf and ezfftb
        ezfftf, & ! a simplified real periodic forward transform
        ezfftb, & ! a simplified real periodic backward transform
        sinti, & ! initialize sint
        sint, & ! sine transform of a real odd sequence
        costi, & ! initialize cost
        cost, & ! cosine transform of a real even sequence
        sinqi, & ! initialize sinqf and sinqb
        sinqf, & ! forward sine transform with odd wave numbers
        sinqb, & ! unnormalized inverse of sinqf
        cosqi, & ! initialize cosqf and cosqb
        cosqf, & ! forward cosine transform with odd wave numbers
        cosqb, & ! unnormalized inverse of cosqf
        cffti, & ! initialize cfftf and cfftb
        cfftf, & ! forward transform of a complex periodic sequence
        cfftb ! unnormalized inverse of cfftf

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: FFTpack

    type, public :: FFTpack
        !-------------------------------------------------------
        ! Type components
        !-------------------------------------------------------
        real(wp), allocatable :: workspace(:)
        !-------------------------------------------------------
    contains
        !-------------------------------------------------------
        ! Type-bound procedures
        !-------------------------------------------------------
        procedure, nopass, public :: rffti
        procedure, nopass, public :: rfftf
        procedure, nopass, public :: rfftb
        procedure, nopass, public :: ezffti
        procedure, nopass, public :: ezfftf
        procedure, nopass, public :: ezfftb
        procedure, nopass, public :: sinti
        procedure, nopass, public :: sint
        procedure, nopass, public :: costi
        procedure, nopass, public :: cost
        procedure, nopass, public :: sinqi
        procedure, nopass, public :: sinqf
        procedure, nopass, public :: sinqb
        procedure, nopass, public :: cosqi
        procedure, nopass, public :: cosqf
        procedure, nopass, public :: cosqb
        procedure, nopass, public :: cffti
        procedure, nopass, public :: cfftf
        procedure, nopass, public :: cfftb
        !-------------------------------------------------------
    end type FFTpack

end module type_FFTpack
