!
!     file type_FishpackWorkspace.f90
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  Version 1.1                    *
!     *                                                               *
!     *                      A Package of Fortran                     *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *               for Modeling Geophysical Processes              *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        John Adams, Paul Swarztrauber and Roland Sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! Purpose:
!
! This module is used by all fishpack solvers
! to allocate real and complex work space arrays for FFTpack
!
module type_FishpackWorkspace

    use fishpack_precision, only: &
        wp, & ! Working precision
        ip ! Integer precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: FishpackWorkspace


    ! Declare derived data type
    type, public :: FishpackWorkspace
        !---------------------------------------------------------------
        ! Type components
        !---------------------------------------------------------------
        real (wp),    allocatable :: real_workspace(:)
        complex (wp), allocatable :: complex_workspace(:)
        integer (ip), allocatable :: workspace_indices(:)
        !---------------------------------------------------------------
    contains
        !---------------------------------------------------------------
        ! Type-bound procedures
        !---------------------------------------------------------------
        procedure,         public :: create => create_fishpack_workspace
        procedure,         public :: destroy => destroy_fishpack_workspace
        procedure, nopass, public :: get_block_tridiagonal_workpace_dimensions
        procedure, nopass, public :: get_genbun_workspace_dimensions
        !final                     :: finalize_fishpack_workspace
        !---------------------------------------------------------------
    end type FishpackWorkspace



contains



    subroutine create_fishpack_workspace(this, irwk, icwk, ierror)
        !
        ! Remark:
        !
        ! ierror is set to 20 if the dynamic allocation is unsuccessful
        ! (e.g., this would happen if m,n are too large for the computers memory
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class (FishpackWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: irwk ! required real work space length
        integer (ip),              intent (in)     :: icwk ! required integer work space length
        integer (ip), optional,    intent (in out) :: ierror
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer (ip) :: allocation_status
        !--------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        !
        !==> allocate irwk words of real workspace
        !
        if (irwk > 0) then

            ! Allocate memory
            allocate(this%real_workspace(irwk), stat=allocation_status)
            !
            ! Check if allocation was successful
            if (allocation_status /= 0 ) then
                error stop 'Object of class (FishpackWorkspace): '&
                    //'failed to allocate real_workspace array '&
                    //'in create_fishpack_workspace'
            end if
        end if

        !
        !==> allocate icwk words of complex workspace
        !
        if (icwk > 0) then

            ! Allocate memory
            allocate(this%complex_workspace(icwk), stat=allocation_status)

            ! Check if allocation was successful
            if (allocation_status /= 0 ) then
                error stop 'Object of class (FishpackWorkspace): '&
                    //'failed to allocate complex_workspace array '&
                    //'in create_fishpack_workspace'
            end if
        end if
        
        ! Address error flag
        if (present(ierror)) then
            if (allocation_status /= 0 ) then
                ierror = 20
            else
                ierror = 0
            end if
        end if

    end subroutine create_fishpack_workspace



    subroutine get_block_tridiagonal_workpace_dimensions(n, m, irwk, icwk)
        !
        ! Purpose:
        !
        ! This subroutine computes the real and complex work space
        ! requirements (generous estimate) of blktri for n,m values
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: n
        integer (ip), intent (in)  :: m
        integer (ip), intent (out) :: irwk
        integer (ip), intent (out) :: icwk
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer (ip) :: log2n
        !--------------------------------------------------------------

        !
        !==> Compute nearest integer greater than or equal to
        !    log base 2 of n+1, i.e., log2n is smallest integer
        !    such that 2**log2n >= n+1
        !
        log2n = 1

        do while (n+1 > 2**log2n)
            log2n = log2n+1
        end do

        associate( l => 2**(log2n+1) )

            irwk = (log2n-2)*l+5+max(2*n,6*m)+log2n+2*n
            icwk = ((log2n-2)*l+5+log2n)/2+3*m+n

        end associate

    end subroutine get_block_tridiagonal_workpace_dimensions



    subroutine get_genbun_workspace_dimensions(n, m, irwk)
        !
        ! Purpose:
        !
        ! This subroutine computes the real work space
        ! requirement (generously) of genbun for the current n,m
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: n
        integer (ip), intent (in)  :: m
        integer (ip), intent (out) :: irwk
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer (ip) :: log2n
        !--------------------------------------------------------------

        !
        !==> Compute nearest integer greater than or equal to
        !    log base 2 of n+1, i.e., log2n is smallest integer
        !    such that 2**log2n >= n+1
        !
        log2n = 1

        do while (n+1 > 2**log2n)
            log2n = log2n+1
        end do

        irwk = 4*n + (10 + log2n)*m

    end subroutine get_genbun_workspace_dimensions



    subroutine destroy_fishpack_workspace(this)
        !
        ! Purpose:
        !
        ! This subroutine releases dynamically allocated work space.
        ! It should be called after a fishpack solver has finished
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class (FishpackWorkspace), intent (in out) :: this
        !--------------------------------------------------------------

        !
        !==> Release memory
        !
        if (allocated(this%real_workspace)) then
            deallocate( this%real_workspace )
        end if

        if (allocated(this%complex_workspace)) then
            deallocate( this%complex_workspace )
        end if

        if (allocated(this%workspace_indices)) then
            deallocate( this%workspace_indices )
        end if

    end subroutine destroy_fishpack_workspace



    subroutine finalize_fishpack_workspace(this)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        type (FishpackWorkspace), intent (in out) :: this
        !--------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_fishpack_workspace



end module type_FishpackWorkspace
