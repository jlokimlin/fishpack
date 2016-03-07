!
!     file fish.f
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  version 1.1                    *
!     *                                                               *
!     *                 A Package of Fortran 77 and 90                *
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
! to allocate real and complex work space arrays for FFT
!
module type_FishpackWorkspace

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    ! Everything is private unless stated otherwise
    private
    public :: FishpackWorkspace

    ! Declare derived data type
    type, public :: FishpackWorkspace
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
        real (wp),    allocatable :: rew(:)
        complex (wp), allocatable :: cxw(:)
        !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Class methods
        !---------------------------------------------------------------------------------
        procedure,         public :: create => create_fish_workspace
        procedure,         public :: destroy => destroy_fish_workspace
        procedure, nopass, public :: get_block_tridiagonal_workpace_dimensions
        procedure, nopass, public :: get_genbun_workspace_dimensions
        final                     :: finalize_fish_workspace
        !---------------------------------------------------------------------------------
    end type FishpackWorkspace

contains
    !
    !*****************************************************************************************
    !
    subroutine create_fish_workspace( this, irwk, icwk, ierror )
        ! Remark:
        ! ierror is set to 20 if the dynamic allocation is unsuccessful
        ! (e.g., this would happen if m,n are too large for the computers memory
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (FishpackWorkspace), intent (in out)  :: this
        integer,                   intent(in)       :: irwk ! required real work space length
        integer (ip),              intent (in)      :: icwk ! required integer work space length
        integer,                   intent(in out)   :: ierror
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: allocation_status
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! allocate irwk words of real work space
        if (irwk > 0) then
            allocate(this%rew(irwk), stat = allocation_status)
        end if
        
        ! allocate icwk words of complex work space
        if (icwk > 0) then
            allocate(this%cxw(icwk), stat = allocation_status)
        end if
	
        ! Set error flag
        ierror = 0
        
        !  flag fatal error if allocation fails
        if (allocation_status /= 0 ) then
            ierror = 20
        end if

    end subroutine create_fish_workspace
    !
    !*****************************************************************************************
    !
    subroutine get_block_tridiagonal_workpace_dimensions( n, m, irwk, icwk)
        !
        ! Purpose:
        !
        ! This subroutine computes the real and complex work space
        ! requirements (generous estimate) of blktri for n,m values
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent(in)  :: n,m
        integer (ip), intent(out) :: irwk
        integer (ip), intent(out) :: icwk
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: log2n
        !--------------------------------------------------------------------------------
        !       compute nearest integer greater than or equal to
        !       log base 2 of n+1, i.e., log2n is smallest integer
        !       such that 2**log2n >= n+1
        log2n = 1
        do
            log2n = log2n+1
            if (n+1 <= 2**log2n) exit
        end do

        associate( l => 2**(log2n+1) )
            irwk = (log2n-2)*l+5+max(2*n,6*m)+log2n+2*n
            icwk = ((log2n-2)*l+5+log2n)/2+3*m+n
        end associate

    end subroutine get_block_tridiagonal_workpace_dimensions
    !
    !*****************************************************************************************
    !
    subroutine get_genbun_workspace_dimensions( n, m, irwk )
        !
        ! Purpose:
        !
        ! This subroutine computes the real work space
        ! requirement (generously) of genbun for the current n,m
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (ip), intent(in)  :: n
        integer (ip), intent (in) :: m
        integer (ip), intent(out) :: irwk
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: log2n
        !--------------------------------------------------------------------------------

        !       compute nearest integer greater than or equal to
        !       log base 2 of n+1, i.e., log2n is smallest integer
        !       such that 2**log2n >= n+1
        log2n = 1
        do
            log2n = log2n+1
            if (n+1 <= 2**log2n) exit
        end do

        irwk = 4*n + (10 + log2n)*m

    end subroutine get_genbun_workspace_dimensions
    !
    !*****************************************************************************************
    !
    subroutine destroy_fish_workspace( this )
        !
        ! Purpose:
        !
        ! This subroutine releases dynamically allocated work space.
        ! It should be called after a fishpack solver has finished
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (FishpackWorkspace), intent (in out) :: this
        !--------------------------------------------------------------------------------

        ! Free dynamically allocated real workspace array
        if ( allocated( this%rew ) ) deallocate( this%rew )

        ! Release dynamically allocated complex workspace array
        if ( allocated( this%cxw ) )  deallocate( this%cxw )

    end subroutine destroy_fish_workspace
    !
    !*****************************************************************************************
    !
    subroutine finalize_fish_workspace( this )
        !
        ! Purpose:
        !
        ! This subroutine releases dynamically allocated work space.
        ! It should be called after a fishpack solver has finished
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (FishpackWorkspace), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_fish_workspace
    !
    !*****************************************************************************************
    !
end module type_FishpackWorkspace
