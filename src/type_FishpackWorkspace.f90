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

    type, public :: FishpackWorkspace
        !---------------------------------------------------------------
        ! Type components
        !---------------------------------------------------------------
        real(wp),    allocatable :: real_workspace(:)
        complex(wp), allocatable :: complex_workspace(:)
        integer(ip), allocatable :: workspace_indices(:)
        !---------------------------------------------------------------
    contains
        !---------------------------------------------------------------
        ! Type-bound procedures
        !---------------------------------------------------------------
        procedure, public :: create => create_fishpack_workspace
        procedure, public :: destroy => destroy_fishpack_workspace
        procedure, public :: compute_blktri_workspace_lengths
        procedure, public :: compute_genbun_workspace_lengths
        procedure, public :: initialize_staggered_workspace
        procedure, public :: initialize_centered_workspace
        procedure, public :: initialize_blktri_workspace
        !final             :: finalize_fishpack_workspace
        !---------------------------------------------------------------
    end type FishpackWorkspace

contains

    subroutine create_fishpack_workspace(self, irwk, icwk, iiwk)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(FishpackWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: irwk ! required real work space length
        integer(ip),              intent(in)    :: icwk ! required integer work space length
        integer(ip), optional,    intent(in)    :: iiwk
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: allocation_status
        !--------------------------------------------------------------

        ! Ensure that object is usable
        call self%destroy()

        if (irwk > 0) then
            ! Allocate irwk words of real workspace
            allocate( self%real_workspace(irwk), stat=allocation_status )

            ! Check if real allocation was successful
            if (allocation_status /= 0 ) then
                error stop 'Object of class(FishpackWorkspace): '&
                    //'failed to allocate real_workspace array '&
                    //'in create_fishpack_workspace'
            end if
        end if

        if (icwk > 0) then
            ! Allocate icwk words of complex workspace
            allocate( self%complex_workspace(icwk), stat=allocation_status )

            ! Check if complex allocation was successful
            if (allocation_status /= 0 ) then
                error stop 'Object of class(FishpackWorkspace): '&
                    //'failed to allocate complex_workspace array '&
                    //'in create_fishpack_workspace'
            end if
        end if
        
        ! Address optional argument
        if (present(iiwk)) then
            if (iiwk > 0) then
                ! Allocate iwwk words of integer workspace
                allocate( self%workspace_indices(iiwk), stat=allocation_status )

                ! Check if integer allocation was successful
                if (allocation_status /= 0 ) then
                    error stop 'Object of class(FishpackWorkspace): '&
                        //'failed to allocate workspace_indices array '&
                        //'in create_fishpack_workspace'
                end if
            end if
        end if

    end subroutine create_fishpack_workspace

    subroutine initialize_staggered_workspace(self, n, m)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        class(FishpackWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: n
        integer(ip),              intent(in)    :: m
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip)  :: irwk, icwk
        !-----------------------------------------------

        ! Ensure that object is usable
        call self%destroy()

        ! Get workspace dimensions for genbun
        call self%compute_genbun_workspace_lengths(n, m, irwk)

        ! Adjust workspace for hstcyl, hstplr
        irwk = irwk + 3 * m
        icwk = 0

        ! Allocate memory
        call self%create(irwk, icwk)

    end subroutine initialize_staggered_workspace

    subroutine initialize_centered_workspace(self, n, m)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        class(FishpackWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: n
        integer(ip),              intent(in)    :: m
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip)         :: irwk, icwk
        real(wp), parameter :: TWO = 2.0_wp
        !-----------------------------------------------

        ! Ensure that object is usable
        call self%destroy()

        ! Compute workspace lengths for hwscrt, hwsplr
        irwk = 4*(n+1)+(m+1)*(13+int(log(real(n+1, kind=wp))/log(TWO)))
        icwk = 0

        ! Allocate memory
        call self%create(irwk, icwk)

    end subroutine initialize_centered_workspace

    subroutine initialize_blktri_workspace(self, n, m)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        class(FishpackWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: n
        integer(ip),              intent(in)    :: m
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer(ip) :: irwk, icwk
        !-----------------------------------------------

        ! Ensure that object is usable
        call self%destroy()

        ! compute and allocate real and complex required work space
        call self%compute_blktri_workspace_lengths(n, m, irwk, icwk)

        ! Allocate memory
        call self%create(irwk, icwk)

    end subroutine initialize_blktri_workspace

    pure subroutine compute_blktri_workspace_lengths(self, n, m, irwk, icwk)
        !
        ! Purpose:
        !
        ! This subroutine computes the real and complex work space
        ! requirements (generous estimate) of blktri for n,m values
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(FishpackWorkspace), intent(inout) :: self
        integer(ip),              intent(in)     :: n
        integer(ip),              intent(in)     :: m
        integer(ip),              intent(out)    :: irwk
        integer(ip),              intent(out)    :: icwk
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: log2n
        !--------------------------------------------------------------

        !
        ! Compute nearest integer greater than or equal to
        !    log base 2 of n+1, i.e., log2n is smallest integer
        !    such that 2**log2n >= n+1
        !
        log2n = 1

        do
            if (n+1 <= 2**log2n) exit
            log2n = log2n+1
        end do

        associate( l => 2**(log2n+1) )
            irwk = (log2n-2)*l+5+max(2*n,6*m)+log2n+2*n
            icwk = ((log2n-2)*l+5+log2n)/2+3*m+n
        end associate

    end subroutine compute_blktri_workspace_lengths

    pure subroutine compute_genbun_workspace_lengths(self, n, m, irwk)
        !
        ! Purpose:
        !
        ! This subroutine computes the real work space
        ! requirement (generously) of genbun for the current n,m
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(FishpackWorkspace), intent(inout) :: self
        integer(ip),              intent(in)     :: n
        integer(ip),              intent(in)     :: m
        integer(ip),              intent(out)    :: irwk
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip) :: log2n
        !--------------------------------------------------------------

        !
        ! Compute nearest integer greater than or equal to
        !    log base 2 of n+1, i.e., log2n is smallest integer
        !    such that 2**log2n >= n+1
        !
        log2n = 1

        do
            if (n+1 <= 2**log2n) exit
            log2n = log2n+1
        end do

        irwk = 4*n + (10 + log2n)*m

    end subroutine compute_genbun_workspace_lengths

    subroutine destroy_fishpack_workspace(self)
        !
        ! Purpose:
        !
        ! This subroutine releases dynamically allocated work space.
        ! It should be called after a fishpack solver has finished
        !
        ! Remark:
        !
        ! When intent(out) is used with a derived type, any component
        ! not assigned in a procedure will automatically deallocate on exit.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        class(FishpackWorkspace), intent(out) :: self
        !--------------------------------------------------------------

    end subroutine destroy_fishpack_workspace

    subroutine finalize_fishpack_workspace(self)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        type(FishpackWorkspace), intent(inout) :: self
        !--------------------------------------------------------------

        call self%destroy()

    end subroutine finalize_fishpack_workspace

end module type_FishpackWorkspace
