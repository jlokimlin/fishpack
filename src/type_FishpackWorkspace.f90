!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                         Fishpack                              *
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
! to allocate real and complex workspace arrays
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

    ! Parameters confined to the module
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: LOG_TWO = log(TWO)

    type, public :: FishpackWorkspace
        ! Type components
        real(wp),    allocatable :: real_workspace(:)
        complex(wp), allocatable :: complex_workspace(:)
        integer(ip), allocatable :: workspace_indices(:)
    contains
        ! Type-bound procedures
        procedure, public :: create => create_fishpack_workspace
        procedure, public :: destroy => destroy_fishpack_workspace
        procedure, public, nopass :: compute_blktri_workspace_lengths
        procedure, public, nopass :: compute_genbun_workspace_lengths
        procedure, public :: initialize_staggered_workspace
        procedure, public :: initialize_centered_workspace
    end type FishpackWorkspace

    ! Declare user-defined constructor
    interface FishpackWorkspace
        module procedure fishpack_workspace_constructor
    end interface

contains

    function fishpack_workspace_constructor(irwk, icwk, iiwk) &
        result (return_value)

        ! Dummy arguments
        integer(ip),              intent(in)    :: irwk ! required real workspace length
        integer(ip),              intent(in)    :: icwk ! required integer workspace length
        integer(ip), optional,    intent(in)    :: iiwk
        type(FishpackWorkspace)                 :: return_value

        ! Local variables
        integer(ip) :: iiwk_op

        ! Address optional argument
        if (present(iiwk)) then
            iiwk_op = iiwk
        else
            iiwk_op = -1
        end if

        ! Initialize workspace
        call return_value%create(irwk, icwk, iiwk_op)

    end function fishpack_workspace_constructor

    subroutine create_fishpack_workspace(self, irwk, icwk, iiwk)

        ! Dummy arguments
        class(FishpackWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: irwk ! required real workspace length
        integer(ip),              intent(in)    :: icwk ! required integer workspace length
        integer(ip), optional,    intent(in)    :: iiwk

        ! Local variables
        integer(ip) :: allocation_status

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

        ! Dummy arguments
        class(FishpackWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: n
        integer(ip),              intent(in)    :: m

        ! Local variables
        integer(ip)  :: irwk, icwk

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

        ! Dummy arguments
        class(FishpackWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: n
        integer(ip),              intent(in)    :: m

        ! Local variables
        integer(ip) :: irwk, icwk

        ! Ensure that object is usable
        call self%destroy()

        ! Compute workspace lengths for hwscrt, hwsplr
        irwk = (4 * (n + 1)) + (m + 1) * (13 + int(log(real(n + 1, kind=wp))/LOG_TWO))
        icwk = 0

        ! Allocate memory
        call self%create(irwk, icwk)

    end subroutine initialize_centered_workspace

    ! Purpose:
    !
    ! Compute nearest integer greater than or equal to
    ! logarithm base 2 of (n + 1), i.e.,
    ! log2n is smallest integer such that
    ! (n + 1) < 2**log2n
    !
    pure function get_log2n(n) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: n
        integer(ip)             :: return_value

        associate( log2n => return_value )
            log2n = 1
            do while (n + 1 > 2**log2n)
                log2n = log2n + 1
            end do
        end associate

    end function

    ! Purpose:
    !
    ! This subroutine computes the real and complex workspace
    ! requirements (generous estimate) of blktri for n,m values
    !
    pure subroutine compute_blktri_workspace_lengths(n, m, irwk, icwk)

        ! Dummy arguments
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: m
        integer(ip), intent(out) :: irwk
        integer(ip), intent(out) :: icwk

        ! Local variables
        integer(ip) :: log2n

        log2n = get_log2n(n)

        associate( k => (log2n - 2) * (2**(log2n + 1)) + 5 )
            irwk = k + max(2 * n, 6 * m) + log2n + (2 * n)
            icwk = (k + log2n)/2 + (3 * m) + n
        end associate

    end subroutine compute_blktri_workspace_lengths

    ! Purpose:
    !
    ! This subroutine computes the real workspace
    ! requirement (generously) of genbun for the current n,m
    !
    pure subroutine compute_genbun_workspace_lengths(n, m, irwk)

        ! Dummy arguments
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: m
        integer(ip), intent(out) :: irwk

        ! Local variables
        integer(ip) :: log2n

        log2n = get_log2n(n)
        irwk = 4*n + (10 + log2n)*m

    end subroutine compute_genbun_workspace_lengths

    ! Purpose:
    !
    ! This subroutine releases dynamically allocated workspace.
    ! It should be called after a fishpack solver has finished
    !
    ! Remark:
    !
    ! When intent(out) is used with a derived type, any component
    ! not assigned in a procedure will automatically deallocate on exit.
    !
    subroutine destroy_fishpack_workspace(self)

        ! Dummy arguments
        class(FishpackWorkspace), intent(out) :: self

    end subroutine destroy_fishpack_workspace

end module type_FishpackWorkspace
