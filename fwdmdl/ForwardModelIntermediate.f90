! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelIntermediate
!=============================================================================

  ! This module simply defines a data type passed around forward model routines.
  ! It contains hydrostatic information and so on, and is remembered by the
  ! calling code.

  use MLSCommon, only: R8

  implicit none
  private

  type, public :: ForwardModelIntermediate_T

    ! These ones are for the scan model
    real (r8), dimension(:,:),             pointer :: BasisGph=>NULL()
    real (r8), dimension(:),               pointer :: R=>NULL()
    real (r8), dimension(:,:),             pointer :: RT=>NULL()
    integer :: BelowRef                 ! T. basis at or below refGPH

  end type ForwardModelIntermediate_T

  type, public :: ForwardModelStatus_T
    integer :: Flags                    ! Bits indicate problems.  See B_...
    integer :: Maf                      ! The MAF to process
    logical :: NewScanHydros            ! Scan model needs to recompute hydrostatic
    logical, dimension(:), pointer :: Rows=>NULL() ! Flag to indicate this row has non zeros
  end type ForwardModelStatus_T

  ! Bits set in Flags:
  integer, public, parameter :: B_Metrics = 1 ! Trouble with metrics
  integer, public, parameter :: B_Ptg_Angles = 2*b_metrics ! Patched some
    ! out-of-order pointing angles
  integer, public, parameter :: B_Refraction = 2*b_ptg_angles ! Computation of
    ! refraction coefficient failed or was not as accurate as desired.

  public :: DestroyForwardModelIntermediate

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here
  !---------------------------------------------------------------------------

contains

  ! -----------------------------  DestroyForwardModelIntemediate  -----
  subroutine DestroyForwardModelIntermediate ( ifm )

    use Allocate_Deallocate, only: DEALLOCATE_TEST

    type (ForwardModelIntermediate_T), intent(inout) :: ifm

    ! Exectuable code

    call Deallocate_test ( ifm%basisGPH, 'basisGPH', ModuleName )
    call Deallocate_test ( ifm%RT, 'RT', ModuleName )
    call Deallocate_test ( ifm%R, 'R', ModuleName )

  end subroutine DestroyForwardModelIntermediate

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ForwardModelIntermediate

! $Log$
! Revision 2.4  2003/07/09 00:53:10  vsnyder
! Remove unreferenced variables
!
! Revision 2.3  2002/10/08 17:08:03  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2002/08/22 00:05:42  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.1  2001/10/02 16:51:41  livesey
! Removed fmStat%finished and reordered loops in forward models
!
! Revision 2.0  2001/09/17 20:26:25  livesey
! New forward model
!
! Revision 1.10.2.3  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.10.2.2  2001/09/07 22:50:09  livesey
! Very intermediate
!
! Revision 1.10.2.1  2001/09/07 19:58:49  zvi
! Starting new code developement
!
! Revision 1.10  2001/06/04 22:43:02  livesey
! Added belowRef
!
! Revision 1.9  2001/05/03 23:07:12  livesey
! Added scan model stuff
!
! Revision 1.8  2001/05/01 00:23:00  livesey
! Nullified fmStat%rows by default
!
! Revision 1.7  2001/04/28 07:05:54  livesey
! Tidy up destruction of ForwardModelIntermediate, for cases where
! forward model has not been run.
!
! Revision 1.6  2001/04/23 21:40:46  livesey
! Moved closestInstances into ifm
!
! Revision 1.5  2001/04/19 23:57:00  livesey
! New fmStat
!
! Revision 1.4  2001/04/19 21:04:42  livesey
! Bug fix, sized arrays wrongly
!
! Revision 1.3  2001/04/19 20:29:46  livesey
! Added destroy routine
!
! Revision 1.2  2001/04/10 23:16:14  livesey
! Working version.
!
! Revision 1.1  2001/04/10 22:17:05  livesey
! Renamed module
!
! Revision 1.1  2001/04/10 22:06:26  livesey
! Very first version.  Bare bones only.
!
