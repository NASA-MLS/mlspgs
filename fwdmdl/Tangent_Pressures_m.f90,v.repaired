! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Tangent_Pressures_m

  implicit NONE
  private
  public :: Tangent_Pressures

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

!---------------------------------------------  Tangent_Pressures  -----

  subroutine Tangent_Pressures ( FwdModelConf, Z_grid, No_Tan_Hts, &
    &                            SurfaceTangentIndex, Tan_Press )

    ! Compute SurfaceTangentIndex, and from that No_Tan_Hts.
    ! Compute Tan_Press from fwdModelConf%tangentGrid%surfs and z_grid

    use Allocate_Deallocate, only: Allocate_Test
    use ForwardModelConfig, only: ForwardModelConfig_t
    use MLSCommon, only: RP

  ! Inputs:
    type (forwardModelConfig_T), intent(in) :: FwdModelConf
    real(rp), intent(in) :: Z_Grid(:)

  ! Outputs:
    integer, intent(out) :: No_Tan_Hts          ! Number of tangent heights
    integer, intent(out) :: SurfaceTangentIndex ! First or surface tangent index
    real(rp), allocatable, intent(out) :: Tan_Press(:)  ! Pressures at tangent
                                                ! points in z_grid

    ! Allocate tan_press and compute it from fwdModelConf%tangentGrid%surfs
    ! and z_grid

    if ( associated(FwdModelConf%tangentGrid) ) then
      surfaceTangentIndex = COUNT( fwdModelConf%tangentGrid%surfs < &
                                 & (z_grid(1) - 0.0001_rp)) + 1
    else
      surfaceTangentIndex = 1
    end if
    no_tan_hts = size(z_grid) + surfaceTangentIndex - 1

    call allocate_test ( tan_press, no_tan_hts, 'tan_press', moduleName )

! Compute tan_press from fwdModelConf%tangentGrid%surfs and z_grid

    if ( associated(FwdModelConf%tangentGrid) ) &
      & tan_press(1:surfaceTangentIndex-1) = &
        & fwdModelConf%tangentGrid%surfs(1:surfaceTangentIndex-1,1)
    tan_press(surfaceTangentIndex:no_tan_hts) = z_grid

  end subroutine Tangent_Pressures

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Tangent_Pressures_m

! $Log$
! Revision 2.3  2016/09/03 00:31:06  vsnyder
! Make Tan_Press intent(out) to ensure it gets deallocated first
!
! Revision 2.2  2016/06/03 23:44:44  vsnyder
! Make Tan_Press allocatable instead of a pointer
!
! Revision 2.1  2013/07/11 00:00:01  vsnyder
! Initial commit
!
