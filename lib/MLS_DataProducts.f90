! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MLS_DataProducts
  use MLSCommon, only: r8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_deallocate
  implicit none
  public :: Deallocate_DataProducts, DataProducts_T
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  type DataProducts_T
     character(len=80) :: name, data_type
     character(len=20), dimension(:), pointer :: Dimensions => NULL()
     real(r8), dimension(:), pointer :: FrequencyCoordinates  => NULL()
     real(r8), dimension(:), pointer :: VerticalCoordinates   => NULL()
     real(r8), dimension(:), pointer :: HorizontalCoordinates => NULL()
  end type DataProducts_T

contains

 subroutine Deallocate_DataProducts( DataProducts )
    type( DataProducts_T ), intent(inout) :: DataProducts
    character(len=480) :: msr
    integer :: status

    if (associated(DataProducts%FrequencyCoordinates)) then 
        deallocate(DataProducts%FrequencyCoordinates, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' DataProducts%FrequencyCoordinates in ' // &
           DataProducts%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if (associated(DataProducts%VerticalCoordinates)) then 
        deallocate(DataProducts%VerticalCoordinates, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' DataProducts%VerticalCoordinates in ' // &
           DataProducts%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if (associated(DataProducts%HorizontalCoordinates)) then 
        deallocate(DataProducts%HorizontalCoordinates, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' DataProducts%HorizontalCoordinates in ' // &
           DataProducts%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if (associated(DataProducts%Dimensions)) then 
        deallocate(DataProducts%Dimensions, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' DataProducts%Dimensions in ' // &
           DataProducts%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

 end subroutine Deallocate_DataProducts

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------
end module MLS_DataProducts
