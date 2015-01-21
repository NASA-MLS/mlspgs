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
  use MLSKinds, only: r8
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
    use Allocate_Deallocate, only: Test_Deallocate
    type( DataProducts_T ), intent(inout) :: DataProducts
    integer :: s, status

    if (associated(DataProducts%FrequencyCoordinates)) then 
        s = size(DataProducts%FrequencyCoordinates) * &
          & storage_size(DataProducts%FrequencyCoordinates) / 8
        deallocate(DataProducts%FrequencyCoordinates, stat=status)
        call test_deallocate ( status, ModuleName, &
          & 'DataProducts%FrequencyCoordinates', s )
    endif

    if (associated(DataProducts%VerticalCoordinates)) then 
        s = size(DataProducts%VerticalCoordinates) * &
          & storage_size(DataProducts%VerticalCoordinates) / 8
        deallocate(DataProducts%VerticalCoordinates, stat=status)
        call test_deallocate ( status, ModuleName, &
          & 'DataProducts%VerticalCoordinates', s )
    endif

    if (associated(DataProducts%HorizontalCoordinates)) then 
        s = size(DataProducts%HorizontalCoordinates) * &
          & storage_size(DataProducts%HorizontalCoordinates) / 8
        deallocate(DataProducts%HorizontalCoordinates, stat=status)
        call test_deallocate ( status, ModuleName, &
          & 'DataProducts%HorizontalCoordinates', s )
    endif

    if (associated(DataProducts%Dimensions)) then 
        s = size(DataProducts%Dimensions) * &
          & storage_size(DataProducts%Dimensions) / 8
        deallocate(DataProducts%Dimensions, stat=status)
        call test_deallocate ( status, ModuleName, &
          & 'DataProducts%Dimensions', s )
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
