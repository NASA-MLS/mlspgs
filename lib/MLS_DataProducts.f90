module MLS_DataProducts
  use MLSCommon, only: r8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_deallocate
  implicit none
  public :: Deallocate_DataProducts, DataProducts_T
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"

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
end module MLS_DataProducts
