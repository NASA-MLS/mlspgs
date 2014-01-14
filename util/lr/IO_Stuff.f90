! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module IO_Stuff

  ! Stuff related to I/O: Unit numbers for listing output and table output,
  ! and a subroutine to get an unused unit number

  use, intrinsic :: ISO_Fortran_Env, only: Output_Unit
  implicit NONE
  public

  integer :: List_Unit = Output_Unit ! Defaults to standard output
  integer :: Table_Unit = -1         ! Defaults to no table output

  integer, private, parameter :: Bottom_Unit_Num = 10, Top_Unit_Num = 99

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! ================================================     GET_LUN     =====

  subroutine GET_LUN ( LUN, MSG, BOTTOM, TOP )
    ! Find a Fortran logical unit number that's not in use.
    ! Args
    integer, intent(out)          :: LUN    ! The logical unit number
    logical, intent(in), optional :: MSG ! Print failure message? (default: T)
    integer, intent(in), optional :: BOTTOM ! Where to begin
    integer, intent(in), optional :: TOP    ! Where to end
    ! Internal variables
    logical :: EXIST, OPENED             ! Used to inquire about the unit
    integer :: myBottom
    integer :: myTop
    ! Executable
    myBottom = bottom_unit_num
    myTop    = top_unit_num
    if ( present(Bottom) ) myBottom = Bottom
    if ( present(Top) ) myTop = Top
    do lun = myBottom, myTop
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) return
    end do
    lun = -1
    if ( present(msg) ) then
      if ( .not. msg ) return
    end if
    write(*,*) 'IO_STUFF%GET_LUN-E- Unable to get a logical unit number'
    return
  end subroutine GET_LUN

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module IO_Stuff

! $Log$
