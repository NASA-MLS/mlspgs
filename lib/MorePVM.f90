! Copyright (c 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MorePVM                          ! Additional PVM stuff

  ! This module contains routines for packing and unpacking integers that are
  ! to be interpreted as lit indices or string indices

  implicit none
  private

  public :: PVMPackLitIndex, PVMPackStringIndex, &
    & PVMUnpackLitIndex, PVMUnpackStringIndex

  !------------------------------- RCS Ident Info ------------------------------
  character(LEN=130), private :: Id = & 
       "$Id$"
  character(LEN=*), private, parameter :: ModuleName="$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

contains

  ! ------------------------------------------ PVMPackLitIndex ---
  subroutine PVMPackLitIndex ( index, info )
    use Intrinsic, only: Lit_indices
    use PVMIDL, only: PVMIDLPack
    ! Dummy arguments
    integer, intent(in) :: INDEX        ! lit index
    integer, intent(out) :: INFO        ! Flag from pvm

    ! Executable code
    if ( index == 0 ) then
      call PVMIDLPack ( '', info )
    else
      call PVMPackStringIndex ( lit_indices ( index ), info )
    end if
  end subroutine PVMPackLitIndex

  ! ------------------------------------------ PVMPackStringIndex --
  subroutine PVMPackStringIndex ( index, info )
    use PVMIDL, only: PVMIDLPack
    use String_table, only: Get_string
    ! Dummy arguments
    integer, intent(in) :: INDEX        ! String index
    integer, intent(out) :: INFO        ! Flag from pvm
    ! Local variables
    character (len=4096) :: WORD
    ! Executable code
    if ( index == 0 ) then
      call PVMIDLPack ( '', info )
    else
      call get_string ( index, word, strip=.true., noError=.true. )
      call PVMIDLPack ( trim(word), info )
    end if
  end subroutine PVMPackStringIndex

  ! ------------------------------------------ PVMUnpackLitIndex ---
  subroutine PVMUnpackLitIndex ( index, info )
    use PVMIDL, only: PVMIDLUnpack
    use MoreTree, only: GetLitIndexFromString
    ! Dummy arguments
    integer, intent(out) :: INDEX        ! String index
    integer, intent(out) :: INFO        ! Flag from pvm
    ! Local variables
    character (len=4096) :: WORD
    ! Executable code
    call PVMIDLUnpack ( word, info )
    if ( info == 0 ) then
      if ( len_trim ( word ) > 0 ) then
        index = GetLitIndexFromString ( trim(word) )
      else
        index = 0
      end if
    end if
  end subroutine PVMUnpackLitIndex

  ! ------------------------------------------ PVMUnpackStringIndex ---
  subroutine PVMUnpackStringIndex ( index, info, outWord )
    use PVMIDL, only: PVMIDLUnpack
    use MoreTree, only: GetStringIndexFromString
    ! Dummy arguments
    integer, intent(out) :: INDEX        ! String index
    integer, intent(out) :: INFO        ! Flag from pvm
    character (len=*), intent(out), optional :: OUTWORD
    ! Local variables
    character (len=4096) :: WORD
    ! Executable code
    call PVMIDLUnpack ( word, info )
    ! print *, 'word in PVMUnpackStringIndex: ', trim(word)
    if ( present(outWord) ) outWord=word
    if ( info == 0 ) then
      if ( len_trim ( word ) > 0 ) then
        index = GetStringIndexFromString ( trim(word), caseSensitive=.true. )
      else
        index = 0
      end if
    end if
  end subroutine PVMUnpackStringIndex

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MorePVM

! $Log$
! Revision 2.7  2004/01/22 00:46:35  pwagner
! Enforces case sensitivity in PVMUnpackStringIndex; optional outWord from PVMUnpackStringIndex
!
! Revision 2.6  2003/09/15 23:15:28  vsnyder
! Remove unused USE for PVMErrorMessage
!
! Revision 2.5  2003/07/07 20:21:17  livesey
! Increased length of string literals
!
! Revision 2.4  2002/12/04 21:54:12  livesey
! Minor bug fixes with string stuff
!
! Revision 2.3  2002/10/08 17:42:46  livesey
! Bug fix
!
! Revision 2.2  2002/10/08 00:09:12  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2002/10/05 00:42:04  livesey
! First version
!
