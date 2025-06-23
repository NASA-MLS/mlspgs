! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MorePVM                          ! Additional PVM stuff

  ! This module contains routines for packing and unpacking integers that are
  ! to be interpreted as lit indices or string indices

  implicit none
  private

  public :: PVMPackLitIndex, PVMPackStringIndex, &
    & PVMUnpackLitIndex, PVMUnpackStringIndex

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ------------------------------------------ PVMPackLitIndex ---
  subroutine PVMPackLitIndex ( index, info, msg )
    use Intrinsic, only: Lit_indices
    ! Dummy arguments
    integer, intent(in) :: INDEX                  ! lit index
    integer, intent(out), optional :: INFO        ! Flag from pvm
    character(len=*), intent(in), optional :: MSG ! In case of error

    ! Executable code
    call PVMPackStringIndex ( index, info, msg, lit_indices )
  end subroutine PVMPackLitIndex

  ! ------------------------------------------ PVMPackStringIndex --
  subroutine PVMPackStringIndex ( index, info, msg, array )
    use PVMIDL, only: PVMIDLPack
    use String_table, only: Get_string
    ! Dummy arguments
    integer, intent(in) :: INDEX                  ! String index
    integer, intent(out), optional :: INFO        ! Flag from pvm
    character(len=*), intent(in), optional :: MSG ! In case of error
    integer, intent(in), optional :: ARRAY(:)     ! Use ARRAY(INDEX)
    ! Local variables
    character (len=4096) :: WORD
    ! Executable code
    if ( index == 0 ) then
      call PVMIDLPack ( '', info, msg )
    else
      if ( present(array) ) then
        call get_string ( array(index), word, strip=.true., noError=.true. )
      else
        call get_string ( index, word, strip=.true., noError=.true. )
      end if
      call PVMIDLPack ( trim(word), info, msg )
    end if
  end subroutine PVMPackStringIndex

  ! ------------------------------------------ PVMUnpackLitIndex ---
  subroutine PVMUnpackLitIndex ( index, info, msg )
    use PVMIDL, only: PVMIDLUnpack
    use MoreTree, only: GetLitIndexFromString
    ! Dummy arguments
    integer, intent(inout) :: INDEX                 ! String index
    integer, intent(out), optional :: INFO        ! Flag from pvm
    character(len=*), intent(in), optional :: MSG ! In case of error
    ! Local variables
    character (len=4096) :: WORD

    ! Executable code
    word = ' ' ! Don't initialize on declaration, Intel compiler ignores it
    call PVMIDLUnpack ( word, info, msg )
    if ( info == 0 ) then
      if ( len_trim ( word ) > 0 ) then
        index = GetLitIndexFromString ( trim(word) )
      else
        index = 0
      end if
    end if
  end subroutine PVMUnpackLitIndex

  ! ------------------------------------------ PVMUnpackStringIndex ---
  subroutine PVMUnpackStringIndex ( index, info, outWord, msg )
    use PVMIDL, only: PVMIDLUnpack
    use MoreTree, only: GetStringIndexFromString
    ! Dummy arguments
    integer, intent(inout) :: INDEX                 ! String index
    integer, intent(out), optional :: INFO        ! Flag from pvm
    character (len=*), intent(out), optional :: OUTWORD
    character(len=*), intent(in), optional :: MSG ! In case of error
    ! Local variables
    character (len=4096) :: WORD
    ! Executable code
    call PVMIDLUnpack ( word, info, msg )
    if ( present(outWord) ) outWord=word
    if ( info == 0 ) then
      if ( len_trim ( word ) > 0 ) then
        index = GetStringIndexFromString ( trim(word), caseSensitive=.true. )
      else
        index = 0
      end if
    end if
  end subroutine PVMUnpackStringIndex

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MorePVM

! $Log$
! Revision 2.12  2011/07/12 22:37:17  honghanh
! Fix bug in PVMUnpackLitIndex to initialize variable 'word'
! in the execution code instead of declaration code
!
! Revision 2.11  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.10  2006/11/29 03:05:19  vsnyder
! Remove unused USE names
!
! Revision 2.9  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.8  2004/10/19 22:59:34  vsnyder
! Add optional 'msg' argument and internal error processing
!
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
