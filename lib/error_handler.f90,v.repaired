! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ERROR_HANDLER

  use Machine, only: IO_Error
  use Output_m, only: Output

  implicit NONE
  public

  ! Values of the WHY argument of ERROR_INTRO
  integer, parameter :: Compiler = 1
  integer, parameter :: Fatal = 2
  integer, parameter :: Simple = 3
  integer, parameter :: Warning = 4

  integer, private :: Last_Error = 0

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine ALLOCATE_STAT ( DOING_WHAT, STAT )
  ! Examine the allocation status.  If it failed, print an error message
  ! and stop
    character(len=*), intent(in) :: DOING_WHAT    ! What/where doing allocate?
    integer, intent(in) :: STAT                   ! From ALLOCATE
    if ( stat == 0 ) return
    call io_error ( 'ALLOCATE Failed while allocating:', stat, doing_what )
    stop
  end subroutine ALLOCATE_STAT

  subroutine ERROR_INTRO ( WHY, WHERE )
    integer, intent(in) :: WHY
    integer, intent(in), optional :: WHERE   ! 256*line + column
    character(len=6) :: LINE
    call output ('*** ')
    if ( why > 0 ) then
      select case ( why )
      case ( compiler ); call output ( 'Compiler error' )
      case ( fatal );    call output ( 'Fatal error' )
      case ( simple );   call output ( 'Error' )
      case ( warning );  call output ( 'Warning' )
      end select
      call output (' *** ')
    end if
    if ( present(where) ) then
      call output ('Near column ')
      write ( line, '(i6)' ) mod(where,256)
      call output ( line(index( line, ' ', back=.true. )+1:) )
      call output ( ' of line ' )
      write ( line, '(i6)' ) where / 256
      call output ( line(index( line, ' ', back=.true. )+1:) )
    end if
    call output ( ': ' )
  end subroutine ERROR_INTRO

  subroutine ERROR_WALKBACK ( LINE_NO )
    integer, intent(in) :: LINE_NO
    character(len=6) :: LINE
    if ( last_error /= 0 .and. last_error /= line_no ) then
      call output ( 'Previous error at line ')
      write ( line, '(i6)' ) last_error
      call output ( line(index( line, ' ', back=.true. )+1:), &
                          advance='yes' )
    end if
    last_error = line_no
  end subroutine ERROR_WALKBACK

  subroutine FINISH_ERROR_WALKBACK
    character(len=6) :: LINE
    if ( last_error /= 0 ) then
      call output ( 'Last error at line ')
      write ( line, '(i6)' ) last_error
      call output ( line(index( line, ' ', back=.true. )+1:), &
                          advance='yes' )
    end if
  end subroutine FINISH_ERROR_WALKBACK

  subroutine INIT_ERROR_WALKBACK
    last_error = 0
  end subroutine INIT_ERROR_WALKBACK

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ERROR_HANDLER

! $Log$
! Revision 2.5  2016/08/31 01:42:06  vsnyder
! Don't print Why in Error_intro if it's zero
!
! Revision 2.4  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.3  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2002/10/08 00:09:09  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2000/10/11 18:33:24  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:49  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
