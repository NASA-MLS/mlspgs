module ERROR_HANDLER

  use MACHINE, only: IO_ERROR
  use OUTPUT_M, only: OUTPUT

  implicit NONE
  public

  ! Values of the WHY argument of ERROR_INTRO
  integer, parameter :: COMPILER = 1
  integer, parameter :: FATAL = 2
  integer, parameter :: SIMPLE = 3
  integer, parameter :: WARNING = 4

  integer, private :: LAST_ERROR = 0

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
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
    select case ( why )
    case ( compiler ); call output ( 'Compiler error' )
    case ( fatal );    call output ( 'Fatal error' )
    case ( simple );   call output ( 'Error' )
    case ( warning );  call output ( 'Warning' )
    end select
    call output (' *** ')
    if ( present(where) ) then
      call output ('near column ')
      write ( line, '(i6)' ) mod(where,256)
      call output ( line(index( line, ' ', back=.true. )+1:) )
      call output ( ' of line ' )
      write ( line, '(i6)' ) where / 256
      call output ( line(index( line, ' ', back=.true. )+1:) )
    end if
    call output ( ': ' )
    return
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
    return
  end subroutine ERROR_WALKBACK

  subroutine FINISH_ERROR_WALKBACK
    character(len=6) :: LINE
    if ( last_error /= 0 ) then
      call output ( 'Last error at line ')
      write ( line, '(i6)' ) last_error
      call output ( line(index( line, ' ', back=.true. )+1:), &
                          advance='yes' )
    end if
    return
  end subroutine FINISH_ERROR_WALKBACK

  subroutine INIT_ERROR_WALKBACK
    last_error = 0
  end subroutine INIT_ERROR_WALKBACK
end module ERROR_HANDLER

! $Log$
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
