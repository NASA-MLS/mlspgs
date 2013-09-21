! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module TRACE_M

  implicit none
  private

  public :: Trace_Begin, Trace_End

  interface Trace_Begin
    module procedure Trace_Begin_B, Trace_Begin_C, Trace_Begin_I
  end interface

  character (len=8), save :: PreviousDate = ' '
  logical, save           :: PreviousDebug
  logical, save           :: PreviousVerbose
  logical, parameter      :: DEEBUG = .false.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

! ------------------------------------------------  TRACE_BEGIN_B  -----
  subroutine TRACE_BEGIN_B ( NAME_I, NAME_C, ROOT, INDEX, String, Cond, Advance )
  ! If Name_I <= 0, use Create_String ( Name_C ) to give it a value.
  ! We assume the actual argument is a SAVE variable.  Thereby, if
  ! Name_I is positive, we assume it's the result of entering Name_C,
  ! and it isn't done again.

  ! Print "ENTER NAME with ROOT = <node_id(root)>" with DEPTH dots in
  ! front.  Increment DEPTH.

    use STRING_TABLE, only: CREATE_STRING, HOW_MANY_STRINGS

    integer, intent(inout) :: NAME_I
    character(len=*), intent(in) :: NAME_C
    integer, intent(in), optional :: ROOT
    integer, intent(in), optional :: INDEX
    integer, intent(in), optional :: String ! To display after Name_I
    logical, intent(in), optional :: Cond   ! Print if true, default true
    character(len=*), intent(in), optional :: Advance
    if ( how_many_strings() < 1 ) return
    if ( name_i <= 0 ) name_i = create_string ( trim(name_c) )
    call trace_begin ( name_i, root, index, string, cond, advance )

  end subroutine TRACE_BEGIN_B

! ------------------------------------------------  TRACE_BEGIN_C  -----
  subroutine TRACE_BEGIN_C ( NAME_C, ROOT, INDEX, String, Cond, Advance )
  ! Print "ENTER NAME with ROOT = <node_id(root)>" with DEPTH dots in
  ! front.  Increment DEPTH.

    use STRING_TABLE, only: CREATE_STRING, HOW_MANY_STRINGS

    character(len=*), intent(in) :: NAME_C
    integer, intent(in), optional :: ROOT
    integer, intent(in), optional :: INDEX
    integer, intent(in), optional :: String ! To display after Name_C
    logical, intent(in), optional :: Cond   ! Print if true, default true
    character(len=*), intent(in), optional :: Advance

    integer :: Name_I

    if ( how_many_strings() < 1 ) return
    name_i = create_string ( trim(name_c) )
    call trace_begin ( name_i, root, index, string, cond )

  end subroutine TRACE_BEGIN_C

! ------------------------------------------------  TRACE_BEGIN_I  -----
  subroutine TRACE_BEGIN_I ( NAME, ROOT, INDEX, String, Cond, Advance )
  ! Print "ENTER NAME with ROOT = <node_id(root)>" with DEPTH dots in
  ! front.  Increment DEPTH.

    use CALL_STACK_M, only: PUSH_STACK
    use MLSCOMMON, only: MLSDEBUG, MLSVERBOSE, MLSNAMESAREDEBUG, MLSNAMESAREVERBOSE
    use MLSMESSAGEMODULE, only: MLSMESSAGECALLS
    use OUTPUT_M, only: OUTPUTOPTIONS
    use STRING_TABLE, only: CREATE_STRING, GET_STRING, HOW_MANY_STRINGS

    integer, intent(in) :: NAME
    integer, intent(in), optional :: ROOT
    integer, intent(in), optional :: INDEX
    integer, intent(in), optional :: String ! To display after Name_I
    logical, intent(in), optional :: Cond   ! Print if true, default true
    character(len=*), intent(in), optional :: Advance

    logical :: MyCond

    ! Executable
 
    if ( how_many_strings() < 1 ) return
    call Checkdate

    myCond = .true.
    if ( present(cond) ) myCond = cond

    if ( myCond ) then
      call push_stack ( name, root, index, string, before='Enter ', &
        & where=.true., advance=advance )
    else
      call push_stack ( name, root, index, string )
    end if

    if ( create_string(MLSNamesAreDebug) == name .and. .not. MLSDebug ) then
      PreviousDebug = MLSDebug
      MLSDebug = .true.
    end if
    if ( create_string(MLSNamesAreVerbose) == name .and. .not. MLSVerbose ) then
      PreviousVerbose = MLSVerbose
      MLSVerbose = .true.
    end if

    call get_string ( name, outputOptions%parentName )
    call MLSMessageCalls( 'push', name=outputOptions%parentName )

  end subroutine TRACE_BEGIN_I

! --------------------------------------------------    TRACE_END  -----
  subroutine TRACE_END ( NAME, INDEX, Cond )
  ! Decrement DEPTH.  Print "EXIT NAME" with DEPTH dots in front.

  ! The only reason to provide Name is if you want to check whether stack
  ! pushes and pops match, by setting -Schktr.  You probably don't need to
  ! bother with this if you can see both Trace_Begin and Trace_end, and
  ! they're invoked with the same condition.

    use CALL_STACK_M, only: POP_STACK, STACK_T, STACK_DEPTH, TOP_STACK
    use MLSCOMMON, only: MLSDEBUG, MLSVERBOSE, MLSNAMESAREDEBUG, MLSNAMESAREVERBOSE
    use MLSMESSAGEMODULE, only: MLSMESSAGECALLS
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use OUTPUT_M, only: NEWLINE, OUTPUT, OUTPUTOPTIONS
    use STRING_TABLE, only: CREATE_STRING, DISPLAY_STRING, HOW_MANY_STRINGS !, Get_String
    use TOGGLES, only: SWITCHES

    character(len=*), optional, intent(in) :: NAME ! Checked but taken from stack
    integer, intent(in), optional :: INDEX ! Use value from stack if not present
    logical, intent(in), optional :: Cond  ! Print if true, default true

    integer :: Check = -2
    type(stack_t) :: Frame
    character(32) :: PARENTNAME

    logical :: MyCond

    ! Executable

    if ( how_many_strings() < 1 ) return
    call Checkdate

    if ( check < -1 ) check = switchDetail ( switches, 'chktr' )

    myCond = merge(.true., .false., present(name))
    if ( present(cond) ) myCond = cond
    if ( myCond ) then
      call pop_stack ( 'Exit ', .true., frame=frame, index=index )
    else
      call pop_stack ( frame=frame, index=index, Silent = .true. )
    end if

    if ( check > -1 ) then
      if ( frame%tree < 0 .or. stack_depth() <= 0 ) then
        call output ( 'In Trace_End, stack underflow noticed ' )
        if ( present(name) ) call output ( 'with NAME = ' // trim(name) )
        if ( present(index) ) call output ( index, before=' INDEX = ' )
        call newLine
      else
        if ( present(name) ) then
          if ( frame%text /= create_string(trim(name)) ) then
            call display_string ( frame%text, &
              & before='In Trace_End, name at top of stack = "' )
            call output ( '" but NAME = "' // trim(name) // '"', advance='yes' )
          end if
        end if
        if ( present(index) ) then
          if ( frame%index /= index ) then
            call display_string ( frame%text, &
              & before='In Trace_End, in frame for "' )
            call output ( frame%index, before='", INDEX at top of stack = ' )
            call output ( index, before=' but INDEX argument = ', advance='yes' )
          end if
        end if
      end if
    end if

    if ( create_string(trim(MLSNamesAreDebug)) == frame%text .and. .not. MLSDebug ) then
      PreviousDebug = MLSDebug
      MLSDebug = .true.
    end if
    if ( create_string(trim(MLSNamesAreVerbose)) == frame%text .and. .not. MLSVerbose ) then
      PreviousVerbose = MLSVerbose
      MLSVerbose = .true.
    end if

!   call top_stack ( frame )
!   call get_string ( frame%text, outputOptions%parentName )
    call MLSMessageCalls( 'pop' )
    call MLSMessageCalls( 'top', outputOptions%parentName )  

  end subroutine TRACE_END

!------------------------ private procedures ---------------------
  ! -- CheckDate
  ! Checks for a jump in the date --
  ! So that we don't accidentally drop (multiples of) 24 hours in 
  ! marking times entering and exiting Named process
  subroutine CheckDate
    character (len=8) :: CurrentDate
    ! Executable
    call date_and_time( date=CurrentDate )
    if ( len_trim(PreviousDate) > 0 ) then
      if ( PreviousDate /= CurrentDate .or. DEEBUG ) &
        & call DateLine ( PreviousDate, CurrentDate )
    endif
    PreviousDate = CurrentDate
  end subroutine CheckDate
  
  subroutine DateLine( olddate, newdate )
    use OUTPUT_M, only: BLANKS, NEWLINE, OUTPUT
    character(len=*), intent(in) :: olddate, newdate
    ! Executable
    call blanks( 15, fillChar='-' )
    call output( ' Crossed date boundary from ' )
    call output( olddate )
    call output( ' to ' )
    call output( newdate )
    call blanks( 1 )
    call blanks( 15, fillChar='-' )
    call newline
  end subroutine DateLine
  
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module TRACE_M

! $Log$
! Revision 2.32  2013/09/21 00:21:39  pwagner
! Check if string table initialized first
!
! Revision 2.31  2013/09/12 03:12:13  vsnyder
! Add Advance to Trace_Begin, pass Index through Trace_End to Pop_Stack
!
! Revision 2.30  2013/09/04 02:49:35  vsnyder
! Simplify stack handling in Trace_End
!
! Revision 2.29  2013/08/31 02:26:21  vsnyder
! Better debugging output
!
! Revision 2.28  2013/08/31 02:01:51  vsnyder
! Always trim names sent to string table
!
! Revision 2.27  2013/08/31 01:24:01  vsnyder
! Better debugging output
!
! Revision 2.26  2013/08/30 23:13:20  pwagner
! Prevent unwanted printing during routine trace_end
!
! Revision 2.25  2013/08/30 03:55:49  vsnyder
! Add String argument to Trace_Begin
!
! Revision 2.24  2013/08/23 02:50:47  vsnyder
! Use Call_Stack
!
! Revision 2.23  2013/08/17 03:10:13  vsnyder
! Use Call_Stack
!
! Revision 2.22  2013/06/28 18:06:14  pwagner
! Automatically set parentName
!
! Revision 2.21  2012/05/24 01:35:59  vsnyder
! Clean up output
!
! Revision 2.20  2011/12/13 01:09:00  pwagner
! Automatically turns on MLSVerbose or MLSDebug according to MLSNamesAre..
!
! Revision 2.19  2011/10/14 00:33:17  pwagner
! Prints notice if date boundary crossed
!
! Revision 2.18  2011/03/31 19:55:22  vsnyder
! Cannonball polishing
!
! Revision 2.17  2009/06/23 18:25:44  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.16  2007/08/13 17:38:42  pwagner
! Push named procedures automatically onto new MLSCallStack
!
! Revision 2.15  2007/07/27 00:20:51  vsnyder
! Spiff up printing, work on memory tracking
!
! Revision 2.14  2006/07/28 01:59:42  vsnyder
! Correct bug in memory reporting, plus cannonball polishing
!
! Revision 2.13  2006/07/19 22:26:04  vsnyder
! Report memory size changes in trace_end, plus some cannonball polishing
!
! Revision 2.12  2005/06/22 17:25:51  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.11  2002/10/08 00:09:15  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.10  2001/11/09 23:14:08  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.9  2001/09/13 19:36:50  livesey
! Added optional index arguments
!
! Revision 2.8  2001/05/03 02:13:25  vsnyder
! Trying to print time exclusive of calls isn't working
!
! Revision 2.7  2001/05/01 23:53:40  vsnyder
! Print CPU time exclusive of deeper ones at each end_trace
!
! Revision 2.6  2001/04/25 00:08:26  vsnyder
! Use 'fill' argument of 'output' to get leading zeroes on milliseconds
!
! Revision 2.5  2001/04/24 23:33:56  vsnyder
! Add timing
!
! Revision 2.4  2001/04/24 20:12:47  vsnyder
! Emit blanks before 'Enter' if no tree node supplied
!
! Revision 2.3  2001/04/17 20:56:10  vsnyder
! Emit 'depth' dots in trace_begin even if 'root' is absent
!
! Revision 2.2  2001/03/16 21:01:17  vsnyder
! Make ROOT optional in Trace_Begin
!
! Revision 2.1  2000/10/11 18:33:25  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
