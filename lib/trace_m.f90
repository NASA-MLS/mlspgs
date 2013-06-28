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
  use MLSCOMMON, only: MLSDEBUG, MLSVERBOSE, MLSNAMESAREDEBUG, MLSNAMESAREVERBOSE
  use MLSSTRINGLISTS, only: SWITCHDETAIL
  implicit none
  private
  public :: Trace_Begin, Trace_End
  integer, public, save :: DEPTH   ! Depth in tree.  Used for trace printing.

  integer, parameter :: ClockStackMax = 100
  real :: ClockStack(0:clockStackMax) = 0.0
  double precision :: Memory(0:clockStackMax) = 0.0d0
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
! --------------------------------------------------  TRACE_BEGIN  -----
  subroutine TRACE_BEGIN ( NAME, ROOT, INDEX )
  ! Print "ENTER NAME with ROOT = <node_id(root)>" with DEPTH dots in
  ! front.  Increment DEPTH.

    use ALLOCATE_DEALLOCATE, only: MEMORY_UNITS, NOBYTESALLOCATED
    use MLSMESSAGEMODULE, only: MLSMESSAGECALLS
    use OUTPUT_M, only: OUTPUTOPTIONS, DUMPSIZE, OUTPUT
    use TIME_M, only: TIME_NOW
    use TREE, only: DUMP_TREE_NODE

    character(len=*), intent(in) :: NAME
    integer, intent(in), optional :: ROOT
    integer, intent(in), optional :: INDEX
    integer :: I              ! Loop inductor
    character(len=10) :: Now  ! For Date_and_time
    ! Executable
    call Checkdate
    if ( depth < 0 ) call output ( depth, before='***** Why is depth = ', &
      & after=' negative?', advance='yes' )
    if ( present(root) ) then
      call output ( root, 6 ); call output ( ': ' )
    else
      call output ( '        ' )
    end if
    do i = 1, depth
      call output ( '.' )
    end do
    call output ( 'Enter ' ); call output ( name )
    if ( present(index) ) call output( index, before=' ' )
    call date_and_time ( time=now )
    call output ( ' at ' // now(1:2) // ':' // now(3:4) // ':' // now(5:) )
    if ( present(root) ) then
      call output ( ' with ' );
      call dump_tree_node ( root, 0 )
    else
      call output ( ' ' )
    end if
    if ( depth >= 0 .and. depth < clockStackMax ) then
      call time_now ( clockStack(depth) )
      clockStack(depth+1) = 0.0
      memory(depth) = NoBytesAllocated
    end if
    call dumpsize ( memory_units * nobytesallocated, before = 'Memory ', &
      & advance='yes' )
    depth = depth + 1
    call MLSMessageCalls( 'push', constantName=Name )
    outputOptions%parentName = Name
    if ( switchDetail( MLSNamesAreDebug, Name, options='-fc' ) > -1 .and. &
      & .not. MLSDebug ) then
      PreviousDebug = MLSDebug
      MLSDebug = .true.
    endif
    if ( switchDetail( MLSNamesAreVerbose, Name, options='-fc' ) > -1 .and. &
      & .not. MLSVerbose ) then
      PreviousVerbose = MLSVerbose
      MLSVerbose = .true.
    endif
  end subroutine TRACE_BEGIN
! --------------------------------------------------    TRACE_END  -----
  subroutine TRACE_END ( NAME, INDEX )
  ! Decrement DEPTH.  Print "EXIT NAME" with DEPTH dots in front.

    use ALLOCATE_DEALLOCATE, only: MEMORY_UNITS, NOBYTESALLOCATED
    use MLSMESSAGEMODULE, only: MLSMESSAGECALLS
    use OUTPUT_M, only: OUTPUTOPTIONS, DUMPSIZE, NEWLINE, OUTPUT
    use TIME_M, only: TIME_NOW

    character(len=*), intent(in) :: NAME
    integer, intent(in), optional :: INDEX
    double precision :: Delta ! memory
    integer :: I              ! Loop inductor
    character(len=10) :: Now  ! For Date_and_time
    character(32) :: PARENTNAME
    integer :: Values(8)      ! For Date_and_time
    real :: T                 ! For timing
    character(12) :: Used     ! For timing
    ! Executable
    call Checkdate
    depth = depth - 1
    if ( depth < 0 ) call output ( depth, before='***** Why is depth = ', &
      & after=' negative?', advance='yes' )
    call output ( '        ' )
    do i = 1, depth
      call output ( '.' )
    end do
    call output ( 'Exit ' ); call output ( name )
    if ( present(index) ) call output( index, before=' ' )
    call date_and_time ( time=now, values=values )
    call output ( ' at ' // now(1:2) // ':' // now(3:4) // ':' // now(5:) )
    if ( depth >= 0 .and. depth < clockStackMax ) then
      call time_now ( t )
      clockStack(depth) = t - clockStack(depth)
!     call output ( clockStack(depth) - clockStack(depth+1), &
!       & format='(g10.3)', before=' used ' )
      write ( used, '(g10.3)' ) clockStack(depth)
      call output ( ' used ' // trim(adjustl(used)) // ' cpu' )
      if ( memory(depth) /= noBytesAllocated ) then
        delta = memory_units * (noBytesAllocated-memory(depth))
        if ( abs(delta) < huge(1) ) then
          call dumpSize ( int(delta), before=', Memory changed by ' )
        else
          call dumpSize ( delta, before=', Memory changed by ' )
        end if
        call dumpSize ( memory_units * nobytesallocated, before = ' to ' )
        memory(depth) = nobytesallocated
      end if
    end if
    call newLine
    call MLSMessageCalls( 'pop' )
    call MLSMessageCalls( 'top', parentName )
    outputOptions%parentName = parentName
    if ( switchDetail( MLSNamesAreDebug, Name, options='-fc' ) > -1 ) then
      MLSDebug = PreviousDebug
    endif
    if ( switchDetail( MLSNamesAreVerbose, Name, options='-fc' ) > -1 ) then
      MLSVerbose = PreviousVerbose
    endif
  end subroutine TRACE_END

!------------------------ private procedures ---------------------
  ! -- CheckDate
  ! Checks for a jump in the date --
  ! So that we don't accidentally drop (multiples of) 24 hours in 
  ! marking times entering and exiting q Naamed process
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
