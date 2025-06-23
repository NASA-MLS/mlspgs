! Copyright 2014, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Memory_m
!=============================================================================

! Calculate the memory usage. Return stack, data, or total

  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! Memory_Used            return stack, data, and total (which includes exe)

!     Example:
!     Assume you want to track how your memory usage rises and falls over time
!     The linux system records this info in up-to-date state in the file
!        /proc/process_id/maps
! === (end of toc) ===

! === (start of api) ===
! Memory_Used ( [real *total], [real *stack], [real *data] )

! === (end of api) ===
  public :: Memory_Used

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Memory_Used ( Total, Stack, Data )
    ! Calculate memory usages in kB from /proc/GetPID()/status
    use GetResourceUsage_m, only: GetPID
    ! Args
    integer, intent(out), optional :: Total ! VmSize: ... kB in /proc/pid/status
                                            ! (everything, including exe)
    integer, intent(out), optional :: Stack ! VmStk: ... kB in /proc/pid/status
    integer, intent(out), optional :: Data  ! VmData: ... kB in /proc/pid/status
    ! Internal variables
    character(len=256)          :: Line
    integer                     :: L1, L2 ! Ends of numbers
    integer                     :: Unit ! to read /proc/pid/maps

    ! Executable

    ! Create the file name for /proc/[pid]/status
    write ( line, '(a,i0,a)' ) '/proc/', GetPID(), '/status'
    open ( newUnit = unit, file=trim(line) ) ! let it crash if it fails,
                                             ! which it almost certainly won't
    do
      read ( unit, '(a)', end=9 ) line
      ! Lines have "Name:     amount kB"
      l1 = index(line,':')                   ! First delimiter
      l2 = index(line, 'kB')                 ! Second delimiter
      select case ( line(:l1) )
      case ( 'VmSize:' )
        if ( present(total) ) read ( line(l1+1:l2-1), * ) total
      case ( 'VmData:' )
        if ( present(data) ) read ( line(l1+1:l2-1), * ) data
      case ( 'VmStk:' )
        if ( present(stack) ) read ( line(l1+1:l2-1), * ) stack
      end select
    end do
  9 continue
    close ( unit )

  end subroutine Memory_Used

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Memory_m

!$Log$
!Revision 2.3  2014/09/04 23:30:41  vsnyder
!New interface, use /proc/[pid]/status
!
!Revision 2.2  2014/08/15 23:28:49  pwagner
!Fixed some of the more obvious errors
!
!Revision 2.1  2014/08/05 00:19:02  pwagner
!First commit
!
