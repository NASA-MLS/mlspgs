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
module memory_m
!=============================================================================

! Calculate the memory usage. Return stack, heap, or total

  use IO_stuff, only: get_nLines, read_textFile
  use MLSFinds, only: FindFirst
  use MLSStringLists, only: GetStringElement
  use MLSStrings, only: readNumFromBaseN
  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! memory_used            return stack, heap, and sum

!     Example:
!     Assume you want to track how your memory usage rises and falls over time
!     If you know your process_id, the linux system records this info
!     in up-gto-date state in the file
!        /proc/process_id/maps
! === (end of toc) ===

! === (start of api) ===
! memory_used ( char* process_id, [real stack], [real heap], [real total] )
!
! Note:
! If process_id is not a valid process id, we will return 0
! === (end of api) ===
  public :: memory_used
  logical, parameter   :: countEmpty = .true.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine memory_used ( process_id, stack, heap, total )
    ! Calculate memory used the running process_id
    ! Args
    character(len=*), intent(in)                 :: process_id
    real, intent(out), optional                  :: stack
    real, intent(out), optional                  :: heap
    real, intent(out), optional                  :: total
    ! Internal variables
    character(len=256)                           :: filename
    character(len=256)                           :: line
    character(len=256), dimension(:), pointer    :: lines
    real                                         :: mstart, mfinish
    integer                                      :: n, kHeap, kStack, status
    ! Executable
    if ( present(stack) ) stack = 0.
    if ( present(heap) ) heap = 0.
    if ( present(total) ) total = 0.
    filename = '/proc/' // trim(process_id) // '/maps'
    call get_nLines ( filename, n )
    print *, 'n: ', n
    if ( n < 1 ) return
    allocate( lines(n), stat=status )
    call read_textfile ( FileName, lines )
    kHeap  = FindFirst( lines, '[heap]', 'p' )
    kStack = FindFirst( lines, '[stack]', 'p' )
    print *, 'kHeap: ', kHeap
    print *, 'kStack: ', kStack
    if ( kHeap > 0 ) then
      call getMemBounds ( lines(kHeap), mStart, mFinish )
      print *, 'heap'
      print *, 'mstart:  ', mstart
      print *, 'mfinish: ', mfinish
      print *,  mfinish -   mstart
      if ( present(heap) ) heap = heap + mFinish - mStart
      if ( present(total) ) total = total + mFinish - mStart
    endif
    if ( kStack > 0 ) then
      call getMemBounds ( lines(kStack), mStart, mFinish )
      print *, 'stack'
      print *, 'mstart:  ', mstart
      print *, 'mfinish: ', mfinish
      print *,  mfinish -   mstart
      if ( present(stack) ) stack = stack + mFinish - mStart
      if ( present(total) ) total = total + mFinish - mStart
    endif
    if ( kHeap > 1 ) then
      kHeap = kHeap - 1
      call getMemBounds ( lines(kHeap), mStart, mFinish )
      print *, 'small heap'
      print *, 'mstart:  ', mstart
      print *, 'mfinish: ', mfinish
      print *,  mfinish -   mstart
      if ( present(heap) ) heap = heap + mFinish - mStart
      if ( present(total) ) total = total + mFinish - mStart
    endif
    deallocate( lines, stat=status )
  end subroutine memory_used

  ! --------------- Private Procedures --------------------------
  ! Read memory start, finish values from line
  subroutine getMemBounds ( inLine, nStart, nFinish )
    ! Args
    character(len=*), intent(in)    :: inLine
    real, intent(out)               :: nstart, nfinish 
    ! Internal variables
    character(len=len(inLine))      :: line
    character(len=10)               :: start, finish   
    ! Executable
    call GetStringElement ( inLine, start, 1, countEmpty, inseparator = '-' )
    call GetStringElement ( inLine, line, 2, countEmpty, inseparator = '-' )
    print *, 'start: ', trim(start)
    print *, 'line: ', trim(line)
    call GetStringElement ( line, finish, 1, countEmpty, inseparator = ' ' )
    print *, 'finish: ', trim(finish)
    call readNumFromBaseN ( start, nstart, 16, options='c' )
    call readNumFromBaseN ( finish, nfinish, 16, options='c' )
  end subroutine getMemBounds

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module memory_m

!$Log$
!Revision 2.1  2014/08/05 00:19:02  pwagner
!First commit
!
