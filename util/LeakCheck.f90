program LeakCheck

  ! Within a single file, connect allocate statements or calls to
  ! Allocate_Test with deallocate statements or calls to Deallocate_Test,
  ! and report any mismatches.

  use Sort_M, only: AsortP

  implicit NONE

  integer, parameter :: MaxVars = 10000 ! Maximum number of variables to hold

  character(len=32) :: allocVar(maxVars), deAllocVar(maxVars)
  integer :: AllocLine(maxVars), deAllocLine(maxVars)
  integer :: AllocPVec(maxVars), deAllocPVec(maxVars)

  logical :: Debug = .false., Verbose = .false.

  integer :: I, J, K, L, LineNum = 0, NumAllocs = 0, NumDeallocs = 0

  character(132) :: Line
  character :: What ! A for allocate, D for deallocate

!---------------------------- RCS Module Info ------------------------------
  character (len=*), parameter :: ModuleName="$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------

  i = 0
  do
    i = i + 1
    call getarg ( i, line )
    l = len_trim(line)
    if ( line(1:1) == '-' ) then
      if ( verify(line(2:l),'dv') /= 0 ) then
        write ( *, * ) &
          & 'The only options I understand are d (debug) and v (verbose)'
        stop
      end if
      if ( index(line(2:l),'d') /= 0 ) debug = .true.
      if ( index(line(2:l),'v') /= 0 ) verbose = .true.
    else if ( line /= '' ) then
      write ( *, * ) &
        & 'The only command-line fields I understand are options beginning with -.'
      stop
    else
      exit
    end if
  end do

  do
    read ( *, '(a)', end=9 ) line
    lineNum = lineNum + 1
    line = adjustl(line)
    if ( line(1:1) == '!' ) cycle
    if ( line(1:1) == '&' ) line = adjustl(line(2:))
    do i = 1, len_trim(line)
      ! Lower-case the line
      if ( line(i:i) >= 'A' .and. line(i:i) <= 'Z' ) &
        & line(i:i) = achar(iachar(line(i:i)) + iachar('a') - iachar('A'))
    end do
    
    if ( line(1:5) == 'call ' ) then
      line = adjustl(line(6:))
      if ( line(1:13) == 'allocate_test' ) then
        line = adjustl(line(14:))
        what = 'A'
      else if ( line(1:15) == 'deallocate_test' ) then
        line = adjustl(line(16:))
        what = 'D'
      else
        cycle
      end if
      if ( line(1:1) /= '(' ) then
        write (*, *) 'Line ', lineNum, ' appears to have a syntax error'
        cycle
      end if
      line = adjustl(line(2:))
      i = index(line,',')
      if ( i /= 0 ) line(i:) = ''
    else
      if ( line(1:8) == 'allocate' ) then
        line = adjustl(line(9:))
        what = 'A'
      else if ( line(1:10) == 'deallocate' ) then
        line = adjustl(line(11:))
        what = 'D'
      else
        cycle
      end if
      if ( line(1:1) /= '(' ) then
        write (*, *) 'Line ', lineNum, ' appears to have a syntax error'
        cycle
      end if
      line = adjustl(line(2:))
      i = scan(line,'(,')
      if ( i /= 0 ) line(i:) = ''
    end if

    ! LINE now contains a variable name.  WHAT contains A or D.
    if ( what == 'A' ) then
      numAllocs = numAllocs + 1
      if ( numAllocs >= maxVars ) then
        write (*, *) 'Too many allocated variables (', numAllocs, '(.  Increase MaxVars.'
        stop
      end if
      allocVar(numAllocs) = line
      allocLine(numAllocs) = lineNum
    else
      numDeallocs = numDeallocs + 1
      if ( numDeallocs >= maxVars ) then
        write (*, *) 'Too many deallocated variables (', numDeallocs, '(.  Increase MaxVars.'
        stop
      end if
      deallocVar(numDeallocs) = line
      deallocLine(numDeallocs) = lineNum
    end if

  end do

9 continue

  call asortp ( allocVar, 1, numAllocs, allocPVec )
  call asortp ( deallocVar, 1, numDeallocs, deallocPVec )

  if ( debug ) then
    write ( *, '(a/(i5,1x,a))' ) 'Allocations:', &
      & (allocLine(allocPVec(i)), allocVar(allocPVec(i)), i = 1, numAllocs)
    write ( *, '(a/(i5,1x,a))' ) 'Deallocations:', &
      & (deallocLine(deallocPVec(i)), deallocVar(deallocPVec(i)), i = 1, numDeallocs)
  end if
    
  ! Insert sentinels
  allocVar(numAllocs+1)(1:1) = achar(iachar('z')+1)
  allocLine(numAllocs+1) = huge(0)
  allocPVec(numAllocs+1) = numAllocs+1
  deallocVar(numDeallocs+1)(1:1) = achar(iachar('z')+1)
  deallocLine(numDeallocs+1) = huge(0)
  deallocPVec(numDeallocs+1) = numDeallocs+1

  ! Merge allocations and deallocations.
  i = 1
  j = 1
  do
    if ( allocVar(allocPVec(i)) == deallocVar(deallocPVec(j)) ) then
      if ( i > numAllocs ) exit
      do k = i + 1, numAllocs
        if ( allocVar(allocPVec(k)) /= allocVar(allocPVec(i)) ) exit
      end do
      do l = j + 1, numDeallocs
        if ( deallocVar(deallocPVec(l)) /= deallocVar(deallocPVec(j)) ) exit
      end do
      if ( k - i >1 .or. l - j > 1 ) then
        write ( *, * ) trim(allocVar(allocPVec(i))), ' is allocated and deallocated several times'
        write ( *, * ) ' Allocations are at lines', allocLine(allocPVec(i:k-1))
        write ( *, * ) ' Deallocations are at lines', deallocLine(deallocPVec(j:l-1))
      else if ( verbose ) then
        write ( *, * ) trim(allocVar(allocPVec(i))), &
          & ' is allocated at line ', allocLine(allocPVec(i)), &
          & ' and deallocated at line ', deallocLine(deallocPVec(j))
      end if
      i = k
      j = l
    else if ( allocVar(allocPVec(i)) < deallocVar(deallocPVec(j)) ) then
      write ( *, * ) trim(allocVar(allocPVec(i))), ' is allocated at line ', &
        & allocLine(allocPVec(i)), ' but nowhere deallocated.'
      i = i + 1
    else
      write ( *, * ) trim(deallocVar(deallocPVec(j))), ' is deallocated at line ', &
        & deallocLine(deallocPVec(j)), ' but nowhere allocated.'
      j = j + 1
    end if
  end do

end program LeakCheck

! $Log$
! Revision 1.2  2005/11/19 02:54:40  vsnyder
! Print multiple allocate/deallocate even if equal numbers
!
! Revision 1.1  2005/11/19 02:05:11  vsnyder
! Initial commit
!
