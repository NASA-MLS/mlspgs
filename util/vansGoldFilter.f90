! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program FilterGold
  implicit none

! Filter output from gold brick comparison run
  logical :: AlreadyPrinted
  logical :: ReadyForIt
  character(255) :: Block, LastLine, Line, Quantity
  character(255) :: Hunt = 'MLS-Aura' ! What to look for to identify run
  character(255) :: WhatMaxAbsDiff    ! Block containing MaxMaxAbsDiff
  integer :: I, J
  real :: MaxAbsDiff, MaxAbsVal, MaxMaxAbsDiff = 0
  real :: RelMaxAbsDiff ! = maxAbsDiff / maxAbsVal
  real :: RelMaxAbsDiffPrint ! might be log10(relMaxAbsDiff)
  logical :: DoLog = .false., NeedHead = .true., SawOne = .true., Summary = .false.
  real :: Ref(2) ! Min, Max

!---------------------------- RCS Module Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
!---------------------------------------------------------------------------

  i = 0
  do
    i = i + 1
    call get_command_argument ( i, line )
    if ( line == '' ) exit
    if ( line(1:3) == '-a' ) then
      summary = .false.
    else if  ( line(1:3) == '-l' ) then
      doLog = .true.
    else if ( line(1:3) == '-s' ) then
      summary = .true.
    else if ( line(1:1) == '-' ) then
      call usage
    else
      hunt = line
    end if
  end do
  LastLine = ' '
  ReadyForIt = .false.
  do
    read ( *, '(a)', end=9 ) line
    if ( len_trim(line) < 1 ) cycle
    AlreadyPrinted = .false.
    i = index( line, 'Comparing' )
    if ( i == 0 ) i = min(1,  index ( line, trim(hunt) ))
    if ( i /= 0 ) then
      if ( .not. sawOne .and. ReadyForIt ) then
        call print_line ( 'File is identical to reference' )
        ReadyForIt = .false.
      end if
      if ( maxMaxAbsDiff > 0 ) then
        if ( doLog ) then
          maxMaxAbsDiff = log10(maxMaxAbsDiff)
          write ( *, 3 ) maxMaxAbsDiff, &
            & ' log10(maximum Rel Max Abs Val), in', adjustl(trim(whatMaxAbsDiff))
      3 format (f8.4,a37,t47,a)
        else
          write ( *, 4 ) maxMaxAbsDiff, &
            & ' is maximum Rel Max Abs Val, in', adjustl(trim(whatMaxAbsDiff))
        end if
      4 format (1p,g12.5,a33,t47,a)
      end if
      call print_line ( trim(line(i:)) )
      needHead = .true.
      sawOne = .false.
      maxMaxAbsDiff = 0
      whatMaxAbsDiff = ''
    end if
    ! Do we have special cases?
    if ( index( line, 'Gold brick summary' ) > 0 ) then
      call print_line (  trim(line) )
      call print_line (  ' ' )
      cycle
    elseif ( index( line, '----- Log files' ) > 0 ) then
      call print_line (  trim(line) )
      call print_line (  ' ' )
      cycle
    elseif ( index( line, 'identical to reference' ) > 0 .and. ReadyForIt ) then
      call print_line (  trim(line) )
      call print_line (  ' ' )
      ReadyForIt = .false.
      cycle
    elseif ( index( line, 'times longer than the reference' ) > 0 ) then
      call print_line (  trim(line) )
      cycle
    elseif ( index( line, 'crashed' ) > 0 ) then
      call print_line (  trim(line) )
      cycle
    endif
    ReadyForIt = ReadyForIt .or. index ( line, trim(hunt) ) > 0
    if ( index ( line, 'Block' ) /= 0 ) then
      i = index ( line, '[' )
      j = index ( line, ']' )
      if ( i /= 0 .and. j /= 0 ) then
        block = adjustl(line(i+1:j-1))
        ref = 0
        maxAbsDiff = 0
        maxAbsVal = 1
      end if
    end if
    if ( index ( line, 'Reference' ) /= 0 ) then
      i = index ( line, ':' )
      if ( i /= 0 ) then
        line = line(i+1:)
        i = index ( line, ':' )
        if ( i /= 0 ) line = line(i+1:)
        i = index ( line, ':' )
        if ( i /= 0 ) then
          line(i:i) = ''
          read ( line, * ) ref
          maxAbsVal = maxval(abs(ref))
        else
          maxAbsVal = 1
        end if
      end if
    end if
    i = index ( line, 'differs by at most' )
    if ( i /= 0 ) then
      i = index ( line, 'Quantity' )
      call print_line (  trim(line(i:)) )
      sawOne = .true.
    end if
    if ( index ( line, 'Max. absolute:' ) /= 0 ) then
      i = index ( line, ':' )
      if ( i /= 0 ) then
        if ( needHead ) then
          if ( .not. summary ) then
            if ( doLog ) then
              write ( *, 11 )
           11 format ( ' Log10(Max Diff /' / &
                     & ' Max Abs Val )    Max Diff      Max Abs Val   Jacobian Block' )
            else
              write ( *, 12 )
           12 format ( ' Max Diff /' / &
                     & ' Max Abs Val      Max Diff      Max Abs Val   Jacobian Block' )
            end if
          end if
          needHead = .false.
          sawOne = .true.
        end if
        read ( line(i+1:), * ) maxAbsDiff
        if ( maxAbsVal == 0 ) maxAbsVal = 1
        relMaxAbsDiff = maxAbsDiff / maxAbsVal
        RelMaxAbsDiffPrint = relMaxAbsDiff
        if ( doLog .and. relMaxAbsDiff > 0 ) relMaxAbsDiffPrint = log10(relMaxAbsDiff) 
        write ( line, 13 ) relMaxAbsDiffPrint, maxAbsDiff, maxAbsVal, &
          & trim(block)
     13 format ( 1pg12.5, 3x, 2g14.5, t47, a )
        if ( maxAbsDiff == 0 ) line(1:12) = 'identical'
        if ( .not. summary ) call print_line ( trim(line) )
        if ( relMaxAbsDiff > maxMaxAbsDiff ) then
          maxMaxAbsDiff = relMaxAbsDiff
          whatMaxAbsDiff = block
        end if
      end if
    else
      i = index ( line, 'Radiances' )
      if ( i /= 0 ) call print_line (  trim(line(i-1:)) )
    end if

  end do

9 continue

contains

  subroutine Usage
    call get_command_argument ( 0, line )
    print 3, 'Usage: ', trim(line), ' [ options ] [ hunt ] < input > output'
  3 format ( 3a )
    print 3, ' Options:'
    print 3, '  -a => Show all details, default ', merge('summary','details',summary)
    print 3, '  -s => Show summary'
    print 3, '  -anything else => show usage'
    print 3, ' "hunt" is the text to look for to identify a new run, default "', &
          &  trim(hunt), '"'
    stop
  end subroutine Usage
  
  subroutine print_line ( line )
    character(len=*), intent(in)       :: line
    if ( LastLine /= line ) write ( *, '(a)' ) trim(line)
    if ( len_trim(line) > 0 ) LastLine = Line
  end subroutine print_line
    

end program

! $Log$
! Revision 1.1  2018/04/27 00:08:38  pwagner
! First commit
!
