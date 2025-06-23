! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program TeX_Grammar

  implicit NONE

  ! Convert the productions part of the output of LR to LaTeX

  integer :: J, I
  character(255) :: Line
  integer :: State = 0   ! 0 = before 'T H E   P R O D U C T I O N S'
                         ! 1 = first production
                         ! 2 = after the first production

  ! This subterfuge usually creates a character with value \ even
  ! on systems that treat \ as some kind of an escape character.
  character(*), parameter :: BS1 = '\\'
  character, parameter :: BS = BS1(1:1)

  call output ( '% LR output converted to TeX by TeX_Grammar program', advance='yes' )
  do
    read ( *, '(a)', end=9 ) line
    select case ( state )
    case ( 0 ) ! before 'T H E   P R O D U C T I O N S'
      if ( index(line,'T H E   P R O D U C T I O N S') /= 0 ) state = 1
    case ( 1, 2 )
      if ( line == '' ) cycle
      if ( index(line,'V O C A B U L A R Y') /= 0 ) exit
      line = adjustl(line)
      i = index(line,' ') ! After the production number
      line = adjustl(line(i:)) ! Skip the production number
      i = index(line,'->')
      if ( i == 0 ) cycle ! Not a production
      if ( i /= 1 ) then  ! First production for a symbol; output the LHS
        if ( state == 2 ) call TeX ( 'end{tabbing}', indent=4 )
        call TeX ( 'item', indent=2 )
        call TeX ( 'begin{tabbing}', indent=4 )
        call convert ( trim(line(:i-1)), advance='no', indent=6 )
        call TeX ( '=', advance='no', indent=1 )
      else
        call TeX ( '>', advance='no', indent=6 )
      end if
      i = i + 3 ! Skip ->
      call output ( ' $', advance='no' )
      call TeX ( 'rightarrow$', advance='no' )
      do
        line = adjustl(line(i:))
        if ( line == '' ) exit
        i = index(line,' ' )
        if ( line(1:3) == '=> ' ) then
          call output ( ' $', advance='no' )
          call TeX ( 'Rightarrow$', advance='no' )
        else
          call convert ( line(:i-1), advance='no', indent=1 )
        end if
      end do
      call TeX ( bs )
      state = 2
    end select
  end do

9 continue
  if ( state == 2 ) call TeX ( 'end{tabbing}', indent=4 )

contains

  subroutine Convert ( Text, Advance, Indent )
    ! Replace < and > with $<$ and $>$, then output Text
    character(len=*), intent(in) :: Text
    character(len=*), intent(in), optional :: Advance ! Default 'yes'
    integer, intent(in), optional :: Indent ! Default zero
    integer :: I
    integer :: MyIndent
    myIndent = 0
    if ( present(indent) ) myIndent = indent
    call output ( repeat(' ',myIndent), advance='no' )
    do i = 1, len(text)
      select case ( line(i:i) )
      case ( '<', '>' )
        call output ( '$' //line(i:i) // '$', advance='no' )
      case ( '^' )
        call output ( '$', advance='no' )
        call TeX ( 'wedge', advance='no' )
        call output ( '$', advance='no' )
      case ( '#' )
        call TeX ( '#', advance='no' )
      case ( bs )
        call output ( '$', advance='no' )
        call TeX ( 'backslash', advance='no' )
        call output ( '$', advance='no' )
      case ( '{', '}' )
        call TeX ( line(i:i), advance='no' )
      case default
        call output ( line(i:i), advance='no' )
      end select
    end do
    call output ( '', advance=advance )
  end subroutine Convert

  subroutine Output ( Text, Advance )
    ! Output, with trim if advancing
    character(len=*), intent(in) :: Text
    character(len=*), intent(in), optional :: Advance ! Default 'yes'
    character(len=3) :: MyAdv
    myAdv = 'yes'
    if ( present(advance) ) myAdv = advance
    if ( myAdv == 'yes' ) then
      write ( *, '(a)' ) trim(text)
    else
      write ( *, '(a)', advance='no' ) text
    end if
  end subroutine Output

  subroutine TeX ( Text, Advance, Indent )
    ! Output a TeX command with \ before it
    character(len=*), intent(in) :: Text
    character(len=*), intent(in), optional :: Advance ! Default 'yes'
    integer, intent(in), optional :: Indent ! Default zero
    integer :: MyIndent
    myIndent = 0
    if ( present(indent) ) myIndent = indent
    call output ( repeat(' ',myIndent) // bs // trim(text), advance=advance )
  end subroutine TeX 
  
end program TeX_Grammar

! $Id$

! $Log$
! Revision 1.1  2014/02/28 22:01:25  vsnyder
! Initial commit
!
