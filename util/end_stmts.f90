program END_STMTS

! Read a Fortran program file, in fixed format.  Copy all lines in it,
! except for END statements, to the output unchanged.

! When a SUBROUTINE or FUNCTION statement is found, remember which
! it is, and the name.  THIS PROGRAM ASSUMES THE NAME IS ON THE SAME LINE!

! When an END statement is found, add SUBROUTINE or FUNCTION as remembered,
! and the remembered name.

!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
!-----------------------------------------------------------------------------

  integer :: I1, I2                     ! Begin and end of SUBROUTINE or
                                        ! FUNCTION in LINE_CAP, else zero.
  integer :: I3, I4                     ! Begin and end of name in LINE_IN.
  integer :: IC                         ! Position of ! in LINE_CAP, else 133
  character(len=133) :: LINE_CAP        ! LINE_IN, capitalized
  character(len=133) :: LINE_IN         ! Input line
  character(len=133) :: NAME = ' '      ! Name of procedure
  character(len=10) :: TYPE = ' '       ! SUBROUTINE or FUNCTION

  do
    read ( *, '(a)', end=999 ) line_in
    line_cap = capitalize(line_in)
    ic = index(line_cap, '!')
    if ( ic == 0 ) ic = len(line_cap) + 1

    ! Look for SUBROUTINE or FUNCTION not after !
    i2 = 0
    i1 = index(line_cap, 'SUBROUTINE')
    if ( i1 /= 0 .and. i1 < ic ) then
      i2 = i1 + len('SUBROUTINE') - 1
    else
      i1 = index(line_cap, 'FUNCTION')
      if ( i1 /= 0 .and. i1 < ic ) i2 = i1 + len('FUNCTION') - 1
    end if
    if ( i2 /= 0 ) then ! Find the name after SUBROUTINE or FUNCTION
      do i3 = i2+1, 132
        if ( line_in(i3:i3) /= ' ' ) then
          type = line_in(i1:i2)
          i4 = index(line_in(i3:), ' ') + i3 - 2
          name = line_in(i3:i4)
          exit
        end if
      end do
      write ( *, '(a)' ) trim(line_in)
    else

      ! Look for END statement
      i1 = index(line_cap, 'END')
      if ( i1 /= 0 .and. i1 < ic .and. line_cap(i1+3:ic-1) == ' ' ) then
        write ( *, '(a)' ) '      ' // line_in(i1:i1+2) // ' ' // &
        &                  trim(type) // ' ' // trim(name) // ' ' // &
        &                  line_in(ic:)
      else
        write ( *, '(a)' ) trim(line_in)
      end if
    end if
  end do
999 stop

contains
  character(len=len(in)) function CAPITALIZE ( IN )
    character(len=*), intent(in) :: IN
    integer :: I                        ! Loop inductor
    do i = 1, len(in)
      if ( in(i:i) < 'a' .or. in(i:i) > 'z' ) then
        capitalize(i:i) = in(i:i)
      else
        capitalize(i:i) = char( ichar(in(i:i)) + ichar('A') - ichar('a') )
      end if
    end do
  end function CAPITALIZE
end program END_STMTS
! $Log$
