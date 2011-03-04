! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program Goldbrick_More

  implicit NONE

!---------------------------- RCS Info -------------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

  character(len=128) :: Line
  real :: RefMin, RefMax
  real :: MaxAbs
  integer :: I

  do
    read ( *, '(a)', end=9 ) line
    write ( *, '(a)' ) trim(line)
    if ( index(line,'Reference') /= 0 ) then
      i = index(line,'):')
      line = line(i+2:)
      i = index(line,':')
      line(i:i) = ','
      read ( line, * ) RefMin, RefMax
    end if
    i = index(line,'absolute:')
    if ( i /= 0 ) then
      read ( line(i+9:), * ) maxAbs
      write ( *, '(" ** (Max absolute) / max: ", 1p,3g11.5)' ) &
        & maxabs/max(abs(refmin),abs(refmax))
    end if
  end do

9 continue

end program Goldbrick_More

! $Log$
