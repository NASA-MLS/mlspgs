! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program UnwrapList

  ! Given a comma-separated list, perhaps with Fortran & continuation
  ! at the beginnings and ends of lines, and Fortran ! comments, unwrap
  ! it into one-per-line output.

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  !---------------------------------------------------------------------------

  integer :: I
  character(255) :: Line

  do
    read ( *, '(a)', end=9 ) line
    do
      line = adjustl(line)
      if ( line(1:1) == '&' ) line = adjustl(line(2:))
      if ( line(1:1) == ',' ) line = adjustl(line(2:))
      i = scan(line,' ,&!')
      if ( i < 2 ) exit
      write ( *, '(a)' ) line(1:i-1)
      line(1:i-1) = ''
    end do
  end do
9 continue

end program UnwrapList

! $Log$
