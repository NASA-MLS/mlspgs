!=================================
PROGRAM outputtest ! tests subroutine
!=================================

   USE MLSCommon , ONLY: r4, r8
   USE output_m , ONLY: blanks, output

   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the output subroutine.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it

! Variables

integer :: i
real(r4) :: t1, t2
real(r8) :: d1, d2
character(len=*), parameter :: fmt = '(g8.2)'
! character(len=*), parameter :: fmt = '(1pe9.2)'

do i=1, 99
  t1 = (i-1)*.01
  t2 = t1 + 1.
  call output(i, format='(i3)')
  call output(t1, format=fmt, advance='no')
  call output(t2, format=fmt, advance='no')
  call blanks(2)
  call output(i+1, format='(i3)', advance='yes')
enddo

do i=1, 99
  d1 = (i-1)*.01
  d2 = d1 + 1.
  call output(i, format='(i3)')
  call output(d1, format=fmt, advance='no')
  call output(d2, format=fmt, advance='no')
  call blanks(2)
  call output(i+1, format='(i3)', advance='yes')
enddo

!==================
END PROGRAM outputtest
!==================

! $Log$
