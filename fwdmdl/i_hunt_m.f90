! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!
module I_HUNT_M
  use MLSCommon, only: I4
  implicit NONE
  private
  public :: I_HUNT, HUNT
  interface HUNT; module procedure I_HUNT; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
! A binary search routine with a hunt procedure, to start from last known
! location (if 0 < JLO < N) or from the begining otherwise.
  SUBROUTINE I_HUNT (ELEMENT, ARRAY, N, JLO, JHI )
!
    integer(i4), intent(in) :: N
!
    Integer(i4), intent(in) :: ELEMENT, ARRAY(:)
!
    integer(i4), intent(in out) :: JLO, JHI

    integer :: INC,JM
    logical :: ASCND
!
    if (n < 2) then
      jlo = 1
      jhi = 1
      return
    end if
!
    ascnd = (array(n) > array(1))
!
    if (jlo <= 0 .or. jlo > n) then
      jlo = 1
      jhi = n
    else
!
      inc = 1
!
      if (element >= array(jlo) .eqv. ascnd) then
 10     jhi = jlo+inc
        if (jhi > n) then
          jhi = n
          if (jlo == n) jlo = n-1
        else if (element >= array(jhi) .eqv. ascnd) then
          jlo = jhi
          inc = inc+inc
          goto 10
        end if
      else
        jhi = jlo
 20     jlo = jhi-inc
        if (jlo < 1) then
          jlo = 1
          if (jhi == 1) jhi = 2
        else if (element < array(jlo) .eqv. ascnd) then
          jhi = jlo
          inc = inc+inc
          goto 20
        end if
      end if
    end if
!
    do while (jhi-jlo > 1)
      jm = (jhi + jlo) / 2
      if (element > array(jm) .eqv. ascnd) then
        jlo = jm
      else
        jhi = jm
      end if
    end do
!
  end subroutine I_HUNT
end module I_HUNT_M
! $Log$
! Revision 1.5  2001/06/07 23:39:31  pwagner
! Added Copyright statement
!
! Revision 1.4  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.3  2001/03/09 00:40:32  zvi
! Correcting an error in HUNT routine
!
! Revision 1.2  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:04  Z.Shippony
! Initial conversion to Fortran 90
!
