!
module VELCOR_M
  use STRINGS, only: STRLWR
  use UNITS, only: velcor_unit
  use MLSCommon, only: I4, R4, R8
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------
contains
!---------------------------------------------------------------------

SUBROUTINE velcor(InDir, tmpts, tot_vel_z_lin_val, ier)

CHARACTER (LEN=*), INTENT(IN) :: InDir, Tmpts

INTEGER(i4), INTENT(OUT) :: ier
REAL(r8), INTENT(OUT) :: tot_vel_z_lin_val

REAL(r8) :: vel_z
INTEGER(i4) :: j, ia, io

CHARACTER (LEN=8) :: time_stamp,ts
CHARACTER (LEN=78) :: absfile
!
! Begin code:
!
  ier = 1
  absfile(1:)=' '
  time_stamp(1:)=' '
  j = LEN_TRIM(InDir)
  absfile = InDir(1:j)//'velcor.dat'
  ia = LEN_TRIM(absfile)
  OPEN(velcor_unit,FILE=absfile,STATUS='OLD',action='READ',IOSTAT=io)
  IF(io /= 0) THEN
    PRINT *,'** Error: Velocity correction file: ',absfile(1:ia)
    CALL errmsg(' ',io)
    RETURN
  END IF

  ts(1:)=' '
  ts = tmpts
  CALL strlwr(ts)

  10  time_stamp(1:)=' '
  READ(velcor_unit,*,END=20,IOSTAT=io) time_stamp, vel_z
  IF(io /= 0) THEN
    PRINT *,'** Read error in file: ',absfile(1:ia),CHAR(7)
    RETURN
  END IF

  CALL strlwr(time_stamp)
  IF(time_stamp(1:7) /= ts(1:7)) GO TO 10

  20  CLOSE(velcor_unit,IOSTAT=io)
  IF(time_stamp(1:7) /= ts(1:7)) THEN
    PRINT *,' ** The requested Time_Stamp not found in file ..'
    PRINT *,'     file: '//absfile(1:ia)
    PRINT *,'     Time_Stamp: ',ts
    RETURN
  END IF

!  Define the velocity correction:

  ier = 0
  tot_vel_z_lin_val = vel_z

  RETURN
END SUBROUTINE velcor
end module VELCOR_M
! $Log$
! Revision 1.1 2000/06/09 00:08:14  Z.Shippony
! Initial conversion to Fortran 90
