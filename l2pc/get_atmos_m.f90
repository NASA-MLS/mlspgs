!
module GET_ATMOS_M
  use UNITS, only: model_unit
  use MLSCommon, only: I4, R4, R8
  use L2PCDim, only: NLVL
  use STRINGS, only: LEFTJ, SQZSTR, STRUPR, STRLWR
  Implicit None
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------

SUBROUTINE get_atmos(spectag,time_stamp,mr_f_sps,sps_press,no_pts, &
                     no_phi,InDir,ld,ier)

!  ===============================================================
!  Declaration of variables for sub-program: get_atmos
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: spectag, ld

Integer(i4), INTENT(OUT) :: ier, no_pts, no_phi

Real(r8), INTENT(OUT) :: sps_press(:), mr_f_sps(:,:)

Character (LEN=*), INTENT(IN) :: time_stamp, InDir
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: io
Character (LEN=80) :: fn

! Begin code:

  ier = 0
  fn(1:) = ' '

  write(Fn,900) InDir(1:ld),time_stamp(1:5),spectag
900  format(A,'zpmodel_',A,'_',i7.7,'.data')

! ** DEBUG, use a "cooked up" file

     fn(1:) = ' '
     WRITE(fn,901) time_stamp(1:5),spectag
901  FORMAT('zpmodel_',a,'_',i7.7,'.data')
     fn = '/home/zvi/'//fn

! ** END DEBUG

  CALL strlwr(fn)
  CALL read_zpm(fn,model_unit,no_pts,no_phi,sps_press,mr_f_sps,io)
  IF(io /= 0) THEN
    ier = 1
    PRINT *,'** Error in: get_atmos subroutine ..'
    CALL errmsg(fn,io)
  END IF

  RETURN
END SUBROUTINE get_atmos

!---------------------------------------------------------------

SUBROUTINE read_zpm(fn,uu,no_zeta,no_phi,zeta,valpz,ier)

!  ===============================================================
!  Declaration of variables for sub-program: read_zpm
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: uu
Integer(i4), INTENT(OUT) :: no_zeta, no_phi, ier

Real(r8), INTENT(OUT) :: zeta(:), valpz(:,:)

Character (LEN=*), INTENT(IN) :: fn
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, k, iz, phi_index

Real(r8) :: vec(nlvl)

Character (LEN=80) :: line

! begin code

  ier = 0
  no_phi = -1
  no_zeta = 0
  j = LEN_TRIM(fn)
  CLOSE(uu,IOSTAT=iz)
  OPEN(uu,FILE=fn,STATUS='OLD',action='READ',IOSTAT=iz)
  IF(iz /= 0) THEN
    PRINT *,'* Error: Could not open File:'
    PRINT *,'   ',fn(1:j)
    GO TO 99
  END IF

  READ(uu,'(A)',IOSTAT=iz) line
  IF(iz /= 0) THEN
    PRINT *,'* Error: Could not read File:'
    PRINT *,'   ',fn(1:j)
    GO TO 99
  END IF

  iz = 0
  DO WHILE(iz == 0)

    line(1:) = ' '
    READ(uu,'(A)',IOSTAT=iz) line
    IF(iz /= 0) GO TO 99

    CALL leftj(line)
    CALL strupr(line)

    IF(line(1:17) == 'NO_OF_ZETA_LEVELS') THEN
      i = INDEX(line,'=')
      line(1:i) = ' '
      CALL leftj(line)
      READ(line,*,IOSTAT=iz) j
      IF(iz /= 0) GO TO 99
      IF(j > 0) no_zeta = j
    ELSE IF(line(1:16) == 'NO_OF_PHI_LEVELS') THEN
      i = INDEX(line,'=')
      line(1:i) = ' '
      CALL leftj(line)
      READ(line,*,IOSTAT=iz) j
      IF(iz /= 0) GO TO 99
      IF(j > 0) no_phi = j
    ELSE IF(line(1:11) == 'ZETA_LEVELS') THEN
      IF(no_zeta < 1) THEN
        iz = -6
        line(1:) = ' '
        line = 'NO_OF_ZETA_LEVELS  is missing or misplaced ..'
        GO TO 99
      END IF
      READ(uu,*,IOSTAT=iz) (zeta(j),j=1,no_zeta)
      IF(iz /= 0) GO TO 99
    ELSE IF(line(1:19) == 'VALUES_AT_PHI_INDEX') THEN
      IF(no_zeta*no_phi < 1) THEN
        IF(no_zeta < 1) THEN
          iz = -6
          line(1:) = ' '
          line = 'NO_OF_ZETA_LEVELS  is missing or misplaced ..'
          GO TO 99
        ELSE IF(no_phi < 1) THEN
          iz = -6
          line(1:) = ' '
          line = 'NO_OF_PHI_LEVELS  is missing or misplaced ..'
          GO TO 99
        END IF
      END IF
      i = INDEX(line,'=')
      line(1:i) = ' '
      CALL leftj(line)
      READ(line,*,IOSTAT=iz) j
      IF(iz /= 0) GO TO 99
      k = (no_phi + 1) / 2
      phi_index = j + k
      READ(uu,*,IOSTAT=iz) (vec(i),i=1,no_zeta)
      IF(iz /= 0) GO TO 99
      DO i = 1, no_zeta
        valpz(i,phi_index) = vec(i)
      END DO
    END IF
  END DO

99  CLOSE(uu,IOSTAT=j)

  IF(iz == 0.OR.iz == -1) RETURN

  ier = 1
  PRINT *
  IF(iz < -2) THEN
    i = LEN_TRIM(line)
    PRINT *,'** Error in subroutine Read_Zpm'
    PRINT *,'   ',line(1:i)
  ELSE IF(iz > 0) THEN
    CALL errmsg('From subroutine Read_Zpm',iz)
  END IF

  RETURN
END SUBROUTINE read_zpm

end module GET_ATMOS_M
! $Log$
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
