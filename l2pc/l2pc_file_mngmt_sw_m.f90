!
module L2PC_FILE_MNGMT_SW_M
  use STRINGS, only: LEFTJ, SQZSTR, STRUPR, STRLWR
  use L2PCdim, only: NCH
  use TIME_MOD, only: NOW_CCSDS
  use MLSCommon, only: I4, R4, R8
  use L2PC_FILE_PARAMETERS, only: max_rec_len_key_file, l2pc_header_key1, &
                                  l2pc_header_key2, l2pc_header_key3
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE, L2PC_HEADER_TWO, &
                                  L2PC_HEADER_TRI, l2pc_keys
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains

!---------------------------------------------------------------

SUBROUTINE open_l2pc(filename, l2pc_lu, l2pc_lu_key, cstat, &
                     l2pc_rec_length, ier)

!  ===============================================================
!  Declaration of variables for sub-program: open_l2pc
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(OUT) :: ier
Integer(i4), INTENT(OUT) :: l2pc_lu
Integer(i4), INTENT(OUT) :: l2pc_lu_key

Integer(i4), INTENT(IN)  :: l2pc_rec_length

Character (LEN=*), INTENT(IN) :: filename, cstat
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, lr

Character (LEN=80) :: filekey
Character (LEN=8) :: astat

type (l2pc_keys) :: key_rec
!
!  Begin code:
!
  ier = 0
  l2pc_lu = 10
  l2pc_lu_key = 11
  CLOSE(l2pc_lu,IOSTAT=i)
  CLOSE(l2pc_lu_key,IOSTAT=i)

  astat(1:8)=' '
  astat = cstat
  CALL strupr(astat)
  IF(astat(1:3) /= 'NEW'.AND.astat(1:3) /= 'OLD') THEN
    PRINT *,' ** Error in open_l2pc subroutine..'
    PRINT *,'    Invalid status parameter: ',astat
   ier = 1
   RETURN
  END IF

  IF(astat(1:3) == 'NEW') astat='UNKNOWN'

  OPEN(l2pc_lu,FILE=filename,STATUS=astat,ACCESS='DIRECT', &
       FORM='UNFORMATTED',RECL=l2pc_rec_length,IOSTAT=ier)
  IF(ier /= 0) GO TO 99

! Do a reverse scanned search so as to avoid stepping on a directory
! specifier

  i = LEN_TRIM(filename)
  DO WHILE (filename(i:i) /= '.')
    i = i - 1
  END DO

  filekey(1:)=' '
  filekey=filename(1:i)//'key'
  OPEN(l2pc_lu_key,FILE=filekey,STATUS='UNKNOWN',iostat=lr)
  Close(l2pc_lu_key,STATUS='DELETE',iostat=lr)
  lr = 4 * max_rec_len_key_file
  OPEN(l2pc_lu_key,FILE=filekey,STATUS='NEW',ACCESS='DIRECT', &
       FORM='UNFORMATTED',RECL=lr,IOSTAT=ier)
  IF(ier /= 0) GO TO 99

  !  write the first record in the key file, close it and open it again,
  !  to insure file exists in case of a crash ..

  key_rec.rec_no = l2pc_rec_length
  key_rec.l2pc_key = 'Record length of the l2pc file .'
  WRITE(l2pc_lu_key,REC=1) key_rec

  CLOSE(l2pc_lu_key,IOSTAT=i)
  OPEN(l2pc_lu_key,FILE=filekey,STATUS='OLD',ACCESS='DIRECT', &
       FORM='UNFORMATTED',RECL=lr,IOSTAT=i)

  RETURN

99 Print 900,filekey
900  FORMAT(' ** Error in open_l2pc subroutine..',/, &
            '    Error opening file: ',a)
  CALL errmsg(' ',ier)

  RETURN
END SUBROUTINE open_l2pc

!---------------------------------------------------------------------------

SUBROUTINE position_l2pc(filename, l2pc_lu, l2pc_lu_key, irec, &
                         l2pc_rec_length, ier)

!  ===============================================================
!  Declaration of variables for sub-program: position_l2pc
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Character (LEN=*), INTENT(IN) :: filename

Integer(i4), INTENT(IN)  :: irec

Integer(i4), INTENT(OUT) :: l2pc_lu, l2pc_lu_key, l2pc_rec_length, ier
!  ----------------
!  Local variables:
!  ----------------

type (l2pc_keys) :: key_rec

Integer(i4) :: i, lr, ieropen

Character (LEN=80) :: filekey, tmp
Character (LEN=40) :: akey

  ier = 0
  ieropen = 0
  l2pc_lu = 10
  l2pc_lu_key = l2pc_lu + 1

  CLOSE(l2pc_lu,IOSTAT=i)
  CLOSE(l2pc_lu_key,IOSTAT=i)

!  Get the record length of the l2pc file from the key file that goes
!  with it. Find the file extension:  Do a reverse scanned search so as
!  to avoid stepping on a directory specifier

  filekey(1:)=' '
  i = LEN_TRIM(filename)
  DO WHILE (filename(i:i) /= '.')
    i = i - 1
  END DO

  tmp(1:)=' '
  filekey = filename(1:i)//'key'
  tmp = filekey

  CLOSE(l2pc_lu_key,IOSTAT=i)
  lr = 4 * max_rec_len_key_file
  OPEN(l2pc_lu_key,FILE=filekey,STATUS='UNKNOWN',ACCESS='DIRECT', &
       FORM='UNFORMATTED',RECL=lr,IOSTAT=ieropen)
  IF(ieropen /= 0) GO TO 99

  READ(l2pc_lu_key,REC=1,IOSTAT=ieropen) key_rec
  IF(ieropen /= 0) GO TO 99

  l2pc_rec_length = key_rec.rec_no      ! L2PC record length

  tmp = filename
  CLOSE(l2pc_lu,IOSTAT=i)
  OPEN(l2pc_lu,FILE=filename,STATUS='OLD',ACCESS='DIRECT', &
       FORM='UNFORMATTED',RECL=l2pc_rec_length,IOSTAT=ieropen)
  IF(ieropen /= 0) GO TO 99

!  Position file pointer to the last record (Irec)

  READ(l2pc_lu,REC=irec,IOSTAT=ier) akey
  IF(ier /= 0) GO TO 99

  RETURN

99   i = LEN_TRIM(tmp)
     IF(ieropen /= 0) THEN
       ier = 1
       Print 900,tmp(1:i)
       CALL errmsg(' ',ieropen)
     ELSE IF(ier /= 0) THEN
       Print 910,irec,tmp(1:i)
       CALL errmsg(' ',ier)
     END IF

900  FORMAT(' ** Error in position_l2pc subroutine..',/, &
    '    Error opening file: ',a)

910  FORMAT(' ** Error in position_l2pc subroutine..',/, &
    '    Error reading record: ',i5,' in file: ',a)

  RETURN
END SUBROUTINE position_l2pc

!---------------------------------------------------------------

SUBROUTINE re_write_header1(no_records_per_bin, l2pc_lu, ier)

!  ===============================================================
!  Declaration of variables for sub-program: re_write_header1
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN)  :: l2pc_lu, no_records_per_bin
Integer(i4), INTENT(OUT) :: ier
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, no_bands, no_pointings, no_avail_keys, &
      &        no_sv_components,no_channels_per_band

Character (LEN=40) :: akey

type (l2pc_header_one) :: header1

  READ(l2pc_lu,REC=1,IOSTAT=ier) akey,header1.line1,header1.line2, &
      header1.line3,no_bands,no_channels_per_band,no_pointings, &
      no_avail_keys,no_sv_components,header1.no_coeff_per_component, &
      header1.no_k_records_per_bin,header1.no_mag_fields, &
      header1.no_b_theta_lin_val,header1.b_phi_lin_val, &
      (header1.avail_keys(i),i=1,no_avail_keys), &
      (header1.pointings(i),i=1,no_pointings), &
      (header1.sv_components(i),i=1,no_sv_components), &
      ((header1.sv_rtrvl_by_band(i,j),j=1,no_bands), i=1,no_sv_components), &
      (header1.no_elmnts_per_sv_component(i),i=1,no_sv_components), &
      (header1.sv_component_first_elmnt_index(i),i=1,no_sv_components), &
      (header1.b_fields(i),i=1,max_no_mag_fields), &
      (header1.b_theta_lin_val(i),i=1,max_no_theta_val)
  IF(ier /= 0) GO TO 99

! Compute number of K_STAR records per bin:
!   Total number of records minus X_STAR(1) and I_STAR(6)

  header1.no_k_records_per_bin = no_records_per_bin - 7

  akey = l2pc_header_key1
  WRITE(l2pc_lu,REC=1,IOSTAT=ier) akey,header1.line1,header1.line2, &
      header1.line3,no_bands,no_channels_per_band,no_pointings, &
      no_avail_keys,no_sv_components,header1.no_coeff_per_component, &
      header1.no_k_records_per_bin,header1.no_mag_fields, &
      header1.no_b_theta_lin_val,header1.b_phi_lin_val, &
      (header1.avail_keys(i),i=1,no_avail_keys), &
      (header1.pointings(i),i=1,no_pointings), &
      (header1.sv_components(i),i=1,no_sv_components), &
      ((header1.sv_rtrvl_by_band(i,j),j=1,no_bands), i=1,no_sv_components), &
      (header1.no_elmnts_per_sv_component(i),i=1,no_sv_components), &
      (header1.sv_component_first_elmnt_index(i),i=1,no_sv_components), &
      (header1.b_fields(i),i=1,max_no_mag_fields), &
      (header1.b_theta_lin_val(i),i=1,max_no_theta_val)
  IF(ier /= 0) GO TO 99

  RETURN

99 PRINT *,' I had a problem accessing the header, sorry ..',CHAR(7)
   CALL errmsg(' ',ier)

  RETURN
END SUBROUTINE re_write_header1
!
!---------------------------------------------------------------

SUBROUTINE write_hdr(header1,header2,header3,l2pc_lu,ier)

!  ===============================================================
!  Declaration of variables for sub-program: write_hdr
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: l2pc_lu
Integer(i4), INTENT(OUT) :: Ier

type (l2pc_header_one), INTENT(IN) :: header1
type (l2pc_header_two), INTENT(IN) :: header2
type (l2pc_header_tri), INTENT(IN) :: header3

!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, no_bands, no_pointings, no_avail_keys, &
  &   no_sv_components,no_channels_per_band,no_sv_elmnts,matdim

Character (LEN=40) :: akey

! Structure defintions are here.

no_bands = header1.no_bands
no_pointings = header1.no_pointings
no_avail_keys = header1.no_avail_keys
no_sv_components = header1.no_sv_components
no_channels_per_band = header1.no_channels_per_band

  akey = l2pc_header_key1
  WRITE(l2pc_lu,REC=1,IOSTAT=ier) akey,header1.line1,header1.line2, &
      header1.line3,no_bands,no_channels_per_band,no_pointings, &
      no_avail_keys,no_sv_components,header1.no_coeff_per_component, &
      header1.no_k_records_per_bin,header1.no_mag_fields, &
      header1.no_b_theta_lin_val,header1.b_phi_lin_val, &
      (header1.avail_keys(i),i=1,no_avail_keys), &
      (header1.pointings(i),i=1,no_pointings), &
      (header1.sv_components(i),i=1,no_sv_components), &
      ((header1.sv_rtrvl_by_band(i,j),j=1,no_bands), i=1,no_sv_components), &
      (header1.no_elmnts_per_sv_component(i),i=1,no_sv_components), &
      (header1.sv_component_first_elmnt_index(i),i=1,no_sv_components), &
      (header1.b_fields(i),i=1,max_no_mag_fields), &
      (header1.b_theta_lin_val(i),i=1,max_no_theta_val)
  IF(ier /= 0) GO TO 99

  akey = l2pc_header_key2
  no_sv_elmnts = header2.no_sv_elmnts
  WRITE(l2pc_lu,REC=2,IOSTAT=ier) akey,no_sv_elmnts, &
      (header2.tri_basis_vert_grid(i),i=1,no_sv_elmnts)
  IF(ier /= 0) GO TO 99

  akey = l2pc_header_key3
  matdim = header3.matdim
  WRITE(l2pc_lu,REC=3,IOSTAT=ier) akey,matdim, &
      (header3.second_der_matrix_bands(i),i=1,matdim), &
      (header3.second_der_matrix_namid(i),i=1,matdim)
  IF(ier /= 0) GO TO 99

  RETURN

99 PRINT *,' I had a problem writing the header,sorry ..',CHAR(7)
   CALL errmsg(' ',ier)

  RETURN
END SUBROUTINE write_hdr
!
!---------------------------------------------------------------

SUBROUTINE close_l2pc_xx(l2pc_lu,l2pc_lu_key,keys,rec_nos,nrec)

!  ===============================================================
!  Declaration of variables for sub-program: close_l2pc_xx
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: rec_nos(:)
Integer(i4), INTENT(OUT) :: l2pc_lu
Integer(i4), INTENT(IN OUT) :: l2pc_lu_key
Integer(i4), INTENT(IN) :: nrec

Character (LEN=*), INTENT(IN) :: keys(:)
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: j

type (l2pc_keys) :: key_rec

  CLOSE(l2pc_lu,IOSTAT=j)          ! Close the l2pc file itself

  CALL c_qsort2(nrec,keys,rec_nos) ! Sort keys, maintain rec_nos

!  Start writing the key file, from record #2. (First record has been
!  written already, and is containing the record length of the l2pc file.)

  DO j = 1, nrec
    key_rec.l2pc_key = keys(j)
    key_rec.rec_no = rec_nos(j)
    WRITE(l2pc_lu_key,REC=j+1) key_rec
  END DO

  CLOSE(l2pc_lu_key,IOSTAT=j)      ! Close the l2pc key file

  RETURN
END SUBROUTINE close_l2pc_xx
!
!---------------------------------------------------------------

SUBROUTINE compose_header_lines(filename,no_pfa_ch,pfa_ch,lines,ier)

!  ===============================================================
!  Declaration of variables for sub-program: compose_header_lines
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_pfa_ch
Integer(i4), INTENT(IN) :: pfa_ch(:)
Integer(i4), INTENT(IN OUT) :: ier

Character (LEN=*), INTENT(IN) :: filename
Character (LEN=*), INTENT(OUT) :: lines(:)
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: uars_month, cver, sver, i, j, k, jl, jj, idash, &
               nspfa, spfach(nch), tmpch(nch)

Character (LEN=4) :: bank_status, ax
Character (LEN=12) :: fbc
Character (LEN=16) :: adir
Character (LEN=20) :: cd
Character (LEN=252) :: tmp

!  Begin code:

  lines(1)(1:)=' '
  lines(2)(1:)=' '
  lines(3)(1:)=' '

  CALL analize_filename(filename,uars_month,cver,sver,jj,bank_status,ier)
  IF(ier /= 0) THEN
    ier = 0
    k = LEN_TRIM(filename)
    PRINT *,'** Warning: Non-Standard L2PC filename !!'
    PRINT *,'   Name: ',filename(1:k)
    RETURN
  END IF

!  Get file creation date (today's date & time)

  cd(1:20)=' '
  Call NOW_CCSDS (cd)

!  Define filter's configuration

  fbc(1:12)=' '
  IF(bank_status(1:1) == 'n') THEN
    fbc = 'NORMAL'
  ELSE IF(bank_status(1:1) == 'a') THEN
    fbc = 'AB-NORMAL'
  ELSE
    fbc = 'UNKNOWN'
  END IF

  adir(1:16)=' '
  IF(jj > 0) THEN
    adir = 'Forward/South'
  ELSE
    adir = 'Backward/North'
  END IF

!  Get scalar & pfa channels

  nspfa = no_pfa_ch
  DO i = 1, no_pfa_ch
    k = pfa_ch(i)
    spfach(i) = k
  END DO

  tmp(1:)=' '
  IF(nspfa < 1) THEN
    tmp='No Scalar pfa channels,'
  ELSE
    ax(1:)=' '
    WRITE(ax,'(i2,'','')') spfach(1)
    CALL leftj(ax)
    j = LEN_TRIM(ax)
    tmp='Scalar PFA channels: '//ax(1:j)
    k = 1
    tmpch(k) = spfach(1)
    IF(nspfa > 1) THEN
      i = 1
      idash = 0
      jl = LEN_TRIM(tmp)
      DO WHILE(i < nspfa)
        i = i + 1
        IF(spfach(i)-spfach(i-1) < 2) THEN
          IF(idash == 0) THEN
            idash = 1
            IF(tmp(jl:jl) == ',') jl = jl - 1
            jl = jl + 1
            tmp(jl:jl)='-'
          END IF
          IF(i == nspfa) THEN
            ax(1:)=' '
            WRITE(ax,'(i2,'','')') spfach(i)
            CALL leftj(ax)
            j = LEN_TRIM(ax)
            tmp(jl+1:jl+j)=ax(1:j)
            jl = j + jl
          END IF
        ELSE
          idash = 0
          DO jj = 1, k
            IF(spfach(i-1) == tmpch(jj)) GO TO 100
          END DO
          ax(1:)=' '
          WRITE(ax,'(i2,'','')') spfach(i-1)
          CALL leftj(ax)
          j = LEN_TRIM(ax)
          tmp(jl+1:jl+j)=ax(1:j)
          jl = j + jl
          k = k + 1
          tmpch(k) = spfach(i-1)
  100     DO jj = 1, k
            IF(spfach(i) == tmpch(jj)) GO TO 110
          END DO
          ax(1:)=' '
          WRITE(ax,'(i2,'','')') spfach(i)
          CALL leftj(ax)
          j = LEN_TRIM(ax)
          tmp(jl+1:jl+j)=ax(1:j)
          jl = j + jl
          k = k + 1
          tmpch(k) = spfach(i)
  110     CONTINUE
        END IF
      END DO
    END IF
  END IF

  jl = LEN_TRIM(tmp)
  IF(tmp(jl:jl) == ',') THEN
    tmp(jl:)=' '
    jl = jl - 1
  END IF

  tmp=tmp(1:jl)//',No Magnetic pfa channels'

  jl = LEN_TRIM(tmp)
  IF(tmp(jl:jl) == ',') THEN
    tmp(jl:)=' '
    jl = jl - 1
  END IF

  IF(jl > 80) THEN
    CALL replace_text(tmp,'channels','ch.')
    jl = LEN_TRIM(tmp)
    IF(jl > 80) THEN
      CALL replace_text(tmp,'Magnetic','Mag.')
      jl = LEN_TRIM(tmp)
      IF(jl > 80) THEN
        CALL replace_text(tmp,'Scalar','Sca.')
        jl = LEN_TRIM(tmp)
        IF(jl > 80) THEN
          CALL replace_text(tmp,'ch.',' ')
          j = LEN_TRIM(tmp)
          jl = MIN0(j,80)
        END IF
      END IF
    END IF
  END IF

  lines(1)(1:jl) = tmp(1:jl)

  tmp(1:)=' '
  WRITE(tmp,920) cver,sver,cd
  jl = LEN_TRIM(tmp)
  lines(2)(1:jl) = tmp(1:jl)

  tmp(1:)=' '
  j = LEN_TRIM(adir)
  WRITE(tmp,930) uars_month,adir(1:j),fbc
  jl = LEN_TRIM(tmp)
  lines(3)(1:jl) = tmp(1:jl)

920  FORMAT(' Cal. Version: ',i3.3,'  S/W version: ',i3.3, &
    ' File creation date: ',a)
930  FORMAT(' UARS month: ',i2.2,'  Look Dir: ',a, &
    '  Filter bank configuration: ',a)

  RETURN
END SUBROUTINE compose_header_lines

!---------------------------------------------------------------------

SUBROUTINE analize_filename(fn,uars_month,cver,sver,idir,bank_status,ier)

!  ===============================================================
!  Declaration of variables for sub-program: analize_filename
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Character (LEN=*), INTENT(IN) :: fn

Integer(i4), INTENT(OUT) :: uars_month, idir, cver, sver
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: jj, kk

Character (LEN=1) :: ch
Character (LEN=4) :: bank_status, aver

!  file name looks like: umls_l2_cal_v032_l2pc_07f_n_003.dat
!                Legend: .............xxx......mmd.b.ccc....

!  xxx - S/W  ver      softwear version number
!  mm  - uars_month    (01  to: 10)
!  d   - Flight Direction ('f'=Forward,'b'=Backward)
!  b   - bank_status   ('n'=Normal,'a'=Abnormal)
!  ccc - cal. ver      calibration I.D. version

  ch = '@'
  ier = 1
  bank_status(1:)=' '
  l = 1 + LEN_TRIM(fn)
  DO WHILE(l > 1.AND.ch /= '.')
    l = l - 1
    ch = fn(l:l)
  END DO
  aver(1:) = ' '
  aver(1:3) = fn(l-3:l-1)
  READ(aver,*,IOSTAT=io) cver
  IF(io /= 0) RETURN
  bank_status(1:1) = fn(l-5:l-5)
  ch = fn(l-7:l-7)
  IF(ch == 'f') THEN
    idir = 1
  ELSE IF(ch == 'b') THEN
    idir = -1
  ELSE
    RETURN
  END IF
  aver(1:) = ' '
  aver(1:2) = fn(l-9:l-8)
  READ(aver,*,IOSTAT=io) uars_month
  IF(io /= 0) RETURN
  kk = 18
  jj = -1
  DO WHILE(jj < 0.AND.kk > 16)
    ch = fn(l-kk:l-kk)
    jj = ICHAR(ch)
    IF(jj < 48.OR.jj > 57) THEN
      kk = kk - 1
      jj = -1
    END IF
  END DO
  aver(1:) = ' '
  jj = kk - 15
  aver(1:jj) = fn(l-kk:l-16)
  READ(aver,*,IOSTAT=io) sver
  IF(io /= 0) RETURN
  ier = 0

  RETURN
END SUBROUTINE analize_filename

!---------------------------------------------------------------------

SUBROUTINE replace_text(line,str1,str2)

!  ===============================================================
!  Declaration of variables for sub-program: replace_text
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Character (LEN=*), INTENT(IN OUT) :: line
Character (LEN=*), INTENT(IN) :: str1, str2
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, lk, ll, l1, l2, jt1, jt2
Character (LEN=252) :: tmp1, tmp2

  i = 1
  lk = LEN(line)
  ll = LEN_TRIM(line)
  l1 = LEN_TRIM(str1)
  l2 = LEN_TRIM(str2)
  DO WHILE(i > 0)
    tmp1(1:)=' '
    tmp2(1:)=' '
    i = INDEX(line,str1)
    IF(i > 0) THEN
      IF(i > 1) tmp1(1:i-1)=line(1:i-1)
      IF(l2 > 0) tmp1(i:i+l2-1)=str2(1:l2)
      jt1 = LEN_TRIM(tmp1)
      jt2 = ll-l1-i+1
      tmp2(1:jt2)=line(i+l1:ll)
      line(1:lk)=' '
      IF(jt1 > 0) THEN
        line=tmp1(1:jt1)//tmp2(1:jt2)
      ELSE
        line=tmp2(1:jt2)
      END IF
      ll = LEN_TRIM(line)
    END IF
  END DO

  RETURN
END SUBROUTINE replace_text

end module L2PC_FILE_MNGMT_SW_M
! $Log$
! Revision 1.1 2000/06/09 00:08:14  Z.Shippony
! Initial conversion to Fortran 90
