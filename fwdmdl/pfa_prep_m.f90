!
module PFA_PREP_M
  use L2PCDim, only: NLVL
  use STRINGS, only: STRLWR
  use MLSCommon, only: I4, R4, R8
  use FILTER_SW_M, only: FILTER
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, PFA_SLAB, MAXFILTPTS, &
                                 MAXAITKENPTS, MAXSPS
  use EOS_MDB, only: EOS_MDB_HDR, EOS_MDB_REC,MAX_FREQ, MAX_ZETA, &
                     MAX_TEMP, MAX_NO_LINES
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA
  use GET_BETA_PATH_M, only: GET_BETA_PATH
  use DSIMPSON_MODULE, only: DSIMPS
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
     "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
     "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------

SUBROUTINE pfa_prep(atmospheric,no_atmos,no_pfa_ch,no_filt_pts,pfa_ch,    &
           no_int_frqs,pfa_spectrum,ndx_sps,f_grid_filter,f_grid,freqs,   &
           filter_func,no_tan_hts,no_freqs,ndx_path,z_path,t_path,  &
           beta_path,freq_grid,InDir,ld,primag,ier)

!  ===============================================================
!  Declaration of variables for sub-program: pfa_prep
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_atmos, no_pfa_ch, no_filt_pts, pfa_ch(:), &
                           no_int_frqs(:), no_freqs, no_tan_hts

Integer(i4), INTENT(OUT) :: ier, ld, ndx_sps(:,:)

Real(r8), INTENT(IN) :: freqs(:),freq_grid(:)

Type(path_index), INTENT(IN) :: ndx_path(*)
Type(path_vector), INTENT(IN) :: z_path(*),t_path(*)

Type (atmos_comp), INTENT(IN) :: atmospheric(*)

Type (pfa_slab), INTENT(INOUT) :: pfa_spectrum(6,*)

Real(r8), INTENT(OUT) :: filter_func(:,:), f_grid(:,:), f_grid_filter(:,:)

Type(path_beta), INTENT(OUT) :: beta_path(:,:,:)  ! (sps_i,frq_i,ptg_i)

Character (LEN=*), INTENT(IN) :: InDir, primag

!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Real(r8), PARAMETER :: Tiny = epsilon(freqs(1))
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, k, m, ch_i, sps_ind, sps_i, pb, no_sps, spectag, &
               j4, mch, band, ptg_i, frq_i, spectags(MAXSPS)

Real(r8) :: xlhs, xrhs, df, q, area, frq

Character (LEN=8) :: sps_name

! Type (eos_mdb_hdr) :: MDB_HDR(MAXSPS)
! Type (eos_mdb_rec) :: MDB_REC(MAX_NO_LINES)
Type (eos_mdb_hdr) :: MDB_HDR(02)      ! ** DEBUG
Type (eos_mdb_rec) :: MDB_REC(02)      ! ** DEBUG

Real(r8), DIMENSION(:), ALLOCATABLE :: values, t_power, dbeta_dw, &
                                       dbeta_dn, dbeta_dnu

! begin code:

  ier = 0

!  Read in the EOS Beta database for all species involved

  pb = 0
  spectags(1:MAXSPS) = -1

  DO ch_i = 1, no_pfa_ch

    mch = pfa_ch(ch_i)
    band = (mch + 14) / 15
    no_sps = pfa_spectrum(band,1)%no_sps

sps:DO sps_i = 1, no_sps

      Spectag = pfa_spectrum(band,sps_i)%sps_spectag

!  Build: ndx_sps, the index in the mixing ratio array:

      sps_ind = 1
      DO WHILE(sps_ind < no_atmos .AND. spectag /=  &
            atmospheric(sps_ind)%spectag)
        sps_ind = sps_ind + 1
      END DO

!  Check if sps is in the database

      IF(spectag /= atmospheric(sps_ind)%spectag) THEN
        ier = 1
        WRITE(6,900) spectag
        RETURN
      END IF

!  'sps_ind' locates the index in the mixing ratio array:

      ndx_sps(sps_i,ch_i) = sps_ind

      IF(.NOT.atmospheric(sps_ind)%fwd_calc(band)) CYCLE sps

!  Now check if we loaded this Spectag already ..

      IF(pb > 0) THEN
        DO k = 1, pb
          IF(spectag == spectags(k)) CYCLE sps
        END DO
      END IF

      IF(pb == 9) THEN
        ier = 1
        PRINT *,'** Error in routine: read_EOS_db ..'
        PRINT *,'   This limited version does not support more then &
            &nine species..'
        RETURN
      END IF

      pb = pb + 1
      spectags(pb) = Spectag
      CALL read_eos_db(Spectag, sps_name, mdb_hdr(pb), mdb_rec, Ier)
      IF(ier /= 0) THEN
        PRINT *,'** Error in routine: read_EOS_db, Spectag:',spectag
        PRINT *,'   Called by routine: pfa_prep ..'
        RETURN
      END IF
!
! Now, build the beta arrays along the path
!
      do ptg_i = 1, no_tan_hts
!
        m = ndx_path(ptg_i)%total_number_of_elements
        ALLOCATE(values(m),t_power(m),dbeta_dw(m),dbeta_dn(m), &
                 dbeta_dnu(m),STAT=ier)
        if(ier /= 0) then
          PRINT *,'** Allocation error in routine: pfa_prep ..'
          PRINT *,'   IER =',ier
          RETURN
        endif
!
        do frq_i = 1, no_freqs
!
          Frq = freq_grid(frq_i)
          CALL get_beta_path (Spectag, Frq, m, z_path(ptg_i), t_path(ptg_i), &
       &                mdb_hdr(pb), mdb_rec, values, t_power, dbeta_dw, &
       &                dbeta_dn,dbeta_dnu)
!
          ALLOCATE(beta_path(sps_ind,frq_i,ptg_i)%values(m),    &
       &           beta_path(sps_ind,frq_i,ptg_i)%t_power(m),   &
       &           beta_path(sps_ind,frq_i,ptg_i)%dbeta_dw(m),  &
       &           beta_path(sps_ind,frq_i,ptg_i)%dbeta_dn(m),  &
       &           beta_path(sps_ind,frq_i,ptg_i)%dbeta_dnu(m), &
                   STAT = ier)
          if(ier /= 0) then
            PRINT *,'** Allocation error in routine: pfa_prep ..'
            PRINT *,'   IER =',ier
            RETURN
          endif
!
          beta_path(sps_ind,frq_i,ptg_i)%values(1:m) = values(1:m)
          beta_path(sps_ind,frq_i,ptg_i)%t_power(1:m) = t_power(1:m)
          beta_path(sps_ind,frq_i,ptg_i)%dbeta_dw(1:m) = dbeta_dw(1:m)
          beta_path(sps_ind,frq_i,ptg_i)%dbeta_dn(1:m) = dbeta_dn(1:m)
          beta_path(sps_ind,frq_i,ptg_i)%dbeta_dnu(1:m) = dbeta_dnu(1:m)
!
        end do

        DEALLOCATE(values,t_power,dbeta_dw,dbeta_dn,dbeta_dnu,STAT=i)

      end do

    END DO sps

  END DO

!  Find the species index in the l2pc mixing ratio database:

  m = pb
  pb = -1
  DO ch_i = 1, no_pfa_ch

    mch = pfa_ch(ch_i)
    band = (mch + 14) / 15
    no_sps = pfa_spectrum(band,1)%no_sps

    j4 = 4 * no_int_frqs(ch_i) + 1     ! Total # of Aitken's points

    IF(j4 > 1) THEN

      DO sps_i = 1, no_sps

        spectag = pfa_spectrum(band,sps_i)%sps_spectag
        spectags(sps_i) = spectag

!  Setup the spectrum record (structure), per band:

        j = 1
        k = spectag
        IF(k == 18999 .OR. k == 28964 .OR. k == 28965) j = -1

        IF(j >= 0 .AND. pb /= band) THEN

          k = -sps_i
          IF(sps_i > 1) k = sps_i

        END IF
!
!  Overwrite the pfa_spectrum data with the appropriate EOS Database
!
        j = 0
        k = -1
        DO WHILE(j < m .AND. k < 1)
          j = j + 1
          IF(mdb_hdr(j)%spectag == spectag) k = j
        END DO

        IF(k > 0) THEN
          j = mdb_hdr(k)%no_lines
          pfa_spectrum(band,sps_i)%no_lines = j
        END IF

      END DO

    END IF

    pb = band
    frq = freqs(mch)
    IF(frq < 1.0D0) THEN
      ier = 1
      WRITE(6,905) mch
      RETURN
    END IF

! Set up filter's response function

    q = 0.0
    df = Filter(q,mch,xlhs,xrhs,area,ier,InDir,primag,ld)
    IF(ier /= 0) GO TO 99

    df = (xrhs-xlhs)/(no_filt_pts-1)
    DO j = 1, no_filt_pts
      q = xlhs + (j - 1) * df
      f_grid_filter(j,ch_i) = frq + q
      filter_func(j,ch_i) = Filter(q)
    END DO

!  Normalize the filter's response array:

    CALL dsimps(filter_func(1:,ch_i),df,no_filt_pts,q)
    DO j = 1, no_filt_pts
      filter_func(j,ch_i) = filter_func(j,ch_i) / q
    END DO

!  If needed, Get Aitken's grid points:

    IF(j4 > 1) THEN
      df = (xrhs - xlhs) / (j4 - 1)
      DO j = 1, j4
        q = xlhs + (j - 1) * df
        f_grid(j,ch_i) = frq + q
      END DO
      f_grid(j4+1,ch_i) = frq                 ! *** DEBUG
    END IF

  END DO                         ! On ch_i

 99 IF(ier /= 0) THEN
      PRINT *,' ** Error in Pfa_Prep subroutine **'
      PRINT *,'    After calling subroutine: Get_Filter_Param'
      CALL errmsg(' ',ier)
    END IF

  900  FORMAT(' ** Error in Pfa_Prep subroutine **',/, &
      '    PFA species: ',i7.7,' not among L2PC species database !')
  905  FORMAT(' ** Error in Pfa_Prep subroutine **',/, &
      '    Inconsistant User Input.',/, &
      '    PFA Channel:',i3,' not among the non-PFA channels !')

  RETURN
END SUBROUTINE pfa_prep

!---------------------------------------------------------------------
!  This routine reads the EOS database and returns the structer holding
!  the required Spectag

SUBROUTINE read_eos_db(spectag,sps_name,mdb_hdr,mdb_rec,ier)

!  ===============================================================
!  Declaration of variables for sub-program: read_eos_db
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: spectag
Integer(i4), INTENT(OUT) :: ier

Character (LEN=8), intent (OUT) :: sps_name

type (eos_mdb_hdr), intent (OUT) :: mdb_hdr
type (eos_mdb_rec), intent (OUT) :: mdb_rec(*)
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, k, io, nl, du, iu, ii, jj, kk, no_lines

Integer(i4), save :: init = 0

Integer(i4) :: no_sps = 14
Integer(i4) :: spectags(14) = (/                          &
        32001, 34001, 18003, 44004, 63001, 27001, 51002,  &
        48004, 28001, 52006, 33001, 36001, 97001, 17001/)

Character (LEN=80) :: fhd, fdt, datdir
Character (LEN=8) :: names(14) = (/                                &
           'O2      ','O-18-O  ','H2O     ','N2O     ','HNO3    ', &
           'HCN     ','CLO     ','O3      ','CO      ','HOCL    ', &
           'HO2     ','HCL     ','BR-81-O ','OH      '/)

type (eos_mdb_hdr), save :: mdb_zero_hdr
type (eos_mdb_rec), save :: mdb_zero_rec

! Begin code:

  IF(init < 1) THEN
!
! Initialize mdb_zero_hdr:
!
    init = 5
    mdb_zero_hdr%spectag = 0
    mdb_zero_hdr%no_lines = 0
    mdb_zero_hdr%q_log(1) = 0.0
    mdb_zero_hdr%q_log(2) = 0.0
    mdb_zero_hdr%q_log(3) = 0.0
    DO ii = 1, max_no_lines
      mdb_zero_hdr%n(ii) = 0.0
      mdb_zero_hdr%w(ii) = 0.0
      mdb_zero_hdr%el(ii) = 0.0
      mdb_zero_hdr%n1(ii) = 0.0
      mdb_zero_hdr%n2(ii) = 0.0
      mdb_zero_hdr%v0(ii) = 0.0
      mdb_zero_hdr%ps(ii) = 0.0
      mdb_zero_hdr%delta(ii) = 0.0
      mdb_zero_hdr%gamma(ii) = 0.0
      mdb_zero_hdr%log_i(ii) = 0.0
      mdb_zero_hdr%no_f_grid(ii) = 0.0
      DO jj = 1, max_freq
        mdb_zero_hdr%x_grid(jj,ii) = 0.0
      END DO
    END DO

    DO ii = 1, max_zeta
      mdb_zero_hdr%zeta(ii) = 0.0
    END DO

    DO ii = 1, max_temp
      mdb_zero_hdr%log_temp(ii) = 0.0
    END DO
!
! Initialize mdb_zero_rec:
!
    DO ii = 1,max_zeta
      DO jj = 1, max_temp
        DO kk = 1, max_freq
          mdb_zero_rec%log_beta(ii,jj,kk) = 0.0
          mdb_zero_rec%dlog_beta_dw(ii,jj,kk) = 0.0
          mdb_zero_rec%dlog_beta_dn(ii,jj,kk) = 0.0
          mdb_zero_rec%dlog_beta_dnu0(ii,jj,kk) = 0.0
        END DO
      END DO
    END DO

  END IF

  i = 0
  j = 0
  ier = 0
  sps_name(1:) = ' '
  DO WHILE(j < no_sps .AND. i < 1)
    j = j + 1
    IF(spectags(j) == spectag) THEN
      i = j
      j = 23
    END IF
  END DO

  IF(i < 1) THEN
    io = -1
    GO TO 99
  END IF

  sps_name = names(i)
  CALL strlwr(sps_name)

  datdir(1:) = ' '
! datdir= '/home/zvi/data/'                 ! PC
! DatDir= '/zvi/eos/data/'                  ! SUN
  DatDir= '/user5/zvi/linux/MLS/data/'      ! MLSGATE

  fdt(1:) = ' '
  fhd(1:) = ' '
  j = LEN_TRIM(datdir)
  i = LEN_TRIM(sps_name)
  fdt = datdir(1:j)//sps_name(1:i)//'_eosmdb.dat'
! CALL strlwr(fdt)
  j = LEN_TRIM(fdt)
  i = INDEX(fdt,'.dat')
  fhd = fdt(1:i-1)//'.hdr'

  du = 43
  iu = du + 1
  CLOSE(iu,IOSTAT=io)
  inquire (iolength=k) mdb_hdr
  OPEN(iu,FILE=fhd,FORM='UNFORMATTED',STATUS='OLD',action='READ', &
       RECL=k,ACCESS='DIRECT',IOSTAT=io)
  IF(io /= 0) GO TO 99

  CLOSE(du,IOSTAT=io)
  inquire (iolength=j) mdb_rec(1)
  OPEN(du,FILE=fdt,FORM='UNFORMATTED',STATUS='OLD',action='READ', &
       RECL=j,ACCESS='DIRECT',IOSTAT=io)
  IF(io /= 0) GO TO 99

  mdb_hdr = mdb_zero_hdr                    ! Initialize mdb_hdr
  READ(iu,REC=1,IOSTAT=io) mdb_hdr
  IF(io /= 0) GO TO 99

  no_lines = mdb_hdr%no_lines

  DO nl = 1, no_lines

    mdb_rec(nl) = mdb_zero_rec              ! Initialize mdb_rec
    READ(du,REC=nl,IOSTAT=io) mdb_rec(nl)
    IF(io /= 0) GO TO 99

  END DO

  99   IF(io > 0) THEN
    ier = 1
    CALL errmsg(' ',io)
  END IF

  CLOSE(du,IOSTAT=i)
  CLOSE(iu,IOSTAT=i)

  RETURN
END SUBROUTINE read_eos_db
end module PFA_PREP_M
! $Log$
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
