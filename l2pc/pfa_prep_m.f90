!
module PFA_PREP_M
  use L2PCDim, only: NLVL
  use STRINGS, only: STRLWR
  use MLSCommon, only: I4, R4, R8
  use FILTER_SW_M, only: FILTER
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, PFA_SLAB, MAXFILTPTS, &
                                 MAXAITKENPTS, MAXGEOPHYS, maxrat
  use EOS_MDB, only: EOS_MDB_HDR, EOS_MDB_REC, MAX_NO_LINES, MAX_FREQ, &
                     MAX_ZETA, MAX_TEMP
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

SUBROUTINE pfa_prep(atmospheric,no_atmos,no_pfa_ch,n_lvls,no_filt_pts, &
           pfa_ch,no_int_frqs,pfa_spectrum,z_grid,t_grid,h2o_grid,     &
           ndx_sps,f_grid_filter,f_grid,filter_func,freqs,vel_z,       &
           h2o_corr,InDir,ld,fnd,lf,primag,mdb_hdr,mdb_rec,ier)

!  ===============================================================
!  Declaration of variables for sub-program: pfa_prep
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_atmos, n_lvls, no_pfa_ch, no_filt_pts, &
                           pfa_ch(:), no_int_frqs(:)

Integer(i4), INTENT(OUT) :: ier, ld, ndx_sps(:,:)

Real(r8), INTENT(IN) :: vel_z, z_grid(:), h2o_grid(:), t_grid(:), freqs(:)

Real(r8), INTENT(OUT) :: filter_func(:,:), f_grid(:,:), h2o_corr(:), &
                         f_grid_filter(:,:)

Character (LEN=*), INTENT(IN) :: InDir

Integer(i4), INTENT(IN) :: lf
Character (LEN=*), INTENT(IN) :: fnd

!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Real(r8), PARAMETER :: Tiny = epsilon(vel_z)
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, k, m, ch_i, sps_ind, sps_i, pb, no_sps, spectag, &
               j4, mch, band, spectags(20)

Real(r8) :: cp, ct, theta, sp, xlhs, xrhs, df, q, area, frq

Character (LEN=80) :: pqm_fnd, pqm_fni

Character (LEN=8) :: sps_name
Character (LEN=1) :: primag

type (atmos_comp) :: atmospheric(*)

type (pfa_slab) :: pfa_spectrum(6,*)
type (pfa_slab) :: pfa_spectrum_i(6,20)

type (eos_mdb_hdr), intent (OUT) :: mdb_hdr(*)
type (eos_mdb_rec), intent (OUT) :: mdb_rec(max_no_lines,*)

! begin code:

  ier = 0

!  Create H2O (18003) non-linear part:

  DO k = 1, n_lvls
    theta = 300.0 / t_grid(k)
    ct = theta**7.5
    cp = 31.6 * h2o_grid(k) * ct
    h2o_corr(k) = 1.0 + cp
  END DO

  IF(no_pfa_ch < 1) RETURN

!  Read in the EOS Beta database for all species involved

  pb = 0
  DO k = 1, 20
    spectags(k) = -1
  END DO

  DO ch_i = 1, no_pfa_ch

    mch = pfa_ch(ch_i)
    band = 1 + (mch - 1) / 15
    no_sps = pfa_spectrum(band,1)%no_sps

    DO  sps_i = 1, no_sps

      spectag = pfa_spectrum(band,sps_i)%sps_spectag

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

      IF(.NOT.atmospheric(sps_ind)%fwd_calc(band)) CYCLE

!  Now check if we loaded this Spectag already ..

      IF(pb > 0) THEN
        DO k = 1, pb
          IF(spectag == spectags(k)) GO TO 22
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
      spectags(pb) = spectag
      CALL read_eos_db(spectag,sps_name,mdb_hdr(pb),mdb_rec(1:,pb), &
                       ier)
      IF(ier /= 0) THEN
        PRINT *,'** Error in routine: read_EOS_db, Spectag:',spectag
        PRINT *,'   Called by routine: pfa_prep ..'
        RETURN
      END IF

 22   k = 1

    END DO

  END DO

!  Create PQM Master_database files from the master_database file given:

  pqm_fnd(1:)=' '
  pqm_fni(1:)=' '

!  ** Code for the PC

  k = -1
  pqm_fnd = '/Tamar/config/data/master_pqm_database.dat'
  pqm_fni = '/Tamar/config/data/master_pqm_database.ndx'

!  ** END Code for the PC

!  ** Code for the HIPP, SGI

!     k = lf
!     do while(k.gt.1)
!       k = k - 1
!       if(Fnd(k:k).eq.'/') then
!         Pqm_Fnd = Fnd(1:k)//'master_pqm_database.dat'
!         Pqm_Fni = Fnd(1:k)//'master_pqm_database.ndx'
!         k = -2
!       endif
!     end do

!  ** END Code for the HIPP, SGI

  IF(k >= 0) THEN
    ier =1
    PRINT *,' ** Error in Pfa_Prep subroutine **'
    PRINT *,'   Incorrect master database filename'
    RETURN
  END IF

!  Find the species index in the l2pc mixing ratio database:

  m = pb
  pb = -1
  DO ch_i = 1, no_pfa_ch

    mch = pfa_ch(ch_i)
    band = 1 + (mch - 1) / 15
    no_sps = pfa_spectrum(band,1)%no_sps

! Get total # of Aitken's points

    j4 = 4 * no_int_frqs(ch_i) + 1

    DO sps_i = 1, no_sps

      spectag = pfa_spectrum(band,sps_i)%sps_spectag
      spectags(sps_i) = spectag

!  Setup the spectrum record (structure), per band:

      j = 1
      k = spectag
      IF(k == 18999.OR.k == 28964.OR.k == 28965) j = -1

      IF(j >= 0 .AND. pb /= band) THEN

        k = -sps_i
        IF(sps_i > 1) k = sps_i

! ** DEBUG: Skip any pfa_spectrum data loading (any mdb_sw.f code ...)
!
!       IF(primag == 'p') THEN
!         CALL bin_pqm_intrp(band,k,n_lvls,spectags,pfa_spectrum, &
!              pfa_spectrum_i,vel_z,z_grid,t_grid,InDir,pqm_fnd,pqm_fni,ier)
!       ELSE
!         CALL bin_pqm_intrp(band,k,n_lvls,spectags,pfa_spectrum_i, &
!              pfa_spectrum,vel_z,z_grid,t_grid,InDir,pqm_fnd,pqm_fni,ier)
!       END IF
!
!       IF(ier /= 0) RETURN
!
! ** END DEBUG

        IF(spectag == 18003) THEN

!  Correct H2O (18003) pfa_spectrum entries, with the non-linear part:

          DO k = 1, n_lvls
            cp = h2o_corr(k)
            DO j = 1, maxrat
              sp = pfa_spectrum(band,sps_i)%yy(j,k)
              pfa_spectrum(band,sps_i)%yy(j,k) = sp * cp
            END DO
          END DO

        END IF

      END IF

! *** DEBUG
!     Overwrite the pfa_spectrum data with the appropriate EOS Database

      j = 0
      k = -1
      DO WHILE(j < m .AND. k < 1)
        j = j + 1
        IF(mdb_hdr(j)%spectag == spectag) k = j
      END DO

      IF(k > 0) THEN

        i = 0
        j = mdb_hdr(k)%no_lines
        pfa_spectrum(band,sps_i)%no_lines = j
        DO WHILE(i < j)
          i = i + 1
          pfa_spectrum(band,sps_i)%sps_n(i) = mdb_hdr(k)%n(i)
          pfa_spectrum(band,sps_i)%sps_w(i) = mdb_hdr(k)%w(i)
          pfa_spectrum(band,sps_i)%sps_el(i) = mdb_hdr(k)%el(i)
          pfa_spectrum(band,sps_i)%sps_ps(i) = mdb_hdr(k)%ps(i)
          pfa_spectrum(band,sps_i)%sps_n1(i) = mdb_hdr(k)%n1(i)
          pfa_spectrum(band,sps_i)%sps_n2(i) = mdb_hdr(k)%n2(i)
          pfa_spectrum(band,sps_i)%sps_v0(i) = mdb_hdr(k)%v0(i)
          pfa_spectrum(band,sps_i)%sps_delta(i) = mdb_hdr(k)%delta(i)
          pfa_spectrum(band,sps_i)%sps_gamma(i) = mdb_hdr(k)%gamma(i)
        END DO

      END IF

! *** END DEBUG

    END DO

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

!  Get Aitken's grid points:

    df = (xrhs - xlhs) / (j4 - 1)
    DO j = 1, j4
      q = xlhs + (j - 1) * df
      f_grid(j,ch_i) = frq + q
    END DO

    f_grid(j4+1,ch_i) = frq                 ! *** DEBUG

  END DO                         ! On ch_i

  99  IF(ier /= 0) THEN
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
          mdb_zero_rec%log_beta_intrp(ii,jj,kk) = ' '
        END DO
      END DO
    END DO

  END IF

  i = 0
  j = 0
  ier = 0
  sps_name(1:) = ' '
  DO WHILE(j < no_sps.AND.i < 1)
    j = j + 1
    IF(spectags(j) == spectag) THEN
      i = j
      j = 20
    END IF
  END DO

  IF(i < 1) THEN
    io = -1
    GO TO 99
  END IF

  sps_name = names(i)

  datdir(1:) = ' '
  datdir= '/home/zvi/data/'             ! PC
! DatDir= '/zvi/eos/data/'        ! HIPP, SGI

  fdt(1:) = ' '
  fhd(1:) = ' '
  j = LEN_TRIM(datdir)
  i = LEN_TRIM(sps_name)
  fdt = datdir(1:j)//sps_name(1:i)//'_eosmdb.dat'
  CALL strlwr(fdt)
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
