!
module GET_CS_M
  use L2PCDim, only: NLVL
  use I_HUNT_M, only: HUNT
  use EOS_MDB, only: CS_MNF => MAX_FREQ
  use MDBETA, only: max_no_zeta, MNP => no_t_phi
  use MLSCommon, only: I4, R4, R8
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------

SUBROUTINE get_cs(primag, spectag, n_lvls, p_indx, ncs, mdb_pres, &
                  mdb_temp, mdb_freq, cs, fnd, ch1, ier)

!  ===============================================================
!  Declaration of variables for sub-program: get_cs
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: spectag, ch1, p_indx(:)

Integer(i4), INTENT(IN OUT) :: n_lvls, ncs, ier

Real(r8), INTENT(OUT) :: mdb_freq(:), mdb_pres(:), mdb_temp(:), cs(:,:,:)

Character (LEN=*), INTENT(IN) :: fnd, primag
!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Integer(i4), PARAMETER :: ksps = 20
!  ----------------
!  Local variables:
!  ----------------
Integer(i4), save :: ns = 0
Integer(i4), save :: lf = 38
Integer(i4), save :: spectags(ksps) = 0

Integer(i4) :: js, j, j_obs

! begin code:

  j_obs = n_lvls
  IF(j_obs < 1) THEN
    ns = -1
    n_lvls = -j_obs
    j_obs = n_lvls
  END IF

  IF(ns < 1) THEN
    ns = 1
    js = -1
    spectags(1) = spectag
  ELSE
    j = 0
    DO WHILE(j < ns)
      j = j + 1
      IF(spectags(j) == spectag) THEN
        js = j
        j = ns + 10
      END IF
    END DO
    IF(j < ns+2) THEN
      ns = ns + 1
      spectags(ns) = spectag
      js = ns
    END IF
  END IF

! Read BINARY Master database:

  CALL read_eos_mdb(primag,js,n_lvls,spectags,p_indx,cs_mnf, &
                    ncs,mdb_pres,mdb_temp,mdb_freq,cs,fnd,lf,ch1,ier)
  IF(ier > 0) THEN
    PRINT *,'** Error in "read_eos_mdb" subroutine **'
    CALL errmsg(' ',ier)
  END IF

  RETURN
END SUBROUTINE get_cs

!-----------------------------------------------------------------------
!  Reading the "Binary"  Master DataBase

SUBROUTINE read_eos_mdb(api,js,n_lvls,spectags,p_indx,cs_mnf, &
           ncs,mdb_pres,mdb_temp,mdb_freq,cs,fnd,lf,ch1,ier)

Character (len=*), intent(in) :: fnd, api

Integer(i4), intent(in) :: spectags(*), p_indx(*), n_lvls, cs_mnf, ch1, lf

Integer(i4), intent(in out) :: js
Integer(i4), intent(out)    :: ncs, Ier

Real(r8), INTENT(OUT) :: mdb_freq(:), mdb_pres(:), mdb_temp(:), cs(:,:,:)

! INCLUDE 'inc/mdbeta_org.inc'

 Integer(i4), Parameter :: org_mnf = 90  ! Original Max. number of frequencies
 type ORG_MDB_BETA_REC
   Integer(i4) :: no_freq, spectag
   Real(r8)    :: Freq_p(org_mnf), Freq_i(org_mnf)
   real(r4)    :: beta_p(max_no_zeta,mnp,org_mnf)
   real(r4)    :: beta_i(max_no_zeta,mnp,org_mnf)
 end type ORG_MDB_BETA_REC

INTEGER(i4), PARAMETER :: ksps=20

REAL(r8), save :: eos_press(max_no_zeta), eos_temp(max_no_zeta)
Integer(i4), save :: nsr,n_eos,sname(ksps),erc(ksps)

INTEGER(i4) :: i, j, k, m, io, spec, irec, kcs, shift, spectag
REAL(r8) :: z, t

type (ORG_MDB_BETA_REC) :: org_mdbeta

!
! Begin code:
!
  ier = 0
!
  IF(js < 0) THEN
!
    js = -js
    inquire (iolength=j) org_mdbeta
    CLOSE(lf,IOSTAT=io)
    OPEN(lf,FILE=fnd,FORM='UNFORMATTED',STATUS='OLD',action='READ', RECL=j, &
         ACCESS='DIRECT',IOSTAT=io)
    IF(io /= 0) THEN
      ier = 1
      k = LEN_TRIM(fnd)
      PRINT *,'** File: ',fnd(1:k)
      CALL errmsg(' ',io)
      RETURN
    END IF
!
!  Load the Spectags and their records number in the index file, into memory:
!
    READ(lf,REC=1,IOSTAT=io) org_mdbeta
    IF(io /= 0) THEN
      ier = 1
      k = LEN_TRIM(fnd)
      PRINT *,'** File: ',fnd(1:k)
      CALL errmsg(' ',io)
      RETURN
    END IF
!
    n_eos = org_mdbeta%no_freq
    nsr = org_mdbeta%spectag
!
    DO k = 1, nsr
      t = org_mdbeta%beta_p(k,1,1)
      z = org_mdbeta%beta_p(k,2,1)
      erc(k) = INT(z)
      sname(k) = INT(t)
    END DO
!
    DO k = 1, n_eos
      eos_temp(k) = DBLE(org_mdbeta%beta_p(k,4,1))
      eos_press(k) = DBLE(org_mdbeta%beta_p(k,3,1))
    END DO
!
  END IF
!
  DO i = 1, n_lvls
    j = p_indx(i)
    mdb_temp(i) = eos_temp(j)
    mdb_pres(i) = eos_press(j)
  END DO
!
  spectag = spectags(js)
  CALL hunt(spectag,sname,nsr,k,i)
  DO WHILE(sname(k) /= spectag.AND.k < nsr)
    k = k + 1
  END DO
!
  IF(sname(k) /= spectag) THEN
    ier = 1
    PRINT *,'** Spectag: ',spectag,' Not in the Header'
    RETURN
  END IF
!
  irec = erc(k)
  READ(lf,REC=irec,IOSTAT=io) org_mdbeta
  IF(io /= 0) THEN
    ier = 1
    k = LEN_TRIM(fnd)
    PRINT *,'** Error reading: ',fnd(1:k)
    CALL errmsg(' ',io)
    RETURN
  END IF
!
  spec = org_mdbeta%spectag
  IF(spec /= spectag) THEN
    ier = 1
    k = LEN_TRIM(fnd)
    PRINT *,'** Spectag: ',spectag,' Not in file: ',fnd(1:k)
    RETURN
  END IF
!
  ncs = 0
  kcs = MIN(org_mdbeta%no_freq,cs_mnf)
  shift = ch1 + org_mdbeta%no_freq - org_mnf
!
  DO k = 1, kcs
    mdb_freq(k) = 0.0D0
  END DO
!
  DO WHILE(ncs < kcs)
!
    ncs = ncs + 1
    m = ncs + shift - 1
!
    IF(api == 'p') THEN
!
      mdb_freq(ncs) = org_mdbeta%freq_p(m)
      DO k = 1, n_lvls
        j = p_indx(k)
        DO i = 1, mnp
          cs(k,i,ncs) = DBLE(org_mdbeta%beta_p(j,i,m))
        END DO
      END DO
!
    ELSE
!
      mdb_freq(ncs) = org_mdbeta%freq_i(m)
      DO k = 1, n_lvls
        j = p_indx(k)
        DO i = 1, mnp
          cs(k,i,ncs) = DBLE(org_mdbeta%beta_i(j,i,m))
        END DO
      END DO
!
    END IF
!
  END DO
!
  RETURN
END SUBROUTINE read_eos_mdb
end module GET_CS_M
! $Log$
! Revision 1.1  2000/06/21 21:56:13  zvi
! First version D.P.
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
