
!-----------------------------------------------------------------------

SUBROUTINE bin_pqm_intrp(ib,js,n_lvls,spectags,pfa_spectrum_p,  &
           pfa_spectrum_i,vel_z,p_grid,t_grid,indir,fnd,fni,ier)


INTEGER*4, INTENT(IN OUT)                :: ib
INTEGER*4, INTENT(OUT)                   :: js
INTEGER*4, INTENT(IN)                    :: n_lvls
NO TYPE, INTENT(IN)                      :: spectags
NO TYPE, INTENT(IN OUT)                  :: pfa_spectr
NO TYPE, INTENT(IN OUT)                  :: pfa_spectr
REAL*4, INTENT(IN OUT)                   :: vel_z
NO TYPE, INTENT(IN)                      :: p_grid
NO TYPE, INTENT(IN)                      :: t_grid
CHARACTER (LEN=*), INTENT(IN)            :: indir
CHARACTER (LEN=*), INTENT(IN OUT)        :: fnd
CHARACTER (LEN=*), INTENT(IN OUT)        :: fni
INTEGER*4, INTENT(IN OUT)                :: ier
IMPLICIT NONE

INCLUDE '/hpdsk/config/inc_v51/l2pcdim.inc'
INCLUDE '/hpdsk/config/inc_v51/pqm.inc'
INCLUDE '/hpdsk/config/inc_v51/units.inc'

INTEGER*4
INTEGER*4, PARAMETER :: maxh=400
INTEGER*4, PARAMETER :: mz=35


CHARACTER (LEN=80) :: pvsh
CHARACTER (LEN=8) :: aname

REAL*8 dp,q,dvz,yv(5),velcor_array(5),tdep(5),td

INTEGER*4 its,kv,n_hlvl, spectags(*), i,j,i1,i2,l,nf,ni

REAL*4 t,p,tmp_lo,tmp_md,tmp_hi,y,r,v_lo(maxrat,5),v_md(maxrat,5), v_hi(maxrat,5),p1,p2

REAL*4  p_grid(*),t_grid(*)

REAL*4 dumh(maxh),plog(maxh),tlog(maxh),temp_step,y_lo,y_md,y_hi

INTEGER*4 npl_p,npl_i
REAL*8 freq_p(mz),freq_i(mz)
REAL*4 w_p(mz),n_p(mz),el_p(mz),n1_p(mz),n2_p(mz),log_i_p(mz), gamma_p(mz),delta_p(mz),w_i(mz),n_i(mz),el_i(mz),n1_i(mz), n2_i(mz),log_i_i(m
ps_p(mz),ps_i(mz)

INTEGER*4 nrat_p(maxh,5),nrat_i(maxh,5)

REAL*8 xx_p(maxrat,maxh,5), yy_lo_p(maxrat,maxh,5), yy_md_p(maxrat,maxh,5), yy_hi_p(maxrat,maxh,5)

REAL*8 xx_i(maxrat,maxh,5), yy_lo_i(maxrat,maxh,5), yy_md_i(maxrat,maxh,5), yy_hi_i(maxrat,maxh,5)

REAL*8 yy(maxrat,5)

record/pfa_slab/pfa_spectrum_p(6,*)
record/pfa_slab/pfa_spectrum_i(6,*)

SAVE n_hlvl,plog,tlog,velcor_array,temp_step,nf,ni,its

DATA temp_step/-1.0/
DATA velcor_array/-0.4D0,-0.2D0,0.0D0,0.2D0,0.4D0/

IF(temp_step < 1.0) THEN

  pvsh(1:)=' '
  j = Len_TRIM(InDir)
  pvsh = indir(1:j)//'p_vs_h_md.dat'
  CALL get_pvsh(pvsh,dumh,tlog,plog,n_hlvl,ier)
  IF(ier /= 0) RETURN

  temp_step = 20.0
  IF(js > 0) js = -js
  its = INT(temp_step+0.5)

  nf = mdb_pqm_unit_dat
  ni = mdb_pqm_unit_ndx

END IF

CALL read_bin_pqm(ib,js,spectags,aname,q_log,npl_p,npl_i,  &
    freq_p,w_p,ps_p,n_p,el_p,n1_p,n2_p,log_i_p,gamma_p,delta_p,  &
    freq_i,w_i,ps_i,n_i,el_i,n1_i,n2_i,log_i_i,gamma_i,delta_i,  &
    nrat_p,xx_p,yy_lo_p,yy_md_p,yy_hi_p,nrat_i,xx_i,  &
    yy_lo_i,yy_md_i,yy_hi_i,indir,fnd,fni,nf,ni,ier)
IF(ier /= 0) RETURN

!  Primary ...

pfa_spectrum_p(ib,js).no_lines = npl_p
pfa_spectrum_p(ib,js).sps_name = aname
pfa_spectrum_p(ib,js).sps_spectag = spectags(js)
DO l = 1, npl_p
  pfa_spectrum_p(ib,js).sps_v0(l) = freq_p(l)
  pfa_spectrum_p(ib,js).sps_el(l) = el_p(l)
  pfa_spectrum_p(ib,js).sps_str(l) = log_i_p(l)
  pfa_spectrum_p(ib,js).sps_w(l) = w_p(l)
  pfa_spectrum_p(ib,js).sps_ps(l) = ps_p(l)
  pfa_spectrum_p(ib,js).sps_n(l) = n_p(l)
  pfa_spectrum_p(ib,js).sps_n1(l) = n1_p(l)
  pfa_spectrum_p(ib,js).sps_n2(l) = n2_p(l)
  pfa_spectrum_p(ib,js).sps_gamma(l) = gamma_p(l)
  pfa_spectrum_p(ib,js).sps_delta(l) = delta_p(l)
  DO i = 1, 3
    pfa_spectrum_p(ib,js).sps_part(i,l) = q_log(i)
  END DO
END DO

DO i = 1, nlvl
  pfa_spectrum_p(ib,js).varm(i) = 0.0D0
  pfa_spectrum_p(ib,js).nrat(i) = maxrat
  DO j = 1, maxrat
    pfa_spectrum_p(ib,js).xx(j,i) = 0.0
    pfa_spectrum_p(ib,js).yy(j,i) = 0.0
  END DO
END DO

!  Image ..

pfa_spectrum_i(ib,js).no_lines = npl_i
pfa_spectrum_i(ib,js).sps_name = aname
pfa_spectrum_i(ib,js).sps_spectag = spectags(js)
DO l = 1, npl_i
  pfa_spectrum_i(ib,js).sps_v0(l) = freq_i(l)
  pfa_spectrum_i(ib,js).sps_el(l) = el_i(l)
  pfa_spectrum_i(ib,js).sps_str(l) = log_i_i(l)
  pfa_spectrum_i(ib,js).sps_w(l) = w_i(l)
  pfa_spectrum_i(ib,js).sps_ps(l) = ps_i(l)
  pfa_spectrum_i(ib,js).sps_n(l) = n_i(l)
  pfa_spectrum_i(ib,js).sps_n1(l) = n1_i(l)
  pfa_spectrum_i(ib,js).sps_n2(l) = n2_i(l)
  pfa_spectrum_i(ib,js).sps_gamma(l) = gamma_i(l)
  pfa_spectrum_i(ib,js).sps_delta(l) = delta_i(l)
  DO i = 1, 3
    pfa_spectrum_i(ib,js).sps_part(i,l) = q_log(i)
  END DO
END DO

DO i = 1, nlvl
  pfa_spectrum_i(ib,js).varm(i) = 0.0D0
  pfa_spectrum_i(ib,js).nrat(i) = maxrat
  DO j = 1, maxrat
    pfa_spectrum_i(ib,js).xx(j,i) = 0.0
    pfa_spectrum_i(ib,js).yy(j,i) = 0.0
  END DO
END DO

dvz = DBLE(vel_z)

DO  i = 1, n_lvls

  t = t_grid(i)
  p = p_grid(i)
  CALL search(p,plog,n_hlvl,i1,i2)

  p1 = plog(i1)
  p2 = plog(i2)

  IF(i2 > i1) THEN

    dp = p2 - p1
    q = (tlog(i2) - tlog(i1)) / dp
    tmp_md = tlog(i1) + q * (p - p1)

  ELSE

    tmp_md = tlog(i1)

  END IF

  tmp_lo = tmp_md - temp_step
  tmp_hi = tmp_md + temp_step

!  Primary  ..

  kv = pfa_spectrum_p(ib,js).nrat(i)
  IF(kv < 4) GO TO 40

  IF(i2 > i1) THEN

    DO j = 1, 5
      DO l = 1, kv
        CALL exp_intrp(p1,p2,p,SNGL(yy_lo_p(l,i1,j)),  &
            SNGL(yy_lo_p(l,i2,j)),v_lo(l,j))
        CALL exp_intrp(p1,p2,p,SNGL(yy_md_p(l,i1,j)),  &
            SNGL(yy_md_p(l,i2,j)),v_md(l,j))
        CALL exp_intrp(p1,p2,p,SNGL(yy_hi_p(l,i1,j)),  &
            SNGL(yy_hi_p(l,i2,j)),v_hi(l,j))
      END DO
    END DO

  ELSE

    DO j = 1, 5
      DO l = 1, kv
        v_lo(l,j) = yy_lo_p(l,i1,j)
        v_md(l,j) = yy_md_p(l,i1,j)
        v_hi(l,j) = yy_hi_p(l,i1,j)
      END DO
    END DO

  END IF

  DO j = 1, 5          ! For fixed velocity, interpolate in temp.
    td = 0.0D0
    DO l = 1, kv
      y_lo = v_lo(l,j)
      y_md = v_md(l,j)
      y_hi = v_hi(l,j)
      IF(y_lo*y_md*y_hi <= 0.0) THEN
        pfa_spectrum_p(ib,js).nrat(i) = 0
        GO TO 40
      END IF
      CALL power_intrp(tmp_lo,tmp_md,tmp_hi,t,y_lo,y_md,y_hi,y,r)
      yy(l,j) = y
      td = td + r
    END DO
    tdep(j) = td / kv
  END DO

  CALL d_polint(velcor_array,tdep,5,dvz,td)
  pfa_spectrum_p(ib,js).varm(i) = td

  DO l = 1, kv         ! For fixed coeff. #, interpolate on velocity
    DO j = 1, 5
      yv(j) = yy(l,j)
    END DO
    CALL d_polint(velcor_array,yv,5,dvz,q)
    pfa_spectrum_p(ib,js).yy(l,i) = q
    pfa_spectrum_p(ib,js).xx(l,i) = xx_p(l,i1,1)
  END DO

!  Image  ..

  40     kv = pfa_spectrum_i(ib,js).nrat(i)
  IF(kv < 4) CYCLE

  IF(i2 > i1) THEN

    DO j = 1, 5
      DO l = 1, kv
        CALL exp_intrp(p1,p2,p,SNGL(yy_lo_i(l,i1,j)),  &
            SNGL(yy_lo_i(l,i2,j)),v_lo(l,j))
        CALL exp_intrp(p1,p2,p,SNGL(yy_md_i(l,i1,j)),  &
            SNGL(yy_md_i(l,i2,j)),v_md(l,j))
        CALL exp_intrp(p1,p2,p,SNGL(yy_hi_i(l,i1,j)),  &
            SNGL(yy_hi_i(l,i2,j)),v_hi(l,j))
      END DO
    END DO

  ELSE

    DO j = 1, 5
      DO l = 1, kv
        v_lo(l,j) = yy_lo_i(l,i1,j)
        v_md(l,j) = yy_md_i(l,i1,j)
        v_hi(l,j) = yy_hi_i(l,i1,j)
      END DO
    END DO

  END IF

  DO j = 1, 5          ! For fixed velocity, interpolate in temp.
    td = 0.0D0
    DO l = 1, kv
      y_lo = v_lo(l,j)
      y_md = v_md(l,j)
      y_hi = v_hi(l,j)
      IF(y_lo*y_md*y_hi <= 0.0) THEN
        pfa_spectrum_i(ib,js).nrat(i) = 0
        CYCLE
      END IF
      CALL power_intrp(tmp_lo,tmp_md,tmp_hi,t,y_lo,y_md,y_hi,y,r)
      yy(l,j) = y
      td = td + r
    END DO
    tdep(j) = td / kv
  END DO

  CALL d_polint(velcor_array,tdep,5,dvz,td)
  pfa_spectrum_i(ib,js).varm(i) = td

  DO l = 1, kv         ! For fixed coeff. #, interpolate on velocity
    DO j = 1, 5
      yv(j) = yy(l,j)
    END DO
    CALL d_polint(velcor_array,yv,5,dvz,q)
    pfa_spectrum_i(ib,js).yy(l,i) = q
    pfa_spectrum_i(ib,js).xx(l,i) = xx_i(l,i1,1)
  END DO

END DO

RETURN
END SUBROUTINE bin_pqm_intrp

!-----------------------------------------------------------------------
!  Reading the "Shorter" binary PQM DataBase

SUBROUTINE read_bin_pqm(band,js,spectags,aname,q_log,npl_p,npl_i,  &
    freq_p,w_p,ps_p,n_p,el_p,n1_p,n2_p,log_i_p,gamma_p,delta_p,  &
    freq_i,w_i,ps_i,n_i,el_i,n1_i,n2_i,log_i_i,gamma_i,delta_i,  &
    nrat_p,xx_p,yy_lo_p,yy_md_p,yy_hi_p,nrat_i,xx_i,  &
    yy_lo_i,yy_md_i,yy_hi_i,indir,fnd,fni,nf,ni,ier)


NO TYPE, INTENT(IN OUT)                  :: band
NO TYPE, INTENT(IN OUT)                  :: js
NO TYPE, INTENT(IN)                      :: spectags
CHARACTER (LEN=8), INTENT(OUT)           :: aname
NO TYPE, INTENT(OUT)                     :: q_log
INTEGER*4, INTENT(OUT)                   :: npl_p
INTEGER*4, INTENT(OUT)                   :: npl_i
NO TYPE, INTENT(OUT)                     :: freq_p
NO TYPE, INTENT(OUT)                     :: w_p
NO TYPE, INTENT(OUT)                     :: ps_p
NO TYPE, INTENT(OUT)                     :: n_p
NO TYPE, INTENT(OUT)                     :: el_p
NO TYPE, INTENT(OUT)                     :: n1_p
NO TYPE, INTENT(OUT)                     :: n2_p
NO TYPE, INTENT(OUT)                     :: log_i_p
NO TYPE, INTENT(OUT)                     :: gamma_p
NO TYPE, INTENT(OUT)                     :: delta_p
NO TYPE, INTENT(OUT)                     :: freq_i
NO TYPE, INTENT(OUT)                     :: w_i
NO TYPE, INTENT(OUT)                     :: ps_i
NO TYPE, INTENT(OUT)                     :: n_i
NO TYPE, INTENT(OUT)                     :: el_i
NO TYPE, INTENT(OUT)                     :: n1_i
NO TYPE, INTENT(OUT)                     :: n2_i
NO TYPE, INTENT(OUT)                     :: log_i_i
NO TYPE, INTENT(OUT)                     :: gamma_i
NO TYPE, INTENT(OUT)                     :: delta_i
NO TYPE, INTENT(OUT)                     :: nrat_p
NO TYPE, INTENT(OUT)                     :: xx_p
NO TYPE, INTENT(OUT)                     :: yy_lo_p
NO TYPE, INTENT(OUT)                     :: yy_md_p
NO TYPE, INTENT(OUT)                     :: yy_hi_p
NO TYPE, INTENT(OUT)                     :: nrat_i
NO TYPE, INTENT(OUT)                     :: xx_i
NO TYPE, INTENT(OUT)                     :: yy_lo_i
NO TYPE, INTENT(OUT)                     :: yy_md_i
NO TYPE, INTENT(OUT)                     :: yy_hi_i
CHARACTER (LEN=*), INTENT(IN OUT)        :: indir
CHARACTER (LEN=*), INTENT(OUT)           :: fnd
CHARACTER (LEN=*), INTENT(OUT)           :: fni
INTEGER*4, INTENT(IN OUT)                :: nf
INTEGER*4, INTENT(IN OUT)                :: ni
INTEGER*4, INTENT(OUT)                   :: ier
IMPLICIT NONE

INCLUDE '/hpdsk/config/inc_v51/l2pcdim.inc'
INCLUDE '/hpdsk/config/inc_v51/mdb.inc'



INTEGER*4  spectags(*),band,spectag,js
INTEGER*4

INTEGER*4, PARAMETER :: ksps=20



REAL*8 freq_p(*),freq_i(*)
REAL*4 w_p(*),n_p(*),el_p(*),n1_p(*),n2_p(*),log_i_p(*), gamma_p(*),delta_p(*),w_i(*),n_i(*),el_i(*),n1_i(*), n2_i(*),log_i_i(*),gamma_i(*),
ps_p(*),ps_i(*)

INTEGER*4 nrat_p(maxh,*),nrat_i(maxh,*)

REAL*8 xx_p(mrat,maxh,*), yy_lo_p(mrat,maxh,*), yy_md_p(mrat,maxh,*), yy_hi_p(mrat,maxh,*)

REAL*8 xx_i(mrat,maxh,*), yy_lo_i(mrat,maxh,*), yy_md_i(mrat,maxh,*), yy_hi_i(mrat,maxh,*)

INTEGER*4 na,ivel,i,j,io,ibd,spec,rectot,k,irec,jrec,nb,  &
    strlen,preco(ksps),pname(ksps),msr,sbd(6),src(6)

record/mdb_pqm_rec/pqm(2)

SAVE rectot,msr,preco,pname

ier = 0

IF(js < 0) THEN

  j = js
  js = -j

  j = 2 * sizeof(pqm(1))
  CLOSE(nf,IOSTAT=io)
  OPEN(nf,FILE=fnd,FORM='UNFORMATTED',STATUS='OLD',RECL=j,  &
      readonly,ACCESS='DIRECT',IOSTAT=io)
  IF(io /= 0) THEN
    ier = 1
    k = strlen(fnd)
    PRINT *,'** File: ',fnd(1:k)
    CALL errmsg(' ',io)
    GO TO 99
  END IF

  CLOSE(ni,IOSTAT=io)
  i = MAX(2*ksps+2,14)
  j = i * sizeof(i)
  OPEN(ni,FILE=fni,FORM='UNFORMATTED',STATUS='OLD',RECL=j,  &
      readonly,ACCESS='DIRECT',IOSTAT=io)
  IF(io /= 0) THEN
    ier = 1
    k = strlen(fni)
    PRINT *,'** File: ',fni(1:k)
    CALL errmsg(' ',io)
    GO TO 99
  END IF

  READ(ni,REC=1) rectot,msr,(pname(i),i=1,msr),(preco(j),j=1,msr)

END IF

spectag = spectags(js)
CALL i_search(spectag,pname,msr,k,i)
DO WHILE(pname(k) /= spectag.AND.k < msr)
  k = k + 1
END DO

IF(pname(k) /= spectag) THEN
  ier = 1
  k = strlen(fni)
  PRINT *,'** Spectag: ',spectag,' Not in Header: ',fni(1:k)
  GO TO 99
ELSE
  jrec = preco(k)
  READ(ni,REC=jrec) spec,nb,(sbd(i),i=1,nb),(src(j),j=1,nb)
  IF(spec /= spectag) THEN
    ier = 1
    k = strlen(fni)
    PRINT *,'** Spectag: ',spectag,' Not in file: ',fni(1:k)
    GO TO 99
  END IF
END IF

i = 1
DO WHILE(i < nb.AND.sbd(i) /= band)
  i = i + 1
END DO

IF(sbd(i) /= band) THEN
  ier = 1
  k = strlen(fni)
  PRINT *,'** For Spectag: ',spectag
  PRINT *,'** Band: ',band,'  Not found in: ',fni(1:k)
  GO TO 99
END IF

irec = src(i)
READ(nf,REC=irec,IOSTAT=io) pqm(1),pqm(2)
IF(io /= 0) THEN
  ier = 1
  k = strlen(fnd)
  PRINT *,'** Error reading: ',fnd(1:k)
  CALL errmsg(' ',io)
  GO TO 99
END IF

spec = pqm(1).spectag
IF(spec /= spectag) THEN
  ier = 1
  k = strlen(fnd)
  PRINT *,'** Spectag: ',spectag,'  Not found in: ',fnd(1:k)
  GO TO 99
END IF

ibd = pqm(1).band
IF(ibd /= band) THEN
  ier = 1
  k = strlen(fnd)
  PRINT *,'** For Spectag: ',spectag
  PRINT *,'** Band: ',band,'  Not found in: ',fnd(1:k)
  GO TO 99
END IF

aname = pqm(1).NAMe
npl_p = pqm(1).no_cat_lines
npl_i = pqm(2).no_cat_lines
q_log(1) = pqm(1).qlog(1)
q_log(2) = pqm(1).qlog(2)
q_log(3) = pqm(1).qlog(3)

IF(npl_p > 0) THEN
  DO i = 1, npl_p
    w_p(i) = pqm(1).w(i)
    ps_p(i) = pqm(1).ps(i)
    n_p(i) = pqm(1).n(i)
    el_p(i) = pqm(1).el(i)
    n1_p(i) = pqm(1).n1(i)
    n2_p(i) = pqm(1).n2(i)
    freq_p(i) = pqm(1).freq(i)
    log_i_p(i) = pqm(1).log_i(i)
    gamma_p(i) = pqm(1).gamma(i)
    delta_p(i) = pqm(1).delta(i)
  END DO
END IF

IF(npl_i > 0) THEN
  DO i = 1, npl_i
    w_i(i) = pqm(2).w(i)
    ps_i(i) = pqm(2).ps(i)
    n_i(i) = pqm(2).n(i)
    el_i(i) = pqm(2).el(i)
    n1_i(i) = pqm(2).n1(i)
    n2_i(i) = pqm(2).n2(i)
    freq_i(i) = pqm(2).freq(i)
    log_i_i(i) = pqm(2).log_i(i)
    gamma_i(i) = pqm(2).gamma(i)
    delta_i(i) = pqm(2).delta(i)
  END DO
END IF

DO i = 1, maxh
  DO ivel = 1, 5
    nrat_p(i,ivel) = 0
    nrat_i(i,ivel) = 0
    DO k = 1, mrat
      xx_p(k,i,ivel)    = 0.0D0
      xx_i(k,i,ivel)    = 0.0D0
      yy_lo_p(k,i,ivel) = 0.0D0
      yy_lo_i(k,i,ivel) = 0.0D0
      yy_md_p(k,i,ivel) = 0.0D0
      yy_md_i(k,i,ivel) = 0.0D0
      yy_hi_p(k,i,ivel) = 0.0D0
      yy_hi_i(k,i,ivel) = 0.0D0
    END DO
  END DO
END DO

DO i = 1, maxh
  DO ivel = 1, 5
    na = pqm(1).nrat(i,ivel)
    IF(na > 0) THEN
      nrat_p(i,ivel) = na
      DO k = 1, na
        xx_p(k,i,ivel)    = pqm(1).xx(k,i,ivel)
        yy_lo_p(k,i,ivel) = pqm(1).yy_lo(k,i,ivel)
        yy_md_p(k,i,ivel) = pqm(1).yy_md(k,i,ivel)
        yy_hi_p(k,i,ivel) = pqm(1).yy_hi(k,i,ivel)
      END DO
    END IF
    na = pqm(2).nrat(i,ivel)
    IF(na > 0) THEN
      nrat_i(i,ivel) = na
      DO k = 1, na
        xx_i(k,i,ivel)    = pqm(2).xx(k,i,ivel)
        yy_lo_i(k,i,ivel) = pqm(2).yy_lo(k,i,ivel)
        yy_md_i(k,i,ivel) = pqm(2).yy_md(k,i,ivel)
        yy_hi_i(k,i,ivel) = pqm(2).yy_hi(k,i,ivel)
      END DO
    END IF
  END DO
END DO

RETURN

99   CLOSE(nf,IOSTAT=io)
CLOSE(ni,IOSTAT=io)

RETURN
END SUBROUTINE read_bin_pqm

!---------------------------------------------------------------------
! Reading any p_vs_h file:

SUBROUTINE get_pvsh(pvsh,h_grid,t_grid,p_grid,no_hts,ier)


CHARACTER (LEN=*), INTENT(OUT)           :: pvsh
NO TYPE, INTENT(OUT)                     :: h_grid
NO TYPE, INTENT(OUT)                     :: t_grid
NO TYPE, INTENT(OUT)                     :: p_grid
INTEGER*4, INTENT(OUT)                   :: no_hts
INTEGER*4, INTENT(OUT)                   :: ier
IMPLICIT NONE

INCLUDE '/hpdsk/config/inc_v51/units.inc'





REAL*4 t_grid(*),p_grid(*),h_grid(*)

CHARACTER (LEN=8) :: ax
INTEGER*4 io,uph
REAL*4 h,t,p,plog,sgn

ier = 0
uph = p_vs_h_unit
CLOSE(uph,IOSTAT=io)
OPEN(uph,FILE=pvsh,STATUS='OLD',IOSTAT=io)
IF(io /= 0) THEN
  ier = 1
  PRINT *,' ** Error in subroutine: Get_PvsH'
  PRINT *,'    Couldn''T OPEN FILE: ',pvsh
  RETURN
END IF

no_hts = 0
sgn = 1.0
READ(uph,'(A)',IOSTAT=io) ax
IF(io == 0) READ(uph,'(A)',IOSTAT=io) ax
DO WHILE(io == 0)
  READ(uph,*,IOSTAT=io,END=10) h,t,p,plog
  IF(io /= 0) THEN
    ier = 1
    PRINT *,' ** Error in subroutine: Get_PvsH'
    PRINT *,'    Error while reading file: ',pvsh
    GO TO 10
  END IF
  IF(no_hts == 0.AND.plog > 0.0) sgn = -1.0
  no_hts = no_hts + 1
  h_grid(no_hts) = h
  t_grid(no_hts) = t
  p_grid(no_hts) = plog * sgn
END DO

10  CLOSE(uph,IOSTAT=io)

RETURN
END SUBROUTINE get_pvsh

!---------------------------------------------------------------------
! Exponential interpolation subroutine

SUBROUTINE exp_intrp(h1,h2,h,a1,a2,a)


REAL*4, INTENT(IN OUT)                   :: h1
REAL*4, INTENT(IN OUT)                   :: h2
REAL*4, INTENT(IN OUT)                   :: h
REAL*4, INTENT(IN)                       :: a1
REAL*4, INTENT(IN)                       :: a2
REAL*4, INTENT(OUT)                      :: a
IMPLICIT NONE

REAL*4  r

IF(ABS(h2-h1) < 1.0E-8) THEN
  a = 0.5 * (a1 + a2)
  RETURN
ELSE IF(ABS(h-h1) < 5.0E-9) THEN
  a = a1
  RETURN
ELSE IF(ABS(h-h2) < 5.0E-9) THEN
  a = a2
  RETURN
END IF

r = -1.0
IF(ABS(a1) > 1.0E-11) r = a2 / a1

IF(r > 0.0) THEN

  a = a1 * EXP(LOG(r)*(h-h1)/(h2-h1))

ELSE

! This is a sign change transition revert to linear interpolation

  a = a1 + (a2 - a1) * (h - h1) / (h2 - h1)

END IF

RETURN
END SUBROUTINE exp_intrp

!--------------------------------------------------------------------
! Power interpolation subroutine

SUBROUTINE power_intrp(t0,t1,t2,t,cs0,cs1,cs2,cs,tp)

x0 = LOG(t0)
y0 = LOG(cs0)

x1 = LOG(t1)
y1 = LOG(cs1)

x2 = LOG(t2)
y2 = LOG(cs2)

tp0 = (y1 - y0) / (x1 - x0)
tp2 = (y2 - y1) / (x2 - x1)

IF(t <= t0) THEN
  tp = tp0
  IF(ABS(tp) > 15.0) tp = 15.0 * SIGN(1.0,tp0)
  cs = cs0 * ((t/t0)**tp)
ELSE IF(t >= t2) THEN
  tp = tp2
  IF(ABS(tp) > 15.0) tp = 15.0 * SIGN(1.0,tp2)
  cs = cs2 * ((t/t2)**tp)
ELSE
  x = LOG(t)
  r = (x - x0) / (x2 - x0)
  tp = (1.0 - r) * tp0 + r * tp2
  IF(ABS(tp) > 15.0) tp = 15.0 * SIGN(1.0,tp2)
  cs = cs1 * ((t/t1)**tp)
END IF

RETURN
END SUBROUTINE power_intrp

!-----------------------------------------------------------------------
!  Integer sorting routine (Quick Sort) for MAIN array & attached array

SUBROUTINE isort2(n,iary,jary)



INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: iary
INTEGER, INTENT(IN OUT)                  :: jary
INTEGER*4 iary(*),jary(*)

l = n / 2  +  1
ir = n

DO WHILE(ir > 1)
  IF(l > 1) THEN
    l = l - 1
    ia = iary(l)
    ja = jary(l)
  ELSE
    ia = iary(ir)
    ja = jary(ir)
    iary(ir) = iary(1)
    jary(ir) = jary(1)
    ir = ir - 1
    IF(ir == 1) THEN
      iary(1) = ia
      jary(1) = ja
      RETURN
    END IF
  END IF
  i = l
  j = l + l
  DO WHILE(j <= ir)
    IF(j < ir) THEN
      IF(iary(j) < iary(j+1)) j = j + 1
    END IF
    IF(ia < iary(j)) THEN
      iary(i) = iary(j)
      jary(i) = jary(j)
      i = j
      j = j + j
    ELSE
      j = ir + 1
    END IF
  END DO
  iary(i) = ia
  jary(i) = ja
END DO

RETURN
END SUBROUTINE isort2

!-----------------------------------------------------------------------
!  Sorting routine (Quick Sort) for TWO arrays, (D.P. & integer)

SUBROUTINE di_sort2(n,dary,iary)


INTEGER, INTENT(IN)                      :: n
REAL, INTENT(IN OUT)                     :: dary
INTEGER, INTENT(IN OUT)                  :: iary
REAL*8 dary(*),da
INTEGER*4 iary(*)

l = n / 2  +  1
ir = n

DO WHILE(ir > 1)
  IF(l > 1) THEN
    l = l - 1
    da = dary(l)
    ia = iary(l)
  ELSE
    da = dary(ir)
    ia = iary(ir)
    dary(ir) = dary(1)
    iary(ir) = iary(1)
    ir = ir - 1
    IF(ir == 1) THEN
      dary(1) = da
      iary(1) = ia
      RETURN
    END IF
  END IF
  i = l
  j = l + l
  DO WHILE(j <= ir)
    IF(j < ir) THEN
      IF(dary(j) < dary(j+1)) j = j + 1
    END IF
    IF(da < dary(j)) THEN
      dary(i) = dary(j)
      iary(i) = iary(j)
      i = j
      j = j + j
    ELSE
      j = ir + 1
    END IF
  END DO
  dary(i) = da
  iary(i) = ia
END DO

RETURN
END SUBROUTINE di_sort2
