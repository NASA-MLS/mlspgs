Module GET_BETA_PATH_M
  use MLSCommon, only: R8, RP, IP
  use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
  use SpectroscopyCatalog_m, only: CATALOG_T, LINES
  use CREATE_BETA_M, only: CREATE_BETA
  Implicit NONE
  private
  PUBLIC :: get_beta_path

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
  "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
! This is a generic form of get coarse beta path. We really don't need
! separate versions of these.
!
 SUBROUTINE get_beta_path(frq,p_path,t_path,Catalog,gl_slabs,path_inds,      &
        &   beta_path,gl_slabs_m,t_path_m,gl_slabs_p,t_path_p,dbeta_dt_path, &
        &   dbeta_dw_path,dbeta_dn_path,dbeta_dv_path)

!  ===============================================================
!  Declaration of variables for sub-program: get_coarse_beta_path
!  ===============================================================
! Inputs:
!
  REAL(r8), INTENT(in) :: frq ! frequency in MHz
  REAL(rp), INTENT(in) :: t_path(:) ! path temperatures
  REAL(rp), INTENT(in) :: p_path(:) ! path pressures in hPa!
  Type(Catalog_T), INTENT(IN) :: Catalog(:)
  Type (slabs_struct), DIMENSION(:,:), POINTER :: gl_slabs
  INTEGER(ip), INTENT(in) :: path_inds(:) ! indicies for reading gl_slabs
!
! Another clumsy feature of f90
! Optional inputs:
!
  TYPE(slabs_struct), OPTIONAL, POINTER :: gl_slabs_m(:,:) ! reduced
!                               strength data for t_path_m
  REAL(rp), OPTIONAL, INTENT(in) :: t_path_m(:) ! path temperatures for
!                                                 gl_slabs_m
  TYPE(slabs_struct), OPTIONAL, POINTER :: gl_slabs_p(:,:) ! reduced
!                               strength data for t_path_p
  REAL(rp), OPTIONAL, INTENT(in) :: t_path_p(:) ! path temperatures for
!
! outputs
!
  REAL(rp), INTENT(out) :: beta_path(:,:) ! path beta for each species
!
! Optional output:
!
  REAL(rp), OPTIONAL, INTENT(out) :: dbeta_dt_path(:,:) ! t dep.
  REAL(rp), OPTIONAL, INTENT(out) :: dbeta_dw_path(:,:) ! line width
  REAL(rp), OPTIONAL, INTENT(out) :: dbeta_dn_path(:,:) ! line width t dep.
  REAL(rp), OPTIONAL, INTENT(out) :: dbeta_dv_path(:,:) ! line position
!
! Local varibles..
!
  Integer(ip) :: n_sps, n_path, i, j, k, m, nl, no_of_lines, Spectag
  REAL(rp) :: bb, vp, v0, vm, vn2, t, tm, tp, bp, bm
  REAL(rp), allocatable, dimension(:) :: LineWidth
!
! begin the code
!
  n_sps = Size(Catalog)
  n_path = SIZE(path_inds)
!
! no derivative call
!
  DO i = 1, n_sps
    Spectag = Catalog(i)%Spec_Tag
    no_of_lines = gl_slabs(1,i)%no_lines
    Allocate(LineWidth(no_of_lines))
    do k = 1, no_of_lines
      m = Catalog(i)%Lines(k)
      LineWidth(k) = Lines(m)%W
    end do
    DO j = 1, n_path
      k = path_inds(j)
      CALL create_beta(Spectag,Catalog(i)%continuum,p_path(k),t_path(k), &
        &  Frq,no_of_lines,LineWidth,gl_slabs(k,i)%v0s,gl_slabs(k,i)%x1, &
        &  gl_slabs(k,i)%y,gl_slabs(k,i)%yi,gl_slabs(k,i)%slabs1,bb,     &
        &  gl_slabs(k,i)%dslabs1_dv0,DBETA_DW=v0,DBETA_DN=vp,DBETA_DV=vm)
      beta_path(j,i) = bb
      if(PRESENT(dbeta_dw_path)) dbeta_dw_path(j,i) = v0
      if(PRESENT(dbeta_dn_path)) dbeta_dn_path(j,i) = vp
      if(PRESENT(dbeta_dv_path)) dbeta_dv_path(j,i) = vm
    END DO
    DEAllocate(LineWidth)
  ENDDO
!
  IF(PRESENT(dbeta_dt_path)) THEN
!
    DO i = 1 , n_sps
      Spectag = Catalog(i)%Spec_Tag
      no_of_lines = gl_slabs_m(1,i)%no_lines
      Allocate(LineWidth(no_of_lines))
      do k = 1, no_of_lines
        m = Catalog(i)%Lines(k)
        LineWidth(k) = Lines(m)%W
      end do
      DO j = 1 , n_path
        k = path_inds(j)
        tm = t_path_m(k)
        CALL create_beta(Spectag,Catalog(i)%continuum,p_path(k),tm,Frq,   &
        &    no_of_lines,LineWidth,gl_slabs_m(k,i)%v0s,gl_slabs_m(k,i)%x1,&
        &    gl_slabs_m(k,i)%y,gl_slabs_m(k,i)%yi,gl_slabs_m(k,i)%slabs1, &
        &    bm,gl_slabs_m(k,i)%dslabs1_dv0)
        tp = t_path_p(k)
        CALL create_beta(Spectag,Catalog(i)%continuum,p_path(k),tp,Frq,   &
        &    no_of_lines,LineWidth,gl_slabs_m(k,i)%v0s,gl_slabs_p(k,i)%x1,&
        &    gl_slabs_p(k,i)%y,gl_slabs_p(k,i)%yi,gl_slabs_p(k,i)%slabs1, &
        &    bp,gl_slabs_p(k,i)%dslabs1_dv0)
        t  = t_path(k)
        bb = beta_path(j,i)
        if ( bp > 0.0 .and. bb > 0.0 .and. bm > 0.0 ) then
          vp = Log(bp/bb)/Log(tp/t)        ! Estimate over [temp+10,temp]
          v0 = Log(bp/bm)/Log(tp/tm)       ! Estimate over [temp+10,temp-10]
          vm = Log(bb/bm)/Log(t/tm)        ! Estimate over [temp,temp-10]
        else if ( bp > 0.0 .and. bb > 0.0) then
          vp = Log(bp/bb)/Log(tp/t)        ! Estimate over [temp+10,temp]
          vm = vp
          v0 = vp
        else if ( bm > 0.0 .and. bb > 0.0) then
          vm = Log(bb/bm)/Log(t/tm)        ! Estimate over [temp,temp-10]
          vp = vm
          v0 = vm
        else if ( bm > 0.0 .and. bp > 0.0) then
          v0 = Log(bp/bm)/Log(tp/tm)       ! Estimate over [temp+10,temp-10]
          vp = v0
          vm = v0
        else
          vp = 0.0
          v0 = 0.0
          vm = 0.0
        endif
        dbeta_dt_path(j,i) = (vp + 2.0 * v0 + vm) / 4.0  ! Weighted Average
      ENDDO
      DeAllocate(LineWidth)
    ENDDO
  ENDIF
!
 END SUBROUTINE get_beta_path

!----------------------------------------------------------------------
End module GET_BETA_PATH_M
! $Log$
! Revision 2.4  2001/11/10 00:46:40  zvi
! Adding the EXTINCTION capabilitis
!
! Revision 2.3  2001/11/07 22:24:45  zvi
! Further modification for the t-power computations
!
! Revision 2.2  2001/11/07 21:13:48  livesey
! Fixed bug with log(0.0/0.0) for molecues with no lines
! or continua
!
! Revision 2.1  2001/10/16 15:07:18  zvi
! Continuum parameters are now part of Catalog
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.22.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.8  2001/03/09 02:26:11  vsnyder
! More work on deallocation
!
! Revision 1.7  2001/03/09 02:11:28  vsnyder
! Repair deallocating
!
! Revision 1.6  2001/03/05 21:37:20  zvi
! New filter format
!
! Revision 1.1 2001/02/01 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
