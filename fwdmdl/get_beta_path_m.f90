Module GET_BETA_PATH_M
  use MLSCommon, only: R8, RP, IP
  use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
  use SpectroscopyCatalog_m, only: CATALOG_T, LINES
  use CREATE_BETA_M, only: CREATE_BETA
  Implicit NONE
  private
  PUBLIC :: get_beta_path

! *** Beta group type declaration:
  type, public :: beta_group_T
    integer :: n_elements
    integer, pointer  :: cat_index(:)
    real(rp), pointer :: ratio(:)
  end type beta_group_T

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
! This is a generic form of get coarse beta path. We really don't need
! separate versions of these.
!
 SUBROUTINE get_beta_path(frq,p_path,t_path,Catalog,beta_group,gl_slabs, &
        &   path_inds,beta_path,gl_slabs_m,t_path_m,gl_slabs_p,t_path_p, &
        &   dbeta_dt_path,dbeta_dw_path,dbeta_dn_path,dbeta_dv_path)

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

  type (beta_group_T), dimension(:), pointer :: beta_group
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
  Integer(ip) :: i, j, k, m, n, nl, ib, nbe, Spectag, no_of_lines, &
            &    no_mol, n_path
  REAL(rp) :: ratio, bb, vp, v0, vm, t, tm, tp, bp, bm
  REAL(rp), allocatable, dimension(:) :: LineWidth
!
! begin the code
!
  no_mol = Size(beta_group)
  n_path = SIZE(path_inds)
!
  beta_path = 0.0
  if(PRESENT(dbeta_dw_path)) dbeta_dw_path = 0.0
  if(PRESENT(dbeta_dn_path)) dbeta_dn_path = 0.0
  if(PRESENT(dbeta_dv_path)) dbeta_dv_path = 0.0
!
  DO i = 1, no_mol
    nbe = beta_group(i)%n_elements
    do n = 1, nbe
      ratio = beta_group(i)%ratio(n)
      ib = beta_group(i)%cat_index(n)
      Spectag = Catalog(ib)%Spec_Tag
      no_of_lines = gl_slabs(1,ib)%no_lines
      Allocate(LineWidth(no_of_lines))
      do k = 1, no_of_lines
        m = Catalog(ib)%Lines(k)
        LineWidth(k) = Lines(m)%W
      end do
      DO j = 1, n_path
        k = path_inds(j)
        CALL create_beta(Spectag,Catalog(ib)%continuum,p_path(k),t_path(k), &
          &  Frq,no_of_lines,LineWidth,gl_slabs(k,ib)%v0s,gl_slabs(k,ib)%x1,&
          &  gl_slabs(k,ib)%y,gl_slabs(k,ib)%yi,gl_slabs(k,ib)%slabs1,bb,   &
          &  gl_slabs(k,ib)%dslabs1_dv0,DBETA_DW=v0,DBETA_DN=vp,DBETA_DV=vm)
        beta_path(j,i) = beta_path(j,i) + ratio * bb
        if(PRESENT(dbeta_dw_path)) &
          &  dbeta_dw_path(j,i) = dbeta_dw_path(j,i) + ratio * v0
        if(PRESENT(dbeta_dn_path)) &
          &  dbeta_dn_path(j,i) = dbeta_dn_path(j,i) + ratio * vp
        if(PRESENT(dbeta_dv_path)) &
          &  dbeta_dv_path(j,i) = dbeta_dv_path(j,i) + ratio * vm
      END DO
      DEAllocate(LineWidth)
    end do
  END DO
!
  IF(PRESENT(dbeta_dt_path)) THEN

    dbeta_dt_path = 0.0
!
    DO i = 1, no_mol
      bm = 0.0_rp
      bp = 0.0_rp
      nbe = beta_group(i)%n_elements
      do n = 1, nbe
        ratio = beta_group(i)%ratio(n)
        ib = beta_group(i)%cat_index(n)
        Spectag = Catalog(ib)%Spec_Tag
        no_of_lines = gl_slabs_m(1,ib)%no_lines
        Allocate(LineWidth(no_of_lines))
        do k = 1, no_of_lines
          m = Catalog(ib)%Lines(k)
          LineWidth(k) = Lines(m)%W
        end do
        DO j = 1 , n_path
          k = path_inds(j)
          tm = t_path_m(k)
          CALL create_beta(Spectag,Catalog(ib)%continuum,p_path(k),tm,Frq,    &
          &    no_of_lines,LineWidth,gl_slabs_m(k,ib)%v0s,gl_slabs_m(k,ib)%x1,&
          &    gl_slabs_m(k,ib)%y,gl_slabs_m(k,ib)%yi,gl_slabs_m(k,ib)%slabs1,&
          &    vm,gl_slabs_m(k,ib)%dslabs1_dv0)
          bm = bm + vm * ratio
          tp = t_path_p(k)
          CALL create_beta(Spectag,Catalog(ib)%continuum,p_path(k),tp,Frq,    &
          &    no_of_lines,LineWidth,gl_slabs_m(k,ib)%v0s,gl_slabs_p(k,ib)%x1,&
          &    gl_slabs_p(k,ib)%y,gl_slabs_p(k,ib)%yi,gl_slabs_p(k,ib)%slabs1,&
          &    vp,gl_slabs_p(k,ib)%dslabs1_dv0)
          bp = bp + vp * ratio
        END DO
        DeAllocate(LineWidth)
      end do
      DO j = 1 , n_path
        k = path_inds(j)
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
      END DO
    ENDDO
  ENDIF
!
 END SUBROUTINE get_beta_path

!----------------------------------------------------------------------
End module GET_BETA_PATH_M
! $Log$
! Revision 2.5  2001/11/15 01:22:01  zvi
! Remove Extiction debug
!
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
