! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GET_BETA_PATH_M

  use MLSCommon, only: RP, R8
  use RHIFromH2O, only: RHIFromH2O_Factor
  implicit NONE
  private
  public :: get_beta_path

! *** Beta group type declaration:
  type, public :: beta_group_T
    integer :: n_elements
    integer, pointer  :: cat_index(:)
    real(rp), pointer :: ratio(:)
  end type beta_group_T

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
! This is a generic form of get coarse beta path. We really don't need
! separate versions of these.

  ! ----------------------------------------------  Get_Beta_Path  -----
  subroutine Get_Beta_Path ( frq, p_path, t_path, z_path, Catalog, beta_group, gl_slabs, &
        & path_inds, beta_path, gl_slabs_m, t_path_m, gl_slabs_p, t_path_p, &
        & dbeta_dt_path, dbeta_dw_path, dbeta_dn_path, dbeta_dv_path, ICON )

    use MLSCommon, only: R8, RP, IP
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use SpectroscopyCatalog_m, only: CATALOG_T, LINES
    use CREATE_BETA_M, only: CREATE_BETA
    use Molecules, only: SP_H2O

! Inputs:

    real(r8), intent(in) :: Frq ! frequency in MHz
    real(rp), intent(in) :: T_path(:) ! path temperatures
    real(rp), intent(in) :: P_path(:) ! path pressures in hPa!
    real(rp), intent(in) :: Z_path(:) ! =-log(p_path)
    type(catalog_t), intent(in) :: Catalog(:)
    type (slabs_struct), dimension(:,:) :: Gl_slabs
    integer(ip), intent(in) :: Path_inds(:) ! indicies for reading gl_slabs

    type (beta_group_T), dimension(:) :: beta_group

! Optional inputs.  GL_SLABS_* are pointers because the caller need not
! allocate them if DBETA_D*_PATH aren't allocated.  They would be
! INTENT(IN) if we could say so.

    type(slabs_struct), pointer :: gl_slabs_m(:,:) ! reduced
!                               strength data for t_path_m
    real(rp) :: t_path_m(:) ! path temperatures for gl_slabs_m
    type(slabs_struct), pointer :: gl_slabs_p(:,:) ! reduced
!                               strength data for t_path_p
    real(rp) :: t_path_p(:) ! path temperatures for gl_slabs_p

! outputs

    real(rp), intent(out) :: beta_path(:,:) ! path beta for each species

! Optional outputs.  We use ASSOCIATED instead of PRESENT so that the
! caller doesn't need multiple branches.  These would be INTENT(OUT) if
! we could say so.

    real(rp), pointer :: dbeta_dt_path(:,:) ! t dep.
    real(rp), pointer :: dbeta_dw_path(:,:) ! line width
    real(rp), pointer :: dbeta_dn_path(:,:) ! line width t dep.
    real(rp), pointer :: dbeta_dv_path(:,:) ! line position

! Local variables..

    integer(ip) :: i, j, k, n, ib, Spectag, no_of_lines, &
              &    no_mol, n_path
    real(rp) :: ratio, bb, vp, v0, vm, t, tm, tp, bp, bm
    real(rp), allocatable, dimension(:) :: LineWidth
    real(rp), dimension(size(path_inds)) :: betam, betap

    integer  :: ICON
    real(rp) :: P, Vapor_P

! begin the code

    no_mol = size(beta_group)
    n_path = size(path_inds)

    beta_path = 0.0
    if ( associated(dbeta_dw_path) ) dbeta_dw_path(1:n_path,:) = 0.0
    if ( associated(dbeta_dn_path) ) dbeta_dn_path(1:n_path,:) = 0.0
    if ( associated(dbeta_dv_path) ) dbeta_dv_path(1:n_path,:) = 0.0

    do i = 1, no_mol
      do n = 1, beta_group(i)%n_elements
        ratio = beta_group(i)%ratio(n)
        ib = beta_group(i)%cat_index(n)
        Spectag = Catalog(ib)%Spec_Tag
         
        do j = 1, n_path
          k = path_inds(j)

          ! mask 100%RH below 100mb
          IF(Spectag .EQ. SP_H2O .AND. ICON .EQ.-1 .and. p_path(k).GE. 100.)THEN
            ratio = RHIFromH2O_Factor (t_path(k), z_path(k), 0, .true.)*100._r8
            ! optional 0 will return ratio as parts per 1, as Bill uses here.
          ENDIF                                 

          call create_beta ( Spectag, Catalog(ib)%continuum, p_path(k), t_path(k), &
            &  Frq, Lines(Catalog(ib)%Lines)%W, gl_slabs(k,ib), bb,     &
            &  DBETA_DW=v0, DBETA_DN=vp, DBETA_DV=vm )
          beta_path(j,i) = beta_path(j,i) + ratio * bb
          if ( associated(dbeta_dw_path)) &
            &  dbeta_dw_path(j,i) = dbeta_dw_path(j,i) + ratio * v0
          if ( associated(dbeta_dn_path)) &
            &  dbeta_dn_path(j,i) = dbeta_dn_path(j,i) + ratio * vp
          if ( associated(dbeta_dv_path)) &
            &  dbeta_dv_path(j,i) = dbeta_dv_path(j,i) + ratio * vm
        end do
      end do
    end do

    if ( associated(dbeta_dt_path) ) then

      dbeta_dt_path(1:n_path,:) = 0.0

      do i = 1, no_mol
        betam = 0.0
        betap = 0.0
        do n = 1, beta_group(i)%n_elements
          ratio = beta_group(i)%ratio(n)
          ib = beta_group(i)%cat_index(n)
          Spectag = Catalog(ib)%Spec_Tag
          no_of_lines = gl_slabs_m(1,ib)%no_lines
          allocate ( LineWidth(no_of_lines) )
          do k = 1, no_of_lines
            LineWidth(k) = Lines(Catalog(ib)%Lines(k))%W
          end do
          do j = 1 , n_path
            k = path_inds(j)
            tm = t_path_m(k)
            call create_beta ( Spectag, Catalog(ib)%continuum, p_path(k), tm, Frq, &
            &    LineWidth, gl_slabs_m(k,ib), vm )
            betam(j) = betam(j) + ratio * vm
            tp = t_path_p(k)
            call create_beta ( Spectag, Catalog(ib)%continuum, p_path(k), tp, Frq, &
            &    LineWidth, gl_slabs_p(k,ib), vp )
            betap(j) = betap(j) + ratio * vp
          end do
          DeAllocate ( LineWidth )
        end do
        do j = 1 , n_path
          k = path_inds(j)
          t  = t_path(k)
          tm = t_path_m(k)
          tp = t_path_p(k)
          bm = betam(j)
          bp = betap(j)
          bb = beta_path(j,i)
          if ( bp > 0.0 .and. bb > 0.0 .and. bm > 0.0 ) then
            vp = Log(bp/bb)/Log(tp/t)        ! Estimate over [temp+10,temp]
            v0 = Log(bp/bm)/Log(tp/tm)       ! Estimate over [temp+10,temp-10]
            vm = Log(bb/bm)/Log(t/tm)        ! Estimate over [temp,temp-10]
          else if ( bp > 0.0 .and. bb > 0.0 ) then
            vp = Log(bp/bb)/Log(tp/t)        ! Estimate over [temp+10,temp]
            vm = vp
            v0 = vp
          else if ( bm > 0.0 .and. bb > 0.0 ) then
            vm = Log(bb/bm)/Log(t/tm)        ! Estimate over [temp,temp-10]
            vp = vm
            v0 = vm
          else if ( bm > 0.0 .and. bp > 0.0 ) then
            v0 = Log(bp/bm)/Log(tp/tm)       ! Estimate over [temp+10,temp-10]
            vp = v0
            vm = v0
          else
            vp = 0.0
            v0 = 0.0
            vm = 0.0
          end if
          dbeta_dt_path(j,i) = (vp + 2.0 * v0 + vm) / 4.0  ! Weighted Average
        end do
      end do
    end if

  end subroutine Get_Beta_Path

!----------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module GET_BETA_PATH_M

! $Log$
! Revision 2.13  2003/01/30 00:17:42  jonathan
! add z_path to get_beta_path & use Paul's RHIFromH2O to compute VMR from RHi
!
! Revision 2.12  2003/01/14 21:49:33  jonathan
! option for saturation below 100mb
!
! Revision 2.11  2003/01/08 00:17:29  vsnyder
! Use "associated" instead of "present" to control optional computations
!
! Revision 2.10  2002/12/13 02:06:51  vsnyder
! Use a SLABS structure for the slabs quantities
!
! Revision 2.9  2002/10/08 17:08:03  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.8  2002/09/12 23:00:04  vsnyder
! Cosmetic changes, move USEs from module scope to procedure scope
!
! Revision 2.7  2001/12/23 23:30:42  zvi
! Fixing a bug in the dbeta_dt computations
!
! Revision 2.6  2001/12/14 23:43:15  zvi
! Modification for Grouping concept
!
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
