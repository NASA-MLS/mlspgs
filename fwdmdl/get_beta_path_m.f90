! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GET_BETA_PATH_M

  use MLSCommon, only: RP, R8

  implicit NONE
  private
  public :: Get_Beta_Path, Get_Beta_Path_Scalar, Get_Beta_Path_Polarized
  public :: Get_Beta_Path_Cloud

  interface Get_Beta_Path
    module procedure Get_Beta_Path_Scalar, Get_Beta_Path_Polarized
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

! ---------------------------------------  Get_Beta_Path_Scalar  -----
  subroutine Get_Beta_Path_Scalar ( frq, p_path, t_path, &
        & Catalog, beta_group, polarized, gl_slabs, path_inds,     &
        & beta_path, gl_slabs_m, t_path_m, gl_slabs_p, t_path_p,   &
        & dbeta_dt_path, dbeta_dw_path, dbeta_dn_path, dbeta_dv_path)

    use CREATE_BETA_M, only: CREATE_BETA
    use Get_Species_Data_m, only: Beta_Group_T
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: R8, RP, IP
    use Molecules, only: SP_H2O
    use Physics, only: H_OVER_K
    use SpectroscopyCatalog_m, only: CATALOG_T, LINES

! Inputs:

    real(r8), intent(in) :: Frq         ! frequency in MHz
    real(rp), intent(in) :: T_path(:)   ! path temperatures
    real(rp), intent(in) :: P_path(:)   ! path pressures in hPa!
    type(catalog_t), intent(in) :: Catalog(:)
    type (slabs_struct), dimension(:,:) :: Gl_slabs
    integer(ip), intent(in) :: Path_inds(:) ! indicies for reading gl_slabs

    type (beta_group_T), dimension(:) :: beta_group

    logical, intent(in) :: Polarized    ! "Don't work on Zeeman-split lines"

! Optional inputs.  GL_SLABS_* are pointers because the caller need not
! allocate them if DBETA_D*_PATH aren't allocated.  They would be
! INTENT(IN) if we could say so.

    type(slabs_struct), pointer :: gl_slabs_m(:,:) ! reduced
!                               strength data for t_path_m
    real(rp), intent(in) :: t_path_m(:) ! path temperatures for gl_slabs_m
    type(slabs_struct), pointer :: gl_slabs_p(:,:) ! reduced
!                               strength data for t_path_p
    real(rp), intent(in) :: t_path_p(:) ! path temperatures for gl_slabs_p

! outputs

    real(rp), intent(out) :: beta_path(:,:) ! path beta for each specie

! Optional outputs.  We use ASSOCIATED instead of PRESENT so that the
! caller doesn't need multiple branches.  These would be INTENT(OUT) if
! we could say so.

    real(rp), pointer :: dbeta_dt_path(:,:) ! t dep.
    real(rp), pointer :: dbeta_dw_path(:,:) ! line width
    real(rp), pointer :: dbeta_dn_path(:,:) ! line width t dep.
    real(rp), pointer :: dbeta_dv_path(:,:) ! line position

! Local variables..

    integer(ip) :: I, J, K, N, IB, Spectag, No_of_lines, &
              &    No_mol, N_path
    real(rp) :: Ratio, BB, VP, V0, VM, T, TM, TP, BP, BM
    real(rp), allocatable, dimension(:) :: LineWidth
    real(rp), dimension(size(path_inds)) :: betam, betap
    real(rp), dimension(size(path_inds)) :: tanh1, tanh1_p, tanh1_m

! begin the code

    no_mol = size(beta_group)
    n_path = size(path_inds)

    beta_path = 0.0
! Note that expa only depends on temperature.
! compute path hyperbolic tangent
    tanh1 = tanh(0.5_rp * h_over_k * frq / t_path(path_inds))

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

          call create_beta ( spectag, catalog(ib)%continuum, p_path(k),   &
            & t_path(k), Frq, lines(catalog(ib)%lines)%w, gl_slabs(k,ib), &
            & tanh1(j), bb, polarized .and. catalog(ib)%polarized,        &
            & DBETA_DW=v0, DBETA_DN=vp, DBETA_DV=vm )

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
      tanh1_p = tanh(0.5_rp * h_over_k * frq / t_path_p(path_inds))
      tanh1_m = tanh(0.5_rp * h_over_k * frq / t_path_m(path_inds))

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
            call create_beta ( spectag, catalog(ib)%continuum, p_path(k), &
            &  tm, frq, linewidth, gl_slabs_m(k,ib),                      &
            &  tanh1_m(j), vm, polarized .and. catalog(ib)%polarized )
            betam(j) = betam(j) + ratio * vm
            tp = t_path_p(k)
            call create_beta ( Spectag, Catalog(ib)%continuum, p_path(k), &
            &  tp, Frq, LineWidth, gl_slabs_p(k,ib),                      &
            &  tanh1_p(j), vp, polarized .and. catalog(ib)%polarized )
            betap(j) = betap(j) + ratio * vp
          end do
          deallocate ( LineWidth )
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

  end subroutine Get_Beta_Path_Scalar

  ! ------------------------------------  Get_Beta_Path_Polarized  -----
  subroutine Get_Beta_Path_Polarized ( Frq, H, &
        & Catalog, Beta_group, GL_slabs, Path_inds, Beta_path )

    use Get_Species_Data_m, only: Beta_Group_T
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: R8, RP, IP
    use O2_Abs_CS_m, only: O2_Abs_CS
    use SpectroscopyCatalog_m, only: CATALOG_T

! Inputs:

    real(r8), intent(in) :: Frq ! frequency in MHz
    real(rp), intent(in) :: H(:)      ! Magnetic field component in instrument
                                      ! polarization on the path
    type(catalog_t), intent(in) :: Catalog(:)
    type (slabs_struct), dimension(:,:), intent(in) :: GL_slabs
    integer(ip), intent(in) :: Path_inds(:) ! indicies for reading gl_slabs
    type (beta_group_T), dimension(:), intent(in) :: Beta_group

! outputs

    complex(rp), intent(out) :: beta_path(-1:,:,:) ! path beta for each species
    ! beta_path(-1,:,:) is Sigma_m, beta_path(0,:,:) is Pi,
    ! beta_path(+1,:,:) is Sigma_p

! Local variables..

    integer(ip) :: I, IB, J, K, M, N, N_PATH
    real(rp) :: RATIO ! Isotope ratio, not mixing ratio
    complex(rp) :: Sigma_m, Pi, Sigma_p

! begin the code

    n_path = size(path_inds)

    beta_path = 0.0

    do i = 1, size(beta_group)
      do n = 1, beta_group(i)%n_elements
        ratio = beta_group(i)%ratio(n)
        ib = beta_group(i)%cat_index(n)

        do j = 1, n_path
          k = path_inds(j)

          call o2_abs_cs ( frq, (/ ( -1, m=1,size(catalog(ib)%polarized) ) /), &
            & h(k), gl_slabs(k,ib), catalog(ib)%polarized,        &
            & catalog(ib)%continuum(1), catalog(ib)%continuum(3), &
            & sigma_p, pi, sigma_m )
          beta_path(-1,j,i) = beta_path(-1,j,i) + ratio * sigma_m
          beta_path( 0,j,i) = beta_path( 0,j,i) + ratio * pi
          beta_path(+1,j,i) = beta_path(+1,j,i) + ratio * sigma_p
        end do ! j
      end do ! n
    end do ! i
  end subroutine Get_Beta_Path_Polarized

  ! ----------------------------------------  Get_Beta_Path_Cloud  -----
  subroutine Get_Beta_Path_Cloud ( Frq, p_path, t_path,                 &
        & beta_group, path_inds, beta_path_cloud,                       &
        & beta_path_w0,  beta_path_phh,                                 & 
        & IPSD, WC, fwdModelConf  )
    use ForwardModelConfig, only: FORWARDMODELCONFIG_T
    use Cloud_extinction, only: get_beta_cloud
    use Get_Species_Data_m, only: Beta_Group_T
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: R8, RP, IP

! Inputs:

    real(r8), intent(in) :: Frq ! frequency in MHz
    real(rp), intent(in) :: T_path(:)   ! path temperatures
    real(rp), intent(in) :: P_path(:)   ! path pressures in hPa!
    
    integer(ip), intent(in) :: Path_inds(:) ! indicies for reading gl_slabs

    type (beta_group_T), dimension(:) :: beta_group
    type (ForwardModelConfig_T) ,intent(in) :: FWDMODELCONF

    INTEGER, intent(in) :: IPSD(:)
    REAL(rp), intent(in)  :: WC(:,:)
    REAL(rp) :: W0       ! SINGLE SCATTERING ALBEDO
    REAL(rp) :: PHH(fwdModelConf%num_scattering_angles)   ! PHASE FUNCTION

! outputs

    real(rp), intent(out) :: beta_path_cloud(:) ! cloud extinction
    real(rp), intent(out) :: beta_path_w0(:)    ! single scattering albedo
    real(rp), intent(out) :: beta_path_phh(:,:)   ! phase function

! Optional outputs.  We use ASSOCIATED instead of PRESENT so that the
! caller doesn't need multiple branches.  These would be INTENT(OUT) if
! we could say so.

! Local variables..

    INTEGER :: NC, NU, NUA, NAB, NR, N
    logical :: Incl_Cld
    integer(ip) :: i, j, k, n_path
    real(rp) :: cld_ext, RHI

! begin the code

          Incl_Cld  = fwdModelConf%Incl_Cld
          NU  = fwdModelConf%NUM_SCATTERING_ANGLES
          NUA = fwdModelConf%NUM_AZIMUTH_ANGLES
          NAB = fwdModelConf%NUM_AB_TERMS
          NR  = fwdModelConf%NUM_SIZE_BINS
          N   = fwdModelConf%no_cloud_species

    n_path = size(path_inds)

    beta_path_cloud = 0.0
    beta_path_w0    = 0.0
    beta_path_phh   = 0.0

        do j = 1, n_path
          k = path_inds(j)

          call get_beta_cloud (Frq, t_path(k),                       &
                          &  WC(:,k), IPSD(k), NC, NU, NUA, NAB, NR, &
                          &  cld_ext, W0, PHH                )      

            beta_path_cloud(j) = beta_path_cloud(j) + cld_ext 
            beta_path_w0(j) = beta_path_w0(j) + W0 
            beta_path_phh(j,:) = beta_path_phh(j,:) + PHH(:) 
            
         end do

  end subroutine Get_Beta_Path_Cloud

!-----------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module GET_BETA_PATH_M

! $Log$
! Revision 2.30  2003/05/09 19:24:38  vsnyder
! Cosmetic change
!
! Revision 2.29  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.25.2.9  2003/03/24 21:50:44  jonathan
! remove unused 'pressure' in call to cloud_extinction
!
! Revision 2.25.2.8  2003/03/22 04:03:04  vsnyder
! Move Beta_Group_T and Dump_Beta_Group from get_beta_path to Get_Species_Data
!
! Revision 2.25.2.7  2003/03/12 21:35:44  vsnyder
! Add Dump_Beta_Group and generic Dump for it
!
! Revision 2.25.2.6  2003/03/01 03:16:15  vsnyder
! Use 'polarized' to specify size of quantum numbers array
!
! Revision 2.25.2.5  2003/02/27 23:19:23  vsnyder
! Add polarized stuff.  Remove unused z_path_c argument of get_beta_path.
! Remove unused p_path argument to get_beta_cloud.  Cosmetics.
!
! Revision 2.25.2.3  2003/02/14 00:21:42  jonathan
! add singl. scat. albedo W0, ph funct PHH
!
! Revision 2.25.2.2  2003/02/13 22:26:30  jonathan
! changes dimension for beta_path_cloud
!
! Revision 2.25.2.1  2003/02/13 17:34:27  bill
! uses new slabs_sw
!
! Revision 2.25  2003/02/11 00:48:18  jonathan
! changes made after adding get_beta_path_cloud
!
! Revision 2.24  2003/02/07 01:57:19  vsnyder
! Delete USE RHIFromH2O because it's not used
!
! Revision 2.23  2003/02/07 01:08:34  jonathan
! remove ICON option for compute super saturation
!
! Revision 2.22  2003/02/06 22:12:49  jonathan
! fix bug
!
! Revision 2.21  2003/02/06 00:20:16  jonathan
! Add in many stuff to deal with clouds CloudIce, iwc_path, etc
!
! Revision 2.20  2003/02/04 22:03:33  jonathan
! ICON now equal to 0 as default
!
! Revision 2.19  2003/02/04 21:46:27  jonathan
! add ICON options for super saturation and dry cases
!
! Revision 2.18  2003/02/03 22:56:58  vsnyder
! Add Get_bata_path_polarized
!
! Revision 2.17  2003/01/31 18:45:09  jonathan
! use cld_ext only if Incl_Cld is ture
!
! Revision 2.16  2003/01/31 17:53:48  jonathan
! change z_path to z_path_c in passing to get_beta_path
!
! Revision 2.15  2003/01/31 17:16:08  jonathan
! add Inc_Cld, and cld_ext
!
! Revision 2.14  2003/01/30 17:43:04  jonathan
! remove RHtoEV
!
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
