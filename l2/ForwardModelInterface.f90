! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelInterface
!=============================================================================

! Set up the forward model.  Interface from the retrieve step to the
! forward model.

!??? Do we want a forward model database ???

  use Expr_M, only: EXPR
  use Init_Tables_Module, only: field_first, field_last
  use Init_Tables_Module, only: F_ATMOS_DER, F_DO_CONV, F_DO_FREQ_AVG, &
    & F_FREQUENCY, F_SPECT_DER, F_TEMP_DER
  use Lexer_Core, only: Print_Source
  use MatrixModule_1, only: Matrix_Database_T, Matrix_T
  use MLSCommon, only: R8
  use MoreTree, only: Get_Boolean, Get_Field_ID
  use Output_M, only: Output
  use String_Table, only: Display_String
  use Toggles, only: Gen, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Node_ID, Nsons, Source_Ref, Subtree
  use Tree_Types, only: N_named
  use VectorsModule, only: Vector_T

  !??? The next USE statement is Temporary for l2load:
  use L2_TEST_STRUCTURES_M, only: FWD_MDL_CONFIG, FWD_MDL_INFO, &
    & TEMPORARY_FWD_MDL_INFO

  implicit NONE
  private
  public :: ForwardModel, ForwardModelGlobalSetup, ForwardModelInfo_T, &
    & ForwardModelSetup

  type ForwardModelInfo_T
    logical :: Atmos_Der      ! Do atmospheric derivatives
    logical :: Do_Conv        ! Do convolution
    logical :: Do_Freq_Avg    ! Do Frequency averaging
    real(r8) :: The_Freq      ! Frequency to use if .not. do_freq_avg
    logical :: Spect_Der      ! Do spectroscopy derivatives
    logical :: Temp_Der       ! Do temperature derivatives
  end type ForwardModelInfo_T

  integer :: Error            ! Error level -- 0 = OK

  !---------------------------- RCS Ident Info -------------------------------
  character (len=130), private :: Id = &
    & "$Id$"
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  ! ------------------------------------  ForwardModelGlobalSetup  -----
  subroutine ForwardModelGlobalSetup ( Root, ForwardModelInfo )
  ! Process the forwardModel specification to produce ForwardModelInfo.

    integer :: Root                     ! of the forwardModel specification.
                                        ! Indexes either a "named" or
                                        ! "spec_args" vertex.
    type(forwardModelInfo_T), intent(inout) :: ForwardModelInfo

    integer :: Field                    ! Field index -- f_something
    integer :: I                        ! Subscript and loop inductor.
    integer :: Key                      ! Indexes the spec_args vertex.
    integer :: Name                     ! sub_rosa of label of specification,
                                        ! if any, else zero.
    integer :: Son                      ! Some subtree of root.
    integer :: Type                     ! Type of value returned by EXPR
    integer :: Units(2)                 ! Units of value returned by EXPR
    double precision :: Value(2)        ! Value returned by EXPR

    ! Error message codes

    error = 0
    if ( toggle(gen) ) call trace_begin ( "ForwardModelGlobalSetup", root )
    if ( node_id(root) == n_named ) then
      name = subtree(1, root)
      key = subtree(2, root)
    else
      name = 0
      key = root
    end if

    ! "Key" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.

    do i = 2, nsons(key)
      son = subtree(i,key)
      field = get_field_id(son)
      select case ( field )
      case ( f_atmos_der )
        forwardModelInfo%atmos_der = get_boolean(son)
      case ( f_do_conv )
        forwardModelInfo%do_conv = get_boolean(son)
      case ( f_do_freq_avg )
        forwardModelInfo%do_freq_avg = get_boolean(son)
      case ( f_frequency )
        call expr ( subtree(2,son), units, value, type )
        forwardModelInfo%the_freq = value(1)
      case ( f_spect_der )
        forwardModelInfo%spect_der = get_boolean(son)
      case ( f_temp_der )
        forwardModelInfo%temp_der = get_boolean(son)
      case default
        ! Shouldn't get here if the type checker worked
      end select
    end do ! i = 2, nsons(key)
    if ( toggle(gen) ) call trace_end ( "ForwardModelGlobalSetup" )
  end subroutine ForwardModelGlobalSetup

  ! ------------------------------------------  ForwardModelSetup  -----
  subroutine ForwardModelSetup ( Root, VectorDatabase, MatrixDatabase, &
    &                            ForwardModelInfo )
  ! Process the forwardModel specification to produce ForwardModelInfo.

    integer :: Root                     ! of the forwardModel specification.
                                        ! Indexes either a "named" or
                                        ! "spec_args" vertex.
    type(vector_T), dimension(:), intent(inout), target :: VectorDatabase
    type(matrix_Database_T), dimension(:), pointer :: MatrixDatabase
    type(forwardModelInfo_T), intent(inout) :: ForwardModelInfo

    integer :: Field                    ! Field index -- f_something
    logical :: Got(field_first:field_last)   ! "Got this field already"
    integer :: I                        ! Subscript and loop inductor.
    integer :: Key                      ! Indexes the spec_args vertex.
    integer :: Name                     ! sub_rosa of label of specification,
                                        ! if any, else zero.
    integer :: Son                      ! Some subtree of root.

    ! Error message codes

    error = 0
    if ( toggle(gen) ) call trace_begin ( "ForwardModelSetup", root )
    if ( node_id(root) == n_named ) then
      name = subtree(1, root)
      key = subtree(2, root)
    else
      name = 0
      key = root
    end if

    ! "Key" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.

    got = .false.
    do i = 2, nsons(key)
      son = subtree(i,key)
      field = get_field_id(son)
      got(field) = .true.
      select case ( field )
      case default
        ! Shouldn't get here if the type checker worked
      end select
    end do ! i = 2, nsons(key)
    if ( toggle(gen) ) call trace_end ( "ForwardModelSetup" )

  end subroutine ForwardModelSetup

  ! -----------------------------------------------  ForwardModel  -----
! subroutine ForwardModel ( FwdModelInfo, FwdModelExtra, FwdModelIn, &
!   &                       Jacobian, RowBlock, FwdModelOut )
 Subroutine ForwardModel ( FwdModelInfo, FwdModelExtra, FwdModelIn, &
    &                      Jacobian, RowBlock, FwdModelOut, FMC, FMI, TFMI )

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

  use GL6P, only: NG
  use MLSCommon, only: I4, R4, R8
  use L2_TEST_STRUCTURES_M
  use L2PC_FILE_PARAMETERS, only: mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT
  use L2PC_PFA_STRUCTURES, only: K_MATRIX_INFO
  use L2PCdim, only: Nlvl, N2lvl, NSPS, Nptg, NCH, MNP => max_no_phi, &
                     MNM => max_no_mmaf
  use ELLIPSE, only: PHI_TAN, ROC
  use L2_LOAD_M, only: L2_LOAD
  use PTG_FRQ_LOAD_M, only: PTG_FRQ_LOAD
  use COMP_PATH_ENTITIES_M, only: COMP_PATH_ENTITIES
  use REFRACTION_M, only: REFRACTION_CORRECTION
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
                             PATH_DERIVATIVE
  use HYDROSTATIC_MODEL_M, only: HYDROSTATIC_MODEL
  use GET_CHI_ANGLES_M, only: GET_CHI_ANGLES
  use GET_BETA_PATH_M, only: GET_BETA_PATH
  use GEOC_GEOD_CONV_M, only: GEOC_GEOD_CONV
  use RAD_TRAN_M, only: RAD_TRAN
! use RAD_TRAN_WD_M, only: RAD_TRAN_WD
  use FREQ_AVG_M, only: FREQ_AVG
  use CONVOLVE_ALL_M, only: CONVOLVE_ALL
  use NO_CONV_AT_ALL_M, only: NO_CONV_AT_ALL
  use D_HUNT_M, only: HUNT

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

 type(forwardModelInfo_T), intent(in) :: FwdModelInfo ! From ForwardModelSetup
 type(vector_T), intent(in) :: FwdModelExtra, FwdModelIn ! ???
 type(matrix_T), intent(inout), optional :: Jacobian
 integer, intent(in), optional :: RowBlock          ! With which block of
 ! rows of F and Jacobian are we computing? All of them if absent.
 type(vector_T), intent(inout), optional :: FwdModelOut  ! Radiances, etc.

!??? Begin temporary stuff to start up the forward model
  type(fwd_mdl_config), optional :: FMC
  type(fwd_mdl_info), dimension(:), pointer, optional :: FMI
  type(temporary_fwd_mdl_info), dimension(:), pointer, optional :: TFMI
!??? End of temporary stuff to start up the forward model

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

Integer(i4) :: i, j, k, kk, kz, ht_i, mnz, no_tan_hts, ch, Spectag, &
               m, prev_npf, ier, mmaf, si, ptg_i, &
               frq_i, io, klo, jj, l, n, brkpt, no_ele, mid, ilo, ihi, &
               k_info_count, gl_count, ld

Integer(i4) :: ch1, ch2, no_pfa_ch, pfa_ch(2)

Type(path_index)  :: ndx_path(Nptg,mnm)
Type(path_vector) :: z_path(Nptg,mnm),t_path(Nptg,mnm),h_path(Nptg,mnm),  &
                     dhdz_path(Nptg,mnm), spsfunc_path(Nsps,Nptg,mnm),    &
                     n_path(Nptg,mnm),phi_path(Nptg,mnm)

Type(path_derivative) :: dh_dt_path(Nptg,mnm)

Real(r8) :: thbs(10),elev_offset
Real(r8) :: t_script(N2lvl),ref_corr(N2lvl,Nptg),tau(N2lvl), &
            tan_dh_dt(Nlvl,mnm,mxco)

Real(r8) :: dx_dt(Nptg,mxco), d2x_dxdt(Nptg,mxco)

Real(r8) :: h_glgrid(ngt,mnm), t_glgrid(ngt,mnm), z_glgrid(ngt/2)
Real(r8) :: dh_dt_glgrid(ngt,mnm,mxco), dhdz_glgrid(ngt,mnp)

Real(r8) :: ptg_angles(Nptg,mnm), center_angle
Real(r8) :: tan_hts(Nptg,mnm), tan_temp(Nptg,mnm)


! Real(r4) :: K_TEMP(Nch,Nptg,mxco,mnp)
! Real(r4) :: K_ATMOS(Nch,Nptg,mxco,mnp,Nsps)
! Real(r4) :: K_SPECT_DW(Nch,Nptg,mxco,mnp,Nsps),  &
!             K_SPECT_DN(Nch,Nptg,mxco,mnp,Nsps),  &
!             K_SPECT_DNU(Nch,Nptg,mxco,mnp,Nsps)

! ** DEBUG, memory limitations force us to have up to 2 channels
!           only (Replacing Nch by: 2)

Real(r4) :: K_TEMP(02,Nptg,mxco,mnp)
Real(r4) :: K_ATMOS(02,Nptg,mxco,mnp,Nsps)
Real(r4) :: K_SPECT_DW(02,Nptg,mxco,mnp,Nsps),  &
            K_SPECT_DN(02,Nptg,mxco,mnp,Nsps),  &
            K_SPECT_DNU(02,Nptg,mxco,mnp,Nsps)

real(r8) :: I_STAR_ALL(Nch,Nptg)

real(r4) :: K_STAR_ALL(02,20,mxco,mnp,Nptg)      ! 02 should be: Nch
Type(k_matrix_info) :: k_star_info(20)

Type(path_derivative) :: k_temp_frq, k_atmos_frq(Nsps), &
                         k_spect_dw_frq(Nsps), k_spect_dn_frq(Nsps), &
                         k_spect_dnu_frq(Nsps)
!
Type(path_beta), DIMENSION(:,:), POINTER :: beta_path

Real(r8) :: Radiances(Nptg,Nch)
Real(r8) :: e_rad, Zeta, Frq, h_tan, Rad, geoc_lat, r
!
Character (LEN=01) :: CA
Character (LEN=08) :: Name
Character (LEN=16) :: Vname
Character (LEN=80) :: Line
Character (LEN=40) :: Ax, Dtm1, Dtm2

Real(r8), DIMENSION(:), ALLOCATABLE :: RadV, F_grid

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

  ch1 = FMC%Channels_range(1)
  ch2 = FMC%Channels_range(2)
  no_pfa_ch = min(2,ch2-ch1+1)
  do i = 1, no_pfa_ch
    pfa_ch(i) = ch1 + i - 1
  end do

  elev_offset = 0.0                         ! Zero elev_offset in any case

!
! Load the "Frequency gridding by pointing" file ("Bill's" file..)
!
  Call ptg_frq_load(FMC, FMI, Ier)
  if(ier /= 0) goto 99
!
! Convert GeoDetic Latitude to GeoCentric Latitude, and convert both to
! Radians (instead of Degrees). Also compute the effective earth radius.
!
  mmaf = 3                     ! Do only this mmaf (middle phi)
  phi_tan = FMC%phi_tan_mmaf(mmaf)
  Call geoc_geod_conv(TFMI%beta_inc,phi_tan,geoc_lat,E_rad)
!
! Compute the hydrostatic_model on the GL-Grid for all mmaf(s):
!
  thbs(1:) = 0.0
  si = FMI%Surface_index
  thbs(1:si-1) = FMI%Tan_hts_below_surface(1:si-1)
  Call hydrostatic_model(si,FMC%N_lvls,TFMI%no_t,FMC%no_mmaf,FMC%t_indx, &
       FMC%no_tan_hts,geoc_lat,TFMI%Href,TFMI%Zref,FMI%z_grid,thbs, &
       TFMI%t_zeta_basis, TFMI%t_coeff, z_glgrid, h_glgrid, t_glgrid, &
       dhdz_glgrid,dh_dt_glgrid,FMI%tan_press,tan_hts,tan_temp,tan_dh_dt, &
       gl_count, Ier)
  IF(ier /= 0) goto 99
!
  Zeta = -1.666667
  no_tan_hts = FMC%no_tan_hts
  Call Hunt(Zeta,FMI%tan_press,no_tan_hts,jj,i)
  IF(ABS(Zeta-FMI%tan_press(i)) < ABS(Zeta-FMI%tan_press(jj))) jj = i
!
! Compute all path entities for all mmafs and tanget pointings
!
  Call comp_path_entities(FMC%n_lvls,TFMI%no_t,gl_count,ndx_path, &
       z_glgrid,t_glgrid,h_glgrid,dhdz_glgrid,dh_dt_glgrid,        &
       TFMI%atmospheric,TFMI%f_zeta_basis,TFMI%mr_f,            &
       TFMI%no_coeffs_f,tan_hts,no_tan_hts,FMI%n_sps,             &
       TFMI%no_phi_f,TFMI%f_phi_basis,z_path,h_path,t_path,      &
       phi_path,n_path,dhdz_path,dh_dt_path,TFMI%no_phi_t,        &
       TFMI%t_phi_basis,spsfunc_path,TFMI%is_f_log,FMC%no_mmaf,  &
       FMC%phi_tan_mmaf,Ier)
  IF(ier /= 0) goto 99
!
! **********************  MAIN Mmaf Loop *******************
!
! DO l = 1, FMC%no_mmaf
! DO l = mmaf, mmaf                 ! ** DEBUG, only one mmaf
!
    phi_tan = FMC%phi_tan_mmaf(l)
!
    TFMI%t_phi_basis(1:TFMI%no_phi_t) = &
                    TFMI%T_PHI_BASIS_COPY(1:TFMI%no_phi_t) + phi_tan

    DO j = 1, FMI%n_sps
      k = TFMI%no_phi_f(j)
      TFMI%f_phi_basis(1:k,j) = TFMI%F_PHI_BASIS_COPY(1:k,j) + phi_tan
    end do

    k = FMI%mfi + 2
    do j = 1, FMI%no_spectro
      FMI%spectroscopic(j)%PHI_BASIS(1:k) = &
     &                TFMI%S_PHI_BASIS_COPY(1:k,j) + phi_tan
    end do
!
! Compute the ptg_angles (chi) for Antenna convolution, also the derivatives
! of chi w.r.t to T and other parameters
!
    Call get_chi_angles(ndx_path(1:,l),n_path(1:,l),FMI%tan_press,         &
   &     tan_hts(1:,l),tan_temp(1:,l),phi_tan,RoC,TFMI%h_obs,elev_offset, &
   &     tan_dh_dt(1:,l,1:),no_tan_hts,TFMI%no_t,TFMI%t_zeta_basis,si,   &
   &     center_angle,ptg_angles(1:,l),dx_dt,d2x_dxdt,ier)
    IF(ier /= 0) goto 99
!
! Compute the refraction correction scaling matrix for this mmaf:
!
    Call refraction_correction(no_tan_hts, tan_hts(1:,l), h_path(1:,l), &
   &                n_path(1:,l), ndx_path(1:,l), E_rad, ref_corr)
!
    prev_npf = -1
    Radiances(1:Nptg,1:Nch) = 0.0
!
! **********************  MAIN Pointing Loop *******************
!
    DO ptg_i = 1, no_tan_hts-1
!
      k = ptg_i
      h_tan = tan_hts(k,l)
      kk = FMI%no_ptg_frq(k)
!
      if(kk /= prev_npf) then
!
        prev_npf = kk
        DEALLOCATE(k_temp_frq%values,STAT=i)
!
        DEALLOCATE(RadV,F_grid,STAT=i)
        ALLOCATE(RadV(kk),F_grid(kk),STAT=ier)
        IF(ier /= 0) then
          Print *,'** ALLOCATE Error: RadV or F_grid arrays, STAT =',ier
          goto 99
        endif
!
        do j = 1, FMI%n_sps
          DEALLOCATE(k_atmos_frq(j)%values,STAT=i)
          DEALLOCATE(k_spect_dw_frq(j)%values,STAT=i)
          DEALLOCATE(k_spect_dn_frq(j)%values,STAT=i)
          DEALLOCATE(k_spect_dnu_frq(j)%values,STAT=i)
        end do
!
        ALLOCATE(k_temp_frq%values(kk,TFMI%no_t,TFMI%no_phi_t),STAT=ier)
        IF(ier /= 0) then
          Print *,'** ALLOCATE Error: k_temp_frq array, STAT =',ier
          goto 99
        endif
!
        do j = 1, FMI%n_sps
          m = max(1,TFMI%no_phi_f(j))
          i = max(1,TFMI%no_coeffs_f(j))
          ALLOCATE(k_atmos_frq(j)%values(kk,i,m),STAT=ier)
          IF(ier /= 0) then
            Print *,'** ALLOCATE Error: k_atmos_frq, STAT =',ier
            goto 99
          endif
        end do
!
        do m = 1, FMI%n_sps
          j = FMI%spect_atmos(m)
          if(j < 1) CYCLE
          if(.not. FMI%spectroscopic(j)%DER_CALC(FMI%band)) CYCLE
          Vname(1:) = ' '
          Spectag = FMI%spectroscopic(j)%Spectag
          DO
            if(FMI%spectroscopic(j)%Spectag /= Spectag) EXIT
            n = FMI%spectroscopic(j)%no_phi_values
            i = FMI%spectroscopic(j)%no_zeta_values
            CA = FMI%spectroscopic(j)%type
            select case ( CA )
              case ( 'W' )
                Vname = 'k_spect_dw_frq'
                ALLOCATE(k_spect_dw_frq(m)%values(kk,i,n),STAT=ier)
              case ( 'N' )
                Vname = 'k_spect_dn_frq'
                ALLOCATE(k_spect_dn_frq(m)%values(kk,i,n),STAT=ier)
              case ( 'V' )
                Vname = 'k_spect_dnu_frq'
                ALLOCATE(k_spect_dnu_frq(m)%values(kk,i,n),STAT=ier)
              case default
                Ier = -99
                Print *,'** Unknown Spectroscopic element !'
            end select
            IF(ier /= 0) then
              Print *,'** ALLOCATE Error: ',Vname,', STAT =',ier
              goto 99
            ENDIF
            j = j + 1
            if(j > 3 * FMI%n_sps) EXIT
          END DO
        end do

      endif            ! On DEALLOCATE/ALLOCATE cycle
!
! Compute the beta's along the path, for this tanget hight and this mmaf:
!
      no_ele = ndx_path(ptg_i,l)%total_number_of_elements
      Call get_beta_path(ptg_i,FMI%pfa_spectrum,no_ele,FMI%no_ptg_frq, &
     &     FMI%ptg_frq_grid,z_path(ptg_i,l),t_path(ptg_i,l),beta_path,ier)
      IF(ier /= 0) goto 99
!
      k_temp_frq%values = 0.0
      do j = 1, FMI%n_sps
        k_atmos_frq(j)%values = 0.0
        k_spect_dw_frq(j)%values = 0.0
        k_spect_dn_frq(j)%values = 0.0
        k_spect_dnu_frq(j)%values = 0.0
      end do
!
      RadV(1:kk) = 0.0
      F_grid(1:kk) = FMI%ptg_frq_grid(k)%values(1:kk)
!
      do frq_i = 1, kk
!
        Frq = F_grid(frq_i)
!
        Call Rad_Tran(Frq, FMC%N_lvls, h_tan, FMI%n_sps, ndx_path(k,l),  &
       &    z_path(k,l), h_path(k,l), t_path(k,l), phi_path(k,l),&
       &    dHdz_path(k,l), TFMI%earth_ref, beta_path(1:,frq_i),      &
       &    spsfunc_path(1:,k,l), ref_corr(1:,k), TFMI%s_temp, brkpt, &
       &    no_ele, mid, ilo, ihi, t_script, tau, Rad, Ier)
        IF(ier /= 0) goto 99
!
        RadV(frq_i) = Rad
!
! Now, Compute the radiances derivatives:
!
!       CALL Rad_Tran_WD(frq_i,FMI%band,Frq,FMC%N_lvls,FMI%n_sps, &
!      &     FMC%temp_der,FMC%atmos_der,FMC%spect_der,            &
!      &     z_path(k,l),h_path(k,l),t_path(k,l),phi_path(k,l),   &
!      &     dHdz_path(k,l),TFMI%atmospheric,beta_path(1:,frq_i),&
!      &     spsfunc_path(1:,k,l),TFMI%t_zeta_basis,  &
!      &     TFMI%f_zeta_basis,TFMI%no_coeffs_f,   &
!      &     TFMI%mr_f,TFMI%no_t,ref_corr(1:,k),TFMI%no_phi_f,       &
!      &     TFMI%f_phi_basis,TFMI%no_phi_t,TFMI%t_phi_basis,  &
!      &     dh_dt_path(k,l),FMI%spect_atmos,         &
!      &     FMI%spectroscopic,k_temp_frq,k_atmos_frq,k_spect_dw_frq,   &
!      &     k_spect_dn_frq,k_spect_dnu_frq,TFMI%is_f_log,brkpt,       &
!      &     no_ele,mid,ilo,ihi,t_script,tau,ier)
!       IF(ier /= 0) goto 99
!
      end do
!
! Frequency Average the radiances with the appropriate filter shapes
!
      do i = 1, no_pfa_ch
        ch = pfa_ch(i)
        if(FMC%do_frqavg) then
          Call Freq_Avg(F_grid,FMI%F_grid_filter(1:,i),  &
         &     FMI%Filter_func(1:,i),RadV,kk,FMI%no_filt_pts, &
         &     Radiances(ptg_i,ch))
        else
          Radiances(ptg_i,ch) = RadV(1)
        endif
      end do
!
      if(FMC%temp_der) then
!
! Frequency Average the temperature derivatives with the appropriate
! filter shapes
!
        RadV(1:kk) = 0.0
        do i = 1, no_pfa_ch
!         ch = pfa_ch(i)
          ch = i               ! ** DEBUG, memory limitations on MLSGATE
          do j = 1, TFMI%no_phi_t
            do k = 1, TFMI%no_t
              if(FMC%do_frqavg) then
                RadV(1:kk) = k_temp_frq%values(1:kk,k,j)
                Call Freq_Avg(F_grid,FMI%F_grid_filter(1:,i), &
               &              FMI%Filter_func(1:,i),          &
               &              RadV,kk,FMI%no_filt_pts,r)
              else
                r = k_temp_frq%values(1,k,j)
              endif
              k_temp(ch,ptg_i,k,j) = r
            end do
          end do
        end do
!
      endif

      if(FMC%atmos_der) then
!
! Frequency Average the atmospheric derivatives with the appropriate
! filter shapes
!
        do i = 1, no_pfa_ch
!         ch = pfa_ch(i)
          ch = i                 ! ** DEBUG, memory limitations on MLSGATE
          do j = 1, FMI%n_sps
            if(TFMI%atmospheric(j)%der_calc(FMI%band)) THEN
              RadV(1:kk) = 0.0
              do k = 1, TFMI%no_phi_f(j)
                do n = 1, TFMI%no_coeffs_f(j)
                  if(FMC%do_frqavg) then
                    RadV(1:kk) = k_atmos_frq(j)%values(1:kk,n,k)
                    Call Freq_Avg(F_grid,FMI%F_grid_filter(1:,i), &
                   &              FMI%Filter_func(1:,i),          &
                   &              RadV,kk,FMI%no_filt_pts,r)
                  else
                    r = k_atmos_frq(j)%values(1,n,k)
                  endif
                  k_atmos(ch,ptg_i,n,k,j) = r
                end do
              end do
            endif
          end do
        end do
!
      endif

      if(FMC%spect_der) then
!
! Frequency Average the spectroscopic derivatives with the appropriate
! filter shapes
!
        do i = 1, no_pfa_ch
!         ch = pfa_ch(i)
          ch = i                 ! ** DEBUG, memory limitations on MLSGATE
          do m = 1, FMI%n_sps
            j = FMI%spect_atmos(m)
            if(.not.  FMI%spectroscopic(j)%DER_CALC(FMI%band)) CYCLE
            Spectag = FMI%spectroscopic(j)%Spectag
            DO
              if(FMI%spectroscopic(j)%Spectag /= Spectag) EXIT
              RadV(1:kk) = 0.0
              CA = FMI%spectroscopic(j)%type
              do k = 1, FMI%spectroscopic(j)%no_phi_values
                do n = 1, FMI%spectroscopic(j)%no_zeta_values
                  select case ( CA )
                    case ( 'W' )
                      RadV(1:kk) = k_spect_dw_frq(m)%values(1:kk,n,k)
                    case ( 'N' )
                      RadV(1:kk) = k_spect_dn_frq(m)%values(1:kk,n,k)
                    case ( 'V' )
                      RadV(1:kk) = k_spect_dnu_frq(m)%values(1:kk,n,k)
                  end select
                  if(FMC%do_frqavg) then
                      Call Freq_Avg(F_grid,FMI%F_grid_filter(1:,i), &
                   &              FMI%Filter_func(1:,i),&
                   &              RadV,kk,FMI%no_filt_pts,r)
                  else
                    r = RadV(1)
                  endif
                  select case ( CA )
                    case ( 'W' )
                      k_spect_dw(ch,ptg_i,n,k,j) = r
                    case ( 'N' )
                      k_spect_dn(ch,ptg_i,n,k,j) = r
                    case ( 'V' )
                      k_spect_dnu(ch,ptg_i,n,k,j) = r
                  end select
                end do
              end do
              j = j + 1
              if(j > 3 * FMI%n_sps) EXIT
            END DO
          end do
        end do
!
      endif
!
    END DO              ! Pointing Loop
!
! Complete the radiances's last location, also  complete k_temp last
! location as well as k_atmos last location and k_spect_d? last location:
!
    kk = no_tan_hts
    do i = 1, no_pfa_ch
      ch = pfa_ch(i)
      Radiances(kk,ch) = Radiances(kk-1,ch)
      if(FMC%temp_der) then
        k_temp(i,kk,1:TFMI%no_t,1:TFMI%no_phi_t) = &
     &              k_temp(i,kk-1,1:TFMI%no_t,1:TFMI%no_phi_t)
      endif
      if(FMC%atmos_der) then
        do m = 1, FMI%n_sps
          if(TFMI%atmospheric(m)%der_calc(FMI%band)) then
            k = TFMI%no_phi_f(m)
            n = TFMI%no_coeffs_f(m)
            k_atmos(i,kk,1:n,1:k,m)=k_atmos(i,kk-1,1:n,1:k,m)
          endif
        end do
      endif
      if(FMC%spect_der) then
        do m = 1, FMI%n_sps
          j = FMI%spect_atmos(m)
          if(.not.  FMI%spectroscopic(j)%DER_CALC(FMI%band)) CYCLE
          Spectag =  FMI%spectroscopic(j)%Spectag
          DO
            if(FMI%spectroscopic(j)%Spectag /= Spectag) EXIT
            k = FMI%spectroscopic(j)%no_phi_values
            n = FMI%spectroscopic(j)%no_zeta_values
            k_spect_dw(i,kk,1:n,1:k,j)=k_spect_dw(i,kk-1,1:n,1:k,j)
            k_spect_dn(i,kk,1:n,1:k,j)=k_spect_dn(i,kk-1,1:n,1:k,j)
            k_spect_dnu(i,kk,1:n,1:k,j)=k_spect_dnu(i,kk-1,1:n,1:k,j)
            j = j + 1
            if(j > 3 * FMI%n_sps) EXIT
          END DO
        end do
      endif
    end do
!
!  Here comes the Convolution code
!
    DO i = 1, no_pfa_ch
!
      ch = pfa_ch(i)
!
      if(FMC%do_conv) then
!
        Call convolve_all(TFMI%ptg_press,TFMI%atmospheric,FMI%n_sps,   &
       &     FMC%temp_der,FMC%atmos_der,FMC%spect_der,                   &
       &     FMI%tan_press,ptg_angles(1:,l),tan_temp(1:,l), &
       &     dx_dt, d2x_dxdt,FMI%band,center_angle,FMI%fft_pts,          &
       &     Radiances(1:,ch),k_temp(i,1:,1:,1:),k_atmos(i,1:,1:,1:,1:), &
       &     k_spect_dw(i,1:,1:,1:,1:),k_spect_dn(i,1:,1:,1:,1:),    &
       &     k_spect_dnu(i,1:,1:,1:,1:),FMI%spect_atmos,no_tan_hts,  &
       &     k_info_count,i_star_all(i,1:),k_star_all(i,1:,1:,1:,1:), &
       &     k_star_info,TFMI%no_t,TFMI%no_phi_t,TFMI%no_phi_f,   &
       &     FMI%spectroscopic,TFMI%t_zeta_basis,FMI%Xlamda,FMI%Aaap,&
       &     FMI%D1Aaap,FMI%D2Aaap,FMI%Ias,ier)
        IF(ier /= 0) goto 99
!
      else
!
        Call no_conv_at_all(TFMI%ptg_press,FMI%n_sps,FMI%tan_press, &
       &     FMI%band,FMC%temp_der,FMC%atmos_der,FMC%spect_der,      &
       &     Radiances(1:,ch),k_temp(i,1:,1:,1:),                    &
       &     k_atmos(i,1:,1:,1:,1:),k_spect_dw(i,1:,1:,1:,1:),       &
       &     k_spect_dn(i,1:,1:,1:,1:),k_spect_dnu(i,1:,1:,1:,1:),   &
       &     FMI%spect_atmos, no_tan_hts,k_info_count,               &
       &     i_star_all(i,1:), k_star_all(i,1:,1:,1:,1:),            &
       &     k_star_info,TFMI%no_t,TFMI%no_phi_t,                  &
       &     TFMI%no_phi_f,TFMI%t_zeta_basis,TFMI%atmospheric,    &
       &     FMI%spectroscopic)
!
      endif
!
    END DO
!
! END DO                ! Mmaf Loop
!
  if(FMC%temp_der) DEALLOCATE(k_temp_frq%values,STAT=i)
  do j = 1, FMI%n_sps
    if(FMC%atmos_der) DEALLOCATE(k_atmos_frq(j)%values,STAT=i)
    if(FMC%spect_der) then
      DEALLOCATE(k_spect_dw_frq(j)%values,STAT=i)
      DEALLOCATE(k_spect_dn_frq(j)%values,STAT=i)
      DEALLOCATE(k_spect_dnu_frq(j)%values,STAT=i)
    endif
  end do
!
! *** DEBUG Print
!
    if(FMC%do_conv) then
      Print *,'Convolution: ON'
    else
      Print *,'Convolution: OFF'
    endif
!
    if(FMC%Zfrq > 0.0) then
      Frq = FMC%Zfrq
      write(*,901) Frq
901   format(' Frequency Averaging: OFF',/,  &
        &    ' (All computations done at Frq =',f12.4,')')
    else
      Print *,'Frequency Averaging: ON'
    endif
    Print *
!
    tau(1:Nptg) = 0.0
    kk = TFMI%ptg_press%no_lin_values
    tau(1:kk) = dble(TFMI%ptg_press%lin_val(1:kk))
    Call Hunt(Zeta,tau,kk,klo,j)
    IF(ABS(Zeta-tau(j)) < ABS(Zeta-tau(klo))) klo=j
!
    do i = 1, no_pfa_ch
      ch = pfa_ch(i)
      write(*,903) ch,char(92),kk
      write(*,905) (i_star_all(i,k),k=1,kk)
    end do
903 format('ch',i2.2,'_avg_conv_pfa_rad',a1,i2.2)
905 format(4(2x,1pg15.8))
!
    if(.not. ANY((/FMC%temp_der,FMC%atmos_der,FMC%spect_der/))) goto 99
!
    ch = 1
    tau(1:) = 0.0
    do i = 1, k_info_count
      Print *
      Name(1:) = ' '
      Name = k_star_info(i)%name
      if(Name == 'PTAN') CYCLE
      kz = k_star_info(i)%first_dim_index
      mnz = k_star_info(i)%no_zeta_basis
      ht_i = k_star_info(i)%no_phi_basis
      l = LEN_TRIM(Name)
      if(Name(l-1:l) == '_W' .or.  &
     &   Name(l-1:l) == '_N' .or.  &
     &   Name(l-1:l) == '_V' ) then
        Print *,Name
        r = SUM(k_star_all(ch,kz,1:mnz,1:ht_i,klo))
        Print *,'  Sum over all zeta & phi coeff:',sngl(r)
      else
        if(Name == 'TEMP') then
          write(6,913) 'dI_dT',char(92),ht_i
        else
          write(6,913) 'dI_d'//Name(1:l),char(92),ht_i
        endif
        tau(1:) = 0.0
        tau(1:mnz) = k_star_info(i)%zeta_basis(1:mnz)
        Call Hunt(Zeta,tau,mnz,m,j)
        IF(ABS(Zeta-tau(j)) < ABS(Zeta-tau(m))) m=j
        Print *,(k_star_all(ch,kz,m,kk,klo),kk=1,ht_i)
      endif
    end do
913 format(a,a1,i2.2)
!
 99  if(io /= 0) then
       Call ErrMsg(Line,io)
     endif

   Return
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

  end subroutine ForwardModel

! =====     Private Procedures     =====================================
  ! ----------------------------------------------  AnnounceError  -----
  subroutine AnnounceError ( Code, Where, FieldIndex )
    integer, intent(in) :: Code       ! Index of error message
    integer, intent(in) :: Where      ! Where in the tree did the error occur?
    integer, intent(in) :: FieldIndex ! f_...

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref ( where ) )
    call output ( ' ForwardModelSetup complained: ' )
    select case ( code )
    end select
  end subroutine AnnounceError
end module ForwardModelInterface

! $Log$
! Revision 2.7  2001/03/08 19:22:12  zvi
! New ForwardModelInterface with Zvi's code in it ..
!
! Revision 2.6  2001/03/08 03:23:45  vsnyder
! More stuff to work with L2_Load
!
! Revision 2.5  2001/03/08 00:42:09  vsnyder
! Add temporary stuff to use with L2_Load
!
! Revision 2.4  2001/03/07 23:59:52  vsnyder
! Add stuff for SIDS.
!
! Revision 2.3  2001/02/21 00:07:57  vsnyder
! Periodic commit.  Still needs a lot of work.
!
! Revision 2.2  2001/02/08 00:56:11  vsnyder
! Periodic commit.  Still needs a lot of work.
!
! Revision 2.1  2001/02/07 00:52:27  vsnyder
! Initial commit
!
