! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GET_BETA_PATH_M

  implicit NONE
  private
  public :: Create_Beta, Create_Beta_Path
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
  subroutine Get_Beta_Path_Scalar ( frq, p_path, t_path, tanh_path, &
        & beta_group, polarized, gl_slabs, path_inds,               &
        & beta_path, gl_slabs_m, t_path_m, gl_slabs_p, t_path_p,    &
        & t_der_path_flags, dTanh_dT, dBeta_dt_path,                &
        & beta_n_path, dBeta_dw_path, dBeta_dn_path, dBeta_dv_path )

    use Get_Species_Data_m, only: Beta_Group_T
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: R8, RP, IP
    use Physics, only: H_OVER_K
    use SpectroscopyCatalog_m, only: CATALOG_T, LINES

! Inputs:

    real(r8), intent(in) :: Frq          ! frequency in MHz
    real(rp), intent(in) :: P_path(:)    ! path pressures in hPa!
    real(rp), intent(in) :: T_path(:)    ! path temperatures
    real(rp), intent(in) :: Tanh_path(:) ! tanh(0.5*h_over_k*frq / t_path)
    type (slabs_struct), dimension(:,:) :: Gl_slabs
    integer(ip), intent(in) :: Path_inds(:) ! indicies for reading gl_slabs

    type (beta_group_T), dimension(:) :: beta_group

    logical, intent(in) :: Polarized     ! "Don't work on Zeeman-split lines"

! Optional inputs.  GL_SLABS_* are pointers because the caller need not
! allocate them if DBETA_D*_PATH aren't allocated.  They would be
! INTENT(IN) if we could say so.

    type(slabs_struct), pointer :: gl_slabs_m(:,:) ! reduced strength data
!                               for t_path_m
    real(rp), pointer :: t_path_m(:) ! path temperatures for gl_slabs_m
    type(slabs_struct), pointer :: gl_slabs_p(:,:) ! reduced
!                               strength data for t_path_p
    real(rp), pointer :: t_path_p(:) ! path temperatures for gl_slabs_p
    logical, pointer :: t_der_path_flags(:)     ! indicies where temperature
!                               derivatives are needed. Only useful for
!                               subsetting.
    real(rp), pointer :: dTanh_dT(:)    ! dTanh( (-h nu) / (k T) ) / dT on path

! outputs

    real(rp), intent(out) :: beta_path(:,:) ! path beta for each specie

! Optional outputs.  We use ASSOCIATED instead of PRESENT so that the
! caller doesn't need multiple branches.  These would be INTENT(OUT) if
! we could say so.

!{ The variable {\tt dBeta\_dT\_path} isn't really
!  $\frac{\partial\beta}{\partial T}$.  First, we evaluate $\beta$ at
!  $T+\delta T$ and $T-\delta T$.  These are called {\tt bp} and {\tt bm}
!  below.  We assume $\beta$ has the form
!
!  $\beta_0 \left(\frac{T}{T_0}\right) ^n$.  Taking logarithms, we have $\ln
!  (\beta / \beta_0) = n ( \ln ( T / T_0 ) )$.  Assuming $n$ doesn't change
!  too quickly as a function of $T$, you can solve for $n$ from estimates of
!  $\beta(T)$ and $\beta(T+\delta T)$ (you need $\beta$ at two temperatures
!  because you don't know $\beta_0$ and $T_0$).  Actually, we use three
!  estimates for $\frac{\partial\beta}{\partial T}$, \emph{viz}.
!  $(\beta(T+\delta T) - \beta(T-\delta T)) / ( 2 \delta T)$,
!  $(\beta(T+\delta T) - \beta(T) ) / \delta T$, and $(\beta(T) - \beta(T -
!  \delta T)) / \delta T$, from which we compute three estimates for $n$. 
!  It's the value of $n$ that's returned in {\tt dBeta\_dT\_path}, not
!  $\frac{\partial\beta}{\partial T}$.  $\frac{\partial\beta}{\partial T}$ is
!  actually assembled in {\tt dRad\_tran\_dT}.

    real(rp), pointer :: dBeta_dT_path(:,:) ! Temperature
    real(rp), pointer :: beta_n_path(:,:)   ! n in beta = beta_0 (T/T_0)^n
    real(rp), pointer :: dBeta_dw_path(:,:) ! line width
    real(rp), pointer :: dBeta_dn_path(:,:) ! line width t dep.
    real(rp), pointer :: dBeta_dv_path(:,:) ! line position

! Local variables.

    type(catalog_t), pointer :: Catalog
    real(rp), pointer :: dBdn(:), dBdv(:), dBdw(:) ! slices of dBeta_d*_path
    integer(ip) :: I, J, K, N, IB, No_of_lines, &
              &    No_mol, N_path
    real(rp) :: BB, BP, BM
    real(rp) :: T
    real(rp), dimension(size(path_inds)) :: betam, betap
    real(rp), dimension(size(path_inds)) :: tanh1_p, tanh1_m, TM, TP

! begin the code

    no_mol = size(beta_group)
    n_path = size(path_inds)

    nullify ( dBdn, dBdv, dBdw )

    ! Determine size of the LineWidths array.
    no_of_lines = 0
    do i = 1, no_mol
      do n = 1, beta_group(i)%n_elements
        ib = beta_group(i)%cat_index(n)
        catalog => gl_slabs(1,ib)%catalog ! GL_slabs(:,ib) all have the same catalog
        no_of_lines = max(no_of_lines,size(catalog%lines))
      end do
    end do

    do i = 1, no_mol
      if ( associated(dBeta_dn_path) ) then
        dBdn => dBeta_dn_path(:,i)
        dBdn = 0.0_rp
      end if
      if ( associated(dBeta_dv_path) ) then
        dBdv => dBeta_dv_path(:,i)
        dBdv = 0.0_rp
      end if
      if ( associated(dBeta_dw_path) ) then
        dBdw => dBeta_dw_path(:,i)
        dBdw = 0.0_rp
      end if
      beta_path(:,i) = 0.0_rp
      do n = 1, beta_group(i)%n_elements
        ib = beta_group(i)%cat_index(n)
        catalog => gl_slabs(1,ib)%catalog ! GL_slabs(:,ib) all have the same catalog
        no_of_lines = size(catalog%lines)

        ! The polarized and nonpolarized calls are split so that Create_Beta_Path
        ! can go to the "don't bother to test polarized" loop.  Create_Beta_Path
        ! is the inner loop of the forward model.  It also avoids construction
        ! of an array temp (polarized .and. catalog%polarized).
        if ( polarized ) then
          call create_beta_path ( path_inds, p_path, t_path, frq,                &
            & beta_group(i)%ratio(n), gl_slabs(:,ib), tanh_path, beta_path(:,i), &
            & catalog%polarized, dBeta_dw=dBdw, dBeta_dn=dBdn, dBeta_dv=dBdv )
        else
          call create_beta_path ( path_inds, p_path, t_path, frq,                &
            & beta_group(i)%ratio(n), gl_slabs(:,ib), tanh_path, beta_path(:,i), &
            & dBeta_dw=dBdw, dBeta_dn=dBdn, dBeta_dv=dBdv )
        end if
      end do
    end do

    if ( associated(beta_n_path) ) then

      beta_n_path(1:n_path,:) = 0.0
      tm = t_path_m(path_inds)
      tp = t_path_p(path_inds)
      ! compute path hyperbolic tangents
      tanh1_p = tanh(0.5_rp * h_over_k * frq / tp)
      tanh1_m = tanh(0.5_rp * h_over_k * frq / tm)

      do i = 1, no_mol
        betam = 0.0
        betap = 0.0
        do n = 1, beta_group(i)%n_elements
          ib = beta_group(i)%cat_index(n)
          catalog => gl_slabs(1,ib)%catalog ! GL_slabs(:,ib) all have the same catalog
          no_of_lines = size(catalog%lines)
          if ( polarized ) then
            call create_beta_path ( path_inds, p_path, tm, frq,             &
              & beta_group(i)%ratio(n), gl_slabs_m(:,ib), tanh1_m, betam,   &
              & catalog%polarized, t_der_path_flags,                        &
              & dBeta_dw=null(), dBeta_dn=null(), dBeta_dv=null() )
            call create_beta_path ( path_inds, p_path, tp, frq,             &
              & beta_group(i)%ratio(n), gl_slabs_p(:,ib), tanh1_p, betap,   &
              & catalog%polarized, t_der_path_flags,                        &
              & dBeta_dw=null(), dBeta_dn=null(), dBeta_dv=null() )
          else
            call create_beta_path ( path_inds, p_path, tm, frq,                &
              & beta_group(i)%ratio(n), gl_slabs_m(:,ib), tanh1_m, betam,      &
              & path_flags=t_der_path_flags, dBeta_dw=null(), dBeta_dn=null(), &
              & dBeta_dv=null() )
            call create_beta_path ( path_inds, p_path, tp, frq,                &
              & beta_group(i)%ratio(n), gl_slabs_p(:,ib), tanh1_p, betap,      &
              & path_flags=t_der_path_flags, dBeta_dw=null(), dBeta_dn=null(), &
              & dBeta_dv=null() )
          end if
        end do ! n
        do j = 1 , n_path
          k = path_inds(j)
          t  = t_path(j)
          bm = betam(j)
          bp = betap(j)
          bb = beta_path(j,i)
          n = 0
          if ( bm > 0.0 ) n = 1
          if ( bb > 0.0 ) n = n + 2
          if ( bp > 0.0 ) n = n + 4
          select case ( n )
          case ( 7 ) ! bp > 0.0 .and. bb > 0.0 .and. bm > 0.0
            beta_n_path(j,i) =  0.25 * ( &           ! Weighted average of
              &       Log(bp/bb)/Log(tp(j)/t) + &    ! Estimate over [temp+10,temp]
              & 2.0 * Log(bp/bm)/Log(tp(j)/tm(j))+ & ! Estimate over [temp+10,temp-10]
              &       Log(bb/bm)/Log(t/tm(j)) )      ! Estimate over [temp,temp-10]
          case ( 6 ) ! bp > 0.0 .and. bb > 0.0
            beta_n_path(j,i) = Log(bp/bb)/Log(tp(j)/t) ! Estimate over [temp+10,temp]
          case ( 5 ) ! bm > 0.0 .and. bp > 0.0
            beta_n_path(j,i) = Log(bp/bm)/Log(tp(j)/tm(j)) ! Estimate over [temp+10,temp-10]
          case ( 3 ) ! bm > 0.0 .and. bb > 0.0
            beta_n_path(j,i) = Log(bb/bm)/Log(t/tm(j)) ! Estimate over [temp,temp-10]
          case default
            beta_n_path(j,i) = 0.0
          end select
        end do ! j
      end do ! i
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

!{ The variable {\tt Beta\_Path} isn't really $\beta$.  It lacks a factor of
!  $\tanh\left( \frac{h \nu}{2 k T}\right)$.  This is put in after we
!  compute the weighted average over species, saving as many multiplies as
!  there are species.

    complex(rp), intent(out) :: Beta_path(-1:,:,:) ! path beta for each species
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
  subroutine Get_Beta_Path_Cloud ( Frq, p_path, t_path,  tt_path,       &
        & beta_group, path_inds, beta_path_cloud,                       &
        & w0_path, tt_path_c,                                           & 
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
    real(rp), intent(in) :: tt_path(:,:)   ! scating source func on gl grids

    integer(ip), intent(in) :: Path_inds(:) ! indicies for reading gl_slabs

    type (beta_group_T), dimension(:) :: beta_group
    type (ForwardModelConfig_T) ,intent(in) :: FWDMODELCONF

    INTEGER, intent(in) :: IPSD(:)
    REAL(rp), intent(in)  :: WC(:,:)
    REAL(rp) :: W0       ! SINGLE SCATTERING ALBEDO
    REAL(rp) :: PHH(fwdModelConf%num_scattering_angles)   ! PHASE FUNCTION

! outputs

    real(rp), intent(out) :: beta_path_cloud(:) ! cloud extinction
    real(rp), intent(out) :: w0_path(:)         ! single scattering albedo
    real(rp), intent(out) :: tt_path_c(:)       ! scattering source func coarse grids

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

    beta_path_cloud = 0.0_rp
    w0_path         = 0.0_rp
    tt_path_c       = 0.0_rp

        do j = 1, n_path
          k = path_inds(j)

          call get_beta_cloud (Frq, t_path(k),                       &
                          &  WC(:,k), IPSD(k), NC, NU, NUA, NAB, NR, &
                          &  cld_ext, W0, PHH                )      

            beta_path_cloud(j) = beta_path_cloud(j) + cld_ext 
            w0_path(j)         = w0_path(j) + W0 
            tt_path_c(j)       = tt_path_c(j)       + tt_path(k,1)            

         end do

  end subroutine Get_Beta_Path_Cloud

! ----------------------------------------------  Create_beta  ---------

  subroutine Create_beta ( pressure, Temp, Fgr, slabs_0, tanh1, &
         &                 Beta_Value, Polarized, dTanh_dT,     &
         &                 dBeta_dT, dBeta_dw, dBeta_dn, dBeta_dv )

!  For a given frequency and height, compute beta_value function.
!  This routine should be called for primary and image separately.

!  If you change Create_Beta, change Create_Beta_Path the same way!
!  Create_Beta and Create_Beta_Path should do the same thing.

!  The reason for the existence of Create_Beta_Path is that by moving the
!  loop for the path down into Create_Beta a substantial improvement in
!  running time was achieved.  Create_Beta is in the inner loop of the
!  forward model.

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: RP, R8, IP
    use Molecules, only: L_N2, L_Extinction, L_O2
    use SLABS_SW_M, only: DVOIGT_SPECTRAL, VOIGT_LORENTZ, &
      & SLABS_DT, SLABS_LINES, SLABSWINT_DT, SLABSWINT_LINES
    use SpectroscopyCatalog_m, only: Catalog_T, Lines

! Inputs:
    real(rp), intent(in) :: pressure   ! pressure in hPa
    real(rp), intent(in) :: temp       ! temperature in K
    real(r8), intent(in) :: fgr        ! frequency in MHz
    type(slabs_struct), intent(in) :: slabs_0 ! contains, among others:

!    catalog        ! Pointer to spectroscopy catalog
!    v0s(:)         ! pressure shifted line centers
!    x1(:)          ! Doppler width
!    y(:)           ! ratio Pressure to Doppler widths
!    yi(:)          ! Interference coefficients
!    expa(:)        ! exponential argument / frequency (not used)
!    slabs1(:)      ! strengths
!    dslabs1_dv0(:) ! strength derivative wrt line position

    real(rp), intent(in) :: tanh1      ! tanh(frq*expa/2)

    ! "Don't do this line" -- same size as slabs_0%catalog%lines
    logical, intent(in), optional :: Polarized(:)

! optional inputs for temperature derivatives
    ! (h frq) / (k T^2) ( tanh( (h frq) / (k T) ) - 1/tanh( (h frq) / (k T) ) ):
    real(rp), intent(in), optional :: dTanh_dT ! 1/tanh d/dT tanh
! outputs
    real(rp), intent(out) :: beta_value
! optional outputs
    real(rp), optional, intent(out) :: DBETA_DT ! Temperature derivative
    real(rp), optional, intent(out) :: DBETA_DW ! line width derivative
    real(rp), optional, intent(out) :: DBETA_DN ! temperature dependence deriv
    real(rp), optional, intent(out) :: DBETA_DV ! line position derivative

! -----     Local variables     ----------------------------------------

    type(catalog_t), pointer :: Catalog
    real(rp), pointer :: Cont(:)    ! Continuum parameters
    integer(ip) :: LN_I             ! Line index
    integer(ip) :: NL               ! no of lines
    integer(ip) :: molecule         ! molecule id

    real(rp) :: bv, dNu, dT, dw, dn, ds, dbdw, dbdn, dbdv

!----------------------------------------------------------------------------

    if ( present(dBeta_dw) .or. present(dBeta_dn) .or. present(dBeta_dv) ) then

      dbdw = 0.0_rp
      dbdn = 0.0_rp
      dbdv = 0.0_rp

    end if

!  Setup absorption coefficients function
!  Now get the beta_value:

    catalog => slabs_0%catalog
    cont => catalog%continuum
    nl = size(catalog%lines)
    select case ( catalog%molecule )
    case ( l_n2 ) ! ...........................................  Dry Air

      if ( present(dBeta_dT) ) then
        call abs_cs_n2_cont_dT ( cont, temp, pressure, fgr, beta_value, dBeta_dT )
      else
        beta_value = abs_cs_n2_cont(cont,Temp,Pressure,Fgr)
      end if

    case ( l_extinction ) ! ................................  Extinction

      beta_value = 1.0_rp
      if ( present(dBeta_dT) ) dBeta_dT = 0.0_rp
      return

    case ( l_o2 ) ! ................................................  O2

      if ( present(dBeta_dT) ) then
        call abs_cs_o2_cont_dT ( cont, temp, pressure, fgr, beta_value, dBeta_dT )
      else
        beta_value = abs_cs_o2_cont(cont,Temp,Pressure,Fgr)
      end if

    case default ! ..............................................  Other

      if ( present(dBeta_dT) ) then
        call abs_cs_cont_dT ( cont, temp, pressure, fgr, beta_value, dBeta_dT )
      else
        beta_value = abs_cs_cont(cont,Temp,Pressure,Fgr)
      end if

    end select

    if ( nl < 1 ) return

    if ( present(dBeta_dw) .or. present(dBeta_dn) .or. present(dBeta_dv) ) then

      do ln_i = 1, nl

        if ( present(polarized) ) then
          if ( polarized(ln_i) ) cycle
        end if

        dNu = Fgr - slabs_0%v0s(ln_i)

        if ( abs(slabs_0%y(ln_i))+0.666666_rp*abs(slabs_0%x1(ln_i)*dNu) &
        & > 100.0_rp ) then
          call Voigt_Lorentz ( dNu, slabs_0%v0s(ln_i), slabs_0%x1(ln_i),         &
            &  slabs_0%yi(ln_i), slabs_0%y(ln_i), lines(catalog%lines(ln_i))%w,  &
            &  Temp, tanh1, slabs_0%slabs1(ln_i), bv, slabs_0%dslabs1_dv0(ln_i), &
            &  dw, dn, ds )
        else
          call DVoigt_Spectral ( dNu, slabs_0%v0s(ln_i), slabs_0%x1(ln_i),       &
            &  slabs_0%yi(ln_i), slabs_0%y(ln_i), lines(catalog%lines(ln_i))%w,  &
            &  Temp, tanh1, slabs_0%slabs1(ln_i), bv, slabs_0%dslabs1_dv0(ln_i), &
            &  dw, dn, ds )
        end if

        beta_value = beta_value + bv
        dbdw = dbdw + dw
        dbdn = dbdn + dn
        dbdv = dbdv + ds

      end do

      if ( present(dBeta_dw)) dBeta_dw = dbdw
      if ( present(dBeta_dn)) dBeta_dn = dbdn
      if ( present(dBeta_dv)) dBeta_dv = dbdv

    else if ( present(polarized) ) then  ! No spectroscopy derivatives required

      ! This is split into polarized and unpolarized loops, to avoid
      ! having "if ( present(polarized) )" inside the loop.  This is
      ! the inner loop of the forward model.

      if ( .not. present(dBeta_dT) ) then
        if ( maxval(ABS(slabs_0%yi)) < 1.0e-06_rp ) then
          beta_value = slabs_lines ( Fgr, slabs_0, tanh1, polarized )
        else
          beta_value = slabswint_lines ( Fgr, slabs_0, tanh1, polarized )
        end if
      else ! Temperature derivatives needed
        if ( maxval(ABS(slabs_0%yi)) < 1.0e-06_rp ) then
          do ln_i = 1, nl
            if ( polarized(ln_i) ) cycle
            call slabs_dT ( &
              & Fgr, lines(catalog%lines(ln_i))%v0, slabs_0%v0s(ln_i),          &
              & slabs_0%x1(ln_i), tanh1, slabs_0%slabs1(ln_i), slabs_0%y(ln_i), &
              & slabs_0%dv0s_dT(ln_i), slabs_0%dx1_dT(ln_i), dTanh_dT,          &
              & slabs_0%dslabs1_dT(ln_i), slabs_0%dy_dT(ln_i), bv, dT )
            beta_value = beta_value + bv
            dBeta_dT = dBeta_dT + dT
          end do
        else
          do ln_i = 1, nl
            if ( polarized(ln_i) ) cycle
            call slabswint_dT ( &
              & Fgr, lines(catalog%lines(ln_i))%v0, slabs_0%v0s(ln_i),          &
              & slabs_0%x1(ln_i), tanh1, slabs_0%slabs1(ln_i), slabs_0%y(ln_i), &
              & slabs_0%yi(ln_i), slabs_0%dv0s_dT(ln_i), slabs_0%dx1_dT(ln_i),  &
              & dTanh_dT, slabs_0%dslabs1_dT(ln_i), slabs_0%dy_dT(ln_i),        &
              & slabs_0%dyi_dT(ln_i), bv, dT )
            beta_value = beta_value + bv
            dBeta_dT = dBeta_dT + dT
          end do
        end if
      end if

    else ! Still no spectroscopy derivatives required

      if ( .not. present(dBeta_dT) ) then
        if ( maxval(ABS(slabs_0%yi)) < 1.0e-06_rp ) then
          beta_value = slabs_lines ( Fgr, slabs_0, tanh1 )
        else
          beta_value = slabswint_lines ( Fgr, slabs_0, tanh1 )
        end if
      else ! Temperature derivatives needed
        if ( maxval(ABS(slabs_0%yi)) < 1.0e-06_rp ) then
          do ln_i = 1, nl
            call slabs_dT ( &
              & Fgr, lines(catalog%lines(ln_i))%v0, slabs_0%v0s(ln_i),          &
              & slabs_0%x1(ln_i), tanh1, slabs_0%slabs1(ln_i), slabs_0%y(ln_i), &
              & slabs_0%dv0s_dT(ln_i), slabs_0%dx1_dT(ln_i), dTanh_dT,          &
              & slabs_0%dslabs1_dT(ln_i), slabs_0%dy_dT(ln_i), bv, dT )
            beta_value = beta_value + bv
            dBeta_dT = dBeta_dT + dT
          end do
        else
          do ln_i = 1, nl
            call slabswint_dT ( &
              & Fgr, lines(catalog%lines(ln_i))%v0, slabs_0%v0s(ln_i),          &
              & slabs_0%x1(ln_i), tanh1, slabs_0%slabs1(ln_i), slabs_0%y(ln_i), &
              & slabs_0%yi(ln_i), slabs_0%dv0s_dT(ln_i), slabs_0%dx1_dT(ln_i),  &
              & dTanh_dT, slabs_0%dslabs1_dT(ln_i), slabs_0%dy_dT(ln_i),        &
              & slabs_0%dyi_dT(ln_i), bv, dT )
            beta_value = beta_value + bv
            dBeta_dT = dBeta_dT + dT
          end do
        end if
      end if

    end if

  end Subroutine Create_beta

! -----------------------------------------  Create_beta_path  ---------

  subroutine Create_beta_path ( Path_inds, Pressure, Temp, Fgr, Ratio, &
         &   Slabs_0, Tanh1, Beta_value, Polarized, Path_flags,        &
         &   Slabs_p, Tanh1_p, Slabs_m, Tanh1_m,                       &
         &   dBeta_dw, dBeta_dn, dBeta_dv  )

!  For a given frequency and height, compute beta_value function.
!  This routine should be called for primary and image separately.

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: RP, R8, IP
    use Molecules, only: L_N2, L_Extinction, L_O2
    use SLABS_SW_M, only: DVOIGT_SPECTRAL, VOIGT_LORENTZ, &
      & SLABS_LINES, SLABSWINT_LINES
    use SpectroscopyCatalog_m, only: LINES

! Inputs:
    integer, intent(in) :: Path_inds(:)! Which Pressures to use
    real(rp), intent(in) :: Pressure(:)! pressure in hPa on the find path grid
    real(rp), intent(in) :: Temp(:)    ! temperature in K along the path
    real(r8), intent(in) :: Fgr        ! frequency in MHz
    real(rp), intent(in) :: Ratio      ! Isotope ratio
    type(slabs_struct), intent(in) :: Slabs_0(:) ! contains, among others:

!    Catalog        ! Pointer to catalog
!    v0s(:)         ! pressure shifted line centers
!    x1(:)          ! Doppler width
!    y(:)           ! ratio Pressure to Doppler widths
!    yi(:)          ! Interference coefficients
!    expa(:)        ! exponential argument / frequency (not used)
!    slabs1(:)      ! strengths
!    dslabs1_dv0(:) ! strength derivative wrt line position

    real(rp), intent(in) :: Tanh1(:)   ! tanh(frq*expa/2)
    ! "Don't do this line" -- same size as slabs_0(k)%catalog%lines
    logical, intent(in), optional :: Polarized(:)
    logical, intent(in), optional :: Path_Flags(:) ! to do on fine path -- default true

! optional inputs for temperature derivatives
    type(slabs_struct), intent(in), optional :: Slabs_p(:), Slabs_m(:)
    real(rp), intent(in), optional :: Tanh1_p(:), Tanh1_m(:) ! tanh(frq*expa/2)
! outputs
    real(rp), intent(inout) :: Beta_value(:)
! optional outputs
    real(rp), pointer :: dBeta_dw(:) ! line width derivative
    real(rp), pointer :: dBeta_dn(:) ! temperature dependence deriv
    real(rp), pointer :: dBeta_dv(:) ! line position derivative

! -----     Local variables     ----------------------------------------

    real(rp), pointer :: Cont(:) ! continuum parameters
    integer :: J, K   ! Subscript, loop inductor
    integer :: LN_I   ! Line index
    integer :: NL     ! no of lines

    real(rp) :: dNu, tp, bp, tm, bm, bv, dw, dn, ds, dbdw, dbdn, dbdv

!----------------------------------------------------------------------------

    nl = size(slabs_0(1)%v0s) ! All of the slabs have the same number of lines

    do j = 1, size(path_inds)
      k = path_inds(j)
      if ( present(path_flags) ) then
        if ( .not. path_flags(k) ) cycle
      end if

      tp = Temp(j) + 10.0_rp
      tm = Temp(j) - 10.0_rp

!  Setup absorption coefficients function
!  Now get the beta_value:

      cont => slabs_0(k)%catalog%continuum
      select case ( slabs_0(k)%catalog%molecule )
      case ( l_n2 ) ! ...........................................  Dry Air
! This nominally gets multiplied by "ratio**2" but in practice this
! function is for all isotopics forms of N2 hence the ratio is one.

        beta_value(j) = beta_value(j) + abs_cs_n2_cont(cont,Temp(j),Pressure(k),Fgr)

      case ( l_extinction ) ! ................................  Extinction

        beta_value(j) = beta_value(j) + ratio
        cycle

      case ( l_o2 ) ! ................................................  O2

        beta_value(j) = beta_value(j) + ratio * abs_cs_o2_cont(cont,Temp(j),Pressure(k),Fgr)

      case default ! ..............................................  Other

        beta_value(j) = beta_value(j) + ratio * abs_cs_cont(cont,Temp(j),Pressure(k),Fgr)

      end select

      if ( nl < 1 ) cycle

      if ( associated(dBeta_dw) .or. associated(dBeta_dn) .or. associated(dBeta_dv) ) then

        if ( present(polarized) ) then
          do ln_i = 1, nl

            if ( polarized(ln_i) ) cycle

            dNu = Fgr - slabs_0(k)%v0s(ln_i)

            if ( abs(slabs_0(k)%y(ln_i))+0.666666_rp*abs(slabs_0(k)%x1(ln_i)*dNu) &
            & > 100.0_rp ) then
              call Voigt_Lorentz ( dNu, slabs_0(k)%v0s(ln_i), slabs_0(k)%x1(ln_i), &
                &  slabs_0(k)%yi(ln_i), slabs_0(k)%y(ln_i), &
                &  lines(slabs_0(k)%catalog%lines(ln_i))%w, temp(j), &
                &  tanh1(j), slabs_0(k)%slabs1(ln_i), bv, slabs_0(k)%dslabs1_dv0(ln_i), &
                &  dw, dn, ds )
            else
              call DVoigt_Spectral ( dNu, slabs_0(k)%v0s(ln_i), slabs_0(k)%x1(ln_i), &
                &  slabs_0(k)%yi(ln_i), slabs_0(k)%y(ln_i), &
                &  lines(slabs_0(k)%catalog%lines(ln_i))%w, temp(j), &
                &  tanh1(j), slabs_0(k)%slabs1(ln_i), bv, slabs_0(k)%dslabs1_dv0(ln_i), &
                &  dw, dn, ds )
            end if

            beta_value(j) = beta_value(j) + ratio * bv
            dbdw = dbdw + dw
            dbdn = dbdn + dn
            dbdv = dbdv + ds

          end do

        else

          do ln_i = 1, nl

            dNu = Fgr - slabs_0(k)%v0s(ln_i)

            if ( abs(slabs_0(k)%y(ln_i))+0.666666_rp*abs(slabs_0(k)%x1(ln_i)*dNu) &
            & > 100.0_rp ) then
              call Voigt_Lorentz ( dNu, slabs_0(k)%v0s(ln_i), slabs_0(k)%x1(ln_i), &
                &  slabs_0(k)%yi(ln_i), slabs_0(k)%y(ln_i), &
                &  lines(slabs_0(k)%catalog%lines(ln_i))%w, temp(j), &
                &  tanh1(j), slabs_0(k)%slabs1(ln_i), bv, slabs_0(k)%dslabs1_dv0(ln_i), &
                &  dw, dn, ds )
            else
              call DVoigt_Spectral ( dNu, slabs_0(k)%v0s(ln_i), slabs_0(k)%x1(ln_i), &
                &  slabs_0(k)%yi(ln_i), slabs_0(k)%y(ln_i), &
                &  lines(slabs_0(k)%catalog%lines(ln_i))%w, temp(j), &
                &  tanh1(j), slabs_0(k)%slabs1(ln_i), bv, slabs_0(k)%dslabs1_dv0(ln_i), &
                &  dw, dn, ds )
            end if

            beta_value(j) = beta_value(j) + ratio * bv
            dbdw = dbdw + dw
            dbdn = dbdn + dn
            dbdv = dbdv + ds

          end do

        end if

        if ( associated(dBeta_dw)) dBeta_dw(j) = ratio * dbdw
        if ( associated(dBeta_dn)) dBeta_dn(j) = ratio * dbdn
        if ( associated(dBeta_dv)) dBeta_dv(j) = ratio * dbdv

      else                ! No derivatives required

        ! Add in sum of betas for all the lines
        if ( maxval(ABS(slabs_0(k)%yi)) < 1.0e-06_rp ) then
          beta_value(j) = beta_value(j) + &
            & ratio * slabs_lines ( Fgr, slabs_0(k), tanh1(j), polarized )
        else
          beta_value(j) = beta_value(j) + &
            & ratio * slabswint_lines ( Fgr, slabs_0(k), tanh1(j), polarized )
        end if

      end if
    end do ! j = 1, size(path_inds)

  end Subroutine Create_beta_path

!     =====     Private Procedures     =================================

  ! ----------------------------------------------  Abs_CS_Cont  -----

  ! Compute the general continuum contribution
  pure function Abs_CS_Cont ( Cont, Temperature, Pressure, Frequency ) &
    & result(Abs_CS_Cont_r)
  ! real(rp) function Abs_CS_Cont ( Cont, Temperature, Pressure, Frequency )
    use MLSCommon, only: RP

    real(rp), intent(in) :: CONT(:)     ! continuum parameters
    real(rp), intent(in) :: TEMPERATURE ! in Kelvin
    real(rp), intent(in) :: PRESSURE    ! in mbar
    real(rp), intent(in) :: FREQUENCY   ! in MegaHertz
    real(rp) :: Abs_CS_Cont_r

    Abs_CS_Cont_r = cont(1) * pressure * pressure * frequency * frequency * &
      & ( (300.0_rp / temperature)**cont(2) )

  end function Abs_CS_Cont

  ! -------------------------------------------  Abs_CS_Cont_dT  -----

  ! Compute the general continuum contribution and its temperature derivative
  subroutine Abs_CS_Cont_dT ( Cont, Temperature, Pressure, Frequency, &
    & Beta, dBeta_dT )
    use MLSCommon, only: RP

    real(rp), intent(in) :: CONT(:)     ! continuum parameters
    real(rp), intent(in) :: TEMPERATURE ! in Kelvin
    real(rp), intent(in) :: PRESSURE    ! in mbar
    real(rp), intent(in) :: FREQUENCY   ! in MegaHertz
    real(rp), intent(out) :: Beta, dBeta_dT

    real(rp) :: Onedt ! 1/T

!{ Let $\theta = \frac{300}T$.  Then the general continuum contribution is
!  $\beta = c_1 p^2 \nu^2 \theta^{c_2}$.  Noticing that
!  $\frac{\partial \theta}{\partial T} = -\frac{\theta}T$, we have
!  $\frac{\partial \beta}{\partial T} = -\beta \frac{c_2}T$.

    onedt = 1.0 / temperature
    beta = cont(1) * pressure * pressure * frequency * frequency * &
      & ( (300.0_rp * onedt)**cont(2) )

    dBeta_dT = -beta * cont(2) * onedt

  end subroutine Abs_CS_Cont_dT

  ! -------------------------------------------  Abs_CS_N2_Cont  -----

  ! Compute the N2 continuum contribution
  pure function Abs_CS_N2_Cont ( Cont, Temperature, Pressure, Frequency ) &
    & result(Abs_CS_N2_Cont_r)
  ! real(rp) Function Abs_CS_N2_cont ( Cont, Temperature, Pressure, Frequency )
    use MLSCommon, only: RP

    real(rp), intent(in) :: CONT(:)     ! continuum parameters
    real(rp), intent(in) :: TEMPERATURE ! in Kelvin
    real(rp), intent(in) :: PRESSURE    ! in mbar
    real(rp), intent(in) :: FREQUENCY   ! in MegaHertz
    real(rp) :: Abs_CS_N2_Cont_r

    REAL(rp) :: THETA, FSQR, FSXT

    theta = 300.0_rp / temperature
    fsqr = frequency * frequency
    fsxt = fsqr * theta
    Abs_CS_N2_Cont_r = pressure * pressure * fsqr * (theta**cont(2)) * &
                   & ( cont(1) * exp(-cont(3) * fsxt) + &
                   &   cont(4) * exp(-cont(5) * fsxt) * &
                   & (cont(6)**2 + fsqr))

  end function Abs_CS_N2_Cont

  ! ----------------------------------------  Abs_CS_N2_Cont_dT  -----

  ! Compute the N2 continuum contribution and its temperature derivative
  subroutine Abs_CS_N2_Cont_dT ( Cont, Temperature, Pressure, Frequency, &
    & Beta, dBeta_dT )

    use MLSCommon, only: RP

    real(rp), intent(in) :: CONT(:)     ! continuum parameters
    real(rp), intent(in) :: TEMPERATURE ! in Kelvin
    real(rp), intent(in) :: PRESSURE    ! in mbar
    real(rp), intent(in) :: FREQUENCY   ! in MegaHertz
    real(rp), intent(out) :: Beta, dBeta_dT

    real(rp) :: E1, E2, F, FSQR, FSXT, OneDT, THETA

!{ Let $\theta = \frac{300}T$ and $f = p^2 \nu^2 \theta^{c_2}$. Then the N2
!  continuum contribution to $\beta$ is\\
!  $\beta = f ( c_1 e^{-c_3 \nu^2 \theta} +
!   c_4 e^{-c_5  \nu^2 \theta} (c_6^2 + \nu^2) )$.  Noticing that
!  $\frac{\partial \theta}{\partial T} = -\frac{\theta}T$, we have\\
!  $\frac{\partial \beta}{\partial T} = \frac1T \left ( -\beta c_2 +
!   f \nu^2 \theta \left( c_3 c_1 e^{-c_3 \nu^2 \theta} +
!    c_5 c_4 (c_6^2 + \nu^2) e^{-c_5  \nu^2 \theta} \right) \right)$.

    onedt = 1.0 / temperature
    theta = 300.0_rp * onedt
    fsqr = frequency * frequency
    fsxt = fsqr * theta
    f = pressure * pressure * fsqr * (theta**cont(2))
    e1 = cont(1) * exp(-cont(3) * fsxt)
    e2 = cont(4) * exp(-cont(5) * fsxt) * (cont(6)**2 + fsqr)
    beta = f * ( e1 + e2 )

    dBeta_dT = onedt * ( f * fsxt * ( cont(3) * e1 + cont(5) * e2 ) &
      &                  - beta * cont(2) )

  end subroutine Abs_CS_N2_Cont_dT

  ! -------------------------------------------  Abs_CS_O2_Cont  -----

  ! Compute the O2 continuum contribution
  pure function Abs_CS_O2_Cont ( Cont, Temperature, Pressure, Frequency ) &
    & result(Abs_CS_O2_Cont_r)
  ! real(rp) Function ABS_CS_O2_CONT ( Cont, Temperature, Pressure, Frequency )
    use MLSCommon, only: RP

    real(rp), intent(in) :: CONT(:)     ! continuum parameters
    real(rp), intent(in) :: TEMPERATURE ! in Kelvin
    real(rp), intent(in) :: PRESSURE    ! in mbar
    real(rp), intent(in) :: FREQUENCY   ! in MegaHertz
    real(rp) :: Abs_CS_O2_Cont_r

    real(rp) :: THETA, FSQR

    theta = log ( 300.0_rp / temperature )
    fsqr = frequency * frequency
    Abs_CS_O2_Cont_r = cont(1) * pressure * pressure * fsqr * exp(theta*cont(2)) &
                   & / (fsqr + (cont(3) * pressure * exp(theta*cont(4)) )**2 )

  end function Abs_CS_O2_Cont

  ! ----------------------------------------  Abs_CS_O2_Cont_dT  -----

  ! Compute the O2 continuum contribution
  subroutine Abs_CS_O2_Cont_dT ( Cont, Temperature, Pressure, Frequency, &
      & Beta, dBeta_dT )

    use MLSCommon, only: RP

    real(rp), intent(in) :: CONT(:)     ! continuum parameters
    real(rp), intent(in) :: TEMPERATURE ! in Kelvin
    real(rp), intent(in) :: PRESSURE    ! in mbar
    real(rp), intent(in) :: FREQUENCY   ! in MegaHertz
    real(rp), intent(out) :: Beta, dBeta_dT

    real(rp) :: D, F, FSQR, Onedt, THETA

!{ Let $\theta = \frac{300}T$, $f = (c_3 p \theta^{c_4})^2$ and
!  $D = \frac1{\nu^2 + f}$.
!  Then the O2 continuum contribution to beta
!  is $\beta = c_1 p^2 \nu^2 \theta^{c_2} D$. Noticing that
!  $\frac{\partial \theta}{\partial T} = -\frac{\theta}T$, we have
!  $\frac{\partial \beta}{\partial T} = \frac{\beta}T
!   ( 2 c_4 f D - c_2)$.

    onedt = 1.0 / temperature
    theta = log( 300.0_rp * onedt )
    fsqr = frequency * frequency
    f = ( cont(3) * pressure * exp(cont(4)*theta) )**2
    d = 1.0 / (fsqr + f)
    beta = cont(1) * pressure * pressure * fsqr * exp(cont(2)*theta) * d

    dBeta_dT = beta * onedt * ( 2.0 * cont(4) * f * d - cont(2) )

  end subroutine Abs_CS_O2_Cont_dT

!-----------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module GET_BETA_PATH_M

! $Log$
! Revision 2.54  2004/03/20 04:08:55  vsnyder
! Steps along the way to analytic temperature derivatives
!
! Revision 2.53  2004/03/20 01:15:23  jonathan
!  minor changes
!
! Revision 2.52  2004/03/19 04:07:31  vsnyder
! Fix some blunders re dNu for spectral derivatives
!
! Revision 2.51  2004/03/19 00:47:02  vsnyder
! Use line center instead of pressure-shifted line center in a few places
!
! Revision 2.50  2004/02/27 22:47:32  bill
! fixed bug in n2 continuum calc
!
! Revision 2.49  2004/01/26 21:57:00  vsnyder
! Improve TeXnicalities
!
! Revision 2.48  2003/12/07 19:46:10  jonathan
! update for use in 2D cloud FWM
!
! Revision 2.47  2003/08/20 21:12:39  bill
! fixed tanh1 bug associated with T-ders
!
! Revision 2.46  2003/07/15 17:50:30  vsnyder
! Callers need t_der_path_flags to be a pointer
!
! Revision 2.45  2003/07/14 22:45:09  vsnyder
! Scale BP, BM by isotope ratio in t_power computation
!
! Revision 2.44  2003/07/11 22:43:37  vsnyder
! Multiply the continuum-derived Beta by the isotope ratio
!
! Revision 2.43  2003/07/09 22:47:43  vsnyder
! Make separate branches for the polarized and nonpolarized cases, so
! we don't need to check "if (present(polarized))" inside the loop.  This
! is the inner loop for the full forward model.
!
! Revision 2.42  2003/07/07 19:53:51  vsnyder
! Move newly-public Create_Beta and Create_Beta_Path above the 'Private
! Procedures' comment
!
! Revision 2.41  2003/07/07 19:08:38  vsnyder
! Make Create_Beta and Create_Beta_Path public
!
! Revision 2.40  2003/07/07 16:47:01  pwagner
! Moved declaration of 3 function results inside body to appease NAG
!
! Revision 2.39  2003/07/04 02:47:50  vsnyder
! Move Create_Beta here, add Create_Beta_Path routine
!
! Revision 2.38  2003/06/27 22:09:19  vsnyder
! Move allocation of LineWidths out of loops; simplify exponent calculation
!
! Revision 2.37  2003/06/18 17:23:40  bill
! fixed NAG associated bug
!
! Revision 2.36  2003/06/18 14:44:53  bill
! added subsetting feature for T-ders
!
! Revision 2.35  2003/06/02 22:41:33  vsnyder
! Remove unused symbols
!
! Revision 2.34  2003/05/16 23:51:51  livesey
! Now uses molecules rather than spectags
!
! Revision 2.33  2003/05/15 03:28:52  vsnyder
! Moved some stuff up to FullForwardModel because Get_d_Deltau_pol_dT needs it
!
! Revision 2.32  2003/05/10 00:48:09  vsnyder
! Add TeXnicalities
!
! Revision 2.31  2003/05/09 20:07:07  vsnyder
! Correct kind parameter for H, specify intent for GL_slabs and Beta_group
!
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
