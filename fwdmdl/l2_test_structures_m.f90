module L2_TEST_STRUCTURES_M
  use MLSCommon, only: I4, R8
  use L2PCDim, only: MAX_NO_PHI, NLVL, NPTG
  use L2PC_PFA_STRUCTURES, only: LIMB_PRESS, ATMOS_COMP, &
 &                               SPECTRO_PARAM, PFA_SLAB
  use PATH_ENTITIES_M, only: PATH_VECTOR
  Implicit NONE
!------------------------------------------------------------------------
!  This is the l2_test structures constructions used for
!  storing the data for the l2_test program
!---------------------------- RCS Ident Info -------------------------------
  PRIVATE :: Id, ModuleName
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName = &
    "$RCSfile$"
!_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
!------------------------------------------------------------
! This structure contains the "Temporary FwdMdlInfo (FMI)", to be replaced
! by state vector type as we progressed.

  type TEMPORARY_FWD_MDL_INFO
    Integer(i4) :: No_t
    Integer(i4) :: No_phi_t
    Integer(i4) :: No_Geometric
    Integer(i4), DIMENSION(:), POINTER :: No_phi_f
    Integer(i4), DIMENSION(:), POINTER :: No_coeffs_f
    Type(LIMB_PRESS) :: Ptg_Press
    Type(ATMOS_COMP), DIMENSION(:), POINTER :: Atmospheric
    Real(r8) :: Zref
    Real(r8) :: beta_inc
    Real(r8) :: elev_183
    Real(r8) :: elev_205
    Real(r8) :: earth_ref
    Real(r8) :: s_temp
    Real(r8) :: h_obs
    Real(r8), DIMENSION(:), POINTER :: Href
    Real(r8), DIMENSION(:), POINTER :: t_zeta_basis
    Real(r8), DIMENSION(:), POINTER :: t_phi_basis
    Real(r8), DIMENSION(:), POINTER :: t_phi_basis_copy
    Logical , DIMENSION(:), POINTER :: is_f_log
    Real(r8), DIMENSION(:,:), POINTER :: t_coeff
    Real(r8), DIMENSION(:,:), POINTER :: f_zeta_basis
    Real(r8), DIMENSION(:,:), POINTER :: f_phi_basis
    Real(r8), DIMENSION(:,:), POINTER :: f_phi_basis_copy
    Real(r8), DIMENSION(:,:), POINTER :: s_phi_basis_copy
    Real(r8), DIMENSION(:,:,:), POINTER :: mr_f
  end type TEMPORARY_FWD_MDL_INFO

!------------------------------------------------------------
! This structure contains the "FMI Config", containing various flags
! and parameters
  type FWD_MDL_CONFIG
    Logical :: do_conv
    Logical :: temp_der
    Logical :: atmos_der
    Logical :: spect_der
    Logical :: do_frqavg
    Real(r8) :: Zfrq
    Integer(i4) :: N_lvls
    Integer(i4) :: No_tan_hts
    Integer(i4) :: No_mmaf
    Integer(i4) :: Channels_range(2)
    Integer(i4) :: Sideband
    Integer(i4), DIMENSION(:), POINTER :: P_indx
    Integer(i4), DIMENSION(:), POINTER :: T_indx
    Character(LEN=80) :: Z
    Character(LEN=80) :: B
    Character(LEN=80) :: Indir
    Character(LEN=80) :: Aaap_file
    Real(r8), DIMENSION(:), POINTER :: Phi_tan_mmaf
  end type FWD_MDL_CONFIG

!------------------------------------------------------------
! This structure contains the "Fwd Model Info", containing all the inputs
! the l2_test ! needs, and which is not state-vaector related
  type FWD_MDL_INFO
    Integer(i4) :: Ias
    Integer(i4) :: mfi
    Integer(i4) :: Band
    Integer(i4) :: N_sps
    Integer(i4) :: Fft_pts
    Integer(i4) :: Max_no_zeta
    Integer(i4) :: No_spectro
    Integer(i4) :: Max_no_phi
    Integer(i4) :: NO_SPECT_CAT
    Integer(i4) :: Surface_index
    Integer(i4) :: No_filt_pts
    Integer(i4) :: no_ptg_frq(Nptg)
    Integer(i4), DIMENSION(:), POINTER :: Spect_atmos
    Real(r8) :: Xlamda
    Real(r8), DIMENSION(:), POINTER :: z_grid
    Real(r8), DIMENSION(:), POINTER :: tan_press
    Real(r8), DIMENSION(:), POINTER :: Tan_hts_below_surface
    Real(r8), DIMENSION(:,:), POINTER :: Filter_func
    Real(r8), DIMENSION(:,:), POINTER :: F_grid_filter
    Real(r8), DIMENSION(:,:), POINTER :: Aaap         ! (maxfft,3)
    Real(r8), DIMENSION(:,:), POINTER :: D1Aaap       ! (maxfft,3)
    Character(LEN=8), DIMENSION(:), POINTER :: Species
    Real(r8), DIMENSION(:,:), POINTER :: D2Aaap       ! (maxfft,3)
    Type(PATH_VECTOR) :: Ptg_frq_grid(Nptg)
    Type(PFA_SLAB), DIMENSION(:), POINTER :: Pfa_spectrum
    Type(SPECTRO_PARAM), DIMENSION(:), POINTER :: Spectroscopic
  end type FWD_MDL_INFO
!
end module L2_TEST_STRUCTURES_M
! $Log$
! Revision 1.0  2001/02/21 21:56:15  zvi
