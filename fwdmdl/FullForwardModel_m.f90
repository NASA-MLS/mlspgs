! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module FullForwardModel_m

  ! This module contains the `full' forward model.

  use GLNP, only: NG, NGP1
  implicit NONE
  private
  public :: FullForwardModel

!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile: FullForwardModel_m.f90,v $"
!------------------------------------------------------------------------------

  ! Private parameters:

  integer, parameter :: Max_New = NGP1 ! Maximum new points in coarse path in
    ! addition to ones derived from the preselected zeta grid: the minimum
    ! zeta point plus NG H_Glgrid intersections below the tangent zeta.

  logical, parameter :: NaN_Fill = .false. ! Fill automatic arrays with NaN

  logical, parameter :: Zeta_Only = .true. ! Don't determine intersections of
    ! the LOS with vertical boundaries of a QTM-based grid.  Calculations
    ! along the line of sight are mostly from one intersection with a constant
    ! zeta surface to the next one.  If vertical boundary intersections are
    ! included, the maximum diameter of the QTM would need to be included in
    ! those numbers of elements along the path.

  logical :: Same_Facets ! The size of the union of the parts of the QTM
    ! traversed by the paths of all MIFs in the MAF is less than twice the
    ! length of the longest one.  Thereby, we can copy the subset of the
    ! state vector of the MAF into Grids_Tmp and Grids_F once, instead of
    ! allocating Grids_Tmp and Grids_F at first for the maximum number of
    ! vertices needed for any MIF, and copying the appropriate subset of the
    ! state vector in each pointing.  Sizes of automatic arrays in
    ! FullForwardModelAuto depend upon them.

  ! This is set from the configuration:

  logical, private :: WrongTrapezoidal = .true.
    ! If GL is not used on a panel, the rectangular estimate used to cancel
    ! the singularity at the tangent point is replaced by a trapezoidal
    ! estimate.  Originally, this was done incorrectly, using delta
    ! s ~ ds/dh dh/dz delta z, which is a rectangular quadrature approximation
    ! of delta s.  We have delta s, so we ought to use it.  This flag
    ! indicates whether the incorrect computation ought to be preserved.

  logical, private, parameter :: Dump_ds = .false.
    ! The difference between the trapezoidal methods is that one uses del_s
    ! and the other uses dsdz_c * del_zeta.  If this parameter is true,
    ! these values are written to units 42 and 43.

contains

  ! -------------------------------------------- FullForwardModel -----

  subroutine FullForwardModel ( FwdModelConf, FwdModelIn, FwdModelExtra, &
                             &  FwdModelOut, FmStat, Jacobian, ExtraJacobian, &
                             &  Hessian )

  ! This gets a little bit of data and then computes the sizes for quantities
  ! in the full forward model proper, FullForwardModelAuto.

    use Allocate_Deallocate, only: Deallocate_Test, Test_Allocate, Test_Deallocate
    use Check_QTM_m, only: Check_QTM
    use Compute_Z_PSIG_M, only: Compute_Z_PSIG
    use ForwardModelConfig, only: Dump, ForwardModelConfig_T
    use ForwardModelIntermediate, only: ForwardModelStatus_T
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Get_Magnetic_Field_m, only: Get_Magnetic_Field
    use Get_Species_Data_M, only:  Get_Species_Data
    use HessianModule_1, only: Hessian_T
    use HGridsDatabase, only: HGrid_T
    use Interpolate_MIF_to_Tan_Press_m, only: Get_Lines_of_Sight
    use Intrinsic, only: Lit_Indices, L_ECRtoFOV, L_PhiTan, L_PTan, L_ScECR, &
      & L_TScat, L_VMR, L_Wrong, L_MagneticField
    use Load_SPS_Data_M, only: DestroyGrids_T, Dump, EmptyGrids_T, Grids_T, &
      & Load_One_Item_Grid, Load_SPS_Data
    use MatrixModule_1, only: Matrix_T
    use MLSKinds, only: RP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MLSStringLists, only: SwitchDetail
    use Molecules, only: L_CloudIce
    use MoreMessage, only: MLSMessage
    use Path_Representation_m, only: Facets_and_Vertices_t, Path_t, Union_Paths
    use QTM_Facets_Under_Path_m, only: QTM_Facets_Under_Path
    use Tangent_Pressures_m, only: Tangent_Pressures
    use Toggles, only: Emit, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    USE VectorsModule, ONLY: Vector_T, VectorValue_T, &
         GetVectorQuantityIndexByType

    type(forwardModelConfig_T), intent(inout) :: FwdModelConf
    type(vector_T), intent(in) :: FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian
    type(matrix_T), intent(inout), optional :: ExtraJacobian
    type(hessian_T), intent(inout), optional :: Hessian

    real(rp), allocatable :: Z_PSIG(:)    ! Surfs from Temperature, tangent grid
                                ! and species grids, sans duplicates.
    real(rp), allocatable :: Tan_Press(:) ! Pressures corresponding to Z_PSIG

    type (Facets_and_Vertices_t), allocatable :: F_and_V(:) ! Facets and
                                ! vertices under each path through
                                ! the QTM that is used to integrate the
                                ! radiative-transfer equation.
    type (Facets_and_Vertices_t), allocatable :: F_and_V_MIF(:) ! Facets and
                                ! vertices under each path through
                                ! the QTM that is a MIF line of sight.
    type (Grids_T) :: Grids_tmp ! All the coordinates for TEMP
    type (Grids_T) :: Grids_f   ! All the coordinates for VMR
    type (Grids_T) :: Grids_IWC ! All the coordinates for IWC
    type (Grids_T) :: Grids_mag ! All the coordinates for Magnetic field
    type (Grids_T) :: Grids_n   ! All the spectroscopy(N) coordinates
    type (Grids_T) :: Grids_v   ! All the spectroscopy(V) coordinates
    type (Grids_T) :: Grids_w   ! All the spectroscopy(W) coordinates
    type (HGrid_T), pointer :: QTM_HGrid  ! HGrid that has finest QTM resolution.
    type (Path_t), allocatable :: QTM_Paths(:)! LOS through QTM at tangent
                                ! pressures given by Tan_Press
    type (VectorValue_T), pointer :: CloudIce ! Ice water content
    type (Facets_and_Vertices_t), allocatable :: Path_Union(:) ! Facets and
                                ! vertices under all paths through the QTM.
    type (VectorValue_T), pointer :: ECRtoFOV ! Rotation matrix from ECR to
                                ! field-of-view (minor frame quantity)
    type (VectorValue_T), pointer :: PhiTan   ! Tangent geodAngle component of
                                ! state vector (minor frame quantity)
    type (VectorValue_T), pointer :: PTan     ! Tangent pressure component of
                                ! state vector (minor frame quantity)
    type (Path_T), allocatable :: Q_LOS(:)    ! Line-of-sight, for QTM,
                                ! minor frame quantity
    type (VectorValue_T), pointer :: ScECR_MIF ! Instrument position in ECR
                                ! (minor frame quantity)
    type (VectorValue_T), pointer :: Temp     ! Temperature component of
                                ! state vector
    type (VectorValue_T), pointer :: magneticfield     ! magnetic field component of
                                ! state vector

    character(127) :: ERMSG     ! From allocate

    integer :: Dump_Conf        ! for debugging, from -Sfmconf
    integer :: K                ! Loop inductor and subscript
    integer :: Longest_QTM_path ! Path with greatest number of QTM vertices
                                ! adjacent to the line of sight
    integer :: Max_C            ! Length of longest possible coarse path,
                                ! Z_PSIG & Min Zeta & surface Zeta
    integer :: Max_F            ! Length of longest possible combined coarse &
                                ! fine path (all npf<max_f)
    integer :: MaxVert          ! Number of points in gl-refined vertical grid:
                                ! (NLVL-1) * (NG+1) + 1, i.e., 1 + NG
                                ! per level, except the last,
                                ! where there's no GL space.
    integer :: Me = -1          ! String index for trace
    integer :: Nlvl             ! Number of levels in coarse zeta grid
    integer :: No_Mol           ! Number of molecules
    integer :: No_sv_p_T        ! number of phi basis for temperature
    integer :: No_Tan_Hts       ! Number of tangent heights
    integer :: NoUsedChannels   ! Number of channels used
    integer :: N_T_zeta         ! Number of zetas for temperature
    integer :: Stat             ! From allocate
    integer :: SurfaceTangentIndex  ! Index in tangent grid of earth's
                                    ! surface
    INTEGER :: Sv_T_len         ! Number of t_phi*t_zeta in the window
    integer :: No_sv_p_H        ! number of phi basis for magnetic field
    integer :: N_H_zeta         ! Number of zetas for magnetic field
    integer :: Sv_H_len         ! Number of t_phi*t_zeta in the window
    integer :: S_A   ! Multiplier for atmos derivative sizes, 0 or 1
    integer :: S_H   ! Multiplier for atmos second derivative sizes, 0 or 1
    integer :: S_I   ! Multiplier for ice/cloud sizes, 0 or 1
    integer :: S_LC  ! Multiplier for line center deriv sizes, 0 or 1
    integer :: S_LW  ! Multiplier for line width deriv sizes, 0 or 1
    integer :: S_PFA ! Multiplier for PFA sizes, 0 or 1
    integer :: S_P   ! Multiplier for polarized sizes, 0 or 1
    integer :: S_QTM ! Multiplier for sizes of QTM-related stuff, 0 or 1
    integer :: S_T   ! Multiplier for temp derivative sizes, 0 or 1
    integer :: S_TD  ! Multiplier for temp dependence deriv sizes, 0 or 1
    integer :: S_TG  ! Multiplier for TScat generation sizes, 0 or 1
    integer :: S_TS  ! Multiplier for using TScat tables, 0 or 1
    integer :: S_M   ! Multiplier for magnetic field derivatives, 0 or 1
    logical :: UsingQTM         ! Temperature and all species have QTM hGrids

    ! Flags for various derivatives
    logical :: Atmos_Der, Atmos_Second_Der, PTan_Der, Spect_Der
    logical :: Spect_Der_Center, Spect_Der_Width, Spect_Der_Width_TDep
    ! What severity is not having derivative in "first" state vector?
    INTEGER, PARAMETER :: DerivativeMissingFromState = MLSMSG_Error
    INTEGER :: compute_magfield_derivatives
    
    compute_magfield_derivatives = GetVectorQuantityIndexByType(fwdmodelin, &
         l_magneticfield, noError = .TRUE.)

    IF (compute_magfield_derivatives /= 0) THEN
       fwdmodelconf%hmag_der = .TRUE.
       fwdmodelconf%htheta_der = .TRUE.
       fwdmodelconf%hphi_der = .TRUE.
    ELSE
       fwdmodelconf%hmag_der = .FALSE.
       fwdmodelconf%htheta_der = .FALSE.
       fwdmodelconf%hphi_der = .FALSE.
    ENDIF
    
!    PRINT '(a)','=========begin derivative flag settings======'
!    PRINT *,'atmos der setting ',fwdmodelconf%atmos_der
!    PRINT *,'temperature der setting ',fwdmodelconf%temp_der
!    PRINT *,'h-mag der setting ',fwdmodelconf%hmag_der
!    PRINT *,'h-theta der setting ',fwdmodelconf%htheta_der
!    PRINT *,'h-phi der setting ',fwdmodelconf%hphi_der
!    PRINT '(a)','=========end derivative flag settings========'

    call trace_begin ( me, 'FullForwardModel, MAF=', index=fmstat%maf, &
      & cond=toggle(emit) ) ! set by -f command-line switch

    ! Create the data structures for the species.  Get the
    ! spectroscopy parameters from the state vector.  Get a pointer to the
    ! temperature quantity from the state vector into the configuration.
    ! This has to be done AFTER DeriveFromForwardModelConfig, which is
    ! invoked from ForwardModelWrappers.

    call get_species_data ( fwdModelConf, fwdModelIn, fwdModelExtra )

    ! Get a shorter handle for the temperature quantity
    temp => fwdModelConf%temp%qty
    
    magneticfield => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
         & quantityType=l_magneticfield, config=fwdModelConf)
    n_h_zeta = magneticfield%TEMPLATE%nosurfs
    no_sv_p_h = magneticfield%TEMPLATE%noinstances
    sv_h_len = no_sv_p_h * n_h_zeta

    no_mol = size(fwdModelConf%beta_group)
    noUsedChannels = size(fwdModelConf%channels)

    ! Check whether temperature and all the mixing ratios either all have QTM
    ! HGrid, or none do.  If so, find the QTM HGrid with the finest resolution.
    ! (We currently require that if they're QTM they all have the same grid.)

    call check_QTM ( fwdModelConf, QTM_HGrid, usingQTM )

    ! Compute the preselected integration grid (all surfs from temperature,
    ! tangent grid and species).  Tan_Press is thrown in for free.
    if ( fwdModelConf%generateTScat ) then
      ! Make sure the TScat zeta is in the preselected grid.
      call compute_Z_PSIG ( fwdModelConf, z_psig,                  &
                          & GetQuantityForForwardModel (            &
                          &  fwdModelOut, quantityType=l_TScat,     &
                          &  signal=fwdModelConf%signals(1)%index,  &
                          &  config=fwdModelConf ) )

    else
      call compute_Z_PSIG ( fwdModelConf, z_psig )
    end if
    nlvl = SIZE(z_psig)
    PRINT *,'con fig file ',fwdmodelconf%tangentgrid%surfs
    call tangent_pressures ( fwdModelConf, z_psig, no_tan_hts,   &
                           & surfaceTangentIndex, tan_press )
!    PRINT *,'====== I got to the tangent pressure bit=====',no_tan_hts
!    PRINT *,z_psig
!    PRINT *,tan_press
!    PRINT *,surfacetangentindex
!    PRINT *,'============================================='

    ! Find the phiTan quantity in the state vector.  The phiTan quantity is
    ! only used to compute the instance window for the temperature quantity.
    ! If the temperature grid is QTM, phiTan isn't actually used because the
    ! entire QTM is used, instead of having an instance window.
    phitan => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_phitan, config=fwdModelConf, &
      & instrumentModule=fwdModelConf%signals(1)%instrumentModule )

    ! Compute some derivative flags.  %derivOK needs to be computed BEFORE
    ! loading Grids_f.
    fwdModelConf%beta_group%qty%derivOK = fwdModelConf%moleculeDerivatives .and. &
      & ( fwdModelConf%beta_group%qty%foundInFirst .or. present(extraJacobian) )
    fwdModelConf%temp%derivOK = present ( jacobian ) .and. FwdModelConf%temp_der
    if ( fwdModelConf%temp%derivOK .and. .not. fwdModelConf%temp%foundInFirst ) &
      & call MLSMessage ( DerivativeMissingFromStateFun(), moduleName, &
        & 'With config(%S): Temperature derivative requested but temperature is not in "first" state vector', &
        & datum=fwdModelConf%name )

    pTan => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_ptan, foundInFirst=pTan_der, config=fwdModelConf, &
      & instrumentModule=fwdModelConf%signals(1)%instrumentModule )
    pTan_der = pTan_der .and. present ( jacobian )

    wrongTrapezoidal = fwdModelConf%trapezoid == l_wrong

    ! Compute some sizes
    n_t_zeta = temp%template%noSurfs ! Number of zeta levels for temperature
    n_h_zeta = magneticfield%template%nosurfs
    if ( usingQTM ) then

      ! These computations are here instead of in Both_Sidebands_Setup
      ! because they calculate No_Sv_P_T, which is bound for many arrays
      ! in FullForwardModelAuto.

      ! Find the facets and vertices under all the paths through the QTM.
      ! The largest number of vertices in any of them is the horizontal
      ! extent of grids related to temperature.

      allocate ( QTM_Paths(no_tan_hts), stat=stat, errmsg=ermsg )
      call test_allocate ( stat, moduleName, "QTM_paths", 1, no_tan_hts, &
        & ermsg=ermsg )
      allocate ( F_and_V(no_tan_hts), stat=stat, errmsg=ermsg )
      call test_allocate ( stat, moduleName, "F_and_V", 1, no_tan_hts, &
        & ermsg=ermsg )
      allocate ( F_and_V_MIF(ptan%template%noSurfs), stat=stat, errmsg=ermsg )
      call test_allocate ( stat, moduleName, "F_and_V_MIF", 1, no_tan_hts, &
        & ermsg=ermsg )
      ! Path_Union is allocatable instead of explicit shape so that it can be
      ! moved to F_And_V or F_And_V_MIF if all paths have (approximately) the
      ! same sets of facets and vertices.
      allocate ( path_union(1), stat=stat, errmsg=ermsg )
      call test_allocate ( stat, moduleName, "Path_Union", 1, 1, &
        & ermsg=ermsg )

      ! Get Q_LOS ( C, U ) vectors for all MIFs for the MAF
      ! The third column of ECRtoFOV is an unit vector in the direction of
      ! the line of sight from the instrument.
      ECRtoFOV => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_ECRtoFOV, config=fwdModelConf )
      allocate ( Q_LOS(ECRtoFOV%template%noSurfs), stat=stat, errmsg=ermsg )
      call test_allocate ( stat, moduleName, "Q_LOS", 1, 1, &
        & ermsg=ermsg )
      scECR_MIF => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_scECR, config=fwdModelConf )
      do k = 1, size(Q_LOS)
        Q_LOS(k)%lines(1,1)%xyz = ScECR_MIF%value3(1:3,k,fmstat%maf) ! C vec
        Q_LOS(k)%lines(2,1)%xyz = ECRtoFOV%value3(7:9,k,fmstat%maf)  ! U vec
      end do

      call get_lines_of_sight ( fmStat%maf, pTan, tan_press, Q_LOS, &
                              & QTM_paths )

      ! Compute the maximum horizontal extent of arrays related to temperature.
      ! Temperature is extracted from the state vector and put onto a grid that
      ! has a horizontal extent equal to the maximum number of vertices adjacent
      ! to any path through the QTM.  All the other temperature-related
      ! quantities have the same horizontal grid.

      ! If MIFs in a MAF might cover very different sets of facets of the QTM,
      ! we need to know the longest one because components of Grids_Tmp
      ! and Grids_F provide sizes of automatic arrays in FullForwardModelAuto,
      ! and then re-create Grids_Tmp and Grids_F for each pointing.
      ! Otherwise, assume that each MIF of a MAF covers approximately the same
      ! set of facets.
      no_sv_p_T = 0
      longest_QTM_path = 1
      do k = 1, no_tan_hts ! size(QTM_paths), size(F_and_V)
        call QTM_Facets_Under_Path ( QTM_paths(k), QTM_HGrid%QTM_tree, &
          & F_and_V(k) )
        if ( size(f_and_v(k)%vertices) > no_sv_p_T ) then
          no_sv_p_T = size(f_and_v(k)%vertices)
          longest_QTM_path = k
        end if
        call union_paths ( path_union(1), f_and_v(k) )
      end do
      same_facets = size(path_union(1)%vertices) <= 2 * no_sv_p_T
      if ( same_facets ) then
        call move_alloc ( path_union, f_and_v )
        longest_QTM_path = 1
        no_sv_p_T = size(F_and_V(1)%vertices)
      end if

      ! Copy the temperature from the state vector via the configuration into
      ! Grids_tmp.  Automatic extents of arrays in FullForwardModelAuto
      ! depend upon components of Grids_Tmp.  If Same_Facets, these values will
      ! be used for all paths, else new ones will be loaded for each path.

      call load_one_item_grid ( grids_tmp, temp, fmStat%maf, phitan, fwdModelConf, &
        & setDerivFlags=.true., subset=f_and_v(longest_QTM_path)%vertices, &
        & short=.not. same_facets )

      ! Copy the mixing ratios from the state vector via the configuration
      ! into Grids_f.  Automatic extents of arrays in FullForwardModelAuto
      ! depend upon components of Grids_F.  If Same_Facets, these values will
      ! be used for all paths, else new ones will be loaded for each path.

      call load_sps_data ( FwdModelConf, phitan, fmStat%maf, grids_f, &
        & subset=f_and_v(longest_QTM_path)%vertices, short=.not. same_facets )

      ! Calculate facets and vertices under MIF paths.  These are used for
      ! the hydrostatic equilibrium calculation that is necessary to compute
      ! MIF-based Chi angles, and their derivatives.

      do k = 1, ptan%template%noSurfs
        call QTM_Facets_Under_Path ( Q_LOS(k), QTM_HGrid%QTM_tree, &
          & F_and_V_MIF(k) )
        call union_paths ( path_union(1), F_and_V_MIF(k) )
      end do
      if ( size(path_union(1)%vertices) <= 2 * size(F_and_V_MIF(k)%vertices) ) &
        & call move_alloc ( path_union, F_and_V_MIF )

    else ! not QTM

      same_facets = .true. ! Used to decide whether to re-load Grids_Tmp and
                           ! Grids_F for each path, even in the non-QTM case.

      ! Copy the temperature quantity from the state vector via the
      ! configuration into Grids_tmp.  Automatic extents of arrays in
      ! FullForwardModelAuto depend upon components of Grids_Tmp.

      CALL load_one_item_grid ( grids_tmp, temp, fmStat%maf, phitan, &
           & fwdModelConf, setDerivFlags=.TRUE. )
      
      CALL load_one_item_grid ( grids_mag, magneticfield, fmStat%maf, phitan, &
        & fwdModelConf, setDerivFlags=.true. )

      ! Copy the mixing ratios from the state vector via the configuration into
      ! Grids_f.  Automatic extents of arrays in FullForwardModelAuto depend
      ! upon components of Grids_F.

      call load_sps_data ( FwdModelConf, phitan, fmStat%maf, grids_f )

      allocate ( QTM_Paths(0), stat=stat, errmsg=ermsg )
      call test_allocate ( stat, moduleName, "QTM_paths", 1, no_tan_hts, &
        & ermsg=ermsg )
      ! One element of F_and_V is needed because routines in Convolve_All_m
      ! test whether F_and_V%Vertices is allocated.
      allocate ( F_and_V(1), stat=stat, errmsg=ermsg )
      call test_allocate ( stat, moduleName, "F_and_V", 1, 1, &
        & ermsg=ermsg )
      no_sv_p_T = grids_tmp%l_p(1)  ! size of temperature's horizontal grid ==
                                    ! windowFinish - windowStart + 1.
      no_sv_p_H = grids_mag%l_p(1)  ! size of magnetic horizontal grid ==
                                    ! windowFinish - windowStart + 1.
    end if ! not QTM

    sv_t_len = no_sv_p_T * n_t_zeta ! Number of temperature values
    sv_h_len = no_sv_p_H * n_h_zeta ! Number of magnetic values

    spect_der = present ( jacobian ) .and. FwdModelConf%spect_der
    spect_der_center = spect_der .and. size(fwdModelConf%lineCenter) > 0
    spect_der_width = spect_der .and. size(fwdModelConf%lineWidth) > 0
    spect_der_width_TDep = spect_der .and. size(fwdModelConf%lineWidth_TDep) > 0

    atmos_der = present ( jacobian ) .and. FwdModelConf%atmos_der
    atmos_second_der = present (hessian ) .and. FwdModelConf%atmos_second_der

    ! Determine where we can compute derivatives.  We know for sure we have
    ! a place to put them if fwdModelConf%Beta_group%qty%foundInFirst, but
    ! we also have a place to put them otherwise if ExtraJacobian is present
    ! and there's a column there for it.

    if ( atmos_der ) then
      do k = 1, no_mol
        if ( fwdModelConf%moleculeDerivatives(k) ) then
          if ( .not. fwdModelConf%beta_group(k)%qty%derivOK ) then
            call MLSMessage ( DerivativeMissingFromStateFun(), moduleName, &
              & 'With config(%S): ' // &
              & '%S derivative requested but %S is not in "first" state vector, ' // &
              & 'or ExtraJacobian is not present', &
              & datum=(/ fwdModelConf%name, lit_indices(k), lit_indices(k) /) )
          end if
        else if ( .not. usingQTM ) then
          ! Turn off deriv_flags where we don't want molecule derivatives,
          ! so we don't have to look at fwdModelConf%moleculeDerivatives
          grids_f%deriv_flags(grids_f%l_v(k-1)+1:grids_f%l_v(k)) = .false.
        end if
      end do                        ! Loop over major molecules
    end if                          ! Want derivatives for atmos

    if ( spect_der_center .and. &
       & any(.not. fwdModelConf%lineCenter%qty%foundInFirst) ) then
       ! If the vector quantity for the desired molecule is not in the
       ! "first" state vector, check whether there's an ExtraJacobian.
       if ( .not. present(extraJacobian) ) &
         & call MLSMessage ( DerivativeMissingFromStateFun(), moduleName, &
           & 'With config(%S): ' // &
           & 'derivative requested but line center is not in "first" state ' // &
           & 'vector, or ExtraJacobian is not present', &
           & datum=(/ fwdModelConf%name /) )
    end if
    if ( spect_der_width .and. &
       & any(.not. fwdModelConf%lineWidth%qty%foundInFirst) ) then
       ! If the vector quantity for the desired molecule is not in the
       ! "first" state vector, check whether there's an ExtraJacobian.
       if ( .not. present(extraJacobian) ) &
         & call MLSMessage ( DerivativeMissingFromStateFun(), moduleName, &
           & 'With config(%S): ' // &
           & 'derivative requested but line width is not in "first" state ' // &
           & 'vector, or ExtraJacobian is not present', &
           & datum=(/ fwdModelConf%name /) )
    end if
    if ( spect_der_width_TDep .and. &
       & any(.not. fwdModelConf%lineWidth_TDep%qty%foundInFirst) ) then
       ! If the vector quantity for the desired molecule is not in the
       ! "first" state vector, check whether there's an ExtraJacobian.
       if ( .not. present(extraJacobian) ) &
         & call MLSMessage ( DerivativeMissingFromStateFun(), moduleName, &
           & 'With config(%S): ' // &
           & 'derivative requested but line width temperature dependence is ' // &
           & 'not in "first" state vector, or ExtraJacobian is not present', &
           & datum=(/ fwdModelConf%name /) )
    end if

    call Get_Magnetic_Field ( fwdModelConf, fwdModelIn, fwdModelExtra, &
         & fmStat, phiTan, grids_mag )

    nullify ( cloudIce )
    if ( FwdModelConf%incl_cld .or. FwdModelConf%useTScat ) then
      cloudIce => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra,  &
        & quantityType=l_vmr, molecule=l_cloudIce, config=fwdModelConf,    &
        & noError=.true. )
    end if
    if ( associated(cloudIce) ) then
      call load_one_item_grid ( grids_IWC, cloudIce, fmStat%maf, phitan, &
        & fwdModelConf, .false., .false. )
    else
      call emptyGrids_t ( grids_IWC ) ! Allocate components with zero size
    end if

    ! Compute some sizes
    maxVert = (nlvl-1) * ngp1 + 1
    if ( usingQTM .and. .not. zeta_only ) then
      ! We don't know the maximum diameter of the QTM HGrid, so assume a ray
      ! could cut through all facets at all altitudes.
      max_c = 2*max(nlvl,QTM_hGrid%QTM_tree%n_facets) + max_new + 1
      max_f = ( max_c - 1 ) * ngp1 + 1
    else ! not QTM, or QTM but only with constant-zeta surface intersections
      max_c = 2*nlvl + max_new + 1 ! Maximum coarse path length
      max_f = 2 * maxVert + ngp1   ! Maximum fine path, including minimum Zeta panel
    end if

    ! Allocate and fill spectroscopy derivative grids.  They'll be empty
    ! if fwdModelConf%line* has size zero.
    call load_sps_data ( FwdModelConf, phitan, fmStat%maf, grids_v, &
      & qtyStuffIn=fwdModelConf%lineCenter%qty )
    call load_sps_data ( FwdModelConf, phitan, fmStat%maf, grids_w, &
      & qtyStuffIn=fwdModelConf%lineWidth%qty )
    call load_sps_data ( FwdModelConf, phitan, fmStat%maf, grids_n, &
      & qtyStuffIn=fwdModelConf%lineWidth_TDep%qty )

    if ( switchDetail(switches,'Grids') > -1 .and. .not. usingQTM ) then
      ! dump the grids
      call dump ( grids_f, "Grids_f", details=9 )
      call dump ( grids_tmp, "Grids_tmp", details=9 )
      if ( size(fwdModelConf%lineCenter) > 0 ) call dump ( grids_v )
      if ( size(fwdModelConf%lineWidth) > 0 ) call dump ( grids_w )
      if ( size(fwdModelConf%lineWidth_TDep) > 0 ) call dump ( grids_n )
    end if

    s_a = merge(1,0,atmos_der)
    s_h = merge(1,0,atmos_second_der)
    s_i = merge(1,0,FwdModelConf%incl_cld)
    s_lc = merge(1,0,spect_der_center)
    s_lw = merge(1,0,spect_der_width)
    s_pfa = merge(1,0,FwdModelConf%anyPFA(1) .or. FwdModelConf%anyPFA(2))
    s_p = merge(1,0,FwdModelConf%polarized)
    s_QTM = merge(1,0,usingQTM)
    s_t = merge(1,0,fwdModelConf%temp%derivOK)
    s_td = merge(1,0,spect_der_width_TDep)
    s_tg = merge(1,0,FwdModelConf%GenerateTScat)
    s_ts = MERGE(1,0,FwdModelConf%useTScat)
    s_m = MERGE(1,0,FwdModelConf%hmag_der .or. FwdModelConf%htheta_der .or. FwdModelConf%hphi_der)

    dump_conf = switchDetail(switches,'fmconf')
    if ( dump_conf > -1 ) call dump( FwdModelConf, details=dump_conf  )

    call FullForwardModelAuto ( FwdModelConf, FwdModelIn, FwdModelExtra,       &
                              & FwdModelOut, FmStat, z_psig, tan_press,        &
                              & grids_tmp, grids_f, grids_mag, grids_iwc,      &
                              & grids_n, grids_v, grids_w,                     &
                              & ptan, phitan, temp,                            &
!Q                              & QTM_HGrid, &
                              & Q_LOS, QTM_Paths, ScECR_MIF,                   &
                              & F_and_V, F_and_V_MIF, no_mol, noUsedChannels,  &
                              & no_sv_p_T, n_t_zeta, sv_t_len, &
                              & no_sv_p_h, n_h_zeta, sv_h_len, &
                              & nlvl, no_tan_hts,&
                              & surfaceTangentIndex,                           &
                              & max_c, maxVert, max_f, ptan_der,               &
                              & s_t, s_a, s_h, s_lc, s_lw, s_td, s_p, s_pfa,   &
                              & s_i, s_tg, s_ts, s_m, s_QTM,                        &
                              ! Optional:
                              & Jacobian, ExtraJacobian, Hessian )

    call destroygrids_t ( grids_f )
    call destroygrids_t ( grids_iwc )
    call destroygrids_t ( grids_mag )
    call destroygrids_t ( grids_tmp )
    call destroygrids_t ( grids_n )
    call destroygrids_t ( grids_w )
    call destroygrids_t ( grids_v )

    ! Allocated in tangent_pressures:
    call deallocate_test ( tan_press,    'tan_press',    moduleName )

    deallocate ( QTM_Paths, stat=stat, errmsg=ermsg )
    call test_deallocate ( stat, moduleName, "QTM_Paths", ermsg=ermsg )
    deallocate ( F_and_V, stat=stat, errmsg=ermsg )
    call test_deallocate ( stat, moduleName, "F_and_V", ermsg=ermsg )
    if ( allocated(F_and_V_MIF) ) then
      deallocate ( F_and_V_MIF, stat=stat, errmsg=ermsg )
      call test_deallocate ( stat, moduleName, "F_and_V_MIF", ermsg=ermsg )
    end if

    call trace_end ( 'FullForwardModel, MAF=', fmStat%maf, cond=toggle(emit) )

  contains

    integer function DerivativeMissingFromStateFun() result( SEVERITY )
      ! Default severity (probably ERROR) can be reduced to a warning 
      ! by adding -Sdmiss to command line
      ! Args
      severity = merge( min( MLSMSG_Warning, DerivativeMissingFromState ), &
                      & DerivativeMissingFromState, &
                      & switchDetail(switches,'dmiss') > -1 )
    end function DerivativeMissingFromStateFun

  end subroutine FullForwardModel

  ! ---------------------------------------- FullForwardModelAuto  -----

  subroutine FullForwardModelAuto ( FwdModelConf, FwdModelIn, FwdModelExtra,   &
                             & FwdModelOut, FmStat, z_psig, tan_press,         &
                             & grids_tmp,  grids_f, grids_mag, grids_iwc,      &
                             & grids_n, grids_v, grids_w, ptan,                &
                             & phitan, temp,                                   &
!Q                             & QTM_HGrid, &
                             & Q_LOS, QTM_Paths, ScECR_MIF,                    &
                             & F_and_V, F_and_V_MIF, no_mol, noUsedChannels,   &
                             & No_sv_p_T, n_t_zeta, sv_t_len, &
                             & no_sv_p_h, n_h_zeta, sv_h_len, &
                             & nlvl, no_tan_hts,&
                             & surfaceTangentIndex,                            &
                             & max_c, maxVert, max_f, ptan_der,                &
                             & s_t, s_a, s_h, s_lc, s_lw, s_td, s_p, s_pfa,    &
                             & s_i, s_tg, s_ts, s_m, s_QTM,                         &
                             ! Optional:
                             & Jacobian, ExtraJacobian, Hessian )

  ! This is the full radiative transfer forward model, the workhorse code.
  ! It's called FullForwardModelAuto because most of the variable-size
  ! work arrays are automatic arrays, instead of being explicitly allocated
  ! and deallocated.  Their extents are sufficient for the longest path
  ! through the atmosphere on which the radiative-transfer equation is
  ! integrated.  If they were allocated in One_Pointing with an extent
  ! appropriate to each path, much more time would be spent in the storage
  ! manager.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Comp_Sps_Path_Sparse_m, only: Comp_Sps_Path_Sparse
    use Compute_GL_Grid_M, only: Compute_GL_Grid
    use Constants, only: Deg2Rad, Rad2Deg
    use Convolution_m, only: Convolution, Convolution_Setup
    use Dump_0, only: Dump
    use FilterShapes_M, only: DACSFilterShapes, FilterShapes
    use ForwardModelConfig, only: Beta_Group_T, Channels_T, &
      & ForwardModelConfig_T, LineCenter, LineWidth, LineWidth_TDep
    use ForwardModelIntermediate, only: ForwardModelStatus_T, B_Refraction
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_V_Geod
    use Geometry, only: Get_R_EQ
    use Get_IEEE_NaN_m, only: Fill_IEEE_NaN
    use HessianModule_1, only: Hessian_T
    use HighOutput, only: HeadLine, OutputNamedValue
    use Intrinsic, only: L_A, L_Radiance, L_TScat, L_VMR
    use Load_SPS_Data_M, only: DestroyGrids_T, Dump, Grids_T
    use MatrixModule_1, only: Matrix_T
    use Metrics_3D_m, only: S_QTM_t
    use MLSKinds, only: R4, R8, RP, RV
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSNumerics, only: Hunt, InterpolateValues, PureHunt
    use MLSSignals_M, only: AreSignalsSuperset, Dump, &
      & GetNameOfSignal, MatchSignal, &
      & Radiometers, Signal_T
    use MLSStringLists, only: SwitchDetail
    use Molecules, only: L_H2O, L_N2O, L_O3
    use MoreMessage, only: MLSMessage
    use Output_M, only: NewLine, Output
    use Path_Contrib_M, only: Get_GL_Inds
    use Path_Representation_m, only: Facets_and_Vertices_t, Path_t
    use Physics, only: SpeedOfLight
    use PointingGrid_M, only: PointingGrids, PointingGrid_T, &
      & Dump_Pointing_Grid
    use Radius_of_Curvature_m, only: Radius_of_Curvature_Normal
    use Read_Mie_m, only: IWC_S, T_S
    use Slabs_sw_m, only: AllocateSLABS, DestroyCompleteSLABS, SLABS_Struct
    use Sparse_Eta_m, only: Sparse_Eta_t
    use Sparse_m, only: Sparse_t
    use Tau_m, only: Destroy_Tau, Dump, Tau_t
    use Toggles, only: Emit, Levels, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use TScat_Support_M, only: Print_TScat_Deriv, Print_TScat_Detail, &
      & TScat_Detail_Heading, TScat_Gen_Setup
    use Two_D_Hydrostatic_M, only: Two_D_Hydrostatic
    use VectorsModule, only: CloneVectorQuantity, DestroyVectorQuantityValue, &
      & GetVectorQuantityByType, Vector_T, VectorValue_T

    type(forwardModelConfig_T), intent(inout) :: FwdModelConf
    type(vector_T), intent(in) :: FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    real(rp), intent(in) :: Z_PSIG(:)       ! Surfs from Temperature, tangent
                                            ! grid and species grids, sans
                                            ! duplicates.
    real(rp), intent(in) :: Tan_Press(:)    ! Pressures corresponding to Z_PSIG
    type (Grids_T), intent(inout) :: Grids_tmp ! All the coordinates for TEMP,
                                            ! inout for cloud model to change
                                            ! supersat
    type (Grids_T), intent(inout) :: Grids_f   ! All the coordinates for VMR,
                                            ! inout for cloud model to change
                                            ! supersat
    type (Grids_T), intent(in) :: Grids_IWC ! All the coordinates for IWC
    type (Grids_T), intent(in) :: Grids_Mag ! All the coordinates for Magnetic field
    type (Grids_T), intent(in) :: Grids_n   ! All the spectroscopy(N) coordinates
    type (Grids_T), intent(in) :: Grids_v   ! All the spectroscopy(V) coordinates
    type (Grids_T), intent(in) :: Grids_w   ! All the spectroscopy(W) coordinates
    type (VectorValue_T), pointer :: PhiTan ! Tangent geodAngle component of
                                            ! state vector (minor frame quantity)
    type (VectorValue_T), intent(in) :: PTan ! Tangent pressure component of
                                            ! state vector (minor frame quantity)
    type (VectorValue_T), intent(in) :: Temp ! Temperature component of state vector
    type (Path_t), allocatable :: Q_LOS(:)  ! Minor-frame lines of sight for QTM
    type (Path_t), intent(inout) :: QTM_Paths(:) ! LOS through QTM
    type (VectorValue_T), intent(in), pointer :: ScECR_MIF ! Instrument
                                            ! position in ECR (MIF quantity)
    type (Facets_and_Vertices_t), intent(in) :: F_and_V(:) ! of QTM under paths
                                            ! used to integrate the radiative-transfer
                                            ! equation
    type (Facets_and_Vertices_t), intent(in), allocatable :: F_and_V_MIF(:) ! of
                                            ! QTM under MIF paths.  Used for
                                            ! hydrostatic equilibrium
                                            ! calculation that is used to
                                            ! calculate MIF-based Chi angles.
    type(QTM_tree_t) :: QTM ! for QTM case; only
    integer, intent(in) :: No_Mol           ! Number of molecules
    integer, intent(in) :: NoUsedChannels   ! Number of channels used
    integer, intent(in) :: No_sv_p_T        ! number of phi basis for temperature
    integer, intent(in) :: N_T_zeta         ! Number of zetas for temperature
    integer, intent(in) :: Sv_T_len         ! Number of t_phi*t_zeta in the window
    integer, intent(in) :: No_sv_p_h        ! number of phi basis for magnetic field
    integer, intent(in) :: N_h_zeta         ! Number of zetas for magnetic field
    integer, intent(in) :: Sv_h_len         ! Number of h_phi*h_zeta in the window
    integer, intent(in) :: Nlvl             ! Number of levels in coarse zeta grid
    integer, intent(in) :: No_Tan_Hts       ! Number of tangent heights
    integer, intent(in) :: SurfaceTangentIndex ! Index in tangent grid of
                                            ! earth's surface
    integer, intent(in) :: Max_C            ! Length of longest possible coarse path
    integer, intent(in) :: MaxVert          ! Number of points in gl-refined vertical grid:
                                            ! (NLVL-1) * (NG+1) + 1, i.e., 1 + NG
                                            ! per level, except the last,
                                            ! where there's no GL space.
    integer, intent(in) :: Max_F            ! Length of longest possible path (all npf<max_f)
    logical, intent(in) :: PTan_Der
    integer, intent(in) :: S_A  ! Multiplier for atmos derivative sizes, 0 or 1
    integer, intent(in) :: S_H  ! Multiplier for atmos second derivative sizes, 0 or 1
    integer, intent(in) :: S_I  ! Multiplier for ice/cloud sizes, 0 or 1
    integer, intent(in) :: S_LC ! Multiplier for line center deriv sizes, 0 or 1
    integer, intent(in) :: S_LW ! Multiplier for line width deriv sizes, 0 or 1
    integer, intent(in) :: S_P  ! Multiplier for polarized sizes, 0 or 1
    integer, intent(in) :: S_PFA ! Multiplier for PFA sizes, 0 or 1
    integer, intent(in) :: S_QTM ! Multiplier for sizes of QTM-related stuff, 0 or 1
    integer, intent(in) :: S_T  ! Multiplier for temp derivative sizes, 0 or 1
    integer, intent(in) :: S_TD ! Multiplier for temp dependence deriv sizes, 0 or 1
    integer, intent(in) :: S_TG ! Multiplier for TScat generation sizes, 0 or 1
    integer, intent(in) :: S_TS ! Multiplier for using TScat tables, 0 or 1
    integer, intent(in) :: S_M  ! Multiplier for magnetic field, 0 or 1

    type(matrix_T),  intent(inout), optional :: Jacobian
    type(matrix_T),  intent(inout), optional :: ExtraJacobian ! This is used
                                ! if a quantity with respect to which we are
                                ! computing derivatives is found in the extra
                                ! vector.  Usually, this means we are retrieving
                                ! a MIF quantity, with profile derivatives
                                ! computed here.
    type(hessian_T), intent(inout), optional :: Hessian

    ! Now define local variables, group by type and then
    ! Alphabetically

    integer :: Channel            ! A Loop counter
    integer :: Dump_Rad_Pol       ! Dump intermediate values in Rad_Tran_Pol
    integer :: Frq_Avg_Sel        ! Summarizes combinations of PFA, LBL,
                                  ! Frequency averaging and derivatives.
                                  ! See Frequency_Average below.
    integer :: H2O_Ind            ! Grids_f%S_Ind(L_H2O) = index of H2O in Grids_f
    integer :: IER                ! Status flag from allocates
    integer :: I                  ! Loop index and other uses.
    integer :: Inst               ! Instance of temperature nearest to MAF
    integer :: J                  ! Loop index and other uses.
    integer :: K                  ! Loop index and other uses.
    integer :: MAF                ! MAF under consideration
    integer :: Me = -1            ! String index for trace
    integer :: Me_Hydro = -1      ! String index for trace
    integer :: Me_PointingLoop = -1 ! String index for trace
    integer :: Me_Sideband = -1   ! String index for trace
    integer :: Me_SidebandLoop = -1 ! String index for trace
    integer :: MIF                ! MIF number for tan_press(ptg_i) and rotation
                                  ! matrix from ECR to IFOVPP (R1A) for magnetic
                                  ! field
    integer :: N_More             ! Effective size of More_*_Path
    integer :: N_Phi              ! Number of phi values for one species
    integer :: NGLMax             ! NGL if all panels need GL
    integer :: NoFreqs            ! Number of frequencies for a pointing
    integer :: NoUsedDACS         ! Number of different DACS in this run.
    integer :: NPC                ! Length of coarse path
    integer :: NPF                ! Length of a combined coarse & fine path
    integer :: NZ_IF              ! Effective size of Z_GLgrid and cohorts
    integer :: NZ_IG              ! Effective size of Zetas, size(z_psig) <=
                                  ! nz_ig <= size(z_psig)+2
    integer :: Prev_Mie_Frq_Ind   ! From a previous execution of One_Frequency
    integer :: PTG_I              ! Loop counter for the pointings
    integer :: Sideband           ! Either zero or from firstSignal
    integer :: SX                 ! 1 = LSB, 2 = USB
    integer :: Tan_Ind_C          ! Index of tangent point in coarse zeta ref grid
    integer :: Tan_Ind_F          ! Index of tangent point in fine zeta ref grid
    integer :: Tan_Pt_C           ! Index of tangent point in coarse path
    integer :: Tan_Pt_F           ! Index of tangent point in fine path
    integer :: ThisSideband       ! Loop counter for sidebands, -1 = LSB, +1 = USB
    integer :: WindowFinish       ! End of temperature `window'
    integer :: WindowStart        ! Start of temperature `window'

    integer :: Nspec              ! No of species for cloud model

    logical :: Any_Der            ! temp_der .or. atmos_der .or. spect_der .or. magnetic_der
    ! logical :: Clean              ! Used for dumping
    character(len=4) :: Clean     ! Used for dumping
    logical :: Do_More_Points     ! Do intersections of path at zetas < zeta(tan)
    logical :: Do_Zmin            ! "Do minimum Zeta calculation"
    logical :: Dump_Tscat         ! For debugging, from -Sdsct
    logical, parameter :: PFAFalse = .false.
    logical, parameter :: PFATrue = .true.
    logical :: Print_Frq          ! For debugging, from -Sffrq
    integer :: Print_Deltau_Pol   ! For debugging, from -Sdeltaupol
    logical :: Print_Incopt       ! For debugging, from -Sincp
    logical :: Print_IncRad       ! For debugging, from -Sincr
    integer :: Print_Mag          ! For debugging, from -Smag#
    logical :: Print_More_Points  ! Print if Do_More_Points finds more, from -SZMOR
    integer :: Print_Path         ! Nicer format than Print_Incopt, for few molecules
    integer :: Print_Pol_Rad      ! For debugging, from -Spolr[n]
                                  ! n>0 => print pointing number
                                  ! n>1 => print frequency
    integer :: Print_Rad          ! For debugging, from -Srad
    logical :: Print_Seez         ! For debugging, from -Sseez
    logical :: Print_TauL         ! For debugging, from -Staul
    logical :: Print_TauP         ! For debugging, from -Staup
    logical :: Print_TScat        ! For debugging, from -STScat

    logical :: UsingQTM

    logical :: temp_der, atmos_der, spect_der ! Flags for various derivatives
    logical :: hmag_der, htheta_der, hphi_der ! Flags for various derivatives
    logical :: atmos_second_der   ! Flag for atmos second derivatives
    logical :: Spect_Der_Center, Spect_Der_Width, Spect_Der_Width_TDep

    character (len=127) :: SigName      ! Name of a Signal

    integer, target :: C_Inds_B(max_c)  ! Base array for C_INDS
    integer, target :: CG_Inds_B(max_c) ! Base array for CG_INDS
    integer :: F_Inds(max_c*ng)         ! Indices of fine grid within combined
                                        ! coarse & fine path.
    integer, target :: GL_Inds_B(max_c*ng) ! Base array for GL_INDS
    integer :: Grids(no_tan_hts)        ! Indices in ptgGrid for each tangent
    integer :: IPSD(s_i*max_f)
    integer :: Vert_Inds(max_f)         ! Height indices of fine path in
                                        ! H_Glgrid etc.
    integer :: Where_GL(max_c)          ! If K = Where_GL(i) /= 0 it indicates
                                        ! GL_Inds(k:k+ng-1) are GL_Inds
                                        ! for panel(i:i+1) in the coarse
                                        ! path.

    real(r8) :: WC(s_i*fwdModelConf%no_cloud_species, max_f)
    real(r8) :: Scat_ang(s_i*fwdModelConf%num_scattering_angles)

    real(rp), target :: Alpha_Path(max_f) ! Absorption coefficients on composite
                                      ! coarse & fine path
    real(rp), pointer :: Alpha_Path_C(:)  ! coarse grid part of Alpha_Path
    real(rp) :: Alpha_Path_F(max_c*ng)    ! absorption coefficients for GL
    real(rp) :: B(max_c)              ! Planck radiation function
    real(rp) :: Beta_Path_Cloud_c(s_i*max_c) ! Beta on path coarse
    real(rp) :: Beta_c_e_Path_c(s_ts*max_c)  ! Beta_c_e on coarse path
    real(rp) :: Beta_c_s_Path_c(s_ts*max_c)  ! Beta_c_s on coarse path
    real(r8), target :: ChannelCenters(noUsedChannels) ! for PFA or non-frequency-averaging
    real(rp) :: DAlpha_DT_Path(max_f)   ! dAlpha/dT on coarse & fine grid
    real(rp) :: DAlpha_DT_Path_F(max_f) ! dAlpha/dT on fine grid
    real(rp) :: Del_S(max_c)          ! Integration lengths along coarse path,
                                      !  (2:npc-1).  Before the tangent point,
                                      !  Del_s(i) is the path length [i-1:i].
                                      !  After the tangent point, Del_s(i) is
                                      !  the path length [i:i+1].  Thereby,
                                      !  IncOptDepth(i) is \int_{i-1}^i Alpha
                                      !  before the tangent point, and
                                      !  \int_i^{i+1} Alpha after the tangent
                                      !  point.
    real(rp) :: Del_Zeta(max_c)       ! Integration lengths in Zeta coords
                                      !  along coarse path, (2:npc-1).  Same
                                      !  indexing as for Del_S.
    real(rp) :: dBeta_c_a_dIWC_Path_C(s_ts*max_c)  ! on coarse path
    real(rp) :: dBeta_c_s_dIWC_Path_C(s_ts*max_c)  ! on coarse path
    real(rp) :: dBeta_c_a_dT_Path_C(s_ts*max_c)    ! on coarse path
    real(rp) :: dBeta_c_s_dT_Path_C(s_ts*max_c)    ! on coarse path
    real(rp) :: DHDZ_Path(max_f)      ! dH/dZ on fine path (1:npf)
    real(rp) :: DHDZ_GW_Path(max_f)   ! dH/dZ * GW on fine path (1:npf)
    real(rp) :: DSDH_Path(max_f)      ! dS/dH on fine path (1:tan_pt_f-1,tan_pt_f+ngp1+1:npf)
    real(rp) :: DSDZ_C(merge(max_c,0,WrongTrapezoidal)) ! ds/dH * dH/dZ on coarse path (1:tan_pt_c-1,tan_pt_c+2:npc)
    real(rp) :: DSDZ_GW_Path(max_f)   ! ds/dH * dH/dZ * GW on path
    real(rp) :: DTanh_DT_C(max_c)     ! 1/tanh1_c d/dT tanh1_c
    real(rp) :: DTanh_DT_F(max_f)     ! 1/tanh1_f d/dT tanh1_f
    real(rp) :: dTScat_df(s_ts*max_c,size(grids_f%values))   ! on coarse path w.r.t. SVE
    real(rp) :: dTScat_dT(s_ts*max_c,sv_t_len) ! on coarse path w.r.t. SVE
    real(rp) :: H_Path(max_f)         ! Heights on path (km)
    real(rp) :: H_Path_C(max_c)       ! H_Path on coarse grid (km)
    real(rp) :: IncOptDepth(max_c)    ! Incremental Optical depth on coarse grid
    real(rp) :: More_H_Path(max_new), More_Phi_Path(max_new), More_Z_Path(max_new)
    real(rp) :: N_Path_C(max_c)       ! Refractive index - 1 on coarse path
    real(rp) :: N_Path_F(max_f)       ! Refractive index - 1 on fine path
    real(rp) :: N_Path_Inst           ! Refractive index - 1 at the instrument,
                                      ! 0.0 for satellite, n_path_c(1) for QTM
    real(rp), target :: Phi_Path(max_f) ! Phi's on fine path, Radians
    real(rp), pointer :: Phi_Path_c(:)  ! Phi on the coarse path
    real(rp) :: P_Path(max_f)         ! Pressure on path
    real(rp) :: Ptg_Angles(no_tan_hts)
    real(rp) :: Ref_Corr(max_c)       ! Refraction correction
    real(rp) :: Tan_dH_dT(s_t*sv_t_len) ! dH/dT at Tangent, used for convolution
    real(rp) :: Tanh1_C(max_c)        ! tanh(0.5 h nu / k T)
    real(rp) :: Tanh1_F(max_f)        ! tanh1 on fine grid
    real(rp) :: T_Path(max_f)         ! Temperatures on entire path
    real(rp) :: T_Path_c(max_c)       ! T_Path on coarse grid
    real(rp) :: T_Path_f(max_f)       ! T_Path on GL points
    real(rp) :: TT_Path_c(max(s_i,s_ts)*max_c)  ! TScat on path coarse
    real(rp) :: W0_Path_c(max(s_i,s_ts)*max_c)  ! w0 on path coarse
    real(rp) :: Z_Coarse(max_c)       ! Z_PSIG & Z_min & surface zeta on path
    real(rp) :: Z_GLGrid(maxvert)     ! Zeta on initial glGrid surfs
    real(rp) :: Z_Path(max_f)         ! Zeta on fine grid path tangent grid and
                                      ! species grids, sans duplicates.
    real(rp) :: Beta_Path_C(max_c,no_mol)    ! on coarse path
    real(rp) :: Beta_Path_F(max_c*ng,no_mol) ! for GL points of path only
    real(rp) :: d2_Delta_dF2(max_c,size(grids_f%values),s_h*size(grids_f%values))
                                      ! Incremental opacity second derivative.
                                      ! Path x SVE x SVE.
    real(rp) :: d_T_SCR_dT(max_c,s_t*sv_t_len)     ! D Delta_B in some notes
                                      ! coarse path x state-vector-components
    real(rp) :: d2X_dXdT(no_tan_hts,s_t*sv_t_len)    ! (No_tan_hts, nz*np)
    real(rp), target :: dAlpha_dF_Path(max_f,s_a*no_mol) ! on composite coarse &
                                                     ! fine path
    real(rp), pointer :: dAlpha_DF_Path_C(:,:)       ! on coarse path
    real(rp) :: dAlpha_dF_Path_F(max_f,s_a*no_mol)   ! GL points only
    real(rp) :: d2Alpha_dF2_Path_C(max_c,s_h*no_mol) ! on coarse path
    real(rp) :: d2Alpha_dF2_Path_F(max_f,s_h*no_mol) ! on GL path
    real(rp) :: dB_dF(s_ts*max(s_a,s_t)*max_c)       ! dB / d one f on the path, for TScat
    real(rp) :: dBeta_dF_Path_c(max_c,count(grids_f%where_dBeta_df /= 0))
    real(rp) :: dBeta_dF_Path_f(max_f,count(grids_f%where_dBeta_df /= 0))
    real(rp) :: dBeta_dIWC_Path_c(max_c,s_tg*no_mol) ! dBeta_dIWC on coarse grid
    real(rp) :: dBeta_dIWC_Path_f(max_f,s_tg*no_mol) ! dBeta_dIWC on fine grid
    real(rp) :: dBeta_dN_Path_c(max_c,s_td*size(fwdModelConf%lineWidth_TDep)) ! dBeta_dn on coarse grid
    real(rp) :: dBeta_dN_Path_f(max_f,s_td*size(fwdModelConf%lineWidth_TDep)) ! dBeta_dn on fine grid
    real(rp) :: dBeta_dT_Path_c(max_c,s_t*no_mol)  ! dBeta_dT on coarse grid
    real(rp) :: dBeta_dT_Path_f(max_f,s_t*no_mol)  ! dBeta_dT on fine grid
    real(rp) :: dBeta_dV_Path_c(max_c,s_lc*size(fwdModelConf%lineCenter)) ! dBeta_dv on coarse grid
    real(rp) :: dBeta_dV_Path_f(max_f,s_lc*size(fwdModelConf%lineCenter)) ! dBeta_dv on fine grid
    real(rp) :: dBeta_dW_Path_c(max_c,s_lw*size(fwdModelConf%lineWidth)) ! dBeta_dw on coarse grid
    real(rp) :: dBeta_dW_Path_f(max_f,s_lw*size(fwdModelConf%lineWidth)) ! dBeta_dw on fine grid
    type(sparse_t) :: dH_dT_Path               ! dH/dT on path X (Zeta * Phi)
    real(rp) :: dHdZ_GLGrid(maxVert,no_sv_p_T) ! dH/dZ on glGrid surfs
    real(rp) :: dX_dT(no_tan_hts,s_t*sv_t_len) ! (No_tan_hts, nz*np)
    type (Sparse_Eta_t) :: Eta_ZP_N(size(grids_n%values)) ! Eta_z x Eta_p for N
    type (Sparse_Eta_t) :: Eta_ZP_V(size(grids_v%values)) ! Eta_z x Eta_p for V
    type (Sparse_Eta_t) :: Eta_ZP_W(size(grids_w%values)) ! Eta_z x Eta_p for W
    real(rp) :: H_GLGrid(maxVert,no_sv_p_T)    ! H on glGrid surfs (km)
    real(rp) :: IWC_Path(max_f,max(max(s_i,s_ts),s_ts)) ! IWC on path
    real(rp), target :: Mag_Path(s_p*max_f,4)  ! Magnetic field on path
    real(rp), target :: Rad_Avg_Path(max_c,s_pfa*noUsedChannels) ! Freq. Avgd.
                                               ! LBL radiance along the path
    real(rp) :: R_Eq  ! Radius of equivalent circular Earth tangent to the
                      ! surface of the Earth reference ellipsoid at the surface
                      ! location of the tangent point, in the plane defined by
                      ! the line of sight and the center of the Earth.
    real(rp) :: Radiances(noUsedChannels,no_tan_hts) ! (noChans,Nptg)
!     real(rp) :: RADIANCES_DIFF(noUsedChannels,no_tan_hts) ! (noChans,Nptg)     ! IGOR
    real(rp) :: Spect_N_Path(max_f,size(fwdModelConf%lineWidth_TDep)) ! Line Width Temperature Dependence
    real(rp) :: Spect_V_Path(max_f,size(fwdModelConf%lineCenter)) ! Line Center
    real(rp) :: Spect_W_Path(max_f,size(fwdModelConf%lineWidth)) ! Line Width
    real(rp) :: SPS_Path(max_f,no_mol)         ! species on path
    real(rp) :: SPS_Path_C(max_c,no_mol)       ! species on coarse path
    real(rp) :: SPS_Path_F(max_f,no_mol)       ! species on GL path
    real(rp) :: TT_Path_F(max_f,max(s_i,s_ts)) ! TScat on entire path
    real(rp) :: T_GLGrid(maxVert,no_sv_p_T)    ! Temp on glGrid surfs
    real(rp) :: T_Script_PFA(max_c,s_pfa*noUsedChannels) ! Delta_B in some notes
    real(r8) :: VMRArray(s_i*n_t_zeta,no_mol)  ! The VMRs for H2O, O3, N2

    ! Temporary space for DACS radiances if we're doing frequency averaging
    ! and there are any LBL molecules and any derivatives are calculated
    real(rp) :: DACsStaging2(lbound(fwdModelConf%DACsStaging,1): &
      &                      ubound(fwdModelConf%DACsStaging,1), &
      & merge(1,0,fwdModelConf%do_freq_avg .and. &
      &           any(fwdModelConf%anyLBL((fwdModelConf%sidebandStart+3)/2:    &
      &                                   (fwdModelConf%sidebandStop+3)/2))) * &
      &   MAX(s_t,s_a,s_m,s_lc,s_lw,s_td) *      & ! merge(1,0,any_der)
      &   MAX(sv_t_len,sv_h_len, SIZE(grids_f%values), &
      &       size(grids_w%values), size(grids_n%values), size(grids_v%values)), &
      & size(fwdModelConf%usedDACSSignals) )

    real(rp) :: DACsStaging3( lbound(DACsStaging2,1):ubound(DACsStaging2,1), &
                            & lbound(DACsStaging2,2):ubound(DACsStaging2,2), &
                            & lbound(DACsStaging2,2):ubound(DACsStaging2,2), &
                            & size(fwdModelConf%usedDACSSignals)                )

    real(rp), target :: DH_DT_GLGrid(maxVert,n_t_zeta,s_t*no_sv_p_T)

    complex(rp) :: D_Rad_Pol_dF(2,2,s_p*s_a*size(grids_f%values)) ! From mcrt_der
    complex(rp) :: D_Rad_Pol_dT(2,2,s_p*s_t*sv_t_len) ! From mcrt_der
    complex(rp) :: D_Rad_Pol_dH(2,2,s_p*s_m*sv_h_len) ! From mcrt_der
    complex(rp) :: D_Rad_Pol_dTheta(2,2,s_p*s_m*sv_h_len) ! From mcrt_der
    complex(rp) :: D_Rad_Pol_dPhi(2,2,s_p*s_m*sv_h_len) ! From mcrt_der
    complex(rp), target :: Alpha_Path_Polarized(-1:1,s_p*max_f)
    complex(rp), pointer :: Alpha_Path_Polarized_C(:,:)
    complex(rp) :: Alpha_Path_Polarized_F(-1:1,s_p*max_f)
    complex(rp), target :: dAlpha_dT_Polarized_Path(-1:1,s_p*s_t*max_f)
    complex(rp), pointer :: dAlpha_dT_Polarized_Path_C(:,:)
    complex(rp) :: dAlpha_dT_Polarized_Path_F(-1:1,s_p*s_t*max_f)
    complex(rp), target :: dAlpha_dH_Polarized_Path(-1:1,s_p*s_m*max_f)
    complex(rp), pointer :: dAlpha_dH_Polarized_Path_C(:,:)
    complex(rp) :: dAlpha_dH_Polarized_Path_F(-1:1,s_p*s_m*max_f)
    complex(rp) :: Beta_Path_Polarized(-1:1,s_p*max_c,no_mol)
    complex(rp) :: Beta_Path_Polarized_F(-1:1,s_p*max_f,no_mol)
    complex(rp) :: dBeta_dT_Polarized_Path_C(-1:1,s_p*s_t*max_c,no_mol)
    complex(rp) :: dBeta_dT_Polarized_Path_F(-1:1,s_p*s_t*2*max_f,no_mol)
    complex(rp) :: dBeta_dH_Polarized_Path_C(-1:1,s_p*s_m*max_c,no_mol)
    complex(rp) :: dBeta_dH_Polarized_Path_F(-1:1,s_p*s_m*2*max_f,no_mol)
    complex(rp) :: DE_DF(2,2,s_p*s_a*max_c,size(grids_f%values)) ! DE/Df in Michael's notes
    complex(rp) :: DE_DT(2,2,s_p*s_t*max_c,sv_t_len) ! DE/DT in Michael's notes
    complex(rp) :: DE_DHMAG(2,2,s_p*s_m*max_c,sv_h_len) ! DE/DH in Michael's notes
    complex(rp) :: DE_DTHETA(2,2,s_p*s_m*max_c,sv_h_len) ! DE/DTheta in Michael's notes
    complex(rp) :: DE_DPHI(2,2,s_p*s_m*max_c,sv_h_len) ! DE/DPhi in Michael's notes
    complex(rp) :: DelTau_Pol(2,2,s_p*max_c) ! E in Michael's notes
!   complex(rp) :: dIncOptDepth_Pol_dT(2,2,s_p*s_t*max_c) ! D Incoptdepth_Pol / DT
!   complex(rp) :: GL_Delta_Polarized(-1:1,s_p*max_f)
    complex(rp) :: IncOptDepth_Pol(2,2,s_p*max_c)
    complex(rp) :: IncOptDepth_Pol_test(2,2,s_p*max_c)
    complex(rp) :: Prod_Pol(2,2,s_p*max_c) ! P in Michael's notes
    complex(rp) :: Tau_Pol(2,2,s_p*max_c)  ! Tau in Michael's notes

    real(rp) :: Est_ScGeocAlt(no_tan_hts)  ! Est S/C geocentric altitude /m
    real(rp) :: Est_LOS_Vel(no_tan_hts)    ! Est S/C line-of-sight velocity M/S
    real(rp) :: Tan_D2H_DHDT(s_t*sv_t_len) ! for convolution
    real(rp) :: Tan_Phi(no_tan_hts)        ! Radians
    real(rp) :: DDHIDHIDTL0(maxVert,n_t_zeta,s_t*no_sv_p_T)
    ! Frequency-averaged derivatives
    ! Channels x pointings x grid values == frequencies x surfaces x instances x molecules:
    real(r4) :: K_Atmos(noUsedChannels,no_tan_hts,s_a*size(grids_f%values))
    real(r4) :: H_Atmos(noUsedChannels,no_tan_hts,s_h*size(grids_f%values),s_h*size(grids_f%values))
    ! Channels x pointings x grid values == frequencies x surfaces x instances x molecules:
    real(r4) :: K_Spect_DN(noUsedChannels,no_tan_hts,s_td*size(grids_n%values))
    real(r4) :: K_Spect_DV(noUsedChannels,no_tan_hts,s_lc*size(grids_v%values))
    real(r4) :: K_Spect_DW(noUsedChannels,no_tan_hts,s_lw*size(grids_w%values))
    real(r4) :: K_Temp(noUsedChannels,no_tan_hts,s_t*sv_t_len)
    REAL(r4) :: K_Magfield(noUsedChannels,no_tan_hts,s_m*sv_h_len*3)

    logical :: Do_GL(max_c)       ! GL indicator.  Before the tangent point,
                                  ! Do_GL(i) means the panel (i-1:i) on the
                                  ! coarse path needs GL.  After the tangent
                                  ! point, Do_GL(i) means panel(i:i+1) needs GL.
    logical :: T_der_Path_flags(s_t*max_f) ! a flag that tells where an
      ! absorption coefficient is needed for a temperature derivative.
      ! Only useful when subsetting temperature derivatives.

    integer, pointer :: C_Inds(:)        ! Indices of coarse grid in combined
                                         ! coarse & fine path
    integer, pointer :: LineCenter_IX(:) ! Where are line center offsets?
    integer, pointer :: LineWidth_IX(:)  ! Where are line width offsets?
    integer, pointer :: LineWidth_TDep_IX(:)  ! Where are line width TDep offsets?
    integer, pointer :: UsedDACSsignals(:) ! Indices in FwdModelConf of signals
                                         ! for our dacs

    real(rp) :: E_Rflty       ! Earth reflectivity at given tan. point
    real(rp) :: H_Surf        ! Height above earth surface of first (usually
                              ! zeta=-3) surface
    real(rp) :: Min_Zeta      ! Minimum zeta along the path
    real(rp) :: Min_Phi       ! Phi at which minimum zeta occurs
    real(rp) :: ROT(3,3)      ! ECR-to-FOV rotation matrix
    real(rp) :: Vel_Rel       ! Vel_z / c

    REAL(rp) :: st(max_f)     ! sine of theta
    REAL(rp) :: sp(max_f)     ! sine of phi
    real(rp) :: cp(max_f)     ! cosine of phi
    
    real(rp), pointer :: CT(:)           ! Cos(Theta), where theta
      ! is the angle between the line of sight and magnetic field vectors.
    real(r8), pointer :: Frequencies(:)  ! Frequencies to compute for
    real(rp), pointer :: H(:)            ! Magnetic field on path, in
                                         ! IFOVPP
    real(rp), allocatable :: RadV(:)     ! Radiances for 1 pointing on
                                         ! Freq_Grid
    real(rp), pointer :: STCP(:)         ! Sin(Theta) Cos(Phi) where
      ! theta is as for CT and phi (for this purpose only) is the angle
      ! between the plane defined by the line of sight and the magnetic
      ! field vector, and the "instrument field of view plane polarized"
      ! (IFOVPP) X axis.
    REAL(rp), POINTER :: STSP(:)         ! Sin(Theta) Sin(Phi)
    real(rp), pointer :: DACsStaging(:,:) ! Temporary space for DACS radiances

    ! Cloud arrays have the same grids as temperature
    real(rp) :: Cext_Path(max_f,s_i)    ! Cloud extinction on path
    real(rp) :: Salb_Path(max_f,s_i)    ! Single Scattering Albedo on path
    real(rp) :: Tscat_Path(max_f,s_i*fwdModelConf%num_scattering_angles) ! TScat on path

    ! Beta on path coarse, which is needed to compute d_W0/dVMR
    real(rp), allocatable :: Inc_Rad_Path(:,:)   ! Incremental radiance along the path
    real(rp), allocatable :: K_Atmos_Frq(:,:)    ! dI/dVMR, ptg.frq X vmr-SV
    real(rp), allocatable :: K_Spect_DN_Frq(:,:) ! ****
    real(rp), allocatable :: K_Spect_DV_Frq(:,:) ! ****
    real(rp), allocatable :: K_Spect_DW_Frq(:,:) ! ****
    REAL(rp), ALLOCATABLE :: K_Temp_Frq(:,:)     ! dI/dT, ptg.frq X T-SV
    REAL(rp), ALLOCATABLE :: K_MagField_Frq(:,:) ! dI/dh, ptg.frq X H-SV X 3
    real(rp), allocatable :: T_Script_LBL(:,:)   ! Delta_B in some notes, Path X Frq
    real(rp), allocatable :: H_Atmos_Frq(:,:,:)  ! d2I/(dVMRq dVMRr), ptg.frq X vmr-SVq X vmr-SVr

    ! Used only to schlep from Convolution_Setup to Convolution
    real(rp) :: Chi_MIF_Tan(ptan%template%nosurfs) ! Chi angle to MIF tangent
    real(rp) :: DH_DZ_OUT(ptan%template%nosurfs)
    real(rp) :: DX_DH_OUT(ptan%template%nosurfs)
    real(rp) :: DXDT_Surface(1,s_t*sv_t_len)
    real(rp) :: DXDT_TAN(ptan%template%nosurfs,s_t*sv_t_len)
    real(rv), pointer :: L1BMIF_TAI(:,:)   ! MIF Times
    real(rv), pointer :: MIFDeadTime(:,:)  ! Not collecting data
    real(rp) :: Surf_Angle(1)

!  The 'all_radiometers grid file' approach variables declaration:

    real(rp) :: Max_ch_freq_grid, Min_ch_freq_grid

! *** Beta & Molecules grouping variables:
    type (Beta_Group_T), pointer :: Beta_Group(:) ! from FwdModelConf%Beta_Group

    ! Tangent point coordinates on GL zeta grid.  Used only for QTM.
    type (H_V_Geod) :: Tan_Pt_Geod(max_f)

    ! Interpolation factors for Temperature, IWC, and Temperature X IWC (not
    ! geometric interpolation factors) to put Mie Beta_c_a and Beta_c_s into
    ! path
    type (Sparse_Eta_t) :: Eta_IWC_Path_c, Eta_T_IWC_Path_c, Eta_T_Path_c

! Channel information from the signals database as specified by fwdModelConf
    type (Channels_T), pointer, dimension(:) :: Channels

    type (Grids_T) :: Grids_Tscat ! All the coordinates for scattering source function
    type (Grids_T) :: Grids_Salb  ! All the coordinates for single scattering albedo
    type (Grids_T) :: Grids_Cext  ! All the coordinates for cloud extinction

    type (PointingGrid_T), pointer :: WhichPointingGrid ! Pointing grids for one signal

    type (Signal_T), pointer :: FirstSignal        ! The first signal we're dealing with

    type (Slabs_Struct), allocatable :: GL_SLABS(:,:) ! Freq. independent single-
                                  ! line absorption data for the combined coarse
                                  ! & fine path; subsets selected by or c_inds
                                  ! or GL_inds.
    type (S_QTM_t), allocatable :: S(:)

    type (Tau_T) :: Tau_LBL, Tau_PFA

    ! Mixing-ratio derivatives of incremental opacity, to schlep from
    ! dRad_tran_df to Get_d_Deltau_Pol_df
    type (sparse_t) :: d_delta_df(s_p*size(fwdModelConf%beta_group))
    ! Interpolation coefficients from Freq X Zeta X Horizontal (phi or QTM)
    ! basis to path for all species
    type (Sparse_Eta_t) :: Eta_fzp(size(fwdModelConf%beta_group))
    ! Interpolation from horizontal (phi or QTM) basis to path for each VMR.
    type (Sparse_Eta_t) :: Eta_p(size(fwdModelConf%beta_group))
    ! Interpolation coefficients from horizontal (phi or QTM) basis to path
    ! for temperature
    type (Sparse_Eta_t) :: Eta_p_T
    ! Interpolation coefficients from zeta X horizontal (phi or QTM) basis to
    ! path for temperature; only needed for temperature derivatives
    type (Sparse_Eta_t) :: Eta_zp_T
    ! Interpolation coefficients from zeta X horizontal (phi or QTM) basis to
    ! path for temperature; only needed for temperature derivatives
    ! Interpolation coefficients from zeta X horizontal (phi or QTM) basis to
    ! path for magnetic field; only needed for magnetic field derivatives
    type (Sparse_Eta_t) :: Eta_zp_H
    ! Interpolation coefficients from zeta X phi basis to path for all species
    type (Sparse_Eta_t) :: Eta_zp(size(fwdModelConf%beta_group))
    ! Interpolation coefficients from zeta for each VMR to path GL zeta.
    type (Sparse_Eta_t) :: Eta_z(size(fwdModelConf%beta_group))
    type (Sparse_Eta_t) :: Eta_zzT ! Interpolation coefficients from
                                   ! temperature Zeta to GL zeta.
    type (VectorValue_T), pointer :: BoundaryPressure
    type (VectorValue_T), pointer :: EarthRefl     ! Earth reflectivity
    type (VectorValue_T), pointer :: ECRtoFOV      ! Rotation matrices
    type (VectorValue_T), pointer :: F             ! An arbitrary species
    type (VectorValue_T), pointer :: GPH           ! Geopotential height
    type (VectorValue_T), pointer :: IWC           ! IWC at scattering points, for TScat gen.
    type (VectorValue_T), pointer :: LOSVel        ! Line of sight velocity m/s
    type (VectorValue_T), pointer :: OrbIncline    ! Orbital inclination
    type (VectorValue_T), pointer :: RefGPH        ! Reference geopotential height
    type (VectorValue_T), pointer :: ScatteringAngles ! for TScat computation
    type (VectorValue_T), pointer :: SCGeocAlt     ! S/C geocentric altitude /m
    type (VectorValue_T), pointer :: SpaceRadiance ! Emission from space
    type (VectorValue_T), pointer :: SurfaceHeight ! km above mean sea level
    type (VectorValue_T), pointer :: ThisRadiance  ! A radiance vector quantity

    ! Intermediate minor-frame tangent quantities calculated from LOS, for QTM
    type (VectorValue_T) :: Q_EarthRadC_sq
    type (VectorValue_T), target :: Q_PhiTan
    type (VectorValue_T) :: Q_TanHt

    ! Coarse-zeta grid tangent quantities interpolated from intermediate
    ! minor-frame tangent quantities calculated from LOS, for QTM.
    real(rp) :: EarthRadC_sq(no_tan_hts) ! Square of minor axis of plane-projected
                                         ! ellipse, meters^2
    real(rp) :: TanHt(no_tan_hts)        ! Tangent heights, meters

! Scattering source function for each temperature surface
    type (VectorValue_T) :: scat_src
    type (VectorValue_T) :: scat_alb
    type (VectorValue_T) :: cld_ext

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!                                                          !!!!!
    !!!!! If an array has an extent of No_SV_P_T or SV_T_Len, in   !!!!!
    !!!!! the QTM case it needs to be subscripted using the index  !!!!!
    !!!!! of a vertex adjacent to the path, in the extracted       !!!!!
    !!!!! subset of the entire state vector.  The mapping from QTM !!!!!
    !!!!! indices to path-adjacent indices is QTM%Path_Vertices(). !!!!!
    !!!!! The inverse mapping, to put derivatives into the correct !!!!!
    !!!!! places in the Jacobian, is F_and_V%Vertices.             !!!!!
    !!!!!                                                          !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Executable code --------------------------------------------------------
    ! ------------------------------------------------------------------------

    call trace_begin ( me, 'FullForwardModelAuto, MAF=', index=fmstat%maf, &
      & cond=toggle(emit) )

    alpha_path_polarized_c(-1:,1:) => alpha_path_polarized ( -1:1, 1 :: ngp1 )
    dAlpha_dT_polarized_path_c(-1:,1:) => dAlpha_dT_polarized_path ( -1:1, 1 :: ngp1 )

    dAlpha_dH_polarized_path_c(-1:,1:) => dAlpha_dH_polarized_path ( -1:1, 1 :: ngp1 )

    ! Get pointer to Alpha_Path_C, the coarse-path part of Alpha_Path
    alpha_path_c => alpha_path ( 1::ngp1 )

    ! Get pointer to dAlpha_df_Path_C, the coarse-path part of dAlpha_df_Path
    dAlpha_df_Path_C => dAlpha_df_Path(1::ngp1,:)

    ! Fill REAL local variables with NaN if requested
    if ( NaN_Fill ) then
      ! Scalar RP:
      call fill_IEEE_NaN ( e_rflty, h_surf, min_zeta, min_phi, vel_rel, &
        & max_ch_freq_grid, min_ch_freq_grid )
      ! Rank 1 RP:
      call fill_IEEE_NaN ( Alpha_path_c, Alpha_path_f, b, &
        & beta_path_cloud_c, beta_c_e_path_c, beta_c_s_path_c, &
        & dAlpha_dT_path_f, del_s, del_zeta, dBeta_c_a_dIWC_path_C, &
        & dBeta_c_s_dIWC_path_C, dBeta_c_a_dT_path_C, dBeta_c_s_dT_path_C, &
        & dhdz_path, dhdz_gw_path, dsdh_path, dsdz_c, dsdz_gw_path, dTanh_dT_c, &
        & dTanh_dT_f, h_path, h_path_c, n_path_f, phi_path, &
        & ptg_angles, ref_corr )
      call fill_IEEE_NaN ( tan_dh_dT, tanh1_c, tanh1_f, t_path, t_path_c, &
        & t_path_f, tt_path_c, w0_path_c, z_coarse, z_glgrid, z_path )
      call fill_IEEE_NaN ( earthradc_sq, est_scgeocalt, est_los_vel, &
        & tan_d2h_dhdT, tan_phi )
      call fill_IEEE_NaN ( dh_dz_out, dx_dh_out )
      call fill_IEEE_NaN ( chi_MIF_Tan, dB_df, surf_angle )
      ! Rank 1 R4:
      ! Rank 1 R8:
      call fill_IEEE_NaN ( Scat_ang, channelCenters )
      ! Rank 2 RP:
      call fill_IEEE_NaN ( dTScat_df, dTScat_dT, beta_path_f, &
        & d_t_scr_dT, d2x_dxdT, dAlpha_df_path, dAlpha_df_path_f, &
        & d2Alpha_df2_path_c, d2Alpha_df2_path_f, dBeta_df_path_c, &
        & dBeta_df_path_f, dBeta_dIWC_path_c, dBeta_dIWC_path_f, &
        & dBeta_dn_path_c, dBeta_dn_path_f, dBeta_dT_path_c, dBeta_dT_path_f, &
        & dBeta_dv_path_c, dBeta_dv_path_f, dBeta_dw_path_c, dBeta_dw_path_f )
      call fill_IEEE_NaN ( dhdz_glgrid, dx_dT, &
        & h_glgrid, IWC_path, mag_path, rad_avg_path, radiances, &
        & spect_n_path, spect_v_path, spect_w_path, sps_path, sps_path_c, &
        & sps_path_f )
      call fill_IEEE_NaN ( tt_path_f, t_glgrid, t_script_PFA, rot )
      call fill_IEEE_NaN (  cext_path, salb_path, TScat_path )
      call fill_IEEE_NaN ( dxdT_surface, dxdT_tan )
      ! Rank 2 R8:
      call fill_IEEE_NaN ( wc, vmrArray )
      ! Rank 3 RP:
      call fill_IEEE_NaN ( d2_delta_df2, DACsStaging2, dh_dT_glgrid, &
        & ddhidhidTl0 ) 
      ! Rank 3 R4:
      call fill_IEEE_NaN ( k_atmos, k_spect_dn, k_spect_dv, k_spect_dw, &
        & k_temp, k_magfield )
      ! Rank 3 R8:
      ! Rank 4 RP:
      call fill_IEEE_NaN ( DACsStaging3 )
      ! Rank 4 R4:
      call fill_IEEE_NaN ( h_atmos )
      ! Rank 4 R8:
    end if
    ! Set flags from command-line switches
    clean = ' '
    if ( switchDetail(switches, 'clean') > -1 ) clean = 'c'
    do_zmin = switchDetail(switches, 'dozm') > -1 ! Do minimum zeta only if requested
    dump_rad_pol = 0 ! Dump rad_tran_pol intermediates
    if ( switchDetail(switches, 'dpri') > -1 ) dump_rad_pol = 1 ! Dump if overflow
    if ( switchDetail(switches, 'Dpri') > -1 ) dump_rad_pol = 2 ! Dump all
    if ( switchDetail(switches, 'DPRI') > -1 ) dump_rad_pol = 3 ! Dump first and stop
    dump_rad_pol = switchDetail(switches, 'dpri') + 1
    dump_TScat = switchDetail(switches, 'dsct' ) > -1
    print_Frq = switchDetail(switches, 'ffrq' ) > -1
    print_Deltau_Pol = switchDetail(switches, 'deltaupol' )
    print_Incopt = switchDetail(switches, 'incp' ) > -1
    print_IncRad = switchDetail(switches, 'incr' ) > -1
    print_Mag = switchDetail(switches, 'mag')
    print_Pol_Rad = switchDetail(switches, 'polr')
    print_path = switchDetail(switches, 'path')
    print_Rad = switchDetail(switches, 'rad')
    print_Seez = switchDetail(switches, 'seez') > -1
    print_TauL = switchDetail(switches, 'taul') > -1
    print_TauP = switchDetail(switches, 'taup') > -1
    print_TScat = switchDetail(switches, 'TScat' ) > -1
    Print_TScat_Deriv = switchDetail(switches, 'dtsct' )
    print_TScat_detail = switchDetail(switches, 'psct' )
    print_more_points = switchDetail(switches, 'ZMOR' ) > -1
    do_more_points = switchDetail(switches, 'zmor') > -1

    usingQTM = s_QTM /= 0

    ! Nullify all our pointers that are allocated because the first thing
    ! Allocate_Test does is ask if they're associated.  If we don't nullify
    ! them, they're undefined, i.e., junk that might be mistaken for
    ! associated.

    nullify ( frequencies ) ! Pointer instead of allocatable because it might
                            ! be associated with channel centers instead of
                            ! being allocated

    ! Nullify pointers that are used to control whether calculations get done
    nullify ( linecenter_ix, linewidth_ix, linewidth_tdep_ix )

    d2_delta_df2 = 0.0

    ! Put zeros into H_Atmos so that the elements that aren't filled will
    ! be defined.  It has zero size if second derivatives aren't requested.

    h_atmos = 0.0

    h2o_ind = grids_f%s_ind(l_h2o)

    CALL both_sidebands_setup
    
    PRINT '(a)','=========begin derivative flag settings======'
    PRINT *,'atmos der setting ',atmos_der
    PRINT *,'temperature der setting ',temp_der
    PRINT *,'h-mag der setting ',hmag_der
    PRINT *,'h-theta der setting ',htheta_der
    PRINT *,'h-phi der setting ',hphi_der
    PRINT '(a)','=========end derivative flag settings========'
    
    if ( FwdModelConf%generateTScat ) then
      call TScat_Gen_Setup ( fwdModelConf, fwdModelIn, fwdModelExtra,     &
      &                      fwdModelOut, noUsedChannels, temp, sideband, &
      &                      IWC, scatteringAngles )
    else if ( FwdModelConf%incl_cld ) then
      call cloud_setup
    end if

    ! Compute hydrostatic grid -----------------------------------------------
    ! If this is moved into BOTH_SIDEBANDS_SETUP the run time increases by
    ! a large unexplainable amount, at least in LF95.

    call Trace_Begin ( me_Hydro, 'ForwardModel.Hydrostatic', &
      & cond=toggle(emit) .and. levels(emit) > 0 )

    ! Get interpolation coefficients from temperature's zetas to Z_GLGrid.
    ! Temperature is coherent and stacked, so Eta_zzT is applicable everywhere.

    call eta_zzt%eta_1d ( Grids_tmp%zet_basis, z_glgrid, what=Grids_tmp%qty(1), &
                        & sorted=.true. )

    ! Insert into bill's 2d hydrostatic equation.
    ! The phi inputs for this subprogram are the orbit plane projected
    ! geodetic locations of the temperature phi basis -- not necessarily
    ! the tangent phi's, which might be somewhat different.

    ! For 2D, temperature's windowStart:windowFinish are correct here.
    ! RefGPH and temperature have the same horizontal basis.

    ! For QTM, WindowStart:WindowFinish is the entire grid.  The hydrostatic
    ! calculation is not done here.  Rather, it's done for each pointing.  We
    ! do the calculations only for serial numbers in the path description. 
    ! This entails splitting the Metrics-3D calculation into two parts: First,
    ! identify the facets that the path crosses and gather the vertex serial
    ! numbers.  Then do the hydrostatic calculation for those vertices. Then
    ! finish the vertical part of the Metrics-3D calculation.

    if ( .not. usingQTM ) then ! For QTM, we need a different t_glgrid,
                               ! h_glgrid, etc. for each pointing, so it's
                               ! done later.
!????? Maybe the hydrostatic could be done for QTM here if all the paths cross
!????? the same set of facets.
      if ( temp_der ) then
        call two_d_hydrostatic ( Grids_tmp, refGPH%template%surfs(1,1), &
          &  refGPH%values(1,windowStart:windowFinish), z_glgrid, &
          &  t_glgrid, h_glgrid, dhdz_glgrid, eta_zzT, &
          &  dHidTlm=dh_dt_glgrid, ddHdHdTl0=ddhidhidtl0 )
      else
        call two_d_hydrostatic ( Grids_tmp, refGPH%template%surfs(1,1), &
          &  refGPH%values(1,windowStart:windowFinish), z_glgrid, &
          &  t_glgrid, h_glgrid, dhdz_glgrid, eta_zzT )
      end if
    end if

    call Trace_End ( 'ForwardModel.Hydrostatic', &
      & cond=toggle(emit) .and. levels(emit) > 0 )

    ! Allocate space for the interpolator components for each quantity, and the
    ! components for each frequency-dependent quantity.  We create them here
    ! with enough rows for the longest path.  The number of columns is constant.

    ! Interpolation coefficients for Phi for temperature
    call eta_p_T%create ( max_f, no_sv_p_T, 2*max_f, what=grids_tmp%mol(1) )
    ! Interpolation coefficients for Zeta X Phi for temperature
    IF ( temp_der ) THEN
      call eta_zp_T%create ( max_f, sv_t_len, 2*max_f, what=grids_tmp%mol(1) )
      CALL dh_dt_path%create ( max_f, sv_t_len, 2*max_f, what=grids_tmp%mol(1))
    END IF
    CALL eta_zp_H%create ( max_f, sv_h_len, 2*max_f, what=grids_mag%mol(1) )
    do i = 1, size(beta_group)
      ! For QTM, everything has the same horizontal grid.  The values of
      ! grids_f%l_p are computed from the longest path through the QTM.
      n_phi = grids_f%l_p(i) - grids_f%l_p(i-1)
      ! Interpolation coefficients for Zeta for all VMR
      call eta_z(i)%create ( max_f, grids_f%l_z(i) - grids_f%l_z(i-1), &
                           & 2*max_f, what=grids_f%mol(i) )
      ! Interpolation coefficients for Phi for all VMR
      call eta_p(i)%create ( max_f, n_phi, 2*max_f, what=grids_f%mol(i) )
      ! Interpolation coefficients for Zeta X Phi for all VMR
      call eta_zp(i)%create ( max_f, &
                            & ( grids_f%l_z(i) - grids_f%l_z(i-1) ) * n_phi, &
                            & 2*max_f, what=grids_f%mol(i) )
      ! Interpolation coefficients for Frequency X Zeta X Phi
      call eta_fzp(i)%create ( max_f, &
                             & grids_f%l_v(i) - grids_f%l_v(i-1), &
                             & 2*max_f, what=grids_f%mol(i) )
      ! Representation of d_delta_df for Get_d_Deltau_Pol_df (polarized)
      if ( s_p > 0 ) call d_delta_df(i)%create ( max_f, &
                            & ( grids_f%l_z(i) - grids_f%l_z(i-1) ) * n_phi, &
                            & 2*max_f, what=grids_f%mol(i) )
    end do
    if ( fwdModelConf%useTScat ) then
      ! Create them here instead of in One_Pointing so as not to need to
      ! do it for each pointing.
      if ( .not. allocated(t_s) ) call MLSMessage ( MLSMSG_Error, &
        & moduleName, "UseTScat requested but no Mie tables loaded" )
      call eta_IWC_path_c%create ( max_c, size(IWC_s)*size(T_s), 2*max_c )
      call eta_T_path_c%create ( max_c, size(IWC_s)*size(T_s), 2*max_c )
      call eta_T_IWC_path_c%create ( max_c, size(IWC_s)*size(T_s), 2*max_c )
    end if
    ! Create empty interpolation coefficients for spectroscopy derivatives
    if ( size(fwdModelConf%lineCenter) > 0 ) then
      do i = 1, size(grids_v%values)
        call eta_zp_v(i)%create ( max_f, &
                                & ( grids_v%l_z(i)-grids_v%l_z(i-1) ) * &
                                & ( grids_v%l_p(i)-grids_v%l_p(i-1) ), &
                                & 2*max_f, what=grids_v%mol(i) )
      end do
    end if
    if ( size(fwdModelConf%lineWidth) > 0 ) then
      do i = 1, size(grids_w%values)
        call eta_zp_v(i)%create ( max_f, &
                                & ( grids_w%l_z(i)-grids_w%l_z(i-1) ) * &
                                & ( grids_w%l_p(i)-grids_w%l_p(i-1) ), &
                                & 2*max_f, what=grids_v%mol(i) )
      end do
    end if
    if ( size(fwdModelConf%lineWidth_TDep) > 0 ) then
      do i = 1, size(grids_n%values)
        call eta_zp_v(i)%create ( max_f, &
                                & ( grids_n%l_z(i)-grids_n%l_z(i-1) ) * &
                                & ( grids_n%l_p(i)-grids_n%l_p(i-1) ), &
                                & 2*max_f, what=grids_v%mol(i) )
      end do
    end if

!      PRINT '(a)','t rows'
!      PRINT *,size(eta_zp_t%rows)
!      PRINT *,eta_zp_t%rows
!      PRINT '(a)','t cols'
!      PRINT *,size(eta_zp_t%cols)
!      PRINT *,eta_zp_t%cols
!      PRINT '(a)','t lbnd'
!      PRINT *,eta_zp_t%lbnd
!      PRINT '(a)','t ubnd'
!      PRINT *,eta_zp_t%ubnd
!      PRINT '(a)','t e'
!      PRINT *,size(eta_zp_t%e)
!      PRINT *,eta_zp_t%e
      


    
    ! Loop over sidebands ----------------------------------------------------
    call Trace_Begin ( me_SidebandLoop, 'ForwardModel.SidebandLoop', &
      & cond=toggle(emit) .and. levels(emit) > 0 )
    do thisSideband = fwdModelConf%sidebandStart, fwdModelConf%sidebandStop, 2

      call Trace_Begin ( me_Sideband, 'ForwardModel.Sideband ', &
        & index=thisSideband, cond=toggle(emit) .and. levels(emit) > 1 )

      sx = ( thisSideband + 3 ) / 2 ! [-1,+1] => [1,2]

      ! Work out Frequency averaging / LBL / PFA / Derivatives steering.
      ! See the decision table in Frequency_Average for the meaning of
      ! values of Frq_Avg_Sel.
      frq_avg_sel = 0
      if ( fwdModelConf%do_freq_avg .or. &
        &  fwdModelConf%anyPFA(sx) .and. .not. &
        &  fwdModelConf%anyLBL(sx) )  frq_avg_sel = ior(frq_avg_sel, 8)
      if ( fwdModelConf%anyPFA(sx) )  frq_avg_sel = ior(frq_avg_sel, 4)
      if ( fwdModelConf%anyLBL(sx) )  frq_avg_sel = ior(frq_avg_sel, 2)
      if ( any_der )                  frq_avg_sel = ior(frq_avg_sel, 1)

      ! Now, allocate gl_slabs arrays
      call allocateSlabs ( gl_slabs, max_f, &
        & fwdModelConf%catalog(thisSideband,:fwdModelConf%cat_size(sx)), &
        & moduleName, temp_der )

      call frequency_setup_1 ( tan_press, grids )
      
      ! Loop over pointings --------------------------------------------------
      call Trace_Begin ( me_PointingLoop, 'ForwardModel.PointingLoop', &
        & cond=toggle(emit) .and. levels(emit) > 2 )

      if ( switchDetail( switches, 'nfmpt', options='-fc' ) > -1 ) then
       ! no op; maybe print
        if ( toggle(emit) .and. levels(emit) > 2 ) &
          & call output( '(Skipping loop of pointings)', advance='yes' )
      else if ( .not. FwdModelConf%generateTScat ) then

        do ptg_i = 1, no_tan_hts

          vel_rel = est_los_vel(ptg_i) / speedOfLight

          ! If we're doing frequency averaging, get the frequencies we need for
          ! this pointing.

          if ( associated(whichPointingGrid) ) &
            & call frequency_setup_2 ( (1.0_rp - vel_rel) * &
            & whichPointingGrid%oneGrid(grids(ptg_i))%frequencies )

          ! Do the ray tracing
          if ( .not. usingQTM ) then
            ! Calculate R_Eq in the orbit plane projected ellipse, for which
            ! Earthradc_sq has the same value for all tangent points in the
            ! non-QTM case because they're all in the same plane.
            r_eq = get_r_eq ( tan_phi(ptg_i), earthradc_sq(ptg_i) )
            call one_pointing ( ptg_i, vel_rel, r_eq, tan_phi=tan_phi(ptg_i), &
              &                 tan_press=tan_press(ptg_i), &
              &                 est_scGeocAlt=est_scGeocAlt(ptg_i) )
          else
            ! For the QTM case, R_Eq might have different values for each
            ! tangent point in the QTM.
            r_eq = radius_of_curvature_normal ( QTM_paths(ptg_i)%lines(2,1), &
                                              & tan_pt_geod(ptg_i)%surf_ECR(.false.) )
            call one_pointing ( ptg_i, vel_rel, r_eq,               &
              &                 tan_loc=tan_pt_geod(ptg_i),         &
              &                 qtm=QTM,                            &
              &                 tan_press=tan_press(ptg_i),         &
              &                 est_scGeocAlt=est_scGeocAlt(ptg_i), &
              &                 path=QTM_paths(ptg_i), s=s,         &
              &                 f_and_v=F_and_V(min(ptg_i,size(F_And_V))) )
          end if
        end do ! ptg_i

      else
        ! Do ray tracing at specified scattering angles
        call generate_TScat ( fwdModelConf )
      end if
      
      call Trace_End ( 'ForwardModel.PointingLoop', &
        & cond=toggle(emit) .and. levels(emit) > 2 )

      if ( .not. fwdModelConf%generateTScat ) then
        IF ( thisSideband == fwdModelConf%sidebandStart ) THEN
           ! Same for both sidebands, so only need to do it once
          if ( .not. usingQTM ) then
            call convolution_setup ( dh_dz_out, dx_dh_out, dxdT_surface, &
              & dxdT_tan, Q_EarthRadC_sq, Est_ScGeocAlt, &
              & FwdModelConf, &
              & FwdModelExtra, FwdModelIn, Grids_f, Grids_tmp, &
              & L1BMIF_TAI, MAF, MIFDeadTime, &
              & PhiTan, PTan, RefGPH, SCGeocAlt, &
              & Surf_Angle, Chi_MIF_Tan, Tan_Phi, Tan_Press, &
              & WindowFinish, WindowStart )
          else ! QTM
            call convolution_setup ( dh_dz_out, dx_dh_out, dxdT_surface, &
              & dxdT_tan, Q_EarthRadC_sq, Est_ScGeocAlt, &
              & FwdModelConf, &
              & FwdModelExtra, FwdModelIn, Grids_f, Grids_tmp, &
              & L1BMIF_TAI, MAF, MIFDeadTime, &
              & PhiTan, PTan, RefGPH, SCGeocAlt, &
              & Surf_Angle, Chi_MIF_Tan, Tan_Phi, Tan_Press, &
              & WindowFinish, WindowStart, &
              & n_path_c(1), scECR_MIF, ECRtoFOV )
          end if
       END IF
!       PRINT '(a)','About to convolve this'
!       PRINT *,SIZE(k_magfield,1),SIZE(k_magfield,2),SIZE(k_magfield,3)
!       DO i = 1, 70
!         PRINT *,k_magfield(1,i,1),k_magfield(1,i,2),k_magfield(1,i,3)
!       end do
        call convolution & ! or interpolate to ptan
          ( dh_dz_out, dx_dh_out, dx_dT, dxdT_surface, &
          & dxdT_tan, d2x_dxdT, F_and_V, &
          & Q_EarthRadC_sq, Est_ScGeocAlt, FirstSignal, FmStat, &
          & FwdModelConf, FwdModelExtra, FwdModelIn, FwdModelOut, &
          & Grids_f, Grids_n, Grids_tmp, Grids_v, Grids_w, Grids_mag, &
          & H_Atmos, K_Atmos, K_Spect_DN, K_Spect_DV, K_Spect_DW, &
          & K_Temp, k_magfield, L1BMIF_TAI, MIFDeadTime, PTan, PTan_Der, &
          & Ptg_Angles, Radiances, Sideband, S_T, Surf_Angle, &
          & Chi_MIF_Tan, Tan_Phi, Temp, ThisSideband, &
          & Jacobian, ExtraJacobian, Hessian )
!        PRINT '(a)','after convoltion'
!        PRINT *,k_magfield_frq
!        PRINT '(a)','after convoltion T'
!        PRINT *,k_temp
      end if

      ! Frequency averaging and any LBL:
      if ( .not. associated(frequencies,channelCenters) ) &
        & call deallocate_test ( frequencies, 'frequencies',       moduleName )

      ! Deallocate maxNoPtgFreqs-sized stuff
      call deallocate_test ( inc_rad_path,     'Inc_Rad_Path',     moduleName )
      call deallocate_test ( radv,             'RadV',             moduleName )
      if ( FwdModelConf%anyLBL(sx) ) then
        call deallocate_test ( t_script_LBL,   'T_Script_LBL',     moduleName )
        call destroy_tau ( tau_LBL,            'Tau_LBL',          moduleName )
      end if

      call DestroyCompleteSlabs ( gl_slabs )
      if ( temp_der ) &
        & call deallocate_test ( k_temp_frq,   'k_temp_frq',       moduleName )
      if ( hmag_der .or. htheta_der .or. hphi_der) &
        & call deallocate_test ( k_magfield_frq, 'k_magfield_frq', moduleName )

      call deallocate_test ( k_atmos_frq,  'k_atmos_frq',          moduleName )
      call deallocate_test ( h_atmos_frq,  'h_atmos_frq',          moduleName )

      call deallocate_test ( k_spect_dw_frq, 'k_spect_dw_frq',     moduleName )
      call deallocate_test ( k_spect_dn_frq, 'k_spect_dn_frq',     moduleName )
      call deallocate_test ( k_spect_dv_frq, 'k_spect_dv_frq',     moduleName )

      call trace_end ( 'ForwardModel.Sideband ', index=thisSideband, &
        & cond=toggle(emit) .and. levels(emit) > 1 )

    end do            ! End of loop over sidebands -------------------------

    call Trace_End ( 'ForwardModel.SidebandLoop', &
      & cond=toggle(emit) .and. levels(emit) > 0 )

    if ( FwdModelConf%anyPFA(1) .or. FwdModelConf%anyPFA(2) ) &
      & call destroy_tau ( tau_PFA, "Tau_PFA", moduleName )

    !  **** DEBUG Printing cycle ...

! *** Create *seez* file for "nasty" purposes:

    if ( print_Seez ) call dump_print_code

! *** End of include

    call GetNameOfSignal ( firstSignal, sigName )
    i = index(sigName, '.B')
    j = index(sigName(i+2:), '.' )
    if ( j /= 0 ) sigName(i+j+1:) = ''

    if ( print_Rad > -1 ) then
      print '(6a)', 'Signal: ', trim(sigName), &
        & ' Convolution: ', trim(merge('ON ','OFF',FwdModelConf%do_conv)), &
        & ' Frequency Averaging: ', trim(merge('ON ','OFF', &
          & FwdModelConf%do_freq_avg .and. any(fwdModelConf%anyLBL)))

      k = ptan%template%noSurfs
      print "( /'ptan\ ',i3.3)", k
      print "( 4(3x, f11.7) )", Ptan%values(1:k,maf)

      do i = 1, noUsedChannels
        print "(/, 'ch', i3.3, '_pfa_rad\ ', i3.3 )", channels(i)%used, k
        thisRadiance => GetQuantityForForwardModel (fwdModelOut,  &
          & quantityType=merge(l_tscat,l_radiance,fwdModelConf%GenerateTScat), &
          & signal=fwdModelConf%signals(channels(i)%signal)%index, &
          & sideband=sideband, config=fwdModelConf )
        channel = channels(i)%used - channels(i)%origin + 1
!         j = thisRadiance%template%noChans
!         This caused a bogus bounds-check violation with -check bounds, or
!         undefined variable violation with -check undef, from ifort 17.0.0.098
!         print "( 4(2x, 1pg15.8) )", &
!           & thisRadiance%values(channel:channel+j*(k-1):j, maf)
        print "( 1p, 4g17.8 )", thisRadiance%value3(channel,:,maf)
      end do

    end if

    ! **** End of Printing cycle ...

    call destroygrids_t ( grids_tscat )
    call destroygrids_t ( grids_salb )
    call destroygrids_t ( grids_cext )

    if ( FwdModelConf%incl_cld ) then
      call DestroyVectorQuantityValue ( scat_src, forWhom='scat_src' )
      call DestroyVectorQuantityValue ( scat_alb, forWhom='scat_alb' )
      call DestroyVectorQuantityValue ( cld_ext, forWhom='cld_ext' )
    end if

    ! If hGrids were QTM, destroy intermediate minor-frame quantities
    if ( usingQTM ) then
      call DestroyVectorQuantityValue ( Q_EarthRadC_sq, forWhom='Q_EarthRadC_sq' )
      call DestroyVectorQuantityValue ( Q_PhiTan, forWhom='Q_PhiTan' )
      call DestroyVectorQuantityValue ( Q_TanHt, forWhom='Q_TanHt' )
    end if

    call trace_end ( 'FullForwardModelAuto, MAF=', fmStat%maf, &
      & cond=toggle(emit) )

  contains

  ! .............................................  Announce_Error  .....
    subroutine Announce_Error ( Message )
    ! Announce Message using MLSMessage.  Include the configuration name
    ! in the message.
      use MLSMessageModule, only: MLSMSG_Error
      use MoreMessage, only: MLSMessage
      character(len=*), intent(in) :: Message
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "With config(%S): " // message, datum=fwdModelConf%name )
    end subroutine Announce_Error

  ! .......................................  Both_Sidebands_Setup  .....
    subroutine Both_Sidebands_Setup
    ! Setup stuff done for both sidebands, other than output angles
    ! for convolution and stuff for clouds.

      use Geometry, only: Orbit_Plane_Minor_Axis_sq
      use Interpolate_MIF_to_Tan_Press_m, only: Interpolate_MIF_to_Tan_Press
      use Intrinsic, only: L_EarthRefl, L_ECRtoFOV, L_GPH,  &
        & L_LOSVel, L_SurfaceHeight, L_OrbitInclination, L_RefGPH, &
        & L_ScGeocAlt, L_SpaceRadiance
      use ManipulateVectorQuantities, only: DoHGridsMatch
      use Tangent_Quantities_m, only: Tangent_Quantities

      integer :: Me = -1  ! String index for trace
      integer :: SigInd   ! Signal index, loop counter

      call trace_begin ( me, 'ForwardModel.Both_Sidebands_Setup', &
        & cond=toggle(emit)  .and. levels(emit) > 0  )

      fmStat%flags = 0   ! Assume no errors

      MAF = fmStat%maf

      if ( usingQTM ) then
        ! Compute minor-frame tangent quantities from the minor-frame
        ! lines-of-sight (Q_LOS) provided by level 1.
        call tangent_quantities ( MAF, fwdModelConf, fwdModelIn, fwdModelExtra, &
                                & Q_LOS, Q_EarthRadC_sq, Q_PhiTan, Q_TanHt )
        phiTan => Q_PhiTan ! Use the calculated one, not the one in FwdModelIn
      end if

      temp_der = fwdModelConf%temp%derivOK
      hmag_der = fwdmodelconf%hmag_der
      htheta_der = fwdmodelconf%htheta_der
      hphi_der = fwdmodelconf%hphi_der
      atmos_der = present ( jacobian ) .and. FwdModelConf%atmos_der
      atmos_second_der = present ( hessian ) .and. FwdModelConf%atmos_second_der

      spect_der = present ( jacobian ) .and. FwdModelConf%spect_der
      spect_der_center = spect_der .and. size(fwdModelConf%lineCenter) > 0
      spect_der_width = spect_der .and. size(fwdModelConf%lineWidth) > 0
      spect_der_width_TDep = spect_der .and. size(fwdModelConf%lineWidth_TDep) > 0
      spect_der = spect_der .and. &
        & ( spect_der_center .or. spect_der_width .or. spect_der_width_TDep )

      any_der = temp_der .OR. atmos_der .OR. spect_der .OR. hmag_der .OR. &
           htheta_der .or. hphi_der

      ! Work out what we've been asked to do -----------------------------------

      beta_group => fwdModelConf%beta_group
      channels => fwdModelConf%channels
      DACsStaging => fwdModelConf%DACsStaging
      usedDACSSignals => fwdModelConf%usedDACSSignals
      noUsedDacs = size(usedDACSSignals)

      ! Identify the vector quantities we're going to need.
      ! The key is to identify the signal we'll be working with first
      firstSignal => fwdModelConf%signals(1) ! Config has verified that signals
        ! are all for same radiometer (actually LO), module and sideband
      sideband = merge ( 0, firstSignal%sideband, fwdModelConf%forceFoldedOutput )

      ! Start sorting out stuff from state vector ------------------------------

      ! Identify the appropriate state vector components
      ! VMRS are in Beta_group%qty, gotten by get_species_data.
      ! First the minor frame quantities:
      losVel => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_losVel, instrumentModule=firstSignal%instrumentModule, &
        & config=fwdModelConf )
      scGeocAlt => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_scGeocAlt, config=fwdModelConf )

      ! Now the others:
      earthRefl => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_earthRefl, config=fwdModelConf )
      refGPH => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_refGPH, config=fwdModelConf )
      spaceRadiance => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_spaceRadiance, config=fwdModelConf )
      surfaceHeight => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_surfaceHeight, config=fwdModelConf, noError=.true. )
      if ( FwdModelConf%polarized .or. usingQTM ) then
        ECRtoFOV => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_ECRtoFOV, config=fwdModelConf )
      end if
      if ( FwdModelConf%incl_cld ) then
        gph => GetQuantityForForwardModel (fwdModelIn, fwdModelExtra, &
          & quantityType=l_gph, config=fwdModelConf )
      end if
      if ( FwdModelConf%generateTScat .or. .not. usingQTM ) then
        orbIncline => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_orbitInclination, config=fwdModelConf )
      end if

      ! Check that RefGPH and Temp have the same hGrid.  This is not checked in
      ! Construct or when the config is created.
      if ( .not. doHGridsMatch ( refGPH, temp ) ) call Announce_Error ( &
        & 'Different horizontal grids for refGPH and temperature' )

      if ( .not. fwdModelConf%generateTScat ) then
        ! Check that we have radiances for the channels that are used
        do sigInd = 1, size(fwdModelConf%signals)
          ! This just emits an error message and stops if we don't have a radiance.
          ! We don't use the vector quantity -- at least not right away.  We get
          ! it again later.
          thisRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
            & signal=fwdModelConf%signals(sigInd)%index, sideband=sideband )
        end do
      end if

      windowStart = grids_tmp%windowStart(1)
      windowFinish = grids_tmp%windowFinish(1)

  ! Compute reference Gauss Legendre (GL) grid ------------------------------

      call compute_GL_grid ( z_psig, z_glgrid )

      ! interpolate Tan_Phi, ScGeocAlt and LOSVel from MIFs to pointings
      if ( .not. FwdModelConf%generateTScat ) then
        if ( .not. usingQTM ) then
          ! At pressure levels PTan, interpolate
          ! phitan,  scgeocalt,     losvel, to
          ! tan_phi, est_scgeocalt, est_los_vel
          ! at pressure levels given by Tan_Press
          call interpolate_MIF_to_tan_press ( nlvl, maf, ptan, phitan, &
            & scgeocalt, losvel, tan_press, &
            & tan_phi, est_scgeocalt, est_los_vel )

          ! Create minor-frame values, all the same:
          ! Q_earthRadC_sq is an array of orbit-plane-projected Earth radii,
          ! all the same, used for convolution.
          call cloneVectorQuantity ( Q_earthRadC_sq, phitan )
          Q_EarthRadC_sq%values = &
            & orbit_plane_minor_axis_sq ( orbIncline%values(1,MAF) * deg2rad )
          ! Create coarse-grid tangent-point values, all the same
          earthradc_sq = &
            & orbit_plane_minor_axis_sq ( orbIncline%values(1,MAF) * deg2rad )
        else
          ! At pressure levels PTan, interpolate
          ! phitan,  scgeocalt,     losvel,      Q_EarthRadC_sq, Q_TanHt to
          ! tan_phi, est_scgeocalt, est_los_vel, earthRadC_sq,   tanHt
          ! at pressure levels given by Tan_Press
          call interpolate_MIF_to_tan_press ( nlvl, maf, ptan, phitan, &
            & scgeocalt, losvel, tan_press, &
            & tan_phi, est_scgeocalt, est_los_vel, &
            & Q_EarthRadC_sq, Q_TanHt, &
            & earthRadC_sq,   tanHt, tan_pt_Geod  )
        end if
      else ! Generating TScat
        ! Compute the square of the minor axis of the orbit plane projected
        ! Earth ellipse.
        earthradc_sq = orbit_plane_minor_axis_sq ( orbIncline%values(1,maf) * deg2rad )
      end if

      ! Now, allocate other variables we're going to need later --------

      if ( FwdModelConf%anyPFA(1) .or. FwdModelConf%anyPFA(2) ) then
        call allocate_test ( tau_PFA%tau, max_c, noUsedChannels, 'Tau_PFA%Tau', &
          & moduleName )
        call allocate_test ( tau_PFA%i_stop, noUsedChannels, 'Tau_PFA%I_Stop', &
          & moduleName )
      end if

      call trace_end ( 'ForwardModel.Both_Sidebands_Setup', &
        & cond=toggle(emit) .and. levels(emit) > 0  )

    end subroutine Both_Sidebands_Setup

  ! ................................................  Cloud_Setup  .....
    subroutine Cloud_Setup

      use Intrinsic, only: L_BOUNDARYPRESSURE, L_CLEAR
      use Load_Sps_Data_m, only: Modify_values_for_supersat
      use ManipulateVectorQuantities, only: FindOneClosestInstance
      use VectorsModule, only: CreateVectorValue

      integer :: Ispec                        ! Species index in cloud model
      integer :: J                            ! Loop inductor, subscript
      type (VectorValue_T), pointer :: VMR    ! Quantity

      ! Find the instance in temp that is closest to the current MAF
      inst = FindOneClosestInstance ( temp, ptan, fmStat%maf )

      ! checking done in ForwardModelSupport%ConstructForwardModelConfig
      nspec = no_mol ! Will be at least 3 if l_n2o is included, because
                     ! l_h2o and l_o3 are required
      vmrarray = 0.0_r8

      do j = 1, nspec      ! Loop over species

        if ( fwdModelConf%molecules(j) == l_h2o ) then
          ispec = 1
        else if ( fwdModelConf%molecules(j) == l_o3 ) then
          ispec = 2
        else if ( fwdModelConf%molecules(j) == l_n2o ) then
          ispec = 3
        else
          cycle
        end if

        vmr => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_vmr, molecule=fwdModelConf%molecules(j) )

        call InterpolateValues ( vmr%template%surfs(:,1),    &    ! Old X
                               & vmr%values(:,inst),         &    ! Old Y
                               & temp%template%surfs(:,1),   &    ! New X
                               & vmrArray(:,ispec),          &    ! New Y
                               & 'Linear', extrapolate='Clamp' )  ! Options

      end do ! End of Loop over species

      scat_src%template = temp%template
      scat_src%template%noInstances = fwdModelConf%num_scattering_angles
      scat_alb%template = temp%template
      scat_alb%template%noInstances = 2
      cld_ext%template = temp%template
      cld_ext%template%noInstances = 2

      call createVectorValue ( scat_src, 'scat_src' )
      call createVectorValue ( scat_alb, 'scat_alb' )
      call createVectorValue (  cld_ext, 'cld_ext' )

      if ( FwdModelConf%i_saturation /= l_clear ) then
        boundaryPressure => GetQuantityForForwardModel ( fwdModelIn, &
          & fwdModelExtra, quantityType=l_boundaryPressure, config=fwdModelConf )
        ! modify h2o mixing ratio if a special supersaturation is requested
        call modify_values_for_supersat ( fwdModelConf, grids_f, h2o_ind, &
          & grids_tmp, boundaryPressure )
      end if

    end subroutine Cloud_Setup

  ! ............................................  Dump_Print_Code  .....
    subroutine Dump_Print_Code
      include "dump_print_code.f9h"
    end subroutine Dump_Print_Code

  ! ..........................................  Frequency_Average  .....
    subroutine Frequency_Average ( Ptg_i )

      ! Here we either frequency average to get the unconvolved radiances, or
      ! we just store what we have as we're using monochromatic channels

      use Freq_Avg_m, only: Freq_Avg, Freq_Avg_DACS
      use SCRT_dn_m, only: SCRT_PFA

      integer, intent(in) :: Ptg_i ! Pointing index

      integer :: C, ShapeInd
      integer :: Me = -1           ! String index for trace

      call trace_begin ( me, 'ForwardModel.Frequency_Average', &
        & cond=toggle(emit) .and. levels(emit) > 4 )

! Conditions:  8 Frequency averaging?   -  N N - - N N Y Y Y Y
!  - means     4 PFA?                   N  N N Y Y Y Y N N Y Y
! "macht nix"  2 LBL?                   N  Y Y N N Y Y Y Y Y Y
!              1 Derivatives?           -  N Y N Y N Y N Y N Y
!                ---------------------------------------------
!                Frq_Avg_Sel value     0,1 2 3 4 5 6 7
!                                      8,9     C D     A B E F
! ============================================================
! Actions:
! Impossible                            x
! Radiances = RadV                         x x x x
! K = K_frq                                  x   x
! Frq Avg path integrated rad                          x x x
! Combine total path radiances                     x x
! Frq Avg path integrated LBL deriv                      x
! Frq Avg rad Along Path                                     x
! Combine radiances along path                             x x
! Combine LBL and PFA derivs                         x       x

      ! Do DACs stuff for all DACs channels first
      select case ( frq_avg_sel )
      case ( 10 : 15 )
        do c = 1, noUsedDACS
          shapeInd = MatchSignal ( dacsFilterShapes%filter%signal, &
            & fwdModelConf%signals(usedDacsSignals(c)), sideband = thisSideband )
          call Freq_Avg_DACS ( frequencies, DACSFilterShapes(shapeInd), &
            & RadV(:noFreqs), DACsStaging(:,c) )
        end do
      end select

      select case ( frq_avg_sel )
      case ( 0, 1, 8, 9 )  ! Impossible: No LBL or PFA
      case ( 2, 3, 4, 5, 12, 13 ) ! Just copy radiance. PFA + DACS - LBL impossible
        radiances(:,ptg_i) = radV(:noFreqs)
      case ( 6, 7 ) ! PFA + monochromatic LBL
        do c = 1, noUsedChannels
          if ( channels(c)%dacs == 0 ) then
            ! Combine LBL and PFA Tau's to get radiances.
            call SCRT_PFA (c, tau_LBL, tau_PFA, t_script_pfa, radiances(c:c,ptg_i) )
          else
            radiances(c,ptg_i) = DACsStaging(channels(c)%used,channels(c)%dacs)
          end if
        end do
      case ( 10, 11, 14 )   ! Frq Avg path integrated radiance
        ! Now go through channel by channel
        do c = 1, noUsedChannels
          if ( channels(c)%dacs == 0 ) then
            if ( fwdModelConf%anyPFA(sx) ) then ! frq_avg_sel = 14
            ! Combine LBL and PFA Tau's to get radiances.
            ! It's OK to combine the Tau's before doing the frequency
            ! averaging because the filter function is normalized.
              call SCRT_PFA (c, tau_LBL, tau_PFA, t_script_pfa, radV(:noFreqs) )
            end if
            shapeInd = channels(c)%shapeInds(sx)
            call freq_Avg ( frequencies, &
              &   filterShapes(shapeInd)%filterGrid,  &
              &   filterShapes(shapeInd)%filterShape, &
              &   radV(:noFreqs), radiances(c,ptg_i) )
          else
            radiances(c,ptg_i) = DACsStaging(channels(c)%used,channels(c)%dacs)
          end if
        end do
      case ( 15 )   ! Frq Avg rad along path
        ! For every channel, we've frequency averaged the incremental radiance
        ! at every point along the path, giving Rad_Avg_Path for every channel
        ! and every point along the path.  Rad_Avg_Path was multiplied by
        ! Tau_PFA in One_Frequency to combine LBL and PFA contributions.
        do c = 1, noUsedChannels
          if ( channels(c)%dacs == 0 ) then
            radiances(c,ptg_i) = radV(c) ! Computed in One_Frequency
          else
            radiances(c,ptg_i) = DACsStaging(channels(c)%used,channels(c)%dacs)
          end if
        end do
      end select

      if ( any_der ) then
        call frequency_average_derivatives ( ptg_i, combine=iand(frq_avg_sel,7) == 7 )
        if ( atmos_second_der ) then
          call frequency_avg_second_derivs ( ptg_i, frq_avg_sel == 15 )
        end if

      end if


      call trace_end ( 'ForwardModel.Frequency_Average', &
        & cond=toggle(emit) .and. levels(emit) > 4 )
    end subroutine Frequency_Average

  ! ...............................  Frequency_Average_Derivative  .....
    subroutine Frequency_Average_Derivative ( Grids, K_Frq, K, Mol, Combine )
      ! Frequency average or simply copy K_Frq to give K, the final
      ! Jacobian.

      use Freq_Avg_m, only: Freq_Avg, Freq_Avg_DACS

      type(grids_T), intent(in) :: Grids
      real(rp), intent(in) :: K_FRQ(:,:)    ! To be averaged  Frq X Grid
      real(r4), intent(inout) :: K(:,:)     ! Averaged        Chan X Grid
      integer, intent(in) :: Mol            ! Which molecule
      logical, intent(in) :: Combine        ! "Combine LBL and PFA"

      integer :: C, ShapeInd
      integer :: Me = -1 ! String index for trace
      real(rp) :: R      ! Frequency-averaged value
      integer :: SV_I    ! State-vector index

      call trace_begin ( me, 'ForwardModel.Frequency_Average_Derivative=', &
           & index=mol, cond=toggle(emit) .AND. levels(emit) > 6 )

!      PRINT '(a)','inside the frequency derivative'
!      PRINT *,SIZE(k_frq,1),SIZE(k_frq,2)
!      PRINT *,SIZE(k,1),SIZE(k,2)

      if ( combine ) then
        ! Simply add newly-computed PFA derivatives in K_frq to
        ! previously-averaged or monochromatic LBL derivatives in K.
        ! Remember that for PFA, the frequency dimension has extent
        ! noUsedChannels, not maxNoPtgFreqs.  The first dimension of K_frq is
        ! at least noUsedChannels, so we're guaranteed this will fit.
        do c = 1, noUsedChannels
          do sv_i = grids%l_v(mol-1)+1, grids%l_v(mol)
            k(c,sv_i) = k(c,sv_i) + k_frq(c,sv_i)
          end do
        end do
        call trace_end ( 'ForwardModel.Frequency_Average_Derivative=', &
          & index=mol, cond=toggle(emit) .and. levels(emit) > 6 )
        return
      end if

      ! Only possible values for Frq_avg_sel here are 3, 5, 11, 13, 15.
      ! See Frequency_Average for definition of Frq_avg_sel.

      select case ( frq_avg_sel )
      case ( 3, 5, 7, 13 ) ! Not frequency averaging, or PFA alone, or
                      ! Monochromatic LBL + PFA + Derivatives after
                      ! first call to One_Frequency; copy.
                      ! Shouldn't get here after both LBL and PFA passes
                      ! through One_Frequency -- should be handled by
                      ! combine == .true.
        do sv_i = grids%l_v(mol-1)+1, grids%l_v(mol)
          if ( grids%deriv_flags(sv_i) ) then
            k(:,sv_i) = min ( max ( k_frq(:,sv_i), &
                          &   real(-huge(0.0_r4), rp ) ), &
                          &   real( huge(0.0_r4), rp ) )
          else
            k(:,sv_i) = 0.0
          end if
        end do

      case ( 11, 15 ) ! See Frequency_Average.  LBL and PFA have already
                      ! been combined.  It's OK to frequency average the
                      ! PFA contribution again because the filter function
                      ! is normalized.
        ! Do DACs stuff for all DACs channels first
        do c = 1, noUsedDACS
          shapeInd = MatchSignal ( dacsFilterShapes%filter%signal, &
            & fwdModelConf%signals(usedDacsSignals(c)), sideband = thisSideband )
          do sv_i = grids%l_v(mol-1)+1, grids%l_v(mol)
            call Freq_Avg_DACS ( frequencies, DACSFilterShapes(shapeInd), &
              & k_frq(:,sv_i), DACsStaging2(:,sv_i,c) )
          end do                  ! Surface loop X Instance loop
        end do
        ! Now go through channel by channel
        do c = 1, noUsedChannels
          shapeInd = channels(c)%shapeInds(sx)
          do sv_i = grids%l_v(mol-1)+1, grids%l_v(mol)
            if ( grids%deriv_flags(sv_i) ) then
              if ( channels(c)%dacs == 0 ) then
                call Freq_Avg ( frequencies,            &
                  & FilterShapes(shapeInd)%FilterGrid,  &
                  & FilterShapes(shapeInd)%FilterShape, &
                  & k_frq(:,sv_i), r )
              else
                r = DACsStaging2 ( channels(c)%used, sv_i, channels(c)%dacs )
              end if
            else
              r = 0.0
            end if
            k(c,sv_i) = r
          end do                ! Grid loop
        end do                  ! Channel loop
      case default              ! Impossible
      end select                ! Frequency averaging or not

      call trace_end ( 'ForwardModel.Frequency_Average_Derivative=', &
        & index=mol, cond=toggle(emit) .and. levels(emit) > 6 )

    end subroutine Frequency_Average_Derivative

  ! ..............................  Frequency_Average_Derivatives  .....
    subroutine Frequency_Average_Derivatives ( Ptg_i, Combine )
      integer, intent(in) :: Ptg_i          ! Pointing index
      logical, intent(in) :: Combine        ! "Combine LBL and PFA"

      ! Frequency Average the temperature derivatives with the appropriate
      ! filter shapes

      integer :: K       ! Loop inductor
      integer :: Me = -1 ! String index for trace
      integer :: UB      ! Upper bound for first dimension of k_..._frq

      call trace_begin ( me, 'ForwardModel.Frequency_Average_Derivatives', &
        & cond=toggle(emit) .and. levels(emit) > 5 )

      ub = noFreqs
      if ( combine ) ub = max(ub,noUsedChannels)

      if ( temp_der ) call frequency_average_derivative ( grids_tmp, &
        &               k_temp_frq(:ub,:), k_temp(:,ptg_i,:), 1, combine )

      IF ( hmag_der .OR. htheta_der .OR. hphi_der) &
           CALL frequency_average_derivative ( grids_mag, &
           &               k_magfield_frq(:ub,:), k_magfield(:,ptg_i,:), &
           & 1, combine )
      
      ! Frequency Average the atmospheric derivatives with the appropriate
      ! filter shapes
      if ( atmos_der ) then
        do k = 1, no_mol
          if ( fwdModelConf%moleculeDerivatives(k) ) then
            call frequency_average_derivative ( grids_f, &
              & k_atmos_frq(:ub,:), k_atmos(:,ptg_i,:), k, combine )
          else
            k_atmos(:,ptg_i,grids_f%l_v(k-1)+1:grids_f%l_v(k)) = 0.0
          end if
        end do                        ! Loop over major molecules
      end if                          ! Want derivatives for atmos

      ! Frequency Average the spectroscopic derivatives with the appropriate
      ! filter shapes

      do k = 1, size(fwdModelConf%lineCenter)
        call frequency_average_derivative &
          & ( grids_v, k_spect_dv_frq(:ub,:), k_spect_dv(:,ptg_i,:), k, &
          & combine )
      end do
      do k = 1, size(fwdModelConf%lineWidth)
        call frequency_average_derivative &
          & ( grids_w, k_spect_dw_frq(:ub,:), k_spect_dw(:,ptg_i,:), k, &
          & combine )
      end do
      DO k = 1, SIZE(fwdModelConf%lineWidth_TDep)
        call frequency_average_derivative &
          & ( grids_n, k_spect_dn_frq(:ub,:), k_spect_dn(:,ptg_i,:), k, &
          & combine )
      END DO

!      IF (hmag_der) THEN
!         k_magfield(:,ptg_i,:) = 5.0
!      endif

      call trace_end ( 'ForwardModel.Frequency_Average_Derivatives', &
        & cond=toggle(emit) .and. levels(emit) > 5 )

    end subroutine Frequency_Average_Derivatives

  ! .....................................  Frequency_Avg_Rad_Path  .....
    subroutine Frequency_Avg_Path ( Frequencies, Path_Freq, Path_Chan )
      ! PFA + LBL + Derivatives, and maybe or maybe not frequency averaging.
      ! Frq_Avg_Sel = 7 (no frequency averaging) or 15 (frequency averaging).
      ! For every channel, frequency average or copy Path_Freq at every
      ! point along the path, giving Path_Chan for every channel and
      ! every point along the path.

      use Freq_Avg_m, only: Freq_Avg

      real(r8), intent(in) :: Frequencies(:)
      real(rp), intent(in) :: Path_Freq(:,:)  ! Path X Frequencies
      real(rp), intent(out) :: Path_Chan(:,:) ! Path X Channels

      integer :: C, P, ShapeInd

      if ( fwdModelConf%do_freq_avg ) then ! Frequency average
        do c = 1, noUsedChannels
          shapeInd = channels(c)%shapeInds(sx)
          if ( channels(c)%dacs == 0 ) then
            do p = 1, size(path_freq,1)
              call Freq_Avg ( frequencies,            &
                      & FilterShapes(shapeInd)%FilterGrid,  &
                      & FilterShapes(shapeInd)%FilterShape, &
                      & path_freq(p,:), path_chan(p,c) )
            end do
          end if
        end do
      else ! Copy
        path_chan = path_freq
      end if
    end subroutine Frequency_Avg_Path

  ! ...............................  Frequency_Avg_Second_Derivative  .....
    subroutine Frequency_Avg_Second_Derivative ( Grids, H_Frq, H, Mol1, Mol2, Combine )

      ! Frequency average or simply copy H_Frq to give H, the final HESSIAN

      use Freq_Avg_m, only: Freq_Avg, Freq_Avg_DACS

      type(grids_T), intent(in) :: Grids
      real(rp), intent(in) :: H_FRQ(:,:,:)    ! To be averaged  Frq X Grid X Grid
      real(r4), intent(inout) :: H(:,:,:)     ! Averaged        Chan X Grid X Grid
      integer, intent(in) :: Mol1, Mol2       ! Which molecules
      logical, intent(in) :: Combine          ! "Combine LBL and PFA"

      integer :: C, ShapeInd
      real(rp) :: H_AVG      ! Frequency-averaged value
      integer :: SV_Q, SV_R  ! State-vector indexes

      if ( combine ) then
        ! Simply add newly-computed PFA derivatives in K_frq to
        ! previously-averaged LBL derivatives in K.  Remember that
        ! for PFA, the frequency dimension has extent noUsedChannels,
        ! not maxNoPtgFreqs.  The first dimension of K_frq is at least
        ! noUsedChannels, so we're guaranteed this will fit.
        do c = 1, noUsedChannels
          do sv_q = grids%l_v(mol1-1)+1, grids%l_v(mol1)
            do sv_r = grids%l_v(mol2-1)+1, grids%l_v(mol2)
              h(c,sv_q,sv_r) = h(c,sv_q,sv_r) + h_frq(c,sv_q,sv_r)
            end do
          end do
        end do
        return
      end if

      ! Only possible values for Frq_avg_sel here are 3, 7, 11, 13, 15.
      ! See Frequency_Average for definition of Frq_avg_sel.

      select case ( frq_avg_sel )
      case ( 11, 15 ) ! See Frequency_Average.
        ! Do DACs stuff for all DACs channels first
        do c = 1, noUsedDACS
          shapeInd = MatchSignal ( dacsFilterShapes%filter%signal, &
            & fwdModelConf%signals(usedDacsSignals(c)), sideband = thisSideband )

          do sv_q = grids%l_v(mol1-1)+1, grids%l_v(mol1)
            do sv_r = grids%l_v(mol2-1)+1, grids%l_v(mol2)
              call Freq_Avg_DACS ( frequencies, DACSFilterShapes(shapeInd), h_frq(:,sv_q,sv_r), DACsStaging3(:,sv_q,sv_r,c) )
            end do
          end do                  ! Surface loop X Instance loop

        end do

       ! Now go through channel by channel
        do c = 1, noUsedChannels
          shapeInd = channels(c)%shapeInds(sx)
          do sv_q = grids%l_v(mol1-1)+1, grids%l_v(mol1)
            do sv_r = grids%l_v(mol2-1)+1, grids%l_v(mol2)
              if ( grids%deriv_flags(sv_q) .and. grids%deriv_flags(sv_r) ) then
                if ( channels(c)%dacs == 0 ) then
                  call Freq_Avg ( frequencies,            &
                    & FilterShapes(shapeInd)%FilterGrid,  &
                    & FilterShapes(shapeInd)%FilterShape, &
                    & h_frq(:,sv_q,sv_r), h_avg )
                else
                  h_avg = DACsStaging3 ( channels(c)%used, sv_q, sv_r, channels(c)%dacs )
                end if
              else
                h_avg = 0.0
              end if
              h(c,sv_q,sv_r) = h_avg
            end do              ! Grid loop 2
          end do                ! Grid loop 1
        end do                  ! Channel loop
      case ( 3, 7, 13 )         ! Not frequency averaging, or PFA alone; copy.
        do sv_q = grids%l_v(mol1-1)+1, grids%l_v(mol1)
          do sv_r = grids%l_v(mol2-1)+1, grids%l_v(mol2)
            if ( grids%deriv_flags(sv_q) .and. grids%deriv_flags(sv_r) ) then
              h(:,sv_q,sv_r) = min ( max ( h_frq(:,sv_q,sv_r), &
                            &   real(-huge(0.0_r4), rp ) ), &
                            &   real( huge(0.0_r4), rp ) )
            else
              h(:,sv_q,sv_r) = 0.0
            end if
          end do
        end do
      case default             ! Impossible
      end select               ! Frequency averaging or not
    end subroutine Frequency_Avg_Second_Derivative

  ! ................................  Frequency_Avg_Second_Derivs  .....
    subroutine Frequency_Avg_Second_Derivs ( Ptg_i, Combine )

      integer, intent(in) :: Ptg_i          ! Pointing index
      logical, intent(in) :: Combine        ! "Combine LBL and PFA"

      integer :: K, KK  ! Loop inductor
      integer :: UB     ! Upper bound for first dimension of h_..._frq

      ub = noFreqs
      if ( combine ) ub = max(ub,noUsedChannels)

      ! Frequency Average the atmospheric (ONLY) SECOND derivatives
      ! with the appropriate filter shapes
      if ( atmos_second_der ) then

        do k = 1, no_mol
          do kk = 1, no_mol

            if ( fwdModelConf%moleculeSecondDerivatives(k) .and. &
               & fwdModelConf%moleculeSecondDerivatives(kk) ) then
              call frequency_avg_second_derivative ( &
                &  grids_f, h_atmos_frq(:ub,:,:), h_atmos(:,ptg_i,:,:), k, kk, combine )
            else
              h_atmos(:, ptg_i, grids_f%l_v(k-1)+1:grids_f%l_v(k), &
                              & grids_f%l_v(kk-1)+1:grids_f%l_v(kk)) = 0.0
            end if

          end do                      ! Loop over major molecules
        end do                        ! Loop over major molecules

      end if                          ! Want derivatives for atmos

    end subroutine Frequency_Avg_Second_Derivs

  ! ..........................................  Frequency_Setup_1  .....
    subroutine Frequency_Setup_1 ( Tan_Press, Grids )

      ! Work out which pointing frequency grid we're going to need

      ! Code splits into two sections, one for when we're doing frequency
      ! averaging, and one when we're not.

      real(rp), intent(in) :: Tan_Press(:)
      integer, intent(out) :: Grids(:)

      integer :: I, K, Ptg_i, ShapeInd
      integer :: MaxNoPtgFreqs     ! Used for sizing arrays
      integer :: MinSuperset       ! Min. value of superset > 0
      real(rp) :: R1, R2           ! real variables for various uses
      integer :: Superset          ! Output from AreSignalsSuperset
      integer :: SV_Dim            ! Second dimension of K_*_FRQ, or zero
      integer, parameter :: DumpingDetails = -1 ! Used in case of trouble

      nullify ( whichPointingGrid )

      if ( fwdModelConf%do_freq_avg .and. fwdModelConf%anyLBL(sx) ) then

        minSuperset = huge(0)
        do i = 1, size(pointingGrids)
          superset = AreSignalsSuperset ( pointingGrids(i)%signals, &
            & fwdModelConf%signals, sideband=thisSideband )
          if ( superset >= 0 .and. superset <= minSuperset ) then
            minSuperset = superset
            whichPointingGrid => pointingGrids(i)
          end if
        end do
        if ( .not. associated(whichPointingGrid) ) then
          call HeadLine ( 'Trouble in FullForwardModel line-by-line' )
          call output( 'ForwardModel signals', advance='yes' )
          call Dump( fwdModelConf%signals, Details=DumpingDetails )
          do i=1, size(pointingGrids)
            call outputNamedValue ( 'ptg grid num', i )
            call Dump_Pointing_Grid(  pointingGrids(i), Details=DumpingDetails )
          enddo
          call output( 'Consider the above carefully. You may consider', advance='yes' )
          call output( '(a) Merging 2 or more pointing grids', advance='yes' )
          call output( '(b) Splitting the fwdmdl signals among separate models', advance='yes' )
          call Announce_Error ( &
               & "No matching pointing frequency grids." )
        endif

        ! Now we've identified the pointing grids.  Locate the tangent grid
        ! within it.
        call Hunt ( whichPointingGrid%oneGrid%height, &
                      & tan_press, grids, allowTopValue=.TRUE., nearest=.TRUE. )
        ! Work out the maximum number of frequencies
        maxNoPtgFreqs = 0
        do ptg_i = 1, no_tan_hts
          maxNoPtgFreqs = max ( maxNoPtgFreqs, &
            & Size(whichPointingGrid%oneGrid(grids(ptg_i))%frequencies) )
        end do

        min_ch_freq_grid =  huge(min_ch_freq_grid)
        max_ch_freq_grid = -huge(min_ch_freq_grid)
        do i = 1, noUsedChannels
          shapeInd = channels(i)%shapeInds(sx)
          if ( channels(i)%dacs == 0 ) then
            k = Size(FilterShapes(shapeInd)%FilterGrid)
            r1 = FilterShapes(shapeInd)%FilterGrid(1)
            r2 = FilterShapes(shapeInd)%FilterGrid(k)
          else
            k = Size(dacsFilterShapes(shapeInd)%filter%FilterGrid)
            r1 = dacsFilterShapes(shapeInd)%filter%FilterGrid(1)
            r2 = dacsFilterShapes(shapeInd)%filter%FilterGrid(k)
          end if
          min_ch_freq_grid = MIN(r1, r2, min_ch_freq_grid)
          max_ch_freq_grid = MAX(r1, r2, max_ch_freq_grid)
        end do

        if ( FwdModelConf%anyPFA(sx) ) &
          & call get_channel_centers ( thisSideband, channelCenters )

      else ! ----------------------------- Not frequency averaging -----

        noFreqs = noUsedChannels
        maxNoPtgFreqs = noUsedChannels
        call get_channel_centers ( thisSideband, channelCenters )
        frequencies => channelCenters

      end if ! ----------------- Either frequency averaging or not -----

      call allocate_test ( RadV, maxNoPtgFreqs, 'RadV', moduleName )

      if ( FwdModelConf%anyLBL(sx) ) then
        call allocate_test ( T_Script_LBL, max_c, maxNoPtgFreqs, 'T_Script_LBL', moduleName )
        call allocate_test ( tau_LBL%tau,  max_c, maxNoPtgFreqs, 'Tau_LBL%Tau', &
          & moduleName )
        call allocate_test ( tau_LBL%i_stop, maxNoPtgFreqs, 'Tau_LBL%I_Stop', &
          & moduleName )
      end if

      call allocate_test ( inc_rad_path, max_c, &
        & max(maxNoPtgFreqs,noUsedChannels), 'Inc_Rad_Path', moduleName )

      call allocate_test ( k_temp_frq, max(maxNoPtgFreqs,noUsedChannels), &
                         & merge(sv_t_len,0,temp_der), 'k_temp_frq', &
                         & moduleName )

      ! Allocate temporary derivative arrays with zero size if they're
      ! not needed.  We can't simply leave them unallocated since they're
      ! used as actual arguments to One_Frequency.  They are allocated with
      ! zero extent in second dimension if they're not needed, because the
      ! first dimension is subscripted.

      ! SV_Dim is calculated as it is here because if the various flags
      ! indicating requests for derivatives are false the VALUES component
      ! of the Grids_t structures might not be associated.  In those cases,
      ! invoking the MERGE intrinsic function might require accessing an
      ! undefined or disassociated pointer.
      sv_dim = 0

!      IF (hmag_der .OR. htheta_der .OR. hphi_der) sv_dim = SIZE(grids_mag%values)
      call allocate_test( k_magfield_frq, max(maxNoPtgFreqs,noUsedChannels), &
           & MERGE(3*sv_h_len,0,hmag_der), 'k_magfield_frq', &
           & moduleName )
      sv_dim = 0
      
      if ( atmos_der ) sv_dim = size(grids_f%values)
      call allocate_test ( k_atmos_frq, max(maxNoPtgFreqs,noUsedChannels), &
                         & sv_dim, 'k_atmos_frq', moduleName )

      sv_dim = 0
      if ( atmos_second_der ) sv_dim = size(grids_f%values)
      call allocate_test ( h_atmos_frq, max(maxNoPtgFreqs,noUsedChannels), &
                         & sv_dim, sv_dim, 'h_atmos_frq', moduleName )

      sv_dim = 0
      if ( spect_der_width ) sv_dim = size(grids_w%values)
      call allocate_test ( k_spect_dw_frq, max(maxNoPtgFreqs,noUsedChannels), &
                         & sv_dim, 'k_spect_dw_frq', moduleName )

      sv_dim = 0
      if ( spect_der_width_TDep ) sv_dim = size(grids_n%values)
      call allocate_test ( k_spect_dn_frq, max(maxNoPtgFreqs,noUsedChannels), &
                         & sv_dim, 'k_spect_dn_frq', moduleName )

      sv_dim = 0
      if ( spect_der_center ) sv_dim = size(grids_v%values)
      call allocate_test ( k_spect_dv_frq, max(maxNoPtgFreqs,noUsedChannels), &
                         & sv_dim, 'k_spect_dv_frq', moduleName )

    end subroutine Frequency_Setup_1

  ! ..........................................  Frequency_Setup_2  .....
    subroutine Frequency_Setup_2 ( GridFrequencies )

      ! Work out what frequencies we're using for
      ! frequency averaging case for this pointing

      real(r8), intent(in) :: GridFrequencies(:) ! from PointingGrids
      ! Include the VELOCITY shift correction in GridFrequencies!

      integer :: J, K, L, M

      j = -1
      k = SIZE(GridFrequencies)
      call purehunt ( min_ch_freq_grid, GridFrequencies, k, j, l )
      call purehunt ( max_ch_freq_grid, GridFrequencies, k, l, m )
      noFreqs = m - j + 1
      call allocate_test ( frequencies, noFreqs, "frequencies", moduleName )

      frequencies = GridFrequencies(j:m)

      if ( print_frq ) call dump ( frequencies, 'Frequencies' )

    end subroutine Frequency_Setup_2

  ! .............................................  Generate_TScat  .....
    include 'Generate_TScat.f9h'

  ! ........................................  Get_Channel_Centers  .....
    subroutine Get_Channel_Centers ( ThisSideband, ChannelCenters )

      integer, intent(in) :: ThisSideband ! +1 for USB, -1 for LSB
      real(r8), intent(out) :: ChannelCenters(:)

      integer :: Channel           ! Loop inductor and subscript
      integer :: Sig               ! Subscript for fwdModelConf%signals

      select case ( thisSideband )
      case ( -1, +1 ) ! OK
      case ( 0 )
        call Announce_Error ( &
          & 'Folded signal requested in non-frequency-averaged forward model' )
      case default
        call Announce_Error ( 'Bad value of signal%sideband' )
      end select

      do channel = 1, noUsedChannels
        sig = channels(channel)%signal
        channelCenters(channel) = firstSignal%lo + thisSideband * &
          & ( fwdModelConf%signals(sig)%centerFrequency + &
          &   fwdModelConf%signals(sig)%direction * &
          &   fwdModelConf%signals(sig)%frequencies(channels(channel)%used) )
      end do

    end subroutine Get_Channel_Centers

  ! ..............................................  One_Frequency  .....
    subroutine One_Frequency ( Ptg_i, Frq_i, Alpha_Path_c, Beta_Path_c,     &
      & C_Inds, Del_S, Del_Zeta, Do_GL, Frq, H_Path_C, Tan_Ht, IncOptDepth, &
      & P_Path, PFA, Ref_Corr, Sps_Path, Tau, T_Path_c, T_Script, Tanh1_c,  &
      & TT_Path_c, W0_Path_c, Z_Path, I_Start, I_End, Inc_Rad_Path, RadV,   &
      & dAlpha_dT_Path, H_Atmos_Frq, K_Atmos_Frq, K_Spect_dN_Frq,           &
      & K_Spect_dV_Frq, K_Spect_dW_Frq, K_Temp_Frq, K_Magfield_Frq )

      ! Having arguments instead of using host association serves two
      ! purposes:  The array sizes are implicit, so we don't need explicitly
      ! to mention them, and the pointer attribute gets stripped during the
      ! trip through the CALL statement -- hopefully thereby helping optimizers.
      use Comp_Eta_And_Sps_m, only: Comp_Sps
      use Comp_Sps_Path_Sparse_m, only: Comp_Sps_Path_Sparse
      use CS_ExpMat_m, only: CS_ExpMat
      use Do_T_Script_m, only: Two_D_T_Script, Two_D_T_Script_Cloud
      use Dump_0, only: Dump_2x2xn
      use D_T_Script_dTNP_m, only: DT_Script_dT_Sparse
      use Dump_Path_m, only: Dump_Path, SPS_List
      use Get_Beta_Path_m, only: Get_Beta_Path, Get_Beta_Path_Cloud, &
        & Get_Beta_Path_PFA, Get_Beta_Path_Polarized
      use Get_dAlpha_dF_m, only: Get_dAlpha_dF
      USE Get_D_Deltau_Pol_m, ONLY: Get_d_Deltau_Pol_dF, &
        & Get_d_Deltau_Pol_dT, Get_d_Deltau_Pol_dH, Get_d_Deltau_Pol_dTheta, &
        & Get_d_Deltau_Pol_dPhi
      use Hessians_m, only: d2Rad_Tran_dF2, Get_d2Alpha_dF2
      use Interpolate_Mie_m, only: Interpolate_Mie
      use Load_Sps_data_m, only:  Load_One_Item_Grid
      use L2PC_m, only: L2PC_T
      use MCRT_m, only: MCRT_Der
      USE Opacity_m, ONLY: Opacity, Test_Opacity
      use Path_Contrib_m, only: Path_Contrib
      use Physics, only: H_Over_K
      use Rad_Tran_m, only: Rad_Tran_Pol, dRad_Tran_dF, dRad_Tran_dT, &
        & dRad_Tran_dX
      use ScatSourceFunc, only: T_Scat, Interp_TScat
      use Tau_m, only: Get_Tau
      use TScat_Support_m, only: Get_dB_dT, Get_TScat, Get_TScat_Setup, &
        & Get_TScat_Teardown, Mie_Freq_Index

      integer, intent(in) :: Ptg_i        ! Pointing index
      integer, intent(in) :: Frq_i        ! Frequency loop index
      real(rp), intent(out) :: Alpha_Path_c(:) ! \sum Beta_Path * mixing ratio
      real(rp), intent(out) :: Beta_Path_c(:,:) ! path x species
      integer, intent(in) :: C_Inds(:)    ! Selectors from complete path to coarse path
      real(rp), intent(in) :: Del_S(:)    ! Integration lengths along path
      real(rp), intent(in) :: Del_Zeta(:) ! Integration lengths in Zeta coords
      logical, intent(out) :: Do_GL(:)    ! Where to do GL correction
      real(r8), intent(in) :: Frq         ! The frequency
      real(rp), intent(in) :: H_Path_C(:) ! Heights on coarse path
      real(rp), intent(in) :: Tan_Ht      ! Geometric tangent height, km,
                                          ! from equivalent Earth center
      real(rp), intent(out) :: IncOptDepth(:)  ! Incremental optical depth
      real(rp), intent(in) :: P_Path(:)   ! Pressures along complete path
      logical, intent(in) :: PFA          ! Are we doing PFA or not?
      real(rp), intent(in) :: Ref_Corr(:) ! Refraction correction
      real(rp), intent(inout) :: Sps_Path(:,:) ! Species on path
      type(tau_t), intent(inout) :: Tau   ! Optical depth, inout so as not to
                                          ! undefine components' association status
      real(rp), intent(in) :: T_Path_c(:) ! Temperature on coarse path
      real(rp), intent(out) :: T_Script(:)! Planck function, Delta_B in some notes
      real(rp), intent(out) :: TT_Path_c(:) ! TScat on coarse path
      real(rp), intent(out) :: Tanh1_c(:) ! tanh(frqhk/t_path_c)
      real(rp), intent(out) :: W0_Path_c(:) ! w0 on coarse path
      real(rp), intent(in) :: Z_Path(:)   ! -Log10(Pressures) along complete path
      integer, intent(in) :: I_Start      ! Start of coarse path integration
      integer, intent(inout) :: I_End     ! End of coarse path integration
      ! INC_RAD_Path is (out) if .not. PFA, and (inout) if PFA
      real(rp), intent(inout) :: Inc_Rad_Path(:) ! Incremental radiance along the path
      real(rp), intent(out) :: RadV       ! Radiance
      real(rp), intent(out), target :: dAlpha_dT_Path(:) ! dAlpha/dT on
                                          ! composite coarse & fine path
      real(rp), intent(out) :: H_Atmos_Frq(:,:) ! d2I/(dVMRq dVMRr), ptg.frq X vmr-SVq X vmr-SVr
      real(rp), intent(out) :: K_Atmos_Frq(:)  ! dI/dVMR, ptg.frq X vmr-SV
      real(rp), intent(out) :: K_Spect_dN_Frq(:) ! ****
      real(rp), intent(out) :: K_Spect_dV_Frq(:) ! ****
      real(rp), intent(out) :: K_Spect_dW_Frq(:) ! ****
      REAL(rp), INTENT(out) :: K_Temp_Frq(:)   ! dI/dT, ptg.frq X T-SV
      REAL(rp), INTENT(out) :: k_magfield_frq(:) ! Di/Dmagfield X h-sv * 3

      integer, pointer :: CG_Inds(:) ! Indices on coarse grid where GL needed
      real(rp), pointer :: dAlpha_dT_Path_c(:)
      real(r8) :: FrqHK              ! 0.5 * Frq * H_Over_K
      integer, pointer :: GL_Inds(:) ! Indices of GL points within combined
                                     ! coarse & fine path that are GL points for
                                     ! panels needing GL -- subset of f_inds.
      integer :: I                   ! Loop inductor and subscript
      integer :: I_Stop              ! Stop path integration before I_End
      integer :: J, L                ! Loop inductor and subscript
      integer :: Me = -1             ! String index for trace
      integer :: Mie_Frq_Index       ! Index of Frq in F_s in Read_Mie
      integer :: NCG                 ! Number of panels needing GL = Size(cg_inds)
      integer :: NGL                 ! Total # of GL points = Size(gl_inds)
      integer :: P_Stop              ! Where to stop in polarized case
      logical :: PFA_or_not_pol      ! PFA .or. .not. fwdModelConf%polarized
      real(rp) :: PhiWindowRadians
      integer :: RadInL2PC           ! Which TScat radiance in L2PC to use
      complex(rp) :: Rad_Pol(2,2)    ! polarized radiance output of mcrt for one freq and pointing
        ! (-1,:,:) are Sigma_-, (0,:,:) are Pi, (+1,:,:) are Sigma_+

      type (Vector_T) :: dX, TScat   ! Clones of Column, Row vectors from TScat L2PC

      type(L2PC_t), pointer :: L2PC  ! The selected TScat L2PC

      ! Cloud stuff
      logical, parameter :: Cld_Fine = .false.
      real(r8), parameter :: TOL = 1.D-190

      call Trace_Begin ( me, 'ForwardModel.One_Frequency=', index=frq_i, &
        & cond=toggle(emit) .and. levels(emit) > 4 )

      pfa_or_not_pol = pfa .or. .not. fwdModelConf%polarized

      do_gl = .false.

      ! Set up path quantities --------------------------------------

      phi_path_c => phi_path(1:npf:ngp1)
      dAlpha_dT_Path_c => dAlpha_dT_Path(1:npf:ngp1)

      ! Compute the sps_path for this Frequency for all frequency-dependent sps
      call comp_sps_path_sparse ( grids_f, frq, eta_zp, eta_fzp, &
                                & sps_path(1:npf,:), firstSignal%lo, thisSideband )

      associate ( sps_path_x => sps_path(1:npf:ngp1,:) )
        sps_path_c(i_start:i_end,:) = sps_path_x(i_start:i_end,:)
      end associate
      sps_path_c(:i_start-1,:) = 0.0
      sps_path_c(i_end+1:npc,:) = 0.0

      if ( pfa ) then
        call get_beta_path_PFA ( frq, frq_i, z_path, c_inds, t_path_c, &
          & beta_group, sx, vel_rel, sps_path, beta_path_c,            &
          & t_der_path_flags, dBeta_dT_path_c, dBeta_dw_path_c,        &
          & dBeta_dn_path_c, dBeta_dv_path_c, dBeta_dIWC_path_c )
      else
        frqhk = 0.5_r8 * frq * h_over_k    ! h nu / 2 k
        tanh1_c = tanh( frqhk / t_path_c ) ! tanh ( h nu / 2 k T )
        ! dTanh_dT = -h nu / (2 k T**2) 1/tanh1 d(tanh1)/dT
        if ( temp_der ) dTanh_dT_c(:npc) = &
            & frqhk / t_path_c**2 * ( tanh1_c - 1.0_rp / tanh1_c )
        call get_beta_path ( Frq, firstSignal%lo, p_path, t_path_c,      &
          &  tanh1_c, beta_group, sx, fwdModelConf%polarized, gl_slabs,  &
          &  c_inds, beta_path_c, t_der_path_flags,                      &
          &  dTanh_dT_c, vel_rel, dBeta_dT_path_c, dBeta_dw_path_c,      &
          &  dBeta_dn_path_c, dBeta_dv_path_c, dBeta_df_path_c,          &
          &  grids_f%where_dBeta_df, sps_path(:npf,:) )
!        PRINT '(a)','====print out the betas======'
!        PRINT *,SIZE(beta_path_c)
!        PRINT *,size(dbeta_df_path_c)
!        PRINT *,size(dbeta_dt_path_c)
      end if

      if ( FwdModelConf%incl_cld .and. .not. pfa ) then
        ! Compute Scattering source function based on temperature profile
        ! at all angles U for each temperature layer assuming a plane
        ! parallel atmosphere.

        if ( ptg_i == 1 ) then
        ! ??? Can this work?  On all pointings after the first one, ???
        ! ??? it uses the result for the first frequency.           ???
        ! ??? We have a different frequency grid for each pointing, ???
        ! ??? so is this the wrong idea in the first place?         ???

          call T_scat ( temp%values(:,inst), Frq, GPH%values(:,inst),      &
             & 10.0**(-temp%template%surfs), vmrArray(:,1), vmrArray(:,2), &
             & vmrArray(:,3),fwdModelConf%num_scattering_angles,           &
             & fwdModelConf%num_azimuth_angles,                            &
             & fwdModelConf%num_ab_terms, fwdModelConf%num_size_bins,      &
             & fwdModelConf%no_cloud_species,                              &
             & scat_src%values, scat_alb%values, cld_ext%values, Scat_ang )

        end if

        call load_one_item_grid ( grids_tscat, scat_src, maf, phitan, &
          & fwdModelConf, SetDerivFlags=.false., SetTscatFlag=.true. )

        call comp_sps ( grids_tscat, tan_pt_f, z_path(1:npf), &
                      & phi_path(1:npf), tscat_path(1:npf,:) )

        ! project Tscat onto LOS
        call interp_tscat ( tscat_path(1:npf,:), Scat_ang(:), &
          & phi_path(1:npf), tt_path_f(1:npf,:) )

        if ( .not. cld_fine ) then  ! interpolate onto gl grids along the LOS

          scat_alb%template = temp%template
          cld_ext%template  = temp%template

          call load_one_item_grid ( grids_salb,  scat_alb, maf, phitan, fwdModelConf, .false.)
          call load_one_item_grid ( grids_cext,  cld_ext,  maf, phitan, fwdModelConf, .false.)

          where ( abs(grids_salb%values) < TOL ) grids_salb%values = 0.0
          where ( abs(grids_cext%values) < TOL ) grids_cext%values = 0.0

          call comp_sps ( Grids_salb, tan_pt_f, z_path(1:npf), &
            &  phi_path(1:npf), salb_path(1:npf,:) )

          call comp_sps ( Grids_cext, tan_pt_f, z_path(1:npf), &
            &  phi_path(1:npf), cext_path(1:npf,:) )

          beta_path_cloud_c(1:npc) = cext_path(1:npf:ngp1,1)
          w0_path_c(1:npc) = salb_path(1:npf:ngp1,1)
          tt_path_c(1:npc) = tt_path_f(1:npf:ngp1,1)

        else

          ! cld_fine              re-compute cext and w0 along the LOS
          call get_beta_path_cloud ( Frq, t_path(1:npf), tt_path_f(1:npf,:), &
            &  c_inds, beta_path_cloud_c(1:npc), w0_path_c, tt_path_c,     &
            &  IPSD(1:npf),  WC(:,1:npf), fwdModelConf )

        end if

        do j = 1, npc ! Don't trust compilers to fuse loops
          alpha_path_c(j) = dot_product( sps_path_c(j,:), &
                                &        beta_path_c(j,:) )     &
                                & + beta_path_cloud_c(j)
          incoptdepth(j) = alpha_path_c(j) * del_s(j)
        end do

        ! Needed to compute inc_rad_path and by rad_tran_pol
        ! Don't restrict this to (i_start:i_end) because the end points
        ! at 1 and npc are special.
        call two_d_t_script_cloud ( t_path_c, tt_path_c, w0_path_c, &
          & spaceRadiance%values(1,1), frq, t_script, B(:npc) )

      else ! Not full cloud model

        !{ {\tt incoptdepth} is $\Delta \delta_{s\rightarrow s-1} =
        !  \int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
        !    \frac{\text{d}h}{\text{d}\zeta} \text{d} \zeta$ where
        !  $G(\zeta)$ is {\tt Alpha_path_c}, which is approximated
        !  here by the rectangle rule, \emph{viz.}
        !  $\Delta \delta_{i\rightarrow i-1} \approx G(\zeta_i) \delta s_i$.

        do j = i_start, i_end ! Don't trust compilers to fuse loops
          alpha_path_c(j) = dot_product( sps_path_c(j,:), &
                                       & beta_path_c(j,:) )
          incoptdepth(j) = alpha_path_c(j) * del_s(j)
        end do
        incoptdepth(:i_start-1) = 0.0  ! in case not integrating full path
        incoptdepth(i_end+1:npc) = 0.0 ! in case not integrating full path

        if ( fwdModelConf%useTScat .and. .not. pfa ) then
          include 'Using_TScat.f9h'
        else ! Not TScat model

          !{ Compute $\Delta B_{ij}$ for $j$ = Frq_i.  See page 42 in
          !  19 August 2004 ATBD JPL D-18130.
          !  T_Script and B needed to compute inc_rad_path by rad_tran_pol.
          !  Don't restrict this to (i_start:i_end) because the end points
          !  at 1 and npc are special.
          call two_d_t_script ( t_path_c(i_start:i_end), &
            & spaceRadiance%values(1,1), frq, &
            & t_script(i_start:i_end), &
            & B(i_start:i_end) )

        end if

      end if ! end of check cld

      IF ( .NOT. fwdModelConf%polarized ) THEN
        ! Determine where to use Gauss-Legendre for scalar instead of a trapezoid.

        call path_contrib ( incoptdepth, tan_pt_c, i_start, i_end, &
          &                 e_rflty, fwdModelConf%tolerance, do_gl )

      else ! Extra stuff for polarized case
           ! Can't be doing TScat, so process the whole path

        call get_beta_path_polarized ( frq, h, beta_group%lbl(sx), gl_slabs, &
             & c_inds, beta_path_polarized, dBeta_dT_polarized_path_c, &
             & dbeta_dh_polarized_path_c, ptg_i)

        ! We put an explicit extent of -1:1 for the first dimension in
        ! the hope a clever compiler will do better optimization with
        ! a constant extent.
        ! Add contributions from nonpolarized molecules 1/4 1/2 1/4
        ! to Alpha here.

!        PRINT '(a)','====print out the polarized betas======'
!        PRINT *,frq,h
!        PRINT *,npc
!        PRINT *,SIZE(beta_path_polarized)
!        PRINT *,size(dbeta_dt_polarized_path_c)
!        PRINT *,SIZE(dbeta_dh_polarized_path_c,1), &
!             & SIZE(dbeta_dh_polarized_path_c,2), &
        !             & SIZE(dbeta_dh_polarized_path_c,3)
!        if (ptg_i == 1) THEN
!        DO j = 1, npc
!           PRINT '(a)','species 1'
!           PRINT *,dbeta_dh_polarized_path_c(-1:1,j,1)
!           PRINT '(a)','species 2'
!           PRINT *,dbeta_dh_polarized_path_c(-1:1,j,2)
!           PRINT '(a)','species 3'
!           PRINT *,j
!           PRINT *,dbeta_dh_polarized_path_c(-1:1,j,3)
!           PRINT '(a)','species 4'
!           PRINT *,dbeta_dh_polarized_path_c(-1:1,j,4)
!        END DO
!        endif
        do j = 1, npc
          alpha_path_polarized_c(-1:1,j) = &
            & matmul( beta_path_polarized(-1:1,j,:), sps_path_c(j,:) ) *  &
            & tanh1_c(j) + 0.25 * alpha_path_c(j)
          alpha_path_polarized_c(0,j) = alpha_path_polarized_c(0,j) + &
            & 0.25 * alpha_path_c(j)
        end do

        ! Turn sigma-, pi, sigma+ into 2X2 matrix incoptdepth_pol
        call opacity ( ct(1:npc), stcp(1:npc), stsp(1:npc), &
             & alpha_path_polarized_c(:,1:npc), incoptdepth_pol(:,:,1:npc) )

!        call test_opacity ( ct(1:npc), stcp(1:npc), stsp(1:npc), &
!             & alpha_path_polarized_c(:,1:npc), &
!             & incoptdepth_pol_test(:,:,1:npc) )

!        IF (ptg_i == 30) then
!           PRINT '(a)','print results of my opacity routine'
!           DO i = 1,npc
!              PRINT *,i,incoptdepth_pol(:,:,i)
!              PRINT *,i,incoptdepth_pol_test(:,:,i)
!           end do
!     ENDIF
!        incoptdepth_pol = incoptdepth_pol_test

!         IF (ptg_i == 1 .and. frq_i == 50) THEN
!            PRINT *,'Some stuff after opacity'
!            PRINT *,'this frequency ',frq
!            PRINT *,npc
!            PRINT *,'tangent height ',tan_ht
!            PRINT *,'heights on path ',h_path_c(1:npc)
!            PRINT *,'temperature on path ',t_path_c(1:npc)
!            PRINT *,'size of del_s ',SIZE(del_s)
!            PRINT *,'path lengths ',del_s(1:npc)
!            PRINT *,'magnetic field ',h            
!            PRINT *,'alpha path polaized ',alpha_path_polarized(1,1:npc)     
!            PRINT *,'incoptdepth_pol ',incoptdepth_pol(1,1,1:npc)           
!        ENDIF
       

        ! We don't add unpolarized incremental optical depth to diagonal
        ! of polarized incremental optical depth because we added the
        ! scalar Alpha_path to the sigma-, pi and sigma+ parts of
        ! Alpha_Path_Polarized_c above.  If we did add it here, we would
        ! need 0.5 factors to scale unpolarized "power absorption" to
        ! get "field absorption"

        do j = 2, npc-1
          incoptdepth_pol(1,1,j) = - incoptdepth_pol(1,1,j) * del_s(j)
          incoptdepth_pol(2,1,j) = - incoptdepth_pol(2,1,j) * del_s(j)
          incoptdepth_pol(1,2,j) = - incoptdepth_pol(1,2,j) * del_s(j)
          incoptdepth_pol(2,2,j) = - incoptdepth_pol(2,2,j) * del_s(j)
        end do

        if ( print_incopt ) call dump_2x2xn ( incoptdepth_pol, 'IncoptDepth_Pol' )

        do j = 2, npc-1
          ! deltau_pol = exp(incoptdepth_pol)
          call cs_expmat ( incoptdepth_pol(:,:,j), deltau_pol(:,:,j) )
        end do

        if ( print_deltau_pol > -1 ) call dump_2x2xn ( deltau_pol, 'Deltau_Pol' )

        ! Determine where to do GL
        call path_contrib ( deltau_pol(:,:,1:npc), tan_pt_c, e_rflty, &
             & fwdModelConf%tolerance, do_gl )

      end if

      !{ We want $\Delta \delta_{s\rightarrow s-1} = \int_{\zeta_i}^{\zeta_{i-1}}
      ! G(\zeta) \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
      ! \text{d} \zeta$ where $G(\zeta)$ is {\tt Alpha_path_c}.  We start with a
      ! rectangular estimate $\Delta \delta_{s\rightarrow s-1} \approx\alpha(s_i)
      ! \, \delta s_i$.  Where this is sufficiently accurate, it is replaced with
      ! a trapezoidal estimate $(\alpha(s_i) + \alpha(s_{i-1}))\, \delta s_i$. 
      ! Where GL is needed, we compute $\Delta \delta_{s\rightarrow s-1} =
      ! G(\zeta_i)\, \delta s_i + \int_{\zeta_i}^{\zeta_{i-1}}
      ! \left(G(\zeta)-G(\zeta_i)\right) \frac{\text{d}s}{\text{d}h}
      ! \frac{\text{d}h}{\text{d}\zeta} \,\text{d} \zeta$.  This cancels the
      ! singularity in $\frac{\text{d}s}{\text{d}h}$ at the tangent point,
      ! assuming $\left(G(\zeta)-G(\zeta_i) \right) \frac{\text{d}s}{\text{d}h}
      ! \rightarrow 0$ at the tangent point.
      !
      ! At one time, for the trapezoidal estimate, we inscrutably used
      ! dsdz_c~*~del_z as an approximation for del_s (compare to computation of
      ! rectangular estimate above) but we don't do that any longer because
      ! $\frac{\text{d}s}{\text{d}h}$ is singular at the tangent, therefore
      ! $\frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}\, \delta z
      ! \neq \delta s$. We maybe could have used sum(dsdz_gw_path) as a GL
      ! estimate of $\int \frac{\text{d}s}{\text{d}h}
      ! \frac{\text{d}h}{\text{d}\zeta}\,{\text{d}\zeta}$, but del_s is simpler.

      if ( WrongTrapezoidal ) then
        do j = i_start+1, min(i_end,tan_pt_c)
          if ( .not. do_gl(j) ) &
            & incoptdepth(j) = incoptdepth(j) + &
              & ( alpha_path_c(j-1) - alpha_path_c(j) ) * dsdz_c(j-1) * del_zeta(j)
        end do
        do j = tan_pt_c+1, i_end-1
          if ( .not. do_gl(j) ) &
            & incoptdepth(j) = incoptdepth(j) + &
              & ( alpha_path_c(j+1) - alpha_path_c(j) ) * dsdz_c(j+1) * del_zeta(j)
        end do
      else
        ! Before the tangent point, del_s(j) is the path length from J-1 to J
        do j = i_start+1, min(i_end,tan_pt_c)
          if ( .not. do_gl(j) ) &
            & incoptdepth(j) = &
              & ( alpha_path_c(j-1) + alpha_path_c(j) ) * 0.5 * del_s(j)
        end do
        ! After the tangent point, del_s(j) is the path length from J to J+1
        do j = tan_pt_c+1, i_end-1
          if ( .not. do_gl(j) ) &
            & incoptdepth(j) = &
              & ( alpha_path_c(j+1) + alpha_path_c(j) ) * 0.5 * del_s(j)
        end do
      end if

      if ( dump_ds ) then
        write ( 42, '(5(a,i0))' ) &
          & 'dsdz_c(', i_start, ':', min(i_end,tan_pt_c)-1, ') * del_zeta(', &
          & i_start+1, ':', min(i_end,tan_pt_c), ') \', &
          & min(i_end,tan_pt_c)-1 - i_start + 1
        write ( 42, '(1p5g15.6)' ) &
          & dsdz_c(i_start:min(i_end,tan_pt_c)-1)*del_zeta(i_start+1:min(i_end,tan_pt_c))
        write ( 42, '(5(a,i0))' ) &
          & 'dsdz_c(', tan_pt_c+2,':',i_end,') * del_zeta(', &
          & tan_pt_c+1, ':', i_end-1, ') \', i_end - tan_pt_c - 1
        write ( 42, '(1p5g15.6)' ) dsdz_c(tan_pt_c+2:i_end)*del_zeta(tan_pt_c+1:i_end-1)
        write ( 43, '(2(a,i0),a,i0)' ) &
          & 'del_s(', i_start+1, ':', min(i_end,tan_pt_c), ') \', &
          & min(i_end,tan_pt_c) - i_start
        write ( 43, '(1p5g15.6)' ) 0.5*del_s(i_start+1:min(i_end,tan_pt_c))
        write ( 43, '(2(a,i0),a,i0)' ) &
          & 'del_s(', tan_pt_c+1,':',i_end-1,') \', i_end - tan_pt_c - 1
        write ( 43, '(1p5g15.6)' ) 0.5*del_s(tan_pt_c+1:i_end-1)
      end if

      ! Get indices for GL points only for panels that need GL, then copy
      ! temperature and mixing ratios only for those points to T_Path_f and
      ! Sps_Path_f
      do_gl(1:npc:npc-1) = .false.
      call get_GL_inds ( do_gl, tan_pt_c, gl_inds_b, ngl, cg_inds_b, ncg, where_gl )
      cg_inds => cg_inds_b(:ncg)
      gl_inds => gl_inds_b(:ngl)
      ! ngl is ng * count(do_gl)
      t_path_f(:ngl) = t_path(gl_inds)
      sps_path_f(:ngl,:) = sps_path(gl_inds,:)

      if ( pfa ) then

        call get_beta_path_PFA ( frq, frq_i, z_path, gl_inds, t_path_f(:ngl), &
          & beta_group, sx, vel_rel, sps_path, beta_path_f(:ngl,:),           &
          & t_der_path_flags, dBeta_dT_path_f, dBeta_dw_path_f,               &
          & dBeta_dn_path_f, dBeta_dv_path_f, dBeta_dIWC_path_f )

      else

        tanh1_f(1:ngl) = tanh( frqhk / t_path_f(:ngl) )
        ! dTanh_dT = -h nu / (2 k T**2) 1/tanh1 d(tanh1)/dT
        if ( temp_der ) &
          & dTanh_dT_f(1:ngl) = frqhk / t_path_f(1:ngl)**2 * &
            & ( tanh1_f(1:ngl) - 1.0_rp / tanh1_f(1:ngl) )

        ! The derivatives that get_Beta_path computes depend upon which
        ! derivative arrays are allocated, not which ones are present.
        ! This avoids having multiple paths through the code, each with a
        ! different set of optional arguments.

        call get_beta_path ( Frq, firstSignal%lo, p_path, t_path_f(:ngl),     &
          & tanh1_f(1:ngl), beta_group, sx, fwdModelConf%polarized, gl_slabs, &
          & gl_inds, beta_path_f(:ngl,:), t_der_path_flags, dTanh_dT_f,       &
          & vel_rel, dBeta_dT_path_f, dBeta_dw_path_f, dBeta_dn_path_f,       &
          & dBeta_dv_path_f, dBeta_df_path_f, grids_f%where_dBeta_df,         &
          & sps_path )

      end if ! .not. pfa

      do j = 1, ngl ! loop around dot_product instead of doing sum(a*b,2)
                    ! to avoid path-length array temps
        alpha_path_f(j) = dot_product( sps_path_f(j,:), beta_path_f(j,:) )
        alpha_path(gl_inds(j)) = alpha_path_f(j)
      end do

      if ( print_incopt ) then
        call sps_list ( fwdModelConf )
        call output ( tan_pt_c, before='Tan_pt_c = ' )
        call output ( tan_pt_f, before=' Tan_pt_f = ', advance='yes' )
        call dump ( beta_path_c(i_start:i_end,:), name="Beta_Path_C", lbound=i_start )
        call dump ( sps_path_c(i_start:i_end,:), name="SPS_Path_C", lbound=i_start )
        call dump ( alpha_path_c(i_start:i_end), name="Alpha_Path_C", lbound=i_start )
        call dump ( incoptdepth(i_start+1:i_end-1), name="Incoptdepth", lbound=i_start+1 )
        call dump ( del_s(i_start+1:i_end-1), name="Del_s", lbound=i_start+1 )
        call dump ( gl_inds(:ngl), name="GL_Inds -- fine path indices needing GL integrands" )
        call dump ( cg_inds, name="CG_Inds -- coarse path indices needing GL" )
      end if

      if ( .not. fwdModelConf%polarized ) then

      ! Compute SCALAR radiative transfer --------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TScat computation wants only incoptdepth(i_start+1:i_end-1),
!!!! but this breaks the gold brick
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call get_tau ( frq_i, gl_inds, cg_inds, i_start, e_rflty,  &
          & del_zeta, alpha_path, ref_corr(:i_end), incoptdepth(:i_end), &
          & tan_pt_c, dsdz_gw_path, tau )

        i_stop = tau%i_stop(frq_i) ! total_opacity(i_stop) < underflow for exp

        radV = 0.0
        if ( .not. pfa .or. .not. fwdModelConf%anyLBL(sx) ) then
          ! Not doing PFA, or doing PFA but haven't done LBL.
          ! Get incremental radiance and radiance from Tau and T_Script.
          ! Don't clobber them if doing PFA and already did LBL.  If
          ! doing LBL, inc_rad_path will be frequency averaged to give
          ! Rad_Avg_Path.
          do j = i_start, i_stop
            inc_rad_path(j) = t_script(j) * tau%tau(j,frq_i)
            radV = radV + inc_rad_path(j)
          end do ! j
        end if

        if ( pfa .and. iand(frq_avg_sel,7) == 7 ) then ! See Frequency_Average.
          ! Doing PFA and did LBL and need derivatives.  Multiply Inc_Rad_Path
          ! by Tau to combine LBL and PFA.  Then sum to give RadV.
          ! Inc_Rad_Path is channel(-averaged) LBL radiance. Remember, when
          ! doing PFA, Frq_I is a channel number, and tau is tau_PFA.  See
          ! wvs-026 and wvs-027.
          do j = i_start, i_stop
            inc_rad_path(j) = inc_rad_path(j) * tau%tau(j,frq_i)
            radV = radV + inc_rad_path(j)
          end do ! j
        end if

        inc_rad_path(:i_start-1) = 0.0
        inc_rad_path(i_stop+1:) = 0.0

        if ( print_incRad ) then
          call dump ( tau%tau(i_start:i_stop,frq_i), name="Tau", lbound=i_start )
          call dump ( inc_rad_path(i_start:i_stop), name="Inc_Rad_Path", lbound=i_start )
          call dump ( t_script(i_start:i_stop), name="T_Script", lbound=i_start )
          call output ( frq_i, before="RadV(" )
          call output ( radV, before=") = ", advance="yes" )
        end if

        if ( print_path > -1 ) then
          call dump_path ( print_path, fwdModelConf, &
            & i_start, i_stop, i_end, phi_path_c, z_path(1:npf:ngp1),       &
            & sps_path_c(:npc,:), beta_path_c, alpha_path_c, incoptdepth,   &
            & inc_rad_path, t_path_c, frq_i, tau, frq )
            if ( print_path > 0 ) stop
        end if

      else ! Polarized model; can't combine with PFA or TScat

      ! Compute POLARIZED radiative transfer -----------------------------

        i_stop = npc ! needed by drad_tran_df.  We don't compute a total opacity
                     ! black out point for polarized like we do for scalar.

        ! Get the corrections to integrals for layers that need GL for
        ! the polarized species.
        call get_beta_path_polarized ( frq, h,beta_group%lbl(sx), gl_slabs, &
             & gl_inds, beta_path_polarized_f, dBeta_dT_polarized_path_f, &
             & dbeta_dh_polarized_path_f, ptg_i)

        ! We put an explicit extent of -1:1 for the first dimension in
        ! the hope a clever compiler will do better optimization with
        ! a constant extent.
        ! Add contributions from nonpolarized molecules 1/4 1/2 1/4
        ! to Alpha here.
        DO j = 1, ngl
          alpha_path_polarized_f(-1:1,j) = &
            & matmul( beta_path_polarized_f(-1:1,j,:), &
            &         sps_path_f(j,:) ) * tanh1_f(j) &
            & + 0.25 * Alpha_path_f(j)
          alpha_path_polarized_f(0,j) = alpha_path_polarized_f(0,j) +               &
            & 0.25 * alpha_path_f(j)
          alpha_path_polarized(-1:1,gl_inds(j)) = alpha_path_polarized_f(-1:1,j)
       END DO
       
!        IF (ptg_i == 1 .and. frq_i == 50) THEN
!           PRINT *,'Some stuff before the radiative transfer'
!           PRINT *,'tan pt c ',tan_pt_c
!           PRINT *,'e_rflty ',e_rflty
!           PRINT *,'del_zeta ',del_zeta
!           PRINT *,'ref_corr ',ref_corr
!           PRINT *,'dsdz_gw_path ',dsdz_gw_path
!           PRINT *,'alpha_path_polarized ',alpha_path_polarized(-1,1:ngl)
!            PRINT *,'incoptdepth_pol ',incoptdepth_pol(1,1,1:npc)           
!            PRINT *,'deltau_pol ',deltau_pol(1,1,1:npc)
!            PRINT *,'p stop ',p_stop
!        ENDIF
           
        call rad_tran_pol ( tan_pt_c, gl_inds, cg_inds, e_rflty, del_zeta,   &
          & alpha_path_polarized(:,1:npf), ref_corr,                         &
          & incoptdepth_pol(:,:,1:npc), deltau_pol(:,:,1:npc), dsdz_gw_path, &
          & ct, stcp, stsp, t_script, dump_rad_pol, prod_pol(:,:,1:npc),     &
          & tau_pol(:,:,1:npc), rad_pol, p_stop )

        IF (ptg_i == 1) THEN
!           PRINT *,'Some stuff after the radiative transfer'
!            PRINT *,'incoptdepth_pol ',incoptdepth_pol(1,1,1:npc)           
!            PRINT *,'deltau_pol ',deltau_pol(1,1,1:npc)
            PRINT '(a,f12.4)','freq and Radiance ',frq
            PRINT *,rad_pol
!            PRINT *,'p stop ',p_stop
        ENDIF

        if ( p_stop < 0 ) then ! exp(incoptdepth_pol(:,:,-p_stop)) failed
          call output ( 'Exp(incoptdepth_pol(:,:,' )
          call output ( -p_stop )
          call output ( ') failed.  Value is', advance='yes' )
          call dump ( incoptdepth_pol(:,:,-p_stop), options='c' ) ! clean=.TRUE.
          call output ( thisSideband, before='thisSideband = ' )
          call output ( ptg_i, before=', ptg_i = ' )
          call output ( frq_i, before=', frq_i = ', advance='true' )
          call dump ( t_path_c(:npc), name='T_Path' )
          call dump ( ref_corr(:npc), name='Ref_Corr' )
          call dump ( n_path_c(:npc), name='N_Path' )
          call dump ( h_path_c(:npc), name='H_Path' )
          call Announce_Error ( 'exp(incoptdepth_pol) failed' )
        end if

        print_pol_rad = -1
!        PRINT *,'Print pol rad flag ',print_pol_rad
!        if ( print_pol_rad > -1 ) then
        if ( print_pol_rad > -1 .and. ptg_i == 1) then
          PRINT *,'Print tangent height ',ptg_i,tan_ht
          PRINT *,'Print frequency ',frq_i,frq
          call output ( 'Radiance' )
          do j = 1, 2
            do i = 1, 2
              call output ( i, before=' (' )
              call output ( j, before=',' )
              call output ( real(rad_pol(i,j)), before=' (' )
              call output ( aimag(rad_pol(i,j)), before=',', after=')' )
            end do
          end do
          if ( print_pol_rad > 0 ) call output ( ptg_i, before=' Pointing ' )
          if ( print_pol_rad > 1 ) call output ( frq_i, before=' Frequency ' )
          call newLine
        end if

        if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
          radV = real(rad_pol(1,1))
        else
          radV = real(rad_pol(2,2))
        end if

      end if

      ! Compute derivatives if needed ----------------------------------

      if ( atmos_der ) then

        ! On coarse grid points
        call get_dAlpha_df ( sps_path_c(:npc,:), beta_path_c(:npc,:), &
          &                  dBeta_df_path_c(:npc,:), Grids_f,      &
          &                  dAlpha_df_path_c(:npc,:) )

        ! On fine grid points for panels where GL is needed.  Data for these
        ! have been extracted from the composite path and stored contiguously.
        call get_dAlpha_df ( sps_path_f(:ngl,:), beta_path_f(:ngl,:), &
          &                  dBeta_df_path_f(:ngl,:), Grids_f,      &
          &                  dAlpha_df_path_f(:ngl,:) )

        ! Put dAlpha_df_path_f into their correct places in dAlpha_df_path
        dAlpha_df_path(gl_inds,:) = dAlpha_df_path_f(:ngl,:)

        if ( fwdModelConf%useTScat ) then ! TScat and not polarized

          call drad_tran_df ( gl_inds, del_zeta, Grids_f, eta_fzp, do_gl,   &
            &  del_s, ref_corr, dsdz_gw_path, inc_rad_path, dAlpha_df_path, &
            &  i_start, tan_pt_c, i_stop, k_atmos_frq, dB_df,               &
            &  tau=tau%tau(:,frq_i), alpha_path_c=alpha_path_c,             &
            &  beta_c_e=beta_c_e_path_c(:npc),                              &
            &  dBeta_c_a_dIWC=dBeta_c_a_dIWC_path_c(:npc),                  &
            &  dBeta_c_s_dIWC=dBeta_c_s_dIWC_path_c(:npc),                  &
            &  dTScat_df=dTScat_df, w0=w0_path_c(:npc) )

        else if ( fwdModelConf%polarized .or. atmos_second_der ) then
          ! ( polarized or Hessians ) and not TScat
          ! print *, 'Calling drad_tran_df polarized'

          call drad_tran_df ( gl_inds, del_zeta, Grids_f, eta_fzp, do_gl,   &
            &  del_s, ref_corr, dsdz_gw_path, inc_rad_path, dAlpha_df_path, &
            &  i_start, tan_pt_c, i_stop, k_atmos_frq, dB_df, d_delta_df )

        else ! ! not TScat and not polarized

          call drad_tran_df ( gl_inds, del_zeta, Grids_f, eta_fzp, do_gl,   &
            &  del_s, ref_corr, dsdz_gw_path, inc_rad_path, dAlpha_df_path, &
            &  i_start, tan_pt_c, i_stop, k_atmos_frq, dB_df )

        end if

        if ( atmos_second_der ) then

          call get_d2Alpha_df2 ( sps_path_c(:npc,:), beta_path_c(:npc,:), &
            &                  dBeta_df_path_c(:npc,:), Grids_f,      &
            &                  d2Alpha_df2_path_c(:npc,:) )

          call get_d2Alpha_df2 ( sps_path_f(:ngl,:), beta_path_f(:ngl,:), &
            &                  dBeta_df_path_f(:ngl,:), Grids_f,      &
            &                  d2Alpha_df2_path_f(:ngl,:) )


          call d2rad_tran_df2 ( gl_inds, del_zeta, Grids_f, eta_fzp, do_gl, &
            &  del_s, ref_corr, dsdz_gw_path, inc_rad_path,                 &
            &  d2Alpha_df2_path_c(:npc,:), d2Alpha_df2_path_f, i_start,     &
            &  tan_pt_c, i_stop, d_delta_df, d2_delta_df2, h_atmos_frq )

        end if ! atmos_second_der

        if ( .not. pfa_or_not_pol ) then ! polarized and not PFA

          ! This block of code is returning faulty results.
          ! Unfortunately it was coded in a way that does not parallel
          ! how vmr derivatives are coded for  non-polarized forward models
          ! which makes debugging even harder.
!           print *, 'VMR derivatives for polarized radiance.'
!           print *, shape(de_df(:,:,1:p_stop,:))
          ! VMR derivatives for polarized radiance.
          ! Compute DE / Df from D Incoptdepth_pol / Df and put it
          ! into DE_DF.
          call Get_D_Deltau_Pol_DF ( ct, stcp, stsp, Grids_f, tan_pt_c,        &
            &  beta_path_polarized(:,1:p_stop,:), tanh1_c(:npc),               &
            &  eta_fzp, sps_path, del_s, incoptdepth_pol(:,:,1:p_stop),        &
            &  ref_corr(1:p_stop), d_delta_df, de_df(:,:,1:p_stop,:) )

          ! Compute D radiance / Df from Tau, Prod, T_Script
          ! and DE / Df.
!           print *, 'Compute D radiance / Df.'
!           print *, shape(de_df(:,:,1:npc,:))
          call mcrt_der ( t_script, sqrt(e_rflty),             &
            & deltau_pol(:,:,1:npc), de_df(:,:,1:npc,:),       &
            & prod_pol(:,:,1:npc), tau_pol(:,:,1:npc), p_stop, &
            & tan_pt_c, d_rad_pol_df )

          if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
            ! print *, 'l_a'
            k_atmos_frq(:) = real(d_rad_pol_df(1,1,:))
          else
            ! print *, 'not l_a'
            k_atmos_frq(:) = real(d_rad_pol_df(2,2,:))
          end if
!           if ( any(abs(k_atmos_frq(:)) > 0.) ) &
!             & print *, 'non-zero'

        end if ! polarized and not PFA

      end if ! atmos_der

      
      if ( temp_der ) then

        dAlpha_dT_path_c(:npc) = sum( sps_path_c(:npc,:) *  &
                                      dBeta_dT_path_c(1:npc,:),dim=2 )
        dAlpha_dT_path_f(:ngl) = sum( sps_path_f(:ngl,:) * &
                                      dBeta_dT_path_f(1:ngl,:),dim=2 )
        ! Put dAlpha_dT_path_f, values for GL points only, back into the
        ! composite path
        dAlpha_dT_path(gl_inds) = dAlpha_dT_path_f(:ngl)

        ! get d Delta B / d T * d T / eta
        if ( fwdModelConf%useTScat ) then

          call get_dB_dT ( alpha_path_c, B(:npc), t_path_c(:npc), frq,    &
                         & beta_c_e_path_c(:npc), dAlpha_dT_Path_c(:npc), &
                         & dBeta_c_a_dT_path_c(:npc),                     &
                         & dBeta_c_s_dT_path_c(:npc), tt_path_c(:npc),    &
                         & dTScat_dT(:npc,:), w0_path_c(:npc),            &
                         & eta_zp_T, i_start, i_stop, d_t_scr_dt(1:npc,:) )

        else

          call dT_script_dT_sparse ( t_path_c(:npc), B(:npc), eta_zp_T,         &
                                   & i_start, i_stop, frq, d_t_scr_dt(1:npc,:), &
                                   & skip=ngp1 )

        end if

        if ( pfa_or_not_pol ) then

          call drad_tran_dt ( gl_inds, del_zeta, h_path_c, dh_dt_path, &
            & alpha_path(1:npf), dAlpha_dT_path(:npf), eta_zp_T,       &
            & del_s, ref_corr, tan_ht, tan_pt_f, do_gl, h_path(:npf),  &
            & t_path(:npf), dsdh_path, dhdz_gw_path, dsdz_gw_path,     &
            & d_t_scr_dt(1:npc,:), tau%tau(:npc,frq_i), inc_rad_path,  &
            & i_start, tan_pt_c, i_stop,  grids_tmp%deriv_flags,       &
            & pfa .and. iand(frq_avg_sel,7) == 7, k_temp_frq )

            ! pfa .and. iand(frq_avg_sel,7) == 7 means doing PFA and did LBL
            ! and need derivatives.

        else ! pol and not pfa

          ! Temperature derivatives for polarized radiance
          ! Compute DE / DT from D Incoptdepth_Pol / DT and put
          ! into DE_DT.

          dAlpha_dT_polarized_path_c(:,1:npc) = 0.0
          dAlpha_dT_polarized_path_f(:,1:ngl) = 0.0
          do j = 1, no_mol
            do l = -1, 1
              dAlpha_dT_polarized_path_c(l,1:npc) = &
            & dAlpha_dT_polarized_path_c(l,1:npc) + &
                & sps_path_c(:npc,j) * dBeta_dT_polarized_path_c(l,1:npc,j)
              dAlpha_dT_polarized_path_f(l,1:ngl) = &
            & dAlpha_dT_polarized_path_f(l,1:ngl) + &
                & sps_path_f(:ngl,j) * dBeta_dT_polarized_path_f(l,1:ngl,j)
            end do ! l
          end do ! j
          DO l = -1, 1
            dAlpha_dT_polarized_path_c(l,1:npc) =               &
              & dAlpha_dT_polarized_path_c(l,1:npc) * tanh1_c + &
              & alpha_path_polarized_c(l,:npc) * dTanh_dT_c(:npc)
            dAlpha_dT_polarized_path_f(l,1:ngl) =                     &
              & dAlpha_dT_polarized_path_f(l,1:ngl) * tanh1_f(:ngl) + &
              & alpha_path_polarized_f(l,:ngl) * dTanh_dT_f(:ngl)
            ! Put GL for panels needing it back into composite path
            dAlpha_dT_polarized_path(l,gl_inds) = &
              & dAlpha_dT_polarized_path_f(l,1:ngl)
         END DO
         IF (hmag_der) THEN
            
          dAlpha_dH_polarized_path_c(:,1:npc) = 0.0
          dAlpha_dH_polarized_path_f(:,1:ngl) = 0.0
          do j = 1, no_mol
            do l = -1, 1
              dAlpha_dH_polarized_path_c(l,1:npc) = &
            & dAlpha_dH_polarized_path_c(l,1:npc) + &
            & sps_path_c(:npc,j) * tanh1_c(1:npc) * &
            & dBeta_dH_polarized_path_c(l,1:npc,j)
              dAlpha_dH_polarized_path_f(l,1:ngl) = &
            & dAlpha_dH_polarized_path_f(l,1:ngl) + &
            & sps_path_f(:ngl,j) * tanh1_f(1:ngl) * &
            & dBeta_dH_polarized_path_f(l,1:ngl,j)
            end do ! l
          end do ! j
          DO l = -1, 1
            ! Put GL for panels needing it back into composite path
            dAlpha_dH_polarized_path(l,gl_inds) = &
              & dAlpha_dH_polarized_path_f(l,1:ngl)
          END DO

          call get_d_deltau_pol_dH ( ct, stcp, stsp,           &
               & alpha_path_polarized(:,1:npf),                          &
               & dAlpha_dH_polarized_path(:,1:npf), eta_zp_h, p_stop,    &
               & del_s, gl_inds, del_zeta, do_gl(1:p_stop),              &
               & dsdh_path, dhdz_gw_path, dsdz_gw_path,                  &
               & incoptdepth_pol(:,:,1:p_stop), ref_corr(1:p_stop),      &
               & grids_mag%deriv_flags, ptg_i, de_dhmag(:,:,1:p_stop,:) )
!           IF (ptg_i == 1) then
!             PRINT '(a)','Compute de_dh'
!             PRINT *,tan_pt_c,SIZE(de_dhmag,1),SIZE(de_dhmag,2),SIZE(de_dhmag,3),SIZE(de_dhmag,4),npc
!             PRINT '(a)','element 1'
!             DO l = 1, npc
!                PRINT *,de_dhmag(1,1,1:npc,1)
!             end do
!             PRINT '(a)','element 2'
!             DO l = 1, npc
!                PRINT *,de_dhmag(1,1,1:npc,2)
!                end do
!             PRINT '(a)','element 3'
!             DO l = 1, npc
!                PRINT *,de_dhmag(1,1,1:npc,3)
!                end do
!          endif
          call mcrt_der ( t_script, sqrt(e_rflty),             &
            & deltau_pol(:,:,1:npc), de_dhmag(:,:,1:npc,:),       &
            & prod_pol(:,:,1:npc), tau_pol(:,:,1:npc), p_stop, &
            & tan_pt_c, d_rad_pol_dh )
!          IF (ptg_i == 1) THEN
!            PRINT '(a)','Compute d_rad_pol_dh'
!            PRINT *,ptg_i
!            PRINT *,SIZE(d_rad_pol_dh,1),SIZE(d_rad_pol_dh,2),SIZE(d_rad_pol_dh,3)
!             PRINT '(a)','element 1'
!            PRINT *,real(d_rad_pol_dh(1,1,1))
!             PRINT '(a)','element 2'
!            PRINT *,real(d_rad_pol_dh(1,1,2))
!             PRINT '(a)','element 3'
!            PRINT *,real(d_rad_pol_dh(1,1,3))
!         ENDIF
!          IF (ptg_i == 55) THEN
!            PRINT '(a)','Compute d_rad_pol_dh'
!            PRINT *,ptg_i
!            PRINT *,SIZE(d_rad_pol_dh,1),SIZE(d_rad_pol_dh,2),SIZE(d_rad_pol_dh,3)
!             PRINT '(a)','element 1'
!            PRINT *,real(d_rad_pol_dh(1,1,1))
!             PRINT '(a)','element 2'
!            PRINT *,real(d_rad_pol_dh(1,1,2))
!             PRINT '(a)','element 3'
!            PRINT *,real(d_rad_pol_dh(1,1,3))
!         ENDIF
            
           if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
              k_magfield_frq(1:sv_h_len) = REAL(d_rad_pol_dh(1,1,:))
!             PRINT *,k_magfield_frq(:)
          else
             k_magfield_frq(1:sv_h_len) = REAL(d_rad_pol_dh(2,2,:))
          end if
       ENDIF
       
         IF (htheta_der) THEN
            
          CALL get_d_deltau_pol_dtheta ( st, ct, sp, cp, stcp, stsp,     &
               & alpha_path_polarized(:,1:npf), eta_zp_h, p_stop,        &
               & del_s, gl_inds, del_zeta, do_gl(1:p_stop),              &
               & dsdh_path, dhdz_gw_path, dsdz_gw_path,                  &
               & incoptdepth_pol(:,:,1:p_stop), ref_corr(1:p_stop),      &
               & grids_mag%deriv_flags, ptg_i, de_dtheta(:,:,1:p_stop,:) )
!           IF (ptg_i == 1) then
!             PRINT '(a)','Compute de_dh'
!             PRINT *,tan_pt_c,SIZE(de_dhmag,1),SIZE(de_dhmag,2),SIZE(de_dhmag,3),SIZE(de_dhmag,4),npc
!             PRINT '(a)','element 1'
!             DO l = 1, npc
!                PRINT *,de_dhmag(1,1,1:npc,1)
!             end do
!             PRINT '(a)','element 2'
!             DO l = 1, npc
!                PRINT *,de_dhmag(1,1,1:npc,2)
!                end do
!             PRINT '(a)','element 3'
!             DO l = 1, npc
!                PRINT *,de_dhmag(1,1,1:npc,3)
!                end do
!          endif
          call mcrt_der ( t_script, sqrt(e_rflty),             &
            & deltau_pol(:,:,1:npc), de_dtheta(:,:,1:npc,:),       &
            & prod_pol(:,:,1:npc), tau_pol(:,:,1:npc), p_stop, &
            & tan_pt_c, d_rad_pol_dtheta )
!          IF (ptg_i == 1) THEN
!            PRINT '(a)','Compute d_rad_pol_dh'
!            PRINT *,ptg_i
!            PRINT *,SIZE(d_rad_pol_dh,1),SIZE(d_rad_pol_dh,2),SIZE(d_rad_pol_dh,3)
!             PRINT '(a)','element 1'
!            PRINT *,real(d_rad_pol_dh(1,1,1))
!             PRINT '(a)','element 2'
!            PRINT *,real(d_rad_pol_dh(1,1,2))
!             PRINT '(a)','element 3'
!            PRINT *,real(d_rad_pol_dh(1,1,3))
!         ENDIF
!          IF (ptg_i == 55) THEN
!            PRINT '(a)','Compute d_rad_pol_dh'
!            PRINT *,ptg_i
!            PRINT *,SIZE(d_rad_pol_dh,1),SIZE(d_rad_pol_dh,2),SIZE(d_rad_pol_dh,3)
!             PRINT '(a)','element 1'
!            PRINT *,real(d_rad_pol_dh(1,1,1))
!             PRINT '(a)','element 2'
!            PRINT *,real(d_rad_pol_dh(1,1,2))
!             PRINT '(a)','element 3'
!            PRINT *,real(d_rad_pol_dh(1,1,3))
!         ENDIF
            
           if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
                k_magfield_frq(sv_h_len+1:2*sv_h_len) = &
              & REAL(d_rad_pol_dtheta(1,1,:))
!             PRINT *,k_magfield_frq(:)
          else
             k_magfield_frq(sv_h_len+1:2*sv_h_len) = &
                  & REAL(d_rad_pol_dtheta(2,2,:))
          end if
       ENDIF

       IF (hphi_der) THEN
          CALL get_d_deltau_pol_dPhi ( stcp, stsp,     &
               & alpha_path_polarized(:,1:npf), eta_zp_h, p_stop,        &
               & del_s, gl_inds, del_zeta, do_gl(1:p_stop),              &
               & dsdh_path, dhdz_gw_path, dsdz_gw_path,                  &
               & incoptdepth_pol(:,:,1:p_stop), ref_corr(1:p_stop),      &
               & grids_mag%deriv_flags, ptg_i, de_dphi(:,:,1:p_stop,:) )
          call mcrt_der ( t_script, sqrt(e_rflty),             &
            & deltau_pol(:,:,1:npc), de_dphi(:,:,1:npc,:),       &
            & prod_pol(:,:,1:npc), tau_pol(:,:,1:npc), p_stop, &
            & tan_pt_c, d_rad_pol_dphi )
            
           if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
                k_magfield_frq(2*sv_h_len+1:3*sv_h_len) = &
              & REAL(d_rad_pol_dphi(1,1,:))
          else
             k_magfield_frq(2*sv_h_len+1:3*sv_h_len) = &
                  & REAL(d_rad_pol_dphi(2,2,:))
          end if
       ENDIF
       
          call get_d_deltau_pol_dT ( ct, stcp, stsp, tan_pt_c,            &
            & t_path(:npf), alpha_path_polarized(:,1:npf),                &
            & dAlpha_dT_path(:npf), dAlpha_dT_polarized_path,             &
            & eta_zp_T, p_stop, del_s, gl_inds, del_zeta,                 &
            & do_gl(1:p_stop), dsdh_path, dhdz_gw_path, dsdz_gw_path,     &
            & incoptdepth_pol(:,:,1:p_stop), ref_corr(1:p_stop),          &
            & h_path(:npf), dh_dt_path, tan_ht, tan_pt_f,                 &
            & grids_tmp%deriv_flags, de_dt(:,:,1:p_stop,:) )

          ! Compute D radiance / DT from Tau, Prod, T_Script, D_T_Scr_dT
          ! and DE / DT.

          call mcrt_der ( t_script, sqrt(e_rflty),             &
            & deltau_pol(:,:,1:npc), de_dt(:,:,1:npc,:),       &
            & prod_pol(:,:,1:npc), tau_pol(:,:,1:npc), p_stop, &
            & tan_pt_c, d_rad_pol_dt, d_t_scr_dt(1:npc,:) )

          if ( radiometers(firstSignal%radiometer)%polarization == l_a ) then
            k_temp_frq(:) = real(d_rad_pol_dt(1,1,:))
          else
            k_temp_frq(:) = real(d_rad_pol_dt(2,2,:))
          end if

        end if ! pol and not pfa

      end if ! temp_der

      if ( spect_der ) then

        if ( fwdModelConf%polarized ) then
          call Announce_Error ( &
            & "Spectroscopic derivatives for Polarized species not implemented yet." )
        else if ( pfa ) then
          call Announce_Error ( &
            & "Spectroscopic derivatives for PFA not implemented yet." )
        else

          ! Spectroscopic derivative  wrt: W
          if ( spect_der_width ) &
            & call dRad_Tran_dX ( gl_inds, del_zeta, grids_w, eta_zp_w,    &
              &  sps_path, fwdModelConf%lineWidth%beta(sx),                &
              &  dBeta_dw_path_c, dBeta_dw_path_f, do_gl, del_s, ref_corr, &
              &  dhdz_gw_path, inc_rad_path, tan_pt_c, i_stop, k_spect_dw_frq )

          ! Spectroscopic derivative  wrt: N
          if ( spect_der_width_TDep ) &
            & call dRad_Tran_dX ( gl_inds, del_zeta, grids_n, eta_zp_n,    &
              &  sps_path, fwdModelConf%lineWidth_tDep%beta(sx),           &
              &  dBeta_dn_path_c, dBeta_dn_path_f, do_gl, del_s, ref_corr, &
              &  dhdz_gw_path, inc_rad_path, tan_pt_c, i_stop, k_spect_dn_frq )

          ! Spectroscopic derivative  wrt: Nu0
          if ( spect_der_center ) &
            & call dRad_Tran_dX ( gl_inds, del_zeta, grids_v, eta_zp_v,    &
              &  sps_path, fwdModelConf%lineCenter%beta(sx),               &
              &  dBeta_dv_path_c, dBeta_dv_path_f, do_gl, del_s, ref_corr, &
              &  dhdz_gw_path, inc_rad_path, tan_pt_c, i_stop, k_spect_dv_frq )

        end if

      end if

      call Trace_End ( 'ForwardModel.One_Frequency=', index=frq_i, &
        & cond=toggle(emit) .and. levels(emit) > 4 )

    end subroutine One_Frequency

  ! ...............................................  One_Pointing  .....
    subroutine One_Pointing ( Ptg_i, Vel_Rel, R_Eq, Tan_Phi, Tan_Loc, QTM,  &
      & Tan_Press, Est_SCGeocAlt, Path, S, F_and_V, Scat_Zeta, Scat_Phi, &
      & Scat_Ht, Xi, Scat_Index, Scat_Tan_Ht, Forward, Rev, Which )

      use Comp_Eta_And_Sps_m, only: Comp_Sps
      use Comp_Eta_Docalc_Sparse_m, only: Comp_Eta
      use Comp_Sps_Path_Sparse_m, only: Comp_1_Sps_Path_Sparse_No_Frq
      use Generate_QTM_m, only: QTM_Tree_t
      use Get_Chi_Angles_M, only: Get_Chi_Angles
      use GLNP, only: GW
      use Load_SPS_Data_M, only: Load_One_Item_Grid, Load_SPS_Data
      use Metrics_m, only: Height_Metrics, More_Metrics, Tangent_Metrics
      use Metrics_3D_m, only: Metrics_3D_QTM, Horizontal
      use Min_Zeta_m, only: Lower_Path_Crossings
      use Path_Representation_m, only: Path_t
      use Phi_Refractive_Correction_M, only: Phi_Refractive_Correction
      use QTM_Tangent_Metrics_m, only: QTM_Tangent_Metrics
      use Read_Mie_M, only: IWC_S, T_S
      use Refraction_M, only: Comp_Refcor, MaxRefraction, Refractive_Index
      use Slabs_sw_M, only: Get_GL_SLABS_Arrays
      use TScat_Support_m, only: Find_Scattering_Point

      integer, intent(in) :: Ptg_i     ! Tangent height pointing index
      real(rp), intent(in) :: Vel_Rel  ! Vel_z / speedOfLight
       real(rp), intent(in) :: R_EQ    ! Equivalent Earth Radius at true surface
      real(rp), intent(in), optional :: Tan_Phi     ! orbit angle at tangent, radians
      class(h_v_geod), intent(in), optional :: Tan_Loc   ! (lon,Lat,ht)
                                       ! degrees and km, for QTM case
      type(QTM_tree_t), intent(inout), optional :: QTM ! for QTM case; only
        ! QTM%Path_Vertices is changed here.
      real(rp), intent(in), optional :: Tan_Press     ! hPa, not zeta
      real(rp), intent(in), optional :: Est_SCGeocAlt ! Est S/C geocentric
        ! altitude /m, used only to compute chi angles for convolution
      type(Path_t), intent(inout), optional :: Path    ! Path from
        ! instrument to tangent.  Path%Lines(1,1) is (probably) the instrument
        ! location.  Path%Lines(1,2) is a vector in the direction from the
        ! instrument to the tangent. Path%Lines(1,2) is the tangent or
        ! reflection point.  Path%Lines(2,2) is a vector in the direction of the
        ! ray after the tangent, which is either the same as Path%Lines(2,1), or
        ! the direction of the reflected ray.  Present iff UsingQTM.
      type(s_QTM_t), allocatable, intent(out), optional :: S(:) ! Positions on
        ! the line of sight are Path%Lines(:,1) + s%s*Path%Lines(:,2).
        ! Present if UsingQTM.
      type(facets_and_vertices_t), intent(in), optional :: F_and_V ! Indices
                                       ! of facets and vertices under Path.
                                       ! Present only for QTM.

      ! The following are all present iff fwdModelConf%generateTScat
      real(rp), intent(in), optional :: Scat_Zeta   ! of scattering point
      real(rp), intent(in), optional :: Scat_Phi    ! of scattering point
      real(rp), intent(in), optional :: Scat_Ht     ! To check we hit the
                                                    ! scattering point
      real(rp), intent(in), optional :: Xi          ! Scattering angle
      integer, intent(out), optional :: Scat_Index  ! in coarse grid
      real(rp), intent(in), optional :: Scat_Tan_Ht ! Tangent height above earth
                                                    ! geometric surface, km, for
                                                    ! subsurface rays and
                                                    ! generateTScat
      logical, intent(in), optional :: Forward      ! For subsurface rays
                                                    ! and Generate_TScat
      logical, intent(in), optional :: Rev          ! Reverse the path
      integer, intent(in), optional :: Which        ! Which TScat call got us here?
                                                    ! Negative if we're trying to
                                                    ! trace an Earth-intersecting ray.

      integer :: Frq_I          ! Frequency loop index
      integer :: I              ! Do index
      integer :: I_start, I_end ! Boundaries of coarse path to use
      integer :: J              ! Do index
      integer :: Me = -1        ! String index for trace
      integer :: Me_Etc = -1    ! String index for trace
      real(rp) :: REQ_S  ! Equivalent Earth Radius at height reference surface
      real(rp) :: Tan_Ht ! Geometric tangent height, km, from equivalent Earth center
      real(rp) :: Tan_Ht_S ! Tangent height above 1 bar reference surface, km

      call Trace_Begin ( me, 'ForwardModel.Pointing=', index=ptg_i, &
        & cond=toggle(emit) .and. levels(emit) > 3 )

      ! Assume it's not an earth-intersecting ray and min zeta is at
      ! the tangent point
      e_rflty = 1.0

      call Trace_Begin ( me_Etc, 'ForwardModel.Metrics_Etc', &
        & cond=toggle(emit) .and. levels(emit) > 4 )

      if ( .not. same_facets ) then ! Same_Facets is false only if UsingQTM
        call load_one_item_grid ( grids_tmp, temp, fmStat%maf, phitan, &
          & fwdModelConf, setDerivFlags=.true., subset=F_and_V%vertices )

        call load_sps_data ( FwdModelConf, phitan, fmStat%maf, grids_f, &
          & subset=F_and_V%vertices )
      end if

      ! Compute the index in the pressure grids where the tangent is,
      ! assuming Zeta_Only if QTM.
      tan_ind_c = max(1,ptg_i-surfaceTangentIndex+1) ! On coarse grid
      tan_ind_f = (tan_ind_c-1) * ngp1 + 1           ! On Z_GLgrid

      nz_ig = nlvl
      nz_if = (nz_ig-1) * ngp1 + 1                ! On Z_GLgrid
      tan_pt_f = nz_if + 1 - tan_ind_f            ! fine path tangent index
      tan_pt_c = (tan_pt_f + Ng) / Ngp1           ! coarse path tangent index
      npc = 2 * tan_pt_c
      npf = 2 * tan_pt_f + ng

      ! Coarse path extraction indices, assuming Zeta_Only if QTM
      c_inds => c_inds_b(:npc)
      c_inds = (/ (i*Ngp1-Ng,i=1,npc) /) ! 1:npf:ngp1
      ! And fine path extraction indices, assuming Zeta_Only if QTM
      do_gl(1:npc:npc-1) = .false.; do_gl(2:npc-1) = .true.
      call get_gl_inds ( do_gl(:npc), tan_pt_c, f_inds, nglMax )

      ! Zetas on the coarse path, assuming Zeta_Only if QTM
      z_coarse(:tan_pt_c) = z_psig(nlvl:tan_ind_c:-1)
      z_coarse(tan_pt_c+1:npc) = z_psig(tan_ind_c:nlvl)

      ! Compute the height of the pressure reference surface and the
      ! tangent height above that surface.
      
      if ( .not. usingQTM ) then
         IF ( ASSOCIATED(surfaceHeight) ) THEN

          call tangent_metrics ( tan_phi, grids_tmp%phi_basis,      &
            &                    h_glgrid, tan_ind_f,               & ! in
            &                    h_surf, tan_ht_s,                  & ! output
            &                    surf_height=surfaceHeight%values(1,:) ) ! optional
        else if ( ptg_i < surfaceTangentIndex ) then
          call tangent_metrics ( tan_phi, grids_tmp%phi_basis,      &
            &                    h_glgrid, tan_ind_f,               & ! in
            &                    h_surf, tan_ht_s,                  & ! output
            &                    z_ref=z_psig(1),                   & ! optional
            &                    Tan_Press=tan_press,               & ! optional
            &                    Surf_Temp=temp%values(1,windowstart:windowfinish) )
        else
!            PRINT *,'I got here 3'
          call tangent_metrics ( tan_phi, grids_tmp%phi_basis,      &
            &                    h_glgrid, tan_ind_f,               & ! in
            &                    h_surf, tan_ht_s )                   ! output
        end if

        if ( present(scat_tan_ht) ) tan_ht_s = scat_tan_ht + h_surf

        ! Get H_Path and Phi_Path on the fine grid.
!             PRINT *,'I got here 4'
       call Height_Metrics ( tan_phi, tan_ind_f, grids_tmp%phi_basis, & ! in
          &  h_glgrid, h_surf, tan_ht_s, z_glgrid(:nz_if), r_eq,       & ! in
          &  req_s, vert_inds(1:npf), h_path(1:npf), phi_path(1:npf),  & ! out
          &  forward=forward )                                           ! opt
        tan_ht = tan_ht_s + req_s
        phi_path_c => phi_path(1:npf:ngp1)
        if (ptg_i == 1) then
        PRINT *,'tangent index and ptg_i ',surfaceTangentIndex,ptg_i
        PRINT *,'surface heights ',h_surf
!        PRINT *,'surf heights%values ',surfaceHeight%values(1,:)
        PRINT *,'equivalent earth radius ',r_eq
        PRINT *,'equivalent earth radius at surface ',req_s
        PRINT *,'tangent height ',tan_ht
        PRINT *,'tangent height surface ',tan_ht_s
        PRINT *,'size of path ',npf
!        PRINT *,'path heights ',h_path(1:npf)
        endif

      else ! QTM

        ! Put the inverse of F_and_V%Vertices in QTM%Path_Vertices.  The values
        ! therein are used to index elements of arrays that only exist
        ! adjacent to the path.
        QTM%path_vertices = 0 ! If this turns out to be expensive, keep an
                              ! "old F_and_V" and use it to put zeros here.
        do i_start = 1, size(F_and_V%vertices) ! I_Start is a temp here
          QTM%path_vertices(F_and_V%vertices(i_start)) = i_start
        end do

        ! Interpolate temperatures that are adjacent to the path onto T_GLgrid
        ! and compute H_GLgrid etc. at those places.  The second extent of
        ! T_GLgrid etc. is the maximum number of vertices adjacent to any path.
        ! The hydrostatic calculation for the non-QTM case is done before the
        ! pointing loop.  For the QTM case, each pointing might cross a different
        ! set of vertices, which would require T_GLgrid etc. to be larger,
        ! perhaps as large as the entire QTM.
        if ( temp_der ) then
          call two_d_hydrostatic ( grids_tmp, refGPH%template%surfs(1,1), &
            &  refGPH%values(1,:), z_glgrid, &
            &  t_glgrid, h_glgrid, dhdz_glgrid, eta_zzT, &
            &  dHidTlm=dh_dt_glgrid, ddHdHdTl0=ddhidhidtl0, &
            &  vertices=F_and_V%vertices )
        else
          call two_d_hydrostatic ( grids_tmp, refGPH%template%surfs(1,1), &
            &  refGPH%values(1,:), z_glgrid, &
            &  t_glgrid, h_glgrid, dhdz_glgrid, eta_zzT, &
            &  vertices=F_and_V%vertices )
        end if

        if ( associated(surfaceHeight) ) then
          call QTM_tangent_metrics ( tan_loc, QTM, h_glgrid, tan_ind_f, & ! in
            &                 h_surf, tan_ht_s, &                     ! output
            &                 surf_height=surfaceHeight%values(1,:) ) ! optional
        else if ( ptg_i < surfaceTangentIndex ) then
          call QTM_tangent_metrics ( tan_loc, QTM, h_glgrid, tan_ind_f, & ! in
            &                 h_surf, tan_ht_s, &                     ! output
            &                 Tan_Press=tan_press,               &    ! optional
            &                 Surf_Temp=temp%values(1,windowstart:windowfinish) )
        else
          call QTM_tangent_metrics ( tan_loc, QTM, h_glgrid, tan_ind_f, & ! in
            &                 h_surf, tan_ht_s )                      ! output
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!                                                         !!!!!
        !!!!! At this point Tan_Loc%v (geodetic tangent height        !!!!!
        !!!!! in meters) should equal 1000 * ( h_surf + tan_ht_s )    !!!!!
        !!!!!                                                         !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call path%get_path_ready ! Calculate tangent or reflection, and
                                 ! continuation of the path thereafter.
                                 ! The reflection is from the surface of the
                                 ! Earth, not the reference pressure surface.
        call metrics_3d_QTM ( path, QTM, h=h_glgrid, s=s, &
          & tangent_index=tan_ind_f, pad=NG, f_and_v=F_and_V, &
          & eta_p=eta_p, eta_p_t=eta_p_T, which=horizontal )
        ! ECR coordinates of points on the line-of-sight are
        ! Path%Lines(1,1) + S%s(:tan_ind_f) * Path%Lines(2,1) from the
        ! instrument to the tangent, and
        ! Path%Lines(1,2) + S%s(tan_ind_f+ngp1:) * Path%Lines(2,2) from the
        ! tangent onward. All values of |S%face| should be Top_Face, with
        ! S(.)%Face < 0 if the intersection is outside the QTM.  S(.)%H_ind is
        ! the subscript of z_psig, and the first subscript of H_GLGrid, i.e.,
        ! the index of the zeta surface in the fine zeta grid.  If it's zero,
        ! S(.) is an Earth-reflecting point below H_GLGrid.

      end if ! QTM

      ! Handle Earth-intersecting ray.  It is assumed to reflect from the
      ! reference height surface (h_glgrid(1,:)) so we didn't add a new
      ! zeta at the reflection point, so we don't need to recompute
      ! heights from hydrostatic principles.
      if ( tan_ht_s <= 0.0 ) e_rflty = earthRefl%values(1,1)

      ! If we're integrating from a scattering point, find where it is.
      if ( present(scat_zeta) .and. .not. usingQTM ) then
        call Find_Scattering_Point ( Scat_Zeta, Scat_Phi, Scat_Ht, Xi, &
                                   & Scat_Index, Tan_Phi, Npc, Npf, &
                                   & Print_Incopt, Print_IncRad, &
                                   & Print_Path, Req_S, PhiTan, Tan_Ht, &
                                   & FwdModelConf, MAF, Z_Coarse, Vert_Inds, &
                                   & H_Path, Phi_Path, Tan_Pt_C, Tan_Pt_F, &
                                   & Forward, Rev, Which )
        if ( scat_index < 1 ) then
          call trace_end ( 'ForwardModel.Metrics_Etc', &
            & cond=toggle(emit) .and. levels(emit) > 4 )
          call trace_end ( 'ForwardModel.Pointing=', index=ptg_i, &
            & cond=toggle(emit) .and. levels(emit) > 3 )
          return ! No ray to trace
        end if
        i_start = 1
        i_end = scat_index
      else
        i_start = 1
        i_end = npc
      end if

      ! Look for path crossings at zetas below the tangent point.
      ! These can only happen if the minimum zeta isn't at the tangent,
      ! and the ray isn't an earth-intersecting ray.
      if ( tan_ht_s > 0.0 .and. .not. present(QTM) ) call Lower_Path_Crossings ( &
        & Tan_Ht_s, Tan_Ht, Tan_Phi, Tan_Ind_f, grids_tmp%Phi_Basis, NZ_IF, &
        & Z_Coarse, Z_GLGrid, H_GLGrid, Vert_Inds, Req_S, H_Surf, Phi_Path, &
        & H_Path, NPC, NPF, T_GLGrid, Do_More_Points, Do_Zmin, &
        & Tan_Pt_C, Tan_Pt_F, More_Z_Path, More_H_Path, More_Phi_Path, &
        & N_More, Print_More_Points, Ptg_i )

      ! Compute Gauss Legendre (GL) zeta and pressure grids on the path
      call compute_GL_grid ( z_coarse(:tan_pt_c), z_path(:tan_pt_f), &
        &                    p_path(:tan_pt_f) )
      p_path(tan_pt_f+1:tan_pt_f+ng) = p_path(tan_pt_f) ! So they're not undefined
      z_path(tan_pt_f+1:tan_pt_f+ng) = z_path(tan_pt_f) ! So they're not undefined
      ! z_path, p_path, and t_path all the same within the zero-thickness
      ! tangent layer
      call compute_GL_grid ( z_coarse(tan_pt_c+1:npc), z_path(tan_pt_f+ngp1:npf), &
        &                    p_path(tan_pt_f+ngp1:npf) )

      ! The 0.5 factor is to compensate for the GL weights adding up to 2.0
      ! because the GL abscissae are on -1..+1.  This is a place where attention
      ! will be needed if present(QTM) .and. .not. Zeta_Only
      del_zeta(1:npc:npc-1) = 0.0_rp ! First and last ones
      del_zeta(2:tan_pt_c) = 0.5_rp * ( z_coarse(1:tan_pt_c-1) - &
        &                               z_coarse(2:tan_pt_c) )
      del_zeta(tan_pt_c+1:npc-1) = 0.5_rp * ( z_coarse(tan_pt_c+2:npc) - &
        &                                     z_coarse(tan_pt_c+1:npc-1) )

      ! Do phi refractive correction
      if ( FwdModelConf%refract .and. .not. usingQTM ) then
        ! Get t_path (and dhdz_path, which we don't need yet) so we can
        ! calculate the refractive index.  We don't do this for QTM because the
        ! horizontal coordinate isn't Phi, which is what's adjusted here.
        ! T_path and dhdz_path will be gotten again later, with slightly
        ! different phi.
        call more_metrics ( tan_ind_f, tan_pt_f, grids_tmp,            &
          &  vert_inds(1:npf), t_glgrid, dhdz_glgrid, phi_path(1:npf), &
          &  eta_p_T, t_path(1:npf), dhdz_path(1:npf) )
        ! Compute refractive index on the path.
        if ( h2o_ind > 0 ) then
          ! Compute eta_fzp (Zeta & Phi only) for water.
          call comp_eta ( grids_f, tan_pt_f, z_path(1:npf), &
                        & eta_z(h2o_ind), phi_path(1:npf),  &
                        & eta_p(h2o_ind), eta_fzp(h2o_ind), &
                        & n=h2o_ind, skip=ngp1 )
          call comp_1_sps_path_sparse_no_frq ( grids_f, eta_fzp(h2o_ind), &
                                             & sps_path(1:npf,h2o_ind),   &
                                             & n=h2o_ind )
          call refractive_index ( p_path(1:npf), t_path(1:npf), n_path_f(1:npf), &
            &  h2o_path=sps_path(1:npf, h2o_ind) )
        else
          call refractive_index ( p_path(1:npf), t_path(1:npf), n_path_f(1:npf) )
          n_path_c(:npc) = n_path_f(1:npf:ngp1)
        end if
        ! Do the refractive correction.  Use t_path to store the correction,
        ! since we're going to recompute t_path right away.
        ! The correction is zero at the tangent points and the zero-thickness
        ! layer between the tangent points.
        n_path_f(1:npf) = min ( n_path_f(1:npf), MaxRefraction )
        call phi_refractive_correction ( tan_pt_f, n_path_f(1:npf), &
          & h_path(1:npf), t_path(1:npf) )
        phi_path(:tan_pt_f) = phi_path(:tan_pt_f) - t_path(:tan_pt_f)
        phi_path(tan_pt_f+ngp1:npf) = phi_path(tan_pt_f+ngp1:npf) + &
                                    & t_path(tan_pt_f+ngp1:npf)
      end if

      ! Get other metrics-related quantities: t_path, dhdz_path, dh_dt_path...
      ! This also works for QTM because Eta_p_T and Eta_zp_T also work for QTM.
      if ( temp_der ) then
        call more_metrics ( tan_ind_f, tan_pt_f, Grids_tmp,              &
          &  vert_inds(1:npf), t_glgrid, dhdz_glgrid, phi_path(1:npf),   &
          &  eta_p_T, t_path(1:npf), dhdz_path(1:npf),                   &
          !  Stuff for temperature derivatives:
          &  ddHidHidTl0 = ddhidhidtl0, dHidTlm = dh_dt_glgrid,          &
          &  Z_Ref=z_glgrid,                                             &
          &  ddHtdHtdTl0 = tan_d2h_dhdt, dHitdTlm = dh_dt_path,          &
          &  dHtdTl0 = tan_dh_dt, eta_zp = eta_zp_T,                     &
          &  t_der_path_flags = t_der_path_flags )
      else
        call more_metrics ( tan_ind_f, tan_pt_f, Grids_tmp ,             &
          &  vert_inds(1:npf), t_glgrid, dhdz_glgrid, phi_path(1:npf),   &
          &  eta_p_T, t_path(1:npf), dhdz_path(1:npf) )
      end if

      h_path_c(1:npc) = h_path(1:npf:ngp1)
      t_path_c(1:npc) = t_path(1:npf:ngp1)

      ! Compute Eta_ZP and Eta_FZP for all species
      call comp_eta ( grids_f, tan_pt_f, z_path(1:npf), eta_z, &
                    & phi_path(1:npf), eta_p, eta_zp, eta_fzp, &
                    & skip=ngp1 )
      ! Compute sps_path for all those with no frequency component, especially
      ! to get WATER (H2O) contribution for refraction calculations.
      call comp_sps_path_sparse ( Grids_f, eta_zp, sps_path(1:npf,:), eta_fzp )

      ! Compute the refractive index - 1
      if ( h2o_ind > 0 ) then
        ! Even if we did the refractive correction we need to do this,
        ! because the refractive correction changes phi_path, which
        ! changes sps_path.

        call refractive_index ( p_path(1:npf:ngp1), &
          &  t_path_c(1:npc), n_path_c(1:npc),  &
          &  h2o_path=sps_path(1:npf:ngp1, h2o_ind) )
      else if ( .not. FwdModelConf%refract ) then
        ! If we didn't do the refractive correction, we haven't yet
        ! computed the refractive index.
        call refractive_index ( p_path(1:npf), t_path(1:npf), n_path_f(1:npf) )
        n_path_c(:npc) = n_path_f(1:npf:ngp1)
      end if

      n_path_c(1:npc) = min ( n_path_c(1:npc), MaxRefraction )

      call comp_eta ( grids_mag, tan_pt_f, z_path(1:npf), &
           & phi_path(1:npf), eta_zp_h)

!      IF (ptg_i == 1) THEN
!         PRINT '(a)','eta_zp_h%rows'
!         PRINT *,SIZE(eta_zp_h%rows)
!         PRINT *,eta_zp_h%rows
!         PRINT '(a)','eta_zp_h%cols'
!         PRINT *,SIZE(eta_zp_h%cols)
!         PRINT *,eta_zp_h%cols
!         PRINT '(a)','eta_zp_h%lbnd'
!         PRINT *,SIZE(eta_zp_h%lbnd)
!         PRINT '(a)','eta_zp_h%ubnd'
!         PRINT *,SIZE(eta_zp_h%ubnd)
!         PRINT '(a)','eta_zp_h%e'
!         PRINT *,SIZE(eta_zp_h%e)
!         PRINT *,eta_zp_h%e
!      endif         
      
                    
      ! Special path quantities for spectroscopy derivatives
      if ( size(fwdModelConf%lineCenter) > 0 ) then
        call comp_sps ( grids_v, tan_pt_f, z_path(1:npf), &
          & phi_path(1:npf), eta_zp_v, spect_v_path(1:npf,:) )
        lineCenter_ix => beta_group%lbl(sx)%spect_der_ix(lineCenter)
      end if
      if ( size(fwdModelConf%lineWidth) > 0 ) then
        call comp_sps ( grids_w, tan_pt_f, z_path(1:npf), &
          & phi_path(1:npf), eta_zp_w, spect_w_path(1:npf,:) )
        lineWidth_ix => beta_group%lbl(sx)%spect_der_ix(lineWidth)
      end if
      if ( size(fwdModelConf%lineWidth_TDep) > 0 ) then
        call comp_sps ( grids_n, tan_pt_f, z_path(1:npf), &
          & phi_path(1:npf), eta_zp_n, spect_n_path(1:npf,:) )
        lineWidth_TDep_ix => beta_group%lbl(sx)%spect_der_ix(lineWidth_TDep)
      end if

      ! Special path quantities for cloud model
      if ( fwdModelConf%Incl_Cld .or. fwdModelConf%useTScat ) then ! s_i == 1 or s_ts == 1 here

        call comp_sps ( grids_IWC, tan_pt_f, z_path(1:npf), phi_path(1:npf), &
                      & iwc_path(1:npf,:) )

        if ( fwdModelConf%Incl_Cld ) then
          !set some cloud parameters to zero
          WC(1,1:npf)=iwc_path(1:npf,1)
          WC(2,1:npf)=0.
          IPSD(1:npf)=1000
        end if
      end if

      if ( FwdModelConf%polarized ) then

        BLOCK
          use Magnetic_Field_On_Path_m, only: Magnetic_Field_On_Path

          ! Get the rotation matrix for the magnetic field.  Use the
          ! matrix for the MIF having ptan nearest to tan_press.  They are
          ! nearly identical anyway.
          mif = minloc(abs(tan_press - &
              &            ptan%values(:ptan%template%nosurfs,maf)),1)

          ! Dimensions of ECRtoFOV%value3 are ( chans, surfs,
          ! instances*cross angles ). Chans is actually 3x3, so we need to
          ! reform it.
          rot = reshape ( ECRtoFOV%value3(1:9,mif,maf), [ 3, 3 ] )

          call magnetic_field_on_path ( grids_mag, tan_pt_f, &
            & h_path(1:npf), r_eq, z_path(1:npf), phi_path(1:npf), &
            & rot, mag_path(1:npf,1:4) )
        END BLOCK
!        if ( tan_pt_f == 5 ) THEN
!        PRINT '(a)','The magnetic path values'
!        PRINT *,'npf ',npf,tan_pt_f
!        PRINT *,mag_path(1:npf,3)
!        PRINT *,mag_path(1:npf,1)
!        PRINT *,mag_path(1:npf,2)
!        PRINT *,mag_path(1:npf,4)
!        PRINT '(a)','ecr to fov matrix'
!        PRINT *,rot
!        PRINT *,'zeta, phi, values'
!        PRINT *,grids_mag%zet_basis(:)
!        PRINT *,grids_mag%phi_basis(:)
!        PRINT *,grids_mag%values(:)
!        endif
        ct => mag_path(1:npf,3)   ! cos(theta)
        stcp => mag_path(1:npf,1) ! sin(theta) cos(phi)
        stsp => mag_path(1:npf,2) ! sin(theta) sin(phi)
        h => mag_path(1:npf,4)    ! magnitude of magnetic field

! to compute the theta and phi derivatives we need more trig functions
! for the jacobian
        IF (htheta_der .OR. hphi_der) THEN
!           PRINT '(a)','enterring trig module'
! we will assume 0 <= theta <= 180 therefore sin(theta) >= 0.0
           st(1:npf) = SQRT(1.0 - mag_path(1:npf,3)**2)
           sp(1:npf) = mag_path(1:npf,2) / st(1:npf)
           cp(1:npf) = mag_path(1:npf,1) / st(1:npf)
           where(st < 0.001_rp) 
              sp = mag_path(:,2)
              cp = mag_path(:,1)
           END where
        ENDIF   
!         if ( tan_pt_f == 5 ) THEN
!        PRINT '(a)','New trig values'
!        PRINT *,'npf ',npf,tan_pt_f
!        PRINT *,st(1:npf),ct(1:npf),sp(1:npf),cp(1:npf),stsp(1:npf),stcp(1:npf)
!        endif
        if ( print_Mag > -1 ) then
          call output ( maf, before='ECR to FOV matrix for MAF ' )
          call output ( mif, before=' and MIF ', advance='yes' )
          call dump ( rot, '', options=clean )
          call dump ( h_path(1:npf), 'Path height (km)', options=clean )
          call dump ( h, 'H', options=clean )
          call dump ( ct, 'Cos(theta)', options=clean )
          call dump ( stcp, 'Sin(theta) Cos(phi)', options=clean )
          call dump ( stsp, 'Sin(theta) Sin(phi)', options=clean )
          if ( iand(print_Mag,1) > 0 ) &
            & call dump ( mag_path(1:npf,1:3), 'Mag_Path', options=clean )
        end if

      end if ! polarized

      if ( .not. usingQTM ) then
        if ( present(est_scGeocAlt) .and. .not. fwdModelConf%generateTScat .and. &
           & thisSideband == fwdModelConf%sidebandStart ) then
          ! Compute the pointing angles.  These are needed for antenna
          ! convolution, not for ray tracing.  They are the same for both
          ! sidebands, so there's no need to compute them twice.  We can't easily
          ! move these computations into the convolution code because they need
          ! Tan_Ht_s and Req_s, which are gotten from Tangent_Metrics and
          ! Height_Metrics.
          n_path_inst = 0.0
          if ( temp_der ) then
            ! Est_SCgeocAlt is in meters, but Get_Chi_Angles wants it in km.
            call get_chi_angles ( 0.001_rp*est_scGeocAlt, n_path_c(tan_pt_c), &
               & n_path_inst, tan_ht_s, tan_phi, req_s, 0.0_rp, &
               & ptg_angles(ptg_i), tan_dh_dt, tan_d2h_dhdt, dx_dt(ptg_i,:), &
               & d2x_dxdt(ptg_i,:) )
          else
            ! Est_SCgeocAlt is in meters, but Get_Chi_Angles wants it in km.
            call get_chi_angles ( 0.001_rp*est_scGeocAlt, n_path_c(tan_pt_c), &
               & n_path_inst, tan_ht_s, tan_phi, req_s, 0.0_rp, &
               & ptg_angles(ptg_i) )
          end if
        end if
      else ! Using QTM
        block
          use Get_Chi_Angles_3D_m, only: Get_Chi_Angles
          if ( .not. fwdModelConf%generateTScat .and. &
             & thisSideband == fwdModelConf%sidebandStart ) then
            n_path_inst = n_path_c(1)
            if ( temp_der ) then
              call get_chi_angles ( n_path_c(tan_pt_c), n_path_inst, &
                & Q_LOS(k)%lines(1,1)%xyz, Q_LOS(k)%lines(2,1)%xyz, & ! ECR & LOS
                & tan_ht_s, req_s, 0.0_rp, & ! Elev_Offset
                & ptg_angles(ptg_i), &
                & tan_dh_dt, tan_d2h_dhdt, dx_dt(ptg_i,:), d2x_dxdt(ptg_i,:) )
            else
              call get_chi_angles ( n_path_c(tan_pt_c), n_path_inst, &
                & Q_LOS(k)%lines(1,1)%xyz, Q_LOS(k)%lines(2,1)%xyz, & ! ECR & LOS
                & tan_ht_s, req_s, 0.0_rp, & ! Elev_Offset
                & ptg_angles(ptg_i) )
            end if
          end if
        end block
      end if
      ! Compute refractive correction and Del_s
      n_path_c(1:npc) = n_path_c(1:npc) + 1.0_rp

      call comp_refcor ( tan_pt_c, h_path_c(:npc), n_path_c(:npc), &
                    &    tan_ht, del_s(:npc), ref_corr(:npc), ier )

      if ( ier /= 0 ) fmStat%flags = ior(fmStat%flags,b_refraction)

      !{ Since $s = \sqrt{h^2-h_t^2}$, $\frac{ds}{dh} = \frac{h}s$. We
      ! need {\tt dsdh_path} on the fine grid for Gauss-Legendre or
      ! Gauss-Lobatto quadrature, and on the coarse grid except at the
      ! tangent point for trapezoidal quadrature and Gauss-Lobatto
      ! quadrature, so compute it everywhere except at the tangent point
      ! and in the zero-thickness tangent layer, where it has a pole. 
      ! Besides, it's probably faster not to use a vector subscript to
      ! restrict it to the fine grid.

      dsdh_path(:tan_pt_f-1) = h_path(:tan_pt_f-1) / &
        & ( sqrt(h_path(:tan_pt_f-1)**2 - tan_ht**2 ) )
      dsdh_path(tan_pt_f:tan_pt_f+ngp1) = 0.0
      dsdh_path(tan_pt_f+ngp1+1:npf) = h_path(tan_pt_f+ngp1+1:npf) / &
        & ( sqrt(h_path(tan_pt_f+ngp1+1:npf)**2 - tan_ht**2 ) )

      do i = 1, 2
        do j = 1, nglMax, ng ! Avoid a temp for (/ ( gw, j = 1, nglMax/ng ) /)
          dhdz_gw_path(f_inds(j:j+ng-1)) = dhdz_path(f_inds(j:j+ng-1)) * gw
        end do

        dsdz_gw_path(f_inds(:nglMax)) = dsdh_path(f_inds(:nglMax)) * &
          & dhdz_gw_path(f_inds(:nglMax))

        ! We need dsdz = ds/dh * dh/dz, not multiplied by GW, for
        ! trapezoidal quadrature on the coarse grid.
        if ( WrongTrapezoidal ) &
          & dsdz_c(:npc) = dsdh_path(1:npf:ngp1) * dhdz_path(1:npf:ngp1)

        ! Multiply dhdz_path by ds / ( sum( ds/dh dh/dz gw ) d zeta ) =
        ! del_s / ( sum (dsdz_gw_path) * del_zeta ), which ought to be 1.0

        if ( i == 2 .or. .not. fwdModelConf%do_path_norm ) exit
        do j = 1, npc-2 ! do_gl is always false at the ends
          dhdz_path(f_inds((j-1)*ng+1:j*ng)) = &
            & dhdz_path(f_inds((j-1)*ng+1:j*ng)) * del_s(j+1) / &
            & (sum(dsdz_gw_path(f_inds((j-1)*ng+1:j*ng)))*del_zeta(j+1))
        end do
      end do

      ! Compute ALL the slabs_prep entities over the path's GL grid for this
      ! pointing & mmaf:
      if ( temp_der ) then
        call get_gl_slabs_arrays ( p_path(1:npf), t_path(1:npf), &
          &  vel_rel, gl_slabs(1:npf,:), fwdModelConf%Do_1D, &
          &  spect_v_path, lineCenter_ix, &
          &  spect_w_path, lineWidth_ix, &
          &  spect_n_path, lineWidth_TDep_ix, &
          &  t_der_path_flags(1:npf) )
      else
        call get_gl_slabs_arrays ( p_path(1:npf), t_path(1:npf), &
          &  vel_rel, gl_slabs(1:npf,:), fwdModelConf%Do_1D, &
          &  spect_v_path, lineCenter_ix, &
          &  spect_w_path, lineWidth_ix, &
          &  spect_n_path, lineWidth_TDep_ix )
      end if

      if ( fwdModelConf%useTScat ) then
        ! Get interpolating coefficients for Mie tables onto path
        ! Make sure Mie tables have been loaded.
        if ( .not. allocated(t_s) ) call MLSMessage ( MLSMSG_Error, &
          & moduleName, "UseTScat requested but no Mie tables loaded" )
        ! Get interpolation coefficients for IWC
        call eta_IWC_path_c%eta ( iwc_s, iwc_path(1:npf:ngp1,1) )
        ! Get interpolation coefficients for Temperature
        call eta_T_path_c%eta ( t_s, t_path_c )
        ! Get interpolation coefficients for Temperature X IWC
        call eta_T_IWC_path_c%eta ( eta_T_path_c, eta_IWC_path_c )

        ! Indicate we don't have Mie tables for this path
        prev_Mie_frq_ind = -1
      end if

      call trace_end ( 'ForwardModel.Metrics_Etc', &
        & cond=toggle(emit) .and. levels(emit) > 4  )

!{{\bfseries Notations}
!
! $N_f$ is the number of points in the pointing frequency grid.\\
! $n$ is an index in the pointing frequency grid.\\
! $N_p$ is the number of points in the line-of-sight path.\\
! $i$ is an index in the line-of-sight path.  $c$ is a channel index.\\
! $s$ indicates a strong-line (LBL) result.
!   $w$ indicates a weak-line (PFA) result.\\
! $\delta I^\sigma_{iq} = \Delta B_{iq} \tau^\sigma_{iq}$ is the incremental
!  radiance contribution at the $i^{\text{th}}$ point along the line-of-sight
!  path, for $\sigma$ either $s$  or $w$ and $q$ either $c$ or $n$.
!
! There are four possible combinations of LBL, PFA and
! frequency-averaging.

!{{\bfseries Frequency averaged, LBL only}
!      \begin{equation*}\begin{split}
!       I^s_c = \sum_{n=1}^{N_f} \phi_{nc} \Delta\nu_{nc}
!             \sum_{i=1}^{N_p} \delta I^s_{in}
!       \text{\phantom{xxxx}}
! %
!       \frac{\partial I^s_n}{\partial x_k} =&
!        \sum_{i=1}^{N_p}
!         \tau^s_{in} \frac{\partial \Delta B_{in}}{\partial x_k}
!         - \delta I^s_{in}
!          \sum_{j=1}^i \frac{\partial \delta^s_{jn}}{\partial x_k}\\
! %
!       \frac{\partial I_c}{\partial x_k} =&
!       \frac{\partial I^s_c}{\partial x_k} =
!        \sum_{n=1}^{N_f} \phi_{nc} \Delta\nu_{nc}
!        \frac{\partial I^s_n}{\partial x_k}
!      \end{split}\end{equation*}

!{{\bfseries Monochromatic, LBL only}
!      \begin{equation*}
!       I^s_c = \sum_{i=1}^{N_p} \delta I^s_{ic}
!       \text{\phantom{xxxx}}
! %
!       \frac{\partial I^s_c}{\partial x_k} =
!        \sum_{i=1}^{N_p}
!         \tau^s_{ic} \frac{\partial \Delta B_{ic}}{\partial x_k}
!         - \delta I^s_{ic}
!          \sum_{j=1}^i \frac{\partial \delta^s_{jc}}{\partial x_k}
!       \text{\phantom{xxxx}}
! %
!       \frac{\partial I_c}{\partial x_k} =
!       \frac{\partial I^s_c}{\partial x_k}
!      \end{equation*}

!{{\bfseries Frequency averaged, PFA only}
!      \begin{equation*}
!       I_c = \sum_{i=1}^{N_p} \delta I^w_{ic}
!       \text{\phantom{xxxx}}
! %
!       \frac{\partial I_c}{\partial x_k} =
!        \sum_{i=1}^{N_p}
!         \tau^w_{ic} \frac{\partial \Delta B_{ic}}{\partial x_k}
!           - \delta I^w_{ic}
!            \sum_{j=1}^i \frac{\partial \delta^w_{jc}}{\partial x_k}
!       \text{\phantom{xxxx}}
!      \end{equation*}

!{{\bfseries Frequency averaged, LBL and PFA}
!      \begin{equation*}
! %
!       \overline{\delta I^s_{ic}} =
!        \sum_{n=1}^{N_f} \phi_{nc} \Delta\nu_{nc} \delta I^s_{in}
!       \text{\phantom{xxxx}}
! %
!       I_{ic} = \overline{\delta I^s_{ic}} \tau^w_{ic}
!       \text{\phantom{xxxx}}
! %
!       I_c = \sum_{i=1}^{N_p} I_{ic}
!       \text{\phantom{xxxx}}
! %
!       \frac{\partial I_c}{\partial x_k} =
!       \frac{\partial I^s_c}{\partial x_k}
!        - \sum_{i=1}^{N_p} I_{ic}
!          \sum_{j=1}^i \frac{\partial \delta^w_{jc}}{\partial x_k}
!      \end{equation*}

!{{\bfseries Program variables}
!
! \begin{tabular}{llll}
! $\Delta B_{in}$ is {\tt T_Script_LBL} &
! $\Delta B_{ic}$ is {\tt T_Script_PFA} &
! $\tau^s_{in}$ is {\tt Tau_LBL} &
! $\tau^w_{ic}$ is {\tt Tau_PFA}
! \\
! $\delta I^\sigma_{iq}$ is {\tt Inc_Rad_Path} &
! $\overline{\delta I^s_{ic}}$ is {\tt Rad_Avg_Path} &
! $I_{ic}$ is also {\tt Rad_Avg_Path} &
! \\
! $\sum_{i=1}^{N_p} \delta I^\sigma_{iq}$ is {\tt RadV} &
! $I_c$ or $I^s_c$ is {\tt Radiances} &
! $\frac{\partial I^\sigma_q}{\partial x_k}$ is {\tt K_}$x${\tt_FRQ} &
! $\frac{\partial I_c}{\partial x_k}$ is {\tt K_}$x$
! \\
      ! \end{tabular}
      
      if ( FwdModelConf%anyLBL(sx) ) then
        DO frq_i = 1, SIZE(frequencies)
          call one_frequency ( ptg_i, frq_i, alpha_path_c(:npc),              &
            & beta_path_c(:npc,:), c_inds(:npc), del_s(:npc), del_zeta(:npc), &
            & do_GL(:npc), frequencies(frq_i), h_path_c, tan_ht,              &
            & incoptdepth(:npc), p_path(:npf), pfaFalse, ref_corr(:npc),      &
            & sps_path(:npf,:), tau_lbl, t_path_c(:npc),                      &
            & t_script_lbl(:npc,frq_i), tanh1_c(:npc), tt_path_c(:s_i*npc),   &
            & w0_path_c(:max(s_i,s_ts)*npc), z_path(:npf), i_start, i_end,    &
            & inc_rad_path(:,frq_i), RadV(frq_i), dAlpha_dT_path(:npf),       &
            & H_Atmos_Frq(frq_i,:,:), K_Atmos_Frq(frq_i,:),                   &
            & K_Spect_dN_Frq(frq_i,:), K_Spect_dV_Frq(frq_i,:),               &
            & K_Spect_dW_Frq(frq_i,:), K_Temp_Frq(frq_i,:),                   &
            & K_Magfield_frq(frq_i,:) )
       END DO
        if ( print_TauL ) then
          call output ( thisSideband, before='Sideband ' )
          call output ( ptg_i, before=' Pointing ' )
          call dump ( tau_lbl, noFreqs, ' Tau_LBL:', scat_index )
          if ( present(scat_index) ) then
            call dump ( t_script_lbl(scat_index:npc-1,:noFreqs), 'T_Script_LBL' )
          else
            call dump ( t_script_lbl(2:npc-1,:noFreqs), 'T_Script_LBL' )
          end if
        end if
      end if

      ! Handle PFA molecules
      if ( FwdModelConf%anyPFA(sx) ) then
        if ( iand(frq_avg_sel,7) == 7 ) then ! LBL + PFA + Derivs
          ! For every channel, frequency average the incremental radiance at
          ! every point along the path, giving Rad_Avg_Path for every channel
          ! and every point along the path.
          call frequency_avg_path ( frequencies, &
            & inc_rad_path(:npc,:), rad_avg_path(:npc,:) )
          ! Multiply by Tau_PFA to combine PFA contribution in One_Frequency.
            call frequency_average_derivatives ( ptg_i, combine=.false. )
            if ( atmos_second_der ) call frequency_avg_second_derivs ( ptg_i, .false. )
        end if

        do frq_i = 1, size(channelCenters)
          call one_frequency ( ptg_i, frq_i, alpha_path_c(:npc),                &
            & beta_path_c(:npc,:), c_inds(:npc), del_s(:npc), del_zeta(:npc),   &
            & do_GL(:npc),  channelCenters(frq_i), h_path_c, tan_ht,            &
            & incoptdepth(:npc), p_path(:npf), pfaTrue, ref_corr(:npc),         &
            & sps_path(:npf,:), tau_pfa, t_path_c(:npc),                        &
            & t_script_pfa(:npc,frq_i), tanh1_c(:npc), tt_path_c(:s_i*npc),     &
            & w0_path_c(:max(s_i,s_ts)*npc), z_path(:npf), i_start, i_end,      &
            & rad_avg_path(:,frq_i), RadV(frq_i), dAlpha_dT_path(:npf),         &
            & H_Atmos_Frq(frq_i,:,:), K_Atmos_Frq(frq_i,:),                     &
            & K_Spect_dN_Frq(frq_i,:), K_Spect_dV_Frq(frq_i,:),                 &
            & K_Spect_dW_Frq(frq_i,:), K_Temp_Frq(frq_i,:), K_Magfield_frq(frq_i,:) )
        end do
        if ( print_TauP ) then
          call output ( thisSideband, before='Sideband ' )
          call output ( ptg_i, before=' Pointing ' )
          call dump ( tau_pfa, noUsedChannels, ' Tau_PFA:' )
          call dump ( t_script_pfa(2:npc-1,:), 'T_Script_PFA' )
        end if

      end if ! FwdModelConf%anyPFA(sx)

      if ( print_Rad > 0 ) then
        print '("Pointing ", i0)', ptg_i
        call dump ( frequencies, name='frequencies', options='c' )
        call dump ( radV, name='radV', options='c' )
      end if

      ! Frequency average, or store, or combine LBL and PFA
      call frequency_average ( ptg_i )

      ! If we're doing frequency averaging, there's a different frequency
      ! grid for each pointing, but we don't need to deallocate it here
      ! because the allocate_test in frequency_setup_2 will deallocate it.

      call trace_end ( 'ForwardModel.Pointing=', index=ptg_i, &
        & cond=toggle(emit) .and. levels(emit) > 3 )

      ! End of pointing loop -------------------------------------------------
    end subroutine One_Pointing

  ! ..............................................  One_TScat_Ray  .....
    include 'One_TScat_Ray.f9h'

  end subroutine FullForwardModelAuto

! ------------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id: FullForwardModel_m.f90,v 2.411 2023/07/06 18:34:10 pwagner Exp $"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module FullForwardModel_m

! $Log: FullForwardModel_m.f90,v $
! Revision 2.411  2023/07/06 18:34:10  pwagner
! Commented-out most uncondtional printing
!
! Revision 2.410  2023/06/23 20:45:07  pwagner
! In middle of debugging pol fwdmdl
!
! Revision 2.409  2020/08/28 21:41:58  vsnyder
! Set up to calculate chi angles for QTM
!
! Revision 2.408  2020/05/05 23:57:56  vsnyder
! Made F_and_V_MIF dummy argument of FullForwardModelAuto allocatable. Intel
! ifort 19 gets a seg fault at the call if the actual argument is not allocated,
! which it isn't if the run isn't QTM.
!
! Revision 2.407  2020/04/22 01:59:27  vsnyder
! Move TScat stuff to includes. Some work on QTM chi angles
!
! Revision 2.406  2020/02/07 01:11:08  pwagner
! Offers advice and sympathy if fwdmdl config fails superset test
!
! Revision 2.405  2019/10/07 20:05:53  vsnyder
! Get WrongTrapezoidal from the config
!
! Revision 2.404  2019/09/05 17:10:33  pwagner
! Dont try to deallocate F_and_V_MIF unless allocated already
!
! Revision 2.403  2019/06/24 23:28:16  pwagner
! Updated to reflect TA-01-143
!
! Revision 2.402  2018/10/30 23:14:36  vsnyder
! Make sure RadV always has a value in One_Frequency
!
! Revision 2.401  2018/10/26 22:05:47  vsnyder
! Don't do the trapezoidal update beyond I_End if it's before the tangent.
!
! Revision 2.400  2018/09/12 22:51:16  vsnyder
! Changed name of dRad_Tran_dX_Sparse to dRad_Tran_dX
!
! Revision 2.399  2018/09/12 22:07:42  vsnyder
! Convert interpolators for full cloud forward model from dense to sparse.
! Convert interpolators for spectroscopy derivatives from dense to sparse.
! Use Comp_Sps.  use dRad_Tran_dX_Sparse instead of dRad_Tran_dX.  Use
! Comp_Eta instead of Comp_Eta_DoCalc_Sparse (because we no longer calculate
! Do_Calc -- which Comp_Eta_DoCalc_Sparse didn't calculate anyway).
!
! Revision 2.398  2018/09/05 20:56:10  vsnyder
! Use sparse interpolation for magnetic field on path
!
! Revision 2.397  2018/08/28 22:17:53  vsnyder
! Add WrongTrapezoidal with the value .true. to indicate the incorrect
! trapezoid rule is used to calculate incOptDepth.  Some changes because of
! rearranged argument lists in Comp_Sps_Path_Sparse_m and
! Comp_Eta_DoCalc_Sparse_m.
!
! Revision 2.396  2018/08/15 01:18:50  vsnyder
! Get S_QTM_t from Metrics_2D_m instead of QTM_Interpolation_Weights_3D_m.
! Eliminate dSdZ_C.  Eliminate Do_Clac_T and Do_Calc_T_1.  Get T_Der_Path_Flags
! from More_Metrics.
!
! Revision 2.395  2018/05/24 03:24:36  vsnyder
! Use sparse representation for dh_dt_path
!
! Revision 2.394  2018/05/17 02:15:45  vsnyder
! Use sparse instead of dense interpolation
!
! Revision 2.393  2018/05/15 03:26:25  vsnyder
! Change Mie tables from pointer to allocatable
!
! Revision 2.392  2018/05/14 23:40:58  vsnyder
! Change to sparse eta representation
!
! Revision 2.390  2017/11/01 00:11:32  vsnyder
! Comment out a bunch of testing stuff, soon to disappear entierly
!
! Revision 2.389  2017/11/01 00:09:32  vsnyder
! Remove recently-added ScECR because it's not needed
!
! Revision 2.388  2017/10/31 23:49:35  vsnyder
! Make Coefficients a parameterized type
!
! Revision 2.387  2017/10/31 17:36:46  vsnyder
! Change QTM path from type MIFLOS vector quantity to Path_T
!
! Revision 2.386  2017/09/20 01:18:49  vsnyder
! Change name of Comp_Eta_DoCalc_Sparse to Comp_Eta_DoCalc_List.  Change
! names of Eta_p and Eta_z to Eta_p_List and Eta_z_list.  Get some state
! vector quantities only if they're needed.
!
! Revision 2.385  2017/08/10 00:15:53  vsnyder
! Add Tan_Pt_C in calls to Get_d_Deltau_Pol_df.  Do trapezoidal update using
! Del_s instead of Del_z because ds/dh is singular at the tangent.  More
! stuff getting ready for QTM.
!
! Revision 2.384  2017/03/31 00:49:43  vsnyder
! Use F_and_V to map to Jacobian for QTM.  Make F_and_V scalar in
! One_Pointing.  Cosmetic changes.
!
! Revision 2.383  2017/03/24 00:10:43  vsnyder
! Change the first extent of nz_d_delta_df back to max_f because the number
! of nonzeroes in a column of d_delta_df is not limited unless the vertical
! grids are the same for all species.
!
! Revision 2.382  2017/03/20 23:25:40  vsnyder
! Double row dimension of Nz_d_Delta_df
!
! Revision 2.381  2017/03/11 00:58:13  vsnyder
! Many changes leading to 3D/QTM model
!
! Revision 2.380  2017/01/14 02:58:48  vsnyder
! Inching toward 3D QTM forward model
!
! Revision 2.379  2016/12/02 02:04:50  vsnyder
! Use 'P' Eta list for Eta_ZZ
!
! Revision 2.378  2016/11/23 21:35:13  vsnyder
! Compute Eta_ZZ to interpolate from Temperature's zeta basis to the GL
! zeta grid, and use it for hydrostatic calculations.
!
! Revision 2.377  2016/11/23 00:14:36  vsnyder
! Use types from Indexed_Values_m.  Some cannonball polishing.
!
! Revision 2.376  2016/11/17 01:45:26  vsnyder
! Use Comp_Sps_Path_No_Frq to get H2O for phi refractive correction.  Some
! work on QTM also.
!
! Revision 2.375  2016/11/14 21:10:47  vsnyder
! Change scat zeta tolerance test from relative to absolute
!
! Revision 2.374  2016/11/14 19:17:12  vsnyder
! Change scat zeta tolerance to 'half the bits match'
!
! Revision 2.373  2016/11/12 01:42:32  vsnyder
! Put the inverse of F_and_V%Vertices into QTM_Tree%Path_Vertices.  Replace
! Facets argument to Metrics_3D with F_and_V.
!
! Revision 2.372  2016/11/11 02:06:27  vsnyder
! For QTM, get LOS from MIF LOS before FullForwardModelAuto.  Use it to
! calculate facets and vertices on all paths.  Use maximum number of vertices
! on any path as the horizontal extent for temperature-related arrays.  Pass
! Vertices array to Two_D_Hydrostatic.  Don't call Two_D_Hydrostatic before
! the loop over sidebands.  Pass RefGPH values as m instead of km because
! Two_D_Hydrostatic does the conversion now.  Other stuff for inching toward
! 3D QTM model.
!
! Revision 2.371  2016/11/09 00:38:10  vsnyder
! Move Interpolate_MIF_to_Tan_Press to a separate model.  Get list of facets
! and pass it to Metrics_3D.  Other stuff for inching toward 3D model.
!
! Revision 2.370  2016/11/03 19:11:47  vsnyder
! Inching toward 3D forward model
!
! Revision 2.369  2016/10/25 22:27:39  vsnyder
! Inching toward 3D QTM
!
! Revision 2.368  2016/10/24 22:20:32  vsnyder
! A bunch of stuff for QTM etc.
!
! Revision 2.367  2016/06/03 23:47:23  vsnyder
! Make Tan_Press allocatable instead of a pointer.  Provide minor-frame
! quantities for "orbit" inclination, minor axis length, phiTan, and tangent
! height if the hGrid for temperature is QTM.  Ensure either temperature and
! all molecules are QTM, or none of them are.
!
! Revision 2.366  2016/05/02 23:32:52  vsnyder
! Get temperature quantity from FwdModelConf, remove some unused stuff
!
! Revision 2.365  2016/04/21 02:00:12  vsnyder
! Move Convolution and Convolution_Setup to Convolution_m
!
! Revision 2.364  2016/03/25 02:02:37  vsnyder
! Add a dump for polarized incremental optical depth just before its
! exponential is attempted to be computed.
!
! Revision 2.363  2016/02/25 00:57:58  vsnyder
! Include bounds of t_path_c in call to dt_script_dt
!
! Revision 2.362  2016/01/23 02:55:24  vsnyder
! Add printing for polarized radiance
!
! Revision 2.361  2015/12/08 23:23:42  vsnyder
! Put bound (:i_end) on Ref_Corr in call to Get_Tau
!
! Revision 2.360  2015/12/08 19:13:08  vsnyder
! Define pointer association status of Phi_Path_C immediately after Height_Metrics
!
! Revision 2.359  2015/10/28 00:34:13  vsnyder
! Add more magnetic field-related dumps
!
! Revision 2.358  2015/09/22 23:37:26  vsnyder
! Add 3D magnetic field
!
! Revision 2.357  2015/08/25 17:23:05  vsnyder
! Compute PhiWindow in radians for TScat
!
! Revision 2.356  2015/05/28 23:22:44  vsnyder
! Use height above geoid for interpolating the magnetic field to the path.
!
! Revision 2.355  2015/05/01 02:08:36  vsnyder
! Interpolate in height or zeta for magnetic field
!
! Revision 2.354  2015/04/11 01:26:03  vsnyder
! Add more dumps, add (km) to comments about h_path and h_glgrid
!
! Revision 2.353  2015/03/28 02:16:24  vsnyder
! Compute the viewing azimuth to determine whether cross-track viewing is
! taking place.  Use Orbit_Plane_Minor_Axis_sq from Geometry.  Use 3-d
! Value field for ECRtoFOV array.  Use Norm2 to normalize Mag_Path.
!
! Revision 2.352  2014/09/05 20:49:32  vsnyder
! Avoid reference to undefined array elements
!
! Revision 2.351  2014/08/01 01:06:05  vsnyder
! Add code to fill real arrays with sNaN if the private parameter NaN_Fill
! is true.  Set p_path to zero between the tangent points, so it's not
! undefined.  It's not actually used for anything productive, but if the
! compiler has an option to fill with sNaN, and those elements are referenced
! (and the results subsequently not used), a pointless trap occurs.
!
! Revision 2.350  2014/07/18 23:15:44  pwagner
! Aimed for consistency in names passed to allocate_test
!
! Revision 2.349  2014/01/11 01:28:53  vsnyder
! Decruftification
!
! Revision 2.348  2013/08/31 02:30:03  vsnyder
! Improve tracing
!
! Revision 2.347  2013/08/30 03:56:23  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.346  2013/08/17 02:56:19  vsnyder
! Regularize trace usage
!
! Revision 2.345  2013/08/08 02:35:25  vsnyder
! Set Qty_Stuff%derivOK if ExtraJacobian is present
!
! Revision 2.344  2013/08/02 01:24:06  vsnyder
! Add ExtraJacobian to compute derivatives not in state vector
!
! Revision 2.343  2013/07/13 00:06:20  vsnyder
! Move computation of tangent pressures from Compute_Z_PSIG to Tangent_Pressures.
! Remove declarations for unused symbols.
!
! Revision 2.342  2013/06/12 02:35:22  vsnyder
! Make Z_psig allocatable instead of pointer
!
! Revision 2.341  2013/05/22 00:19:10  vsnyder
! Remove unreferenced USE names
!
! Revision 2.340  2013/05/18 00:34:43  vsnyder
! Insert NG fine-grid (GL) points between tangent points, thereby
! regularizing coarse-grid spacing, and reducing significantly the need
! to use c_inds to extract coarse-grid points from the composite grid.
!
! Revision 2.339  2013/04/09 18:23:56  pwagner
! Fixed error in failing to use DerivativeMissingFromStateFun
!
! Revision 2.338  2013/02/04 22:06:27  pwagner
! Added dmiss switch to downgrade severity when qty missing from state
!
! Revision 2.337  2013/01/23 21:24:26  vsnyder
! Use |sin| in TScat phase convolution
!
! Revision 2.336  2012/08/08 20:06:37  vsnyder
! Use CreateVectorValue, DestroyVectorQuantityValue in full cloud model.
! Exchange subscript ordering for vmrArray to improve locality.
!
! Revision 2.335  2012/07/07 00:14:33  vsnyder
! Shorten some comments to avoid gripes about long lines
!
! Revision 2.334  2012/07/06 21:30:41  yanovsky
! Make the module ready for making calls to Convolve_Radiance_Normalization
! and Convolve_Temperature_Deriv_Normalization that compute normalized
! Temperature derivatives.  These calls are currently commented.
!
! Revision 2.333  2012/06/15 23:33:13  vsnyder
! Include sin(theta) factor in TScat phase convolution.  Improve dumps.
!
! Revision 2.332  2012/02/13 23:20:12  pwagner
! DerivativeMissingFromState eases switch from strict to lenient
!
! Revision 2.331  2012/02/10 23:51:33  vsnyder
! Move DeriveFromForwardModelConfig to Retrieval module
!
! Revision 2.330  2011/11/10 23:23:32  vsnyder
! Correct computation of 'combine' argument in frequency averaging
!
! Revision 2.329  2011/11/09 00:29:48  vsnyder
! Monochromatic PFA, more debugging
!
! Revision 2.327  2011/08/12 18:59:57  vsnyder
! Add Do_Calc_Fzp into one call to Comp_Sps_Path_Frq, because it is no longer
! optional, and it is actually needed where it wasn't specified.  Add checks
! that temperature and mixing ratio quantities for which derivatives are
! requested appear in the first state vector.
!
! Revision 2.326  2011/07/29 01:55:18  vsnyder
! Use CloudIce instead of Cloud_A and Cloud_S.  Only one IWC, not IWC_A and
! IWC_S.  More dumps.
!
! Revision 2.325  2011/07/21 20:48:38  honghanh
! Fix an array-index-out-of-bound bug by fixing the declaration fo DBeta_DF_Path_C
! and DBeta_DF_Path_F
!
! Revision 2.324  2011/07/08 18:19:35  yanovsky
! Use get_d2Alpha_df2
!
! Revision 2.323  2011/06/24 23:15:53  pwagner
! Fixed erroneous declaration for D2_DELTA_DF2 that gave non-Hessian runs excess memory footprint
!
! Revision 2.322  2011/06/02 22:31:48  yanovsky
! Add D2_DELTA_DF2 for computations of analytical Hessians in logarithmic basis
!
! Revision 2.321  2011/05/09 17:46:38  pwagner
! Converted to using switchDetail
!
! Revision 2.320  2011/03/31 19:51:03  vsnyder
! Allow 'scattering point not in path' to be a warning
!
! Revision 2.319  2011/03/25 20:46:59  vsnyder
! Delete declarations of unused objects
!
! Revision 2.318  2011/03/23 23:50:42  vsnyder
! Finishing -- hopefully -- TScat derivatives
!
! Revision 2.317  2011/03/23 23:45:32  vsnyder
! This log entry is bogus.  Check in again to get the right one.
! FOV_Convolve_m.f90
!
! Revision 2.316  2011/03/11 03:09:08  vsnyder
! Use Get_dAlpha_df
!
! Revision 2.315  2011/03/04 03:44:40  vsnyder
! Remove declarations for unused stuff
!
! Revision 2.314  2011/02/12 03:57:40  vsnyder
! Add mixing-ratio dependence for H2O derivatives
!
! Revision 2.313  2011/01/29 00:52:32  vsnyder
! Allow PFA without frequency averaging
!
! Revision 2.312  2011/01/28 19:18:44  vsnyder
! Lots of stuff for TScat
!
! Revision 2.311  2010/11/05 23:03:18  pwagner
! Lots of new NaN checks--should make them optional
!
! Revision 2.310  2010/11/03 18:43:51  vsnyder
! Initialize eta_fzp and H_Atmos.  Frequency average LBL second derivatives.
!
! Revision 2.309  2010/09/25 01:08:39  vsnyder
! Cannonball polishing
!
! Revision 2.308  2010/08/28 00:03:12  vsnyder
! Shortened some overly long names.  Corrected allocation for H_Atmos_Frq.
! Moved some TScat stuff to TScat_Support.  Some cannonball polishing.
!
! Revision 2.307  2010/08/27 06:02:39  yanovsky
! Major additions to support for computations of atmospheric second
! derivatives (atmospheric Hessians). FullForwardModel subroutine
! now has dummy variable Hessian.  Also, added
! Frequency_Average_Second_Derivative and
! Frequency_Average_Second_Derivatives subroutines.
!
! Revision 2.306  2010/08/19 02:14:03  vsnyder
! Substantial changes, the most significant being to put the body of the
! frequency loop in a subroutine of its own called One_Frequency, delete
! the subroutine Frequency_Loop, put the frequency loop around the calls
! to One_Frequency, and pass in sections of some arrays that used to be
! accessed by host association, and subscripted with the frequency index.
! Also, some other stuff inching toward TScat processing.
!
! Revision 2.305  2010/06/12 01:30:54  vsnyder
! Give frequency dimension to Alpha_path_c and Beta_path_c.  Make
! Frequency_Loop the body of the frequency loop instead of doing the
! loop in it.  Pass in lower-dimensional sections of several arrays,
! with the section selected according to the frequency loop index.
!
! Revision 2.304  2010/06/07 23:23:53  vsnyder
! Numerous changes inching toward using TScat tables, some of which will
! almost certainly prove to have been dead ends.
!
! Revision 2.303  2010/05/18 18:13:35  honghanh
! Change both_sidebands_setup to get GPH quantity with GetQuantityForForwardModel as the rest of the quantity
!
! Revision 2.302  2010/03/24 20:50:02  vsnyder
! Add more checks to TScat generation code
!
! Revision 2.301  2010/02/05 03:29:19  vsnyder
! Remove unused stuff
!
! Revision 2.300  2010/02/02 01:39:06  vsnyder
! Move USEs closer to where the used stuff is referenced.  Don't compute
! theta if it trying to do so will cause an exception.  Finish applying the
! product rule correctly for TScat derivatives.
!
! Revision 2.299  2010/01/23 01:30:49  vsnyder
! Move USEs into closer proximity to references.  Get cloudIce as vmr if it
! isn't available as a quantityType.  Create TScat stuff with zero size if
! not doing TScat table generation.  Get index of Cloud_A if there is one.
! Provide that Betas for Cloud_A and Cloud_S depend upon IWC mixing ratio.
! Calculate TScat derivatives (still not complete).  Don't compute TScat
! or derivatives where there's no IWC.  Don't use the logIWC quantity type.
! During TScat table computation, interpolate P to pointing angles rather
! than interpolating incident radiance to angles where T is tabluated.
!
! Revision 2.298  2009/12/15 03:30:04  vsnyder
! TScat radiance calculation working again -- P interpreted to xi instead
! of radiance interpolated to theta.
! Exchange order of subscripts of "radiances" to be consistent with k_...
! and to get better locality in frequency loop
!
! Revision 2.297  2009/11/17 23:45:32  vsnyder
! Add R_eq, R_sc arguments to FOV_Convolve_Setup, incomplete TScat stuff
!
! Revision 2.296  2009/09/25 02:45:06  vsnyder
! TScat computation appears to be working
!
! Revision 2.295  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.294  2009/06/16 17:38:17  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.293  2009/06/13 01:15:03  vsnyder
! Intermediate commit of extensive changes for TScat calculation
!
! Revision 2.292  2009/05/14 00:46:02  pwagner
! Gets Deg2Rad from Constants now
!
! Revision 2.291  2009/01/21 01:00:39  pwagner
! Compatible with hastily committed select_nz_list
!
! Revision 2.290  2008/10/03 16:27:14  livesey
! Pushed down LO to support EXTINCTIONV2
!
! Revision 2.289  2008/06/26 00:28:03  vsnyder
! Simplify call to Compute_Z_PSIG.  Simplify some stuff.  Correct a bug:
! If H-PHI iteration doesn't converge and no new zeta is added, z_ig isn't
! defined.  Sneak up a little bit on TScat computation.
!
! Revision 2.288  2008/05/20 00:23:21  vsnyder
! Much rearranging to prepare for TScat calculation.
! Use pointing-interpolated LOS velocity instead of MIF(1) for Doppler
! correction for pointing frequency grid.
! Send vel/C instead of vel to slabs_prep.
!
! Revision 2.287  2008/02/29 01:59:39  vsnyder
! Added a separate H2O continuum routine
!
! Revision 2.286  2007/11/08 02:02:19  vsnyder
! Change name of path_dsdh to dsdh_path for consistency with other names.
! Add do_path_norm switch.  Rearrange printing if switches set to show
! desires to add points but not to add them.
!
! Revision 2.285  2007/07/31 23:49:21  vsnyder
! Add an argument metrics needs for H/Phi failure recovery
!
! Revision 2.284  2007/07/11 22:26:32  vsnyder
! More dumps, change some error handling
!
! Revision 2.283  2007/06/29 19:33:59  vsnyder
! Put the pointing loop body into an internal subroutine
! Make ForwardModelIntermediate_t private to ScanModelModule
!
! Revision 2.282  2007/06/26 01:05:02  vsnyder
! Use column-sparse eta
!
! Revision 2.281  2007/06/08 22:05:33  vsnyder
! Faster d_delta_df = 0, metrics stuff
!
! Revision 2.280  2007/02/01 02:53:47  vsnyder
! Stuff for min zeta and more intersections, plus cannonball polishing
!
! Revision 2.279  2007/01/20 01:08:08  vsnyder
! Decrufting
!
! Revision 2.278  2007/01/19 02:38:53  vsnyder
! Include water in phi refractive correction
!
! Revision 2.277  2007/01/18 00:27:10  vsnyder
! Split Pure_Metrics into Tangent_Metrics and Height_Metrics, insert Earth
! intersection into Zeta grid.
!
! Revision 2.276  2006/12/21 22:59:00  vsnyder
! Get rid of some unused variables
!
! Revision 2.275  2006/12/21 01:34:53  vsnyder
! Finally implemented minimum Zeta
!
! Revision 2.274  2006/12/20 21:22:16  vsnyder
! Split metrics into pure H-Phi calculation, and everything else, in
! preparation for inserting the minimum-Zeta point into the path.
!
! Revision 2.273  2006/12/19 02:53:15  vsnyder
! Change some names, send max coarse path from FullForwardModel to
! FullForwardModelAuto instead of using 2*NLVL, get rid of STATUS
! argument from Metrics, reference H_Path from Earth center instead
! of surface.
!
! Revision 2.272  2006/12/13 02:32:02  vsnyder
! Drag the tangent point around instead of assuming it's the middle one
!
! Revision 2.271  2006/12/08 23:57:08  vsnyder
! Revise earth-intersecting metrics
!
! Revision 2.270  2006/12/04 21:17:28  vsnyder
! Reorganize FullForwardModel to use automatic arrays instead of allocating
! pointer arrays.  Requires testing for zero size instead of testing for
! associated in several subsidiary procedures.
!
! Revision 2.267  2006/08/25 19:42:31  vsnyder
! Recommited to get correct comment into the log: Compute earthradc_sq for
! metrics, since we both need it, more accurate tracing, cannonball polishing.
!
! Revision 2.266  2006/08/25 19:37:07  vsnyder
! Earthradc_sq, etc.
!
! Revision 2.265  2006/07/06 23:16:19  vsnyder
! Need to allocate dh_dt_glgrid even if no temperature derivatives
!
! Revision 2.264  2006/06/29 19:34:43  vsnyder
! Changes due to metrics handling only one tangent
!
! Revision 2.263  2006/06/19 19:26:56  vsnyder
! OOPS, out of bounds subscript possible in dsdh_path
!
! Revision 2.262  2006/06/16 20:32:30  vsnyder
! Define NGP1 in glnp
!
! Revision 2.261  2006/06/16 00:49:10  vsnyder
! Improved non-GL correction, add TeXnicalities
!
! Revision 2.260  2006/04/25 23:25:36  vsnyder
! Revise DACS filter shape data structure
!
! Revision 2.259  2006/04/21 22:25:20  vsnyder
! Cannonball polishing
!
! Revision 2.258  2006/04/11 18:34:37  vsnyder
! Add tanh(h nu / k T) to Get_D_Deltau_Pol.  Use DACS frequencies in
! Frequency_Setup_1.
!
! Revision 2.257  2006/03/06 20:44:06  vsnyder
! Avoid appearance of out-of-bounds subscript during frequency averaging
!
! Revision 2.256  2006/02/08 21:38:18  vsnyder
! Add relative humidity (RHi) calculation
!
! Revision 2.255  2006/02/08 01:02:01  vsnyder
! More stuff for spectroscopy derivatives
!
! Revision 2.254  2006/01/05 00:03:52  vsnyder
! Implement refractive correction for Phi
!
! Revision 2.253  2005/12/10 01:51:54  vsnyder
! Cannonball polishing
!
! Revision 2.252  2005/12/07 19:43:32  vsnyder
! Mistakenly committed needing Phi_Refractive_Correction
!
! Revision 2.251  2005/12/07 01:30:04  vsnyder
! More on getting correct size for RadV
!
! Revision 2.250  2005/12/07 00:35:04  vsnyder
! Allocate RadV with correct size for PFA and no LBL
!
! Revision 2.249  2005/11/21 22:57:41  vsnyder
! PFA derivatives stuff, plug some memory leaks
!
! Revision 2.248  2005/11/05 03:38:13  vsnyder
! Frequency_Average_Derivative doesn't need Tau, cannonball polishing
!
! Revision 2.247  2005/11/03 03:57:45  vsnyder
! Don't try to look at filter shapes for DACS
!
! Revision 2.246  2005/11/02 21:38:48  vsnyder
! Repair a broken tracing message
!
! Revision 2.245  2005/11/01 23:02:08  vsnyder
! PFA Derivatives, use precomputed ShapeInds
!
! Revision 2.244  2005/09/17 01:02:38  vsnyder
! Out of bounds subscript
!
! Revision 2.243  2005/09/17 00:49:53  vsnyder
! Revise arrays for spectroscopic derivatives, plus some cannonball polishing
!
! Revision 2.242  2005/09/03 01:21:33  vsnyder
! Spectral parameter offsets stuff
!
! Revision 2.241  2005/08/03 18:03:42  vsnyder
! Scan averaging, some spectroscopy derivative stuff
!
! Revision 2.240  2005/07/08 19:40:51  vsnyder
! OOPS, forgot to nullify Rad_FFT
!
! Revision 2.239  2005/07/08 00:12:11  vsnyder
! Get Rad_FFT from Convolve_Radiance to Convolve_Temperature_Deriv
!
! Revision 2.238  2005/07/06 02:17:20  vsnyder
! Revisions for spectral parameter derivatives
!
! Revision 2.237  2005/06/09 02:34:16  vsnyder
! Move stuff from l2pc_pfa_structures to slabs_sw_m
!
! Revision 2.236  2005/06/03 01:59:25  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades
!
! Revision 2.235  2005/04/19 20:16:45  livesey
! Added a couple of nullifys
!
! Revision 2.234  2005/03/28 20:26:25  vsnyder
! Lots of PFA stuff
!
! Revision 2.233  2005/03/10 00:28:09  pwagner
! Better patching of ptg_angles avoids list-out-of-order in Hunt
!
! Revision 2.232  2005/02/17 02:35:29  vsnyder
! Do PFA on fine path if necessary
!
! Revision 2.231  2005/02/16 23:16:49  vsnyder
! Revise data structures for split-sideband PFA
!
! Revision 2.230  2004/12/28 00:28:02  vsnyder
! Remove unreferenced use name
!
! Revision 2.229  2004/12/13 20:37:23  vsnyder
! Moved stuff from get_species_data to ForwardModelConfig, some cannonball polishing
!
! Revision 2.228  2004/11/05 19:38:39  vsnyder
! Got rid of DerivedFromForwardModel component of config
!
! Revision 2.227  2004/11/01 20:26:35  vsnyder
! Reorganization of representation for molecules and Beta groups; PFA may be broken for now
!
! Revision 2.226  2004/10/13 01:08:27  vsnyder
! Moved some checking to ForwardModelSupport
!
! Revision 2.225  2004/10/07 23:26:09  vsnyder
! Changes in Beta_Group structure
!
! Revision 2.224  2004/10/06 21:27:41  vsnyder
! More work on PFA -- seems to work now
!
! Revision 2.223  2004/09/14 17:20:45  bill
! add mif dependent los vel correction-wgr
!
! Revision 2.222  2004/09/04 01:50:30  vsnyder
! get_Beta_path_m.f90
!
! Revision 2.221  2004/09/01 01:48:27  vsnyder
! Status flags, more work on PFA
!
! Revision 2.220  2004/08/06 22:40:17  livesey
! Better patch for ptg_angles
!
! Revision 2.219  2004/08/06 01:40:39  livesey
! Upgraded the precision of the ptg dump.
!
! Revision 2.218  2004/08/06 01:24:55  livesey
! Typo!
!
! Revision 2.217  2004/08/06 01:24:22  livesey
! Minor bug fix in ptg_angles dump
!
! Revision 2.216  2004/08/05 20:53:50  vsnyder
! More PFA preparations, some cannonball polishing
!
! Revision 2.215  2004/08/05 20:24:00  livesey
! Added ptg switch to dump ptg_angles
!
! Revision 2.214  2004/08/03 22:07:10  vsnyder
! Inching further toward PFA
!
! Revision 2.213  2004/07/30 19:53:18  livesey
! Bug fix, forbid extrapolation in Estimate_tan_phi
!
! Revision 2.212  2004/07/08 21:00:23  vsnyder
! Inching toward PFA
!
! Revision 2.211  2004/06/10 00:59:56  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.210  2004/05/27 23:23:50  pwagner
! named parameter clean= to dump procedures
!
! Revision 2.209  2004/05/17 22:05:11  livesey
! A change to refraction and to k_atmos handling to avoid explosions.
!
! Revision 2.208  2004/04/24 02:27:05  vsnyder
! Cosmetic changes
!
! Revision 2.207  2004/04/19 21:01:37  vsnyder
! Put size of gl_slabs in call to get_gl_slabs_arrays
!
! Revision 2.206  2004/04/17 00:37:00  vsnyder
! Analytic temperature derivatives
!
! Revision 2.205  2004/04/05 21:08:42  jonathan
! delet temp_prof
!
! Revision 2.204  2004/04/02 00:59:24  vsnyder
! Get catalog from slabs structure
!
! Revision 2.203  2004/03/31 20:35:26  jonathan
! bug fix in handling clouds
!
! Revision 2.202  2004/03/30 02:26:17  livesey
! Bug fix in Jonathan's w0 handling
!
! Revision 2.201  2004/03/30 02:00:36  vsnyder
! Remove USE for unreferenced symbol.  Don't try to fill tpath_m and
! tpath_p if they're not allocated.
!
! Revision 2.200  2004/03/27 03:35:27  vsnyder
! Add pointer to catalog in slabs_struct.  Use it so as not to need to drag
! line centers and line widths around.  Write slabs_lines and slabswint_lines
! to get sum of Beta over all lines; put slabs_struct instead of its components
! in the calling sequence.
!
! Revision 2.199  2004/03/20 04:05:50  vsnyder
! Moved SpeedOfLight from units to physics
!
! Revision 2.198  2004/03/20 01:15:16  jonathan
! add in scattering correction term in two t_script
!
! Revision 2.196  2004/03/01 19:22:14  jonathan
! following the changes made to load_one_item
!
! Revision 2.195  2004/02/14 00:23:48  vsnyder
! New DACS convolution algorithm
!
! Revision 2.194  2004/02/05 23:30:01  livesey
! Finally implemented code to do correct handing of sideband fraction in
! single sideband radiometers.
!
! Revision 2.193  2004/02/03 02:48:35  vsnyder
! Progress (hopefully) on polarized temperature derivatives.
! Implement DACs frequency convolution.
!
! Revision 2.192  2004/01/23 19:13:42  jonathan
! add an extra-flag for tscat in load_one_item
!
! Revision 2.191  2003/12/08 21:38:33  jonathan
! some minor changes
!
! Revision 2.190  2003/12/08 17:52:02  jonathan
! update for 2d cldfwm
!
! Revision 2.189  2003/12/07 19:45:46  jonathan
! update for 2D cloud FWM
!
! Revision 2.188  2003/12/01 17:24:05  jonathan
! add scat_alb
!
! Revision 2.187  2003/11/24 22:10:14  vsnyder
! Multiply Beta_path_polarized_f by tanh1_f if derivatives needed
!
! Revision 2.186  2003/11/20 17:33:50  pwagner
! Nullify some things otherwise left unassociated
!
! Revision 2.185  2003/11/19 22:21:34  jonathan
! interpolate scat_src to tscat_path
!
! Revision 2.184  2003/11/17 18:04:15  jonathan
! scat_src output from T_scat is correct
!
! Revision 2.183  2003/11/12 00:10:55  jonathan
! some changes due to cloud construction
!
! Revision 2.182  2003/11/07 03:18:49  vsnyder
! Cosmetic changes
!
! Revision 2.181  2003/11/06 01:13:14  bill
! fixed DACS freq extrapolation problem
!
! Revision 2.180  2003/11/05 19:26:01  jonathan
! add stuff for use later in cld model
!
! Revision 2.179  2003/11/04 02:49:13  vsnyder
! Calculate coarse-path indices where GL is needed
!
! Revision 2.178  2003/11/04 01:57:15  vsnyder
! Move trapezoid correction to a better place, cosmetics
!
! Revision 2.177  2003/11/03 23:15:16  vsnyder
! Get rid of path_ds_dh procedure -- a one-liner used in one place
!
! Revision 2.176  2003/11/01 03:02:57  vsnyder
! Compute del_zeta, ds_dz_gw and dh_dz_gw for [d]rad_tran*; change
! indices_c to c_inds for consistency with usage in rad_tran_m.
!
! Revision 2.175  2003/10/30 20:37:00  vsnyder
! Compute del_zeta here for *rad_tran_*
!
! Revision 2.174  2003/10/29 00:43:51  livesey
! Added support for the forceFoldedOutput option
!
! Revision 2.173  2003/10/28 22:05:53  jonathan
! add l_gph for use in cloud model
!
! Revision 2.172  2003/10/09 19:30:18  vsnyder
! Cosmetic changes
!
! Revision 2.171  2003/09/24 02:53:48  vsnyder
! Cosmetic changes
!
! Revision 2.170  2003/09/09 00:04:27  vsnyder
! Supply E and Sqrt_Earth_Rflty to mcrt_der
!
! Revision 2.169  2003/08/16 01:14:58  vsnyder
! Use radiometers%polarization to choose 1,1 or 2,2 element
!
! Revision 2.168  2003/08/15 20:29:26  vsnyder
! Implement polarized VMR derivatives
!
! Revision 2.167  2003/08/15 18:50:21  vsnyder
! Preparing the way for polarized vmr derivatives
!
! Revision 2.166  2003/08/12 23:07:32  vsnyder
! Futzing with comments
!
! Revision 2.165  2003/08/12 21:58:37  vsnyder
! Use trapezoid instead of rectangle to integrate non-GL panels
!
! Revision 2.164  2003/08/12 18:22:10  michael
! Contribution of scalar molecules to polarized absorption is now added to
! coarse grid Alpha_path_polarized (1/4 1/2 1/4) instead of to diagonal of
! incoptdepth_pol.  This make Alpha_path_polarized correct for later use
! in gl corrections.
!
! Revision 2.163  2003/07/15 23:07:21  vsnyder
! Simplify Freq_Avg
!
! Revision 2.162  2003/07/15 22:10:09  livesey
! Added support for hybridModel
!
! Revision 2.161  2003/07/15 18:16:48  livesey
! Catalog now split by sideband, also changed no_ele to max_ele in
! allocates
!
! Revision 2.160  2003/07/09 23:40:13  vsnyder
! Use new AllocateSlabs routine
!
! Revision 2.159  2003/07/09 22:27:42  vsnyder
! More futzing
!
! Revision 2.158  2003/07/09 20:23:50  vsnyder
! Futzing
!
! Revision 2.157  2003/07/09 20:14:19  livesey
! Anticipative bug fix commented out.
!
! Revision 2.156  2003/07/04 03:40:13  vsnyder
! Simplify dump in case exp(incoptdepth_pol) fails
!
! Revision 2.155  2003/07/04 02:50:15  vsnyder
! Simplify interface to Get_GL_Slabs_Arrays, correct a blunder introduced around July 17
!
! Revision 2.154  2003/06/27 23:43:34  vsnyder
! Remove unreferenced USE names
!
! Revision 2.153  2003/06/27 22:09:48  vsnyder
! Check status from rad_tran_pol and report an error if overflow occurred
!
! Revision 2.152  2003/06/27 00:59:53  vsnyder
! Simplify interface to Get_Species_Data
!
! Revision 2.151  2003/06/25 02:41:37  vsnyder
! Futzing
!
! Revision 2.150  2003/06/18 22:26:41  vsnyder
! Restored the check to the wrong place at 2.149
!
! Revision 2.149  2003/06/18 19:29:30  vsnyder
! Replace a check inadvertently deleted in rev 2.146
!
! Revision 2.148  2003/06/18 17:23:16  bill
! fixed NAG associated bug
!
! Revision 2.147  2003/06/18 14:54:04  bill
! added subsetting feature for T-ders
!
! Revision 2.146  2003/06/18 01:59:20  vsnyder
! Move checking that all signals in config are for same radiometer, module
! and sideband, and computation of sidebandStart and sidebandStop, to
! ForwardModelSupport.  Remove vector quantity validation, as that's now
! done in Construct.  Futzing.
!
! Revision 2.145  2003/06/13 00:00:27  vsnyder
! Move multiplication of Beta_path by tanh into FullForwardModel
!
! Revision 2.144  2003/06/10 15:06:54  bill
! fixed polarized t-derivs
!
! Revision 2.143  2003/06/09 20:52:37  vsnyder
! More work on polarized derivatives
!
! Revision 2.142  2003/06/09 17:38:47  livesey
! Use GetQuantityForForwardModel in more places
!
! Revision 2.141  2003/05/29 16:37:38  livesey
! Renamed sideband fraction
!
! Revision 2.140  2003/05/26 01:42:50  michael
! Two temporary fixes only relevant to the polarized model.
! Added a bug-fix removing scalar contribution to magnetic GL.
! I don't know why this works but it makes agreement with scalar model
! much better.  I also added an undocumented switch -Scrosspol to allow
! the use of the (2,2) element of the radiance tensor rather than the (1,1).
! This code should be removed when we set up the l2cf to do this properly.
!
! Revision 2.139  2003/05/22 20:01:17  vsnyder
! Cosmetic changes
!
! Revision 2.138  2003/05/22 04:03:41  livesey
! Now handles elevation offset as a channel by channel quantity
!
! Revision 2.137  2003/05/20 00:05:39  vsnyder
! Move some stuff to subroutines
!
! Revision 2.136  2003/05/15 03:29:44  vsnyder
! Implement polarized model's temperature derivatives
!
! Revision 2.135  2003/05/09 19:26:36  vsnyder
! Expect T+DT instead of T and DT separately in Get_GL_Slabs_Arrays,
! initial stuff for temperature derivatives of polarized radiance.
!
! Revision 2.134  2003/05/05 23:00:24  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.126.2.40  2003/04/24 21:57:05  vsnyder
! Check FwdModelConf%incl_cld instead of associated(cloudIce)
!
! Revision 2.126.2.39  2003/04/24 19:01:19  vsnyder
! Compute MIF correctly
!
! Revision 2.126.2.38  2003/04/16 20:44:20  vsnyder
! Move allocate for scat_src%values out of the loops
!
! Revision 2.126.2.37  2003/04/16 19:32:02  vsnyder
! Move working storage for T_Scat into T_Scat
!
! Revision 2.126.2.36  2003/04/15 23:45:45  vsnyder
! Correct choice of mag field coordinate indices, futzing
!
! Revision 2.126.2.35  2003/04/15 17:55:19  jonathan
! FullForwardModel_m.f90 change scat_src as non-pointer
!
! Revision 2.126.2.34  2003/04/15 14:18:22  jonathan
! interpolate scat_src%values to LOS
!
! Revision 2.126.2.33  2003/04/10 16:16:06  jonathan
! fix bug
!
! Revision 2.126.2.32  2003/04/08 21:44:23  jonathan
! remove if condition for Modify_values_for_supersat
!
! Revision 2.126.2.31  2003/04/07 17:13:52  jonathan
! modified cal to T_scat
!
! Revision 2.126.2.30  2003/04/04 20:43:53  jonathan
! update call to T_scat
!
! Revision 2.126.2.29  2003/03/28 18:14:19  jonathan
! add new scaterring source module
!
! Revision 2.126.2.28  2003/03/27 23:32:19  vsnyder
! Houst a reshape out of a loop
!
! Revision 2.126.2.27  2003/03/27 00:50:35  vsnyder
! Get ECRtoFOV, use it to rotate mag_path
!
! Revision 2.126.2.26  2003/03/22 04:19:13  vsnyder
! Cosmetic changes
!
! Revision 2.126.2.25  2003/03/22 04:03:04  vsnyder
! Move Beta_Group_T and Dump_Beta_Group from get_Beta_path to Get_Species_Data
!
! Revision 2.126.2.24  2003/03/22 02:48:53  vsnyder
! Interpolate the magnetic field onto the path
!
! Revision 2.126.2.23  2003/03/21 02:49:55  vsnyder
! Use an array of pointers to quantities instead of searching several times
! using GetQuantityForForwardModel.  Move stuff around.  Cosmetic changes.
!
! Revision 2.126.2.22  2003/03/20 19:21:05  vsnyder
! More futzing with grids_t and stuff that uses it
!
! Revision 2.126.2.21  2003/03/20 01:42:25  vsnyder
! Revise Grids_T structure
!
! Revision 2.126.2.19  2003/03/07 23:17:05  vsnyder
! Use ASSOCIATED instead of PRESENT for "optional" arguments GET_CHI_OUT
!
! Revision 2.126.2.18  2003/03/06 22:47:29  vsnyder
! Polarized radiative transfer sort-of works
!
! Revision 2.126.2.17  2003/03/05 03:40:07  vsnyder
! Do tanh correction for polarized, more polarized work, simplifications
!
! Revision 2.126.2.16  2003/03/04 20:24:44  dwu
! add a temporay fix for the tangent height crossover problem
!
! Revision 2.126.2.15  2003/03/01 03:19:28  vsnyder
! More work on polarized radiative transfer.  Still doesn't work, but at
! least it doesn't break the nonpolarized case.
!
! Revision 2.126.2.14  2003/02/28 00:15:59  vsnyder
! Print 3 digits in channel number in 'rad' diagnostic output
!
! Revision 2.126.2.13  2003/02/27 23:35:01  vsnyder
! Revise polarized processing
!
! Revision 2.126.2.12  2003/02/25 00:53:09  jonathan
! add grids_tscat
!
! Revision 2.126.2.11  2003/02/24 23:40:31  jonathan
! change input for get_Beta_path_cloud
!
! Revision 2.126.2.10  2003/02/24 23:26:43  jonathan
! change temp_supersat to temp_prof
!
! Revision 2.126.2.9  2003/02/22 00:49:26  vsnyder
! Delete moleculesPol, moleculeDerivativesPol, add Polarized to ForwardModelConfig
!
! Revision 2.126.2.8  2003/02/21 21:34:32  michael
! Made index for dumping dacs radiances be 1:129 for dacs channels 0:128
!
! Revision 2.126.2.7  2003/02/15 00:29:11  vsnyder
! Don't exp(incoptdepth) in mcrt -- it's done here
!
! Revision 2.126.2.6  2003/02/14 23:31:14  vsnyder
! Delete unreferenced names
!
! Revision 2.126.2.5  2003/02/14 23:27:31  vsnyder
! Move stuff to Get_Species_Data
!
! Revision 2.126.2.4  2003/02/14 03:53:41  vsnyder
! Initial commit of polarized model
!
! Revision 2.126.2.3  2003/02/14 00:21:42  jonathan
! add singl. scat. albedo W0, ph funct PHH
!
! Revision 2.126.2.2  2003/02/13 22:26:18  jonathan
! changes dimension for Beta_path_cloud also delocate it
!
! Revision 2.126.2.1  2003/02/13 17:35:01  bill
! fixes gl_ind bug and interfaces to get_Beta
!
! Revision 2.126  2003/02/11 00:48:18  jonathan
! changes made after adding get_Beta_path_cloud
!
! Revision 2.125  2003/02/08 01:03:00  livesey
! Bug fix in call to rad_tran
!
! Revision 2.124  2003/02/07 01:07:41  jonathan
! add in option to compute dry and super-saturation case in load_sps
!
! Revision 2.123  2003/02/06 22:04:25  vsnyder
! Add f_moleculesPol, f_moleculeDerivativesPol, delete f_polarized
!
! Revision 2.122  2003/02/06 19:15:47  jonathan
!  fix a bug
!
! Revision 2.121  2003/02/06 19:06:49  jonathan
! add eta_iwc, eta_iwc_zp, do_calc_iwc, do_cala_iwc_zp and also not passing
! through comp_eta_docalc and comp_sps_path_frq if fwdModelConf%Incl_Cld is
! false
!
! Revision 2.120  2003/02/06 05:55:47  livesey
! Fix to sort of fix Jonathan's cloud ice stuff.
!
! Revision 2.119  2003/02/06 00:20:08  jonathan
! Add in many stuff to deal with clouds CloudIce, iwc_path, etc
!
! Revision 2.118  2003/02/03 23:18:43  vsnyder
! Squash a bug in deallocating Beta_path_polarized
!
! Revision 2.117  2003/02/03 22:58:17  vsnyder
! Plug a memory leak, delete gl_ndx, some polarized stuff
!
! Revision 2.116  2003/02/03 19:00:36  bill
! changed interface to rad tran to speed up program
!
! Revision 2.115  2003/02/01 02:33:22  vsnyder
! Plug a bunch of memory leaks
!
! Revision 2.114  2003/01/31 17:53:39  jonathan
! change z_path to z_path_c in passing to get_Beta_path
!
! Revision 2.113  2003/01/31 17:15:49  jonathan
! add Inc_Cld to get_Beta_path
!
! Revision 2.112  2003/01/31 01:53:01  vsnyder
! Move array temps to arrays explicitly allocated outside the loop
!
! Revision 2.111  2003/01/30 00:16:35  jonathan
! add z_path to get_Beta_path
!
! Revision 2.110  2003/01/26 04:42:42  livesey
! Added profiles/angle options for phiWindow
!
! Revision 2.109  2003/01/21 18:20:32  vsnyder
! Put dimensions back onto actual arguments to path_contrib
!
! Revision 2.108  2003/01/18 03:36:09  vsnyder
! Undo ill-advised cosmetic changes -- that weren't cosmetic
!
! Revision 2.106  2003/01/16 23:13:03  livesey
! Added MaxRefraction stuff
!
! Revision 2.105  2003/01/16 18:04:01  jonathan
! add Do_1D option to get_gl_slabs_arrays
!
! Revision 2.104  2003/01/14 21:48:58  jonathan
! add i_saturation
!
! Revision 2.103  2003/01/10 21:55:26  vsnyder
! Move SpeedOfLight from Geometry ot Units
!
! Revision 2.102  2003/01/08 00:16:39  vsnyder
! Use "associated" instead of "present" to control optional computations.
! Cosmetic changes, too.
!
! Revision 2.101  2002/12/12 01:12:47  vsnyder
! Let InvalidQuantity have a length > 1
!
! Revision 2.100  2002/11/15 01:33:08  livesey
! Added allLinesForRadiometer functionality.
!
! Revision 2.99  2002/11/13 17:07:44  livesey
! Passes FwdModelExtra into convolve/no_conv
!
! Revision 2.98  2002/10/26 00:13:35  livesey
! Made the warning about lines less common as it checks for continuum
! aswell.
!
! Revision 2.97  2002/10/25 01:12:43  livesey
! Added an array to accumulate the stuff such as one_tan_ht etc.
! Also put in but commented out some code to dump it.
!
! Revision 2.96  2002/10/10 19:38:22  vsnyder
! Mostly cosmetic, some performance improvements
!
! Revision 2.95  2002/10/10 01:46:50  livesey
! Whoops, typo fix
!
! Revision 2.94  2002/10/10 01:28:06  livesey
! Bug fix in k_temp windowing
!
! Revision 2.93  2002/10/08 17:08:03  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.92  2002/10/04 23:49:50  vsnyder
! More cosmetic changes
!
! Revision 2.91  2002/10/02 23:06:42  vsnyder
! Add 'seez' switch, instead of uncommenting include, to get Zvi's debug print.
! Get SpeedOfLight from Geometry.  Cosmetic changes.
!
! Revision 2.90  2002/09/26 18:02:49  livesey
! Now uses GetQuantityForForwardModel.
!
! Revision 2.89  2002/09/10 17:05:38  livesey
! New update arguments to convolve/noconvole
!
! Revision 2.88  2002/09/07 02:18:50  vsnyder
! More cosmetic changes
!
! Revision 2.87  2002/09/06 20:23:22  livesey
! Merged in bug fixes with change from Van
!
! Revision 2.86  2002/09/06 18:19:04  vsnyder
! Cosmetic changes
!
! Revision 2.85  2002/09/05 21:00:38  vsnyder
! Get rid of some auxiliary variables
!
! Revision 2.84  2002/08/26 20:02:22  livesey
! Checks whether the jacobian is present when setting internal
! temp/atmos_der
!
! Revision 2.83  2002/08/22 23:13:03  livesey
! New intermediate frequency based frq_bases
!
! Revision 2.82  2002/08/20 23:33:23  livesey
! Fixed bug with handling of extinction
!
! Revision 2.81  2002/08/20 22:37:04  livesey
! Moved uses inside routine
!
! Revision 2.80  2002/08/02 00:12:07  bill
! just testing
!
! Revision 2.79  2002/07/31 23:27:54  bill
! added feature to user supply grid points to the psig and tangent grids
!
! Revision 2.78  2002/07/31 00:03:45  livesey
! Embarassing bug fix, was seeking wrong vector quantity.
!
! Revision 2.77  2002/07/30 20:03:59  livesey
! More sideband fixes
!
! Revision 2.76  2002/07/30 20:03:19  livesey
! Fixed bug which had it using the wrong sideband ratio
!
! Revision 2.75  2002/07/29 21:41:30  bill
! no changes, just debugging
!
! Revision 2.74  2002/07/23 22:26:31  livesey
! Added ptan_der, set k_atmos to zero on creation
!
! Revision 2.73  2002/07/19 23:35:49  bill
! fixed undefined surf angle
!
! Revision 2.72  2002/07/11 20:51:20  bill
! fixed bug regarding req
!
! Revision 2.71  2002/07/09 17:37:10  livesey
! Fixed more DACs problems
!
! Revision 2.70  2002/07/08 17:45:37  zvi
! Make sure spect_der is turned off for now
!
! Revision 2.69  2002/07/05 07:52:45  zvi
! Coor. switch (phi,z) -> (z,phi)
!
! Revision 2.68  2002/06/28 21:41:36  livesey
! Repeat of bug fix with atmos_der being deallocated in wrong place.
!
! Revision 2.67  2002/06/28 11:06:46  zvi
! Now computing dI/dPtan using chain rule
!
! Revision 2.66  2002/06/26 19:58:48  livesey
! Bug fix with DAC channel numbering
!
! Revision 2.65  2002/06/24 21:11:24  zvi
! Adding Grids_tmp stracture and modifying calling sequences
!
! Revision 2.61  2002/06/17 17:12:15  bill
! fixed yet another bug--wgr
!
! Revision 2.60  2002/06/17 16:30:52  bill
! inc zvis changes--wgr
!
! Revision 2.58  2002/06/13 22:40:38  bill
! some variable name changes--wgr
!
! Revision 2.57  2002/06/11 22:20:10  bill
! work in progress--wgr
!
! Revision 2.56  2002/06/07 23:22:58  bill
! add debug switch--wgr
!
! Revision 2.55  2002/06/07 23:01:25  bill
! work in progress
!
! Revision 2.54  2002/06/07 04:50:03  bill
! fixes and improvements--wgr
!
! Revision 2.53  2002/06/05 17:20:28  livesey
! Fixed tan_temp
!
! Revision 2.52  2002/06/04 23:06:47  livesey
! On the way to having phiTan
!
! Revision 2.51  2002/06/04 10:27:59  zvi
! Encorporate deriv. flag into convolution
!
! Revision 2.50  2002/05/28 17:09:14  livesey
! Removed print statement
!
! Revision 2.49  2002/05/24 17:10:57  livesey
! Fixed bug with my_catalog(?)%lines not being associated for parent
! species.
!
! Revision 2.48  2002/05/23 21:01:11  livesey
! No, that was the wrong thing to do.  Think a bit more.
!
! Revision 2.47  2002/05/23 20:55:20  livesey
! Put more checking around case where a molecule has no lines.
!
! Revision 2.46  2002/05/22 19:44:51  zvi
! Fix a bug in the mol. index loop
!
! Revision 2.45  2002/05/17 22:13:20  livesey
! Bug fix for case where channels start at zero.
!
! Revision 2.44  2002/05/14 22:40:45  livesey
! Bug fix in change in line gathering.  Never got to run it mercifully!
!
! Revision 2.43  2002/05/14 22:32:45  livesey
! Added single sideband stuff.  Also skip line gathering for parent
! molecules.
!
! Revision 2.42  2002/05/14 00:19:10  livesey
! Minor bug fixes
!
! Revision 2.41  2002/05/10 16:18:45  livesey
! Code for dealing with new channel shape information
!
! Revision 2.40  2002/05/08 08:53:42  zvi
! All radiometers grid concept implemented
!
! Revision 2.39  2002/05/03 23:29:18  livesey
! Added direction and split sideband ratio stuff
!
! Revision 2.38  2002/02/22 00:52:06  bill
! fixed units error for light speed--wgr
!
! Revision 2.37  2002/02/20 22:19:45  zvi
! Reversing the subset logic ..
!
! Revision 2.36  2002/02/16 10:32:16  zvi
! Fixing small bug..
!
! Revision 2.35  2002/02/16 06:49:59  zvi
! Retain deriv flag code ..
!
! Revision 2.32  2002/02/14 19:05:01  bill
! Fixed no spectral avg bug--wgr
!
! Revision 2.31  2002/02/13 20:35:47  livesey
! Added some nullifies
!
! Revision 2.30  2002/02/08 00:46:05  zvi
! Some cosmetic changes..
!
! Revision 2.29  2002/02/07 00:36:31  zvi
! Fix a bug - phi_tan non defined..
!
! Revision 2.28  2002/02/05 21:54:29  zvi
! Fix a bug ..
!
! Revision 2.27  2002/02/04 22:44:40  zvi
! Fixing some bugs in the automatic grid selection process
!
! Revision 2.26  2002/02/02 11:20:17  zvi
! Code to overwrite the l2cf integration & tanget grids
!
! Revision 2.25  2002/01/30 01:11:18  zvi
! Fix bug in user selectable coeff. code
!
! Revision 2.24  2002/01/27 08:37:45  zvi
! Adding Users selected coefficients for derivatives
!
! Revision 2.22  2002/01/08 01:01:19  livesey
! Some changes to my_catalog and one_tan_height and one_tan_temp
!
! Revision 2.21  2001/12/26 04:05:04  zvi
! Convert phi_tan to Radians
!
! Revision 2.20  2001/12/14 23:43:05  zvi
! Modification for Grouping concept
!
! Revision 2.19  2001/11/25 07:57:12  zvi
! Fixing inconsistency in k_xxx instance loops
!
! Revision 2.18  2001/11/20 10:23:27  zvi
! Fixing window bug & diemsion
!
! Revision 2.17  2001/11/20 01:18:59  zvi
! Fixing Shifting Window bug
!
! Revision 2.16  2001/11/15 01:21:57  zvi
! Extiction debug fix
!
! Revision 2.15  2001/11/10 00:46:40  zvi
! Adding the EXTINCTION capabilitis
!
! Revision 2.14  2001/11/08 21:52:23  jonathan
! add spec_tags to Molecules
!
! Revision 2.13  2001/11/08 09:56:59  zvi
! Fixing a bug..
!
! Revision 2.12  2001/11/08 00:11:29  livesey
! Added extinction stuff
!
! Revision 2.11  2001/11/07 21:19:01  livesey
! Put Zvi's change comments back
!
! Revision 2.10  2001/11/07 21:16:56  livesey
! Now defaults to *not* using a line if no bands listed.
!
! Revision 2.9  2001/11/07 09:58:41  zvi
! More effective code for sps_path calculations
!
! Revision 2.8  2001/11/03 01:33:35  livesey
! Add more informative message if no spectroscopy information available
! for a molecule
!
! Revision 2.7  2001/11/02 10:47:57  zvi
! Implementing frequecy grid
!
! Revision 2.6  2001/10/12 20:40:25  livesey
! Moved sideband ratio check
!
! Revision 2.5  2001/10/09 22:39:08  livesey
! Allow for molecules with zero lines.  This may need attention
! from Bill/Zvi later on.
!
! Revision 2.4  2001/10/02 16:51:41  livesey
! Removed fmStat%finished and reordered loops in forward models
!
! Revision 2.3  2001/09/19 04:38:48  livesey
! Lines per band stuff works now
!
! Revision 2.2  2001/09/18 02:04:38  livesey
! Bug fix with signals/spectroscopy interaction
!
! Revision 2.1  2001/09/18 01:23:19  livesey
! Added band discrimination for lines catalog.  Not tested yet.
!
! Revision 2.0  2001/09/17 20:26:25  livesey
! New forward model
!
! Revision 1.5.2.56  2001/09/14 22:19:39  livesey
! Fixed bug with sv_i in frequency averaging of k_atmos_frq
!
! Revision 1.5.2.55  2001/09/13 19:57:43  livesey
! Fixed a small memory leak
!
! Revision 1.5.2.54  2001/09/13 19:36:27  livesey
! Added some more useful diagnotics/trace statements
!
! Revision 1.5.2.53  2001/09/13 11:18:15  zvi
! Fix temp. derv. bug
!
! Revision 1.5.2.52  2001/09/13 02:03:08  livesey
! Left a print statement in by mistake
!
! Revision 1.5.2.51  2001/09/13 01:48:39  livesey
! Added lots of deallocates, slight problem with temperature derivatives
! and convolution.
!
! Revision 1.5.2.50  2001/09/13 00:42:18  livesey
! Fixed memory leak

! Revision 1.5.2.49  2001/09/12 23:44:18  livesey
! Fixed bug with calling sequence for metrics

! Revision 1.5.2.48  2001/09/12 22:56:05  livesey
! Got rid of print statement.

! Revision 1.5.2.47  2001/09/12 22:46:26  livesey
! Put dimension limit back on dh_dt_path

! Revision 1.5.2.46  2001/09/12 22:38:25  livesey
! Bug fix

! Revision 1.5.2.45  2001/09/12 04:42:32  zvi
! Fixing Conv. bug

! Revision 1.5.2.44  2001/09/12 01:00:34  livesey
! Fixed problem with vmr derivatives

! Revision 1.5.2.43  2001/09/12 00:32:36  livesey
! Fixed printing

! Revision 1.5.2.42  2001/09/12 00:30:00  zvi
! Put printing stmt. back

! Revision 1.5.2.41  2001/09/12 00:12:03  livesey
! Single channel no convolution, derivatives or frequency averaging works

! Revision 1.5.2.40  2001/09/11 21:24:27  livesey
! Interim version

! Revision 1.5.2.39  2001/09/11 20:48:04  livesey
! Added dump and stop statement

! Revision 1.5.2.38  2001/09/11 20:48:22  zvi
! Develop.

! Revision 1.5.2.37  2001/09/11 08:13:33  zvi
! Fixing bugs, adding printing code

! Revision 1.5.2.36  2001/09/11 01:36:54  livesey
! More tidy ups

! Revision 1.5.2.35  2001/09/11 01:27:14  livesey
! It compiles

! Revision 1.5.2.34  2001/09/11 00:50:32  zvi
! Convolution code..done

! Revision 1.5.2.33  2001/09/11 00:01:06  zvi
! adding convolution..incomplete yet

! Revision 1.5.2.32  2001/09/10 23:50:47  livesey
! Added a use statement

! Revision 1.5.2.31  2001/09/10 23:34:51  zvi
! Added freq_avg code..

! Revision 1.5.2.30  2001/09/10 21:07:58  livesey
! Interim

! Revision 1.5.2.29  2001/09/10 21:06:42  livesey
! Interim

! Revision 1.5.2.28  2001/09/10 20:46:42  livesey
! Interim

! Revision 1.5.2.27  2001/09/10 20:24:06  livesey
! More tidying up

! Revision 1.5.2.26  2001/09/10 19:56:36  livesey
! Added call to AllocateOneSlabs

! Revision 1.5.2.25  2001/09/10 19:38:19  livesey
! Tidied up variable definitions a bit.

! Revision 1.5.2.24  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90

! Revision 1.5.2.23  2001/09/09 11:17:43  zvi
! Cleaning up..

! Revision 1.5.2.22  2001/09/09 03:16:48  livesey
! End of working day (ha ha!)

! Revision 1.5.2.21  2001/09/09 03:03:25  livesey
! More

! Revision 1.5.2.20  2001/09/09 03:02:53  zvi
! more definitions ..

! Revision 1.5.2.19  2001/09/09 02:42:52  zvi
! more definitions ..

! Revision 1.5.2.18  2001/09/09 02:16:56  livesey
! More..

! Revision 1.5.2.17  2001/09/09 02:18:03  zvi
! more definirions ..

! Revision 1.5.2.16  2001/09/09 01:57:31  livesey
! Work

! Revision 1.5.2.15  2001/09/09 01:57:41  zvi
! more Allocatios ..

! Revision 1.5.2.14  2001/09/09 01:40:18  livesey
! More work

! Revision 1.5.2.13  2001/09/09 01:41:58  zvi
! Allocatios ..

! Revision 1.5.2.12  2001/09/09 00:43:34  zvi
! more metric work

! Revision 1.5.2.11  2001/09/09 00:09:15  livesey
! Another intermediate

! Revision 1.5.2.10  2001/09/08 23:47:52  livesey
! More updates

! Revision 1.5.2.9  2001/09/08 23:11:23  livesey
! Sent to zvi

! Revision 1.5.2.8  2001/09/08 21:41:14  zvi
! Starting to import codes

! Revision 1.5.2.7  2001/09/07 22:49:10  livesey
! Very intermediate

! Revision 1.5.2.6  2001/09/07 22:34:21  zvi
! change comments..

! Revision 1.5.2.5  2001/09/07 22:18:04  livesey
! More comments

! Revision 1.5.2.4  2001/09/07 20:16:37  livesey
! Changed stuff to lower case

! Revision 1.5.2.3  2001/09/07 20:05:26  livesey
! Cosmetic change again.

! Revision 1.5.2.2  2001/09/07 19:58:49  zvi
! Starting new code developement

! Revision 1.5  2001/07/04 00:34:24  zvi
! Modified & new code(s) for better timing

! Revision 1.4  2001/06/22 02:35:08  zvi
! Fixing some memory leaks..

! Revision 1.3  2001/06/21 15:05:42  livesey
! Gets tolerance from fwdModelConf

! Revision 1.2  2001/06/21 13:07:08  zvi
! Speed enhancement MAJOR update

! Revision 1.1  2001/05/29 22:53:51  livesey
! First version, taken from old ForwardModelInterface.f90
