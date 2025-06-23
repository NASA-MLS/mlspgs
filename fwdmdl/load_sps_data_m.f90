! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Load_SPS_Data_M

  use ForwardModelConfig, only: QtyStuff_T
  use MLSCommon, only: R8, RP

  implicit NONE

  private
  public :: Load_Sps_Data, Load_One_Item_Grid, Load_Grid_From_Vector
  public :: Modify_Values_For_Supersat
  public :: Create_Grids_1, Create_Grids_2, Fill_Grids_1, Fill_Grids_2
  public :: FindInGrid, Get_SPS_Bounds
  public :: EmptyGrids_t, Destroygrids_t, Dump, Dump_Grids

  type, public :: C_t ! "Contents" type, to get an array of pointers
    real(rp), pointer, contiguous :: V1(:) => NULL()      ! Frq * Zeta * Phi * Cross
    real(rp), pointer, contiguous :: V4(:,:,:,:) => NULL() ! Frq X Zeta X Phi X Cross
    logical, pointer, contiguous :: L1(:) => NULL()       ! Frq * Zeta * Phi * Cross
    logical, pointer, contiguous :: L4(:,:,:,:) => NULL() ! Frq X Zeta X Phi X Cross
  end type C_t

  type, public :: Grids_T                 ! Fit all Gridding categories
    type(qtyStuff_t), pointer :: QtyStuff(:) => null()
    integer,  pointer :: l_f(:) => null() ! Last entry in frq_basis per sps
    integer,  pointer :: l_p(:) => null() ! Last entry in phi_basis per sps
    integer,  pointer :: l_v(:) => null() ! Last entry in values per sps
    integer,  pointer :: l_x(:) => null() ! Last entry in cross angles per sps
    integer,  pointer :: l_z(:) => null() ! Last entry in zet_basis per sps
    integer,  pointer :: l_zp(:) => null() ! Last entry in Zeta x Phi per sps;
                                          ! Not the same if l_v if there is an
                                          ! earlier frequency-dependent species
    integer :: P_Len = 0 ! \sum_{i=1}^n (l_z(i)-l_z(i-1))*(l_p(i)-l_p(i-1))*
                         !              (l_x(i)-l_x(i-1))
    integer,  pointer :: windowstart(:) => null()! horizontal starting index
                                          ! from l2gp for 2D, or 1 for QTM
    integer,  pointer :: windowfinish(:) => null()! horizontal ending index
                                          ! from l2gp for 2D per sps, or the
                                          ! number of vertices adjacent to the
                                          ! path for QTM
    integer,  pointer :: mol(:) => null()       ! Qty molecule, l_...
    integer,  pointer :: Z_coord(:) => null()   ! l_geo[cd]Altitude, l_zeta
    logical,  pointer :: Coherent(:) => null()  ! From each Qty template
    logical,  pointer :: Stacked(:) => null()   ! From each Qty template
    ! Which column of dBeta_path_df to use for each molecule.
    ! Molecule beta is not dependent upon mixing ratio if zero.
    ! Only nonzero where derivatives for molecules with mixing-ratio-dependent
    ! beta are selected.
    integer,  pointer :: where_dBeta_df(:) => null()
    integer,  pointer :: qty(:) => null() ! Qty type, l_...
    integer,  pointer :: s_ind(:) => null()     ! Indexed by mol, gives index
    ! in l_f etc if qty=l_vmr.
    logical,  pointer :: lin_log(:) => null()   ! true for log representation
                                                ! basis
    real(r8), pointer :: min_val(:) => null()   ! Minimum value
    real(r8), pointer :: frq_basis(:) => null() ! frq  grid entries for all
                                                ! molecules
    real(rp), pointer :: zet_basis(:) => null() ! zeta grid entries for all
                                                ! molecules
    real(rp), pointer :: phi_basis(:) => null() ! phi  grid entries for all
                                                ! molecules, radians
    real(rp), pointer :: cross_angles(:) => null() ! cross-angles grid entries
                                                ! for all molecules, degrees
    real(rp), pointer :: values(:) => null()    ! species values (eg vmr).
    ! For 2D, this is really a four-dimensional quantity dimensioned frequency
    ! (or 1) X zeta (or 1) X phi (or 1) X Cross, taken in Fortran's column-
    ! major array-element order.  For 3D, this is a three-dimensional quantity
    ! dimensioned frequency (or 1) X zeta (or 1) X QTM index, taken in
    ! Fortran's column-major array-element order.
    logical,  pointer :: deriv_flags(:) => null() ! flags to do derivatives,
                                                ! corresponding to the values
    type(c_t), pointer :: C(:) => null()        ! Contents pointers, associated
                                                ! with parts of Values(:) and
                                                ! Deriv_Flags(:).
  contains
    procedure :: IsQTM
  end type Grids_T

  interface Dump
    module procedure Dump_Grids
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! ----------------------------------------------  Load_SPS_Data  -----

  subroutine Load_Sps_Data ( FwdModelConf, Phitan, MAF, Grids_x, &
                           & QtyStuffIn, Subset, Short )

    use ForwardModelConfig, only: ForwardModelConfig_t, QtyStuff_T
    use Toggles, only: Emit, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use VectorsModule, only: VectorValue_T

    type (forwardModelConfig_T), intent(in) :: fwdModelConf
    type (vectorValue_T), intent(in) ::  Phitan  ! Tangent geodAngle
    integer, intent(in) :: MAF

    type (Grids_T), intent(out) :: Grids_x     ! All the coordinates

    type (qtyStuff_t), intent(in), optional , target :: QtyStuffIn(:)
    integer, intent(in), optional :: Subset(:) ! Subset of horizontal basis,
                                               ! primarily for QTM
    logical, intent(in), optional :: Short     ! Don't do Fill_Grids_2

    ! Local variables:

    logical :: Long      ! true if short is absent, else .not. short
    integer :: Me = -1   ! String index cache for tracing
    integer :: Mol, No_Mol

    type (qtyStuff_t), pointer :: QtyStuff(:)

    ! Begin code:

    call trace_begin ( me, 'Load_Sps_Data', &
      & cond=toggle(emit) .and. levels(emit) > 2 ) ! set by -f command-line switch

    qtyStuff => fwdModelConf%beta_group%qty
    if ( present(qtyStuffIn) ) qtyStuff => qtyStuffIn

    no_mol = size( qtyStuff )

    ! Stuff that doesn't depend on the numbers of values or lengths of bases
    call create_grids_1 ( Grids_x, no_mol )

    do mol = 1, no_mol

      grids_x%qtyStuff(mol) = qtyStuff(mol)

      if ( .not. associated(qtyStuff(mol)%qty) ) then
        grids_x%l_f(mol) = grids_x%l_f(mol-1)
        grids_x%l_p(mol) = grids_x%l_p(mol-1)
        grids_x%l_v(mol) = grids_x%l_v(mol-1)
        grids_x%l_z(mol) = grids_x%l_z(mol-1)
        grids_x%s_ind(mol) = 0
      else
        call fill_grids_1 ( grids_x, mol, maf, phitan, fwdModelConf, subset )
      end if

    end do

    ! Allocate space for the zeta, phi, freq. basis and value components.

    call create_grids_2 ( Grids_x )

    long = .true.
    if ( present(short) ) long = .not. short
    if ( long ) then
      do mol = 1, no_mol
        ! Fill the zeta, phi, freq. basis and value components.
        if ( associated(qtyStuff(mol)%qty) ) &
          & call fill_grids_2 ( grids_x, mol, qtyStuff(mol)%qty, &
            & fwdModelConf%moleculeDerivatives(mol) .and. &
            & qtyStuff(mol)%derivOK, subset=subset )
      end do
    end if

! ** ZEBUG - Simulate qty%values for EXTINCTION, using the N2 function
!  (Some code here ...)
! ** END ZEBUG

    call trace_end ( 'Load_Species_Data', &
      & cond=toggle(emit) .and. levels(emit) > 2 )

  end subroutine Load_Sps_Data

  ! --------------------------------------  Load_Grid_From_Vector  -----
  subroutine Load_Grid_From_Vector ( Grids_X, Vector, Phitan, Maf, &
    & FwdModelConf, SetDerivFlags )

    use ForwardModelConfig, only: ForwardModelConfig_t
    use VectorsModule, only: Vector_T, VectorValue_t

    type(grids_t), intent(out) :: Grids_X
    type(vector_t), intent(in) :: Vector
    type(vectorValue_t), intent(in), optional :: Phitan
    integer, intent(in), optional :: Maf
    type(forwardModelConfig_t), intent(in), optional :: FwdModelConf
    logical, intent(in), optional :: SetDerivFlags(:) ! size(vector%quantities)

    integer :: I_Qty

    call create_grids_1 ( Grids_x, size(vector%quantities) )
    do i_qty = 1, size(vector%quantities)

     grids_x%qtyStuff(i_qty)%qty => vector%quantities(i_qty)

      if ( present(phitan) ) then
        call fill_grids_1 ( grids_x, i_qty, maf, phitan, fwdModelConf )
      else
        ! Use the whole phi space for the window
        call fill_grids_1 ( grids_x, i_qty, maf )
      end if
    end do

    ! Allocate space for the zeta, phi, freq. basis and value components.
    call create_grids_2 ( Grids_x )

    do i_qty = 1, size(vector%quantities)
      ! Fill the zeta, phi, freq. basis and value components.
      call fill_grids_2 ( grids_x, i_qty, vector%quantities(i_qty), setDerivFlags(i_qty) )
    end do

  end subroutine Load_Grid_From_Vector

  ! -----------------------------------------  Load_One_Item_Grid  -----
  subroutine Load_One_Item_Grid ( Grids_X, Qty, Maf, Phitan, FwdModelConf, &
    & SetDerivFlags, SetTscatFlag, Across, Subset, Short )
  ! A simplification of Load_Sps_Data to load a grid that has only one
  ! quantity in it.

    use ForwardModelConfig, only: ForwardModelConfig_t
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use VectorsModule, only: VectorValue_T

    type(grids_t), intent(out) :: Grids_X
    type(vectorValue_t), intent(in), target :: Qty ! Assumes actual argument
                                                   ! is a pointer or target
    integer, intent(in) :: Maf
    type(vectorValue_t), intent(in), optional :: Phitan
    type(forwardModelConfig_t), intent(in), optional :: FwdModelConf
    logical, intent(in), optional :: SetDerivFlags
    logical, intent(in), optional :: SetTscatFlag
    logical, intent(in), optional :: Across ! Viewing angle is not in orbit plane
    integer, intent(in), optional :: Subset(:) ! Subset of horizontal basis,
                                               ! primarily for QTM
    logical, intent(in), optional :: Short     ! Don't do Fill_Grids_2

    type(vectorValue_t) :: QtyStuff
    logical :: Long      ! true if short is absent, else .not. short
    logical :: MyAcross, MyFlag
    integer :: II, No_Ang

    if ( .not. qty%template%stacked .or. .not. qty%template%coherent ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
         & 'Cannot load a quantity that is unstacked or incoherent' )

    myAcross = .false.
    if ( present(across) ) myAcross = across
    myFlag = .false.
    if ( present(SetTscatFlag) ) myFlag = SetTscatFlag

    if ( myFlag ) then

       no_ang = FwdModelConf%num_scattering_angles

       qtyStuff%template = qty%template     ! having same template as temp

       call create_grids_1 ( grids_x, no_ang )


       do ii = 1, no_ang
         grids_x%qtyStuff(ii)%qty => qty
         call fill_grids_1 ( grids_x, ii, maf, phitan, fwdModelConf )
       end do

       call create_grids_2 ( grids_x )

       do ii = 1, no_ang
         ! ??? Can this work?  In Fill_Grids_2, the second subscript   ???
         ! ??? range for qtyStuff%values is grids_x%windowStart:       ???
         ! ??? grids_x%windowFinish, but qtyStuff%values will have a   ???
         ! ??? shape of [ size(qty%values,1), 1 ].  Perhaps in the     ???
         ! ??? call to Fill_Grids_1 above, WS and WF should have been  ???
         ! ??? specified with the value II, and then all of qty%values ???
         ! ??? should be passed to Fill_Grids_2 here.                  ???
         qtyStuff%values => qty%values(:,ii:ii)
         call fill_grids_2 ( grids_x, ii, qtyStuff, setDerivFlags )
       end do

    else

       call create_grids_1 ( grids_x, 1 )

       grids_x%qtyStuff(1)%qty => qty
       if ( present(phitan) ) then
         call fill_grids_1 ( grids_x, 1, maf, phitan, fwdModelConf, subset )
       else if ( myAcross ) then
         call MLSMessage ( MLSMSG_Error, moduleName, &
         & 'Cross-track viewing needs PhiTan quantity' )
       else
         ! Use the whole phi space for the window
         call fill_grids_1 ( grids_x, 1, maf, subset=subset )
       end if

       call create_grids_2 ( grids_x )

       long = .true.
       if ( present(short) ) long = .not. short
       if ( long ) &
         & call fill_grids_2 ( grids_x, 1, qty, setDerivFlags, subset=subset )

    end if

  end subroutine Load_One_Item_Grid

  ! ---------------------------------  Modify_Values_For_Supersat  -----

  subroutine Modify_Values_For_Supersat ( FwdModelConf, &
    & Grids_f, H2O_Ind, Grids_Tmp, BoundaryPressure )

    use Constants, only: Deg2Rad
    use ForwardModelConfig, only: ForwardModelConfig_t
    use Intrinsic, only: L_Clear_110RH_Below_Top, L_Clear_0RH
    use HyperSlabs, only: EssentiallyEqual
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSNumerics, only: Hunt
    use RHIFromH2O, only: RHIFromH2O_Factor
    use VectorsModule, only: VectorValue_T

    type (forwardModelConfig_T) ,intent(in) :: FwdModelConf
    type (grids_T), intent(inout) :: Grids_F   ! All the vmrs
    integer, intent(in) :: H2O_Ind             ! Where is H2O in Grids_f?
    type (grids_T), intent(in) :: Grids_Tmp    ! All the temperatures
    type (vectorValue_T), intent(in) :: BoundaryPressure

    integer :: KF, KZ, KP
    real(r8) :: RHI
    integer :: Supersat_Index
    integer :: WF1, WF2
    integer :: V0               ! One before starting point in Values array
    integer :: Z_index
    integer :: Z0, P, P0
    logical :: Failed

    kf = Grids_f%l_f(h2o_ind) - Grids_f%l_f(h2o_ind-1)
    kz = Grids_f%l_z(h2o_ind) - Grids_f%l_z(h2o_ind-1)
    kp = grids_f%l_p(h2o_ind) - grids_f%l_p(h2o_ind-1)

    v0 = Grids_f%l_v(h2o_ind-1)
    wf1 = Grids_f%windowStart(h2o_ind)
    wf2 = Grids_f%windowFinish(h2o_ind)
    z0 = Grids_f%l_z(h2o_ind-1)
    p0 = Grids_f%l_p(h2o_ind-1)

    ! Here, we have to assert that temperature, h2o and boundarypressure share
    ! the same horizontal and vertical (except boundary pressure) grids.
    failed =  kf /= 1
    failed = failed .or. kf /= grids_tmp%l_f(1) - grids_tmp%l_f(0)
    failed = failed .or. kz /= grids_tmp%l_z(1) - grids_tmp%l_z(0)
    failed = failed .or. &
      & grids_tmp%windowStart(1) /= wf1 .or. grids_tmp%windowFinish(1) /= wf2
    if ( failed ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'When overwriting H2O, it must share dimensions with temperature' )
    else
      ! Do more checks
      failed = .not. all ( EssentiallyEqual ( &
        & grids_f%zet_basis(z0+1:z0+kz), grids_tmp%zet_basis ) )
      failed = failed .or. .not. all ( EssentiallyEqual ( &
        & grids_f%phi_basis(p0+1:p0+kp), grids_tmp%phi_basis ) )
      if ( failed ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'When overwriting H2O, it must share coordinates with temperature' )
    end if
    ! Check the boundary pressure
    failed = wf1 < 1 .or. wf2 > boundaryPressure%template%noInstances
    if ( .not. failed ) then
      failed = .not. all ( EssentiallyEqual ( &
        & boundaryPressure%template%phi(1,wf1:wf2)*Deg2Rad, grids_tmp%phi_basis ) )
    end if
    if ( failed ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'When overwriting H2O, boundary pressure must share coordinates with temperature' )

    ! To get height index h, profile index p the values are in:
    ! h2o: grids_f%values ( v0 + h + kz*p )
    ! temperature: grids_tmp%values ( h + kz*p )
    ! tpPres: boundaryPressure%values ( 1, p )
    ! zetas (-log pressure) are given by grids_tmp%zet_basis

    select case (FwdModelConf%i_saturation)
      case (l_clear_110RH_below_top)
        RHi = 110._r8              ! 110% supersaturation
      case (l_clear_0RH)
        RHi = 1.e-9_r8             ! 0% dry saturation
      case default
        call MLSMessage(MLSMSG_Error, ModuleName,'invalid i_saturation')
    end select

    ! H2O may have log basis functions
    do p = 1, kp
      ! find the index for the top of saturation levels
      call Hunt (Grids_f%zet_basis(z0+1:z0+kz), -log10(boundaryPressure%values(1,p)), &
        & supersat_Index, 1, nearest=.true.)

      if (supersat_Index < 1 .or. supersat_Index > kz) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'the top for supersaturation is out of range' )
      else
        if ( Grids_f%lin_log(h2o_ind) ) then
          do z_index = 1, supersat_Index
            Grids_f%values(v0 + z_index + kz*(p-1)) = &
              & log ( RHIFromH2O_Factor ( grids_tmp%values(z_index + kz*(p-1)), &
              &   grids_f%zet_basis(z0+z_index), 0, .true. ) * RHI )
          end do
        else
          do z_index = 1, supersat_Index
            Grids_f%values(v0 + z_index + kz*(p-1)) = &
              & RHIFromH2O_Factor ( grids_tmp%values(z_index + kz*(p-1)), &
              &   grids_f%zet_basis(z0+z_index), 0, .true. ) * RHI
          end do
        end if   ! linear or log basis
      end if     ! check supersat_index
    end do

  end subroutine Modify_Values_For_Supersat

  ! ---------------------------------------------  Create_Grids_1  -----
  subroutine Create_Grids_1 ( Grids_x, N )
  ! Create the part of the Grids_T structure that doesn't depend on the
  ! number of values or lengths of bases

    use Allocate_Deallocate, only: Allocate_test, Test_Allocate
    use Molecules, only: First_Molecule, Last_Molecule

    type (Grids_T), intent(inout) :: Grids_X
    integer, intent(in) :: N

    integer :: Stat

    allocate ( Grids_x%qtyStuff(n), stat=stat )
    call test_allocate ( stat, moduleName, 'Grids_x%QtyStuff' )
    call allocate_test ( Grids_x%l_f, n, 'Grids_x%l_f', moduleName, lowBound=0 )
    call allocate_test ( Grids_x%l_p, n, 'Grids_x%l_p', moduleName, lowBound=0 )
    call allocate_test ( Grids_x%l_v, n, 'Grids_x%l_v', moduleName, lowBound=0 )
    call allocate_test ( Grids_x%l_x, n, 'Grids_x%l_x', moduleName, lowBound=0 )
    call allocate_test ( Grids_x%l_z, n, 'Grids_x%l_z', moduleName, lowBound=0 )
    call allocate_test ( Grids_x%l_zp, n, 'Grids_x%l_zp', moduleName, lowBound=0 )
    call allocate_test ( Grids_x%mol, n, 'Grids_x%mol', moduleName )
    call allocate_test ( Grids_x%qty, n, 'Grids_x%qty', moduleName )
    call allocate_test ( Grids_x%where_dBeta_df, n, 'Grids_x%where_dBeta_df', ModuleName )
    call allocate_test ( Grids_x%s_ind, last_molecule, &
      & 'Grids_x%S_ind', moduleName, lowBound=first_molecule, fill=0 )
    grids_x%l_f(0) = 0
    grids_x%l_p(0) = 0
    grids_x%l_v(0) = 0
    grids_x%l_x(0) = 0
    grids_x%l_z(0) = 0
    grids_x%l_zp(0) = 0
    grids_x%p_len = 0
    grids_x%mol = 0
    grids_x%qty = 0
    grids_x%where_dBeta_df = 0
    call allocate_test ( Grids_x%windowstart, n, 'Grids_x%windowstart', &
                       & ModuleName )
    call allocate_test ( Grids_x%windowfinish, n, 'Grids_x%windowfinish',&
                       & ModuleName )
    call allocate_test ( Grids_x%lin_log, n, 'lin_log', ModuleName )
    call allocate_test ( Grids_x%min_val, n, 'min_val', ModuleName )
    call allocate_test ( Grids_x%z_coord, n, 'Grids_x%z_coord', ModuleName )
    call allocate_test ( Grids_x%coherent, n, 'Grids_x%coherent', ModuleName )
    call allocate_test ( Grids_x%stacked, n, 'Grids_x%stacked', ModuleName )

    grids_x%min_val = -huge(0.0_r8)

  end subroutine Create_Grids_1

  ! ---------------------------------------------  Create_Grids_2  -----
  subroutine Create_Grids_2 ( Grids_x )
  ! Create the part of the Grids_T structure that depends on the
  ! number of values and lengths of bases

    use Allocate_Deallocate, only: Allocate_Test, Test_Allocate

    type (Grids_T), intent(inout) :: Grids_X

    integer :: N, Stat
    character(127) :: ERMSG

    n = ubound(grids_x%l_z,1)

    call allocate_test ( Grids_x%zet_basis, Grids_x%l_z(n), 'Grids_x%zet_basis', &
                       & ModuleName )
    call allocate_test ( Grids_x%phi_basis, Grids_x%l_p(n), 'Grids_x%phi_basis', &
                       & ModuleName )
    call allocate_test ( Grids_x%frq_basis, Grids_x%l_f(n), 'Grids_x%frq_basis', &
                       & ModuleName )
    call allocate_test ( Grids_x%cross_angles, Grids_x%l_x(n), 'Grids_x%cross_angles', &
                       & ModuleName )
    call allocate_test ( Grids_x%values, Grids_x%l_v(n), 'Grids_x%values', &
                       & ModuleName )
    call allocate_test ( Grids_x%deriv_flags, Grids_x%l_v(n), 'Grids_x%deriv_flags', &
                       & ModuleName )
    allocate ( Grids_x%c(n), stat=stat, errmsg=ermsg )
    call test_allocate ( stat, moduleName, 'Grids_x%C', ermsg=ermsg )

  end subroutine Create_Grids_2

  ! -----------------------------------------------  Fill_Grids_1  -----
  subroutine Fill_Grids_1 ( Grids_x, II, Maf, Phitan, FwdModelConf, Subset )
  ! Fill in the size information for the II'th element of Grids_x

    use ForwardModelConfig, only: ForwardModelConfig_t
    use Intrinsic, only: L_Vmr
    use ManipulateVectorQuantities, only: FindInstanceWindow
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use VectorsModule, only: VectorValue_T

    type(grids_T), intent(inout) :: Grids_x
    integer, intent(in) :: II
    integer, intent(in) :: MAF
    type (vectorValue_T), intent(in), optional :: PHITAN  ! Tangent geodAngle
    type(forwardModelConfig_T), intent(in), optional :: FwdModelConf
    integer, intent(in), optional :: Subset(:) ! Subset of horizontal basis,
                                               ! primarily for QTM

    integer :: KF ! Number of frequencies
    integer :: KP ! Number of horizontal coordinates
    integer :: KX ! Number of cross-track coordinates
    integer :: KZ ! Number of vertical coordinates

    associate ( qty => grids_x%qtyStuff(ii)%qty )
 
      if ( .not. qty%template%stacked .or. .not. qty%template%coherent ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
           & 'Cannot load a quantity that is unstacked or incoherent' )

      if ( present(subset) ) then
        grids_x%windowStart(ii) = 1
        grids_x%windowFinish(ii) = size(subset)
      else if ( qty%template%isQTM() ) then
        grids_x%windowStart(ii) = 1
        grids_x%windowFinish(ii) = size(qty%template%the_hGrid%QTM_tree%geo_in)
      else if ( present(phitan) ) then
        call FindInstanceWindow ( qty, phitan, maf, fwdModelConf%phiWindow, &
          & fwdModelConf%windowUnits, grids_x%windowStart(ii), grids_x%windowFinish(ii) )
      else ! Use the whole phi space for the window
        grids_x%windowStart(ii) = 1
        grids_x%windowFinish(ii) = size(qty%template%phi,2)
      end if
      if ( grids_x%windowStart(ii) > grids_x%windowFinish(ii) ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Why is windowStart greater than windowFinish?' )

      kz = qty%template%noSurfs

      kp = grids_x%windowFinish(ii) - grids_x%windowStart(ii) + 1

      kx = qty%template%noCrossTrack ! kx >= 1 even if there is no cross track

      if ( associated(qty%template%frequencies) ) then
        kf = size(qty%template%frequencies)
      else
        kf = qty%template%noChans ! == 1 if qty%template%frequencyCoordinate == l_none
      end if

      grids_x%l_f(ii) = grids_x%l_f(ii-1) + kf
      grids_x%l_p(ii) = grids_x%l_p(ii-1) + &
                      &   kp * merge(1,kz,qty%template%stacked) * kx
      grids_x%l_x(ii) = grids_x%l_x(ii-1) + kx
      grids_x%l_v(ii) = grids_x%l_v(ii-1) + kf * kz * kp * kx
      grids_x%z_coord(ii) = qty%template%verticalCoordinate
      grids_x%l_z(ii) = grids_x%l_z(ii-1) + kz * merge(1,kx,qty%template%stacked)
      grids_x%l_zp(ii) = grids_x%l_zp(ii-1) + kz * kp
      grids_x%p_len = grids_x%p_len + kz * kp * kx
      grids_x%mol(ii) = qty%template%molecule
      grids_x%qty(ii) = qty%template%quantityType
      grids_x%coherent(ii) = qty%template%coherent
      grids_x%stacked(ii) = qty%template%stacked

      ! Remember positions of molecules in grids_x
      if ( qty%template%quantityType == l_vmr ) &
        & grids_x%s_ind(qty%template%molecule) = ii
        ! Note the ambiguity here as to whether it's extinction or extinctionv2:
        ! the last one wins.
        ! Also note, however, that s_ind(l_*extinction*) are never actually used.
    end associate

  end subroutine Fill_Grids_1

  ! -----------------------------------------------  Fill_Grids_2  -----
  subroutine Fill_Grids_2 ( Grids_x, II, Qty, SetDerivFlags, Phi_Offset, Subset )
  ! Fill the zeta, phi, freq. basis and value components for the II'th
  ! "molecule" in Grids_x.

    use Constants, only: Deg2Rad
    use Geometry, only: To_XYZ, XYZ_to_Geod
    use Intrinsic, only: L_Channel, L_geocAltitude, L_geodAltitude, &
      & L_IntermediateFrequency, L_Zeta
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Molecules, only: L_CloudIce, L_H2O
    use MoreMessage, only: MLSMessage
    use To_Log_Basis_m, only: To_Log_Basis
    use VectorsModule, only: M_FullDerivatives, VectorValue_T

    ! For which molecules do we compute dBeta_df?
    !integer, parameter :: Which_dBeta_df(2) = (/ L_CloudIce, L_H2O /)
    integer, parameter :: Which_dBeta_df(1) = (/ L_CloudIce /)

    type(grids_t), intent(inout) :: Grids_x
    integer, intent(in) :: II
    real(r8), pointer :: P3(:,:,:) ! Rank 3 pointer to grids_x%phi_basis(pp+1:qp)
    type(vectorValue_T), intent(in) :: QTY     ! An arbitrary vector quantity
    logical, intent(in), optional :: SetDerivFlags
    real(rp), intent(in), optional :: Phi_Offset ! Radians
    integer, intent(in), optional :: Subset(:) ! Subset of horizontal basis,
                                               ! primarily for QTM

    integer :: H, I, J, K, KF, KP, KX, KZ, L, N, PF, PP, PV, PX, PZ
    integer :: QF, QP, QV, QX, QZ, WF, WS
    integer :: InstOr1
    logical :: PackFrq ! Need to pack the frequency "dimension"
    real(rp) :: Geod(3)   ! Geodetic coordinates lat (deg), lon(deg), alt (m)
    
    pf = Grids_x%l_f(ii-1)
    pp = Grids_x%l_p(ii-1)
    pv = Grids_x%l_v(ii-1)
    px = grids_x%l_x(ii-1)
    pz = Grids_x%l_z(ii-1)

    qf = Grids_x%l_f(ii)
    qp = Grids_x%l_p(ii)
    qv = Grids_x%l_v(ii)
    qx = Grids_x%l_x(ii)
    qz = Grids_x%l_z(ii)

    kf = qf - pf
    kx = qx - px
    kz = qz - pz
    ws = Grids_x%windowStart(ii)
    wf = Grids_x%windowFinish(ii)
    kp = wf - ws + 1

    ! Associate components of Grids_x%c with parts of Grids_x%Values
    ! and Grids_x%Deriv_Flags.

    Grids_x%c(ii)%v1 => Grids_x%values(pv+1:qv)
    Grids_x%c(ii)%v4(1:kf,1:kz,ws:wf,1:kx) => Grids_x%c(ii)%v1
    Grids_x%c(ii)%l1 => Grids_x%deriv_flags(pv+1:qv)
    Grids_x%c(ii)%l4(1:kf,1:kz,ws:wf,1:kx) => Grids_x%c(ii)%l1

    select case ( qty%template%verticalCoordinate )
    case ( l_geocAltitude ) ! This will be used for interpolation with
                            ! H_Path, which is in geodetic height measured
                            ! from the center of the equivalent circular
                            ! Earth, in kilometers, not meters.  We don't add
                            ! the equivalent circular Earth radius here, not
                            ! least because we don't know it yet.
      n = pz
      ! The order for zet_basis is surfs (fastest), instances, cross angles
      do k = px+1, qx
        h = k
        if ( present(subset) ) h = subset(k)
        do j = ws, wf
          instOr1 = merge(1,j,qty%template%coherent)
          do i = 1, qty%template%noSurfs
            l = merge(1,i,qty%template%stacked)
            n = n + 1
            ! Get geodetic coordinates equivalent to geocentric coordinates
            geod = xyz_to_geod ( to_xyz ( qty%template%geodLat3(l,j,h), &
                                        & qty%template%lon3(l,j,h) ) * &
                                        & qty%template%surfs(i,instOr1) )
            ! Get height above the geoid in kilometers
            Grids_x%zet_basis(n) = geod(3)/1000.0
          end do
        end do
      end do
    case ( l_geodAltitude ) ! This will be used for interpolation with
                            ! H_Path, which is in kilometers, not meters
      Grids_x%zet_basis(pz+1:qz) = qty%template%surfs(1:kz,1) / 1000.0
    case ( l_zeta )
      Grids_x%zet_basis(pz+1:qz) = qty%template%surfs(1:kz,1)
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unexpected vertical coordinate for quantity $S', &
        & datum=[qty%template%name] )
    end select

    p3(1:size(qty%template%phi,1),1:kp,1:kx) => Grids_x%phi_basis(pp+1:qp)

    do i = 1, kx ! kx >= 1 even if there are no crossAngles, and
                 ! qty%template%crossAngles == 0 if there are none.
      p3(:,:,i) = (qty%template%phi(:,ws:wf) + qty%template%crossAngles(i)) * Deg2Rad
    end do

    if ( present(phi_offset) ) &
      & Grids_x%phi_basis(pp+1:qp) = Grids_x%phi_basis(pp+1:qp) + phi_offset

    grids_x%cross_angles(px+1:qx) = qty%template%crossAngles

    if ( associated ( qty%template%frequencies ) ) then
      if ( qty%template%frequencyCoordinate /= l_intermediateFrequency .and. &
         & qty%template%frequencyCoordinate /= l_channel ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unexpected frequency coordinate for quantity $S', &
        & datum=[qty%template%name] )
      Grids_x%frq_basis(pf+1:qf) = qty%template%frequencies
    else
      Grids_x%frq_basis(pf+1:qf) = 0.0
    end if

    packFrq = .false.
    if ( associated(qty%template%frequencies) ) &
      & packFrq = size(qty%template%frequencies) /= qty%template%noChans

    if ( packFrq ) then ! Only the values for which we have frequencies
      if ( present(subset) ) then
        Grids_x%c(ii)%v4 = qty%value4(qty%template%chanInds,1:kz,subset,1:kx)
      else
        Grids_x%c(ii)%v4 = qty%value4(qty%template%chanInds,1:kz,ws:wf,1:kx)
      end if
    else ! All the values
      ! kz = noSurfs * noCrossTrack here
      ! shape(qty%value4) = [ freqs, zetas, instances, cross-angles ]
      if ( present(subset) ) then
        Grids_x%c(ii)%v4 = qty%value4(1:kf,1:kz,subset,1:kx)
      else
        Grids_x%c(ii)%v4 = qty%value4(1:kf,1:kz,ws:wf,1:kx)
      end if
    end if

    !  mixing ratio values are manipulated if it's log basis
    Grids_x%lin_log(ii) = qty%template%logBasis
    if ( qty%template%logBasis ) then
      Grids_x%min_val(ii) = qty%template%minValue
      call to_log_basis ( Grids_x%c(ii)%v1, qty%template%minValue )
    end if

    ! set 'do derivative' flags

    if ( present(setDerivFlags) ) then
      if ( setDerivFlags ) then
        ! We set grids_x%where_dBeta_df here instead of in Fill_Grids_1
        ! because we have setDerivFlags here, but not there.
        if ( any(qty%template%molecule == Which_dBeta_df) ) &
          grids_x%where_dBeta_df(ii) = count(grids_x%where_dBeta_df /= 0) + 1
        if ( associated(qty%mask) ) then
          if ( present(subset) ) then
            Grids_x%c(ii)%l4 = &
              & iand( M_FullDerivatives, ichar(qty%mask4(1:kf,1:kz,subset,1:kx) ) ) == 0
          else
            Grids_x%c(ii)%l4 = &
              & iand( M_FullDerivatives, ichar(qty%mask4(1:kf,1:kz,ws:wf,1:kx) ) ) == 0
          end if
        else
          Grids_x%deriv_flags(pv+1:qv) = .true.
        end if
      else
        grids_x%deriv_flags(pv+1:qv) = .false.
      end if
    else
      grids_x%deriv_flags(pv+1:qv) = .false.
    end if

  end subroutine Fill_Grids_2

  ! -------------------------------------------------  FindInGrid  -----
  integer function FindInGrid ( Grids_x, Phi, Zeta, Sps )
    ! Find the index in l_v for Sps, Phi, Zeta if Sps is present.
    ! Find the index in l_v for Phi, Zeta otherwise, assuming Grids_x is temperature.
    ! Assume no frequency dependence.
    ! Return zero if no data for Sps.

    use Constants, only: Pi
    use MLSKinds, only: RP

    type (Grids_T), intent(in) :: Grids_x
    real(rp), intent(in) :: Phi             ! Radians
    real(rp), intent(in) :: Zeta
    integer, intent(in), optional :: Sps    ! A molecule index

    real(rp), parameter :: Pi2 = 2.0*Pi
    integer :: I_Phi, I_Zeta, P0, P1, V0, V1, Z0, Z1

    ! Get bounds in phi, zeta, values, for sps
    if ( present(sps) ) then
      findInGrid = grids_x%s_ind(sps) ! findInGrid is convenient temp here
      if ( findInGrid == 0 ) return
    else
      findInGrid = 1
    end if
    z0 = grids_x%l_z(findInGrid-1)    ! First zeta - 1 for Sps
    z1 = grids_x%l_z(findInGrid)      ! Last zeta for Sps
    p0 = grids_x%l_p(findInGrid-1)    ! First phi - 1 for Sps
    p1 = grids_x%l_p(findInGrid)      ! Last phi for Sps
    v0 = grids_x%l_v(findInGrid-1)    ! First value - 1 for Sps
    v1 = grids_x%l_v(findInGrid)      ! Last value for Sps

    ! Compute location of (phi,zeta) in values
    i_zeta = minloc(abs(grids_x%zet_basis(z0+1:z1)-zeta),1)
    i_phi = minloc(abs(mod(grids_x%phi_basis(p0+1:p1),pi2)-phi),1)
    findInGrid = v0 + i_zeta + (z1-z0) * (i_phi-1)
    if ( findInGrid > v1 ) then
      findInGrid = 0
    else if ( abs(grids_x%zet_basis(i_zeta+z0)-zeta) > 0.25 * &
       & max(abs(grids_x%zet_basis(min(i_zeta+z0+1,z1)) - &
              &  grids_x%zet_basis(i_zeta+z0)), &
         &   abs(grids_x%zet_basis(i_zeta+z0) - &
              &  grids_x%zet_basis(max(i_zeta+z0-1,z0+1)))) ) then
      findInGrid = 0 ! Zeta too far from grid point
    else if ( abs(mod(grids_x%phi_basis(i_phi+p0),pi2)-phi) > 0.25 * &
       & max(abs(grids_x%phi_basis(min(i_phi+p0+1,p1)) - &
              &  grids_x%phi_basis(i_phi+p0)), &
         &   abs(grids_x%phi_basis(i_phi+p0) - &
              &  grids_x%phi_basis(max(i_phi+p0-1,p0+1)))) ) then
      findInGrid = 0 ! Phi too far from grid point
    end if
  end function FindInGrid

! --------------------------------------------------------  IsQTM  -----
  pure logical function IsQTM ( Grid, I )
    class(grids_t), intent(in) :: Grid
    integer, intent(in) :: I
    IsQTM = grid%qtyStuff(i)%qty%template%isQTM()
  end function IsQTM

  ! -----------------------------------------------  EmptyGrids_t  -----
  subroutine EmptyGrids_t ( Grids_x )
    ! Create a grids structure with all empty grids
    use Allocate_Deallocate, only: Allocate_test, Test_Allocate

    type (Grids_T), intent(inout) :: Grids_x

    integer :: Stat

    grids_x%p_len = 0
    allocate ( Grids_x%qtyStuff(0), stat=stat )
    call test_allocate ( stat, moduleName, 'Grids_x%QtyStuff' )
    call allocate_test(grids_x%l_f,0,'grids_x%l_f',modulename)
    call allocate_test(grids_x%l_z,0,'grids_x%l_z',modulename)
    call allocate_test(grids_x%l_p,0,'grids_x%l_p',modulename)
    call allocate_test(grids_x%l_x,0,'grids_x%l_x',modulename)
    call allocate_test(grids_x%l_v,0,'grids_x%l_v',modulename)
    call allocate_test(grids_x%windowstart,0,'grids_x%windowstart',modulename)
    call allocate_test(grids_x%windowfinish,0,'grids_x%windowfinish',modulename)
    call allocate_test(grids_x%mol,0,'grids_x%mol',modulename)
    call allocate_test(grids_x%where_dBeta_df,0,'grids_x%where_dBeta_df',modulename)
    call allocate_test(grids_x%qty,0,'grids_x%qty',modulename)
    call allocate_test(grids_x%s_ind,0,'grids_x%s_ind',modulename)
    call allocate_test(grids_x%lin_log,0,'grids_x%lin_log',modulename)
    call allocate_test(grids_x%min_val,0,'grids_x%min_val',modulename)
    call allocate_test(grids_x%frq_basis,0,'grids_x%frq_basis',modulename)
    call allocate_test(grids_x%zet_basis,0,'grids_x%zet_basis',modulename)
    call allocate_test(grids_x%phi_basis,0,'grids_x%phi_basis',modulename)
    call allocate_test(grids_x%cross_angles,0,'grids_x%cross_angles',modulename)
    call allocate_test(grids_x%values,0,'grids_x%values',modulename)
    call allocate_test(grids_x%deriv_flags,0,'grids_x%deriv_flags',modulename)
    call allocate_test(grids_x%z_coord,0,'grids_x%z_coord',modulename)
    call allocate_test(grids_x%coherent,0,'grids_x%coherent',modulename)
    call allocate_test(grids_x%stacked,0,'grids_x%stacked',modulename)
    allocate ( Grids_x%c(0), stat=stat )
    call test_allocate ( stat, moduleName, 'Grids_x%C' )

  end subroutine EmptyGrids_t

  ! ---------------------------------------------  Get_SPS_Bounds  -----
  subroutine Get_SPS_Bounds ( Grids_x, SPS, V1, V2 )

    ! Get the bounds in grids_x%values, etc., for SPS

    type (Grids_T), intent(in) :: Grids_x
    integer, intent(in) :: SPS
    integer, intent(out) :: V1, V2

    v1 = grids_x%l_v(sps-1)+1
    v2 = grids_x%l_v(sps)

  end subroutine Get_SPS_Bounds

  ! ---------------------------------------------  DestroyGrids_t  -----
  subroutine DestroyGrids_t ( grids_x )
    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate

    type (Grids_T), intent(inout) :: Grids_x

    character(127) :: ERMSG
    integer :: Stat

    grids_x%p_len = 0
    if ( associated(Grids_x%qtyStuff) ) then
      deallocate ( Grids_x%qtyStuff, stat=stat, errmsg=ermsg )
      call test_deallocate ( stat, moduleName, 'Grids_x%QtyStuff', ermsg=ermsg )
    end if
    call deallocate_test(grids_x%l_f,'Grids_x%l_f',modulename)
    call deallocate_test(grids_x%l_z,'Grids_x%l_z',modulename)
    call deallocate_test(grids_x%l_p,'Grids_x%l_p',modulename)
    call deallocate_test(grids_x%l_x,'Grids_x%l_x',modulename)
    call deallocate_test(grids_x%l_v,'Grids_x%l_v',modulename)
    call deallocate_test(grids_x%windowstart,'Grids_x%windowstart',modulename)
    call deallocate_test(grids_x%windowfinish,'Grids_x%windowfinish',modulename)
    call deallocate_test(Grids_x%mol,'Grids_x%mol',modulename)
    call deallocate_test(Grids_x%where_dBeta_df,'Grids_x%where_dBeta_df',modulename)
    call deallocate_test(Grids_x%qty,'Grids_x%qty',modulename)
    call deallocate_test(Grids_x%s_ind,'Grids_x%s_ind',modulename)
    call deallocate_test(grids_x%lin_log,'Grids_x%lin_log',modulename)
    call deallocate_test(grids_x%min_val,'Grids_x%min_val',modulename)
    call deallocate_test(grids_x%frq_basis,'Grids_x%frq_basis',modulename)
    call deallocate_test(grids_x%zet_basis,'Grids_x%zet_basis',modulename)
    call deallocate_test(grids_x%phi_basis,'Grids_x%phi_basis',modulename)
    call deallocate_test(grids_x%cross_angles,'Grids_x%cross_angles',modulename)
    call deallocate_test(grids_x%values,'Grids_x%values',modulename)
    call deallocate_test(grids_x%deriv_flags,'Grids_x%deriv_flags',modulename)
    call deallocate_test(grids_x%z_coord,'Grids_x%z_coord',modulename)
    call deallocate_test(grids_x%coherent,'Grids_x%coherent',modulename)
    call deallocate_test(grids_x%stacked,'Grids_x%stacked',modulename)
    if ( associated(Grids_x%c) ) then
      deallocate ( Grids_x%c, stat=stat, errmsg=ermsg )
      call test_deallocate ( stat, moduleName, 'Grids_x%C', ermsg=ermsg )
    end if

  end subroutine DestroyGrids_t

  ! -------------------------------------------------  Dump_Grids  -----
  subroutine Dump_Grids ( The_Grid, Name, Details, OneGrid )
  ! Dump The_Grid

    use Constants, only: Rad2Deg
    use Dump_0, only: Dump
    use Intrinsic, only: Lit_Indices
    use Output_M, only: NewLine, Output
    use String_Table, only: Display_String

    type(grids_t), intent(in) :: The_Grid
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: Details ! <= 0 => don't dump bases (default)
                                             ! <= 1 => don't dump values
    integer, intent(in), optional :: OneGrid ! Only dump this grid

    integer :: I, I1, IN, KX, L1, MyDetails, NZG, NZV, W

    myDetails = 0
    if ( present(details) ) myDetails = details

    call output ( 'Dump of Grids_T structure', advance='no' )
    if ( present(name) ) then
      call output ( ' ', advance='no' )
      call output ( name, advance='no' )
    end if
    call output ( the_grid%p_len, before = ', P_Len = ', advance='yes' )
    call output ( 'Molecules:', advance='yes' )
    i1 = 1
    in = size(the_grid%qtyStuff)
    if ( present(oneGrid) ) then
      i1 = oneGrid
      in = oneGrid
    end if
    l1 = max(i1-1,1)
    do i = i1, in
      call output ( i )
      if ( the_grid%qtyStuff(i)%qty%template%name > 0 ) then
        call display_string ( the_grid%qtyStuff(i)%qty%template%name, before=': ' )
      else
        call output ( ': ?' )
      end if
      if ( the_grid%mol(i) > 0 ) then
        call display_string ( lit_indices(the_grid%mol(i)), &
          & before=' molecule: ' )
      end if
      if ( the_grid%qty(i) > 0 ) then
        call display_string ( lit_indices(the_grid%qty(i)), &
          & before=' quantity: ' )
      end if
      call display_string ( lit_indices(the_grid%z_coord(i)), &
        & before=' vertical coordinate: ' )
      if ( the_grid%where_dBeta_df(i) > 0 ) then
        call output ( the_grid%where_dBeta_df(i), before=' dBeta_df at ' )
      end if
      call newLine
    end do
    call dump ( the_grid%l_f(l1:in), 'The_grid%l_f' )
    call dump ( the_grid%l_z(l1:in), 'The_grid%l_z' )
    call dump ( the_grid%l_p(l1:in), 'The_grid%l_p' )
    call dump ( the_grid%l_x(l1:in), 'The_grid%l_x' )
    call dump ( the_grid%l_v(l1:in), 'The_grid%l_v' )
    call dump ( the_grid%l_zp(l1:in), 'The_grid%l_zp' )
    call dump ( the_grid%windowStart(i1:in), 'The_grid%WindowStart' )
    call dump ( the_grid%windowFinish(i1:in), 'The_grid%WindowFinish' )
    call dump ( the_grid%lin_log(i1:in), 'The_grid%Lin_Log' )
    if ( myDetails > 0 ) then
      call dump ( the_grid%min_val(i1:in), 'The_grid%Min_Val' )
      call dump ( the_grid%frq_basis(the_grid%l_f(l1)+1:the_grid%l_f(in)), &
        & 'The_grid%Frq_Basis' )
      call dump ( the_grid%coherent(i1:in), 'The_grid%Coherent' )
      call dump ( the_grid%stacked(i1:in), 'The_grid%Stacked' )
      do i = i1, in
        kx = the_grid%l_x(i) - the_grid%l_x(i-1)
        ! NZ for geolocation fields:
        nzg = merge(1,(the_grid%l_z(i)-the_grid%l_z(i-1))/kx,the_grid%stacked(i))
        ! NZ for values fields:
        nzv = (the_grid%l_z(i)-the_grid%l_z(i-1))/merge(1,kx,the_grid%stacked(i))
        w = the_grid%windowFinish(i) - the_grid%windowStart(i) + 1
        call output ( i, before='The_grid%Phi_Basis(' )
        call dump_what
        call output ( nzg, before=') surfs,insts,cross (' )
        call output ( w, before=',' )
        call output ( the_grid%l_x(i)-the_grid%l_x(i-1), before=',', &
                    & after=') (degrees)', advance='yes' )
        call dump ( rad2deg*the_grid%phi_basis(the_grid%l_p(i-1)+1: &
                                             & the_grid%l_p(i)), &
                  & lbound=the_grid%l_p(i-1)+1 )
        call output ( i, before='The_grid%Cross_Angles(' )
        call dump_what
        call output ( kx, before=') cross (', after=') (degrees)', advance='yes' )
        call dump ( the_grid%cross_angles(the_grid%l_x(i-1)+1:the_grid%l_x(i)), &
                  & lbound=the_grid%l_x(i-1)+1 )
        call output ( i, before='The_grid%Zet_Basis(' )
        call dump_what
        call output ( nzv, before=') surfs,insts,cross (' )
        call output ( w, before=',' )
        call output ( merge(1,kx,the_grid%stacked(i)), before=',', after=')', advance='yes' )
        call dump ( the_grid%zet_basis(the_grid%l_z(i-1)+1:the_grid%l_z(i)), &
          & lbound=the_grid%l_z(i-1)+1 )
        if ( myDetails > 1 ) then
           call output ( i, before='The_grid%Values(' )
           call dump_what
           call output ( the_grid%l_f(i)-the_grid%l_f(i-1), &
             & before=') freqs,surfs,insts,cross (' )
           call output ( nzv, before=',' )
           call output ( w, before=',' )
           call output ( kx, before=',', after=')', advance='yes' )
           call dump ( the_grid%values(the_grid%l_v(i-1)+1:the_grid%l_v(i)), &
             lbound=the_grid%l_v(i-1)+1 )
           call output ( i, before='The_grid%deriv_flags(' )
           call dump_what
           call output ( the_grid%l_f(i)-the_grid%l_f(i-1), &
             & before=') freqs,surfs,insts,cross (' )
           call output ( nzv, before=',' )
           call output ( w, before=',' )
           call output ( kx, before=',', after=')', advance='yes' )
           call dump ( the_grid%deriv_flags(the_grid%l_v(i-1)+1:the_grid%l_v(i)), &
             lbound=the_grid%l_v(i-1)+1 )
        end if
      end do
    end if

  contains
    subroutine Dump_What
      if ( the_grid%qtyStuff(i)%qty%template%name > 0 ) then
        call display_string ( the_grid%qtyStuff(i)%qty%template%name, &
          & before='=' )
      else if ( the_grid%mol(i) > 0 ) then
        call display_string ( lit_indices(the_grid%mol(i)), before='=' )
      else if ( the_grid%qty(i) > 0 ) then
        call display_string ( lit_indices(the_grid%qty(i)), &
          & before='=' )
      end if
    end subroutine Dump_What

  end subroutine Dump_Grids

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module LOAD_SPS_DATA_M

! $Log$
! Revision 2.121  2017/11/03 20:57:45  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.120  2017/03/11 00:57:54  vsnyder
! Make components of C_t contiguous.  Add Subset and Short arguments.  Use
! components of C_t in more places.
!
! Revision 2.119  2017/03/02 00:31:58  vsnyder
! Add %L1, %V1 components of C_t, add Get_SPS_Bounds
!
! Revision 2.118  2017/02/04 02:13:56  vsnyder
! Add type C_t and component C of that type, having components V4 and L4 to
! view sections of Values and Deriv_Flags as rank-4 objects.
! Use To_Log_Basis.
!
! Revision 2.117  2017/01/21 03:09:39  vsnyder
! Spiff the dump
!
! Revision 2.116  2016/11/02 01:30:46  vsnyder
! Remove unused USE name
!
! Revision 2.115  2016/09/21 00:14:20  vsnyder
! Use IsQTM function
!
! Revision 2.114  2016/08/30 20:29:36  vsnyder
! Add IsQTM type-bound function
!
! Revision 2.113  2016/08/23 00:43:11  vsnyder
! Components within or adjacent to the polygon are now within the QTM_Tree_t
! structure instead of the HGrid_t structure.
!
! Revision 2.112  2016/06/03 23:44:05  vsnyder
! Eliminate QTM_Geo component.  Make sure grids_x%qtyStuff(:)%qty is
! associated with the quantity.  Correct some labels in DestroyGrids_t.
! Allow to dump only grid.
!
! Revision 2.111  2016/05/27 01:26:01  vsnyder
! Cannonball polishing
!
! Revision 2.110  2016/05/16 23:22:30  vsnyder
! Check for stacked and coherent always, not just for QTM
!
! Revision 2.108  2016/05/10 00:11:33  vsnyder
! Copy qtyStuff into Grids_x, simplify Create_Grids_1
!
! Revision 2.107  2016/05/10 00:00:32  vsnyder
! Test Grids_x%qtyStuff before deallocating it
!
! Revision 2.106  2016/05/02 23:31:32  vsnyder
! Add QtyStuff, horizontal grid type (phi or QTM), some cannonball polishing
!
! Revision 2.105  2015/10/28 00:32:12  vsnyder
! Check that windowStart is not after windowFinish
!
! Revision 2.104  2015/09/22 23:18:23  vsnyder
! Add a cross-track horizontal coordinate; spiff the dump
!
! Revision 2.103  2015/08/25 18:42:03  vsnyder
! Filling zet_basis and values was still getting subscript bounds errors.
! Hopefully, that's repaired now.  The Grids_t dump is improved.
!
! Revision 2.102  2015/07/27 22:37:22  vsnyder
! Store crossAngles in grids_t instead of conflating with phi.  Use absence
! of PhiTan to indicate the phi window is all available phi's.  Remove the
! assumption that if there are cross angles, there's only one profile.  Do
! not fake using WS and WF for crossAngles bounds.
!
! Revision 2.101  2015/07/08 01:25:35  vsnyder
! Repair? problem with values being shifted
!
! Revision 2.92  2014/08/01 01:03:45  vsnyder
! Eliminate unreferenced USE name
!
! Revision 2.91  2014/07/18 23:16:28  pwagner
! Aimed for consistency in names passed to allocate_test
!
! Revision 2.90  2014/04/22 00:08:27  vsnyder
! Add tracing
!
! Revision 2.89  2013/08/08 02:36:03  vsnyder
! Use derivOK instead of foundInFirst to set deriv_flags
!
! Revision 2.88  2012/08/08 20:08:43  vsnyder
! Insert some comments about potential problems in Load_One_Item_Grid
!
! Revision 2.87  2011/11/11 00:40:17  vsnyder
! Update a comment about extinction
!
! Revision 2.86  2011/08/25 22:37:36  vsnyder
! Delete s_ind dummy argument, since it is never used.  The s_ind component
! of grids_x is used instead.  Move filling grids_x%s_ind to fill_grids_1,
! instead of doing it at some of the calls to fill_grids_1.
!
! Revision 2.85  2011/08/20 02:07:34  vsnyder
! Allocate %s_ind with bounds wide enough to accomodate l_RHi, which is not a molecule
!
! Revision 2.84  2011/08/20 00:45:33  vsnyder
! Remove unused USE statements and declarations for unused variables
!
! Revision 2.83  2011/07/29 01:50:20  vsnyder
! Add the s_ind field.  Make CloudIce a molecule.  Add the FindInGrid function.
!
! Revision 2.82  2011/07/08 20:58:18  yanovsky
! Remove L_H2O from Which_dBeta_df
!
! Revision 2.81  2011/04/28 00:14:24  vsnyder
! Give default value to P_Len component
!
! Revision 2.80  2011/02/12 03:02:04  vsnyder
! Calculate which column of dBeta_df to use
!
! Revision 2.79  2011/02/05 01:17:17  vsnyder
! Add mol component, some cannonball polishing
!
! Revision 2.78  2010/12/07 01:06:24  vsnyder
! Many changes for TScat.  Mostly making irrelevant arguments optional, and
! providing other stuff explicitly.  Also added the qty component.
!
! Revision 2.77  2010/09/25 01:13:35  vsnyder
! Add Load_Grid_From_Vector, Pack_Frq
!
! Revision 2.76  2010/06/07 23:22:53  vsnyder
! Replaced H2O_ind with S_Ind
!
! Revision 2.75  2009/12/22 02:13:38  vsnyder
! Add species names
!
! Revision 2.74  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.73  2009/05/13 20:03:02  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.72  2009/01/16 23:38:32  vsnyder
! Correct problem with absent molecules.  Dump angles in degrees.  Add
! PRINT statement to not_used_here.
!
! Revision 2.71  2008/10/03 16:30:31  livesey
! Added EXTINCTIONV2
!
! Revision 2.70  2008/06/06 22:51:44  pwagner
! EssentiallyEqual moved to MLSFillValues
!
! Revision 2.69  2008/05/20 00:15:07  vsnyder
! Change GRIDS_TMP to INTENT(IN)
!
! Revision 2.68  2006/12/04 21:17:28  vsnyder
! Reorganize FullForwardModel to use automatic arrays instead of allocating
! pointer arrays.  Requires testing for zero size instead of testing for
! associated in several subsidiary procedures.
!
! Revision 2.67  2006/07/21 00:18:09  vsnyder
! Remove unused USEs
!
! Revision 2.66  2006/04/05 21:46:44  vsnyder
! Allow state vector not to include all molecules
!
! Revision 2.65  2006/01/05 03:26:09  vsnyder
! Remove unused variable
!
! Revision 2.64  2005/09/16 23:41:45  vsnyder
! Spiff up a dump
!
! Revision 2.63  2005/08/03 18:04:09  vsnyder
! Some spectroscopy derivative stuff
!
! Revision 2.62  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.61  2004/11/01 20:26:36  vsnyder
! Reorganization of representation for molecules and beta groups; PFA may be broken for now
!
! Revision 2.60  2004/08/03 02:24:21  vsnyder
! Polish up dump a little bit
!
! Revision 2.59  2004/07/08 21:00:23  vsnyder
! Inching toward PFA
!
! Revision 2.58  2004/03/30 00:55:30  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.57  2004/03/20 01:15:29  jonathan
!  minor changes
!
! Revision 2.56  2004/03/01 19:21:46  jonathan
! modify load_one_item for scat source function
!
! Revision 2.55  2004/02/03 02:46:39  vsnyder
! Make ScatFlag argument optional
!
! Revision 2.54  2004/01/23 19:13:10  jonathan
! add SetTscatFlag and related changes
!
! Revision 2.53  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.52  2003/05/16 23:52:26  livesey
! Now uses molecule indices rather than spectags
!
! Revision 2.51  2003/05/16 02:47:12  vsnyder
! Removed USE's for unreferenced symbols
!
! Revision 2.50  2003/05/08 23:42:50  livesey
! Bug fix for requesting derivatives etc.
!
! Revision 2.49  2003/05/06 23:36:59  dwu
! fix a bug in modify_h2o
!
! Revision 2.48  2003/05/06 20:35:34  livesey
! Another bug fix
!
! Revision 2.47  2003/05/06 20:23:21  livesey
! Bug fixes and cosmetic changes, renamed some variables etc.
!
! Revision 2.46  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.42.2.7  2003/03/21 02:48:54  vsnyder
! Get QtyStuff_T from ForwardModelVectorTools, embellish some comments
!
! Revision 2.42.2.6  2003/03/20 19:21:05  vsnyder
! More futzing with grids_t and stuff that uses it
!
! Revision 2.42.2.5  2003/03/20 01:42:26  vsnyder
! Revise Grids_T structure
!
! Revision 2.42.2.4  2003/03/08 00:03:08  vsnyder
! Split some cloud stuff into a separate routine
!
! Revision 2.42.2.3  2003/03/06 21:53:09  vsnyder
! Cosmetic changes
!
! Revision 2.42.2.2  2003/02/14 18:11:43  pwagner
! Non-supersat bug fixed
!
! Revision 2.42.2.1  2003/02/13 18:33:18  dwu
! a clean up
!
! Revision 2.42  2003/02/13 00:42:46  dwu
! fix bugs and add comments for i_saturation
!
! Revision 2.41  2003/02/13 00:07:22  jonathan
! another bug fix
!
! Revision 2.40  2003/02/12 23:50:49  jonathan
! add optional to i_supersat and temp_supersat
!
! Revision 2.39  2003/02/11 00:48:25  jonathan
! fix index bug
!
! Revision 2.38  2003/02/07 18:47:47  vsnyder
! Back to 2.35
!
! Revision 2.35  2003/02/07 03:30:37  vsnyder
! Correct i_supersat test?
!
! Revision 2.34  2003/02/07 02:00:48  vsnyder
! Move some USE statements down
!
! Revision 2.33  2003/02/07 01:07:53  jonathan
! add in option to compute dry and super-saturation case in load_sps
!
! Revision 2.32  2003/02/06 20:00:06  vsnyder
! Make Load_One_Grid a public module procedure instead of internal
!
! Revision 2.31  2003/01/26 04:42:42  livesey
! Added profiles/angle options for phiWindow
!
! Revision 2.30  2002/11/23 02:49:33  vsnyder
! Cosmetic changes
!
! Revision 2.29  2002/10/08 17:08:05  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.28  2002/10/03 05:37:53  livesey
! Minor efficiency improvement in derivative selection
!
! Revision 2.27  2002/10/02 22:42:39  vsnyder
! Move USE statements from module scope to procedure scope.  Make Load_One_Grid
! an internal subroutine of Load_SPS_Data.  Cosmetic changes.
!
! Revision 2.26  2002/09/26 18:02:10  livesey
! Now uses GetQuantityForForwardModel.
!
! Revision 2.25  2002/09/24 21:37:01  livesey
! Added min_val stuff
!
! Revision 2.24  2002/09/05 20:48:59  livesey
! Added moleculeDerivatives info into deriv_flags
!
! Revision 2.23  2002/08/22 23:13:45  livesey
! New frequency basis on IF
!
! Revision 2.22  2002/08/20 22:37:34  livesey
! Minor change in handling of frequency coordinate
!
! Revision 2.21  2002/07/08 17:45:41  zvi
! Updated spectroscopy handling
!
! Revision 2.20  2002/07/05 07:52:49  zvi
! Some cosmetic changes
!
! Revision 2.19  2002/06/19 11:00:34  zvi
! Removing unused variables, some cosmetic changes
!
! Revision 2.17  2002/06/13 22:39:12  bill
! fixed phi window selection--wgr
!
! Revision 2.16  2002/06/04 10:28:03  zvi
! Adding comments, fixing a bug with species ruuning index
!
! Revision 2.15  2002/02/20 22:19:46  zvi
! Reversing the subset logic ..
!
! Revision 2.14  2002/02/16 20:43:54  zvi
! Changing the code for Log of Neg. VMR
!
! Revision 2.13  2002/02/16 10:32:17  zvi
! Guaranties against taking Log(0.0)
!
! Revision 2.12  2002/02/16 06:37:34  zvi
! New code for derivative flags..
!
! Revision 2.11  2002/02/08 00:46:04  zvi
! Fixing a bug in the t_deriv_flag code
!
! Revision 2.10  2002/01/30 01:11:20  zvi
! Fix bug in user selectable coeff. code
!
! Revision 2.9  2002/01/27 08:37:49  zvi
! Adding Users selected coefficients for derivatives
!
! Revision 2.8  2002/01/09 00:30:48  zvi
! Fix a bug with skip_eta_frq
!
! Revision 2.7  2002/01/08 01:02:54  livesey
! Made my_catalog intent(in) rather than pointer, also 'fixed'
! problem with frequency coordinate?
!
! Revision 2.6  2001/12/14 23:43:24  zvi
! Modification for Grouping concept
!
! Revision 2.5  2001/11/15 01:21:59  zvi
! Extiction debug fix
!
! Revision 2.4  2001/11/10 00:46:40  zvi
! Adding the EXTINCTION capabilitis
!
! Revision 2.3  2001/11/08 09:56:59  zvi
! Fixing a bug..
!
! Revision 2.2  2001/11/08 00:10:13  livesey
! Interim version for extinction
!
! Revision 2.1  2001/11/02 10:48:39  zvi
! Implementing frequecy grid
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.1.2.10  2001/09/12 21:38:51  zvi
! Added CVS stuff
!
! Revision 1.1.2.9  2001/09/08 23:11:41  zvi
! Bug fixed..
!
! Revision 1.1.2.8  2001/09/08 22:43:37  zvi
! Bug fixed..
!
! Revision 1.1.2.7  2001/09/08 22:37:11  zvi
! Eliminatin molecule coeff. stuff
!
! Revision 1.1.2.6  2001/09/08 20:22:00  zvi
! Fixing a bug..
!
! Revision 1.1.2.5  2001/09/08 20:19:48  zvi
! Developing code..
!
! Revision 1.1.2.4  2001/09/08 18:30:49  zvi
! Working on developement version
!
! Revision 1.1.2.3  2001/09/07 20:41:42  livesey
! Messing around.
!
! Revision 1.1.2.2  2001/09/07 20:16:28  livesey
! Changed stuff to lower case
!
! Revision 1.1.2.1  2001/09/07 19:56:41  zvi
! New module..
!
! Revision 1.0  2001/09/06 13:07:09  zvi
