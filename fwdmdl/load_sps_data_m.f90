! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module LOAD_SPS_DATA_M

  use MLSCommon, only: R8, RP

  implicit NONE

  private
  public :: Load_Sps_Data, Load_One_Item_Grid, Modify_Values_For_Supersat
  public :: Create_Grids_1, Create_Grids_2, Fill_Grids_1, Fill_Grids_2
  public :: Destroygrids_t, Dump, Dump_Grids

  type, public :: Grids_T                 ! Fit all Gridding categories
    integer,  pointer :: l_f(:) => null() ! Last entry in frq. grid per sps
    integer,  pointer :: l_z(:) => null() ! Last entry in zeta grid per sps
    integer,  pointer :: l_p(:) => null() ! Last entry in phi  grid per sps
    integer,  pointer :: l_v(:) => null() ! Last entry in values grid per sps
    integer :: P_Len ! \sum_i=1^n (l_z(i)-l_z(i-1))*(l_p(i)-l_p(i-1))
    integer,  pointer :: windowstart(:) => null()! horizontal starting index
!                                                  from l2gp
    integer,  pointer :: windowfinish(:) => null()! horizontal ending index
!                                                   from l2gp
    logical,  pointer :: lin_log(:) => null()   ! type of representation basis
    real(r8), pointer :: min_val(:) => null()   ! Minimum value
    real(r8), pointer :: frq_basis(:) => null() ! frq grid entries for all
!                                                 molecules
    real(rp), pointer :: zet_basis(:) => null() ! zeta grid entries for all
!                                                 molecules
    real(rp), pointer :: phi_basis(:) => null() ! phi  grid entries for all
!                                                 molecules
    real(rp), pointer :: values(:) => null()    ! species values (ie vmr). 
!     This is really a three-dimensional quantity dimensioned frequency (or
!     1) X zeta (or 1) X phi (or 1), taken in Fortran's column-major
!     array-element order.
    logical,  pointer :: deriv_flags(:) => null() ! do derivatives flags
!     corresponding to the values component.
  end type Grids_T

  interface Dump
    module procedure Dump_Grids
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
  private :: NOT_USED_HERE 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! ----------------------------------------------  Load_SPS_Data  -----

  subroutine Load_Sps_Data ( FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
       &    radiometer, mol_cat_index, QuantityType, Grids_x,      &
       &    h2o_ind, ext_ind )

    use ForwardModelConfig, only: ForwardModelConfig_t
    use ForwardModelIntermediate, only: ForwardModelStatus_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel, QtyStuff_T
    use Intrinsic, only: L_Phitan
    use ManipulateVectorQuantities, only: FindInstanceWindow
    use Molecules, only: Spec_tags, SP_Extinction, SP_H2O
    use VectorsModule, only: Vector_T, VectorValue_T

    type(forwardModelConfig_T), intent(in) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(forwardModelStatus_t), intent(in) :: FmStat ! Reverse comm. stuff

    integer, intent(in) :: RADIOMETER       
    integer, intent(in) :: MOL_CAT_INDEX(:) 

    integer, intent(in) :: QuantityType      ! L_vmr, L_dv, etc.

    type (Grids_T), intent(out) :: Grids_x   ! All the coordinates

    integer, intent(out), optional :: H2O_IND
    integer, intent(out), optional :: EXT_IND

! Local variables:


    integer :: ii, kk, no_mol, mol, my_ext_ind, my_h2o_ind

    type (VectorValue_T), pointer :: PHITAN  ! Tangent geodAngle component of

    type(qtyStuff_t) :: QtyStuff( size( mol_cat_index ) )

! Begin code:

    no_mol = size( mol_cat_index )

    call create_grids_1 ( Grids_x, no_mol )

    grids_x%min_val = -huge(0.0_r8)

    grids_x%p_len = 0
    my_ext_ind = 0
    my_h2o_ind = 0

    phitan => GetQuantityforForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_phitan, config=fwdModelConf, &
      & instrumentModule=FwdModelConf%signals(1)%instrumentModule )

    do ii = 1, no_mol

      mol = mol_cat_index(ii)
      kk = FwdModelConf%molecules(mol)
      if ( spec_tags(kk) == sp_h2o ) my_h2o_ind = ii        ! memorize h2o index
      if ( spec_tags(kk) == sp_extinction ) my_ext_ind = ii ! memorize extiction ix
      
      qtyStuff(ii)%qty => GetQuantityforForwardModel ( &
        & fwdModelIn, fwdModelExtra, &
        & quantityType=quantityType, molIndex=mol, &
        & radiometer=radiometer, config=fwdModelConf, &
        & foundInFirst=qtyStuff(ii)%foundInFirst )

      call fill_grids_1 ( grids_x, ii, qtyStuff(ii)%qty, phitan, fmStat%maf, &
        &                 fwdModelConf )

    end do

    ! Allocate space for the zeta, phi, freq. basis and value components.

    call create_grids_2 ( Grids_x )

    do ii = 1, no_mol

      ! Fill the zeta, phi, freq. basis and value components.
      call fill_grids_2 ( grids_x, ii, qtyStuff(ii)%qty, &
        & fwdModelConf%moleculeDerivatives(mol) .and. qtyStuff(ii)%foundInFirst )

    end do

! ** ZEBUG - Simulate qty%values for EXTINCTION, using the N2 function
!  (Some code here ...)
! ** END ZEBUG

    if ( present(ext_ind) ) ext_ind = my_ext_ind
    if ( present(h2o_ind) ) h2o_ind = my_h2o_ind

  end subroutine Load_Sps_Data

  ! -----------------------------------------  Load_One_Item_Grid  -----
  subroutine Load_One_Item_Grid ( Grids_X, Qty, Phitan, Maf, FwdModelConf, &
    & SetDerivFlags )
  ! A simplification of Load_Sps_Data to load a grid that has only one
  ! quantity in it.

    use ForwardModelConfig, only: ForwardModelConfig_t
    use VectorsModule, only: VectorValue_T

    type(grids_t), intent(out) :: Grids_X
    type(vectorValue_t), intent(in) :: Qty, Phitan
    integer, intent(in) :: Maf
    type(forwardModelConfig_t), intent(in) :: FwdModelConf
    logical, intent(in) :: SetDerivFlags

    call create_grids_1 ( grids_x, 1 )
    call fill_grids_1 ( grids_x, 1, qty, phitan, maf, fwdModelConf )
    call create_grids_2 ( grids_x )
    call fill_grids_2 ( grids_x, 1, qty, setDerivFlags )

  end subroutine Load_One_Item_Grid

  ! ---------------------------------  Modify_Values_For_Supersat  -----

  subroutine Modify_Values_For_Supersat ( fwdModelConf, &
    & grids_f, h2o_ind, grids_tmp, boundaryPressure )

    use Intrinsic, only: l_clear_110RH_below_top, l_clear_0RH
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSNumerics, only: HUNT
    use RHIFromH2O, only: RHIFromH2O_Factor
    use ForwardModelConfig, only: ForwardModelConfig_t
    use VectorsModule, only: VECTORVALUE_T
    use MLSNumerics, only: ESSENTIALLYEQUAL

    type (ForwardModelConfig_T) ,intent(in) :: FWDMODELCONF
    type (Grids_T), intent(inout) :: GRIDS_F   ! All the vmrs
    integer, intent(in) :: H2O_IND             ! Where is H2O in Grids_f?
    type (Grids_T), intent(inout) :: GRIDS_TMP ! All the temperatures
    type (VectorValue_T), intent(in) :: BOUNDARYPRESSURE

    integer :: KF, KZ, KP
    real(r8) :: RHI 
    integer :: Supersat_Index
    integer :: WF1, WF2
    integer :: V0               ! One before starting point in Values array
    integer :: Z_index
    integer :: Z0, P, P0
    logical :: FAILED

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
    failed = size ( grids_tmp%phi_basis, 1 ) /= boundaryPressure%template%noInstances
    if ( .not. failed ) then
      failed = .not. all ( EssentiallyEqual ( &
        & boundaryPressure%template%phi(1,wf1:wf2), grids_tmp%phi_basis ) )
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

    use Allocate_Deallocate, only: Allocate_test

    type (Grids_T), intent(inout) :: Grids_X
    integer, intent(in) :: N

    call allocate_test ( Grids_x%l_z, n, 'Grids_x%l_z', ModuleName, lowBound=0 )
    call allocate_test ( Grids_x%l_p, n, 'Grids_x%l_p', ModuleName, lowBound=0 )
    call allocate_test ( Grids_x%l_f, n, 'Grids_x%l_f', ModuleName, lowBound=0 )
    call allocate_test ( Grids_x%l_v, n, 'Grids_x%l_v', ModuleName, lowBound=0 )
    Grids_x%l_z(0) = 0
    Grids_x%l_p(0) = 0
    Grids_x%l_f(0) = 0
    Grids_x%l_v(0) = 0
    grids_x%p_len = 0
    call allocate_test ( Grids_x%windowstart, n, 'Grids_x%windowstart', &
                       & ModuleName )
    call allocate_test ( Grids_x%windowfinish, n, 'Grids_x%windowfinish',&
                       & ModuleName )
    call allocate_test ( Grids_x%lin_log, n, 'lin_log', ModuleName )
    call allocate_test ( Grids_x%min_val, n, 'min_val', ModuleName )

  end subroutine Create_Grids_1

  ! ---------------------------------------------  Create_Grids_2  -----
  subroutine Create_Grids_2 ( Grids_x )
  ! Create the part of the Grids_T structure that depends on the
  ! number of values and lengths of bases

    use Allocate_Deallocate, only: Allocate_test

    type (Grids_T), intent(inout) :: Grids_X

    integer :: N
    n = ubound(grids_x%l_z,1)

    call allocate_test ( Grids_x%zet_basis, Grids_x%l_z(n), 'Grids_x%zet_basis', &
                       & ModuleName )
    call allocate_test ( Grids_x%phi_basis, Grids_x%l_p(n), 'Grids_x%phi_basis', &
                       & ModuleName )
    call allocate_test ( Grids_x%frq_basis, Grids_x%l_f(n), 'Grids_x%frq_basis', &
                       & ModuleName )
    call allocate_test ( Grids_x%values, Grids_x%l_v(n), 'Grids_x%values', &
                       & ModuleName )
    call allocate_test ( Grids_x%deriv_flags, Grids_x%l_v(n), 'Grids_x%deriv_flags', &
                       & ModuleName )
  end subroutine Create_Grids_2

  ! -----------------------------------------------  Fill_Grids_1  -----
  subroutine Fill_Grids_1 ( Grids_x, II, Qty, Phitan, Maf, FwdModelConf )
  ! Fill in the size information for the II'th element of Grids_x

    use ForwardModelConfig, only: ForwardModelConfig_t
    use ManipulateVectorQuantities, only: FindInstanceWindow
    use VectorsModule, only: VectorValue_T

    type(grids_T), intent(inout) :: Grids_x
    integer, intent(in) :: II
    type (vectorValue_T), intent(in) :: QTY     ! An arbitrary vector quantity
    type (vectorValue_T), intent(in) :: PHITAN  ! Tangent geodAngle component of
    integer, intent(in) :: MAF
    type(forwardModelConfig_T), intent(in) :: FwdModelConf

    integer :: KF, KP, KZ

    kf = qty%template%noChans ! == 1 if qty%template%frequencyCoordinate == l_none
    call FindInstanceWindow ( qty, phitan, maf, fwdModelConf%phiWindow, &
      & fwdModelConf%windowUnits, grids_x%windowStart(ii), grids_x%windowFinish(ii) )
    kp = grids_x%windowFinish(ii) - grids_x%windowStart(ii) + 1
    kz = qty%template%noSurfs
    grids_x%l_f(ii) = grids_x%l_f(ii-1) + kf
    grids_x%l_p(ii) = grids_x%l_p(ii-1) + kp
    grids_x%l_v(ii) = grids_x%l_v(ii-1) + kz * kp * kf
    grids_x%l_z(ii) = grids_x%l_z(ii-1) + kz
    grids_x%p_len = grids_x%p_len + kz * kp

  end subroutine Fill_Grids_1

  ! -----------------------------------------------  Fill_Grids_2  -----
  subroutine Fill_Grids_2 ( Grids_x, II, Qty, SetDerivFlags )
  ! Fill the zeta, phi, freq. basis and value components for the II'th
  ! "molecule" in Grids_x.

    use Intrinsic, only: L_IntermediateFrequency
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Units, only: Deg2Rad
    use VectorsModule, only: VectorValue_T, M_FullDerivatives

    type(grids_t), intent(inout) :: Grids_x
    integer, intent(in) :: II
    type(vectorValue_T), intent(in) :: QTY     ! An arbitrary vector quantity
    logical, intent(in) :: SetDerivFlags

    integer :: KF, KZ, PF, PP, PV, PZ, QF, QP, QV, QZ, WF1, WF2

    pf = Grids_x%l_f(ii-1)
    pp = Grids_x%l_p(ii-1)
    pv = Grids_x%l_v(ii-1)
    pz = Grids_x%l_z(ii-1)

    qf = Grids_x%l_f(ii)
    qp = Grids_x%l_p(ii)
    qv = Grids_x%l_v(ii)
    qz = Grids_x%l_z(ii)

    kz = qz - pz
    kf = qf - pf
    wf1 = Grids_x%windowStart(ii)
    wf2 = Grids_x%windowFinish(ii)

    Grids_x%zet_basis(pz+1:qz) = qty%template%surfs(1:kz,1)
    Grids_x%phi_basis(pp+1:qp) = qty%template%phi(1,wf1:wf2)*Deg2Rad
    if ( associated ( qty%template%frequencies ) ) then
      if ( qty%template%frequencyCoordinate /= l_intermediateFrequency ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unexpected frequency coordinate for quantity' )
      Grids_x%frq_basis(pf+1:qf) = qty%template%frequencies
    else
      Grids_x%frq_basis(pf+1:qf) = 0.0
    end if

    ! qv - pv == (wf2 - wf1 + 1) * kf * kz here
    Grids_x%values(pv+1:qv) = reshape(qty%values(1:kf*kz,wf1:wf2), &
                                    & (/qv-pv/))
    !  mixing ratio values are manipulated if it's log basis
    Grids_x%lin_log(ii) = qty%template%logBasis
    if ( qty%template%logBasis ) then
      Grids_x%min_val(ii) = qty%template%minValue
      where ( Grids_x%values(pv+1:qv) <= grids_x%min_val(ii) ) &
           & Grids_x%values(pv+1:qv) = grids_x%min_val(ii)
      Grids_x%values(pv+1:qv) = log(Grids_x%values(pv+1:qv))
    end if

    ! set 'do derivative' flags

    if ( setDerivFlags ) then
      if ( associated(qty%mask) ) then
        Grids_x%deriv_flags(pv+1:qv) = reshape(( iand (M_FullDerivatives, &
          & ichar(qty%mask)) == 0),(/qv-pv/))
      else
        Grids_x%deriv_flags(pv+1:qv) = .true.
      end if
    else
      grids_x%deriv_flags(pv+1:qv) = .false.
    end if

  end subroutine Fill_Grids_2

  ! ---------------------------------------------  DestroyGrids_t  -----
  subroutine DestroyGrids_t ( grids_x )
    use Allocate_Deallocate, only: Deallocate_test

    type (Grids_T), intent(inout) :: Grids_x

    call deallocate_test(grids_x%l_f,'grids_x%l_f',modulename)
    call deallocate_test(grids_x%l_z,'grids_x%l_z',modulename)
    call deallocate_test(grids_x%l_p,'grids_x%l_p',modulename)
    call deallocate_test(grids_x%l_v,'grids_x%l_v',modulename)
    call deallocate_test(grids_x%values,'grids_x%values',modulename)
    call deallocate_test(grids_x%lin_log,'grids_x%lin_log',modulename)
    call deallocate_test(grids_x%min_val,'grids_x%min_val',modulename)
    call deallocate_test(grids_x%frq_basis,'grids_x%frq_basis',modulename)
    call deallocate_test(grids_x%zet_basis,'grids_x%zet_basis',modulename)
    call deallocate_test(grids_x%phi_basis,'grids_x%phi_basis',modulename)
    call deallocate_test(grids_x%deriv_flags,'grids_x%deriv_flags',modulename)
    call deallocate_test(grids_x%windowstart,'grids_x%windowstart',modulename)
    call deallocate_test(grids_x%windowfinish,'grids_x%windowfinish',modulename)

  end subroutine Destroygrids_t

  ! -------------------------------------------------  Dump_Grids  -----
  subroutine Dump_Grids ( The_Grid, Name )
  ! Dump The_Grid

    use Dump_0, only: Dump
    use Output_M, only: Output

    type(grids_t), intent(in) :: The_Grid
    character(len=*), intent(in), optional :: Name

    call output ( 'Dump of Grids_T structure ', advance='no' )
    if ( present(name) ) call output ( name, advance='no' )
    call output ( '', advance='yes' )
    call dump ( the_grid%windowStart, 'The_grid%WindowStart' )
    call dump ( the_grid%windowFinish, 'The_grid%WindowFinish' )
    call dump ( the_grid%lin_log, 'The_grid%Lin_Log' )
    call dump ( the_grid%min_val, 'The_grid%Min_Val' )
    call dump ( the_grid%frq_basis, 'The_grid%Frq_Basis' )
    call dump ( the_grid%zet_basis, 'The_grid%Zet_Basis' )
    call dump ( the_grid%phi_basis, 'The_grid%Phi_Basis' )
    call dump ( the_grid%values, 'The_grid%Values' )
    call dump ( the_grid%deriv_flags, 'The_grid%Deriv_Flags' )

  end subroutine Dump_Grids

  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module LOAD_SPS_DATA_M
! $Log$
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
