! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module LOAD_SPS_DATA_M

  use MLSCommon, only: R8, RP

  implicit NONE

  private
  public :: Load_SPS_Data, Load_One_Grid, Destroygrids_t

  type, public :: Grids_T             ! Fit all Gridding categories
    integer,  pointer :: no_f(:) => null()! No. of entries in frq. grid per sps
    integer,  pointer :: no_z(:) => null()! No. of entries in zeta grid per sps
    integer,  pointer :: no_p(:) => null()! No. of entries in phi  grid per sps
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
    real(rp), pointer :: values(:) => null() ! species values (ie vmr) in lvf
    logical,  pointer :: deriv_flags(:) => null() ! do derivatives flags in lvf
  end type Grids_T

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

  subroutine Load_SPS_Data ( FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
       &    radiometer, mol_cat_index, p_len, f_len, h2o_ind, ext_ind,        &
       &    Grids_f, f_len_dw, Grids_dw, f_len_dn, Grids_dn, f_len_dv,        &
       &    Grids_dv, i_supersat, temp_supersat )

    use ForwardModelConfig, only: ForwardModelConfig_t
    use ForwardModelIntermediate, only: ForwardModelStatus_t
    use Intrinsic, only: L_DN, L_DV, L_DW, L_VMR
    use VectorsModule, only: Vector_T

    type(forwardModelConfig_T), intent(in) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(forwardModelStatus_t), intent(in) :: FmStat ! Reverse comm. stuff

    integer, intent(in)  :: RADIOMETER
    integer, intent(in)  :: MOL_CAT_INDEX(:)

    integer, intent(out) :: P_LEN, F_LEN
    integer, intent(out) :: H2O_IND
    integer, intent(out) :: EXT_IND

    integer, optional, intent(out) :: F_LEN_DW, F_LEN_DN, F_LEN_DV

! All the VMR coordinates
    type (Grids_T), intent(out) :: Grids_f

! All the spectroscopy(W) coordinates
    type (Grids_T), optional, intent(out) :: Grids_dw

! All the spectroscopy(N) coordinates
    type (Grids_T), optional, intent(out) :: Grids_dn

! All the spectroscopy(V) coordinates
    type (Grids_T), optional, intent(out) :: Grids_dv

    integer, optional, intent(in)  :: i_supersat     ! Do the suprsaturation calculation?
!-----------------------------------------------------------------------------
! i_supersat indicates different clear and cloudy sky combinations:
!        i_supersat =-1 is for clear-sky radiance limit assuming 110%RHi
!        i_supersat =-2 is for clear-sky radiance limit assuming 0%RHi
!-----------------------------------------------------------------------------
    real(r8), dimension(:), optional, intent(in) :: &
                                    & temp_supersat  ! What temperatures to use for supersaturation
! Begin code:

    call load_one_grid ( FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
       & radiometer, mol_cat_index, f_len, l_vmr, Grids_f, p_len, h2o_ind,&
       & ext_ind, i_supersat, temp_supersat )

! ** When the spectroscopy flags are properly introduced into the database,
!    un-comment the following codes:

!   if ( PRESENT ( Grids_dw ) ) &
!   & Call load_one_grid(FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
!        & radiometer, mol_cat_index, f_len_dw, l_dw, Grids_dw)

!   if ( PRESENT ( Grids_dn ) ) &
!   & Call load_one_grid(FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
!        & radiometer, mol_cat_index, f_len_dn, l_dn, Grids_dn)

!   if ( PRESENT ( Grids_dv ) ) &
!   & Call load_one_grid(FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
!        & radiometer, mol_cat_index, f_len_dv, l_dv, Grids_dv)

  end subroutine Load_Sps_Data

  !  ---------------------------------------------  Load_One_Grid  -----

  subroutine Load_One_Grid ( FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
       &    radiometer, mol_cat_index, f_len, QuantityType, Grids_x, p_len,      &
       &    h2o_ind, ext_ind, i_supersat, temp_supersat )

    use Allocate_Deallocate, only: Allocate_test
    use ForwardModelConfig, only: ForwardModelConfig_t
    use ForwardModelIntermediate, only: ForwardModelStatus_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: L_IntermediateFrequency, L_None, L_Phitan
    use ManipulateVectorQuantities, only: FindInstanceWindow
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSNumerics, only: HUNT
    use Molecules, only: Spec_tags, SP_H2O
    use RHIFromH2O, only: RHIFromH2O_Factor
    use Units, only: Deg2Rad
    use VectorsModule, only: Vector_T, VectorValue_T, M_FullDerivatives

    type(forwardModelConfig_T), intent(in) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(forwardModelStatus_t), intent(in) :: FmStat ! Reverse comm. stuff

    integer, intent(in) :: RADIOMETER       
    integer, intent(in) :: MOL_CAT_INDEX(:) 

    integer, intent(out) :: F_LEN

    integer, intent(in) :: QuantityType      ! L_vmr, L_dv, etc.

    type (Grids_T), intent(out) :: Grids_x   ! All the coordinates

    integer, optional, intent(out) :: P_LEN
    integer, optional, intent(out) :: H2O_IND
    integer, optional, intent(out) :: EXT_IND
    integer, optional, intent(in)  :: i_supersat     ! Do the suprsaturation calculation?
!-----------------------------------------------------------------------------
! i_supersat indicates different clear and cloudy sky combinations:
!        i_supersat =-1 is for clear-sky radiance limit assuming 110%RHi
!        i_supersat =-2 is for clear-sky radiance limit assuming 0%RHi
!-----------------------------------------------------------------------------
    real(r8), dimension(:), optional, intent(in) :: &
                                    & temp_supersat  ! What temperatures to use for supersaturation
! Local variables:

    integer :: i, j, k, l, m, n, r, s, kz, kp, kf
    integer :: n_f_phi, n_f_zet, n_f_frq, no_mol
    integer :: ii, kk, wf1, wf2

    type (VectorValue_T), pointer :: F       ! An arbitrary species
    type (VectorValue_T), pointer :: PHITAN  ! Tangent geodAngle component of

    logical :: FoundInFirst                  ! Flag
    integer :: supersat_Index
    integer :: wf
    integer :: f_index
    integer :: z_index
    integer :: values_index
    real(r8), parameter :: refPRESSURE = 100._r8
    real(r8) :: RHI 


    !******************* LOAD SPECIES DATA ************

    no_mol = size( mol_cat_index )

    call allocate_test ( Grids_x%no_z,no_mol,'Grids_x%no_z',ModuleName )
    call allocate_test ( Grids_x%no_p,no_mol,'Grids_x%no_p',ModuleName )
    call allocate_test ( Grids_x%no_f,no_mol,'Grids_x%no_f',ModuleName )
    call allocate_test ( Grids_x%windowstart,no_mol,'Grids_x%windowstart', &
                       & ModuleName )
    call allocate_test ( Grids_x%windowfinish,no_mol,'Grids_x%windowfinish',&
                       & ModuleName )
    call Allocate_test ( Grids_x%lin_log, no_mol, 'lin_log', ModuleName )
    call Allocate_test ( Grids_x%min_val, no_mol, 'min_val', ModuleName )

    Grids_x%no_z = 0
    Grids_x%no_p = 0
    Grids_x%no_f = 0
    grids_x%min_val = -huge(0.0_r8)

    f_len = 0

    if ( present(p_len) ) p_len=0
    if ( present(ext_ind) ) ext_ind = 0
    if ( present(h2o_ind) ) h2o_ind = 0

    phitan => GetQuantityforForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_phitan, config=fwdModelConf, &
      & instrumentModule=FwdModelConf%signals(1)%instrumentModule )

    do ii = 1, no_mol
      kk = FwdModelConf%molecules(mol_cat_index(ii))
      if ( present(h2o_ind) .and. spec_tags(kk) == sp_h2o ) h2o_ind = ii
      f => GetQuantityforForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=quantityType, molIndex=mol_cat_index(ii), &
        & radiometer=radiometer, config=fwdModelConf )
      kz = f%template%noSurfs
      if ( f%template%frequencyCoordinate == l_none ) then
        kf = 1
      else
        kf = f%template%noChans
      end if
      call FindInstanceWindow ( f, phitan, fmStat%maf, FwdModelConf%phiWindow, &
        & FwdModelConf%windowUnits, wf1, wf2 )
      Grids_x%windowStart(ii) = wf1
      Grids_x%windowFinish(ii) = wf2
      kp = wf2 - wf1 + 1
      Grids_x%no_f(ii) = kf
      Grids_x%no_z(ii) = kz
      Grids_x%no_p(ii) = kp
      if ( present(p_len) ) p_len = p_len + kz * kp
      f_len = f_len + kz * kp * kf
      if ( f%template%logBasis ) then
        Grids_x%lin_log(ii) = .true.
        Grids_x%min_val(ii) = f%template%minValue
      else
        Grids_x%lin_log(ii) = .false.
      end if
    end do

    n_f_zet = sum(Grids_x%no_z)
    n_f_phi = sum(Grids_x%no_p)
    n_f_frq = sum(Grids_x%no_f)

! Allocate space for the zeta, phi & freq. basis componenets

    call allocate_test ( Grids_x%zet_basis,n_f_zet,'Grids_x%zet_basis', &
                       & ModuleName )
    call allocate_test ( Grids_x%phi_basis,n_f_phi,'Grids_x%phi_basis', &
                       & ModuleName )
    call allocate_test ( Grids_x%frq_basis,n_f_frq,'Grids_x%frq_basis', &
                       & ModuleName )
    call allocate_test ( Grids_x%values,f_len,'Grids_x%values', &
                       & ModuleName )
    call allocate_test ( Grids_x%deriv_flags,f_len,'Grids_x%deriv_flags',&
                       & ModuleName )

    j = 1
    l = 1
    s = 1
    f_len = 1
    do ii = 1, no_mol
      i = mol_cat_index(ii)
      f => GetQuantityforForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=quantityType, molIndex=i, radiometer=radiometer, &
        & config=fwdModelConf, foundInFirst=foundInFirst )
      kz = Grids_x%no_z(ii)
      kp = Grids_x%no_p(ii)
      kf = Grids_x%no_f(ii)
      n = l + kz
      m = s + kf
      k = j + kp
      wf1 = Grids_x%windowStart(ii)
      wf2 = Grids_x%windowFinish(ii)
      Grids_x%zet_basis(l:n-1) = f%template%surfs(1:kz,1)
      Grids_x%phi_basis(j:k-1) = f%template%phi(1,wf1:wf2)*Deg2Rad
      if ( associated ( f%template%frequencies ) ) then
        if ( f%template%frequencyCoordinate /= l_intermediateFrequency ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unexpected frequency coordinate for quantity' )
        Grids_x%frq_basis(s:m-1) = f%template%frequencies
      else
        Grids_x%frq_basis(s:m-1) = 0.0
      end if

! ** ZEBUG - Simulate f%values for EXTINCTION, using the N2 function
!  (Some code here ...)
! ** END ZEBUG

      r = f_len + kf * kz * kp
      Grids_x%values(f_len:r-1) = reshape(f%values(1:kf*kz,wf1:wf2), &
                                      & (/kf*kz*kp/))
      if ( Grids_x%lin_log(ii) ) then
        where ( Grids_x%values(f_len:r-1) <= grids_x%min_val(ii) ) &
             & Grids_x%values(f_len:r-1) = grids_x%min_val(ii)
        Grids_x%values(f_len:r-1) = log(Grids_x%values(f_len:r-1))
      end if

! set 'do derivative' flags

      if ( fwdModelConf%moleculeDerivatives(i) .and. foundInFirst ) then
        if ( associated(f%mask) ) then
          Grids_x%deriv_flags(f_len:r-1) = reshape(( iand (M_FullDerivatives,&
            & ichar(f%mask)) == 0),(/kf*kz*kp/))
        else
          Grids_x%deriv_flags(f_len:r-1) = .true.
        end if
      else
        grids_x%deriv_flags(f_len:r-1) = .false.
      end if

      j = k
      l = n
      s = m
      f_len = r

      ! do the following only if i_supersat not equal to 0
!??? The next line is almost surely wrong.  It is inconsistent with the
!??? above comment, and it leads to using RHI without it having a value
!??? about twenty lines down from here.
!     if ( i_supersat .eq. 0 ) then
      if ( i_supersat /= 0 ) then
        select case ( i_supersat )
        case ( -1 )
          RHI=110.0_r8
        case ( -2 )
          RHI=1.0e-9_r8
        end select  

        call Hunt (Grids_x%zet_basis(1:n-1), -log(refPRESSURE), supersat_Index, &
        & l, nearest=.true.)

        if (supersat_Index < 1) then        
           call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Could not find value of Zeta in Grid_T; returned index too small' )

        else if (supersat_Index > size(Grids_x%values) ) then
           call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Could not find value of Zeta in Grid_T; returned index too BIG' )
        else
          do wf = wf1, wf2
            do f_index = 1, kf
              do z_index = 1, supersat_Index
                values_index = f_len - 1 + f_index + kf*(z_index-1) + kf*kz*(wf-wf1)
                Grids_x%values(values_index) = RHIFromH2O_Factor (temp_supersat(z_index), &
                 & grids_x%zet_basis(l+z_index-1), 0, .true.)*RHI
              end do
            end do
          end do
        end if
      end if

    end do

    f_len = f_len - 1

  end subroutine Load_One_Grid

!----------------------------------------------------------------
  subroutine DestroyGrids_t( grids_x )
    use Allocate_Deallocate, only: Deallocate_test

    type(Grids_T), intent(inout) :: Grids_x

    call deallocate_test(grids_x%no_f,'grids_x%no_f',modulename)
    call deallocate_test(grids_x%no_z,'grids_x%no_z',modulename)
    call deallocate_test(grids_x%no_p,'grids_x%no_p',modulename)
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

  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module LOAD_SPS_DATA_M
! $Log$
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
