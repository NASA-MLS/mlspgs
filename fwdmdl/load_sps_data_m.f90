! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module LOAD_SPS_DATA_M
  use MLSCommon, only: R8, RP, IP
  use Units, only: Deg2Rad
  use ForwardModelConfig, only: FORWARDMODELCONFIG_T
  use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
  USE INTRINSIC, ONLY: L_VMR, L_FREQUENCY, L_NONE, L_PHITAN
! USE INTRINSIC, ONLY: L_VMR, L_FREQUENCY, L_NONE, L_PHITAN, L_DW, L_DN, L_DV
  use VectorsModule, only: Vector_T, VectorValue_T, GetVectorQuantityByType, &
                        &  M_FullDerivatives
  use Molecules, only: spec_tags, L_EXTINCTION
  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error

  use SpectroscopyCatalog_m, only: CATALOG_T
  USE manipulatevectorquantities, only: FindInstanceWindow

  implicit none

  Private
  PUBLIC :: load_sps_data, Destroygrids_t

  type, public :: Grids_T             ! Fit all Gridding categories
    integer,  pointer :: no_f(:) => null()! No. of entries in frq. grid per sps
    integer,  pointer :: no_z(:) => null()! No. of entries in zeta grid per sps
    integer,  pointer :: no_p(:) => null()! No. of entries in phi  grid per sps
    INTEGER,  pointer :: windowstart(:) => null()! horizontal starting index
!                                                  from l2gp
    INTEGER,  pointer :: windowfinish(:) => null()! horizontal ending index
!                                                   from l2gp
    LOGICAL,  pointer :: lin_log(:) => null()   ! type of representation basis
    real(r8), pointer :: frq_basis(:) => null() ! frq grid entries for all
!                                                 molecules
    real(rp), pointer :: zet_basis(:) => null() ! zeta grid entries for all
!                                                 molecules
    real(rp), pointer :: phi_basis(:) => null() ! phi  grid entries for all
!                                                 molecules
    REAL(rp), pointer :: values(:) => null() ! species values (ie vmr) in lvf
    LOGICAL,  pointer :: deriv_flags(:) => null() ! do derivatives flags in lvf
  end type Grids_T

!---------------------------- RCS Ident Info -------------------------------
  character (LEN=256) :: Id = &
   "$Id$"
  character (LEN=*), parameter :: ModuleName= &
   "$RCSfile$"
!---------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------

 SUBROUTINE load_sps_data(FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
       &    radiometer, mol_cat_index, p_len, f_len, h2o_ind, ext_ind,     &
       &    Grids_f, f_len_dw, Grids_dw, f_len_dn, Grids_dn, f_len_dv,     &
       &    Grids_dv)

    Type(forwardModelConfig_T), intent(in) :: fwdModelConf
    Type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    Type(forwardModelStatus_t), intent(in) :: FmStat ! Reverse comm. stuff

    Integer, intent(in)  :: RADIOMETER
    Integer, intent(in)  :: MOL_CAT_INDEX(:)

    Integer, intent(out) :: P_LEN, F_LEN
    Integer, intent(out) :: H2O_IND
    Integer, intent(out) :: EXT_IND

    Integer, optional, intent(out) :: F_LEN_DW, F_LEN_DN, F_LEN_DV

! All the VMR coordinates
    Type (Grids_T), intent(out) :: Grids_f

! All the spectroscopy(W) coordinates
    Type (Grids_T), optional, intent(out) :: Grids_dw

! All the spectroscopy(N) coordinates
    Type (Grids_T), optional, intent(out) :: Grids_dn

! All the spectroscopy(V) coordinates
    Type (Grids_T), optional, intent(out) :: Grids_dv
!
! Begin code:
!
    Call load_one_grid(FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
       & radiometer, mol_cat_index, f_len, 'f', Grids_f, p_len, h2o_ind,&
       & ext_ind )
!
! ** When the spectroscopy flags are properly introduced into the database,
!    un-comment the following codes:
!
!   if( PRESENT ( Grids_dw ) ) &
!   & Call load_one_grid(FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
!        & radiometer, mol_cat_index, f_len_dw, 'w', Grids_dw)
!
!   if( PRESENT ( Grids_dn ) ) &
!   & Call load_one_grid(FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
!        & radiometer, mol_cat_index, f_len_dn, 'n', Grids_dn)
!
!   if( PRESENT ( Grids_dv ) ) &
!   & Call load_one_grid(FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
!        & radiometer, mol_cat_index, f_len_dv, 'v', Grids_dv)
!
 END SUBROUTINE load_sps_data
!-------------------------------------------------------------------
!
 SUBROUTINE load_one_grid(FwdModelConf, fwdModelIn, fwdModelExtra, FmStat, &
       &    radiometer, mol_cat_index, f_len, Grid_type, Grids_x, p_len,   &
       &    h2o_ind, ext_ind)

    type(forwardModelConfig_T), intent(in) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(forwardModelStatus_t), intent(in) :: FmStat ! Reverse comm. stuff

    integer, intent(in)  :: RADIOMETER
    integer, intent(in)  :: MOL_CAT_INDEX(:)

    integer, intent(out) :: F_LEN

    integer, optional, intent(out) :: P_LEN
    integer, optional, intent(out) :: H2O_IND
    integer, optional, intent(out) :: EXT_IND

    type (Grids_T), intent(out) :: Grids_x   ! All the coordinates

    Character(LEN=1), intent(in) :: Grid_type
!
!** ZEBUG: When Intrinsic has (l_dw, l_dn & l_dv), get rid of the following 
!   4 lines of code
!
    Integer, parameter :: l_dw = 1
    Integer, parameter :: l_dn = 2
    Integer, parameter :: l_dv = 3

! Local variables:

    Integer :: i,j,k,l,m,n,r,s,kz,kp,kf,n_f_phi,n_f_zet,n_f_frq,no_mol,l_x, &
            &  ii, kk, wf1, wf2

    type (VectorValue_T), pointer :: F             ! An arbitrary species
    type (VectorValue_T), pointer :: PHITAN ! Tangent geodAngle component of

!
    !******************* LOAD SPECIES DATA ************

    no_mol = size( mol_cat_index )

    Call allocate_test ( Grids_x%no_z,no_mol,'Grids_x%no_z',ModuleName )
    Call allocate_test ( Grids_x%no_p,no_mol,'Grids_x%no_p',ModuleName )
    Call allocate_test ( Grids_x%no_f,no_mol,'Grids_x%no_f',ModuleName )
    Call allocate_test ( Grids_x%windowstart,no_mol,'Grids_x%windowstart', &
                       & ModuleName )
    Call allocate_test ( Grids_x%windowfinish,no_mol,'Grids_x%windowfinish',&
                       & ModuleName )
    Call Allocate_test ( Grids_x%lin_log, no_mol, 'lin_log', ModuleName )

    Grids_x%no_z = 0
    Grids_x%no_p = 0
    Grids_x%no_f = 0

    f_len = 0

    if( PRESENT(p_len) ) p_len=0
    if( PRESENT(ext_ind) ) ext_ind = 0
    if( PRESENT(h2o_ind) ) h2o_ind = 0

    phitan => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_phitan, &
      & instrumentModule=FwdModelConf%signals(1)%instrumentModule )

    l_x = l_vmr
    if(Grid_type == 'W' .OR. Grid_type == 'w') l_x = l_dw
    if(Grid_type == 'N' .OR. Grid_type == 'n') l_x = l_dn
    if(Grid_type == 'V' .OR. Grid_type == 'v') l_x = l_dv

    do ii = 1, no_mol
      kk = FwdModelConf%molecules(mol_cat_index(ii))
      if(PRESENT(h2o_ind) .AND. spec_tags(kk) == 18003) h2o_ind = ii
      if ( kk == l_extinction ) then
        if( PRESENT(ext_ind) ) ext_ind = ii
        f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_extinction, radiometer=radiometer)
      else
        f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_x, molecule=kk)
      endif
      kz = f%template%noSurfs
      if ( f%template%frequencyCoordinate == l_none ) then
        kf = 1
      else
        kf = f%template%noChans
      endif
      Call FindInstanceWindow(f,phitan,fmStat%maf,FwdModelConf%phiWindow, &
                            & wf1, wf2)
      Grids_x%windowStart(ii) = wf1
      Grids_x%windowFinish(ii) = wf2
      kp = wf2 - wf1 + 1
      Grids_x%no_f(ii) = kf
      Grids_x%no_z(ii) = kz
      Grids_x%no_p(ii) = kp
      if( PRESENT(p_len) ) p_len = p_len + kz * kp
      f_len = f_len + kz * kp * kf
      if (f%template%logBasis) then
        Grids_x%lin_log(ii) = .TRUE.
      else
        Grids_x%lin_log(ii) = .FALSE.
      endif
   end do
!
    n_f_zet = SUM(Grids_x%no_z)
    n_f_phi = SUM(Grids_x%no_p)
    n_f_frq = SUM(Grids_x%no_f)
!
! Allocate space for the zeta, phi & freq. basis componenets
!
    Call allocate_test ( Grids_x%zet_basis,n_f_zet,'Grids_x%zet_basis', &
                       & ModuleName)
    Call allocate_test ( Grids_x%phi_basis,n_f_phi,'Grids_x%phi_basis', &
                       & ModuleName)
    Call allocate_test ( Grids_x%frq_basis,n_f_frq,'Grids_x%frq_basis', &
                       & ModuleName)
    Call allocate_test ( Grids_x%values,f_len,'Grids_x%values', &
                       & ModuleName)
    Call allocate_test ( Grids_x%deriv_flags,f_len,'Grids_x%deriv_flags',&
                       & ModuleName)
!
    j = 1
    l = 1
    s = 1
    f_len = 1
    do ii = 1, no_mol
      i = mol_cat_index(ii)
      kk = FwdModelConf%molecules(i)
      if ( kk == l_extinction ) then
        f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_extinction, radiometer=radiometer )
      else
        f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_x, molecule=kk )
      endif
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
        Grids_x%frq_basis(s:m-1) = f%template%frequencies
      else
        Grids_x%frq_basis(s:m-1) = 0.0
      endif
!
! ** ZEBUG - Simulate f%values for EXTINCTION, using the N2 function
!  (Some code here ...)
! ** END ZEBUG
!
      r = f_len + kf * kz * kp
      Grids_x%values(f_len:r-1) = RESHAPE(f%values(1:kf*kz,wf1:wf2), &
                                      & (/kf*kz*kp/))
      if (Grids_x%lin_log(ii)) then
        WHERE (Grids_x%values(f_len:r-1) <= 1.0e-16_rp) &
             & Grids_x%values(f_len:r-1) = 1.0e-16_rp
        Grids_x%values(f_len:r-1) = LOG(Grids_x%values(f_len:r-1))
      endif
!
! set 'do derivative' flags
!
      IF (ASSOCIATED(f%mask)) THEN
        Grids_x%deriv_flags(f_len:r-1) = RESHAPE(( iand (M_FullDerivatives,&
                          & ICHAR(f%mask)) == 0),(/kf*kz*kp/))
      ELSE
        Grids_x%deriv_flags(f_len:r-1) = .TRUE.
      ENDIF
!
      j = k
      l = n
      s = m
      f_len = r
!
    end do
!
    f_len = f_len - 1
!
 END subroutine load_one_grid
!
!----------------------------------------------------------------
 Subroutine DestroyGrids_t( grids_x )
!
  TYPE(Grids_T), intent(inout) :: Grids_x
!
  Call deallocate_test(grids_x%no_f,'grids_x%no_f',modulename)
  Call deallocate_test(grids_x%no_z,'grids_x%no_z',modulename)
  Call deallocate_test(grids_x%no_p,'grids_x%no_p',modulename)
  Call deallocate_test(grids_x%values,'grids_x%values',modulename)
  Call deallocate_test(grids_x%lin_log,'grids_x%lin_log',modulename)
  Call deallocate_test(grids_x%frq_basis,'grids_x%frq_basis',modulename)
  Call deallocate_test(grids_x%zet_basis,'grids_x%zet_basis',modulename)
  Call deallocate_test(grids_x%phi_basis,'grids_x%phi_basis',modulename)
  Call deallocate_test(grids_x%deriv_flags,'grids_x%deriv_flags',modulename)
  Call deallocate_test(grids_x%windowstart,'grids_x%windowstart',modulename)
  Call deallocate_test(grids_x%windowfinish,'grids_x%windowfinish',modulename)

 End subroutine Destroygrids_t

end module LOAD_SPS_DATA_M
! $Log$
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
