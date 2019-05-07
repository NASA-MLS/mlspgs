! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_Species_Data_m

  ! Get species data for the molecules in the beta groups.
  ! Get the spectral parameters from the state vector.
  ! Get the temperature quantity.

  implicit NONE
  private
  public :: Get_Species_Data

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Get_Species_Data ( FwdModelConf, FwdModelIn, FwdModelExtra )

    ! Fill in the Beta_Groups' isotope ratios in FwdModelConf.
    ! Get vector quantities.

    use ForwardModelConfig, only: Dump, ForwardModelConfig_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel, &
      GetQtyStuffForForwardModel
    use Intrinsic, only: L_IsotopeRatio, L_LineCenter, L_LineWidth, &
      & L_LineWidth_TDep, L_Temperature, L_VMR
    use MLSStringLists, only: SwitchDetail
    use SpectroscopyCatalog_m, only: Catalog, Dump
    use Toggles, only: Emit, Levels, Switches, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use VectorsModule, only: GetVectorQuantityByType, Vector_t, VectorValue_t

    type(forwardModelConfig_t), intent(inout) :: FwdModelConf ! Fills Beta_Group
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra

    integer :: B         ! Index for beta groups 
    integer, save :: DumpFWM = -2          ! -2 = not called yet, -1 = no dumps,
                                           ! low-order digit: catalog dump level
                                           ! high-order digit: 1 => stop
    type (VectorValue_T), pointer :: F  ! An arbitrary species
    integer :: L         ! Index in spectral parameter data structure
    integer :: M         ! Index for molecules in beta groups, or size thereof
    integer :: Me = -1   ! String index cache for tracing
    integer :: Mol       ! Molecule in a beta group
    integer :: S         ! Sideband index, 1 = LSB, 2 = USB
    integer :: S1, S2    ! Bounds for sideband index S

    call trace_begin ( me, 'Get_Species_Data', &
      & cond=toggle(emit) .and. levels(emit) > 1 ) ! set by -f command-line switch

    if ( dumpFWM < -1 ) dumpFwm = switchDetail(switches,'fwmd')

    s1 = (fwdModelConf%sidebandStart+3)/2; s2 = (fwdModelConf%sidebandStop+3)/2

    if ( associated(fwdModelConf%beta_group) ) then
      ! Get isotope ratios for molecules in a beta group, else 1.0 if not a group
      do b = 1, size(fwdModelConf%beta_group)
        if ( fwdModelConf%beta_group(b)%group ) then ! A molecule group
          ! First LBL molecules' ratios
          do s = s1, s2
            do m = 1, size(fwdModelConf%beta_group(b)%lbl(s)%molecules)
              mol = fwdModelConf%beta_group(b)%lbl(s)%molecules(m)
              f => getQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
                & quantityType=l_isotopeRatio,                             &
                & molecule=mol, noError=.TRUE., config=fwdModelConf )
              fwdModelConf%beta_group(b)%lbl(s)%ratio(m) = &
                & catalog(mol)%defaultIsotopeRatio
              if ( associated ( f ) ) & ! Have an isotope ratio
                & fwdModelConf%beta_group(b)%lbl(s)%ratio(m) = f%values(1,1)
            end do ! m
            if ( associated(fwdModelConf%beta_group(b)%pfa(s)%molecules) ) then
              ! Now PFA molecules' ratios
              do m = 1, size(fwdModelConf%beta_group(b)%pfa(s)%molecules)
                mol = fwdModelConf%beta_group(b)%pfa(s)%molecules(m)
                f => getQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
                  & quantityType=l_isotopeRatio,                             &
                  & molecule=mol, noError=.TRUE., config=fwdModelConf )
                fwdModelConf%beta_group(b)%pfa(s)%ratio(m) = &
                  & catalog(mol)%defaultIsotopeRatio
                if ( associated ( f ) ) & ! Have an isotope ratio
                  & fwdModelConf%beta_group(b)%pfa(s)%ratio(m) = f%values(1,1)
              end do ! m
            end if
          end do ! s
  !     else ! Not a molecule group, but this is handled in ForwardModelSupport
  !       fwdModelConf%beta_group(b)%lbl_ratio(1)     = 1.0
  !       fwdModelConf%beta_group(b)%pfa_ratio(1)     = 1.0
        end if

      ! Get state vector quantity for species
        fwdModelConf%beta_group(b)%qty = GetQtyStuffForForwardModel ( &
          &  fwdModelIn, fwdModelExtra, quantityType=l_vmr, molIndex=b,    &
          &  config=fwdModelConf, radiometer=fwdModelConf%signals(1)%radiometer, &
          &  noError=.false. )
      end do ! b
    ! else
      ! No species, maybe doing a hydrostatic fill
    end if

    ! Get temperature vector quantity
    fwdModelConf%temp = GetQtyStuffForForwardModel ( &
        &  fwdModelIn, fwdModelExtra, quantityType=l_temperature, &
        &  config=fwdModelConf )

    ! Get state vector quantities for spectral parameters
    ! fwdModelConf%line..._ix have indices into these
    if ( associated(fwdModelConf%lineCenter) ) then
      do l = 1, size(fwdModelConf%lineCenter)
        fwdModelConf%lineCenter(l)%qty%qty => getVectorQuantityByType ( &
        &  fwdModelIn, fwdModelExtra, quantityType=l_lineCenter, &
        &  molecule=fwdModelConf%lineCenter(l)%molecule, &
        &  foundInFirst=fwdModelConf%lineCenter(l)%qty%foundInFirst )
        call check_no_frq_coord ( fwdModelConf%lineCenter(l)%qty%qty )
      end do ! l
    end if
    if ( associated(fwdModelConf%lineWidth) ) then
      do l = 1, size(fwdModelConf%lineWidth)
        fwdModelConf%lineWidth(l)%qty%qty => getVectorQuantityByType ( &
        &  fwdModelIn, fwdModelExtra, quantityType=l_lineWidth, &
        &  molecule=fwdModelConf%lineWidth(l)%molecule, &
        &  foundInFirst=fwdModelConf%lineWidth(l)%qty%foundInFirst )
        call check_no_frq_coord ( fwdModelConf%lineWidth(l)%qty%qty )
      end do
    end if
    if ( associated(fwdModelConf%lineWidth_TDep) ) then
      do l = 1, size(fwdModelConf%lineWidth_TDep)
        fwdModelConf%lineWidth_TDep(l)%qty%qty => getVectorQuantityByType ( &
        &  fwdModelIn, fwdModelExtra, quantityType=l_lineWidth_TDep, &
        &  molecule=fwdModelConf%lineWidth_TDep(l)%molecule, &
        &  foundInFirst=fwdModelConf%lineWidth_TDep(l)%qty%foundInFirst )
        call check_no_frq_coord ( fwdModelConf%lineWidth_TDep(l)%qty%qty )
      end do
    end if

    if ( dumpFWM > -1 ) then
      call dump ( fwdModelConf, moduleName(11:len_trim(moduleName)-8) )
      call dump ( fwdModelConf%catalog, details=mod(dumpFwm,10) )
      if ( dumpFWM > 9 ) stop
    end if

    call trace_end ( 'Get_Species_Data', &
      & cond=toggle(emit) .and. levels(emit) > 1 )

  contains

    subroutine Check_No_Frq_Coord ( Qty )
      use Intrinsic, only: L_None
      use MLSMessageModule, only: MLSMSG_Error
      use MoreMessage, only: MLSMessage
      type(vectorValue_t), intent(in) :: Qty
      if ( qty%template%frequencyCoordinate == l_none ) return
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'The %s vector quantity shall not have a frequency coordinate', &
        & (/ qty%template%name /) )
    end subroutine Check_No_Frq_Coord

  end subroutine Get_Species_Data

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Get_Species_Data_m

! $Log$
! Revision 2.40  2016/05/05 23:29:45  vsnyder
! Get temperature the same way it used to be gotten in FullForwardModel
!
! Revision 2.39  2016/05/02 23:33:53  vsnyder
! Get the temperature QtyStuff into FwdModelConf
!
! Revision 2.38  2014/04/22 00:37:18  vsnyder
! Add tracing
!
! Revision 2.37  2013/08/02 01:23:33  vsnyder
! Use GetQtyStuffForForwardModel
!
! Revision 2.36  2013/05/15 03:08:54  vsnyder
! Revise processing of dump switch
!
! Revision 2.35  2012/05/08 01:34:54  vsnyder
! Get default isotope ratio from catalog, not from 1.0
!
! Revision 2.34  2011/05/09 17:48:06  pwagner
! Converted to using switchDetail
!
! Revision 2.33  2011/03/31 19:52:49  vsnyder
! Only get isotope ratios in the specified sidebands
!
! Revision 2.32  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.31  2008/10/01 21:01:54  vsnyder
! Require state vector to include all molecules; 2.30 was a bad idea
!
! Revision 2.30  2006/04/05 21:46:44  vsnyder
! Allow state vector not to include all molecules
!
! Revision 2.29  2005/09/07 23:07:26  vsnyder
! Find spectral parameters in FwdModelIn or FwdModelExtra
!
! Revision 2.28  2005/09/03 01:21:33  vsnyder
! Spectral parameter offsets stuff
!
! Revision 2.27  2005/08/03 18:04:09  vsnyder
! Some spectroscopy derivative stuff
!
! Revision 2.26  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.25  2005/05/05 01:11:49  vsnyder
! Don't inquire the size of fwdModelConf%beta_group(b)%pfa(s)%molecules if
! it's not associated.
!
! Revision 2.24  2005/03/15 19:55:51  vsnyder
! Spiff up a dump
!
! Revision 2.23  2005/02/16 23:16:49  vsnyder
! Revise data structures for split-sideband PFA
!
! Revision 2.22  2004/12/28 00:28:43  vsnyder
! Remove unused procedure and unreferenced use names
!
! Revision 2.21  2004/12/13 20:38:24  vsnyder
! Moved stuff that doesn't depend on state vector to ForwardModelConfig
!
! Revision 2.20  2004/11/05 19:37:23  vsnyder
! Move some stuff to ForwardModelConfig%DeriveFromForwardModel
!
! Revision 2.19  2004/11/04 03:40:54  vsnyder
! Index spectroscopy catalog by molecule instead of searching
!
! Revision 2.18  2004/11/01 20:26:35  vsnyder
! Reorganization of representation for molecules and beta groups; PFA may be broken for now
!
