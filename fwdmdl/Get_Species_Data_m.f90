! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Get_Species_Data_m

  ! Get species data for the molecules in the beta groups.

  implicit NONE
  private
  public :: Get_Species_Data, Destroy_Species_Data

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)), save :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !-----------------------------------------------------------------------

contains

  subroutine Get_Species_Data ( FwdModelConf, FwdModelIn, FwdModelExtra )

    ! Fill in the Beta_Groups' isotope ratios in FwdModelConf.
    ! Get vector quantities.

    use Allocate_Deallocate, only: Allocate_Test
    use ForwardModelConfig, only: Dump, ForwardModelConfig_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: LIT_INDICES, L_ISOTOPERATIO, L_NONE, L_VMR
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
      & MLSMSG_Warning
    use MLSSignals_m, only: GetRadiometerFromSignal
    use SpectroscopyCatalog_m, only: Catalog, Dump, Empty_Cat, Line_t, &
      & Lines, MostLines
    use String_table, only: GET_STRING
    use Toggles, only: Switches
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T

    type(forwardModelConfig_t), intent(inout) :: FwdModelConf ! Fills Beta_Group
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra

    integer :: B         ! Index for beta groups 
    integer, save :: DumpFWM = -1
    type (VectorValue_T), pointer :: F  ! An arbitrary species
    integer :: M         ! Index for molecules in beta groups, or size thereof

    if ( dumpFWM < 0 ) then ! done only once
      if ( index(switches,'fwmg') > 0 ) dumpFWM = 1 ! Dump but don't stop
      if ( index(switches,'fwmG') > 0 ) dumpFWM = 2 ! Dump and stop
    end if

    ! Get isotope ratios for molecules in a beta group, else 1.0 if not a group
    do b = 1, size(fwdModelConf%beta_group)
      if ( fwdModelConf%beta_group(b)%group ) then ! A molecule group
        ! First LBL molecules' ratios
        do m = 1, size(fwdModelConf%beta_group(b)%lbl_molecules)
          f => getQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_isotopeRatio, &
            & molecule=fwdModelConf%beta_group(b)%lbl_molecules(m), &
            & noError=.TRUE., config=fwdModelConf )
          fwdModelConf%beta_group(b)%lbl_ratio(m)     = 1.0
          if ( associated ( f ) ) & ! Have an isotope ratio
            & fwdModelConf%beta_group(b)%lbl_ratio(m) = f%values(1,1)
        end do ! m
        ! Now PFA molecules' ratios
        do m = 1, size(fwdModelConf%beta_group(b)%pfa_molecules)
          f => getQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_isotopeRatio, &
            & molecule=fwdModelConf%beta_group(b)%pfa_molecules(m), &
            & noError=.TRUE., config=fwdModelConf )
          fwdModelConf%beta_group(b)%pfa_ratio(m)     = 1.0
          if ( associated ( f ) ) & ! Have an isotope ratio
            & fwdModelConf%beta_group(b)%pfa_ratio(m) = f%values(1,1)
        end do ! m
!     else ! Not a molecule group, but this is handled in ForwardModelSupport
!       fwdModelConf%beta_group(b)%lbl_ratio(1)     = 1.0
!       fwdModelConf%beta_group(b)%pfa_ratio(1)     = 1.0
      end if
    end do ! b

    ! Get state vector quantities for species
    do b = 1, size(fwdModelConf%beta_group)
      fwdModelConf%beta_group(b)%qty%qty => getQuantityForForwardModel ( &
        &  fwdModelIn, fwdModelExtra, quantityType=l_vmr, molIndex=b,    &
        &  config=fwdModelConf, radiometer=fwdModelConf%signals(1)%radiometer, &
        &  foundInFirst=fwdModelConf%beta_group(b)%qty%foundInFirst )
    end do ! b

    if ( dumpFWM > 0 ) then
      call dump ( fwdModelConf, moduleName )
      call dump ( fwdModelConf%catalog )
      if ( dumpFWM > 1 ) stop
    end if

  end subroutine Get_Species_Data

  ! -------------------------------------------  Destroy_Species_Data  -----
  subroutine Destroy_Species_Data ( FwdModelConf )

  ! Destroy the spectroscopy catalog extract that was allocated by
  ! Get_Species_Data.

    use Allocate_Deallocate, only: Deallocate_Test
    use ForwardModelConfig, only: ForwardModelConfig_t
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error

    type(forwardModelConfig_t), intent(inout) :: FwdModelConf

  end subroutine Destroy_Species_Data

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Get_Species_Data_m

! $Log$
! Revision 2.20  2004/11/05 19:37:23  vsnyder
! Move some stuff to ForwardModelConfig%DeriveFromForwardModel
!
! Revision 2.19  2004/11/04 03:40:54  vsnyder
! Index spectroscopy catalog by molecule instead of searching
!
! Revision 2.18  2004/11/01 20:26:35  vsnyder
! Reorganization of representation for molecules and beta groups; PFA may be broken for now
!
