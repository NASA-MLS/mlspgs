! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Get_Species_Data_m

  ! Get species data for the molecules in the beta groups.

  implicit NONE
  private
  public :: Get_Species_Data

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

    use ForwardModelConfig, only: Dump, ForwardModelConfig_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: L_ISOTOPERATIO, L_VMR
    use SpectroscopyCatalog_m, only: Dump
    use Toggles, only: Switches
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T

    type(forwardModelConfig_t), intent(inout) :: FwdModelConf ! Fills Beta_Group
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra

    integer :: B         ! Index for beta groups 
    integer, save :: DumpFWM = -1
    type (VectorValue_T), pointer :: F  ! An arbitrary species
    integer :: M         ! Index for molecules in beta groups, or size thereof
    integer :: S         ! Sideband index, 1 = LSB, 2 = USB

    if ( dumpFWM < 0 ) then ! done only once
      if ( index(switches,'fwmg') > 0 ) dumpFWM = 1 ! Dump but don't stop
      if ( index(switches,'fwmG') > 0 ) dumpFWM = 2 ! Dump and stop
    end if

    ! Get isotope ratios for molecules in a beta group, else 1.0 if not a group
    do b = 1, size(fwdModelConf%beta_group)
      if ( fwdModelConf%beta_group(b)%group ) then ! A molecule group
        ! First LBL molecules' ratios
        do s = 1, 2
          do m = 1, size(fwdModelConf%beta_group(b)%lbl(s)%molecules)
            f => getQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
              & quantityType=l_isotopeRatio,                             &
              & molecule=fwdModelConf%beta_group(b)%lbl(s)%molecules(m), &
              & noError=.TRUE., config=fwdModelConf )
            fwdModelConf%beta_group(b)%lbl(s)%ratio(m)     = 1.0
            if ( associated ( f ) ) & ! Have an isotope ratio
              & fwdModelConf%beta_group(b)%lbl(s)%ratio(m) = f%values(1,1)
          end do ! m
          ! Now PFA molecules' ratios
          do m = 1, size(fwdModelConf%beta_group(b)%pfa(s)%molecules)
            f => getQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
              & quantityType=l_isotopeRatio,                             &
              & molecule=fwdModelConf%beta_group(b)%pfa(s)%molecules(m), &
              & noError=.TRUE., config=fwdModelConf )
            fwdModelConf%beta_group(b)%pfa(s)%ratio(m)     = 1.0
            if ( associated ( f ) ) & ! Have an isotope ratio
              & fwdModelConf%beta_group(b)%pfa(s)%ratio(m) = f%values(1,1)
          end do ! m
        end do ! s
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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Get_Species_Data_m

! $Log$
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
