! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ForwardModelWrappers

  ! This module contains a wrapper routine for calling the various forward
  ! models we have.
  
  use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T, &
    & FORWARDMODELSTATUS_T
  use VectorsModule, only: VECTOR_T
  use MatrixModule_1, only: MATRIX_T
  use ForwardModelConfig, only: FORWARDMODELCONFIG_T
  use Init_tables_module, only: L_LINEAR, L_SCAN, L_FULL
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error

  use ForwardModelInterface, only: FULLFORWARDMODEL
  use LinearizedForwardModel_m, only: LINEARIZEDFORWARDMODEL
!  use ScanModelModule, only: SCANFORWARDMODEL

  implicit none
  private

  public :: ForwardModel

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------
  
contains ! ============= Public Procedures ==========================

  !----------------------------------------- ForwardModel -----------
  subroutine ForwardModel ( ForwardModelConfig, FwdModelIn, FwdModelExtra, &
    FwdModelOut, Ifm, fmStat, Jacobian )

    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: FORWARDMODELCONFIG
    type(vector_T), intent(in) ::  FWDMODELIN, FwdModelExtra
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN

    ! Executable code
    select case (ForwardModelConfig%fwmType)
    case ( l_full )
      call FullForwardModel ( ForwardModelConfig, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
    case ( l_linear )
      call LinearizedForwardModel ( ForwardModelConfig, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
    case ( l_scan )
      !call ScanForwardModel ( ForwardModelConfig, FwdModelIn, FwdModelExtra, &
      !  FwdModelOut, Ifm, fmStat, Jacobian )
    case default ! Shouldn't get here if parser etc. worked
    end select
  end subroutine ForwardModel

end module ForwardModelWrappers

! $Log$
! Revision 2.3  2001/04/28 17:48:48  livesey
! Removed some unnecessary checks
!
! Revision 2.2  2001/04/26 23:54:26  livesey
! Now uses linear forward model
!
! Revision 2.1  2001/04/26 19:47:41  livesey
! First version
!
