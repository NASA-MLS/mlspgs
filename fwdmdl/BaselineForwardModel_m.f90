! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module BaselineForwardModel_m

  use ForwardModelConfig, only: FORWARDMODELCONFIG_T
  use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T, FORWARDMODELSTATUS_T
  use MLSCommon, only: RP
  use MLSSignals_m, only: SIGNALS, SIGNAL_T
  use VectorsModule, only: VECTOR_T, VECTORVALUE_T, GETVECTORQUANTITYBYTYPE
  use MatrixModule_1, only: MATRIX_T

  ! This module contains a special forward model for baseline related effects.

  implicit none
  private
  public :: BaselineForwardModel

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains ! ======================================== BaselineForwardModel ======

  subroutine BaselineForwardModel ( FwdModelConf, FwdModelIn, FwdModelExtra, &
    & FwdModelOut, oldIFM, fmStat, jacobian )
    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: oldIfm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian

  end subroutine BaselineForwardModel

end module BaselineForwardModel_m
  
! $Log$
! Revision 2.1  2001/10/02 16:52:00  livesey
! Very early version!
!
