! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module HybridForwardModel_m

  ! This is a 'hybrid' forward model which uses the linear forward model
  ! for one sideband, and the non-linear forward model for the other.

  implicit none
  private
  public :: HybridForwardModel

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = & 
    "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  subroutine HybridForwardModel ( fmConf, FwdModelIn, FwdModelExtra,&
    & FwdModelOut, Ifm, fmStat, Jacobian, Vectors )

    ! Import stuff
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T, CLONEVECTOR, ADDTOVECTOR, &
      & DESTROYVECTORINFO
    use ForwardModelConfig, only: FORWARDMODELCONFIG_T
    use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T, &
      & FORWARDMODELINTERMEDIATE_T
    use Intrinsic, only: L_LINEAR, L_FULL
    use FullForwardModel_m, only: FULLFORWARDMODEL
    use LinearizedForwardModel_m, only: LINEARIZEDFORWARDMODEL
    use MatrixModule_1, only: MATRIX_T, COPYMATRIX, ADDTOMATRIX, CLEARMATRIX, &
      & CREATEEMPTYMATRIX, DESTROYMATRIX
    use MLSSignals_M, only: SIGNAL_T

    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: FMCONF
    type(vector_T), intent(in) ::  FWDMODELIN
    type(vector_T), intent(in) ::  FWDMODELEXTRA
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN
    type(vector_T), dimension(:), pointer, optional :: VECTORS

    ! Local variables
    integer :: SIDEBAND                 ! Loop counter
    integer :: SIGNAL                   ! Loop counter
    type(forwardModelConfig_T) :: THISCONFIG
    type(signal_T), target, dimension(1) :: THISSIGNAL
    type(Vector_T), target :: LINEARRADIANCE
    type(Matrix_T), target :: LINEARJACOBIAN
    ! Executable code

    ! Setup our temporary stuff
    if ( present(jacobian) ) then
      call CreateEmptyMatrix ( linearJacobian, 0, &
        & jacobian%row%vec, jacobian%col%vec, &
        & .not. jacobian%row%instFirst, .not. jacobian%col%instFirst, &
        & 'linearJacobian' )
    end if
    call CloneVector ( linearRadiance, fwdModelOut )

    ! First we do the linear model.  We do this for both sidebands,
    ! though the second one we update/overwrite with the full forward model
    ! results.
    do sideband = fmConf%linearSideband, -fmConf%linearSideband, -2*fmConf%linearSideband
      ! i.e. do sideband = 1,-1,-2 or -1,1,2
      do signal = 1, size ( fmConf%signals )
        ! Set it up to do this signal
        thisConfig = fmConf
        thisConfig%fwmType = l_linear
        thisConfig%forceFoldedOutput = .true.
        thisSignal(1) = thisConfig%signals(signal)
        thisSignal%sideband = sideband
        thisConfig%signals => thisSignal
        ! We'll want the linear model to supply all the derivatives it can
        ! as it's not much of a cost and we want them from somewhere.
        thisConfig%moleculeDerivatives = .true.
        thisConfig%forceSidebandFraction = .true.
        ! Invoke it, tried to do this with temporary pointers
        ! to avoid so many calls, but was difficult as I didn't want to
        ! make fwdModelOut a target.
        if ( sideband == fmConf%linearSideband ) then
          if ( present(jacobian) ) then
            call LinearizedForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
              & linearRadiance, ifm, fmStat, linearJacobian, vectors=vectors )
          else
            call LinearizedForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
              & linearRadiance, ifm, fmStat, vectors=vectors )
          end if
        else
          if ( present(jacobian) ) then
            call LinearizedForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
              & fwdModelOut, ifm, fmStat, jacobian, vectors )
          else
            call LinearizedForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
              & fwdModelOut, ifm, fmStat, vectors=vectors )
          end if
        end if
      end do                              ! Signal loop
    end do                                ! Sideband loop

    ! At the end of this linearRadiance/linearJacobian contain the results for the
    ! sideband that's purely linear, and fwdModelOut/Jacobian contain that which
    ! we're going to overwrite with the full model.

    ! Now we do the full forward model
    thisConfig = fmConf
    thisConfig%fwmType = l_full
    thisConfig%forceSidebandFraction = .true.
    thisConfig%signals%sideband = - fmConf%linearSideband
    thisConfig%sidebandStart = - fmConf%linearSideband
    thisConfig%sidebandStop = - fmConf%linearSideband
    call FullForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
      & fwdModelOut, ifm, fmStat, Jacobian )
    ! The sidebands must originally have been marked as folded
    fmConf%signals%sideband = 0

    ! Now add the terms together and tidy up
    call AddToVector ( fwdModelOut, linearRadiance )
    call DestroyVectorInfo ( linearRadiance )

    if ( present(jacobian ) ) then
      call AddToMatrix ( jacobian, linearJacobian )
      call DestroyMatrix ( linearJacobian )
    end if

  end subroutine HybridForwardModel

  ! ----------------------------------------------------------------------------

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HybridForwardModel_m

! $Log$
! Revision 2.4  2003/10/28 23:44:25  livesey
! Various deficiencies and bugs fixed.
!
! Revision 2.3  2003/09/11 23:11:28  livesey
! Now includes the vectors argument to push down into the linear model.
!
! Revision 2.2  2003/08/13 00:47:49  livesey
! Cosmetic change
!
! Revision 2.1  2003/07/15 22:10:15  livesey
! First version
!
