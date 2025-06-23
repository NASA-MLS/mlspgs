! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module HybridForwardModel_m

  ! This is a 'hybrid' forward model which uses the linear forward model
  ! for one sideband, and the non-linear forward model for the other.

  implicit none
  private
  public :: HybridForwardModel

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  subroutine HybridForwardModel ( fmConf, FwdModelIn, FwdModelExtra,&
    & FwdModelOut, fmStat, Jacobian, Vectors )

    ! Import stuff
    use VectorsModule, only: VECTOR_T, CLONEVECTOR, ADDTOVECTOR, &
      & DESTROYVECTORINFO
    use ForwardModelConfig, only: FORWARDMODELCONFIG_T
    use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T
    use Intrinsic, only: L_LINEAR, L_FULL
    use FullForwardModel_m, only: FULLFORWARDMODEL
    use LinearizedForwardModel_m, only: LINEARIZEDFORWARDMODEL
    use MatrixModule_1, only: MATRIX_T, ADDTOMATRIX, CREATEEMPTYMATRIX, &
      & DESTROYMATRIX
    use MLSSignals_M, only: SIGNAL_T

    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: FMCONF
    type(vector_T), intent(in) ::  FWDMODELIN
    type(vector_T), intent(in) ::  FWDMODELEXTRA
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN
    type(vector_T), dimension(:), target, optional :: VECTORS

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
              & linearRadiance, fmStat, linearJacobian, vectors=vectors )
          else
            call LinearizedForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
              & linearRadiance, fmStat, vectors=vectors )
          end if
        else
          if ( present(jacobian) ) then
            call LinearizedForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
              & fwdModelOut, fmStat, jacobian, vectors )
          else
            call LinearizedForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
              & fwdModelOut, fmStat, vectors=vectors )
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
    thisConfig%forceFoldedOutput = .true.
    thisConfig%forceSidebandFraction = .true.
    thisConfig%signals%sideband = - fmConf%linearSideband
    thisConfig%sidebandStart = - fmConf%linearSideband
    thisConfig%sidebandStop = - fmConf%linearSideband
    call FullForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
      & fwdModelOut, fmStat, Jacobian )
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HybridForwardModel_m

! $Log$
! Revision 2.10  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.9  2007/10/06 00:01:01  vsnyder
! Delete unused symbols
!
! Revision 2.8  2007/06/29 19:32:42  vsnyder
! Make ForwardModelIntermediate_t private to ScanModelModule
!
! Revision 2.7  2007/04/03 17:43:35  vsnyder
! Replace pointer attribute on VectorDatabase with target attribute
!
! Revision 2.6  2005/06/03 01:59:43  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades
!
! Revision 2.5  2003/10/29 00:44:09  livesey
! Added forceFoldedOuptut to full forward model call.
!
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
