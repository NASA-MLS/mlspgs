! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelConfig
!=============================================================================

! Set up the forward model configuration, except for actually processing
! the command.


  use MLSCommon, only: R8
  use MLSSignals_M, only: Signal_T
  use VGridsDatabase, only: VGrid_T

  implicit NONE
  private

  ! Public procedures:
  interface Dump
    module procedure Dump_ForwardModelConfigs
  end interface Dump

  public :: AddForwardModelConfigToDatabase, DestroyFWMConfigDatabase, Dump, &
    & StripForwardModelConfigDatabase

  ! Public Types:

  type, public :: ForwardModelConfig_T
    logical :: globalConfig   ! If set is shared between all chunks
    integer :: fwmType        ! l_linear, l_full or l_scan
    logical :: Atmos_Der      ! Do atmospheric derivatives
    logical :: do_Baseline    ! Do a baseline computation
    logical :: Do_Conv        ! Do convolution
    logical :: Do_Freq_Avg    ! Do Frequency averaging
    integer, dimension(:), pointer :: molecules=>NULL() ! Which molecules to consider
    logical, dimension(:), pointer :: moleculeDerivatives=>NULL() ! Want jacobians
    type (Signal_T), dimension(:), pointer :: signals=>NULL()
    integer :: instrumentModule         ! Module for scan model
    logical :: differentialScan         ! Differential scan model
    logical :: LockBins                 ! Use same l2pc bin for whole chunk
    logical :: Spect_Der      ! Do spectroscopy derivatives
    logical :: Temp_Der       ! Do temperature derivatives
    logical :: skipOverlaps   ! Don't calculate for MAFs in overlap regions
    type(vGrid_T), pointer :: integrationGrid=>NULL() ! Zeta grid for integration
    type(vGrid_T), pointer :: tangentGrid=>NULL()     ! Zeta grid for integration
    integer, dimension(:), pointer :: specificQuantities=>NULL() ! Specific quantities to use
    integer :: surfaceTangentIndex  ! Index in Tangentgrid of Earth's surface
    real (r8) :: phiWindow            ! Window size for examining stuff
    real (r8) :: tolerance          ! Accuracy desired when choosing approximations
    ! CloudForwardModel
    logical :: Default_spectroscopy      !
    integer :: no_cloud_species     ! No of Cloud Species '2'
    integer :: no_model_surfs       ! No of Model surfaces '640'
    integer :: NUM_SCATTERING_ANGLES! No of scattering angles '16'
    integer :: NUM_AZIMUTH_ANGLES   ! No of azmuth angles '8'
    integer :: NUM_AB_TERMS         ! No of AB terms '50'
    integer :: NUM_SIZE_BINS        ! No of size bins '40'
    integer :: cloud_der            ! Compute cloud sensitivity in cloud models.
    integer :: cloud_width          ! Flag for cloud horizontal extend.
    integer :: cloud_fov            ! Flag for cloud model field-of-view averaging.
    integer :: NameFragment         ! For e.g. restricting bins in l2pc
  end type ForwardModelConfig_T

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  ! ----------------------------  AddForwardModelConfigToDatabase  -----
  integer function AddForwardModelConfigToDatabase ( Database, Item )

    ! Add a quantity template to a database, or create the database if it
    ! doesn't yet exist

    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
      & MLSMSG_Deallocate, MLSMSG_Error

    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: Database
    type (ForwardModelConfig_T), intent(in) :: Item

    ! Local variables
    type (ForwardModelConfig_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddForwardModelConfigToDatabase = newSize
  end function AddForwardModelConfigToDatabase

  ! --------------------------  StripForwardModelConfigDatabase --------
  subroutine StripForwardModelConfigDatabase ( database )
    ! This routine removes the non-global forward model configs from the database
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
      & MLSMSG_Deallocate, MLSMSG_Error

    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: DATABASE

    ! Local variables
    type (ForwardModelConfig_T), dimension(:), pointer :: TMPDATABASE
    integer :: CONFIG                   ! Loop counter
    integer :: STATUS

    ! Executable code
    ! Clear out dying configs
    do config = 1, size ( database )
      if ( .not. database(config)%globalConfig ) &
        & call DestroyOneForwardModelConfig ( database(config) )
    end do

    ! Create new database in tmp space, pack old one into
    allocate ( tmpDatabase ( count ( database%globalConfig ) ), STAT=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate // 'tmpDatabase' )
    tmpDatabase = pack ( database, database%globalConfig )

    ! Destroy old database, then point to new one
    deallocate ( database, STAT=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate // 'database' )
    database => tmpDatabase
  end subroutine StripForwardModelConfigDatabase

  ! --------------------------  DestroyForwardModelConfigDatabase  -----
  subroutine DestroyFWMConfigDatabase ( Database )

    use Allocate_Deallocate, only: Deallocate_Test
    use MLSMessageModule, only: MLSMessage,  MLSMSG_Deallocate, MLSMSG_Error
    use MLSSignals_M, only: DestroySignal

    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: Database

    ! Local variables
    integer :: Config                   ! Loop counter
    integer :: Signal                   ! Loop counter
    integer :: Status                   ! Flag

    if ( associated(database) ) then
      do config = 1, size(database)
        call DestroyOneForwardModelConfig ( database(config) )
      end do

      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate // "Database" )
    end if
  end subroutine DestroyFWMConfigDatabase

  ! =====     Private Procedures     =====================================

  ! ------------------------------------ DestroyOneForwardModelConfig --
  subroutine DestroyOneForwardModelConfig ( config )
    use Allocate_Deallocate, only: Deallocate_Test
    use MLSSignals_M, only: DestroySignal
    use MLSMessageModule, only: MLSMSG_Deallocate, MLSMSG_Error, MLSMessage

    ! Dummy arguments
    type ( ForwardModelConfig_T), intent(inout) :: config

    ! Local variables
    integer :: SIGNAL                   ! Loop counter
    integer :: STATUS                   ! Flag from allocate etc.

    ! Executable code
    if ( associated(config%signals) ) then
      do signal = 1, size(config%signals)
        call destroySignal ( config%signals(signal), &
          & justChannels=.true. )
      end do
      deallocate ( config%signals, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate // "database%signals" )
    end if
    ! Don't destroy integrationGrid and tangentGrid.  Assume they will
    ! be (or already are) destroyed by destroyVGridDatabase.
    call deallocate_test ( config%molecules, &
      & "config%molecules", moduleName )
    call deallocate_test ( config%moleculeDerivatives, &
      & "config%moleculeDerivatives", moduleName )
    call Deallocate_test ( config%specificQuantities, &
      & "config%specificQuantities", ModuleName )
  end subroutine DestroyOneForwardModelConfig

  ! ------------------------------------  Dump_FowardModelConfigs  -----
  subroutine Dump_ForwardModelConfigs ( Database )

    use Dump_0, only: DUMP
    use Intrinsic, only: Lit_indices
    use MLSSignals_M, only: GetSignalName, MaxSigLen
    use Output_M, only: Output
    use String_Table, only: Display_String

    type (ForwardModelConfig_T), pointer, dimension(:) :: Database

    ! Local variables
    integer :: I, J                          ! Loop counters
    character (len=MaxSigLen) :: SignalName  ! A line of text

    ! executable code
    if ( associated(database) ) then
      do i = 1, size(database)
        call output ( 'FowardModelConfig: ' )
        call output ( i, advance = 'yes' )
        call output ( '  Atmos_der:' )
        call output ( database(i)%atmos_der, advance='yes' )
        call output ( '  Do_conv:' )
        call output ( database(i)%do_conv, advance='yes' )
        call output ( '  Do_Baseline:' )
        call output ( database(i)%do_Baseline, advance='yes' )
        call output ( '  Do_freq_avg:' )
        call output ( database(i)%Default_spectroscopy, advance='yes' )
        call output ( '  Default_spectroscopy:' )
        call output ( database(i)%do_freq_avg, advance='yes' )
        call output ( '  Spect_der:' )
        call output ( database(i)%spect_der, advance='yes' )
        call output ( '  Temp_der:' )
        call output ( database(i)%temp_der, advance='yes' )
        call output ( '  SkipOverlaps:' )
        call output ( database(i)%skipOverlaps, advance='yes' )
        call output ( '  Molecules: ', advance='yes' )
        do j = 1, size(database(i)%molecules)
          call output ( '    ' )
          call display_string(lit_indices(database(i)%molecules(j)))
          if (database(i)%moleculeDerivatives(j)) then
            call output (' compute derivatives', advance='yes')
          else
            call output (' no derivatives', advance='yes')
          end if
        end do
        call output ( '  Signals:', advance='yes')
        do j = 1, size(database(i)%signals)
          call output ( '    ' )
          !call GetSignalName( signal=database(i)%signals(j), signalName)
          !??? Sort this out later!
          ! call output ( signalName//' channelIncluded:', advance='yes')
          call dump ( database(i)%signals(j)%channels )
        end do
      end do
    end if
  end subroutine Dump_ForwardModelConfigs

end module ForwardModelConfig

! $Log$
! Revision 2.9  2002/08/21 23:53:57  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.8  2002/07/17 06:02:13  livesey
! New config elements for hdf5 l2pcs
!
! Revision 2.7  2002/06/12 17:00:49  livesey
! Changed phiWindow to float
!
! Revision 2.6  2002/03/07 17:17:49  livesey
! Removed frqGap
!
! Revision 2.5  2002/02/14 23:01:35  livesey
! Added justChannels in call to destroySignal
!
! Revision 2.4  2002/02/13 00:09:24  livesey
! Added differential Scan model
!
! Revision 2.3  2001/11/15 23:50:11  jonathan
! rename DF_spectroscopy to default_spectroscopy
!
! Revision 2.2  2001/11/15 20:56:38  jonathan
! add df_spectroscopy
!
! Revision 2.1  2001/10/02 20:37:09  livesey
! Added do_baseline
!
! Revision 2.0  2001/09/17 20:26:25  livesey
! New forward model
!
! Revision 1.15  2001/09/04 15:59:01  jonathan
! add cloud_fov, jonathan
!
! Revision 1.14  2001/07/17 22:38:05  jonathan
! add cloud_width, jonathan/paul
!
! Revision 1.13  2001/07/16 22:07:57  jonathan
! change cloud_der to integer-type, jonathan
!
! Revision 1.12  2001/07/06 18:55:13  jonathan
! Modified for cloud model, Paul/Jonathan
!
! Revision 1.11  2001/06/21 15:05:53  livesey
! Added tolerance
!
! Revision 1.10  2001/05/31 23:07:45  livesey
! Added cloud_der
!
! Revision 1.9  2001/05/25 20:26:09  livesey
! Added skipOverlaps option
!
! Revision 1.8  2001/05/14 23:17:35  livesey
! Added frqGap parameter
!
! Revision 1.7  2001/05/03 23:07:02  livesey
! Added scan model stuff
!
! Revision 1.6  2001/05/02 20:30:36  livesey
! Removed frequency from config
!
! Revision 1.5  2001/04/26 02:36:52  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 1.4  2001/04/21 01:08:57  vsnyder
! Deallocate Molecules and MoleculeDerivatives in DestroyFWMConfigDatabase
!
! Revision 1.3  2001/04/12 17:00:08  vsnyder
! Comment out a line with an undefined variable on it
!
! Revision 1.2  2001/04/10 22:17:05  livesey
! Renamed module
!
! Revision 1.1  2001/04/07 01:56:25  vsnyder
! Initial commit
!
