! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelConfig
!=============================================================================

! Set up the forward model configuration, except for actually processing
! the command.


  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dump_0, only: DUMP
  use Intrinsic, only: Lit_indices
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate,&
    & MLSMSG_Error
  use MLSSignals_M, only: DestroySignal, GetSignalName, MaxSigLen, Signal_T
  use Output_M, only: Output
  use String_Table, only: Display_String
  use VGridsDatabase, only: DestroyVGridContents, VGrid_T

  implicit NONE
  private

  ! Public procedures:
  interface Dump
    module procedure Dump_ForwardModelConfigs
  end interface Dump

  public :: AddForwardModelConfigToDatabase, DestroyFWMConfigDatabase, Dump

  ! Public Types:

  type, public :: ForwardModelConfig_T
    integer :: fwmType        ! l_linear, l_full or l_scan
    logical :: Atmos_Der      ! Do atmospheric derivatives
    logical :: Do_Conv        ! Do convolution
    logical :: Do_Freq_Avg    ! Do Frequency averaging
    integer, dimension(:), pointer :: molecules=>NULL() ! Which molecules to consider
    logical, dimension(:), pointer :: moleculeDerivatives=>NULL() ! Want jacobians
    type (Signal_T), dimension(:), pointer :: signals=>NULL()
    logical :: Spect_Der      ! Do spectroscopy derivatives
    logical :: Temp_Der       ! Do temperature derivatives
    type(vGrid_T), pointer :: integrationGrid=>NULL() ! Zeta grid for integration
    type(vGrid_T), pointer :: tangentGrid=>NULL()     ! Zeta grid for integration
    integer :: surfaceTangentIndex  ! Index in Tangentgrid of Earth's surface
    integer :: phiWindow            ! Window size for examining stuff
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

    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: Database
    type (ForwardModelConfig_T), intent(in) :: Item

    ! Local variables
    type (ForwardModelConfig_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddForwardModelConfigToDatabase = newSize
  end function AddForwardModelConfigToDatabase

  ! --------------------------  DestroyForwardModelConfigDatabase  -----
  subroutine DestroyFWMConfigDatabase ( Database )
    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: Database

    ! Local variables
    integer :: Config                   ! Loop counter
    integer :: Signal                   ! Loop counter
    integer :: Status                   ! Flag

    if ( associated(database) ) then
      do config = 1, size(database)
        if ( associated(database(config)%signals) ) then
          do signal = 1, size(database(config)%signals)
            call destroySignal ( database(config)%signals(signal) )
          end do
          deallocate ( database(config)%signals, stat=status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_Deallocate // "database%signals" )
        end if
        ! Don't destroy integrationGrid and tangentGrid.  Assume they will
        ! be (or already are) destroyed by destroyVGridDatabase.
        call deallocate_test ( database(config)%molecules, &
          & "database(config)%molecules", moduleName )
        call deallocate_test ( database(config)%moleculeDerivatives, &
          & "database(config)%moleculeDerivatives", moduleName )
      end do

      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate // "Database" )
    end if
  end subroutine DestroyFWMConfigDatabase

  ! =====     Private Procedures     =====================================
  ! ------------------------------------  DUMP_FOWARDMODELCONFIGS  -----
  subroutine Dump_ForwardModelConfigs ( Database )
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
        call output ( '  Do_freq_avg:' )
        call output ( database(i)%do_freq_avg, advance='yes' )
        call output ( '  Spect_der:' )
        call output ( database(i)%spect_der, advance='yes' )
        call output ( '  Temp_der:' )
        call output ( database(i)%temp_der, advance='yes' )
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
