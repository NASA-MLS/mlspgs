! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module PFADataBase_m

  ! Read PFA data.  Build a database.  Provide for access to it.
  ! Write PFA data.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use HighOutput, only: OutputNamedValue
  use MLSKinds, only: R4
  use MLSCommon, only: FileNameLen
  use MLSSignals_M, only: MaxSigLen, Signal_T
  use Molecules, only: First_Molecule, Last_Molecule
  use VGridsDatabase, only: VGrid_T
  ! use FOR HDF5 INTENTIONALLY LAST TO AVOID LONG LF95 COMPILE TIMES
  use HDF5, only: hid_t

  implicit NONE
  private
  public :: PFAData_t, PFAFile_t
  public :: PFAData, PFAFiles, RK
  public :: AddPFADatumToDatabase
  public :: Destroy_All_PFAData_Arrays, Destroy_PFADataBase, Destroy_PFADatum
  public :: Destroy_PFADatum_Arrays, Destroy_PFAFile, Destroy_PFAFiles
  public :: Dump, Dump_PFADataBase, Dump_PFADatum
  public :: Dump_PFAFileDataBase, Dump_PFAFileDatum, Dump_PFAStructure
  public :: Flush_PFADataBase, GetGroupName, HookTableToFindPFA
  public :: Process_PFA_File, Read_PFADataBase
  public :: Test_And_Fetch_PFA, Write_PFADatum, Write_PFADataBase

  interface Dump
    module procedure Dump_PFADatum, Dump_PFAFileDatum
  end interface Dump

  interface Process_PFA_File
    module procedure Process_PFA_File_datum
    module procedure Process_PFA_File_node, Process_PFA_File_name
  end interface Process_PFA_File

  integer, parameter :: RK = r4 ! Kind of real fields in PFAData_t

  type PFAData_T
    integer :: Name = 0                ! of the pfaData spec
    integer :: Channel                 ! Duh
    integer :: FilterFile = 0          ! Index in string table
    integer(hid_t) :: HDF5_GroupID     ! HDF5 group id if open
    integer :: Molecule                ! Molecule index
    logical :: Open = .false.          ! "HDF5 group is open"
    integer :: Signal                  ! Sub-rosa index of the signal string
    integer :: SignalIndex = 0         ! in Signals database
    integer :: SpectroscopyFile = 0    ! Index in string table
    type(signal_t) :: TheSignal        ! The signal, with channels
                                       ! and sidebands added
    type(vGrid_t) :: TGrid ! shallow copy from VGrids database of Log temperatures
    type(vGrid_t) :: VGrid ! shallow copy from VGrids database of vertical grid
    real(rk) :: Vel_Rel                ! vel_lin / c
    integer :: WhichLines = 0          ! 0 = Lines for channel only,
                                       ! 1 = AllLinesInCatalog,
                                       ! 2 = AllLinesForRadiometer,
    integer :: FileIndex = 0           ! Index in PFAFiles
    integer :: GroupIndex = 0          ! Index in PFAFiles%ix
    real(rk), pointer :: Absorption(:,:) => NULL() ! Ln Absorption data, T x P
    real(rk), pointer :: dAbsDwc(:,:) => NULL()    ! d Ln Absorption / d wc data
    real(rk), pointer :: dAbsDnc(:,:) => NULL()    ! d Ln Absorption / d nc data
    real(rk), pointer :: dAbsDnu(:,:) => NULL()    ! d Ln Absorption / d nu data
  end type PFAData_T

  type(PFAData_t), pointer, save :: PFAData(:) => NULL()

  type PFAFile_T ! For all the PFA tables in one HDF file
    integer :: FileName = 0            ! Sub-rosa index of file name
    character(len=FileNameLen) :: nameString
    integer(hid_t) :: HDF5_FileID      ! HDF5 file ID if open
    logical :: Open = .false.          ! "HDF5 file is open"
    integer, pointer :: Ix(:) => NULL() ! Indices in PFAData
  end type PFAFile_T

  ! All PFA files mentioned in GlobalSettings
  type(PFAFile_t), pointer, save :: PFAFiles(:) => NULL()


  ! Structures for finding PFA data quickly.  There are four levels of tables,
  ! indexed respectively by molecule, signal, sideband and channel.  We do this
  ! instead of having a rank-4 array because some molecules have no PFA data,
  ! and of those that do, some have no data for some signals, and of those that
  ! do, some have no data for one sideband, and then the number of channels is
  ! different for different signals.

  type, public :: BySideband_T
    integer, pointer :: C(:) => NULL() ! Indexed by channel
  end type BySideband_T

  type, public :: BySignal_T
    type(bySideband_T) :: SB(2)        ! Indexed by sideband, 1 = LSB, 2 = USB
  end type BySignal_T

  type, public :: ByMolecule_T
    type(bySignal_t), pointer :: S(:) => NULL()   ! Indexed by signal
  end type ByMolecule_T

  type(byMolecule_t), public, save :: FindPFA(first_molecule:last_molecule)

  integer, parameter, public :: MolNameLen = 31 ! Length for Molecule names in files

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------  AddPFADatumToDatabase  -----
  integer function AddPFADatumToDatabase ( DATABASE, ITEM )

  ! This routine adds a PFA Datum to a database of PFA Data, creating the
  ! database if necessary.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (PFAData_T), dimension(:), pointer :: DATABASE
    type (PFAData_T), intent(in) :: ITEM

    ! Local variables
    type (PFAData_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddPFADatumToDatabase = newSize
  end function AddPFADatumToDatabase

  ! ---------------------------------  Destroy_All_PFAData_Arrays  -----
  subroutine Destroy_All_PFAData_Arrays
    integer :: I ! Loop inductor, subscript
    if ( .not. associated(PFAData) ) return
    do i = 1, ubound(PFAData,1)
      call Destroy_PFADatum_Arrays ( PFAData(i) )
    end do
  end subroutine

  ! ----------------------------------------  Destroy_PFADataBase  -----
  subroutine Destroy_PFADataBase
    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I, M, S, SB
    if ( .not. associated(pfaData) ) return
    do i = 1, ubound(PFAData,1)
      call destroy_PFADatum ( pfaData(i) )
    end do
    call destroy_PFAFiles
    s = size(PFAData) * storage_size(PFAData) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(PFAData(1)), addr)
    deallocate ( PFAData, stat=i )
    call test_deallocate ( i, moduleName, 'PDAData', s, address=addr )
    do m = first_molecule, last_molecule
      if ( associated(findPFA(m)%s) ) then
        do s = lbound(findPFA(m)%s,1), ubound(findPFA(m)%s,1)
          do sb = 1, 2
            call deallocate_test ( findPFA(m)%s(s)%sb(sb)%c, &
              & 'FindPFA(m)%s(s)%sb(sb)%c', moduleName )
          end do
        end do
        s = size(findPFA(m)%s) * storage_size(findPFA(m)%s) / 8
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(findPFA(m)%s(1)), addr)
        deallocate ( findPFA(m)%s, stat=i )
        call test_deallocate ( i , moduleName, 'findPFA(m)%s', s, address=addr )
      end if
    end do
  end subroutine Destroy_PFADataBase

  ! -------------------------------------------  Destroy_PFADatum  -----
  subroutine Destroy_PFADatum ( PFADatum )
    type(PFAData_T), intent(inout) :: PFADatum
    call deallocate_test ( PFADatum%theSignal%channels, 'PFADatum...Channels', &
      & moduleName )
    call Destroy_PFADatum_Arrays ( PFADatum )
  end subroutine Destroy_PFADatum

  ! ------------------------------------  Destroy_PFADatum_Arrays  -----
  subroutine Destroy_PFADatum_Arrays ( PFADatum )
    type(PFAData_T), intent(inout) :: PFADatum
    call deallocate_test ( PFADatum%absorption, 'absorption', moduleName )
    call deallocate_test ( PFADatum%dAbsDwc, 'dAbsDwc', moduleName )
    call deallocate_test ( PFADatum%dAbsDnc, 'dAbsDnc', moduleName )
    call deallocate_test ( PFADatum%dAbsDnu, 'dAbsDnu', moduleName )
  end subroutine Destroy_PFADatum_Arrays

  ! --------------------------------------------  Destroy_PFAFile  -----
  subroutine Destroy_PFAFile ( PFAFile )
  ! Destroy one PFAFile_t object
    type(PFAFile_t), intent(inout) :: PFAFile
    integer :: I
    if ( associated(PFAFile%ix) ) then
      do i = 1, ubound(PFAFile%ix,1)
        call destroy_PFADatum ( PFAData(PFAFile%ix(i)) )
      end do
      call deallocate_test ( PFAFile%ix, 'PFAFile%ix', moduleName )
    end if
  end subroutine Destroy_PFAFile

  ! -------------------------------------------  Destroy_PFAFiles  -----
  subroutine Destroy_PFAFiles
  ! Destroy the elements of the PFAFiles array, then deallocate it.
    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I, S
    if ( associated(PFAFiles) ) then
      do i = 1, ubound(PFAFiles,1)
        call destroy_PFAFile ( PFAFiles(i) )
      end do
      s = size(PFAFiles) * storage_size(PFAFiles) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(PFAFiles(1)), addr)
      deallocate ( PFAFiles, stat=i )
      call test_deallocate ( i, moduleName, 'PFAFiles', s, address=addr )
    end if
  end subroutine Destroy_PFAFiles

  ! -------------------------------------------  Dump_PFADataBase  -----
  subroutine Dump_PFADataBase ( Details )
    use Output_m, only: Output
    integer, intent(in), optional :: Details
    ! Local Variables
    integer            :: I
    real               :: total
    if ( .not. associated(pfaData) ) then
      call output ( 'No PFA Database to dump', advance='yes' )
      return
    end if
    total = 0.
    do i = 1, ubound(PFAData,1)
      call dump_PFADatum ( pfaData(i), details, i )
      if (   associated(pfaData(i)%absorption) ) &
        & total = total + &
        & product(shape(pfaData(i)%absorption))
      if (   associated(pfaData(i)%dAbsDwc) ) &
        & total = total + &
        & product(shape(pfaData(i)%dAbsDwc))
      if (   associated(pfaData(i)%dAbsDnc) ) &
        & total = total + &
        & product(shape(pfaData(i)%dAbsDnc))
      if (   associated(pfaData(i)%dAbsDnu) ) &
        & total = total + &
        & product(shape(pfaData(i)%dAbsDnu))
    end do
    call outputNamedValue( 'PFA Database Total Memory Footprint', &
      & storage_size(PFAData(1)%Vel_Rel) / 8 * total )
  end subroutine Dump_PFADataBase

  ! ----------------------------------------------  Dump_PFADatum  -----
  subroutine Dump_PFADatum ( PFADatum, Details, Index )

    use Dump_0, only: Dump
    use Intrinsic, only: Lit_Indices
    use MLSsignals_M, only: Displaysignalname
    use Physics, only: Speedoflight
    use String_Table, only: Display_String
    use Output_M, only: Newline, Output
    use VgridsDatabase, only: Dump

    type(PFAData_t), intent(in) :: PFADatum
    integer, intent(in), optional :: Details ! <=-2 -> Signal and molecule if data
                                             ! -1 -> Complete summary if data
                                             ! 0 -> Signal and molecule only
                                             ! 1 -> Complete summary
                                             ! 2 -> Betas (default)
                                             ! 3 -> Betas and unnamed grids
                                             ! >3 -> Betas, grids and derivatives
    integer, intent(in), optional :: Index   ! In PFA Database

    integer, parameter          :: CK = kind(speedOfLight)
    real(ck)                    :: C = speedOfLight / 1000.0_ck ! km/s
    integer                     :: MyDetails                          
    real                        :: total
    character(len=*), parameter :: WhichLines(0:2) = &
      (/ 'Channel   ', &
       & 'Radiometer', &
       & 'Catalog   ' /)

    myDetails = 2
    if ( present(details) ) myDetails = details

    if ( myDetails < 0 .and. .not. associated(pfaDatum%absorption) ) return
    total = 0.
    if (   associated(PFADatum%absorption) ) &
      & total = total + &
      & product(shape(PFADatum%absorption))
    if (   associated(PFADatum%dAbsDwc) ) &
      & total = total + &
      & product(shape(PFADatum%dAbsDwc))
    if (   associated(PFADatum%dAbsDnc) ) &
      & total = total + &
      & product(shape(PFADatum%dAbsDnc))
    if (   associated(PFADatum%dAbsDnu) ) &
      & total = total + &
      & product(shape(PFADatum%dAbsDnu))
    call outputNamedValue( 'PFA Datum Total Memory Footprint', &
      & storage_size(PFADatum%Vel_Rel) / 8 * total )

    call output ( 'PFA Datum' )
    if ( present(index) ) call output ( index, before=' ' )
    if ( pfaDatum%name /= 0 ) then
      if ( present(index) ) call output ( ': ' )
      call display_string ( pfaDatum%name )
    end if
    if ( myDetails <= -2 .or. myDetails == 0 ) then
      call display_string ( lit_indices(pfaDatum%molecule), before=' for ' )
      call display_string ( pfaDatum%signal, before=' and ' )
      if ( associated(pfaDatum%absorption) ) call output ( ' has data' )
      call newLine
      return
    end if
    call output ( ':', advance='yes' )
    call output ( ' Molecule: ' )
    call display_string ( lit_indices(pfaDatum%molecule), advance='yes' )

    call output ( ' Specified signal: ' )
    call display_string ( pfaDatum%signal, advance='yes' )
    if ( associated(pfaDatum%vGrid%surfs) ) then ! Something has been read
      call output ( ' Actual signal: ' )
      call output ( pfaDatum%signalIndex, after=': ' )
      if ( pfaDatum%theSignal%name /= 0 ) then
        call display_string ( pfaDatum%theSignal%name )
        call output ( ': ' )
      end if
      call displaySignalName ( pfaDatum%theSignal, advance='yes', &
        & channel=pfaDatum%channel )
      if ( pfaDatum%filterFile /= 0 ) &
        & call display_string ( pfaDatum%filterFile, before=' Filter file: ', &
        & advance='yes' )
      if ( pfaDatum%spectroscopyFile /= 0 ) &
        & call display_string ( pfaDatum%spectroscopyFile, before=' Spectroscopy file: ', &
        & strip=.true., advance='yes' )

      call output ( real(pfaDatum%vel_rel*c,rk), before=' Velocity linearization: ', &
        & after='kms, all lines for ' )

      call output ( trim(whichLines(pfaDatum%whichLines)), advance='yes' )

      if ( pfaDatum%tGrid%name /= 0 ) then
        call display_string ( pfaDatum%tGrid%name, before=' TGrid: ' )
      else if ( myDetails > 2 ) then
        call output ( ' TGrid: ' )
        call dump ( pfaDatum%tGrid )
      end if

      if ( pfaDatum%vGrid%name /= 0 ) then
        call display_string ( pfaDatum%vGrid%name, before=' VGrid: ' )
      else if ( myDetails > 2 ) then
        if ( pfaDatum%tGrid%name /= 0 ) call newLine
        call dump ( pfaDatum%vGrid )
      end if

      if ( pfaDatum%tGrid%name /= 0 .and. pfaDatum%vGrid%name == 0 .and. &
        &  myDetails > 2 &
        & .or. pfaDatum%vGrid%name /= 0 ) &
        & call newLine

      if ( myDetails <= 1 ) return

      if ( associated(pfaDatum%absorption) ) &
        & call dump ( pfaDatum%absorption, name=' ln Absorption' )
      if ( myDetails <= 3 ) return
      if ( associated(pfaDatum%dAbsDwc) ) &
        & call dump ( pfaDatum%dAbsDwc, name=' d ln Absorption / d wc' )
      if ( associated(pfaDatum%dAbsDnc) ) &
        & call dump ( pfaDatum%dAbsDnc, name=' d ln Absorption / d nc' )
      if ( associated(pfaDatum%dAbsDnu) ) &
        &  call dump ( pfaDatum%dAbsDnu, name=' d ln Absorption / d nu' )
    end if

  end subroutine Dump_PFADatum

  ! ---------------------------------------  Dump_PFAFileDataBase  -----
  subroutine Dump_PFAFileDataBase ( Details )
    use Output_m, only: Output
    integer, intent(in), optional :: Details ! See Dump_PFAFileDatum

    integer :: I

    if ( .not. associated(PFAFiles) )then
      call output ( 'No PFA File Database to dump', advance='yes' )
      return
    end if

    do i = 1, ubound(PFAFiles,1)
      call dump ( PFAFiles(i), details )
    end do
  end subroutine Dump_PFAFileDataBase

  ! ------------------------------------------  Dump_PFAFileDatum  -----
  subroutine Dump_PFAFileDatum ( PFAFileDatum, Details )
    use MLSsignals_M, only: Maxsiglen
    use Output_M, only: Blanks, Newline, Output
    use String_Table, only: Display_String
    type(PFAFile_t), intent(in) :: PFAFileDatum
    integer, intent(in), optional :: Details ! -4 (default) => file names,
                                             ! -3 => group names too
                                             ! >=-2 => see Dump_PFADatum

    integer :: G, MyDetails
    character(len=molNameLen+1+maxSigLen) :: GroupName

    myDetails = 0
    if ( present(details) ) myDetails = details

    if ( PFAFileDatum%fileName > 0 ) then
      call display_string ( PFAFileDatum%fileName, before='PFA File ', &
        & strip=.true. )
    else
      call output ( 'PFA file: ', advance='no' )
      call output ( trim(PFAFileDatum%nameString), advance='yes' )
    end if
    if ( myDetails > -4 ) then
      call output ( ' has groups:', advance='yes' )
      if ( associated(PFAFileDatum%ix) ) then
        do g = 1, ubound(PFAFileDatum%ix,1)
          if ( myDetails <= -3 ) then
            call blanks ( 1 )
            call getGroupName ( PFAData(PFAFileDatum%ix(g)), groupName )
            call output ( trim(groupName) )
            if ( myDetails == 2 ) then
              call output ( PFAData(PFAFileDatum%ix(g))%signalIndex, before=' Signal Index ' )
              call output ( PFAData(PFAFileDatum%ix(g))%channel, before=' Channel ' )
              call output ( PFAData(PFAFileDatum%ix(g))%theSignal%sideband, before=' Sideband ' )
            end if
            if ( associated(PFAData(PFAFileDatum%ix(g))%absorption) ) &
              call output ( ' has data' )
            call newLine
          else
            call dump ( PFAData(PFAFileDatum%ix(g)), myDetails, PFAFileDatum%ix(g) )
          end if
        end do
      else
        call output ( ' has no associated groups', advance='yes' )
      end if
    else
      call newLine
    end if

  end subroutine Dump_PFAFileDatum

  ! ------------------------------------------  Dump_PFAStructure  -----
  subroutine Dump_PFAStructure ( Details )
    ! Dump PFA data in Find_PFA order
    use Output_m, only: Output
    integer, intent(in) :: Details ! Passed through to Dump_PFADatum
    integer :: C, M, S, SX ! Channel, Molecule, Signal, Sideband (1..2)
    integer :: P           ! Index in PFAData

    if ( .not. associated(pfaData) ) then
      call output ( 'No PFA Database to dump', advance='yes' )
      return
    end if

    do m = first_molecule, last_molecule
      if ( associated(findPFA(m)%s) ) then
        do s = lbound(findPFA(m)%s,1), ubound(findPFA(m)%s,1)
          do sx = 1, 2
            if ( associated(findPFA(m)%s(s)%sb(sx)%c) ) then
              do c = lbound(findPFA(m)%s(s)%sb(sx)%c,1), ubound(findPFA(m)%s(s)%sb(sx)%c,1)
                p = findPFA(m)%s(s)%sb(sx)%c(c)
                if ( p /= 0 ) call dump ( PFAData(p), details, p )
              end do ! c
            end if
          end do ! sx
        end do ! s
      end if
    end do ! m
  end subroutine Dump_PFAStructure

  ! ------------------------------------------  Flush_PFADatabase  -----
  subroutine Flush_PFADatabase ( Signals, Molecules, Error )
    integer, pointer :: Signals(:)   ! All signals if disassociated
    integer, pointer :: Molecules(:) ! All molecules if disassociated
    integer, intent(out) :: Error    ! 0 => OK, else trouble

    integer :: I, J

    error = 0
    if ( .not. associated(PFAData) ) return ! Nothing to do

    ! For now, flush for all signals
    if ( associated(molecules) ) then
      ! Flush listed molecules for all signals
      do i = 1, ubound(PFAData,1)
        do j = 1, ubound(molecules,1)
          if ( PFAData(i)%molecule == molecules(j) ) then
            call destroy_PFADatum_Arrays ( PFAData(i) )
            exit
          end if
        end do ! j
      end do ! i
    else ! Flush for all molecules and signals
      do i = 1, ubound(PFAData,1)
        call destroy_PFADatum_Arrays ( PFAData(i) )
      end do
    end if

  end subroutine Flush_PFADatabase

  ! -----------------------------------------------  GetGroupName  -----
  subroutine GetGroupName ( PFADatum, GroupName )
  ! Get the group name for a PFA datum
    use Intrinsic, only: Lit_Indices
    use MLSSignals_m, only: GetSignalName
    use String_Table, only: Get_String, String_Length

    type(PFAData_t), intent(in) :: PFADatum
    character(len=*), intent(out) :: GroupName

    integer :: L

    call get_string ( lit_indices(PFADatum%molecule), groupName )
    l = string_length(lit_indices(PFADatum%molecule))
    if ( l > len(groupName) - 1 ) return ! with incomplete name -- sorry
    groupName(l+1:l+1) = '%'
    if ( PFADatum%signalIndex > 0 ) then
      call getSignalName ( PFADatum%signalIndex, groupName(l+2:), &
        & sideband=PFADatum%theSignal%sideband, channel=PFADatum%channel )
    else
      call get_string ( PFADatum%signal, groupName(l+2:) )
    end if

  end subroutine GetGroupName

  ! -----------------------------------------  HookTableToFindPFA  -----
  integer function HookTableToFindPFA ( F, G, PFADatum, Ix, Replace )
  ! Hook a PFA datum to the FindPFA structure, so that it can be found
  ! quickly given its molecule index, signal index, sideband and channel.
  ! Return index of one found one in the table already (if replace is present
  ! and true).  See Test_And_Fetch_PFA.
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, &
      & MLSMSG_Warning
    use MLSSignals_m, only: Signals
    use MoreMessage, only: MLSMessage

    integer, intent(in) :: F  ! Index in PFAFiles; zero if not from a file
                              ! specified in a PFAFile= parameter
    integer, intent(in) :: G  ! Index in PFAFiles(f)%ix; zero if not from
                              ! a file.  -1 if F /= 0 and HookTableToFindPFA is
                              ! to find where to hook it.
    type(PFAData_t), intent(in) :: PFADatum
    integer, intent(in) :: Ix ! Index of PFADatum in PFAData if not zero
    logical, intent(in), optional :: Replace ! Replace old one, default false

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: M, MyG, S, SB, C ! Molecule, MyG, Signal, Sideband, Channel
    logical :: MyReplace
    integer :: STAT
    character(len=molNameLen+1+maxSigLen+1) :: GroupName

    myReplace = .false.
    if ( present(replace) ) myReplace = replace

    m = PFADatum%molecule
    myG = g
    s = PFADatum%signalIndex
    sb = (PFADatum%theSignal%sideband + 3 ) / 2
    c = PFADatum%channel

    if ( m == 0 .or. s == 0 ) then
      hookTableToFindPFA = 0
      return
    end if

    if ( .not. associated(findPFA(m)%s) ) then
      allocate ( findPFA(m)%s(1:size(signals)), stat=stat )
      addr = 0
      if ( stat == 0 .and. size(signals) > 0 ) &
        & addr = transfer(c_loc(findPFA(m)%s(1)), addr)
      call test_allocate ( stat, moduleName, 'FindPFA(m)%s)', &
        & uBounds=[size(signals)], elementSize=storage_size(findPFA(m)%s)/8, &
        & address=addr )
    end if
    if ( .not. associated(findPFA(m)%s(s)%sb(sb)%c) ) &
      & call allocate_test ( findPFA(m)%s(s)%sb(sb)%c, &
        & ubound(signals(s)%frequencies,1),&
        & 'FindPFA(m)%s(s)%sb(sb)%c).', moduleName, &
        & lowBound=lbound(signals(s)%frequencies,1), fill=0 )
    hookTableToFindPFA = findPFA(m)%s(s)%sb(sb)%c(c)
    if ( hookTableToFindPFA /= 0 ) then
      if ( associated(PFAData(findPFA(m)%s(s)%sb(sb)%c(c))%absorption) ) then
        if ( myReplace ) then
            call destroy_PFADatum ( PFAData(hookTableToFindPFA) ) ! Don't leak
            call getGroupName ( PFAData(hookTableToFindPFA), groupName )
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Duplicate PFA data for '//trim(groupName) )
            findPFA(m)%s(s)%sb(sb)%c(c) = Ix
        else
          call getGroupName ( PFADatum, groupName )
          call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Duplicate PFA data for '//trim(groupName) )
        end if
      else
        findPFA(m)%s(s)%sb(sb)%c(c) = Ix
      end if
    else
      findPFA(m)%s(s)%sb(sb)%c(c) = Ix
    end if
    if ( ix == 0 ) return ! Just looking
    PFAData(findPFA(m)%s(s)%sb(sb)%c(c))%fileIndex = f
    if ( f == 0 ) then
      myG = 0
    else if ( myG < 0 ) then
      myG = 0
      if ( .not. associated(PFAFiles) ) then
        myG = -1
      else if ( .not. associated(PFAFiles(f)%ix) ) then
        myG = -1
      end if
      if ( myG < 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'PFA File or Group data structure not allocated' )
      do myG = ubound(PFAFiles(f)%ix,1), 1, -1
        if ( PFAData(PFAFiles(f)%ix(myG))%molecule == PFADatum%molecule .and. &
          &  PFAData(PFAFiles(f)%ix(myG))%signalIndex == PFADatum%signalIndex ) exit
      end do
    end if
    PFAData(findPFA(m)%s(s)%sb(sb)%c(c))%groupIndex = myG

  end function HookTableToFindPFA

  ! -------------------------------------------  Process_PFA_File  -----

  ! Process a PFA file name from the PFAFile parameter in GlobalSettings
  ! Open the HDF5 file, read the Signals and Molecules groups, allocate
  ! the groups, but don't open the groups.  Leave the file open.
  integer function Process_PFA_File_node ( PFAFileIndex, Where )

    integer, intent(in) :: PFAFileIndex ! Sub-rosa index for file name
    integer, intent(in) :: Where  ! tree node index, for error messages

    type(PFAFile_t) :: PFAFileDatum
    ! Open the file
    PFAFileDatum%fileName = PFAFileIndex
    Process_PFA_File_node = Process_PFA_File ( PFAFileDatum, where )

  end function Process_PFA_File_node

  integer function Process_PFA_File_name ( PFAFileName, Where )

    character(len=*), intent(in) :: PFAFileName ! Actual path/name
    integer, intent(in) :: Where  ! tree node index, for error messages

    type(PFAFile_t) :: PFAFileDatum
    ! Open the file
    PFAFileDatum%fileName = 0
    PFAFileDatum%NameString = PFAFileName
    Process_PFA_File_name = Process_PFA_File ( PFAFileDatum, where )

  end function Process_PFA_File_name

  integer function Process_PFA_File_datum ( PFAFileDatum, Where )

    use Allocate_Deallocate, only: Allocate_Test
    use MLSHdf5, only: Loadptrfromhdf5ds
    use MLSMessagemodule, only: MLSmessage, MLSmsg_Error, &
      & MLSMsg_Severity_To_Quit, MLSmsg_Warning
    use MLSSignals_M, only: Signals
    use Moremessage, only: MLSmessage
    use Moretree, only: Getlitindexfromstring, Getstringindexfromstring
    use Parse_Signal_M, only: Parse_Signal
    use Printit_M, only: Set_Config
    ! Hdf5 Intentionally Last To Avoid Long Lf95 Compiles
    use Hdf5, only: H5gclose_F, H5gopen_F
    use Toggles, only: Gen, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Where_At => Where

    type(PFAFile_t) :: PFAFileDatum
    integer, intent(in) :: Where  ! tree node index, for error messages and tracing

    logical, pointer :: Channels(:)
    logical :: Error
    integer :: F, G                    ! Subscripts, loop inductors
    integer(hid_t) :: GroupID          ! From HDF5 open group
    integer :: IPFA                    ! Index in PFA database
    integer :: Me = -1                 ! String index for trace
    character(len=molNameLen), pointer :: MyMolecules(:)
    character(len=maxSigLen), pointer :: MySignalStrings(:) ! From the HDF5
    integer :: MLSMSG_Severity_to_quit_old
    integer :: SB                      ! Sideband, from the signal
    integer, pointer :: SignalIndices(:)
    integer :: STAT                    ! From HDF5 open or read, or allocate
    logical :: Trouble                 ! If molecule or signal doesn't work

    call trace_begin ( me, "Process_PFA_File_datum", where, cond=toggle(gen) )
    ! Open the file

    call OpenPFAFile ( PFAFileDatum, where )

    ! Open the Index group and read the Molecules and Signals data sets
    call h5gOpen_f ( PFAFileDatum%HDF5_fileID, 'Index', groupID, stat )
    if ( stat /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open HDF5 PFA Index group in %s at %w', &
        & (/PFAFileDatum%fileName /), where=where_at(where) )
    nullify ( myMolecules, mySignalStrings )
    call loadPtrFromHDF5DS ( groupID, 'Molecules', myMolecules )
    call loadPtrFromHDF5DS ( groupID, 'Signals', mySignalStrings )
    call h5gClose_f ( groupID, stat )
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close HDF5 PFA index group.' )

    ! Allocate the groups component
    call allocate_test ( PFAFileDatum%ix, size(myMolecules), 'PFAFileDatum%ix', &
      & ModuleName )

    call create_or_expand_PFADatabase ( size(myMolecules), iPFA )

    ! Fill each group's fields, except for the data and group ID.
    ! PFAData(PFAFileDatum%ix(g))%open is default initialized to be .false.
    nullify ( channels, signalIndices )
    do g = 1, size(myMolecules)
      trouble = .false.
      iPFA = iPFA + 1
      PFAData(iPFA)%molecule = getLitIndexFromString ( trim(myMolecules(g)) )
      trouble = PFAData(iPFA)%molecule < 1
      if ( trouble ) then
        if ( PFAFileDatum%fileName /= 0 ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & "The string " // trim(myMolecules(g)) // &
            & " in %s at %w is not a molecule.%n" // &
            & "Molecule list probably changed after PFA tables generated.", &
            & (/PFAFileDatum%fileName /), where=where_at(where) )
        else
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & "The string " // trim(myMolecules(g)) // &
            & " in " // trim(PFAFileDatum%namestring) // &
            & " at %w is not a molecule.%n" // &
            & "Molecule list probably changed after PFA tables generated.", &
            & where=where_at(where) )
        end if
      end if
      PFAData(iPFA)%signal = &
        & getStringIndexFromString ( trim(mySignalStrings(g)) )
      trouble = trouble .or. PFAData(iPFA)%signal < 1
      if ( .not. trouble) &
        & call parse_signal ( trim(mySignalStrings(g)), signalIndices, &
        & channels=channels, sideband=sb )
      if ( .not. associated(signalIndices) ) then
        call MLSMessage ( MLSMSG_Warning, &
        & moduleName, 'Unable to parse signal ' // trim(mySignalStrings(g)) )
        trouble = .true.
      end if
      if ( .not. trouble ) then
        PFAData(iPFA)%channel = lbound(channels,1)
        PFAData(iPFA)%signalIndex = signalIndices(1)
        PFAData(iPFA)%theSignal = signals(signalIndices(1))
        PFAData(iPFA)%theSignal%channels => channels
        PFAData(iPFA)%theSignal%sideband = sb
        PFAFileDatum%ix(g) = iPFA
      else
        PFAData(iPFA)%molecule = 0
        PFAData(iPFA)%signal = 0
        PFAData(iPFA)%channel = 0
        PFAData(iPFA)%signalIndex = 0
        nullify ( PFAData(iPFA)%theSignal%channels )
        PFAData(iPFA)%theSignal%sideband = 0
        PFAFileDatum%ix(g) = 0
      end if
      nullify ( channels ) ! so as not to clobber the one we just stored
    end do

    ! Add the PFAFileDatum to the database
    f = addPFAFileDatumToDatabase ( PFAFiles, PFAFileDatum )

    ! Clean up
    call deallocate_test ( channels, 'Channels', moduleName )
    call deallocate_test ( signalIndices, 'SignalIndices', moduleName )
    call deallocate_test ( myMolecules, 'MyMolecules', moduleName )
    call deallocate_test ( mySignalStrings, 'mySignalStrings', moduleName )

    ! Now hook the tables to the FindPFA structure
    ! Report all errors instead of quitting on the first one
    error = .false.
    MLSMSG_Severity_to_quit_old = MLSMSG_Severity_to_quit
    MLSMSG_Severity_to_quit = MLSMSG_Error + 1
    call set_config ( severity_to_quit = MLSMSG_Severity_to_quit )
    do g = 1, ubound(PFAFileDatum%ix,1)
      error = error .or. &
        & hookTableToFindPFA ( f, g, PFAData(PFAFileDatum%ix(g)), PFAFileDatum%ix(g) ) /= 0
    end do
    MLSMSG_Severity_to_quit = MLSMSG_Severity_to_quit_old ! MLSMSG_Error
    call set_config ( severity_to_quit = MLSMSG_Severity_to_quit )
    if ( error ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Errors while hooking up and finding PFA tables' )

    Process_PFA_File_datum = f
    call trace_end ( "Process_PFA_File_datum", cond=toggle(gen) )

  end function Process_PFA_File_datum

  ! ................................  AddPFAFileDatumToDatabase  .....
  integer function AddPFAFileDatumToDatabase ( DATABASE, ITEM )

  ! This routine adds a PFA File Datum to a database of PFA File Data,
  ! creating the database if necessary.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (PFAFile_T), dimension(:), pointer :: DATABASE
    type (PFAFile_T), intent(in) :: ITEM

    ! Local variables
    type (PFAFile_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddPFAFileDatumToDatabase = newSize
  end function AddPFAFileDatumToDatabase

  ! -------------------------------------------  Read_PFADatabase  -----
  subroutine Read_PFADatabase ( FileNameIndex, FileTypeIndex, TheMolecules, &
    & TheSignals, Where )
  ! Read the PFA data from FileName.  If both TheMolecules and
  ! TheSignals have zero size, all PFA data are read from FileName.
  ! If TheMolecules has zero size but TheSignals does not, the PFA
  ! data for all molecules for each specified signal are read from
  ! FileName.  If TheSignals has zero size but TheMolecules does
  ! not, the PFA data for all signals for each specified molecule are read
  ! from FileName. Otherwise the PFA data for the Cartesian product of
  ! TheMolecules and TheSignals are read from FileName.

    use Intrinsic, only: Lit_Indices
    use MLSmessagemodule, only: MLSmessage, MLSmsg_Error
    use Moremessage, only: MLSmessage
    use Moretree, only: Getstringindexfromstring
    use Parse_Signal_M, only: Get_Individual_Signals
    use Tree, only: Where_At => Where
    ! Hdf5 Intentionally Last To Avoid Long Lf95 Compiles
    use Hdf5, only: H5fclose_F, H5gclose_F

    integer, intent(in) :: FileNameIndex
    integer, intent(in) :: FileTypeIndex ! HDF5 is all we can do
    integer, intent(in) :: TheMolecules(:), TheSignals(:)
    integer, intent(in) :: Where ! tree index, for error messages

    integer, pointer :: AllSignals(:) ! String table indices for AllSignalStrings
    character(len=maxSigLen), pointer :: AllSignalStrings(:) ! from expanding a signal
    integer :: F, G ! Index in PFAFiles, PFAFiles(f)%ix
    integer, pointer :: Groups(:) ! Indices in MyMolecules, MySignals
    integer :: I, IOSTAT, J, K, L
    integer :: NGroups, NPFA

    f = 0
    if ( associated(PFAFiles) ) then
      do f = ubound(PFAFiles,1), 1, -1
        if ( PFAFiles(f)%fileName == fileNameIndex ) exit
      end do
    end if
    if ( f == 0 ) f = process_PFA_File ( fileNameIndex, Where )

    nullify ( allSignalStrings, allSignals )
    ! Expand theSignals to allSignalStrings
    call get_individual_signals ( allSignalStrings, theSignals )

    call allocate_test ( allSignals, size(allSignalStrings), 'AllSignals', &
      & moduleName )
    do i = 1, size(allSignalStrings)
      allSignals(i) = getStringIndexFromString ( allSignalStrings(i) )
    end do

!   if ( fileType == 'HDF5' ) then
      nGroups = ubound(PFAFiles(f)%ix,1)
      nullify ( groups )
      call allocate_test ( groups, nGroups, 'Groups', moduleName )

      ! Decide what to read
      l = max(size(theMolecules) * size(allSignals), &
        &     size(theMolecules), size(allSignals) )
      nPFA = 0
      if ( l == 0 ) then ! Read everything
        nPFA = nGroups
        call allocate_test ( groups, nPFA, 'Groups', moduleName )
        do i = 1, nPFA
          groups(i) = i
        end do
      else if ( size(theMolecules) == 0 ) then
        do i = 1, size(allSignals)
          do k = 1, nGroups
            if ( allSignals(i) == PFAData(PFAFiles(f)%ix(k))%signal ) then
              nPFA = nPFA + 1
              groups(nPFA) = k
            end if
          end do
        end do
      else if ( size(allSignals) == 0 ) then
        do i = 1, size(theMolecules)
          do k = 1, nGroups
            if ( theMolecules(i) == PFAData(PFAFiles(f)%ix(k))%molecule ) then
              nPFA = nPFA + 1
              groups(nPFA) = k
            end if
          end do
        end do
      else ! Both theSignals and theMolecules are present -- get as much
           ! of the Cartesian product of them as possible
        do i = 1, size(theMolecules)
          do j = 1, size(allSignals)
            do k = 1, nGroups
              if ( theMolecules(i) == PFAData(PFAFiles(f)%ix(k))%molecule .and. &
                & allSignals(j) == PFAData(PFAFiles(f)%ix(k))%signal ) then
                nPFA = nPFA + 1
                groups(nPFA) = k
              end if
            end do
          end do
        end do
      end if
      call deallocate_test ( allSignals, 'AllSignals', moduleName )
      call deallocate_test ( allSignalStrings, 'AllSignalStrings', moduleName )
      if ( nPFA == 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'No PFA Data to read from %s at %w.', &
        & (/ fileNameIndex /), where=where_at(where) )

      ! Read the groups
      call openPFAFile ( PFAFiles(f), where )
      do i = 1, nPFA
        g = groups(i)
        call openPFAGroup ( PFAFiles(f), PFAData(PFAFiles(f)%ix(g)), where )
        call read_PFADatum_H5 ( PFAData(PFAFiles(f)%ix(g))%HDF5_groupID, &
          & PFAData(PFAFiles(f)%ix(g)), .true. )
        call h5gClose_f ( PFAData(PFAFiles(f)%ix(g))%HDF5_groupID, iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close HDF5 PFA group %s%%s in %s at %w.', &
          & (/lit_indices(PFAData(PFAFiles(f)%ix(g))%molecule),&
          &   PFAData(PFAFiles(f)%ix(g))%signal,PFAFiles(f)%fileName/), &
          &   where=where_at(where) )
        PFAData(PFAFiles(f)%ix(g))%open = .false.
      end do ! i = 1, nPFA
      call H5FClose_F ( PFAFiles(f)%HDF5_FileID, iostat )
      PFAFiles(f)%open = .false.
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close PFA data file %s at %w.', &
        & (/PFAFiles(f)%fileName/),where=where_at(where) )
      call deallocate_test ( groups, 'Groups', moduleName )
!   else
!     Nothing -- Read_PFAData only allows HDF5 file type
!   end if

  end subroutine Read_PFADatabase

  ! -----------------------------------------  Test_And_Fetch_PFA  -----
  integer function Test_And_Fetch_PFA ( Molecule, Signal, Sideband, Channel, Derivs) &
    &                             result ( PFADatumIx )

    ! Test whether the file and group of PFA data are known.  If not, emit an
    ! error message and crash, or return NULL, depending upon the Crash
    ! parameter here.  See HookTableToFindPFA.
    ! If it is known, test whether the datum is filled.  If not, read the
    ! necessary data sets and add any grids that aren't already in VGrids.
    ! Return a pointer associated with the datum.

    use Intrinsic, only: Lit_Indices
    use MLSHDF5, only: LoadptrfromHDF5ds
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Error
    use Moremessage, only: MLSMessage

    integer, intent(in) :: Molecule ! first_molecule : last_molecule
    integer, intent(in) :: Signal   ! 1 : size(signals)
    integer, intent(in) :: Sideband ! -1 = LSB, +1 = USB
    integer, intent(in) :: Channel  ! depends on signal
    logical, intent(in) :: Derivs   ! "Read dAbs dwc and dAbs dnc"

    logical, parameter :: Crash = .true.

    integer :: SB                   ! Sideband Subscript 1..2
    logical, parameter :: debug = .false.

    PFADatumIx = 0

    sb = ( sideband + 3 ) / 2 ! -1..+1 => 1..2
    if ( debug ) print *, 'Molecule, Signal, Sideband, Channel, Derivs ', &
      & Molecule, Signal, Sideband, Channel, Derivs
    if ( associated(findPFA(molecule)%s) ) then
      if ( associated(findPFA(molecule)%s(signal)%sb(sb)%c ) ) then
        PFADatumIx = findPFA(molecule)%s(signal)%sb(sb)%c(channel)
        if ( PFADatumIx /= 0 ) then
          if ( debug ) print *, 'PFADatumIx ', PFADatumIx
          call readDS ( 'absorption', PFAData(PFADatumIx)%absorption )
          call readDS ( 'dAbsDnu', PFAData(PFADatumIx)%dAbsDnu )
          if ( derivs ) then
            call readDS ( 'dAbsDwc', PFAData(PFADatumIx)%dAbsDwc )
            call readDS ( 'dAbsDnc', PFAData(PFADatumIx)%dAbsDnc )
          end if
          return
        end if
      end if
    end if

    if ( crash ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'No PFA data for signal %g, channel %d, sideband %d and molecule %s.', &
        & (/ signal, channel, sideband, lit_indices(molecule)/) )

  contains
    subroutine ReadDS ( DSName, DSData )
      character(len=*), intent(in) :: DSName
      real(rk), pointer :: DSData(:,:)
      integer :: F, G ! Subscripts
      if ( associated(dsData) ) return
      f = PFAData(PFADatumIx)%fileIndex
      if ( f == 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'PFA Data was not read from a file and therefore cannot be recovered.' )
      g = PFAData(PFADatumIx)%groupIndex
      if ( .not. associated(PFAData(PFADatumIx)%theSignal%channels) .or. &
        &  .not. associated(PFAData(PFADatumIx)%vGrid%surfs) ) then ! nothing there
        call openPFAGroup ( PFAFiles(f), PFAData(PFAFiles(f)%ix(g)), 0 )
        call read_PFADatum_H5 ( PFAData(PFAFiles(f)%ix(g))%HDF5_GroupID, &
          & PFAData(PFADatumIx), derivs )
      end if
      call openPFAGroup ( PFAFiles(f), PFAData(PFAFiles(f)%ix(g)), 0 )
      call loadPtrFromHDF5DS ( PFAData(PFAFiles(f)%ix(g))%HDF5_GroupID, dsName, dsData )
    end subroutine ReadDS
  end function Test_And_Fetch_PFA

  ! ------------------------------------------  Write_PFADatabase  -----
  subroutine Write_PFADatabase ( FileName, FileType )
    use Intrinsic, only: Lit_Indices
    use MLSHDF5, only: SaveasHDF5ds
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_M, only: Maxsiglen
    ! use MLSStrings, only: Capitalize
    use String_Table, only: Get_String
    ! HDF5 Intentionally Last To Avoid Long Lf95 Compiles
    use HDF5, only: H5fcreate_F, H5fclose_F, &
      & H5f_Acc_Trunc_F, H5gclose_F, H5gcreate_F

    character(len=*), intent(in) :: FileName, FileType

    integer :: FileID, GroupID
    integer :: I, IOSTAT
    character(len=molNameLen) :: Molecules(ubound(pfaData,1))
    character(len=maxSigLen) :: SignalText(ubound(pfaData,1))

    if ( .not. associated(pfaData) ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'No PFA Data to write' )
!   if ( capitalize(fileType) == 'HDF5' ) then ! open HDF5 file here
      call H5FCreate_F ( trim(fileName), H5F_ACC_TRUNC_F, fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to open hdf5 PFA file ' // trim(fileName) // ' for output.' )
      ! Make the Index group
      call h5gCreate_f ( fileID, 'Index', groupID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to create hdf5 Index group in ' // trim(fileName) // '.' )
      do i = 1, ubound(pfaData,1)
        call get_string ( lit_indices(pfaData(i)%molecule), molecules(i) )
        call get_string ( pfaData(i)%signal, signalText(i) )
      end do
      call saveAsHDF5DS ( groupID, 'Molecules', molecules )
      call saveAsHDF5DS ( groupID, 'Signals', signalText )
      call h5gClose_F ( groupID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to close hdf5 Index group in ' // trim(fileName) // '.' )
      ! Write the PFA Data tables, one per group
      do i = 1, ubound(pfaData,1)
        if ( associated(pfaData(i)%absorption) ) &
          & call write_PFADatum ( pfaData(i), FileName, FileType, lun=fileID )
      end do
      call H5FClose_F ( fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to close hdf5 PFA file ' // trim(fileName) // '.' )
!   else
!     Nothing -- file type is checked in Write_PFAData
!   end if

  end subroutine Write_PFADatabase

  ! ---------------------------------------------  Write_PDADatum  -----
  subroutine Write_PFADatum ( PFADatum, FileName, FileType, &
    & Lun )

    ! Write the PFADatum on FileName using the format given by FileType

    use MLSHDF5, only: MakeHDF5attribute, SaveasHDF5ds, WritelitindexasHDF5attribute
    use MLSMessagemodule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_M, only: Maxsiglen, Getnameofsignal
    use MLSStringlists, only: Switchdetail
    use Output_M, only: Output
    use String_Table, only: Get_String, String_Length
    use Toggles, only: Switches
    use HDF5, only: H5fcreate_F, H5fclose_F, & ! HDF5 use Intentionally Last
      & H5f_Acc_Trunc_F, H5gclose_F, H5gcreate_F

    type(PFAData_t), intent(in) :: PFADatum
    character(len=*), intent(in) :: FileName, FileType
    integer, intent(in), optional :: Lun ! Don't open a new file if present

    character(len=1023) :: Attrib
    integer :: GroupID
    character(len=molNameLen+1+maxSigLen+1) :: GroupName
    integer :: IOSTAT, MyLun
    character(len=maxSigLen) :: SignalText

    if ( present(lun) ) myLun = lun
!   if ( capitalize(fileType) == 'HDF5' ) then
      if ( .not. present(lun) ) then
        ! Open the HDF file
        call H5FCreate_F ( trim(fileName), H5F_ACC_TRUNC_F, myLun, &
          & iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to open hdf5 PFA file ' // trim(fileName) // ' for output.' )
      end if

      ! Create a group named <molecule>%<signal>
      call getGroupName ( pfaDatum, groupName )
      call h5gCreate_f ( myLun, trim(groupName), groupID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create group for PFA table ' // trim(groupName) )
      if ( switchDetail(switches,'pfaw') > -1 ) then
        call output ( 'Created PFA group ' )
        call output ( trim(groupName) )
        call output ( ' in HDF5 file ' )
        call output ( trim(fileName), advance='yes' )
      end if

      ! Fill the group
      if ( pfaDatum%name /= 0 ) then
        call get_string ( pfaDatum%name, attrib )
        call MakeHDF5Attribute ( groupID, 'name', attrib(:string_length(pfaDatum%name)) )
      end if
      if ( pfaDatum%filterFile /= 0 ) then
        call get_string ( pfaDatum%filterFile, attrib )
        call MakeHDF5Attribute ( groupID, 'filterFile', &
          & attrib(:string_length(pfaDatum%filterFile)) )
      end if
      if ( pfaDatum%spectroscopyFile /= 0 ) then
        call get_string ( pfaDatum%spectroscopyFile, attrib )
        call MakeHDF5Attribute ( groupID, 'spectroscopyFile', &
          & attrib(:string_length(pfaDatum%spectroscopyFile)) )
      end if

      call WriteLitIndexAsHDF5Attribute ( groupID, 'molecule', pfaDatum%molecule )
      call GetNameOfSignal ( pfaDatum%theSignal, signalText )
      call MakeHDF5Attribute ( groupID, 'signal', signalText )
      call MakeHDF5Attribute ( groupID, 'sideband', pfaDatum%theSignal%sideband )
      call MakeHDF5Attribute ( groupID, 'vel_rel', pfaDatum%vel_rel )
      call MakeHDF5Attribute ( groupID, 'whichLines', pfaDatum%whichLines )
      call MakeHDF5Attribute ( groupID, 'nTemps', pfaDatum%tGrid%noSurfs )
      call MakeHDF5Attribute ( groupID, 'tStart', pfaDatum%tGrid%surfs(1,1) )
      call MakeHDF5Attribute ( groupID, 'tStep', &
        & pfaDatum%tGrid%surfs(2,1)-pfaDatum%tGrid%surfs(1,1) )
      call MakeHDF5Attribute ( groupID, 'nPress', pfaDatum%vGrid%noSurfs )
      call MakeHDF5Attribute ( groupID, 'vStart', pfaDatum%vGrid%surfs(1,1) )
      call MakeHDF5Attribute ( groupID, 'vStep', &
        & pfaDatum%vGrid%surfs(2,1)-pfaDatum%vGrid%surfs(1,1) )
      call SaveAsHDF5DS ( groupID, 'absorption', pfaDatum%absorption )
      call SaveAsHDF5DS ( groupID, 'dAbsDwc', pfaDatum%dAbsDwc )
      call SaveAsHDF5DS ( groupID, 'dAbsDnc', pfaDatum%dAbsDnc )
      call SaveAsHDF5DS ( groupID, 'dAbsDnu', pfaDatum%dAbsDnu )

      ! Close the group
      call h5gClose_f ( groupID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close qroup for PFA table ' // trim(groupName) )

      if ( .not. present(lun) ) then
        ! Close the HDF file
        call H5FClose_F ( myLun, iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close hdf5 PFA file ' // trim(fileName) // '.' )
      end if
!   else
!     Nothing -- file type is checked by caller
!   end if

  end subroutine Write_PFADatum

! =====     Private Procedures     =====================================

  ! -------------------------------  Create_or_Expand_PFADatabase  -----
  subroutine Create_or_Expand_PFADatabase ( ToAdd, PrevSize )
    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    integer, intent(in) :: ToAdd
    integer, intent(out) :: PrevSize
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: S, Stat
    type(PFAData_t), pointer :: TempPFAData(:) => NULL()
    if ( associated(PFAData) .and. toAdd > 0 ) then
      tempPFAData => PFAData
      prevSize = ubound(tempPFAData,1)
      allocate ( PFAData(0:prevSize+toAdd), stat=stat )
      addr = 0
      if ( stat == 0 ) then
        if ( size(PFAData) > 0 ) addr = transfer(c_loc(PFAData(0)), addr)
      end if
      call test_allocate ( stat, moduleName, 'PFAData', lBounds=[0], &
        & uBounds=[prevSize+toAdd], elementSize=storage_size(PFAData) / 8, &
        & address=addr )
      pfaData(:prevSize) = tempPFAData
      s = size(tempPFAData) * storage_size(tempPFAData) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(tempPFAData(lbound(tempPFAData,1))), addr)
      deallocate ( tempPFAData, stat=stat )
      call test_deallocate ( stat, moduleName, 'TempPFAData', s, address=addr )
    else
      prevSize = 0
      if ( toAdd > 0 ) then
        allocate ( PFAData(0:toAdd), stat=stat )
        if ( stat == 0 .and. toAdd >= 0 ) addr = transfer(c_loc(PFAData(0)), addr)
        call test_allocate ( stat, moduleName, 'PFAData', lBounds=[0], &
          & uBounds=[toAdd], elementSize=storage_size(PFAData) / 8, address=addr )
      end if
    end if
  end subroutine Create_or_Expand_PFADatabase

  ! ------------------------------------------------  OpenPFAFile  -----
  subroutine OpenPFAFile ( PFAFileDatum, Where )
    use MLSMessageModule, only: MLSMSG_Error, MLSMSG_FileOpen
    use MoreMessage, only: MLSMessage
    use String_Table, only: Get_String
    use Tree, only: Where_At => Where
    use HDF5, only: H5F_ACC_RDONLY_F, H5FOpen_F
    type(PFAFile_t), intent(inout) :: PFAFileDatum
    integer, intent(in) :: Where
    character(len=1023) :: PFAFileName
    integer :: Stat ! from Open
    if ( PFAFileDatum%open ) return
    if ( PFAFileDatum%fileName > 0 ) then
      call get_string ( PFAFileDatum%fileName, PFAFileName, strip=.true. )
      PFAFileDatum%nameString = PFAFileName
    else
      PFAFileName = PFAFileDatum%nameString
    end if
    call h5fopen_f ( trim(PFAFileName), H5F_ACC_RDONLY_F, &
      & PFAFileDatum%HDF5_fileID, stat )
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Fileopen // ' %s (HDF5 PFA data) at %w.', &
          & (/PFAFileDatum%fileName/),where=where_at(where) )
    PFAFileDatum%open = .true.
  end subroutine OpenPFAFile

  ! -----------------------------------------------  OpenPFAGroup  -----
  subroutine OpenPFAGroup ( PFAFileDatum, PFADatum, Where )
    use Intrinsic, only: Lit_Indices
    use String_Table, only: Get_String, String_Length
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MoreMessage, only: MLSMessage
    use Tree, only: Where_At => Where
    use HDF5, only: H5GOpen_f
    type(PFAFile_t), intent(inout) :: PFAFileDatum
    type(PFAData_t), intent(inout) :: PFADatum
    integer, intent(in) :: Where ! For error messages
    integer :: IOSTAT, L
    character(len=molNameLen+maxSigLen+1) :: TheGroup
    call openPFAFile ( PFAFileDatum, where )
    if ( .not. PFADatum%open ) then
      call get_string ( lit_indices(PFADatum%molecule), &
        & theGroup )
      l = string_length(lit_indices(PFADatum%molecule))
      theGroup(l+1:l+1) = '%'
      call get_string ( PFADatum%signal, theGroup(l+2:) )
      l = l + string_length(PFADatum%signal) + 2
      call h5gOpen_f ( PFAFileDatum%HDF5_FileID, theGroup(:l), &
        & PFADatum%HDF5_groupID, iostat )
      if ( iostat /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to open HDF5 PFA group ' // theGroup(:l) // ' in %s at %w.', &
          & (/PFAFileDatum%fileName/),where=where_at(where) )
      PFADatum%open = .true.
    end if
  end subroutine OpenPFAGroup

  ! -------------------------------------------  Read_PFADatum_H5  -----
  subroutine Read_PFADatum_H5 ( GroupID, PFADatum, Derivs, Data )
  ! Read the PFA Datum from HDF5 group GroupId into PFADatum.
  ! Read dAbsDwc and dAbsDnc if and only if Derivs is true.  Add any
  ! grids that aren't already in VGrids.

    use Intrinsic, only: L_Theta, L_Zeta, Lit_Indices
    use MLSHDF5, only: GetHDF5attribute, IsHDF5attributepresent, LoadptrfromHDF5ds, &
      & ReadlitindexfromHDF5attr
    use MLSMessagemodule, only: MLSMSG_Error
    use MLSSignals_M, only: Getsignalname, Maxsiglen, Signals
    use MoreMessage, only: MLSMessage
    use Moretree, only: Getstringindexfromstring
    use Parse_Signal_M, only: Parse_Signal
    ! use VgridsDatabase, only: Addvgridifnecessary, Vgrids ! Rs=kind For Surfs
    use VgridsDatabase, only: Rs ! Rs=kind For Surfs

    integer(hid_t), intent(in) :: GroupId
    type(PFAData_T), target :: PFADatum
    logical, intent(in) :: Derivs
    logical, intent(in), optional :: Data ! "Read the data -- default true"

    logical, pointer :: Channels(:) ! output from Parse_Signal
    integer :: J, K
    character(1023) :: Line ! Text, e.g. filter file name
    logical :: MyData
    integer, pointer :: SignalIndices(:) ! output from Parse_Signal
    character(len=maxSigLen) :: SignalText
    real(rs) :: SurfStep ! for temperature and pressure grids
    type(vGrid_t) :: TGrid, VGrid
!     integer :: GRIDINDEX
    integer :: SB  ! Sideband from input

    myData = .true.
    if ( present(data) ) myData = data

    nullify ( channels, signalIndices )
    if ( isHDF5AttributePresent(groupID, 'filterFile') ) then
      call getHDF5Attribute ( groupID, 'filterFile', line )
      PFADatum%filterFile = getStringIndexFromString ( trim(line), .true. )
    else
      PFADatum%filterFile = 0
    end if
    if ( isHDF5AttributePresent(groupID, 'spectroscopyFile') ) then
      call getHDF5Attribute ( groupID, 'spectroscopyFile', line )
      PFADatum%spectroscopyFile = getStringIndexFromString ( trim(line), .true. )
    else
      PFADatum%spectroscopyFile = 0
    end if
    call ReadLitIndexFromHDF5Attr ( groupID, 'molecule', k )
    if ( k < first_molecule .or. k > last_molecule ) call MLSMessage ( &
      & MLSMSG_Error, moduleName, 'The string %s is not a molecule name.', &
      & (/ lit_indices(k) /) )
    PFADatum%molecule = k
    call getHDF5Attribute ( groupID, 'signal', signalText )
    call getHDF5Attribute ( groupID, 'sideband', sb )
    call parse_signal ( signalText, signalIndices, channels=channels )
    do j = lbound(channels,1), ubound(channels,1)
      if ( channels(j) ) exit
    end do
    PFADatum%channel = j
    PFADatum%signalIndex = signalIndices(1)
    call getSignalName ( PFADatum%signalIndex, signalText, sideband=sb, &
      & channel=PFADatum%channel )
    PFADatum%signal = getStringIndexFromString ( trim(signalText), .true. )
    PFADatum%theSignal = signals(PFADatum%signalIndex)
    PFADatum%theSignal%channels => channels
    PFADatum%theSignal%sideband = sb
    nullify ( channels ) ! so as not to clobber PFADatum%theSignal%channels
      ! in next iteration of the loop
    call getHDF5Attribute ( groupID, 'vel_rel', PFADatum%vel_rel )
    PFADatum%whichLines = 0
    if ( isHDF5AttributePresent(groupID, 'whichLines') ) &
      & call getHDF5Attribute ( groupID, 'whichLines', PFADatum%whichLines )
    tGrid%name = 0
    nullify ( tGrid%surfs, vGrid%surfs )
    tGrid%verticalCoordinate = l_theta
    call getHDF5Attribute ( groupID, 'nTemps', tGrid%noSurfs )
    call allocate_test ( tGrid%surfs, tGrid%noSurfs, 1, &
      & 'tGrid%surfs', moduleName )
    call getHDF5Attribute ( groupID, 'tStart', tGrid%surfs(1,1) )
    call getHDF5Attribute ( groupID, 'tStep', surfStep )
    do j = 2, tGrid%noSurfs
      tGrid%surfs(j,1) = tGrid%surfs(j-1,1) + surfStep
    end do

    ! These lines seemed to have resulted in the VGrids components of
    ! FwdmodelConfs becoming clobbered in sids tests 2005 Aug 9-10
    ! Similarly below
!    gridIndex = addVGridIfNecessary(tGrid, dontDestroy=.true.) ! Needed in case vGrids not yet associated
!     PFADatum%tGrid = vGrids(gridIndex)
     PFADatum%tGrid = tGrid
    vGrid%name = 0
    vGrid%verticalCoordinate = l_zeta
    call getHDF5Attribute ( groupID, 'nPress', vGrid%noSurfs )
    call allocate_test ( vGrid%surfs, vGrid%noSurfs, 1, &
      & 'vGrid%surfs', moduleName )
    call getHDF5Attribute ( groupID, 'vStart', vGrid%surfs(1,1) )
    call getHDF5Attribute ( groupID, 'vStep', surfStep )
    do j = 2, vGrid%noSurfs
      vGrid%surfs(j,1) = vGrid%surfs(j-1,1) + surfStep
    end do
!     gridIndex = addVGridIfNecessary(vGrid, &
!       &                       relErr=vGrid%noSurfs*0.2_rs*epsilon(1.0_rs), &
!       & dontDestroy=.true.)
!     PFADatum%vGrid = vGrids(gridIndex)
     PFADatum%vGrid = vGrid

    if ( myData ) then
      call loadPtrFromHDF5DS ( groupID, 'absorption', PFADatum%absorption )
      call loadPtrFromHDF5DS ( groupID, 'dAbsDnu', PFADatum%dAbsDnu )
      if ( derivs ) then
        call loadPtrFromHDF5DS ( groupID, 'dAbsDwc', PFADatum%dAbsDwc )
        call loadPtrFromHDF5DS ( groupID, 'dAbsDnc', PFADatum%dAbsDnc )
      end if
    end if

    call deallocate_test ( signalIndices, 'SignalIndices', moduleName )

  end subroutine Read_PFADatum_H5

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module PFADataBase_m

! $Log$
! Revision 2.55  2018/04/11 22:25:23  vsnyder
! Remove USE for unused names
!
! Revision 2.54  2017/12/07 02:41:56  vsnyder
! Remove unreferenced use name
!
! Revision 2.53  2017/02/10 01:08:20  pwagner
! Report all PFA errors before quitting
!
! Revision 2.52  2015/03/28 02:00:00  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.51  2014/09/05 20:50:53  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.50  2014/08/06 23:25:48  vsnyder
! Remove USE for BYTES, which is not reference.  Comment out declaration of
! GRIDINDEX, which is only referenced in commented-out code.
!
! Revision 2.49  2014/07/18 23:15:03  pwagner
! Aimed for consistency in names passed to allocate_test
!
! Revision 2.48  2014/01/11 01:28:53  vsnyder
! Decruftification
!
! Revision 2.47  2013/09/24 23:28:17  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.46  2013/08/30 03:56:23  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.45  2013/08/23 02:51:26  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 2.44  2013/06/13 21:05:18  vsnyder
! More cruft removal
!
! Revision 2.43  2013/06/12 02:20:59  vsnyder
! Cruft removal
!
! Revision 2.42  2011/07/23 00:17:32  vsnyder
! More robust response to not finding a molecule string or signal string,
! which might happen if the molecules list or signals database is changed
! after the PFA tables are made.
!
! Revision 2.41  2011/05/09 17:51:35  pwagner
! Converted to using switchDetail
!
! Revision 2.40  2010/01/23 01:05:01  vsnyder
! Cannonball polishing
!
! Revision 2.39  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.38  2008/09/04 19:59:32  vsnyder
! Add PRINT statement in not_used_here
!
! Revision 2.37  2008/05/22 01:07:22  vsnyder
! Cannonball polishing
!
! Revision 2.36  2007/10/03 23:58:26  vsnyder
! Add 'where' for tracing
!
! Revision 2.35  2006/04/22 01:30:46  vsnyder
! Get channel number into signal read by read_PFADatum_H5
!
! Revision 2.34  2006/04/21 22:23:27  vsnyder
! Allow to flush specific molecules, other stuff for updating PFA
!
! Revision 2.32  2006/01/26 03:06:52  vsnyder
! Accumulate all errors before crashing
!
! Revision 2.31  2005/08/11 00:17:50  pwagner
! Comment out  addVGridIfNecessary in Read_PFADatum_H5
!
! Revision 2.30  2005/06/03 22:55:04  vsnyder
! Make the 'details' argument for some dumps more sensible
!
! Revision 2.29  2005/06/03 01:58:53  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! Revise PFA data structures.
!
! Revision 2.28  2005/05/28 03:27:43  vsnyder
! More work on making PFA databases consistent
!
! Revision 2.27  2005/05/27 23:58:07  vsnyder
! Add Flush PFAData
!
! Revision 2.26  2005/05/24 01:54:55  vsnyder
! Bug fixes -- no channel, wrong subscript...
!
! Revision 2.25  2005/05/13 00:21:07  livesey
! More bug fixes.  Things not being written/read correctly.
!
! Revision 2.24  2005/05/02 23:03:16  vsnyder
! Stuff for PFA Cacheing
!
! Revision 2.23  2005/04/19 20:16:27  livesey
! Bug fix for case when no initial vGrids on reading pfa file.
!
! Revision 2.22  2005/04/04 19:53:05  vsnyder
! Make Read_PFADatum_H5 subroutine
!
! Revision 2.21  2005/03/28 20:25:31  vsnyder
! Add WhichLines, SpectroscopyFile metadata
!
! Revision 2.20  2005/03/17 01:32:26  vsnyder
! Put spectroscopy file's string index in PFAData structure
!
! Revision 2.19  2005/03/17 00:00:09  vsnyder
! Spiff up a dump
!
! Revision 2.18  2005/03/03 21:12:36  vsnyder
! Remove UseMolecule from WritePFAData, remove unreferenced symbols
!
! Revision 2.17  2005/02/05 01:39:12  vsnyder
! Handle separate readPFA commands correctly
!
! Revision 2.16  2005/02/05 00:02:12  vsnyder
! Read all of the specified tables
!
! Revision 2.15  2005/01/27 21:19:44  vsnyder
! Remove unformatted, nonscalar molecule
!
! Revision 2.14  2005/01/12 23:59:44  vsnyder
! Use CAPITALIZE to test requested PFA tables
!
! Revision 2.13  2005/01/12 03:17:41  vsnyder
! Read and write PFA data in HDF5
!
! Revision 2.12  2004/12/31 02:41:24  vsnyder
! Working on read/write PFA database
!
! Revision 2.11  2004/12/13 23:59:37  pwagner
! Re-ordered use hdf5 statements to avoid familiar Lahey internal compiler error
!
! Revision 2.10  2004/12/13 20:41:40  vsnyder
! Filled in Write_PFADatabase.  Handle HDF5 in Write_PFADatum.  Some cannonball
! polishing.
!
! Revision 2.9  2004/11/04 03:42:09  vsnyder
! Provide for both LBL_Ratio and PFA_Ratio in beta_group
!
! Revision 2.8  2004/10/06 21:19:50  vsnyder
! Add sorting and comparing, some cannonball polishing
!
! Revision 2.7  2004/09/04 01:50:31  vsnyder
! Got checked in with get_beta_path_m.f90 for some reason
!
! Revision 2.6  2004/09/02 00:50:15  vsnyder
! Replace velLin with vel_cor
!
! Revision 2.5  2004/09/01 00:28:54  vsnyder
! Make kind parameters more abstract, improve some comments
!
! Revision 2.4  2004/07/17 02:28:52  vsnyder
! Add 'details' arguments for dumps
!
! Revision 2.3  2004/06/17 00:18:23  vsnyder
! Added Write_PFADatum
!
! Revision 2.2  2004/06/09 17:53:13  vsnyder
! OOPS -- got the module name wrong in the new file
!
! Revision 2.1  2004/06/09 17:46:43  vsnyder
! Initial commit after splitting from PFAData.f90
!
! Revision 2.3  2004/06/08 19:29:27  vsnyder
! Add file field
!
! Revision 2.2  2004/05/29 02:51:40  vsnyder
! Allow signal string to denote only one signal
!
! Revision 2.1  2004/05/22 02:29:48  vsnyder
! Initial commit
!
