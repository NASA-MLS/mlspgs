! Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

module PFADataBase_m

  ! Read PFA data.  Build a database.  Provide for access to it.
  ! Write PFA data.

  use MLSCommon, only: R4
  use MLSSignals_m, only: MaxSigLen, Signal_T
  use Molecules, only: First_Molecule, Last_Molecule
  use VGridsDatabase, only: VGrid_t
  ! use for HDF5 intentionally last to avoid long LF95 compile times
  use HDF5, only: Hid_t

  implicit NONE
  private
  public :: PFAData_t, PFAFile_t, PFAPointer_t
  public :: PFAData, PFAFiles, RK, SortPFAData
  public :: AddPFADatumToDatabase
  public :: Destroy_All_PFAData_Arrays, Destroy_PFADataBase, Destroy_PFADatum
  public :: Destroy_PFADatum_Arrays, Destroy_PFAFile, Destroy_PFAFiles
  public :: Dump, Dump_PFADataBase, Dump_PFADatum
  public :: Dump_PFAFileDataBase, Dump_PFAFileDatum, GetGroupName
  public :: HookTableToFindPFA
  public :: PFA_By_Molecule, PFADataOrder, PFADataOrderIndexed, Process_PFA_File
  public :: Read_PFADataBase, Read_PFADatum_H5, Sort_PFADataBase
  public :: Test_And_Fetch_PFA, Write_PFADatum, Write_PFADataBase

  interface Dump
    module procedure Dump_PFADatum, Dump_PFAFileDatum
  end interface Dump

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
    integer :: GroupIndex = 0          ! Index in PFAFiles%PFAData
    real(rk), pointer :: Absorption(:,:) => NULL() ! Ln Absorption data, T x P
    real(rk), pointer :: dAbsDwc(:,:) => NULL()    ! d Ln Absorption / d wc data
    real(rk), pointer :: dAbsDnc(:,:) => NULL()    ! d Ln Absorption / d nc data
    real(rk), pointer :: dAbsDnu(:,:) => NULL()    ! d Ln Absorption / d nu data
  end type PFAData_T

  type PFAFile_T ! For all the PFA tables in one HDF file
    integer :: FileName                ! Sub-rosa index of file name
    integer(hid_t) :: HDF5_FileID      ! HDF5 file ID if open
    logical :: Open = .false.          ! "HDF5 file is open"
    type(PFAData_t), pointer :: PFAData(:) => NULL() ! The groups
  end type PFAFile_T

  ! All PFA files mentioned in GlobalSettings
  type(PFAFile_t), pointer :: PFAFiles(:) => NULL()

  type(PFAData_t), pointer, save :: PFAData(:) => NULL()

  type :: PFAPointer_T ! to make an array of pointers to PFAData_t
    type(PFAData_t), pointer :: Datum => NULL()  ! The data
  end type PFAPointer_T

  ! Structures for finding PFA data quickly.  There are four levels of tables,
  ! indexed respectively by molecule, signal, sideband and channel.  We do this
  ! instead of having a rank-4 array because some molecules have no PFA data,
  ! and of those that do, some have no data for some signals, and of those that
  ! do, some have no data for one sideband, and then the number of channels is
  ! different for different signals.

  type, public :: BySideband_T
    type(PFAPointer_t), pointer :: C(:) => NULL() ! Indexed by channel
  end type BySideband_T

  type, public :: BySignal_T
    type(bySideband_T) :: SB(2)        ! Indexed by sideband, 1 = LSB, 2 = USB
  end type BySignal_T

  type, public :: ByMolecule_T
    type(bySignal_t), pointer :: S(:) => NULL()   ! Indexed by signal
  end type ByMolecule_T

  type(byMolecule_t), public, save :: FindPFA(first_molecule:last_molecule)

  ! PFAData(SortPFAData(i)) < PFAData(SortPFAData(j)) if i < j, with "<"
  ! defined by PFADataOrder.
  integer, pointer, save :: SortPFAData(:) => NULL()

  ! PFA_By_Molecule is used for quick access to the PFA data.
  ! PFA_By_Molecule(m) indexes the last element of SortPFAData for
  ! molecule m.  PFA_By_Molecule(first_molecule-1) is zero.  So the
  ! number of PFA for molecule m is PFA_By_Molecule(m) - PFA_By_Molecule(m-1).
  integer, save :: PFA_By_Molecule(first_molecule-1:last_molecule)
  data PFA_By_Molecule(first_molecule-1) /0/

  integer, parameter, public :: MolNameLen = 31 ! Length for Molecule names in files

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------  AddPFADatumToDatabase  -----
  integer function AddPFADatumToDatabase ( DATABASE, ITEM )

  ! This routine adds a PFA Datum to a database of PFA Data, creating the
  ! database if necessary.

    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate, &
      & MLSMSG_Error
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
    integer :: F, G ! Loop inductors, subscripts
    if ( .not. associated(PFAFiles) ) return
    do f = 1, size(PFAFiles)
      if ( associated(PFAFiles(f)%PFAData) ) then
        do g = 1, size(PFAFiles(f)%PFAData)
          call Destroy_PFADatum_Arrays ( PFAFiles(f)%PFAData(g) )
        end do
      end if
    end do
  end subroutine

  ! ----------------------------------------  Destroy_PFADataBase  -----
  subroutine Destroy_PFADataBase
    use Allocate_Deallocate, only: Deallocate_Test
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error
    integer :: I
    if ( .not. associated(pfaData) ) return
    do i = 1, size(pfaData)
      call destroy_PFADatum ( pfaData(i) )
    end do
    deallocate ( PFAData, stat=i )
    if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Deallocate // 'PDAData' )
    call deallocate_test ( SortPFAData, 'SortPFAData', moduleName )
  end subroutine Destroy_PFADataBase

  ! -------------------------------------------  Destroy_PFADatum  -----
  subroutine Destroy_PFADatum ( PFADatum )
    use Allocate_Deallocate, only: Deallocate_Test
    type(PFAData_T), intent(inout) :: PFADatum
    call deallocate_test ( PFADatum%theSignal%channels, 'PFADatum...Channels', &
      & moduleName )
    call Destroy_PFADatum_Arrays ( PFADatum )
  end subroutine Destroy_PFADatum

  ! ------------------------------------  Destroy_PFADatum_Arrays  -----
  subroutine Destroy_PFADatum_Arrays ( PFADatum )
    use Allocate_Deallocate, only: Deallocate_Test
    type(PFAData_T), intent(inout) :: PFADatum
    call deallocate_test ( PFADatum%absorption, 'PFADatum%absorption', moduleName )
    call deallocate_test ( PFADatum%dAbsDwc, 'PFADatum%dAbsDwc', moduleName )
    call deallocate_test ( PFADatum%dAbsDnc, 'PFADatum%dAbsDnc', moduleName )
    call deallocate_test ( PFADatum%dAbsDnu, 'PFADatum%dAbsDnu', moduleName )
  end subroutine Destroy_PFADatum_Arrays

  ! --------------------------------------------  Destroy_PFAFile  -----
  subroutine Destroy_PFAFile ( PFAFile )
  ! Destroy one PFAFile_t object
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error
    type(PFAFile_t), intent(inout) :: PFAFile
    integer :: I
    if ( associated(PFAFile%PFAData) ) then
      do i = 1, size(PFAFile%PFAData)
        call destroy_PFADatum ( PFAFile%PFAData(i) )
      end do
      deallocate ( PFAFile%PFAData, stat=i )
      if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Deallocate // 'PFAFile%PFAData' )
    end if
  end subroutine Destroy_PFAFile

  ! -------------------------------------------  Destroy_PFAFiles  -----
  subroutine Destroy_PFAFiles
  ! Destroy the elements of the PFAFiles array, then deallocate it.
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error
    integer :: I
    if ( associated(PFAFiles) ) then
      do i = 1, size(PFAFiles)
        call destroy_PFAFile ( PFAFiles(i) )
      end do
      deallocate ( PFAFiles, stat=i )
      if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & MLSMSG_Deallocate // 'PFAFiles' )
    end if
  end subroutine Destroy_PFAFiles

  ! -------------------------------------------  Dump_PFADataBase  -----
  subroutine Dump_PFADataBase ( Details )
    use Output_m, only: Output
    integer, intent(in), optional :: Details
    integer :: I, J
    if ( .not. associated(pfaData) ) then
      call output ( 'No PFA Database to dump', advance='yes' )
      return
    end if
    do i = 1, size(pfaData)
      j = i
      if ( associated(sortPFAdata) ) then
        if ( size(sortPFAdata) == size(pfaData) ) j = sortPFAdata(i)
      end if
      call dump_PFADatum ( pfaData(j), details, j )
    end do
  end subroutine Dump_PFADataBase

  ! ----------------------------------------------  Dump_PFADatum  -----
  subroutine Dump_PFADatum ( PFADatum, Details, Index )

    use Dump_0, only: Dump
    use Intrinsic, only: Lit_Indices
    use MLSSignals_m, only: DisplaySignalName
    use Physics, only: SpeedOfLight
    use String_Table, only: Display_String
    use Output_m, only: NewLine, Output
    use VGridsDatabase, only: Dump

    type(PFAData_t), intent(in) :: PFADatum
    integer, intent(in), optional :: Details ! <= 0 -> Signal and molecule only
                                             ! 1 -> Dump complete summary
                                             ! 2 -> Dump Betas (default)
                                             ! 3 -> Dump Betas and unnamed grids
                                             ! >3 -> Dump Betas, grids and derivatives
    integer, intent(in), optional :: Index   ! In PFA Database

    integer, parameter :: CK = kind(speedOfLight)
    real(ck) :: C = speedOfLight / 1000.0_ck ! km/s
    integer :: MyDetails
    character(len=*), parameter :: WhichLines(0:2) = &
      (/ 'Channel   ', &
       & 'Radiometer', &
       & 'Catalog   ' /)

    myDetails = 2
    if ( present(details) ) myDetails = details

    call output ( 'PFA Datum' )
    if ( present(index) ) call output ( index, before=' ' )
    if ( pfaDatum%name /= 0 ) then
      if ( present(index) ) call output ( ': ' )
      call display_string ( pfaDatum%name )
    end if
    if ( myDetails == 0 ) then
      call display_string ( lit_indices(pfaDatum%molecule), before=' for ' )
      call display_string ( pfaDatum%signal, before=' and ', advance='yes' )
      return
    end if
    call newLine
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
      call displaySignalName ( pfaDatum%theSignal, advance='yes' )
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
    integer, intent(in), optional :: Details ! 0 (default) => file names,
                                             ! 1 => group names too
                                             ! 2 => location too
                                             ! >2 => dump group with details-3

    integer :: I

    if ( .not. associated(PFAFiles) )then
      call output ( 'No PFA File Database to dump', advance='yes' )
      return
    end if

    do i = 1, size(PFAFiles)
      call dump ( PFAFiles(i), details )
    end do
  end subroutine Dump_PFAFileDataBase

  ! ------------------------------------------  Dump_PFAFileDatum  -----
  subroutine Dump_PFAFileDatum ( PFAFileDatum, Details )
    use Intrinsic, only: Lit_Indices
    use MLSSignals_m, only: MaxSigLen
    use Output_m, only: Blanks, NewLine, Output
    use String_Table, only: Display_String
    type(PFAFile_t), intent(in) :: PFAFileDatum
    integer, intent(in), optional :: Details ! 0 (default) => file names,
                                             ! >0 => group names too
                                             ! >1 => signal stuff too
                                             ! >2 => dump group with details-3

    integer :: G, MyDetails
    character(len=molNameLen+1+maxSigLen) :: GroupName

    myDetails = 0
    if ( present(details) ) myDetails = details

    call display_string ( PFAFileDatum%fileName, before='PFA File ', strip=.true. )
    if ( myDetails > 0 ) then
      call output ( ' has groups:', advance='yes' )
      if ( associated(PFAFileDatum%PFAData) ) then
        do g = 1, size(PFAFileDatum%PFAData)
          if ( myDetails > 2 ) then
            call dump ( PFAFileDatum%PFAData(g), myDetails-3 )
          else
            call blanks ( 1 )
            call getGroupName ( PFAFileDatum%PFAData(g), groupName )
            call output ( trim(groupName) )
            if ( myDetails == 2 ) then
              call output ( PFAFileDatum%PFAData(g)%signalIndex, before=' Signal Index ' )
              call output ( PFAFileDatum%PFAData(g)%channel, before=' Channel ' )
              call output ( PFAFileDatum%PFAData(g)%theSignal%sideband, before=' Sideband ' )
            end if
            call newLine
          end if
        end do
      else
        call output ( ' has no associated groups', advance='yes' )
      end if
    else
      call newLine
    end if

  end subroutine Dump_PFAFileDatum

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
  logical function HookTableToFindPFA ( F, G, PFADatum, Replace )
  ! Hook a PFA datum to the FindPFA structure, so that it can be found
  ! quickly given its molecule index, signal index, sideband and channel.
  ! Return "Found one in the table already" (if replace is present and true),
  ! so the caller can avoid memory leaks.
  ! See Test_And_Fetch_PFA.
    use Intrinsic, only: Lit_Indices
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error
    use MLSSignals_m, only: MaxSigLen, Signals
    use MoreMessage, only: MLSMessage
    use String_Table, only: Get_String
    integer, intent(in) :: F  ! Index in PFAFiles; zero if not from a file
                              ! specified in a PFAFile= parameter
    integer, intent(in) :: G  ! Index in PFAFiles(f)%PFAData; zero if not from a file
    type(PFAData_t), intent(in), target :: PFADatum
    logical, intent(in), optional :: Replace ! Replace old one, default false

    integer :: M, S, SB, C ! Molecule, Signal, Sideband, Channel
    logical :: MyReplace
    integer :: STAT

    myReplace = .false.
    if ( present(replace) ) myReplace = replace

    m = PFADatum%molecule
    s = PFADatum%signalIndex
    sb = (PFADatum%theSignal%sideband + 3 ) / 2
    c = PFADatum%channel

    if ( .not. associated(findPFA(m)%s) ) then
      allocate ( findPFA(m)%s(1:size(signals)), stat=stat )
      if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate // 'FindPFA(m)%s).' )
    end if
    if ( .not. associated(findPFA(m)%s(s)%sb(sb)%c) ) then
      allocate ( findPFA(m)%s(s)%sb(sb)%c( &
        & lbound(signals(s)%frequencies,1):ubound(signals(s)%frequencies,1)), &
        & stat=stat )
      if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate // 'FindPFA(m)%s(s)%sb(sb)%c).' )
    end if
    hookTableToFindPFA = associated(findPFA(m)%s(s)%sb(sb)%c(c)%datum)
    if ( hookTableToFindPFA ) then
      if ( myReplace ) then
        findPFA(m)%s(s)%sb(sb)%c(c)%datum = PFADatum
      else
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Duplicate PFA data for %s and %s', &
        & (/ lit_indices(PFADatum%molecule), PFADatum%signal /) )
      end if
    else
      findPFA(m)%s(s)%sb(sb)%c(c)%datum => PFADatum
    end if
    findPFA(m)%s(s)%sb(sb)%c(c)%datum%fileIndex = f
    findPFA(m)%s(s)%sb(sb)%c(c)%datum%groupIndex = g

  end function HookTableToFindPFA

  ! ------------------------------------------------- PFADataOrder -----
  pure integer function PFADataOrder ( A, B ) result ( N )
    type (PFAData_t), intent(in) :: A, B
    n = a%molecule - b%molecule
    if ( n /= 0 ) return
    n = a%theSignal%band - b%theSignal%band
    if ( n /= 0 ) return
    n = a%theSignal%instrumentModule - b%theSignal%instrumentModule
    if ( n /= 0 ) return
    n = a%theSignal%radiometer - b%theSignal%radiometer
    if ( n /= 0 ) return
    n = a%theSignal%spectrometer - b%theSignal%spectrometer
    if ( n /= 0 ) return
    n = lbound(a%theSignal%channels,1) - lbound(b%theSignal%channels,1)
    if ( n /= 0 ) return
    n = a%theSignal%sideband - b%theSignal%sideband
  end function PFADataOrder

  ! ------------------------------------------ PFADataOrderIndexed -----
  pure integer function PFADataOrderIndexed ( A, B )
    integer, intent(in) :: A, B
    PFADataOrderIndexed = PFADataOrder ( PFAData(a), PFAData(b) )
  end function PFADataOrderIndexed

  ! -------------------------------------------  Process_PFA_File  -----
  integer function Process_PFA_File ( PFAFileIndex, Where )
  ! Process a PFA file name from the PFAFile parameter in GlobalSettings
  ! Open the HDF5 file, read the Signals and Molecules groups, allocate
  ! the groups, but don't open the groups.  Leave the file open.

    use Allocate_Deallocate, only: Deallocate_Test
    use MLSHDF5, only: LoadPtrFromHDF5DS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
      & MLSMSG_FileOpen
    use MLSSignals_m, only: Signals
    use MoreMessage, only: MLSMessage
    use MoreTree, only: GetLitIndexFromString, GetStringIndexFromString
    use Parse_Signal_m, only: Parse_Signal
    use String_Table, only: Get_String
    ! HDF5 intentionally last to avoid long LF95 compiles
    use HDF5, only: H5F_ACC_RDONLY_F, H5FOpen_F, H5GClose_f, H5GOpen_F

    integer, intent(in) :: PFAFileIndex ! Sub-rosa index for file name
    integer, intent(in) :: Where  ! Source_ref field, for error messages

    integer :: F, G                    ! Subscripts, loop inductors
    logical, pointer :: Channels(:)
    integer(hid_t) :: GroupID          ! From HDF5 open group
    character(len=molNameLen), pointer :: MyMolecules(:)
    character(len=maxSigLen), pointer :: MySignalStrings(:) ! From the HDF5
    type(PFAFile_t) :: PFAFileDatum
    character(len=1023) :: PFAFileName
    integer :: SB                      ! Sideband, from the signal
    integer, pointer :: SignalIndices(:)
    integer :: STAT                    ! From HDF5 open or read, or allocate

    ! Open the file
    PFAFileDatum%fileName = PFAFileIndex
    call get_string ( PFAFileIndex, PFAFileName, strip=.true. )
    call h5fopen_f ( trim(PFAFileName), H5F_ACC_RDONLY_F, &
      & PFAFileDatum%HDF5_fileID, stat )
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Fileopen // ' %s (HDF5 PFA data) at %l.', (/PFAFileIndex,where/) )
    PFAFileDatum%open = .true.

    ! Open the Index group and read the Molecules and Signals data sets
    call h5gOpen_f ( PFAFileDatum%HDF5_fileID, 'Index', groupID, stat )
    if ( stat /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open HDF5 PFA Index group in ' // &
        & trim(PFAFileName) )
    nullify ( myMolecules, mySignalStrings )
    call loadPtrFromHDF5DS ( groupID, 'Molecules', myMolecules )
    call loadPtrFromHDF5DS ( groupID, 'Signals', mySignalStrings )
    call h5gClose_f ( groupID, stat )
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close HDF5 PFA index group.' )

    ! Allocate the groups component
    allocate ( PFAFileDatum%PFAData(size(myMolecules)), stat=stat )
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      MLSMSG_Allocate // 'PFAFileDatum%PFAData' )

    ! Fill each group's fields, except for the data and group ID.
    ! PFAFileDatum%PFAData(g)%open is default initialized to be .false.
    nullify ( channels, signalIndices )
    do g = 1, size(myMolecules)
      PFAFileDatum%PFAData(g)%molecule = &
        & getLitIndexFromString ( trim(myMolecules(g)) )
      PFAFileDatum%PFAData(g)%signal = &
        & getStringIndexFromString ( trim(mySignalStrings(g)) )
      call parse_signal ( trim(mySignalStrings(g)), signalIndices, &
        & channels=channels, sideband=sb )
      if ( .not. associated(signalIndices) ) call MLSMessage ( MLSMSG_Error, &
        & moduleName, 'Unable to parse signal ' // trim(mySignalStrings(g)) )
      PFAFileDatum%PFAData(g)%channel = lbound(channels,1)
      PFAFileDatum%PFAData(g)%signalIndex = signalIndices(1)
      PFAFileDatum%PFAData(g)%theSignal = signals(signalIndices(1))
      PFAFileDatum%PFAData(g)%theSignal%channels => channels
      PFAFileDatum%PFAData(g)%theSignal%sideband = sb
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
    do g = 1, size(PFAFileDatum%PFAData)
      if ( hookTableToFindPFA ( f, g, PFAFileDatum%PFAData(g) ) ) continue
    end do

    process_PFA_File = f

  contains

    ! ................................  AddPFAFileDatumToDatabase  .....
    integer function AddPFAFileDatumToDatabase ( DATABASE, ITEM )

    ! This routine adds a PFA File Datum to a database of PFA File Data,
    ! creating the database if necessary.

      use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate, &
        & MLSMSG_Error
      ! Dummy arguments
      type (PFAFile_T), dimension(:), pointer :: DATABASE
      type (PFAFile_T), intent(in) :: ITEM

      ! Local variables
      type (PFAFile_T), dimension(:), pointer :: tempDatabase

      include "addItemToDatabase.f9h"

      AddPFAFileDatumToDatabase = newSize
    end function AddPFAFileDatumToDatabase

  end function Process_PFA_File

  ! -------------------------------------------  Read_PFADatabase  -----
  subroutine Read_PFADatabase ( FileNameIndex, FileTypeIndex, TheMolecules, &
    & TheSignals, VGrids, Where )
  ! Read the PFA data from FileName.  If both TheMolecules and
  ! TheSignals have zero size, all PFA data are read from FileName.
  ! If TheMolecules has zero size but TheSignals does not, the PFA
  ! data for all molecules for each specified signal are read from
  ! FileName.  If TheSignals has zero size but TheMolecules does
  ! not, the PFA data for all signals for each specified molecule are read
  ! from FileName. Otherwise the PFA data for the Cartesian product of
  ! TheMolecules and TheSignals are read from FileName.

    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    use MLSHDF5, only: LoadPtrFromHDF5DS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
      & MLSMSG_Error, MLSMSG_FileOpen
    use MLSStrings, only: Capitalize
    use MoreMessage, only: MLSMessage
    use MoreTree, only: GetLitIndexFromString, GetStringIndexFromString
    use Parse_Signal_m, only: Get_Individual_Signals
    use String_Table, only: Get_String
    ! HDF5 intentionally last to avoid long LF95 compiles
    use HDF5, only: H5F_ACC_RDONLY_F, H5FOpen_F, H5FClose_F, &
      & H5GClose_F, H5GOpen_f

    integer, intent(in) :: FileNameIndex
    integer, intent(in) :: FileTypeIndex ! HDF5 is all we can do
    integer, intent(in) :: TheMolecules(:), TheSignals(:)
    type(vGrid_t), pointer :: VGrids(:)
    integer, intent(in) :: Where ! Source_ref field, for error messages

    integer, pointer :: AllSignals(:) ! String table indices for AllSignalStrings
    character(len=maxSigLen), pointer :: AllSignalStrings(:) ! from expanding a signal
    integer :: F ! Index in PFAFiles
    integer(hid_t) :: GroupID
    character(len=255) :: FileName ! , FileType
    integer, pointer :: Groups(:) ! Indices in MyMolecules, MySignals
    integer :: I, IOSTAT, IPFA, J, K, L
    integer, pointer :: MyMolecules(:)
    character(len=molNameLen), pointer :: MyMoleculeStrings(:)
    integer, pointer :: MySignals(:)
    character(len=maxSigLen), pointer :: MySignalStrings(:) ! From the HDF5
    integer :: NGroups, NPFA
    character(len=maxSigLen) :: OneSignal
    type(PFAData_t), pointer, save :: TempPFAData(:) => NULL()
    character(len=molNameLen+maxSigLen+1) :: TheGroup

    call get_string ( fileNameIndex, fileName, strip=.true. )
    f = 0
    if ( associated(PFAFiles) ) then
      do f = size(PFAFiles), 1, -1
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
      allSignals(i) = getLitIndexFromString ( allSignalStrings(i) )
    end do
!   if ( fileType == 'HDF5' ) then
      nullify ( groups, myMolecules, myMoleculeStrings, mySignals, mySignalStrings )
      if ( .not. PFAFiles(f)%open ) then
        call h5fopen_f ( trim(fileName), H5F_ACC_RDONLY_F, &
          & PFAFiles(f)%HDF5_FileID, iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Fileopen // ' %s (HDF5 PFA data) at %l.', (/fileNameIndex,where/) )
      end if
      ! Construct the list of group names
      call h5gOpen_f ( PFAFiles(f)%HDF5_FileID, 'Index', groupID, iostat )
      if ( iostat /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to open HDF5 PFA index group in ' // &
          & trim(fileName) )
      call loadPtrFromHDF5DS ( groupID, 'Molecules', myMoleculeStrings )
      call allocate_test ( myMolecules, size(myMoleculeStrings), &
        'MyMolecules', moduleName )
      do i = 1, size(myMoleculeStrings)
        myMolecules(i) = getLitIndexFromString ( myMoleculeStrings(i) )
      end do
      call loadPtrFromHDF5DS ( groupID, 'Signals', mySignalStrings )
      call allocate_test ( mySignals, size(mySignalStrings), &
        'MySignals', moduleName )
      do i = 1, size(mySignalStrings)
        mySignals(i) = getStringIndexFromString ( mySignalStrings(i) )
      end do
      call h5gClose_f ( groupID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close HDF5 PFA index group.' )
      nGroups = size(mySignals)
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
            if ( allSignals(i) == mySignals(k) ) then
              nPFA = nPFA + 1
            end if
          end do
        end do
        call allocate_test ( groups, nPFA, 'Groups', moduleName )
        nPFA = 0
        do i = 1, size(allSignals)
          do k = 1, nGroups
            if ( allSignals(i) == mySignals(k) ) then
              nPFA = nPFA + 1
              groups(nPFA) = k
            end if
          end do
        end do
      else if ( size(allSignals) == 0 ) then
        do i = 1, size(theMolecules)
          do k = 1, nGroups
            if ( theMolecules(i) == myMolecules(k) ) then
              nPFA = nPFA + 1
            end if
          end do
        end do
        call allocate_test ( groups, nPFA, 'Groups', moduleName )
        nPFA = 0
        do i = 1, size(theMolecules)
          do k = 1, nGroups
            if ( theMolecules(i) == myMolecules(k) ) then
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
              if ( myMolecules(k) == theMolecules(i) .and. &
                & mySignals(k) == theSignals(j) ) then
                nPFA = nPFA + 1
              end if
            end do
          end do
        end do
        call allocate_test ( groups, nPFA, 'Groups', moduleName )
        nPFA = 0
        do i = 1, size(theMolecules)
          do j = 1, size(allSignals)
            do k = 1, nGroups
              if ( myMolecules(k) == theMolecules(i) .and. &
                & mySignals(k) == theSignals(j) ) then
                nPFA = nPFA + 1
                groups(nPFA) = k
              end if
            end do
          end do
        end do
      end if
      call deallocate_test ( allSignals, 'AllSignals', moduleName )
      call deallocate_test ( allSignalStrings, 'AllSignalStrings', moduleName )
      call deallocate_test ( myMolecules, 'MyMolecules', moduleName )
      call deallocate_test ( mySignals, 'MySignals', moduleName )
      if ( nPFA == 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'No PFA Data to read from ' // trim(fileName) )
      ! Create or expand the PFADatabase
      if ( associated(PFAData) .and. nPFA > 0 ) then
        tempPFAData => PFAData
        iPFA = size(tempPFAData)
        allocate ( PFAData(iPFA+nPFA), stat=iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & MLSMSG_Allocate // 'PFAData' )
        pfaData(:iPFA) = tempPFAData
        deallocate ( tempPFAData, stat=iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & MLSMSG_DeAllocate // 'TempPFAData' )
      else
        iPFA = 0
        if ( nPFA > 0 ) then
          allocate ( PFAData(nPFA), stat=iostat )
          if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
            & MLSMSG_Allocate // 'PFAData' )
        end if
      end if
      ! Read the groups
      do i = 1, nPFA
        theGroup = trim(myMoleculeStrings(groups(i))) // '%' // &
                   trim(mySignalStrings(groups(i)))
        call h5gOpen_f ( PFAFiles(f)%HDF5_FileID, theGroup, groupID, iostat )
        if ( iostat /= 0 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to open HDF5 PFA group ' // theGroup // ' in ' // &
            & trim(fileName) )
        iPFA = iPFA + 1
        call read_PFADatum_H5 ( groupID, PFAData(iPFA), .true., VGrids )
        call h5gClose_f ( groupID, iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close HDF5 PFA group ' // theGroup // ' in ' // &
            & trim(fileName) )
      end do ! i = 1, nPFA
      call H5FClose_F ( PFAFiles(f)%HDF5_FileID, iostat )
      PFAFiles(f)%open = .false.
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close PFA data file ' // trim(fileName) // '.' )
      call deallocate_test ( groups, 'Groups', moduleName )
      call deallocate_test ( myMoleculeStrings, 'MyMoleculeStrings', moduleName )
      call deallocate_test ( mySignalStrings, 'MySignalStrings', moduleName )
!   else
!     Nothing -- Read_PFAData checks file type
!   end if

    call sort_PFADatabase

  end subroutine Read_PFADatabase

  ! -------------------------------------------  Read_PFADatum_H5  -----
  subroutine Read_PFADatum_H5 ( GroupID, PFADatum, Derivs, VGrids, Data, F, G )
  ! Read the PFA Datum from HDF5 group GroupId into PFADatum.
  ! Read dAbsDwc and dAbsDnc if and only if Derivs is true.  Add any
  ! grids that aren't already in VGrids.

    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    use Intrinsic, only: L_Theta, L_Zeta
    use MLSHDF5, only: GetHDF5Attribute, IsHDF5AttributePresent, LoadPtrFromHDF5DS
    use MLSMessageModule, only: MLSMessage,  MLSMSG_Error
    use MLSSignals_m, only: MaxSigLen, Signals
    use MoreTree, only: GetLitIndexFromString, GetStringIndexFromString
    use Parse_Signal_m, only: Parse_Signal
    use VGridsDatabase, only: AddVGridIfNecessary, RS ! Kind for Surfs

    integer(hid_t), intent(in) :: GroupId
    type(PFAData_T), target :: PFADatum
    logical, intent(in) :: Derivs
    type(vGrid_t), pointer :: VGrids(:)
    logical, intent(in), optional :: Data ! "Read the data -- default true"
    integer, intent(in), optional :: F, G ! Where to hook into PFAFiles,
                                          ! f == zero (default) means don't do it

    logical, pointer :: Channels(:) ! output from Parse_Signal
    integer :: J, K
    character(1023) :: Line ! Text, e.g. filter file name
    character(len=molNameLen) :: Molecule
    logical :: MyData
    integer :: MyF, MyG
    integer, pointer :: SignalIndices(:) ! output from Parse_Signal
    character(len=maxSigLen) :: SignalText
    real(rs) :: SurfStep ! for temperature and pressure grids
    type(vGrid_t) :: TGrid, VGrid
    integer :: GRIDINDEX

    myData = .true.
    if ( present(data) ) myData = data

    myF = 0; myG = 0
    if ( present(f) ) myF = f
    if ( present(g) ) myG = g

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
    call getHDF5Attribute ( groupID, 'molecule', molecule )
    k = getLitIndexFromString ( trim(molecule) )
    if ( k < first_molecule .or. k > last_molecule ) call MLSMessage ( &
      & MLSMSG_Error, moduleName, 'The string ' // &
      & trim(molecule) // ' is not a molecule name.' )
    PFADatum%molecule = k
    call getHDF5Attribute ( groupID, 'signal', signalText )
    PFADatum%signal = getStringIndexFromString ( trim(signalText), .true. )
    call parse_signal ( signalText, signalIndices, &
      & channels=channels )
    PFADatum%signalIndex = signalIndices(1)
    PFADatum%theSignal = signals(PFADatum%signalIndex)
    PFADatum%theSignal%channels => channels
    nullify ( channels ) ! so as not to clobber PFADatum%theSignal%channels
      ! in next iteration of the loop
    call getHDF5Attribute ( groupID, 'sideband', PFADatum%theSignal%sideband )
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
    gridIndex = addVGridIfNecessary(tGrid,vGrids)
    PFADatum%tGrid = vGrids(gridIndex)
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
    PFADatum%vGrid = vGrids(addVGridIfNecessary(vGrid,vGrids, &
      &                       relErr=vGrid%noSurfs*0.2_rs*epsilon(1.0_rs)))

    if ( myData ) then
      call loadPtrFromHDF5DS ( groupID, 'absorption', PFADatum%absorption )
      call loadPtrFromHDF5DS ( groupID, 'dAbsDnu', PFADatum%dAbsDnu )
      if ( derivs ) then
        call loadPtrFromHDF5DS ( groupID, 'dAbsDwc', PFADatum%dAbsDwc )
        call loadPtrFromHDF5DS ( groupID, 'dAbsDnc', PFADatum%dAbsDnc )
      end if
    end if

    call deallocate_test ( signalIndices, 'SignalIndices', moduleName )

    if ( hookTableToFindPFA ( myF, myG, PFADatum, replace=.true. ) ) continue

  end subroutine Read_PFADatum_H5

  ! -------------------------------------------  Sort_PFADatabase  -----
  subroutine Sort_PFADatabase
  ! Create the array SortPFAData such that PFAData(SortPFAData(i)) <
  ! PFAData(SortPFAData(j)) if i < j. "<" is defined by PFADataOrder: First
  ! the molecules, then the signal's fields, with sideband last.  We use the
  ! molecule for a quick search.  Otherwise, the order isn't really
  ! important.  What we're trying to do is make sure that we find stuff in
  ! the same order when we're processing the upper and lower sidebands.

    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    use Sort_m, only: GSORTP

    integer :: I, J, N

    if ( .not. associated(PFAData) ) return ! Can't sort it if it's not there

    n = size(PFAData)

    if ( associated(sortPFAData) ) then
      if ( size(sortPFAData) == n ) return ! only do it once
      call deallocate_test ( sortPFAData, 'SortPFAData', moduleName )
    end if

    call allocate_test ( sortPFAData, n, 'SortPFAData', moduleName )

    ! Now sort the database
    call gsortp ( PFADataOrderIndexed, n, sortPFAData )

    ! Now that it's sorted with molecule as the major key, make
    ! PFA_By_Molecule.
    j = size(SortPFAData)
    do i = last_molecule, first_molecule, -1
      do while ( j > 0 )
        if ( PFAData(sortPFAData(j))%molecule > i ) then
          j = j - 1
        else
      exit
        end if
      end do
      PFA_by_molecule(i) = j
    end do

  end subroutine Sort_PFADatabase

  ! -----------------------------------------  Test_And_Fetch_PFA  -----
  subroutine Test_And_Fetch_PFA ( Molecule, Signal, Sideband, Channel, Derivs, &
    &                             PFADatum, VGrids )

    ! Test whether the file and group of PFA data are known.  If not, emit an
    ! error message and crash, or return NULL, depending upon the Crash
    ! parameter here.  See HookTableToFindPFA.
    ! If it is known, test whether the datum is filled.  If not, read the
    ! necessary data sets and add any grids that aren't already in VGrids.
    ! Return a pointer associated with the datum.

    use Intrinsic, only: Lit_Indices
    use MLSHDF5, only: LoadPtrFromHDF5DS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only: GetSignalName
    use String_Table, only: Get_String
    ! HDF5 intentionally last to avoid long LF95 compiles
    use HDF5, only: H5F_ACC_RDONLY_F, H5FOpen_F, H5GOpen_F

    integer, intent(in) :: Molecule ! first_molecule : last_molecule
    integer, intent(in) :: Signal   ! 1 : size(signals)
    integer, intent(in) :: Sideband ! -1 = LSB, +1 = USB
    integer, intent(in) :: Channel  ! depends on signal
    logical, intent(in) :: Derivs   ! "Read dAbs dwc and dAbs dnc"
    type(PFAData_t), pointer :: PFADatum
    type(vGrid_t), pointer :: VGrids(:)

    logical, parameter :: Crash = .true.

    integer :: I, SB                ! Subscripts
    character(255) :: Msg           ! could also be a group name
    type(PFAPointer_t), pointer :: PFAPointer
    integer :: STAT

    sb = ( sideband + 3 ) / 2 ! -1..+1 => 1..2
    if ( .not. associated(findPFA(molecule)%s) ) go to 9
    if ( .not. associated(findPFA(molecule)%s(signal)%sb(sb)%c ) ) go to 9
    PFAPointer => findPFA(molecule)%s(signal)%sb(sb)%c(channel)
    if ( .not. associated(PFAPointer) ) go to 9
    PFADatum => PFAPointer%datum
    if ( .not. associated(PFADatum) ) go to 9

    call readDS ( 'absorption', PFADatum%absorption )
    call readDS ( 'dAbsDnu', PFADatum%dAbsDnu )
    if ( .not. derivs ) return
    call readDS ( 'dAbsDwc', PFADatum%dAbsDwc )
    call readDS ( 'dAbsDnc', PFADatum%dAbsDnc )

    return

  9 if ( crash ) then
      msg = 'No PFA data for'
      i = len('No PFA data for')
      call getSignalName ( signal, msg(i+2:), sideband=sideband, channel=channel )
      i = len_trim(msg)
      msg(i+2:i+4) = 'and'
      call get_string ( lit_indices(molecule), msg(i+6:) )
      call MLSMessage ( MLSMSG_Error, moduleName, trim(msg) )
    end if
    return

  contains
    subroutine ReadDS ( DSName, DSData )
      character(len=*), intent(in) :: DSName
      real(rk), pointer :: DSData(:,:)
      integer :: F, G ! Subscripts
      if ( associated(dsData) ) return
      f = PFADatum%fileIndex
      if ( f == 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'PFA Data was not read from a file and therefore cannot be recovered.' )
      g = PFADatum%groupIndex
      if ( .not. associated(PFADatum%theSignal%channels) ) then ! nothing there
        if ( .not. PFAFiles(f)%open ) then
          call get_string ( PFAFiles(f)%fileName, msg )
          call h5fopen_f ( trim(msg), H5F_ACC_RDONLY_F, &
            & PFAFiles(f)%HDF5_fileID, stat )
          if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Unable to open HDF5 PFA file ' // trim(msg) )
          PFAFiles(f)%open = .true.
        end if
        if ( .not. PFAFiles(f)%PFAData(g)%open ) then
          call get_string ( lit_indices(molecule), msg )
          i = len_trim(msg)
          msg(i+1:i+1) = '%'
          call getSignalName ( signal, msg(i+2:), sideband=sideband, channel=channel )
          call h5gOpen_f ( PFAFiles(f)%HDF5_fileID, trim(msg), &
            & PFAFiles(f)%PFAData(g)%HDF5_GroupID, stat )
          if ( stat /= 0 ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Unable to open HDF5 PFA group ' // trim(msg) )
          PFAFiles(f)%PFAData(g)%open = .true.
        end if
        call read_PFADatum_H5 ( PFAFiles(f)%PFAData(g)%HDF5_GroupID, PFADatum, &
          & derivs, vGrids )
      end if
      call loadPtrFromHDF5DS ( PFAFiles(f)%PFAData(g)%HDF5_GroupID, dsName, dsData )
    end subroutine ReadDS
  end subroutine Test_And_Fetch_PFA

  ! ------------------------------------------  Write_PFADatabase  -----
  subroutine Write_PFADatabase ( FileName, FileType )
    use Intrinsic, only: Lit_Indices
    use MLSHDF5, only: SaveAsHDF5DS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only: MaxSigLen
!   use MLSStrings, only: Capitalize
    use String_Table, only: Get_String
    ! HDF5 intentionally last to avoid long LF95 compiles
    use HDF5, only: H5FCreate_F, H5FClose_F, &
      & H5F_ACC_TRUNC_F, H5GClose_F, H5GCreate_F

    character(len=*), intent(in) :: FileName, FileType

    integer :: FileID, GroupID
    integer :: I, IOSTAT
    character(len=molNameLen) :: Molecules(size(pfaData))
    character(len=maxSigLen) :: SignalText(size(pfaData))

    if ( .not. associated(pfaData) ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'No PFA Data to write' )
!   if ( capitalize(fileType) == 'HDF5' ) then ! open HDF5 file here
      call H5FCreate_F ( trim(fileName), H5F_ACC_TRUNC_F, fileID, &
        & iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to open hdf5 PFA file ' // trim(fileName) // ' for output.' )
      ! Make the Index group
      call h5gCreate_f ( fileID, 'Index', groupID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to create hdf5 Index group in ' // trim(fileName) // '.' )
      do i = 1, size(pfaData)
        call get_string ( lit_indices(pfaData(i)%molecule), molecules(i) )
        call get_string ( pfaData(i)%signal, signalText(i) )
      end do
      call saveAsHDF5DS ( groupID, 'Molecules', molecules )
      call saveAsHDF5DS ( groupID, 'Signals', signalText )
      call h5gClose_F ( groupID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to close hdf5 Index group in ' // trim(fileName) // '.' )
      ! Write the PFA Data tables, one per group
      do i = 1, size(pfaData)
        call write_PFADatum ( pfaData(i), FileName, FileType, lun=fileID )
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

    use Intrinsic, only: Lit_Indices
    use MLSHDF5, only: MakeHDF5Attribute, SaveAsHDF5DS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only: MaxSigLen
!   use MLSStrings, only: Capitalize
    use Output_m, only: Output
    use String_Table, only: Get_String, String_Length
    use Toggles, only: Switches
    ! HDF5 intentionally last to avoid long LF95 compiles
    use HDF5, only: H5FCreate_F, H5FClose_F, & ! HDF5 USE intentionally last
      & H5F_ACC_TRUNC_F, H5GClose_F, H5GCreate_F

    type(PFAData_t), intent(in) :: PFADatum
    character(len=*), intent(in) :: FileName, FileType
    integer, intent(in), optional :: Lun ! Don't open a new file if present

    character(len=1023) :: Attrib
    integer :: GroupID
    character(len=molNameLen+1+maxSigLen+1) :: GroupName
    integer :: IOSTAT, MyLun
    character(len=molNameLen) :: Molecule
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
      if ( index(switches,'pfaw') /= 0 ) then
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
      call MakeHDF5Attribute ( groupID, 'molecule', molecule )
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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PFADataBase_m

! $Log$
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
