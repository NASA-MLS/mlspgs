! Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

module PFADataBase_m

  ! Read PFA data.  Build a database.  Provide for access to it.
  ! Write PFA data.

  use MLSCommon, only: R4
  use MLSSignals_m, only: MaxSigLen, Signal_T
  use Molecules, only: First_Molecule, Last_Molecule
  use VGridsDatabase, only: VGrid_t

  implicit NONE
  private
  public :: PFAData_t, PFAData, RK, SortPFAData
  public :: AddPFADatumToDatabase
  public :: Destroy_PFADataBase, Dump, Dump_PFADataBase, Dump_PFADatum
  public :: PFA_By_Molecule, PFADataOrder, PFADataOrderIndexed, Sort_PFADataBase
  public :: Read_PFADataBase, Write_PFADatum, Write_PFADataBase

  interface Dump
    module procedure Dump_PFADatum
  end interface Dump

  integer, parameter :: RK = r4 ! Kind of real fields in PFAData_t

  type PFAData_T
    integer :: Name = 0                            ! of the pfaData spec
    integer :: FilterFile = 0                      ! Index in string table
    integer :: Molecule                            ! Molecule index
    character(len=maxSigLen) :: Signal             ! The signal string
    integer :: SignalIndex                         ! in Signals database
    integer :: SpectroscopyFile = 0                ! Index in string table
    type(signal_t) :: TheSignal                    ! The signal, with channels
                                                   ! and sidebands added
    type(vGrid_t) :: TGrid ! shallow copy from VGrids database of Log temperatures
    type(vGrid_t) :: VGrid ! shallow copy from VGrids database of vertical grid
    real(rk) :: Vel_Rel                            ! vel_lin / c
    real(rk), pointer :: Absorption(:,:) => NULL() ! Ln Absorption data, T x P
    real(rk), pointer :: dAbsDwc(:,:) => NULL()    ! d Ln Absorption / d wc data
    real(rk), pointer :: dAbsDnc(:,:) => NULL()    ! d Ln Absorption / d nc data
    real(rk), pointer :: dAbsDnu(:,:) => NULL()    ! d Ln Absorption / d nu data
  end type PFAData_T

  type(PFAData_t), pointer, save :: PFAData(:) => NULL()

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

  ! ----------------------------------------  Destroy_PFADataBase  -----
  subroutine Destroy_PFADataBase
    use Allocate_Deallocate, only: Deallocate_Test
    use MLSMessageModule, only: MLSMessage, MLSMSG_Deallocate, MLSMSG_Error
    integer :: I
    if ( .not. associated(pfaData) ) return
    do i = 1, size(pfaData)
      call deallocate_test ( pfaData(i)%theSignal%channels, 'pfaData...Channels', &
          & moduleName )
      call deallocate_test ( pfaData(i)%absorption, 'pfaData%absorption', moduleName )
      call deallocate_test ( pfaData(i)%dAbsDwc, 'pfaData%dAbsDwc', moduleName )
      call deallocate_test ( pfaData(i)%dAbsDnc, 'pfaData%dAbsDnc', moduleName )
      call deallocate_test ( pfaData(i)%dAbsDnu, 'pfaData%dAbsDnu', moduleName )
    end do
    deallocate ( PFAData, stat=i )
    if ( i /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Deallocate // 'PDAData' )
    call deallocate_test ( SortPFAData, 'SortPFAData', moduleName )
  end subroutine Destroy_PFADataBase

  ! -------------------------------------------  Dump_PFADataBase  -----
  subroutine Dump_PFADataBase ( Details )
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, intent(in), optional :: Details
    integer :: I, J
    if ( .not. associated(pfaData) ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Cannot dump unallocated PFA data base" )
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
    integer, intent(in), optional :: Details ! <= 0 -> Don't dump arrays,
                                             ! 1 -> Dump Betas
                                             ! 2 -> Dump Betas and unnamed grids
                                             ! >2 -> Dump Betas, grids and derivatives
    integer, intent(in), optional :: Index   ! In PFA Database

    integer, parameter :: CK = kind(speedOfLight)
    real(ck) :: C = speedOfLight / 1000.0_ck ! km/s
    integer :: MyDetails

    myDetails = 1
    if ( present(details) ) myDetails = details

    call output ( 'PFA Datum ' )
    if ( present(index) ) call output ( index )
    if ( pfaDatum%name /= 0 ) then
      if ( present(index) ) call output ( ': ' )
      call display_string ( pfaDatum%name )
    end if
    call newLine
    call output ( ' Molecule: ' )
    call display_string ( lit_indices(pfaDatum%molecule), advance='yes' )

    call output ( ' Specified signal: ' )
    call output ( trim(pfaDatum%signal), advance='yes' )
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
      & after='kms', advance='yes' )

    if ( pfaDatum%tGrid%name /= 0 ) then
      call display_string ( pfaDatum%tGrid%name, before=' TGrid: ' )
    else if ( myDetails > 1 ) then
      call output ( ' TGrid: ' )
      call dump ( pfaDatum%tGrid )
    end if

    if ( pfaDatum%vGrid%name /= 0 ) then
      call display_string ( pfaDatum%vGrid%name, before=' VGrid: ' )
    else if ( myDetails > 1 ) then
      if ( pfaDatum%tGrid%name /= 0 ) call newLine
      call dump ( pfaDatum%vGrid )
    end if

    if ( pfaDatum%tGrid%name /= 0 .and. pfaDatum%vGrid%name == 0 .and. &
      &  myDetails > 1 &
      & .or. pfaDatum%vGrid%name /= 0 ) &
      & call newLine

    if ( myDetails <= 0 ) return

    call dump ( pfaDatum%absorption, name=' ln Absorption' )
    if ( myDetails <= 2 ) return
    call dump ( pfaDatum%dAbsDwc, name=' d ln Absorption / d wc' )
    call dump ( pfaDatum%dAbsDnc, name=' d ln Absorption / d nc' )
    call dump ( pfaDatum%dAbsDnu, name=' d ln Absorption / d nu' )

  end subroutine Dump_PFADatum

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

  ! -------------------------------------------  Read_PFADatabase  -----
  subroutine Read_PFADatabase ( FileName, FileType, TheMolecules, &
    & TheSignalStrings, VGrids )
  ! Read the PFA data from FileName.  If both TheMolecules and
  ! TheSignalStrings have zero size, all PFA data are read from FileName. 
  ! If TheMolecules has zero size but TheSignalStrings does not, the PFA
  ! data for all molecules for each specified signal are read from
  ! FileName.  If TheSignalStrings has zero size but TheMolecules does
  ! not, the PFA data for all signals for each specified molecule are read
  ! from FileName. Otherwise the PFA data for the Cartesian product of
  ! TheMolecules and TheSignalStrings are read from FileName.

    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    use Intrinsic, only: L_Theta, L_Zeta
    use MLSHDF5, only: GetHDF5Attribute, IsHDF5AttributePresent, LoadPtrFromHDF5DS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
      & MLSMSG_Error
    use MLSSignals_m, only: Signals
    use MLSStrings, only: Capitalize
    use MoreTree, only: GetLitIndexFromString, GetStringIndexFromString
    use Parse_Signal_m, only: Get_Individual_Signals, Parse_Signal
    use VGridsDatabase, only: AddVGridIfNecessary, RS ! Kind for Surfs
    ! HDF5 intentionally last to avoid long LF95 compiles
    use HDF5, only: H5F_ACC_RDONLY_F, H5FOpen_F, H5FClose_F, &
      & H5GClose_F, H5GOpen_f, Hid_t

    character(len=*), intent(in) :: FileName
    character(len=*), intent(in) :: FileType ! Upper case, HDF5 or UNFORMATTED
    character(len=*), intent(in) :: TheMolecules(:), TheSignalStrings(:)
    type(vGrid_t), pointer :: VGrids(:)

    character(len=maxSigLen), pointer :: AllSignals(:) ! from expanding theSignalStrings
    logical, pointer :: Channels(:) ! output from Parse_Signal
    integer(hid_t) :: FileID, GroupID
    character(len=molNameLen+maxSigLen+1) :: GroupTrial
    character(len=len(groupTrial)), pointer :: Groups(:), GroupsTest(:)
    integer :: I, IOSTAT, IPFA, J, K, L
    character(1023) :: Line ! Text, e.g. filter file name
    character(len=molNameLen) :: Molecule
    character(len=molNameLen), pointer :: MyMolecules(:)
    character(len=maxSigLen), pointer :: MySignalStrings(:) ! From the HDF5
    integer :: NPFA
    integer, pointer :: SignalIndices(:) ! output from Parse_Signal
    real(rs) :: SurfStep ! for temperature and pressure grids
    type(PFAData_t), pointer, save :: TempPFAData(:) => NULL()
    type(vGrid_t) :: TGrid, VGrid

    nullify ( allSignals, channels, signalIndices )
    ! Expand theSignalStrings to allSignals
    call get_individual_signals ( theSignalStrings, allSignals )
!   if ( fileType == 'HDF5' ) then
      nullify ( groups, groupsTest, myMolecules, mySignalStrings )
      call h5fopen_f ( trim(fileName), H5F_ACC_RDONLY_F, fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open PFA data file ' // trim(fileName) // '.' )
      ! Construct the list of group names
      call h5gOpen_f ( fileID, 'Index', groupID, iostat )
      if ( iostat /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to open HDF5 PFA index group in ' // &
          & trim(fileName) )
      call loadPtrFromHDF5DS ( groupID, 'Molecules', myMolecules )
      call loadPtrFromHDF5DS ( groupID, 'Signals', mySignalStrings )
      call h5gClose_f ( groupID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close HDF5 PFA index group.' )
      nPFA = size(mySignalStrings)
      call allocate_test ( groupsTest, nPFA, 'Groups', moduleName )
      do i = 1, nPFA
        groupsTest(i) = trim(myMolecules(i)) // '%' // trim(mySignalStrings(i))
      end do
      ! Decide what to read
      l = max(size(theMolecules) * size(allSignals), &
        &     size(theMolecules), size(allSignals) )
      nPFA = 0
      if ( l == 0 ) then ! Read everything
        groups => groupsTest
        nullify ( groupsTest )
        nPFA = size(groups)
      else if ( size(theMolecules) == 0 ) then
        do i = 1, size(allSignals)
          do k = 1, size(groupsTest)
            if ( capitalize(allSignals(i)) == capitalize(mySignalStrings(k)) ) then
              nPFA = nPFA + 1
            end if
          end do
        end do
        call allocate_test ( groups, nPFA, 'Groups', moduleName )
        nPFA = 0
        do i = 1, size(allSignals)
          do k = 1, size(groupsTest)
            if ( capitalize(allSignals(i)) == capitalize(mySignalStrings(k)) ) then
              nPFA = nPFA + 1
              groups(nPFA) = groupsTest(k)
            end if
          end do
        end do
      else if ( size(allSignals) == 0 ) then
        do i = 1, size(theMolecules)
          do k = 1, size(groupsTest)
            if ( theMolecules(i) == myMolecules(k) ) then
              nPFA = nPFA + 1
            end if
          end do
        end do
        call allocate_test ( groups, nPFA, 'Groups', moduleName )
        nPFA = 0
        do i = 1, size(theMolecules)
          do k = 1, size(groupsTest)
            if ( theMolecules(i) == myMolecules(k) ) then
              nPFA = nPFA + 1
              groups(nPFA) = groupsTest(k)
            end if
          end do
        end do
      else ! Both allSignals and theMolecules are present -- get as much
           ! of the Cartesian product of them as possible
        do i = 1, size(theMolecules)
          do j = 1, size(allSignals)
            groupTrial = trim(theMolecules(i)) // '%' // trim(allSignals(j))
            do k = 1, size(groupsTest)
              if ( capitalize(groupsTest(k)) == capitalize(groupTrial) ) then
                nPFA = nPFA + 1
              end if
            end do
          end do
        end do
        call allocate_test ( groups, nPFA, 'Groups', moduleName )
        nPFA = 0
        do i = 1, size(theMolecules)
          do j = 1, size(allSignals)
            groupTrial = trim(theMolecules(i)) // '%' // trim(allSignals(j))
            do k = 1, size(groupsTest)
              if ( capitalize(groupsTest(k)) == capitalize(groupTrial) ) then
                nPFA = nPFA + 1
                groups(nPFA) = groupTrial
              end if
            end do
          end do
        end do
      end if
      call deallocate_test ( groupsTest, 'GroupsTest', moduleName )
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
        call h5gOpen_f ( fileID, trim(groups(i)), groupID, iostat )
        if ( iostat /= 0 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to open HDF5 PFA group ' // trim(groups(i)) // ' in ' // &
            & trim(fileName) )
        iPFA = iPFA + 1
        if ( isHDF5AttributePresent(groupID, 'filterFile') ) then
          call getHDF5Attribute ( groupID, 'filterFile', line )
          PFAData(iPFA)%filterFile = getStringIndexFromString ( trim(line), .true. )
        else
          PFAData(iPFA)%filterFile = 0
        end if
        call getHDF5Attribute ( groupID, 'molecule', molecule )
        k = getLitIndexFromString ( trim(molecule) )
        if ( k < first_molecule .or. k > last_molecule ) call MLSMessage ( &
          & MLSMSG_Error, moduleName, 'The string ' // &
          & trim(molecule) // ' is not a molecule name.' )
        PFAData(iPFA)%molecule = k
        call getHDF5Attribute ( groupID, 'signal', PFAData(iPFA)%signal )
        call parse_signal ( PFAData(iPFA)%signal, signalIndices, &
          & channels=channels )
        PFAData(iPFA)%signalIndex = signalIndices(1)
        PFAData(iPFA)%theSignal = signals(PFAData(iPFA)%signalIndex)
        PFAData(iPFA)%theSignal%channels => channels
        nullify ( channels ) ! so as not to clobber PFAData(iPFA)%theSignal%channels
          ! in next iteration of the loop
        call getHDF5Attribute ( groupID, 'sideband', PFAData(iPFA)%theSignal%sideband )
        call getHDF5Attribute ( groupID, 'vel_rel', PFAData(iPFA)%vel_rel )
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
        PFAData(iPFA)%tGrid = vGrids(addVGridIfNecessary(tGrid,vGrids))
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
        PFAData(iPFA)%vGrid = vGrids(addVGridIfNecessary(vGrid,vGrids, &
          &                       relErr=vGrid%noSurfs*0.2_rs*epsilon(1.0_rs)))
        call loadPtrFromHDF5DS ( groupID, 'absorption', PFAData(iPFA)%absorption )
        call loadPtrFromHDF5DS ( groupID, 'dAbsDwc', PFAData(iPFA)%dAbsDwc )
        call loadPtrFromHDF5DS ( groupID, 'dAbsDnc', PFAData(iPFA)%dAbsDnc )
        call loadPtrFromHDF5DS ( groupID, 'dAbsDnu', PFAData(iPFA)%dAbsDnu )
        call h5gClose_f ( groupID, iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to close HDF5 PFA group ' // trim(groups(i)) // ' in ' // &
            & trim(fileName) )
      end do ! i = 1, nPFA
      call H5FClose_F ( fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close PFA data file ' // trim(fileName) // '.' )
      call deallocate_test ( groups, 'Groups', moduleName )
      call deallocate_test ( groupsTest, 'GroupsTest', moduleName )
      call deallocate_test ( myMolecules, 'MyMolecules', moduleName )
      call deallocate_test ( mySignalStrings, 'MySignalStrings', moduleName )
      call deallocate_test ( signalIndices, 'SignalIndices', moduleName )
!   else
!     Nothing -- Read_PFAData checks file type
!   end if

    call sort_PFADatabase

  end subroutine Read_PFADatabase

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

  ! ------------------------------------------  Write_PFADatabase  -----
  subroutine Write_PFADatabase ( FileName, FileType )
    use Intrinsic, only: Lit_Indices
    use MLSHDF5, only: SaveAsHDF5DS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
!   use MLSStrings, only: Capitalize
    use String_Table, only: Get_String
    ! HDF5 intentionally last to avoid long LF95 compiles
    use HDF5, only: H5FCreate_F, H5FClose_F, &
      & H5F_ACC_TRUNC_F, H5GClose_F, H5GCreate_F

    character(len=*), intent(in) :: FileName, FileType

    integer :: FileID, GroupID
    integer :: I, IOSTAT
    character(len=molNameLen) :: Molecules(size(pfaData))

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
      call saveAsHDF5DS ( groupID, 'Signals', pfaData%signal )
      do i = 1, size(pfaData)
        call get_string ( lit_indices(pfaData(i)%molecule), molecules(i) )
      end do
      call saveAsHDF5DS ( groupID, 'Molecules', molecules )
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
    character(len=molNameLen+1+len(PFADatum%signal)+1) :: GroupName
    integer :: IOSTAT, MyLun
    character(len=molNameLen) :: Molecule

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
      call get_string ( lit_indices(pfaDatum%molecule), molecule )
      groupName = trim(molecule) // '%' // trim(pfaDatum%signal)
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
      call MakeHDF5Attribute ( groupID, 'molecule', molecule )
      call MakeHDF5Attribute ( groupID, 'signal', pfaDatum%signal )
      call MakeHDF5Attribute ( groupID, 'sideband', pfaDatum%theSignal%sideband )
      call MakeHDF5Attribute ( groupID, 'vel_rel', pfaDatum%vel_rel )
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
