! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PFAData_m

  ! Read the PFA data file(s).  Build a database.  Provide for access to it.
  ! Setup to make the PFA Data tables, as specified by a MakePFA statement.

  implicit NONE
  private
  public :: Get_PFAdata_from_l2cf, Make_PFAData, Read_PFAData, Write_PFAData

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------  Get_PDAdata_from_l2cf  -----
  subroutine Get_PFAdata_from_l2cf ( Root, Name, VGrids, Error )
  ! Process a PFAdata specification from the l2cf.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Expr_m, only: Expr
    use Init_Tables_Module, only: F_Absorption, F_dAbsDnc, F_dAbsDnu, &
      & F_dAbsDwc, F_File, F_Molecules, F_Signal, F_Temperatures, F_VelLin, &
      & F_VGrid, Field_First, Field_Last, L_Zeta
    use Intrinsic, only: PHYQ_Dimensionless, PHYQ_Velocity
    use IO_Stuff, only: Get_Lun
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only: Signals
    use MLSStrings, only: Capitalize
    use Molecules, only: T_Molecule
    use MoreTree, only: Get_Field_ID, GetLitIndexFromString
    use Parse_Signal_m, only: Parse_Signal
    use PFADataBase_m, only: AddPFADatumToDatabase, PFAData, PFAData_T, RK, &
      & Write_PFADatum
    use Physics, only: SpeedOfLight
    use String_Table, only: Create_String, Get_String
    use Tree, only: Decorate, Decoration, Node_Id, NSons, Sub_Rosa, Subtree
    use Tree_Checker, only: Check_Type
    use Tree_Types, only: N_String
    use VGridsDatabase, only: VGrid_t

    integer, intent(in) :: Root            ! of the pfaData subtree in the l2cf
    integer, intent(in) :: Name            ! of the pfaData spec, else zero
    type(vGrid_t), intent(in), target :: VGrids(:) ! database of vgrids
    integer, intent(out) :: Error          ! 0 => OK, else trouble

    ! Error codes
    integer, parameter :: CannotOpen = 1
    integer, parameter :: CannotRead = cannotOpen + 1
    integer, parameter :: DSBSignal = cannotRead + 1
    integer, parameter :: NotMolecule = DSBSignal + 1
    integer, parameter :: NotZeta = notMolecule + 1
    integer, parameter :: ShowSize = notZeta + 1
    integer, parameter :: SignalParse = showSize + 1
    integer, parameter :: TooManyChannels = signalParse + 1
    integer, parameter :: TooManySignals = tooManyChannels + 1
    integer, parameter :: UnsupportedFormat = tooManySignals + 1
    integer, parameter :: WrongFields = unsupportedFormat + 1
    integer, parameter :: WrongSignal = wrongFields + 1
    integer, parameter :: WrongSize = wrongSignal + 1
    integer, parameter :: WrongUnits = wrongSize + 1

    integer :: AbsTree
    integer, parameter :: CK = kind(speedOfLight)
    real(ck) :: C = speedOfLight / 1000.0_ck ! km/s
    logical, pointer :: Channels(:)
    integer :: dAbsDncTree, dAbsDnuTree, dAbsDwcTree
    integer :: Field, FileIndex ! Where in the tree is the filename?
    character(255) :: FileName, FileType ! Formatted(default), Unformatted
    character(255) :: FilterFile ! From PFAData file, not from L2cf
    logical :: Got(field_first:field_last)
    integer :: I, IOStat, J
    character(127) :: Line      ! of formatted PFA data file
    integer :: Lun
    integer :: NArrays, NMolT, NPress, NPressT, NTemps, NTempsT
    integer :: Sideband
    integer, pointer :: SignalIndices(:)
    character(127) :: SignalT
    integer :: Son, Units(2)
    type(pfaData_t) :: PFADatum
    double precision :: Value(2)
    real(rk) :: VelLin

    error = 0
    got = .false.
    nullify ( channels, signalIndices )
    pfaDatum%name = name
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = get_field_id(son)
      got(field) = .true.
      select case ( field )
      case ( f_absorption )
        absTree = son
      case ( f_dAbsDnc )
        dAbsDncTree = son
      case ( f_dAbsDnu )
        dAbsDnuTree = son
      case ( f_dAbsDwc )
        dAbsDwcTree = son
      case ( f_file ) ! "filename" [: "<format>"]
        fileIndex = subtree(2,son)
        fileType = 'UNFORMATTED'
        if ( node_id(fileIndex) /= n_string ) then ! must be n_*colon
          call get_string( sub_rosa(subtree(2,fileIndex)), fileType, strip=.true. )
          fileIndex = subtree(1,fileIndex)
        end if
        call get_string ( sub_rosa(fileIndex), fileName, strip=.true. )
      case ( f_molecules )
        nullify ( pfaDatum%molecules )
        call allocate_test ( pfaDatum%molecules, nsons(son)-1, &
          & 'pfaDatum%molecules', moduleName )
        do j = 2, nsons(son)
          pfaDatum%molecules(j-1) = decoration(subtree(j,son))
        end do
      case ( f_signal )
        call get_string ( sub_rosa(subtree(2,son)), pfaDatum%signal, strip=.true. )
        call parse_signal ( pfaDatum%signal, signalIndices, &
          & tree_index=son, sideband=sideband, channels=channels )
        if ( .not. associated(signalIndices) ) & ! A parse error occurred
          & call announce_error ( subtree(2,son), signalParse, pfaDatum%signal )
        if ( size(signalIndices) > 1 ) &
          & call announce_error ( subtree(2,son), tooManySignals, pfaDatum%signal )
        if ( size(channels) > 1 ) &
          & call announce_error ( subtree(2,son), tooManyChannels, pfaDatum%signal )
        pfaDatum%signalIndex = signalIndices(1)
        pfaDatum%theSignal = signals(pfaDatum%signalIndex)
        pfaDatum%theSignal%channels => channels
        pfaDatum%theSignal%sideband = sideband
        call deallocate_test ( signalIndices, 'SignalIndices', moduleName )
      case ( f_temperatures )
        pfaDatum%tGrid = vGrids(decoration(decoration(subtree(2,son))))
      case ( f_velLin )
        call expr ( subtree(2,son), units, value )
        if ( units(1) /= phyq_velocity ) &
          & call announce_error ( subtree(1,son), wrongUnits, 'Velocity' )
        velLin = value(1) / 1000.0 ! fundamental unit is m/s, fwdmdl wants km/s
      case ( f_vGrid )
        pfaDatum%vGrid = vGrids(decoration(decoration(subtree(2,son))))
        if ( pfaDatum%vGrid%verticalCoordinate /= l_zeta ) &
          & call announce_error ( subtree(1,son), notZeta )
      end select
    end do

    nPress = pfaDatum%vGrid%noSurfs
    nTemps = pfaDatum%tGrid%noSurfs
    nArrays = nPress * nTemps + 1

    ! Nullify shouldn't be needed, but just to be safe....
    nullify ( pfaDatum%absorption, pfaDatum%dAbsDnc, pfaDatum%dAbsDnu, pfaDatum%dAbsDwc )
    call allocate_test ( pfaDatum%absorption, nTemps, nPress, 'pfaDatum%absorption', moduleName )
    call allocate_test ( pfaDatum%dAbsDnc,    nTemps, nPress, 'pfaDatum%dAbsDnc',    moduleName )
    call allocate_test ( pfaDatum%dAbsDnu,    nTemps, nPress, 'pfaDatum%dAbsDnu',    moduleName )
    call allocate_test ( pfaDatum%dAbsDwc,    nTemps, nPress, 'pfaDatum%dAbsDwc',    moduleName )

    ! Got(file) and not any of the others means read PFAdata from file
    if ( got(f_file) .and. .not. any( (/  &
       & got(f_absorption), got(f_dAbsDnc), got(f_dAbsDnu), got(f_dAbsDwc), &
       & got(f_molecules), got(f_velLin) /) ) ) then
      call get_lun ( lun )
      j = scan(fileName,'$%')
      if ( j == 0 ) then
        fileName = trim(fileName) // pfaDatum%signal
      else
        fileName = fileName(:j-1) // trim(pfaDatum%signal) // fileName(j+1:)
      end if
      open ( unit=lun, file=fileName, form=fileType, status='old', iostat=iostat )
      if ( iostat /= 0 ) then
        call announce_error ( fileIndex, cannotOpen, fileName, iostat )
      else
        if ( capitalize(fileType) == 'UNFORMATTED' ) then ! Unformatted
          signalT = ''
          read ( lun, iostat=iostat ) nTempsT, nPressT, nMolT, velLin, &
            & i, signalT(:i)
          if ( iostat /= 0 ) &
            & call announce_error ( fileIndex, cannotRead, fileName, iostat )
          ! Check the signal
          if ( signalT /= pfaDatum%signal ) call announce_error ( &
            fileIndex, wrongSignal, fileName )
          ! Check numbers of temperatures and pressures
          if ( nTemps /= nTempsT ) then
            call announce_error ( fileIndex, showSize, &
              & 'Temperature in L2CF', nTemps )
            call announce_error ( fileIndex, showSize, &
              & 'Temperature in File', nTempsT )
          end if
          if ( nPress /= nPressT ) then
            call announce_error ( fileIndex, showSize, &
              & 'Pressure in L2CF', nPress )
            call announce_error ( fileIndex, showSize, &
              & 'Pressure in File', nPressT )
          end if
          if ( error == 0 ) then
            ! Read the absorption and derivative arrays
            read ( lun, iostat=iostat ) pfaDatum%absorption, &
              & pfaDatum%dAbsDwc, pfaDatum%dAbsDnc, pfaDatum%dAbsDnu
            if ( iostat /= 0 ) &
              & call announce_error ( fileIndex, cannotRead, fileName, iostat )
            ! Read and check the molecules
            call allocate_test ( pfaDatum%molecules, nMolT, 'pfaDatum%molecules', moduleName )
            do i = 1, nMolT
              read ( lun, iostat=iostat ) j, line(:j)
              if ( iostat /= 0 ) &
                & call announce_error ( fileIndex, cannotRead, fileName, iostat )
              pfaDatum%molecules(i) = getLitIndexFromString ( line(:j) )
              if ( .not. check_type ( t_molecule, pfaDatum%molecules(i) ) ) &
                & call announce_error ( fileIndex, notMolecule, line(:j) )
              pfaDatum%filterFile = 0
              read ( lun, iostat=iostat ) j, filterFile(:j)
              if ( iostat > 0 ) &
                & call announce_error ( fileIndex, cannotRead, fileName, iostat )
              if ( iostat == 0 ) pfaDatum%filterFile = &
                & create_string ( filterFile(:j), caseless=.false. )
            end do
          end if
        else if ( capitalize(fileType) == 'HDF' ) then ! HDF
          call announce_error ( fileIndex, unsupportedFormat, fileType )
        else
          call announce_error ( fileIndex, unsupportedFormat, fileType )
        end if
      end if
    else ! Got(file) and not the others means write PFAData to file.
      if ( .not. all( (/ & ! Check that we have all required fields
           & got(f_absorption), got(f_dAbsDnc), got(f_dAbsDnu), got(f_dAbsDwc), &
           & got(f_molecules), got(f_velLin) /) ) ) &
        & call announce_error ( root, wrongFields, stop=.true. )
      if ( nSons(absTree) /= nArrays ) &
        call announce_error ( subtree(1,absTree), wrongSize, 'Absorption', &
          & nArrays )
      if ( nSons(dAbsDncTree) /= nArrays ) &
        call announce_error ( subtree(1,dAbsDncTree), wrongSize, 'd Abs / d nc', &
          & nArrays )
      if ( nSons(dAbsDnuTree) /= nArrays ) &
        call announce_error ( subtree(1,dAbsDnuTree), wrongSize, 'd Abs / d nu', &
          & nArrays )
      if ( nSons(dAbsDwcTree) /= nArrays ) &
        call announce_error ( subtree(1,dAbsDwcTree), wrongSize, 'd Abs / d wc', &
          & nArrays )
      ! Get data from the tree into the data structure
      call store_2d ( absTree, pfaDatum%absorption )
      call store_2d ( dAbsDncTree, pfaDatum%dAbsDnc )
      call store_2d ( dAbsDnuTree, pfaDatum%dAbsDnu )
      call store_2d ( dAbsDwcTree, pfaDatum%dAbsDwc )
      ! Write it?
      if ( got(f_file) ) call write_PFADatum ( pfaDatum, fileName, fileType )
    end if
    PFADatum%vel_rel = velLin / c ! Doppler correction factor
    if ( PFADatum%theSignal%sideband == 0 ) &
      & call announce_error ( fileIndex, DSBSignal )

    if ( error == 0 ) then
      call decorate ( root, addPFADatumToDatabase ( pfaData, pfaDatum ) )
    else
      call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Execution terminated.' )
    end if

  contains

    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( Where, What, String, More, Stop )
      use Machine, only: IO_Error
      use MoreTree, only: StartErrorMessage
      use Output_m, only: Output
      use PFADataBase_m, only: Dump
      integer, intent(in) :: Where             ! Tree node index
      integer, intent(in) :: What              ! Error index
      character(len=*), intent(in), optional :: String  ! For more info
      integer, intent(in), optional :: More    ! For more info
      logical,  intent(in),optional :: Stop    ! Stop via MLSMessage if true
      error = 1
      call startErrorMessage ( where )
      select case ( what )
      case ( cannotOpen )
        call output ( 'Cannot open file ' )
        call output ( trim(string) )
        call output ( more, before='.  IOSTAT = ', after='.', advance='yes' )
        call io_error ( 'Cannot open file ', more, trim(string) )
      case ( cannotRead )
        call output ( 'Cannot read file ' )
        call output ( trim(string) )
        call output ( more, before='.  IOSTAT = ', after='.', advance='yes' )
        call io_error ( 'Cannot read file ', more, trim(string) )
      case ( DSBSignal )
        call output ( 'DSB signals not allowed for PFA', advance='yes' )
        call dump ( pfaDatum )
      case ( notMolecule )
        call output ( 'Symbol ' )
        call output ( trim(string) )
        call output ( ' read from file is not a molecule.', advance='yes' )
      case ( notZeta )
        call output ( 'Vertical coordinate for pressure grid must be Zeta.', &
          & advance='yes' )
      case ( showSize )
        call output ( 'Size of ' )
        call output ( trim(string) )
        call output ( more, before=' = ', after='.', advance='yes' )
      case ( signalParse )
        call output ( 'Unable to parse signal ' )
        call output ( trim(string), advance='yes' )
      case ( tooManyChannels )
        call output ( string )
        call output ( ' Describes more than one channel.', advance='yes' )
      case ( tooManySignals )
        call output ( string )
        call output ( ' Describes more than one signal.', advance='yes' )
      case ( unsupportedFormat )
        call output ( trim(string) )
        call output ( ' is not a supported file format.', advance='yes' )
      case ( wrongFields )
        call output ( 'If file appears, either none of absorption, dAbsDnc, dAbsDnu, dAbsDwc,', &
          advance='yes' )
        call output ( 'molecules or velLin shall appear, or all shall appear.', &
          advance='yes' )
      case ( wrongSignal )
        call output ( 'The signal in file ' )
        call output ( trim(string) )
        call output ( ' is not consistent with the SIGNAL in the PFADATA.', advance='yes' )
      case ( wrongSize )
        call output ( 'Incorrect size for ' )
        call output ( trim(string) )
        call output ( more, before=' -- should be ', after='.', advance='yes' )
      case ( wrongUnits )
        call output ( 'Incorrect units -- should be ' )
        call output ( trim(string) )
        call output ( '.', advance='yes' )
      end select
      if ( present(stop) ) then
        if ( stop ) call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Execution terminated.' )
      end if
    end subroutine Announce_Error

    ! .................................................  ReadLine  .....
    subroutine ReadLine
    ! Read from LUN into LINE with '(a)' format.
    !  Announce an error if it fails.
      read ( lun, '(a)', iostat=iostat ) line
      if ( iostat /= 0 ) &
        & call announce_error ( fileIndex, cannotRead, fileName, iostat, .true. )
    end subroutine ReadLine

    ! .................................................  Store_2d  .....
    subroutine Store_2d ( Where, What )
    ! Store data from Where in the L2CF into an nTemps X nPress
    ! array What

      use MLSCommon, only: R4
      integer, intent(in) :: Where
      real(r4), pointer :: What(:,:)
      integer :: I, J, K

      j = 0
      k = 1
      do i = 2, nsons(where)
        call expr ( subtree(i,where), units, value )
        if ( units(1) /= phyq_dimensionless ) &
          & call announce_error ( subtree(i,where), wrongUnits, &
            & 'Dimensionless' )
        j = j + 1
        if ( j > nTemps ) then
          j = 1
          k = k + 1
        end if
        what(j,k) = value(1)
      end do
    end subroutine Store_2d

  end subroutine Get_PFAdata_from_l2cf

  ! -----------------------------------------------  Make_PFAData  -----
  subroutine Make_PFAData ( Root, VGrids, Error )

    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    use Create_PFAData_m, only: Create_PFAData
    use Expr_m, only: Expr
    use FilterShapes_m, only: FilterShapes
    use Init_Tables_Module, only: F_LOSVEL, F_Molecules, F_Signals, &
      & F_Temperatures, F_VGrid, L_Zeta
    use Intrinsic, only: PHYQ_Velocity
    use MLSCommon, only: RP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
      & MLSMSG_Error
    use MLSSignals_m, only: Signal_T, Signals
    use MoreTree, only: Get_Field_ID
    use Parse_Signal_m, only: Parse_Signal
    use String_Table, only: Get_String
    use Tree, only: Decorate, Decoration, Node_Id, NSons, Sub_Rosa, Subtree
    use Tree_types, only: N_Array
    use VGridsDatabase, only: VGrid_t

    integer, intent(in) :: Root ! of the MakePFA subtree
    type(vGrid_T), intent(in), target :: VGrids(:) ! Both temperature and pressure
    integer, intent(out) :: Error

    logical, pointer :: Channels(:)
    integer :: Field, I, J, K, Son
    real(rp) :: LosVel
    integer, pointer :: Molecules(:)
    type(signal_t), pointer :: MySignals(:), SignalsTemp(:)
    integer :: Sideband
    character(127) :: Signal
    integer, pointer :: SignalIndices(:)
    integer :: Stat
    integer :: TGrid, VGrid ! Indices in vGrids database
    integer :: Units(2) ! of losVel
    double precision :: Values(2) ! of losVel

    ! Error codes
    integer, parameter :: NoFilterShapes = 1
    integer, parameter :: NoFolded = noFilterShapes + 1
    integer, parameter :: NoGroup = noFolded + 1
    integer, parameter :: NotVelocity = noGroup + 1
    integer, parameter :: NotZeta = notVelocity + 1
    integer, parameter :: SignalParse = notZeta + 1

    ! Gather the data from the command
    error = 0
    nullify ( signalIndices )
    allocate ( mySignals(0), stat=stat )
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Allocate // 'Signals' )
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = get_field_id(son)
      select case ( field )
      case ( f_losvel )
        call expr ( subtree(2,son), units, values ) ! can't be a range
        if ( units(1) /= phyq_velocity ) call announce_error ( son, notVelocity )
        losVel = values(1)
      case ( f_molecules )
        nullify ( molecules )
        call allocate_test ( molecules, nsons(son)-1, 'molecules', moduleName )
        do j = 2, nsons(son)
          k = subtree(j,son)
          if ( node_id(k) == n_array ) call announce_error ( son, noGroup )
          molecules(j-1) = decoration(k)
        end do
      case ( f_signals )
        nullify ( channels )
        do j = 2, nsons(son)
          call get_string ( sub_rosa(subtree(j,son)), signal, strip=.true. )
          call parse_signal ( signal, signalIndices, tree_index=son, &
            & sideband=sideband, channels=channels )
          if ( .not. associated(signalIndices) ) & ! A parse error occurred
            & call announce_error ( subtree(j,son), signalParse, signal )
          if ( sideband == 0 ) &
            & call announce_error ( subtree(j,son), noFolded, signal )
          signalsTemp => mySignals
          allocate ( mySignals(size(signalsTemp)+size(signalIndices)), stat=stat )
          if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
            & MLSMSG_Allocate // 'Signals' )
          mySignals(:size(signalsTemp)) = signalsTemp
          mySignals(size(signalsTemp)+1:) = signals(signalIndices)
          mySignals(size(signalsTemp)+1:)%sideband = sideband
          ! Unlike the full forward model, use all of the matching signals.
          do k = size(signalsTemp)+1, size(mySignals)
            mySignals(k)%sideband = sideband
            ! Don't hose channels in database
            nullify ( mySignals(k)%channels )
            call allocate_test ( mySignals(k)%channels, &
              & size(mySignals(k)%frequencies), 'mySignals%channels', &
              & ModuleName ) ! , lowBound=lbound(mySignals(k)%frequencies,1) )
            if ( associated(channels) ) then
              mySignals(k)%channels(1:lbound(channels,1)-1) = .false.
              mySignals(k)%channels(lbound(channels,1):ubound(channels,1)) = &
                & channels
              mySignals(k)%channels(ubound(channels,1)+1:) = .false.
            else
              mySignals(k)%channels = .true.
            end if
          end do ! k
          call deallocate_test ( channels, 'Channels', moduleName )
          call deallocate_test ( signalIndices, 'signalIndices', moduleName )
          deallocate ( signalsTemp, stat=stat )
          if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
            & MLSMSG_DeAllocate // 'SignalsTemp' )
        end do ! j
      case ( f_temperatures )
        tGrid = decoration(decoration(subtree(2,son)))
      case ( f_vGrid )
        vGrid = decoration(decoration(subtree(2,son)))
        if ( vGrids(vGrid)%verticalCoordinate /= l_zeta ) &
          & call announce_error ( subtree(1,son), notZeta )
      end select
    end do ! i
    if ( .not. associated(filterShapes) ) call announce_error ( root, noFilterShapes )
    if ( error /= 0 ) return

    call decorate ( root, &
      & create_PFAData ( molecules, mySignals, vGrids(tGrid), vGrids(vGrid), &
      & losVel, root ) )

  contains

    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( Where, What, String )
      use MoreTree, only: StartErrorMessage
      use Output_m, only: Output
      integer, intent(in) :: Where             ! Tree node index
      integer, intent(in) :: What              ! Error index
      character(len=*), intent(in), optional :: String
      error = 1
      call startErrorMessage ( where )
      select case ( what )
      case ( noFilterShapes )
        call output ( 'Filter shapes are required to compute PFA', advance='yes' )
      case ( noFolded )
        call output ( 'Folded sidaband PFA not allowed: ' )
        call output ( trim(string), advance='yes' )
      case ( noGroup )
        call output ( 'Molecule groups not allowed.', advance='yes' )
      case ( notVelocity )
        call output ( 'Units are not velocity.', advance='yes' )
      case ( notZeta )
        call output ( 'Vertical coordinate for pressure grid must be Zeta.', &
          & advance='yes' )
      case ( signalParse )
        call output ( 'Unable to parse signal ' )
        call output ( trim(string), advance='yes' )
      end select
    end subroutine Announce_Error

  end subroutine Make_PFAData

  ! -----------------------------------------------  Read_PFAData  -----
  subroutine Read_PFAData ( Root, VGrids )

    ! Read PFA data.  If the file= field is a range, the second element
    ! specifies the format.  Default is "UNFORMATTED".

    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    use Init_Tables_Module, only: F_File, F_Molecules, F_Signals
    use PFADataBase_m, only: Read_PFADatabase
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only: MaxSigLen
    use MLSStrings, only: Capitalize
    use MoreTree, only: Get_Field_ID
    use String_Table, only: Get_String
    use Tree, only: Node_Id, NSons, Sub_Rosa, Subtree
    use Tree_Types, only: N_String
    use VGridsDatabase, only: VGrid_t

    integer, intent(in) :: Root ! of the MakePFA subtree, below name if any
    type(vGrid_t), pointer :: VGrids(:) ! database of vgrids

    logical :: Error
    character(255) :: FileName, FileType ! HDF5(default), Unformatted
    integer :: I, J
    character(31), pointer :: Molecules(:)
    character(maxSigLen), pointer :: Signals(:)
    integer :: Son

    integer, parameter :: BadFileType = 1 ! Neither HDF5 nor Unformatted
    integer, parameter :: NoRange = badFileType + 1 ! Range not allowed

    error = .false.
    nullify ( molecules, signals )
    do i = 2, nsons(root)
      son = subtree(i,root)
      select case ( get_field_id(son) )
      case ( f_file )
        j = subtree(2,son)
        fileType = 'HDF5'
        if ( node_id(j) /= n_string ) then ! must be n_*colon
          call get_string ( sub_rosa(subtree(2,j)), fileType, strip=.true. )
          fileType = capitalize(fileType)
          if ( fileType /= 'HDF5' .and. fileType /= 'UNFORMATTED' ) &
            & call announce_error ( subtree(2,j), badFileType )
          j = subtree(1,j)
        end if
        call get_string ( sub_rosa(j), fileName, strip=.true. )
      case ( f_molecules )
        call fillArray ( son, molecules, 'Molecules' )
      case ( f_signals )
        call fillArray ( son, signals, 'Signals' )
      end select
    end do
    if ( .not. associated(molecules) ) &
      & call allocate_test ( molecules, 0, 'Molecules', moduleName )
    if ( .not. associated(signals) ) &
      & call allocate_test ( signals, 0, 'Signals', moduleName )
    if ( error ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Error trying to read PFAData' )
    call read_PFADatabase ( fileName, fileType, molecules, signals, vGrids )
    call deallocate_test ( molecules, 'Molecules', moduleName )
    call deallocate_test ( signals, 'Signals', moduleName )

  contains
    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( Where, What )
      use MoreTree, only: StartErrorMessage
      use Output_m, only: Output
      integer, intent(in) :: Where, What
      error = .true.
      call startErrorMessage ( where )
      select case ( what )
      case ( badFileType )
        call output ( 'Unsupported file type.', advance='yes' )
      case ( noRange )
        call output ( 'Range not allowed.', advance='yes' )
      end select
    end subroutine Announce_Error

    ! ................................................  FillArray  .....
    subroutine FillArray ( Where, Array, ArrayName )
      use MLSStrings, only: Capitalize
      use Tree, only: Node_Kind, Pseudo
      integer, intent(in) :: Where ! in the tree
      character(len=*), pointer :: Array(:)
      character(len=*), intent(in) :: ArrayName
      integer :: Gson, J
      call allocate_test ( array, nsons(son)-1, arrayName, moduleName )
      do j = 2, nsons(where)
        gson = subtree(j,where)
        if ( node_kind(gson) == pseudo ) then
          call get_string ( sub_rosa(gson), array(j-1), strip=.true. )
          array(j-1) = capitalize(array(j-1))
        else
          call announce_error ( where, noRange )
        end if
      end do
    end subroutine FillArray

  end subroutine Read_PFAData

  ! ----------------------------------------------  Write_PFAData  -----
  subroutine Write_PFAData ( Root, Error )

    ! Write PFA data on a file generated from the file= field, the molecule,
    ! and the signal.  If the file= field is a range, the second element
    ! specifies the format.  Default is "UNFORMATTED".

    use Init_Tables_Module, only: F_AllPFA, F_File, F_PFAData
    use MoreTree, only: Get_Boolean, Get_Field_ID
    use PFADataBase_m, only: PFAData, Write_PFADatabase, Write_PFADatum
    use String_Table, only: Get_String
    use Tree, only: Decoration, Node_Id, NSons, Sub_Rosa, Subtree
    use Tree_Types, only: N_String

    integer, intent(in) :: Root ! of the MakePFA subtree, below name if any
    integer, intent(out) :: Error ! 0 => OK

    logical :: AllPFA
    integer :: Field, FileIndex ! Where in the tree is the filename?
    character(255) :: FileName, FileType ! Unformatted(default), ???
    integer :: I
    integer :: PFATree          ! Tree index of PFAData=...
    integer :: Son

    integer, parameter :: AtLeastOne = 1           ! Either allPFA and pfaDATA
    integer, parameter :: NotBoth = atLeastOne + 1 ! not both allPFA and pfaDATA

    allPFA = .false.
    error = 0
    pfaTree = 0
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = get_field_id(son)
      select case ( field )
      case ( f_allPFA )
        allPFA = get_boolean(son)
      case ( f_file )
        fileIndex = subtree(2,son)
        fileType = 'UNFORMATTED'
        if ( node_id(fileIndex) /= n_string ) then ! must be n_*colon
          call get_string( sub_rosa(subtree(2,fileIndex)), fileType, strip=.true. )
          fileIndex = subtree(1,fileIndex)
        end if
        call get_string ( sub_rosa(fileIndex), fileName, strip=.true. )
      case ( f_pfaData )
        pfaTree = son
      end select
    end do

    if ( .not. allPFA .and. pfaTree == 0 ) call announce_error ( root, atLeastOne )
    if ( allPFA .and. pfaTree /= 0 )  call announce_error ( root, notBoth )

    if ( error /= 0 ) return

    if ( allPFA ) then
      call write_PFADataBase ( fileName, fileType )
    else
      do i = 2, nsons(pfaTree)
        call write_PFADatum ( pfaData(decoration(decoration(subtree(i,pfaTree)))), &
          & fileName, fileType, useMolecule=.true. )
      end do
    end if

  contains

    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( Where, What )
      use MoreTree, only: StartErrorMessage
      use Output_m, only: Output
      integer, intent(in) :: Where             ! Tree node index
      integer, intent(in) :: What              ! Error index
      error = 1
      call startErrorMessage ( where )
      select case ( what )
      case ( atLeastOne )
        call output ( 'Either AllPFA or PFAData is required', advance='yes' )
      case ( notBoth )
        call output ( 'Cannot specify both AllPFA and PFAData', advance='yes' )
      end select
    end subroutine Announce_Error

  end subroutine Write_PFAData

! =====     Private Procedures     =====================================

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PFAData_m

! $Log$
! Revision 2.14  2005/01/13 00:00:13  vsnyder
! Delete an unreferenced use name
!
! Revision 2.13  2005/01/12 03:18:22  vsnyder
! Read and write PFAData in HDF5
!
! Revision 2.12  2004/12/31 02:41:56  vsnyder
! Working on read/write PFA database
!
! Revision 2.11  2004/12/13 23:58:48  vsnyder
! Add Make_PFAData; add HDF5 to Write_PFAData
!
! Revision 2.10  2004/09/25 00:29:44  vsnyder
! Don't know how the defective one got committed....
!
! Revision 2.9  2004/09/24 23:44:05  vsnyder
! Don't allow DSB signal
!
! Revision 2.8  2004/09/05 21:14:10  pwagner
! component of PFAData type renamed vel_rel
!
! Revision 2.7  2004/09/02 00:49:38  vsnyder
! Replace velLin with vel_cor
!
! Revision 2.6  2004/07/08 19:33:23  vsnyder
! Set up to read unformatted files
!
! Revision 2.5  2004/06/09 19:58:55  pwagner
! Corrected module name to PFADataBase_m
!
! Revision 2.4  2004/06/09 17:47:10  vsnyder
! Split off PFADataBase to fwdmdl
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
