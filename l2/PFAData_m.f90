! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PFAData_m

  ! Read the PFA data file(s).  Build a database.  Provide for access to it.

  implicit NONE
  private
  public :: Get_PFAdata_from_l2cf

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
    use String_Table, only: Get_String
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
    integer, parameter :: NotMolecule = cannotRead + 1
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
    logical :: Write

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
      case ( f_file )
        fileIndex = subtree(2,son)
        fileType = 'unformatted'
        if ( node_id(fileIndex) /= n_string ) then
          call get_string( sub_rosa(subtree(2,fileIndex)), fileType, strip=.true. )
          fileIndex = subtree(1,fileIndex)
        end if
        call get_string ( sub_rosa(fileIndex), fileName, strip=.true. )
      case ( f_molecules )
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
        pfaDatum%tGrid => vgrids(decoration(decoration(subtree(2,son))))
      case ( f_velLin )
        call expr ( subtree(2,son), units, value )
        if ( units(1) /= phyq_velocity ) &
          & call announce_error ( subtree(1,son), wrongUnits, 'Velocity' )
        velLin = value(1) / 1000.0 ! fundamental unit is m/s, fwdmdl wants km/s
      case ( f_vGrid )
        pfaDatum%vGrid => vgrids(decoration(decoration(subtree(2,son))))
        if ( pfaDatum%vGrid%verticalCoordinate /= l_zeta ) &
          & call announce_error ( subtree(1,son), notZeta )
      end select
    end do

    write = got(f_file) .and. all( (/ &
           & got(f_absorption), got(f_dAbsDnc), got(f_dAbsDnu), got(f_dAbsDwc), &
           & got(f_molecules), got(f_velLin) /) )

    nPress = pfaDatum%vGrid%noSurfs
    nTemps = pfaDatum%tGrid%noSurfs
    nArrays = nPress * nTemps + 1

    call allocate_test ( pfaDatum%absorption, nTemps, nPress, 'pfaDatum%absorption', moduleName )
    call allocate_test ( pfaDatum%dAbsDnc,    nTemps, nPress, 'pfaDatum%dAbsDnc',    moduleName )
    call allocate_test ( pfaDatum%dAbsDnu,    nTemps, nPress, 'pfaDatum%dAbsDnu',    moduleName )
    call allocate_test ( pfaDatum%dAbsDwc,    nTemps, nPress, 'pfaDatum%dAbsDwc',    moduleName )

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
            end do
          end if
        else if ( capitalize(fileType) == 'HDF' ) then ! HDF
          call announce_error ( fileIndex, unsupportedFormat, fileType )
        else
          call announce_error ( fileIndex, unsupportedFormat, fileType )
        end if
      end if
    else
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
    PFADatum%vel_cor = 1.0_ck - velLin / c ! Doppler correction factor

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
      use OUTPUT_M, only: OUTPUT
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

! =====     Private Procedures     =====================================

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PFAData_m

! $Log$
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
