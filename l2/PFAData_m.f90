! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PFAData_m

  ! Read the PFA data file(s).  Build a database.  Provide for access to it.

  use MLSCommon, only: R4
  use MLSSignals_m, only: Signal_T
  use VGridsDatabase, only: VGrid_t

  implicit NONE
  private
  public :: PFAData_t, PFAData
  public :: Get_PFAdata_from_l2cf
  public :: Destroy_PFADataBase, Dump_PFADataBase, Dump
  public :: Write_PFADataBase, Read_PFADataBase

  interface Dump
    module procedure Dump_PFADatum
  end interface Dump

  type PFAData_t
    integer :: Name                                ! of the pfaData spec
    integer, pointer :: Molecules(:) => NULL()     ! Molecule indices
    character(len=127) :: Signal                   ! The signal string
    integer :: SignalIndex                         ! in Signals database
    type(signal_t) :: TheSignal                    ! The signal, with channels
                                                   ! and sidebands added
    type(vGrid_t), pointer :: TGrid => NULL()      ! Log temperatures
    type(vGrid_t), pointer :: VGrid => NULL()      ! vertical grid
    real(r4) :: VelLin                             ! Velocity linearization, km/s
    real(r4), pointer :: Absorption(:,:) => NULL() ! Ln Absorption data
    real(r4), pointer :: dAbsDwc(:,:) => NULL()    ! d Ln Absorption / d wc data
    real(r4), pointer :: dAbsDnc(:,:) => NULL()    ! d Ln Absorption / d nc data
    real(r4), pointer :: dAbsDnu(:,:) => NULL()    ! d Ln Absorption / d nu data
  end type PFAData_t

  type(PFAData_t), pointer,save :: PFAData(:) => NULL()

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
    use Intrinsic, only: PHYQ_Dimensionless, PHYQ_Temperature, PHYQ_Velocity
    use IO_Stuff, only: Get_Lun
    use MLSSignals_m, only: Signals
    use MLSStrings, only: Capitalize
    use MoreTree, only: Get_Field_ID
    use Parse_Signal_m, only: Parse_Signal
    use String_Table, only: Get_String
    use Tree, only: Decorate, Decoration, Node_Id, NSons, Sub_Rosa, Subtree
    use Tree_Types, only: N_String

    integer, intent(in) :: Root            ! of the pfaData subtree in the l2cf
    integer, intent(in) :: Name            ! of the pfaData spec, else zero
    type(vGrid_t), intent(in), target :: VGrids(:) ! database of vgrids
    integer, intent(out) :: Error          ! 0 => OK, else trouble

    ! Error codes
    integer, parameter :: CannotOpen = 1
    integer, parameter :: CannotRead = cannotOpen + 1
    integer, parameter :: NotZeta = cannotRead + 1
    integer, parameter :: SignalParse = notZeta + 1
    integer, parameter :: TooManyChannels = signalParse + 1
    integer, parameter :: TooManySignals = tooManyChannels + 1
    integer, parameter :: WrongFields = tooManySignals + 1
    integer, parameter :: WrongSize = wrongFields + 1
    integer, parameter :: WrongUnits = wrongSize + 1

    integer :: AbsTree
    logical, pointer :: Channels(:)
    integer :: dAbsDncTree, dAbsDnuTree, dAbsDwcTree, Dim
    integer :: Field, FileIndex ! Where in the tree is the filename?
    character(255) :: FileName, FileType ! Formatted(default), Unformatted
    logical :: Got(field_first:field_last)
    integer :: I, IOStat, J, Lun, NArrays, NPress, NTemps, Sideband
    integer, pointer :: SignalIndices(:)
    integer :: Son, Units(2)
    type(pfaData_t) :: PFADatum
    double precision :: Value(2)

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
        fileType = 'formatted'
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
      case ( f_temperatures )
        pfaDatum%tGrid => vgrids(decoration(decoration(subtree(2,son))))
      case ( f_vGrid )
        pfaDatum%vGrid => vgrids(decoration(decoration(subtree(2,son))))
        if ( pfaDatum%vGrid%verticalCoordinate /= l_zeta ) &
          & call announce_error ( subtree(1,son), notZeta )
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
      case ( f_velLin )
        call expr ( subtree(2,son), units, value )
        if ( units(1) /= phyq_velocity ) &
          & call announce_error ( subtree(1,son), wrongUnits, 'Velocity' )
        pfaDatum%velLin = value(1) / 1000.0 ! fundamental unit is m/s, fwdmdl wants km/s
      end select
    end do

    ! Check that we have file and not absorption etc, or vice-versa
    if ( got(f_file) .and. &
         any( (/ &
           & got(f_absorption), got(f_dAbsDnc), got(f_dAbsDnu), got(f_dAbsDwc) /) ) &
    & .or. .not. got(f_file) .and. .not. &
         all( (/ &
           & got(f_absorption), got(f_dAbsDnc), got(f_dAbsDnu), got(f_dAbsDwc) /)) ) &
      call announce_error ( root, wrongFields )

    nPress = pfaDatum%vGrid%noSurfs
    nTemps = pfaDatum%tGrid%noSurfs
    nArrays = nPress * nTemps + 1

    call allocate_test ( pfaDatum%absorption, nTemps, nPress, 'pfaDatum%absorption', moduleName )
    call allocate_test ( pfaDatum%dAbsDnc,    nTemps, nPress, 'pfaDatum%dAbsDnc',    moduleName )
    call allocate_test ( pfaDatum%dAbsDnu,    nTemps, nPress, 'pfaDatum%dAbsDnu',    moduleName )
    call allocate_test ( pfaDatum%dAbsDwc,    nTemps, nPress, 'pfaDatum%dAbsDwc',    moduleName )

    if ( got(f_file) ) then
      call get_lun ( lun )
      fileName = trim(fileName) // pfaDatum%signal
      open ( unit=lun, file=fileName, form=fileType, status='old', iostat=iostat )
      if ( iostat /= 0 ) then
        call announce_error ( fileIndex, cannotOpen, fileName, iostat )
      else
        if ( capitalize(fileType) /= 'UNFORMATTED' ) then
          read ( lun, *, iostat=iostat ) pfaDatum%absorption, pfaDatum%dAbsDnc, &
            & pfaDatum%dAbsDnu, pfaDatum%dAbsDwc
        else
          read ( lun, iostat=iostat ) pfaDatum%absorption, pfaDatum%dAbsDnc, &
            & pfaDatum%dAbsDnu, pfaDatum%dAbsDwc
        end if
        if ( iostat /= 0 ) &
          & call announce_error ( fileIndex, cannotRead, fileName, iostat )
      end if
    else
      ! Check sizes
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
    end if

    if ( error == 0 ) &
      & call decorate ( root, addPFADatumToDatabase ( pfaData, pfaDatum ) )

  contains

    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( Where, What, String, More )
      use Machine, only: IO_Error
      use MoreTree, only: StartErrorMessage
      use OUTPUT_M, only: OUTPUT
      integer, intent(in) :: Where             ! Tree node index
      integer, intent(in) :: What              ! Error index
      character(len=*), intent(in), optional :: String  ! For more info
      integer, intent(in), optional :: More    ! For more info
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
      case ( notZeta )
        call output ( 'Vertical coordinate for pressure grid must be Zeta.', &
          & advance='yes' )
      case ( signalParse )
        call output ( 'Unable to parse signal ' )
        call output ( trim(string), advance='yes' )
      case ( tooManyChannels )
        call output ( string )
        call output ( ' Describes more than one channel.', advance='yes' )
      case ( tooManySignals )
        call output ( string )
        call output ( ' Describes more than one signal.', advance='yes' )
      case ( wrongFields )
        call output ( 'Need either file and not absorption, dAbsDnc, dAbsDnu or dAbsDwc,', &
          advance='yes' )
        call output ( 'or not file and all of absorption, dAbsDnc, dAbsDnu and dAbsDwc,', &
          advance='yes' )
      case ( wrongSize )
        call output ( 'Incorrect size for ' )
        call output ( trim(string) )
        call output ( more, before=' -- should be ', advance='yes' )
      case ( wrongUnits )
        call output ( 'Incorrect units -- should be ' )
        call output ( trim(string), advance='yes' )
      end select
    end subroutine Announce_Error

    ! .................................................  Store_2d  .....
    subroutine Store_2d ( Where, What )
    ! Store data from Where in the L2CF into an nTemps X nPress
    ! array What
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

  ! ----------------------------------------  Destroy_PFADataBase  -----
  subroutine Destroy_PFADataBase
    use Allocate_Deallocate, only: Deallocate_Test
    integer :: I
    if ( .not. associated(pfaData) ) return
    do i = 1, size(pfaData)
      call deallocate_test ( pfaData(i)%molecules, 'pfaData%molecules', moduleName )
      call deallocate_test ( pfaData(i)%theSignal%channels, 'pfaData...Channels', &
          & moduleName )
      call deallocate_test ( pfaData(i)%absorption, 'pfaData%absorption', moduleName )
      call deallocate_test ( pfaData(i)%dAbsDwc, 'pfaData%dAbsDwc', moduleName )
      call deallocate_test ( pfaData(i)%dAbsDnc, 'pfaData%dAbsDnc', moduleName )
      call deallocate_test ( pfaData(i)%dAbsDnu, 'pfaData%dAbsDnu', moduleName )
    end do
  end subroutine Destroy_PFADataBase

  ! -------------------------------------------  Dump_PFADataBase  -----
  subroutine Dump_PFADataBase
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer :: I
    if ( .not. associated(pfaData) ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Cannot dump unallocated PFA data base" )
    do i = 1, size(pfaData)
      call dump_PFADatum ( pfaData(i) )
    end do
  end subroutine Dump_PFADataBase

  ! ----------------------------------------------  Dump_PFADatum  -----
  subroutine Dump_PFADatum ( PFADatum )

    use Dump_0, only: Dump
    use Intrinsic, only: Lit_Indices
    use MLSSignals_m, only: DisplaySignalName
    use String_Table, only: Display_String, String_Length
    use Output_m, only: Blanks, NewLine, Output

    type(PFAData_t), intent(in) :: PFADatum

    integer :: I, L, W
    character(len=*), parameter :: Molecules = ' Molecules:'

    if ( pfaDatum%name /= 0 ) then
      call output ( ' ' )
      call display_string ( pfaDatum%name )
    end if
    call newLine
    call output ( Molecules )
    w = len(Molecules)
    do i = 1, size(pfaDatum%molecules)
      l = string_length(lit_indices(pfaDatum%molecules(i)))
      if ( w + l > 72 ) then
        call newLine
        w = len(molecules)
        call blanks ( w )
      end if
      call blanks ( 1 )
      call display_string ( lit_indices(pfaDatum%molecules(i)) )
      w = w + l + 1
    end do
    call newline

    call output ( ' Specified signal: ' )
    call output ( trim(pfaDatum%signal), advance='yes' )
    call output ( ' Actual signal: ' )
    call output ( pfaDatum%signalIndex, after=': ' )
    if ( pfaDatum%theSignal%name /= 0 ) then
      call display_string ( pfaDatum%theSignal%name )
      call output ( ': ' )
    end if
    call displaySignalName ( pfaDatum%theSignal, advance='yes' )


    call output ( ' TGrid: ' )
    call display_string ( pfaDatum%tGrid%name )

    call output ( ', VGrid: ' )
    call display_string ( pfaDatum%vGrid%name, advance='yes' )

    call output ( pfaDatum%velLin, before=' Velocity Linearization: ', &
      & advance='yes' )

    call dump ( pfaDatum%absorption, name=' ln Absorption' )
    call dump ( pfaDatum%dAbsDwc, name=' d ln Absorption / d wc' )
    call dump ( pfaDatum%dAbsDnc, name=' d ln Absorption / d nc' )
    call dump ( pfaDatum%dAbsDnu, name=' d ln Absorption / d nu' )

  end subroutine Dump_PFADatum

  ! ------------------------------------------  Write_PFADatabase  -----
  subroutine Write_PFADatabase ( FileName )
    character(len=*), intent(in) :: FileName
  end subroutine Write_PFADatabase

  ! -------------------------------------------  Read_PFADatabase  -----
  subroutine Read_PFADatabase ( FileName )
    character(len=*), intent(in) :: FileName
  end subroutine Read_PFADatabase

! =====     Private Procedures     =====================================

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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PFAData_m

! $Log$
! Revision 2.3  2004/06/08 19:29:27  vsnyder
! Add file field
!
! Revision 2.2  2004/05/29 02:51:40  vsnyder
! Allow signal string to denote only one signal
!
! Revision 2.1  2004/05/22 02:29:48  vsnyder
! Initial commit
!
