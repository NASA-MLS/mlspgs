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
    real(r4), pointer :: LnT(:) => NULL()          ! Log temperatures
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
  subroutine Get_PFAdata_from_l2cf ( Root, Name, VGrids )
  ! Process a PFAdata specification from the l2cf.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Expr_m, only: Expr
    use Init_Tables_Module, only: F_Absorption, F_dAbsDnc, F_dAbsDnu, &
      & F_dAbsDwc, F_Molecules, F_Signal, F_Temperatures, F_VelLin, F_VGrid
    use Intrinsic, only: PHYQ_Dimensionless, PHYQ_Temperature, PHYQ_Velocity
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only: Signals
    use MoreTree, only: Get_Field_ID
    use Parse_Signal_m, only: Parse_Signal
    use String_Table, only: Get_String
    use Tree, only: Decorate, Decoration, NSons, Sub_Rosa, Subtree

    integer, intent(in) :: Root            ! of the pfaData subtree in the l2cf
    integer, intent(in) :: Name            ! of the pfaData spec, else zero
    type(vGrid_t), intent(in), target :: VGrids(:) ! database of vgrids

    ! Error codes
    integer, parameter :: SignalParse = 1
    integer, parameter :: TooManyChannels = signalParse + 1
    integer, parameter :: TooManySignals = tooManyChannels + 1
    integer, parameter :: WrongSize = tooManySignals + 1
    integer, parameter :: WrongUnits = wrongSize + 1

    integer :: AbsTree
    logical, pointer :: Channels(:)
    integer :: dAbsDncTree, dAbsDnuTree, dAbsDwcTree, Dim
    logical :: Error
    integer :: I, J, NArrays, NPress, NTemps, Sideband
    integer, pointer :: SignalIndices(:)
    integer :: Son, TemperatureTree, Units(2)
    type(pfaData_t) :: PFADatum
    double precision :: Value(2)

    error = .false.
    nullify ( channels, signalIndices )
    pfaDatum%name = name
    do i = 2, nsons(root)
      son = subtree(i,root)
      select case ( get_field_id(son) )
      case ( f_absorption )
        absTree = son
      case ( f_dAbsDnc )
        dAbsDncTree = son
      case ( f_dAbsDnu )
        dAbsDnuTree = son
      case ( f_dAbsDwc )
        dAbsDwcTree = son
      case ( f_molecules )
        call allocate_test ( pfaDatum%molecules, nsons(son)-1, &
          & 'pfaDatum%molecules', moduleName )
        do j = 2, nsons(son)
          pfaDatum%molecules(j-1) = decoration(subtree(j,son))
        end do
      case ( f_vGrid )
        pfaDatum%vGrid => vgrids(decoration(decoration(subtree(2,son))))
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
        temperatureTree = son
      case ( f_velLin )
        call expr ( subtree(2,son), units, value )
        if ( units(1) /= phyq_velocity .and. units(1) /= phyq_dimensionless ) &
          & call announce_error ( subtree(1,son), wrongUnits, 'Velocity' )
        pfaDatum%velLin = value(1)
      end select
    end do

    ! Check sizes
    nPress = pfaDatum%vGrid%noSurfs
    nTemps = nsons(temperatureTree) - 1
    nArrays = nPress * nTemps + 1
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

    ! Store temperatures
    call allocate_test ( pfaDatum%lnT, nTemps, 'pfaDatum%lnT', moduleName )
    dim = phyq_dimensionless
    do j = 2, nTemps + 1
      call expr ( subtree(j,temperatureTree), units, value )
      pfaDatum%lnT(j-1) = value(1)
      if ( units(1) /= phyq_dimensionless ) then
        if ( units(1) /= phyq_temperature ) then
          call announce_error ( subtree(j,temperatureTree), wrongUnits, &
            & 'Temperature' )
        else
          dim = phyq_temperature
        end if
      end if
    end do
    if ( dim /= phyq_temperature ) &
      & call announce_error ( subtree(1,temperatureTree), wrongUnits, &
        & 'Temperature' )

    call store_2d ( absTree, pfaDatum%absorption, 'pfaDatum%absorption' )
    call store_2d ( dAbsDncTree, pfaDatum%dAbsDnc, 'pfaDatum%dAbsDnc' )
    call store_2d ( dAbsDnuTree, pfaDatum%dAbsDnu, 'pfaDatum%dAbsDnu' )
    call store_2d ( dAbsDwcTree, pfaDatum%dAbsDwc, 'pfaDatum%dAbsDwc' )

    if ( error ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to complete processing of PFAData' )

    call decorate ( root, addPFADatumToDatabase ( pfaData, pfaDatum ) )

  contains

    ! ...........................................  Announce_Error  .....
    subroutine Announce_Error ( Where, What, String, More )
      use MoreTree, only: StartErrorMessage
      use OUTPUT_M, only: OUTPUT
      integer, intent(in) :: Where             ! Tree node index
      integer, intent(in) :: What              ! Error index
      character(len=*), intent(in), optional :: String  ! For more info
      integer, intent(in), optional :: More    ! For more info
      error = .true.
      call startErrorMessage ( where )
      select case ( what )
      case ( signalParse )
        call output ( 'Unable to parse signal ' )
        call output ( trim(string), advance='yes' )
      case ( tooManyChannels )
        call output ( string )
        call output ( ' Describes more than one channel.', advance='yes' )
      case ( tooManySignals )
        call output ( string )
        call output ( ' Describes more than one signal.', advance='yes' )
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
    subroutine Store_2d ( Where, What, Name )
    ! Store data from Where in the L2CF into an nTemps X nPress
    ! array What
      integer, intent(in) :: Where
      real(r4), pointer :: What(:,:)
      character(len=*), intent(in) :: Name
      integer :: I, J, K

      call allocate_test ( what, nTemps, nPress, name, moduleName )
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
      call deallocate_test ( pfaData(i)%lnT, 'pfaData%lnT', moduleName )
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

    call dump ( pfaDatum%lnT, name=' Ln Temperatures' )

    call output ( ' VGrid: ' )
    call display_string ( pfaDatum%vGrid%name, advance='yes' )

    call output ( pfaDatum%velLin, before=' Velocity Linearization: ', &
      & advance='yes' )

    call dump ( pfaDatum%absorption, name=' ln Absorption' )
    call dump ( pfaDatum%dAbsDwc, name=' d ln Absorption / d wc' )
    call dump ( pfaDatum%dAbsDnc, name=' d ln Absorption / d nc' )
    call dump ( pfaDatum%dAbsDnu, name=' d ln Absorption / d nu' )

  end subroutine Dump_PFADatum

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
! Revision 2.2  2004/05/29 02:51:40  vsnyder
! Allow signal string to denote only one signal
!
! Revision 2.1  2004/05/22 02:29:48  vsnyder
! Initial commit
!
