! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PFADataBase_m

  ! Read the PFA data file(s).  Build a database.  Provide for access to it.

  use MLSCommon, only: R4
  use MLSSignals_m, only: Signal_T
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
    integer, pointer :: Molecules(:) => NULL()     ! Molecule indices
    character(len=127) :: Signal                   ! The signal string
    integer :: SignalIndex                         ! in Signals database
    type(signal_t) :: TheSignal                    ! The signal, with channels
                                                   ! and sidebands added
    type(vGrid_t), pointer :: TGrid => NULL()      ! Log temperatures
    type(vGrid_t), pointer :: VGrid => NULL()      ! vertical grid
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
      call deallocate_test ( pfaData(i)%molecules, 'pfaData%molecules', moduleName )
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
    use Output_m, only: Output
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
    use String_Table, only: Display_String, String_Length
    use Output_m, only: Blanks, NewLine, Output

    type(PFAData_t), intent(in) :: PFADatum
    integer, intent(in), optional :: Details ! >0 => Dump arrays, 0 => Don't
    integer, intent(in), optional :: Index   ! In PFA Database

    integer, parameter :: CK = kind(speedOfLight)
    real(ck) :: C = speedOfLight / 1000.0_ck ! km/s
    integer :: I, L, W
    character(len=*), parameter :: Molecules = ' Molecules:'
    integer :: MyDetails

    myDetails = 1
    if ( present(details) ) myDetails = details

    call output ( 'PDA Datum ' )
    if ( present(index) ) call output ( index )
    if ( pfaDatum%name /= 0 ) then
      if ( present(index) ) call output ( ': ' )
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
    if ( pfaDatum%filterFile /= 0 ) &
      & call display_string ( pfaDatum%filterFile, before=' Filter file: ', &
      & advance='yes' )

    call output ( ' TGrid: ' )
    call display_string ( pfaDatum%tGrid%name )

    call output ( ', VGrid: ' )
    call display_string ( pfaDatum%vGrid%name )

    call output ( pfaDatum%vel_rel*c, before=', Velocity linearization: ', &
      & after='kms', advance='yes' )

    if ( myDetails <= 0 ) return

    call dump ( pfaDatum%absorption, name=' ln Absorption' )
    call dump ( pfaDatum%dAbsDwc, name=' d ln Absorption / d wc' )
    call dump ( pfaDatum%dAbsDnc, name=' d ln Absorption / d nc' )
    call dump ( pfaDatum%dAbsDnu, name=' d ln Absorption / d nu' )

  end subroutine Dump_PFADatum

  ! ------------------------------------------------- PFADataOrder -----
  pure integer function PFADataOrder ( A, B ) result ( N )
    type (PFAData_t), intent(in) :: A, B
    integer :: I
    do i = 1, min(size(a%molecules),size(b%molecules))
      n = a%molecules(i) - b%molecules(i)
      if ( n /= 0 ) return
    end do
    n = size(a%molecules) - size(b%molecules)
    if ( n /= 0 ) return
    n = a%theSignal%band - b%theSignal%band
    if ( n /= 0 ) return
    n = a%theSignal%instrumentModule - b%theSignal%instrumentModule
    if ( n /= 0 ) return
    n = a%theSignal%radiometer - b%theSignal%radiometer
    if ( n /= 0 ) return
    n = a%theSignal%spectrometer - b%theSignal%spectrometer
    if ( n /= 0 ) return
    n = a%theSignal%sideband - b%theSignal%sideband
  end function PFADataOrder

  ! ------------------------------------------ PFADataOrderIndexed -----
  pure integer function PFADataOrderIndexed ( A, B )
    integer, intent(in) :: A, B
    PFADataOrderIndexed = PFADataOrder ( PFAData(a), PFAData(b) )
  end function PFADataOrderIndexed

  ! -------------------------------------------  Read_PFADatabase  -----
  subroutine Read_PFADatabase ( FileName )
    character(len=*), intent(in) :: FileName
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
    use Sort_m, only: ISORT, GSORTP

    integer :: I, J, N

    if ( .not. associated(PFAData) ) return ! Can't sort it if it's not there

    n = size(PFAData)

    if ( associated(sortPFAData) ) then
      if ( size(sortPFAData) == n ) return ! only do it once
      call deallocate_test ( sortPFAData, 'SortPFAData', moduleName )
    end if

    call allocate_test ( sortPFAData, n, 'SortPFAData', moduleName )

    ! First sort the Molecules field of every PFA datum.
!   do i = 1, n
!     call isort ( pfaData(i)%molecules, 1, size(pfaData(i)%molecules) )
!   end do

    ! Now sort the database
    call gsortp ( PFADataOrderIndexed, n, sortPFAData )

    ! Now that it's sorted with molecule as the major key, make
    ! PFA_By_Molecule.
    j = size(SortPFAData)
    do i = last_molecule, first_molecule, -1
      do while ( j > 0 )
        if ( PFAData(sortPFAData(j))%molecules(1) > i ) then
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
    use HDF5, only: H5FCREATE_F, H5FClose_F, H5F_ACC_TRUNC_F
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSStrings, only: Capitalize
    character(len=*), intent(in) :: FileName, FileType
    integer :: FileID
    integer :: I, IOSTAT
    if ( capitalize(fileType) == 'HDF5' ) then ! open HDF5 file here
      call H5FCreate_F ( trim(fileName), H5F_ACC_TRUNC_F, fileID, &
        & iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open hdf5 PFA file ' // trim(fileName) // ' for output.' )
      do i = 1, size(pfaData)
        call write_PFADatum ( pfadata(i), FileName, FileType, lun=fileID )
      end do
      call H5FClose_F ( fileID, iostat )
      if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & 'Unable to close hdf5 PFA file ' // trim(fileName) // '.' )
    else ! Open file(s) in Write_PFADatum
      do i = 1, size(pfaData)
        call write_PFADatum ( pfadata(i), FileName, FileType, useMolecule=.true. )
      end do 
    end if
  end subroutine Write_PFADatabase

  ! ---------------------------------------------  Write_PDADatum  -----
  subroutine Write_PFADatum ( PFADatum, FileName, FileType, UseMolecule, Lun )

    ! Write the PFADatum on FileName using the format given by FileType
    ! If FileType is "UNFORMATTED" (case insensitive) and UseMolecule is
    ! absent or present but false, the output file name consists of the
    ! part of FileName before "$" (or all of it if "$" does not appear,
    ! followed by the pfaDatum%signal, followed by the part of FileName
    ! after the "$".  If UseMolecule is present and true, instead of the
    ! signal, the embedded part consists of the (first) molecule name,
    ! followed by an underscore, followed by the signal.

    use HDF5, only: H5FCREATE_F, H5FClose_F, H5F_ACC_TRUNC_F, &
      & H5GCLOSE_F, H5GCREATE_F
    use Intrinsic, only: Lit_Indices
    use IO_Stuff, only: Get_Lun
    use Machine, only: IO_Error
    use MLSHDF5, only: MakeHDF5Attribute, SaveAsHDF5DS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSStrings, only: Capitalize
    use Output_m, only: Output
    use Physics, only: SpeedOfLight
    use String_Table, only: Get_String, String_Length
    use Toggles, only: Switches

    type(PFAData_t), intent(in) :: PFADatum
    character(len=*), intent(in) :: FileName, FileType
    logical, intent(in), optional :: UseMolecule
    integer, intent(in), optional :: Lun ! Don't open a new file if present

    character(len=1023) :: Attrib
    integer, parameter :: CK = kind(speedOfLight)
    real(ck) :: C = speedOfLight / 1000.0_ck ! km/s
    logical DoMolecule
    character(len=255) :: FilterFile
    integer :: GroupID
    character(len=len_trim(pfaDatum%signal)+32) :: GroupName
    integer :: I, IOSTAT, L, MyLun
    character(len=len(fileName)+len_trim(pfaDatum%signal)+31) :: MyFile
    character(len=31) :: Molecule
    character(len=31) :: Molecules(size(pfaDatum%molecules))
    character(len=5) :: What

    doMolecule = .false.
    if ( present(useMolecule) ) doMolecule = useMolecule

    if ( doMolecule ) then
      call get_string ( lit_indices(pfaDatum%molecules(1)), molecule )
      molecule = trim(molecule) // '_'
    else
      molecule = ''
    end if

    if ( present(lun) ) myLun = lun
    if ( capitalize(fileType) == 'UNFORMATTED' ) then
      if ( .not. present(lun) ) then
        i = index(fileName,'$')
        if ( i == 0 ) then
          myFile = trim(fileName) // trim(molecule) // trim(pfaDatum%signal)
        else
          myFile = fileName(:i-1) // trim(molecule) // trim(pfaDatum%signal) // trim(fileName(i+1:))
        end if
        call get_lun ( myLun )
        what = 'open'
        open ( myLun, file=trim(myFile), form='unformatted', iostat=iostat, err=9 )
      end if
      what = 'write'
      write ( myLun, iostat=iostat, err=9 ) pfaDatum%tGrid%noSurfs, pfaDatum%vGrid%noSurfs, &
        & size(pfaDatum%molecules), real(pfaDatum%vel_rel*c,rk), &
        & len_trim(pfaDatum%signal), trim(pfaDatum%signal), &
        & real(pfaDatum%vGrid%surfs(1,1)), &
        & real(pfaDatum%vGrid%surfs(2,1)-pfaDatum%vGrid%surfs(1,1)), &
        & real(pfaDatum%tGrid%surfs(1,1)), &
        & real(pfaDatum%tGrid%surfs(2,1)-pfaDatum%tGrid%surfs(1,1))
      write ( myLun, iostat=iostat, err=9 ) pfaDatum%absorption, pfaDatum%dAbsDwc, &
        & pfaDatum%dAbsDnc, pfaDatum%dAbsDnu
      do i = 1, size(pfaDatum%molecules)
        l = string_length(lit_indices(pfaDatum%molecules(i)))
        call get_string ( lit_indices(pfaDatum%molecules(i)), molecule )
        write ( myLun, iostat=iostat, err=9 ) l, molecule(:l)
      end do
      write ( myLun, iostat=iostat, err=9 ) pfaDatum%tGrid%surfs, pfaDatum%vGrid%surfs
      if ( pfaDatum%filterFile /= 0 ) then
        l = string_length(pfaDatum%filterFile)
        call get_string ( pfaDatum%filterFile, filterFile )
        write ( myLun,iostat=iostat, err=9 ) l, filterFile(:l)
      end if
      what = 'close'
      if ( .not. present(lun) ) close ( myLun, iostat=iostat, err=9 )
    else if ( capitalize(fileType) == 'HDF5' ) then
      if ( .not. present(lun) ) then
        ! Open the HDF file
        call H5FCreate_F ( trim(fileName), H5F_ACC_TRUNC_F, myLun, &
          & iostat )
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to open hdf5 PFA file ' // trim(fileName) // ' for output.' )
      end if

      ! Create a group named <molecule>_<signal>
      call get_string ( lit_indices(pfaDatum%molecules(1)), molecule )
      groupName = trim(molecule) // '_' // trim(pfaDatum%signal)
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
      if ( pfaDatum%name /= 0 ) then
        call get_string ( pfaDatum%filterFile, attrib )
        call MakeHDF5Attribute ( groupID, 'filterFile', &
          & attrib(:string_length(pfaDatum%filterFile)) )
      end if
      do i = 1, size(pfaDatum%molecules)
        call get_string ( pfaDatum%molecules(i), molecules(i) )
      end do
      call MakeHDF5Attribute ( groupID, 'numMolecules', size(molecules) )
      call MakeHDF5Attribute ( groupID, 'molecules', molecules )
      call MakeHDF5Attribute ( groupID, 'signal', pfaDatum%signal )
      call MakeHDF5Attribute ( groupID, 'vel_rel', pfaDatum%vel_rel )
      call MakeHDF5Attribute ( groupID, 'nTemps', pfaDatum%tGrid%noSurfs )
      call MakeHDF5Attribute ( groupID, 'tStart', pfaDatum%tGrid%surfs(1,1) )
      call MakeHDF5Attribute ( groupID, 'tStep', pfaDatum%tGrid%surfs(2,1)-pfaDatum%tGrid%surfs(1,1) )
      call MakeHDF5Attribute ( groupID, 'nPress', pfaDatum%vGrid%noSurfs )
      call MakeHDF5Attribute ( groupID, 'vStart', pfaDatum%vGrid%surfs(1,1) )
      call MakeHDF5Attribute ( groupID, 'vStep', pfaDatum%vGrid%surfs(2,1)-pfaDatum%vGrid%surfs(1,1) )
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
        if ( iostat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
          & 'Unable to close hdf5 PFA file ' // trim(fileName) // '.' )
      end if
    else
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unsupported file format in Write_PFADatum' )
    end if

    return

9   call io_error ( 'Unable to ' // trim(what) // ' output file ', iostat, myFile )
    call MLSMessage ( MLSMSG_Error, moduleName, 'Execution terminated' )

  end subroutine Write_PFADatum

! =====     Private Procedures     =====================================

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PFADataBase_m

! $Log$
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
