! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PFADataBase_m

  ! Read the PFA data file(s).  Build a database.  Provide for access to it.

  use MLSCommon, only: R4
  use MLSSignals_m, only: Signal_T
  use VGridsDatabase, only: VGrid_t

  implicit NONE
  private
  public :: PFAData_t, PFAData
  public :: AddPFADatumToDatabase
  public :: Destroy_PFADataBase, Dump_PFADataBase, Dump
  public :: Write_PFADatum, Write_PFADataBase, Read_PFADataBase

  interface Dump
    module procedure Dump_PFADatum
  end interface Dump

  type PFAData_T
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
  end type PFAData_T

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

  ! -------------------------------------------  Read_PFADatabase  -----
  subroutine Read_PFADatabase ( FileName )
    character(len=*), intent(in) :: FileName
  end subroutine Read_PFADatabase

  ! ------------------------------------------  Write_PFADatabase  -----
  subroutine Write_PFADatabase ( FileName )
    character(len=*), intent(in) :: FileName
  end subroutine Write_PFADatabase

  ! ---------------------------------------------  Write_PDADatum  -----
  subroutine Write_PFADatum ( PFADatum, FileName, FileType )

    ! Write the PFADatum on FileName using the format given by FileType
    ! If FileType is "UNFORMATTED" (case insensitive), the output file
    ! name consists of the part of FileName before "$" (or all of it if
    ! "$" does not appear, followed by the pfaDatum%signal, followed by
    ! the part of FileName after the "$".

    use IO_Stuff, only: Get_Lun
    use Machine, only: IO_Error
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSStrings, only: Capitalize
    use String_Table, only: Get_String, String_Length

    type(PFAData_t), intent(in) :: PFADatum
    character(len=*), intent(in) :: FileName, FileType

    integer :: I, IOSTAT, L, Lun
    character(len=len(fileName)+len_trim(pfaDatum%signal)) :: MyFile
    character(len=31) :: Molecule
    character(len=5) :: What

    if ( capitalize(fileType) == 'UNFORMATTED' ) then
      i = index(fileName,'$')
      if ( i == 0 ) then
        myFile = trim(fileName) // trim(pfaDatum%signal)
      else
        myFile = fileName(:i-1) // trim(pfaDatum%signal) // trim(fileName(i+1:))
      end if
      call get_lun ( lun )
      what = 'open'
      open ( lun, file=trim(myFile), form='unformatted', iostat=iostat, err=9 )
      what = 'write'
      write ( lun, iostat=iostat, err=9 ) pfaDatum%tGrid%noSurfs, pfaDatum%vGrid%noSurfs, &
        & size(pfaDatum%molecules), pfaDatum%velLin, &
        & len_trim(pfaDatum%signal), trim(pfaDatum%signal)
      write ( lun, iostat=iostat, err=9 ) pfaDatum%absorption, pfaDatum%dAbsDwc, &
        & pfaDatum%dAbsDnc, pfaDatum%dAbsDnu
      do i = 1, size(pfaDatum%molecules)
        l = string_length(pfaDatum%molecules(i))
        call get_string ( pfaDatum%molecules(i), molecule )
        write ( lun, iostat=iostat, err=9 ) l, molecule(:l)
      end do
      what = 'close'
      close ( lun, iostat=iostat, err=9 )
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
