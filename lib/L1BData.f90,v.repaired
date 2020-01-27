! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module L1BData

  ! Reading and interacting with Level 1B data (HDF4 or HDF5)

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, &
    & Test_Allocate, Test_DeAllocate
  use HyperSlabs, only: EssentiallyEqual
  use Diff_1, only: Diff, Diff_Fun
  use Dump_Options, only: NameOnEachLine, StatsOnOneLine
  use Dump_0, only: Dump
  use Dump_1, only: Dump
  use HDF, only: DFACC_Rdonly, SFGInfo, SFN2Index, SFSelect, &
    & SFRData_F90, &
    & SFRCData, SFENDACC, DFNT_Char8, DFNT_Int32, DFNT_Float64, &
    & DFNT_Float32
  use HighOutput, only: OutputNamedValue
  use IEEE_Arithmetic, only: IEEE_Is_Finite
  use, Intrinsic :: ISO_C_Binding, only: C_Intptr_T, C_Loc
  use Intrinsic, only: L_HDF
  use Lexer_Core, only: Print_Source
  use Machine, only: Crash_Burn
  use MLSCommon, only: MLSFile_T, &
    & UnDefinedValue, FileNameLen, NameLen
  use MLSFiles, only: FileNotFound, HDFVersion_4, HDFVersion_5, &
    & AddFileToDatabase, GetMLSFileByType, InitializeMLSfile, &
    & MLS_OpenFile, MLS_CloseFile
  use MLSFinds, only: FindFirst
  use MLSHDF5, only: MaxNDSNames, GetAllHDF5DSNames, SaveAsHDF5DS
  use MLSKinds, only: R4, R8
  use MLSMessageModule, only: MLSMSG_Error, &
    & MLSMSG_L1BRead, MLSMSG_Warning, MLSMessage
  use MLSStrings, only: Indexes, Reverse, Streq
  use MLSStringLists, only: NumStringElements, ReplaceSubstring, SwitchDetail
  use Moretree, only: Get_Field_Id
  use Output_M, only: NewLine, Output
  use String_Table, only: Get_String
  use Toggles, only: Gen, Levels, Switches, Toggle
  use Trace_M, only: Trace_Begin, Trace_End
  use Tree, only: Nsons, Sub_Rosa, Subtree, Dump_Tree_Node, Where

  implicit none

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (data types and parameters)
! L1BData_T                       Quantities from an L1B data file
! PRECISIONSUFFIX                 Suffix stuck on end of array name if precision
! NOERROR                         ReadL1BData had no error
! NOCOUNTERMAFINDX                ReadL1BData unable to find counterMAF index
! NOCOUNTERMAFID                  ReadL1BData unable to find counterMAF sds
! NOQUANTITYINDEX                 ReadL1BData unable to find quantity index
! NODATASETID                     ReadL1BData unable to find quantity sds
! NODATASETRANK                   ReadL1BData unable to find quantity rank
! FIRSTMAFNOTFOUND                ReadL1BData unable to requested first MAF
! LASTMAFNOTFOUND                 ReadL1BData unable to requested last MAF
! CANTREADCOUNTERMAF              ReadL1BData unable to read counterMAF
! CANTALLOCATECHARS               ReadL1BData unable to allocate memory
! CANTREAD3DFIELD                 ReadL1BData unable to read sds
! UNKNOWNDATATYPE                 ReadL1BData unable to read data type
! CANTENDCOUNTERMAF               ReadL1BData unable to end access to counterMAF
! CANTENDQUANTITY                 ReadL1BData unable to end access to quantity

!     (subroutines and functions)
! AssembleL1BQtyName              Returns Quantity Name depending on hdfVersion
! ContractL1BData                 Contract an l1bData to just a range of mafs
! ConvertL1BData                  Convert an l1bData from integer-valued to d.p.
! ConvertL1BOADSNames             Convert L1BOA DS names between hdf versions
! CpL1BData                       Duplicate an l1bData
! DeallocateL1BData               Called when an l1bData is finished with
! DupL1BData                      Duplicate an l1bData
! Diff                            Diff two l1b quantities
! Dump                            Print facts about l1brad quantity
! Findl1bdata                     Which file handle contains a given sd name
! FindMaxMAF                      What is the maximum MAF number in the file?
! Getl1BFile                      Which MLSFile contains a given sd name
!                                  (or attribute)
! L1BOAHDF4DSName                 Convert L1BOA DS names from HDF5 to HDF4
! L1BOAHDF4DSName                 Convert L1BOA DS names from HDF4 to HDF5
! L1BOASetup                      From l2cf, open, and save l1boa file info
! L1BRadSetup                     From l2cf, open, and save l1brad file info
! ReadL1BData                     Read all info concerning a l1brad quantity
! ReadL1BAttribute                Read attributes
! === (end of toc) ===

! === (start of api) ===
!     (user-defined types)
! L1BData_T (char L1BName, int FirstMAF, int NoMAFs, int MaxMIFs, int NoAuxInds,
!             *int CounterMAF(:),
!             *char CharField(:,:,:),
!             *r8   DpField(:,:,:),
!             *int  IntField(:,:,:),
!             )

!     (subroutines and functions)
! char*64 AssembleL1BQtyName (char name, int hdfVersion, log isTngtQty,
!                         [char InstrumentName] )
! ConvertL1BData ( l1bData_T l1bDatain , l1bData_T l1bDataOut, int noMAFs,
!   int firstMAF, int lastMAF, int firstChannel, int lastChannel,
!   int firstMIF, int lastMIF )
! ConvertL1BData ( l1bData_T l1bData )
! ConvertL1BOADSNames ( char* HDF4name, char* HDF5Name, char* direction )
! CpL1BData ( MLSFile_t L1BFile1, MLSFile_t L1BFile2, char QuantityName, &
!    [char NewName], [log l2aux] )
! DupL1BData ( l1bData_T l1bData1, l1bData_T l1bData2, int offsetMAF )
! DeallocateL1BData (l1bData_T l1bData)
! Diff (l1bData_T l1bData1, l1bData_T l1bData2, int details, &
!   [char options], [int numDiffs], [int mafStart], [int mafEnd])
! Dump (l1bData_T l1bData, int details)
! int FindL1BData (int files(:), char fieldName, [int hdfVersion])
! int FindMaxMAF (MLSFile_t file, [int minMAF])
! MLSFile_T GetL1BFile (MLSFile_t filedatabase(:), char fieldName, [char options])
! char* L1BOAHDF4DSName ( char* HDF5Name )
! char* L1BOAHDF5DSName ( char* HDF4Name )
! L1boaSetup (int root, MLSFile_t filedatabase(:), int f_file, [int hdfVersion])
! L1bradSetup (int root, MLSFile_t filedatabase(:), int f_file, [int hdfVersion])
! ReadL1BData (int L1FileHandle, char QuantityName, l1bData_T l1bData,
!               int NoMAFs, int Flag, [int FirstMAF], [int LastMAF],
!               [log NeverFail], [int hdfVersion])
! ReadL1BAttribute (int L1FileHandle, value(:), nchar AttributeName,
!               int Flag, [int hdfVersion])
! === (end of api) ===

  private

  public :: L1BData_T, NameLen, PrecisionSuffix, &
    & AllocateL1BData, AssembleL1BQtyName, &
    & CheckForCorruptFileDatabase, ContractL1BData, &
    & ConvertL1BData, ConvertL1BOADSNames, CpL1BData, &
    & DeallocateL1BData, Diff, Dump, DupL1BData, &
    & FindL1BData, FindMaxMAF, GetL1BFile, &
    & L1BOAHDF4DSName, L1BOAHDF5DSName, L1BRadSetup, L1BOASetup, PadL1BData, &
    & ReadL1BAttribute, ReadL1BData

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  interface Diff
    module procedure DiffL1BData
  end interface

  interface Dump
    module procedure DumpL1BData
  end interface

  interface FindMaxMAF
    module procedure FindMaxMAF_arr, FindMaxMAF_sca
  end interface

  interface ReadL1BAttribute
    module procedure ReadL1BAttribute_intarr1, &
      & ReadL1BAttribute_dblarr1
  end interface

  interface ReadL1BData
    module procedure ReadL1BData_FileHandle
    module procedure ReadL1BData_MLSFile
  end interface

  ! Parameters
  ! suffix of sd precision; check against 'grep -i precision l1/OutputL1B.f90'
  character  (len=*), parameter :: PRECISIONSUFFIX = ' precision'

  ! If TRUE, treats l1brad and l2aux files the same
  ! Among other things, this means forgiving absence of counterMAF
  ! so we may read l2aux and l1brad files alike
  logical, parameter            :: JUSTLIKEL2AUX = .true.

  ! Assume l1b files w/o explicit hdfVersion field are this
  ! 4 corresponds to hdf4, 5 to hdf5 in L2GP, L2AUX, etc.
  integer, parameter :: L1BDEFAULT_HDFVERSION = HDFVERSION_5

  ! Unless the following is true, a subgroup named 'tp' will be
  ! created under each instrument Module for hdf5 versions
  ! (An early idea we later repented of)
  logical, parameter :: DROPTPSUBGROUP = .true.

  ! No MAF number can ever be this big
  integer, parameter :: BIGGESTMAFCTR = huge(0)/2

  ! This data type is used to store quantities from an L1B data file.
  ! (This next param may need to be increased)
  integer, parameter :: MaxCharFieldLen = 128  ! max char field length

  type L1BData_T
    character (len=namelen) :: L1BName ! Name of field in file
    character (len=16) :: data_type     ! 'character', 'double', or 'integer'
    integer :: FirstMAF                 ! First major frame read (usu. 0)
    integer :: FirstMAFCtr              ! depends on date
    integer :: LastMAFCtr               ! depends on date
    integer :: NoMAFs                   ! # of MAFs read
    integer :: MaxMIFs                  ! Max # of MIFs/MAF in SD array
    integer :: NoAuxInds                ! # of auxilliary indices
    integer :: TrueRank                 ! # necessary indices (e.g., 1 for MAFs)
    character (len=16) :: NameInst = ' '

    integer, dimension(:), pointer :: CounterMAF => NULL() ! dimensioned (noMAFs)
    character(len=MaxCharFieldLen), &
      &        dimension(:,:,:), pointer :: CharField => NULL()
    real(r8),  dimension(:,:,:), pointer :: DpField => NULL()
    integer,   dimension(:,:,:), pointer :: IntField => NULL()
    ! all the above dimensioned (noAuxInds,maxMIFs,noMAFs)
    ! logical :: mustPad                  ! Gaps in counterMAF
  contains
    final :: DeallocateL1BData
  end type L1BData_T

  integer :: Error            ! Error level -- 0 = OK
  logical, parameter           :: DEEBUG = .FALSE.

  ! Error flags returned from ReadL1BData when NeverFail=TRUE
  integer, public, parameter :: NOERROR =            0
  integer, public, parameter :: NOCOUNTERMAFINDX =   NOERROR + 1
  integer, public, parameter :: NOCOUNTERMAFID =     NOCOUNTERMAFINDX + 1
  integer, public, parameter :: NOQUANTITYINDEX =    NOCOUNTERMAFID + 1
  integer, public, parameter :: NODATASETID =        NOQUANTITYINDEX + 1
  integer, public, parameter :: NODATASETRANK =      NODATASETID + 1
  integer, public, parameter :: FIRSTMAFNOTFOUND =   NODATASETRANK + 1
  integer, public, parameter :: LASTMAFNOTFOUND =    FIRSTMAFNOTFOUND + 1
  integer, public, parameter :: CANTREADCOUNTERMAF = LASTMAFNOTFOUND + 1
  integer, public, parameter :: CANTREADSCALARDS   = CANTREADCOUNTERMAF + 1
  integer, public, parameter :: CANTALLOCATECHARS =  CANTREADSCALARDS + 1
  integer, public, parameter :: CANTREAD3DFIELD =    CANTALLOCATECHARS + 1
  integer, public, parameter :: UNKNOWNDATATYPE =    CANTREAD3DFIELD + 1
  integer, public, parameter :: CANTENDCOUNTERMAF =  UNKNOWNDATATYPE + 1
  integer, public, parameter :: CANTENDQUANTITY =    CANTENDCOUNTERMAF + 1

  ! We plan to use this 2-d array to convert hdf5-formatted l1boa ds names
  ! to the older hdf4 format.
  include 'l1boa_dsnames.f9h'
contains ! ============================ MODULE PROCEDURES =======================

  ! --------------------------------------------  AllocateL1BData  -----
  subroutine AllocateL1BData ( l1bData, &
    & indims, trueRank, datatype, &
    & L1bDataSibling )
    ! Allocate arrays in l1bData type to match one of
    !   indims(1:3)       explicit sizes
    !   L1bDataSibling    L1bDataSibling%dims(1:3)
    ! Their values will be initialized as undefined value or as ' '
    type( L1BData_T ), intent(out) :: L1bData
    integer, dimension(3), optional,  intent(in) :: indims
    integer, optional, intent(in) :: trueRank
    character(len=*), optional, intent(in) :: datatype
    type( L1BData_T ), optional, intent(in) :: L1bDataSibling
    ! Internal variables
    integer :: noMAFs
    character(len=16) :: myDataType
    integer :: rank
    integer, dimension(3) :: dims

    ! Executable code
    if ( present(datatype) ) then
      myDataType = datatype
      rank = trueRank
      dims = indims
      noMAFs = dims(3)
    else if( present(L1bDataSibling) ) then
      myDataType = L1bDataSibling%data_type
      rank = L1bDataSibling%trueRank
      noMAFs = L1bDataSibling%noMAFs
      if ( associated(L1bDataSibling%charField)) then
        dims = shape(L1bDataSibling%charField)
      else if ( associated(L1bDataSibling%intField)) then
        dims = shape(L1bDataSibling%intField)
      else if ( associated(L1bDataSibling%dpField)) then
        dims = shape(L1bDataSibling%dpField)
      else
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'allocatel1bData was passed L1bDataSibling w/o allocating it' )
      end if
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'allocatel1bData must be passed either datatype or L1bDataSibling' )
    end if
    if ( present(indims) ) then
      dims = indims
      noMAFs = dims(3)
    end if
    ! print *, 'rank in allocatel1bdata ', rank
    ! print *, 'dims in allocatel1bdata ', dims
    ! print *, 'noMAFs in allocatel1bdata ', noMAFs
    if ( any(dims < 1) ) then
      call output('NoMAFs ', advance='no')
      call output(NoMAFs , advance='yes')
      call output('dims ', advance='no')
      call output(dims , advance='yes')
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'allocatel1bData encountered illegal dims' )
    end if
    call allocate_test ( l1bData%counterMAF, noMAFs, 'l1bData%counterMAF', &
      & ModuleName )
    l1bData%counterMAF = int(undefinedValue)
    l1bdata%data_type = myDataType
    l1bdata%noMAFs = noMAFs
    select case( myDataType(1:1) )
    case ( 'c' ) ! character
      call allocate_test ( l1bData%charField, dims(1), dims(2), dims(3), &
        & 'l1bData%charField', ModuleName )
      l1bData%charField = ''
    case ( 'i' ) ! integer
      call allocate_test ( l1bData%intField, dims(1), dims(2), dims(3), &
        & 'l1bData%intField', ModuleName )
      l1bData%intField = int(undefinedValue)
    case ( 'd', 'r' ) ! double
      call allocate_test ( l1bData%dpField, dims(1), dims(2), dims(3), &
        & 'l1bData%dpField', ModuleName )
      l1bData%dpField = undefinedValue
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'datatype not recognized in allocating l1bdata field' )
    end select
    if ( .not. present(L1bDataSibling) ) return
    l1bData%L1BName                = l1bDataSibling%L1BName
    l1bData%FirstMAF               = l1bDataSibling%firstMAF
    l1bData%FirstMAFCtr            = l1bDataSibling%firstMAFCtr
    l1bData%LastMAFCtr             = l1bDataSibling%LastMAFCtr
    l1bData%MaxMIFs                = l1bDataSibling%MaxMIFs
    l1bData%NoAuxInds              = l1bDataSibling%NoAuxInds
    l1bData%TrueRank               = l1bDataSibling%TrueRank
    l1bData%NameInst               = l1bDataSibling%NameInst
  end subroutine AllocateL1BData

  ! -----------------------------------------  AssembleL1BQtyName  -----
  function AssembleL1BQtyName ( name, hdfVersion, isTngtQty, &
    & InstrumentName ) &
    & result(QtyName)
    use MLSStrings, only: IsDigits
    ! Returns a QtyName to be found in the L1b file
    ! If given InstrumentName, name should be a fragment:
    ! e.g., name='VelECI' and InstrumentName='sc'
    ! Without InstrumentName, name should be complete (in hdf4-style)
    ! e.g., name='scVelECI' or name='R1A:118.B1F:PT.S0.FB25-1'
    ! and isTngtQty is simply ignored
    !
    ! we should rewrite it using the L1BOADSNames array.
    character(len=*), intent(in)  :: name        ! bare name; e.g. SolarZenith
    integer, intent(in)           :: hdfVersion  ! which QtyName must conform to
    logical, intent(in)           :: isTngtQty   ! T or F (ignored)
    character(len=*), intent(in), optional :: InstrumentName ! e.g. THz
    ! logical, intent(in), optional :: dont_compress_name ! Won't use
    character(len=Namelen)        :: QtyName

    ! Private
    character(len=Namelen)        :: HDF4Name
    character(len=Namelen)        :: HDF5Name
    integer                       :: indx
    ! Executable
    ! Are we a radiance?
    if ( name(1:1) == 'R' .and. IsDigits(name(2:2) ) ) then
      QtyName = Name
      return
    elseif ( name == 'counterMAF' ) then
      QtyName = Name
      return
    endif
    ! Apparently we are an l1boa quantity
    ! Were we given an instrument name?
    if ( present(InstrumentName) ) then
      HDF5Name = trim(InstrumentName) // '/' // Name
    else
      HDF5Name = Name
    endif
    ! print *, 'HDF5Name: ', trim(HDF5Name)
    if ( hdfVersion == HDFVersion_4 ) then
      QtyName = L1BOAHDF4DSName( HDF5Name )
      ! print *, 'QtyName: ', trim(QtyName)
      if ( QtyName == 'unknown' ) QtyName = Name
    else
      HDF4Name = Name
      QtyName = L1BOAHDF5DSName( HDF4Name )
      if ( QtyName == 'unknown' ) QtyName = HDF5Name
    endif
  end function AssembleL1BQtyName

  function AssembleL1BQtyName_old ( name, hdfVersion, isTngtQty, &
    & InstrumentName, dont_compress_name ) &
    & result(QtyName)
    use MLSStrings, only: CompressString, Lowercase
    ! Returns a QtyName to be found in the L1b file
    ! If given InstrumentName, name should be a fragment:
    ! e.g., name='VelECI' and InstrumentName='sc'
    ! Without InstrumentName, name should be complete (in hdf4-style)
    ! e.g., name='scVelECI' or name='R1A:118.B1F:PT.S0.FB25-1'
    ! and isTngtQty is simply ignored
    !
    ! This is flaky and difficult to maintain as currently written;
    ! we should rewrite it using the L1BOADSNames array.
    character(len=*), intent(in) :: name        ! bare name; e.g. SolarZenith
    integer, intent(in)          :: hdfVersion  ! which QtyName must conform to
    logical, intent(in)          :: isTngtQty   ! T or F
    character(len=*), intent(in), optional :: InstrumentName ! e.g. THz
    logical, intent(in), optional :: dont_compress_name
    character(len=namelen)       :: QtyName

    ! Private
    character(len=1), dimension(HDFVersion_4:HDFVersion_5), parameter :: &
      heads =       (/ ' ', '/' /), &
      instr_tails = (/ '.', '/' /), &
      tp_tails =    (/ ' ', '/' /)
    logical, parameter     :: DEEBUG = .false.
    character(len=1)       :: head
    character(len=1)       :: instr_tail
    character(len=1)       :: tp_tail
    character(len=16)      :: my_instrument
    integer                :: p
    character(len=namelen) :: the_rest
    logical                :: is_a_signal
    logical                :: compress
    ! Executable code
    compress = .true.
    if ( present(dont_compress_name) ) compress = .not. dont_compress_name
    is_a_signal = (NumStringElements(trim(name), .true., ':') > 2)
    if ( DEEBUG ) then
      print *, 'name: ', trim(name)
      print *, 'hdfVersion: ', hdfVersion
      print *, 'isTngtQty: ', isTngtQty
      print *, 'is a signal: ', is_a_signal
      if ( present(InstrumentName) ) print *, 'Instrument: ', trim(InstrumentName)
    end if
    my_instrument = ''
    if ( present(InstrumentName) ) then
      my_instrument = instrumentName
      if ( my_instrument == 'SC' ) my_instrument = 'sc'
    end if
    head = heads(hdfVersion)
    instr_tail = instr_tails(hdfVersion)
    tp_tail = tp_tails(hdfVersion)
    if ( hdfVersion == HDFVERSION_4 .and. my_instrument == 'sc' ) instr_tail= ''
    if ( DEEBUG ) then
      print *, 'my_instrument: ', trim(my_instrument)
    endif
    if ( present(InstrumentName) ) then
      QtyName = head // trim(my_instrument) // instr_tail
    else if ( is_a_signal .or. hdfVersion /= HDFVERSION_5 ) then
      QtyName = head
    else
      ! Need only to convert complete hdf4-name to hdf5-name
      ! This means we must parse hdf4-name fully, however
      QtyName = head
      ! Do we match the pattern 'Instrument.tpItem'
      if ( streq(name, '*.tp*', options='-wc' ) ) then
        p = index( name, '.tp' )
        my_instrument = name(:p-1)
        the_rest = name(p+3:)
      ! Is there an instrument prefixed to name?
      elseif ( name(1:2) == 'sc' ) then
        my_instrument = 'sc'
        the_rest = name(3:)
      else
        my_instrument = ' '
        the_rest = name
      end if
      ! 2nd--Does the_rest start with 'tp'?
      if ( my_instrument == ' ' ) then
        QtyName = '/' // trim(name)
      else if ( the_rest(1:2) /= 'tp' ) then
        QtyName = '/' // trim(my_instrument) // &
          & '/' // trim(the_rest)
      else if ( DROPTPSUBGROUP .and. HDFVersion /= HDFVersion_4 ) then
        QtyName = '/' // trim(my_instrument) // &
          & '/' // trim(the_rest(3:))
      else
        QtyName = '/' // trim(my_instrument) // &
          & '/' // 'tp' // &
          & '/' // trim(the_rest(3:))
      end if
      ! We're done
      if ( DEEBUG ) then
        print *, 'converted name: ', trim(QtyName)
      end if
      call restoreHDF4Dot ( QtyName )
      return
    end if
    if ( isTngtQty .and. Lowercase(my_Instrument) /= 'sc' ) then
      if ( DROPTPSUBGROUP .and. hdfVersion == HDFVERSION_5 ) then
        QtyName = trim(QtyName)
      else
        QtyName = trim(QtyName) // 'tp' // tp_tail
      end if
    end if
    QtyName = trim(QtyName) // trim(name)
    if ( compress ) QtyName = CompressString(QtyName)
    if ( DEEBUG ) then
      print *, 'Before restoring dot: ', trim(QtyName)
    end if
    call restoreHDF4Dot ( QtyName )
    ! No name should begin with a dot
    if ( QtyName(1:1) == '.' ) QtyName = QtyName(2:)
    if ( DEEBUG ) then
      print *, 'more converted name: ', trim(QtyName)
    end if
  contains
    subroutine restoreHDF4Dot ( name )
      ! If hdf4, make certain the "." is part of the dataset name, e.g. 
      !               GHz.tpGeodAngle
      ! If hdf5, make certain "/" is not followed by "tp", e.g. 
      !               not /GHz/tpGeodAngle but /GHz/GeodAngle
      ! Dummy arg
      character(len=*), intent(inout)     :: name
      ! Internal variable
      integer                             :: i
      character(len=len(name))            :: temp
      ! Executable
      if ( hdfVersion == HDFVERSION_5 ) then
        i = index(name, '/tp')
        if ( i > 0 ) name = name(:i) // name(i+3:)
      else
        if ( index(name, 'tptp' ) > 0 ) then
          temp = name
          call ReplaceSubString ( temp, name, 'tptp', 'tp' )
        endif
        if ( index(name, '.tp') > 0 ) return
        if ( index(name, 'tp') < 1 ) return
        ! OK, we apparently have an undotted tp in name
        ! so we'll simply insert a "." before it
        ! call output( 'Restored the dot in ' // trim(name), advance='yes' )
        i = index( name, 'tp' )
        name = name(:i-1) // '.' // name(i:)
      endif
    end subroutine restoreHDF4Dot
  end function AssembleL1BQtyName_old

  ! --------------------------------------------  CheckForCorruptFileDatabase  -----
  subroutine CheckForCorruptFileDatabase ( FileDatabase )
    ! Check whether filedatabase has become corrupted
    type (MLSFile_T), dimension(:), pointer, optional ::     FILEDATABASE
    ! Internal variables
    character(len=1024) :: DSNames
    type (MLSFile_T), pointer :: L1BFile

    ! Executable code
    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    call output ( 'Checking file ' // trim(L1BFile%name) // ' for corrupt database', &
                & advance='yes' )
    if ( L1BFile%HDFVersion == HDFVersion_5 ) then
      call GetAllHDF5DSNames ( L1BFile, DSNames )
      call Dump ( DSNames, 'DSNames from MLSFile type' )
      call newLine
      !call GetAllHDF5DSNames ( L1BFile%name, '/', DSNames )
      !call Dump ( DSNames, 'DSNames from MLSFile name' )
      !call newLine
      ! call crash_burn
    else
      call output ( 'Unable to check HDF4 file', advance='yes' )
    end if
  end subroutine CheckForCorruptFileDatabase

  ! --------------------------------------------  ContractL1BData  -----
  subroutine ContractL1BData ( L1BDataIn, L1BDataOut, noMAFs, &
    & firstMAF, lastMAF, firstChannel, lastChannel, firstMIF, lastMIF )
    ! Contract full l1bdataIn to just those mafs
    ! beginning with firstMAF (if present) and ending with lastMAF
    ! Remember 1st maf of full l1bData is numbered 0
    ! (I disapprove, but nobody ever asks my opinion)

    ! Perhaps related to this source of confusion, is the following:
    ! What do we do when the lastMAF blows past the end of l1bDataIn?
    ! For now, we'll try to keep noMAFs correct, but scoot everything down
    ! Dummy arguments
    type(L1BData_T), intent(in)  :: L1BDataIn
    type(L1BData_T), intent(out) :: L1BDataOut
    integer, optional, intent(in):: firstMAF ! Remember, 1st of full set is 0
    integer, optional, intent(in):: lastMAF
    integer, optional, intent(in):: firstChannel
    integer, optional, intent(in):: lastChannel
    integer, optional, intent(in):: firstMIF
    integer, optional, intent(in):: lastMIF
    integer, intent(out)         :: NoMAFs
    ! Internal variables
    integer, dimension(3) :: dims
    integer :: maf
    integer :: mafOffSet
    integer :: myFirstMAF
    integer :: myLastMAF
    integer :: rank
    logical, parameter :: DEEBug = .false.
    ! Executable
    myFirstMAF = 0
    if ( present(firstMAF) ) then
      call output( 'Should not print, but who knows?', advance='yes' )
      myFirstMAF = firstMAF
    end if
    myLastMAF = l1bDataIn%NoMAFs - l1bDataIn%firstMAF - 1
    if ( present(lastMAF) ) myLastMAF = lastMAF
    if ( associated(l1bDataIn%charField)) then
      dims = shape(l1bDataIn%charField)
    else if ( associated(l1bDataIn%intField)) then
      dims = shape(l1bDataIn%intField)
    else if ( associated(l1bDataIn%dpField)) then
      dims = shape(l1bDataIn%dpField)
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Contractl1bData was passed l1bDataIn w/o allocating it' )
    end if
    noMAFs = myLastMAF - myFirstMAF + 1
    dims(3) = noMAFs
    rank = l1bDataIn%trueRank
    l1bDataOut%noMAFs = noMAFs
    l1bDataOut%firstMAF = myFirstMAF
    if ( present(firstChannel) .and. present(lastChannel) ) &
      & dims(1) = lastChannel - firstChannel + 1
    if ( present(firstMIF) .and. present(lastMIF) ) &
      & dims(1) = lastMIF - firstMIF + 1
    if ( DEEBug ) print *, 'Preparing to contract ', trim(l1bDataIn%L1BName)
    if ( DEEBug ) print *, 'rank ', rank
    if ( DEEBug ) print *, 'NoMAFs ', NoMAFs
    if ( DEEBug ) print *, 'l1b(in)NoMAFs ', l1bDataIn%NoMAFs
    if ( DEEBug ) print *, 'dims ', dims
    if ( DEEBug ) print *, 'min(counterMAF) ', myminval(l1bDataIn%counterMAF)
    if ( DEEBug ) print *, 'min(dpField) ',minval(l1bDataIn%dpField(1,1,:))
    if ( DEEBug ) print *, 'max(dpField) ',maxval(l1bDataIn%dpField(1,1,:))
    call allocateL1BData ( l1bDataOut, dims, L1bDataSibling=l1bDataIn )
    mafOffSet = myFirstMAF
    if ( mafOffSet+NoMAFs > size(l1bDataIn%counterMAF) ) &
      & mafOffSet = size(l1bDataIn%counterMAF) - NoMAFs
    do maf=1, NoMAFs
      l1bDataOut%counterMAF(maf) = l1bDataIn%counterMAF(mafOffSet+maf)
      call cpField(l1bDataIn, mafOffSet+maf, l1bdataOut, maf)
    end do
    l1bDataOut%firstMAFCtr = myminval(l1bDataOut%counterMAF)
    l1bDataOut%lastMAFCtr = maxval(l1bDataOut%counterMAF)
  end subroutine ContractL1BData

  ! --------------------------------------------  ConvertL1BData  -----
  subroutine ConvertL1BData ( L1BData )
    ! Convert an l1bdatafrom integer type to d.p.
    type(L1BData_T), intent(inout)  :: L1BData
    ! Internal variables
    integer, dimension(3) :: dims
    integer :: rank
    logical, parameter :: DEEBug = .false.
    ! Executable
    dims = shape(l1bData%intField)
    call allocate_test ( l1bData%dpField, dims(1), dims(2), dims(3), &
      & 'l1bData%dpField', ModuleName )
    l1bData%dpField = l1bData%intField
    call deallocate_test ( l1bData%intField, 'l1bData%intField', ModuleName )
  end subroutine ConvertL1BData

  ! --------------------------------------------  ConvertL1BOADSNames  -----
  subroutine ConvertL1BOADSNames ( HDF4name, HDF5Name, direction )
    ! Convert an l1boa name from how it appears in an hdf5-formatted file
    ! to the name we used in the older hdf4 format or vice versa.
    ! Useful only in retrievals based on sids radiances
    ! direction             conversion
    ! ---------             ----------
    !   5to4            HDF5Name -> HDF4Name
    !   4to5            HDF4Name -> HDF5Name
    character(len=*), intent(inout)     :: HDF5name
    character(len=*), intent(inout)     :: HDF4name
    character(len=*), intent(in)        :: direction ! '5to4' or '4to5'
    ! Internal variables
    character                           :: goal
    integer                             :: i
    ! Executable
    goal = Reverse(trim(direction)) ! '4' or '5'
    select case( goal )
    case ( '4' )
      HDF4Name = 'unknown'
      i = FindFirst( L1BOADSNames(2,:), HDF5Name )
      ! print *, 'HDF5Name, i', trim(HDF5Name), i
      if ( i > 0 ) HDF4Name = L1BOADSNames(1,i)
    case ( '5' )
      HDF5Name = 'unknown'
      i = FindFirst( L1BOADSNames(1,:), HDF4Name )
      if ( i > 0 ) HDF5Name = L1BOADSNames(2,i)
    case default
    end select
  end subroutine ConvertL1BOADSNames

  ! ----------------------------------------  CpL1BData  -----
  subroutine CpL1BData ( L1BFile1, L1BFile2, QuantityName, NewName, l2aux )

    ! Dummy arguments
    type(MLSFile_T), pointer                :: L1BFile1
    type(MLSFile_T), pointer                :: L1BFile2
    character(len=*), intent(in)            :: quantityName ! Name of SD to cp
    character(len=*), optional, intent(in)  :: newName
    logical, optional, intent(in)           :: l2aux

    ! Local variables
    type(l1bdata_t)                :: L1BDATA ! Result
    integer                        :: Me = -1 ! String index for trace
    integer                        :: flag
    integer                        :: noMAFs
    character(len=256)             :: rename
    integer                        :: status
    ! Executable code
    call trace_begin ( me, "CpL1BData", &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    rename = QuantityName
    if ( present(NewName) ) rename = NewName
    if ( .not. associated(L1BFile1) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'null pointer passed to readl1bdata--probably cant find ' // &
        & trim(QuantityName) )
    call mls_openFile ( L1BFile2, status )
    call ReadL1BData_MLSFile ( L1BFile1, QuantityName, L1BData, noMAFs, flag, &
      & l2aux=l2aux  )
    if ( associated(L1BData%CharField) ) then
      call SaveAsHDF5DS( L1BFile2%FileID%f_id, rename, &
        & L1BData%CharField(:,:,1) )
    endif
    if ( associated(L1BData%DpField) ) then
      call SaveAsHDF5DS( L1BFile2%FileID%f_id, rename, L1BData%DpField )
    endif
    if ( associated(L1BData%IntField) ) then
      call SaveAsHDF5DS( L1BFile2%FileID%f_id, rename, L1BData%IntField )
    endif
    call deallocateL1BData ( L1BData )
    call mls_CloseFile( L1BFile2 )
    call trace_end ( "CpL1BData", &
      & cond=toggle(gen) .and. levels(gen) > 1 )
  end subroutine CpL1BData

  ! --------------------------------------------------  DupL1BData  -----
  subroutine DupL1BData ( l1bData1, l1bData2, offsetMAF )
    ! cp all components from l1bdata1 to l1bdata2
    ! increment counterMAF by offsetMAF (if present)
    type( L1BData_T ), intent(in)  :: L1bData1
    type( L1BData_T ), intent(out) :: L1bData2
    integer, optional, intent(in)  :: offsetMAF

    ! Executable code
    call allocatel1bdata ( l1bData2, L1bDataSibling=l1bData1 )
    ! print *, 'l1bData1%counterMAF(1): ', l1bData1%counterMAF
    if ( .not. present(offsetMAF) ) then
      l1bData2%counterMAF = l1bData1%counterMAF
    else
      l1bData2%counterMAF = l1bData1%counterMAF + offsetMAF
      l1bData2%FirstMAFCtr = l1bData1%FirstMAFCtr + offsetMAF
      l1bData2%LastMAFCtr = l1bData1%LastMAFCtr + offsetMAF
    end if
    ! print *, 'l1bData2%counterMAF(1): ', l1bData2%counterMAF
    if ( associated(l1bdata1%charField) ) l1bdata2%charField=l1bdata1%charField
    if ( associated(l1bdata1%intField) ) l1bdata2%intField=l1bdata1%intField
    if ( associated(l1bdata1%dpField) ) l1bdata2%dpField=l1bdata1%dpField
  end subroutine DupL1BData

  ! ------------------------------------------  DeallocateL1BData  -----
  subroutine DeallocateL1BData ( l1bData )
    ! This should be called when finished with an l1bData.
    type( L1BData_T ), intent(inout) :: L1bData

    ! Executable code
    call deallocate_test ( l1bData%counterMAF, 'l1bData%counterMAF', ModuleName )
    call deallocate_test ( l1bData%charField, 'l1bData%charField', ModuleName )
    call deallocate_test ( l1bData%intField, 'l1bData%intField', ModuleName )
    call deallocate_test ( l1bData%dpField, 'l1bData%dpField', ModuleName )
  end subroutine DeallocateL1BData

  ! ------------------------------------------------  DiffL1BData  -----
  subroutine DiffL1BData ( l1bData1, l1bData2, &
    & details, options, numDiffs, mafStart, mafEnd, l1bValues1, l1bValues2, &
    & Period )
  use MLSStrings, only: Asciify, IsAllAscii
    ! Diff two l1brad quantities
    type( L1BData_T ), intent(inout) :: L1bData1
    type( L1BData_T ), intent(inout) :: L1bData2
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even counterMAF
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1
    character(len=*), intent(in), optional :: options
    ! logical, intent(in), optional :: silent  ! don't print anything
    integer, intent(out), optional :: numDiffs  ! how many diffs
    integer, intent(in), optional :: mafStart, mafEnd
    real(r8), dimension(:), optional :: l1bValues1
    real(r8), dimension(:), optional :: l1bValues2
    real(r8), optional, intent(in)   :: Period
    ! If options contains 'r' or 's', print much less
    ! if options contains 'd' don't bother with essentially equal and so on
    ! Local variables
    integer :: i,j,k
    logical :: hideAssocStatus
    logical :: l1b1NotFinite
    logical :: l1b2NotFinite
    integer :: MYDETAILS
    logical :: myDirect
    logical :: myPeriodic
    integer :: mafStart1, mafStart2, mafEnd1, mafEnd2
    integer :: myNumDiffs
    real(r8) :: myPeriod
    logical :: mySilent
    logical :: prntAssocStatus  ! Whether to remark on association status
                                !  of multidimensional arrays
    logical, parameter :: DEBUG = .false.
    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details
    hideAssocStatus = .false.
    if ( present(options) ) then
      hideAssocStatus = any(indexes(options, (/'s','r','@'/)) > 0 )
    end if
    prntAssocStatus = .not. hideAssocStatus
    if ( present(mafStart) ) then
      mafStart1 = mafStart
      mafStart2 = mafStart
    else
      mafStart1 = 1
      mafStart2 = 1
    end if
    if ( present(mafEnd) ) then
      mafEnd1 = mafEnd
      mafEnd2 = mafEnd
    else
      mafEnd1 = L1bData1%NoMAFs
      mafEnd2 = L1bData1%NoMAFs
    end if

    mySilent = .false.
    if ( present(options) ) mySilent = ( index( options, 'h' ) > 0 )
    myDirect = .false.
    if ( present(options) ) myDirect = ( index( options, 'd' ) > 0 )
    myPeriod = 360._r8
    if ( present(Period) ) myPeriod = Period
    myPeriodic = .false.
    if ( present(options) ) myPeriodic = ( index( options, 'p' ) > 0 )
    nameOnEachLine = ' '
    if ( DEBUG ) then
      call outputNamedValue( 'options', options )
      call outputNamedValue( 'myDetails', myDetails )
      call outputNamedValue( 'mySilent', mysilent )
    end if
    ! if ( mySilent ) call suspendOutput

    myNumDiffs = 0
    if ( trim(L1bData1%NameInst) /= trim(L1bData2%NameInst) ) then
      call output(trim(L1bData1%NameInst), advance='yes')
      call output(trim(L1bData2%NameInst), advance='yes')
      myNumDiffs = myNumDiffs + 1
    end if
    if ( trim(L1bData1%L1BName) /= trim(L1bData2%L1BName) ) then
      call output('L1B rad quantity (1) Name = ', advance='no')
      call output(trim(L1bData1%L1BName), advance='yes')
      call output('L1B rad quantity (2) Name = ', advance='no')
      call output(trim(L1bData2%L1BName), advance='yes')
      myNumDiffs = myNumDiffs + 1
    end if
    if ( myDetails < -1 ) then
      call doneHere
      return
    end if
    if ( L1bData1%FirstMAF /= L1bData2%FirstMAF ) then
      call output(' (1) First major frame read = ', advance='no')
      call output(L1bData1%FirstMAF, advance='yes')
      call output(' (2) First major frame read = ', advance='no')
      call output(L1bData2%FirstMAF, advance='yes')
      myNumDiffs = myNumDiffs + 1
    end if
    if ( L1bData1%NoMAFs /= L1bData2%NoMAFs ) then
      call output(' (1) Num of MAFs read = ', advance='no')
      call output(L1bData1%NoMAFs, advance='yes')
      call output(' (2) Num of MAFs read = ', advance='no')
      call output(L1bData2%NoMAFs, advance='yes')
      myNumDiffs = myNumDiffs + 1
    end if
    if ( L1bData1%MaxMIFs /= L1bData2%MaxMIFs ) then
      call output(' (1) Max # of MIFs/MAF in SD array = ', advance='no')
      call output(L1bData1%MaxMIFs, advance='yes')
      call output(' (2) Max # of MIFs/MAF in SD array = ', advance='no')
      call output(L1bData2%MaxMIFs, advance='yes')
      myNumDiffs = myNumDiffs + 1
    end if
    if ( L1bData1%NoAuxInds /= L1bData2%NoAuxInds ) then
      call output(' (1) Num of auxilliary indices = ', advance='no')
      call output(L1bData1%NoAuxInds, advance='yes')
      call output(' (2) Num of auxilliary indices = ', advance='no')
      call output(L1bData2%NoAuxInds, advance='yes')
      myNumDiffs = myNumDiffs + 1
    end if
    if ( L1bData1%FirstMAFCtr /= L1bData2%FirstMAFCtr ) then
      call output(' (1) First major frame counter = ', advance='no')
      call output(L1bData1%FirstMAFCtr, advance='yes')
      call output(' (2) First major frame counter = ', advance='no')
      call output(L1bData2%FirstMAFCtr, advance='yes')
      myNumDiffs = myNumDiffs + 1
    end if
    if ( L1bData1%LastMAFCtr /= L1bData2%LastMAFCtr ) then
      call output(' (1) Last major frame counter = ', advance='no')
      call output(L1bData1%LastMAFCtr, advance='yes')
      call output(' (2) Last major frame counter = ', advance='no')
      call output(L1bData2%LastMAFCtr, advance='yes')
      myNumDiffs = myNumDiffs + 1
    end if
    if ( L1bData1%FirstMAFCtr < L1bData2%FirstMAFCtr ) then
      mafStart1 = mafStart1 + L1bData2%FirstMAFCtr - L1bData1%FirstMAFCtr
    else if ( L1bData1%FirstMAFCtr > L1bData2%FirstMAFCtr ) then
      mafStart2 = mafStart2 + L1bData1%FirstMAFCtr - L1bData2%FirstMAFCtr
    end if
    if ( L1bData1%LastMAFCtr < L1bData2%LastMAFCtr ) then
      mafEnd2 = mafEnd2 - ( L1bData2%LastMAFCtr - L1bData1%LastMAFCtr )
    else if ( L1bData1%LastMAFCtr > L1bData2%LastMAFCtr ) then
      mafEnd1 = mafEnd1 - ( L1bData1%LastMAFCtr - L1bData2%LastMAFCtr )
    end if
    if ( myDetails < 0 ) then
      call doneHere
      return
    end if
    if ( associated(l1bData1%counterMAF) .and. &
      & associated(l1bData2%counterMAF) ) then
      if ( any(l1bData1%counterMAF(1:mafEnd2) /= l1bData2%counterMAF(1:mafEnd2))) then
        call dump ( l1bData1%counterMAF(1:mafEnd2) - l1bData2%counterMAF(1:mafEnd2), &
          & 'l1bData%counterMAF (diff)' )
         myNumDiffs = myNumDiffs + 1
       end if
    else
      if ( prntAssocStatus ) &
        & call output('(CounterMAF arrays not associated)', advance='yes')
    end if

    if ( myDetails < 1 ) then
      if ( DEBUG ) call output( ' Done here', advance='yes' )
      call doneHere
      return
    end if
    if ( statsOnOneLine ) nameOnEachLine = L1BData1%L1BName
    if ( associated(l1bData1%charField) .and. &
      & associated(l1bData2%charField)) then
      if ( any(l1bData1%charField /= l1bData2%charField) ) then
        myNumDiffs = myNumDiffs + count(l1bData1%charField /= l1bData2%charField)
        ! where ( .not. isAllAscii(l1bData1%CharField) )
        !  l1bData1%CharField = asciify( l1bData1%CharField, 'snip')
        ! end where
        do k=1, size(l1bData1%CharField, 3)
          do j=1, size(l1bData1%CharField, 2)
            do i=1, size(l1bData1%CharField, 1)
              if ( .not. isAllAscii(l1bData1%CharField(i,j,k) ) ) &
              l1bData1%CharField(i,j,k) = asciify( l1bData1%CharField(i,j,k), 'snip')
            end do
          end do
        end do
        ! where ( .not. isAllAscii(l1bData2%CharField) )
        !  l1bData2%CharField = asciify( l1bData2%CharField, 'snip')
        ! end where
        do k=1, size(l1bData1%CharField, 3)
          do j=1, size(l1bData1%CharField, 2)
            do i=1, size(l1bData1%CharField, 1)
              if ( .not. isAllAscii(l1bData2%CharField(i,j,k) ) ) &
              l1bData2%CharField(i,j,k) = asciify( l1bData2%CharField(i,j,k), 'snip')
            end do
          end do
        end do
        call dump ( l1bData1%CharField, 'l1bData1%CharField' )
        call dump ( l1bData2%CharField, 'l1bData2%CharField' )
      end if
    else
      if ( prntAssocStatus ) &
        & call output('(CharField arrays not associated)', advance='yes')
    end if

    if ( associated(l1bData1%intField) &
      & .and. associated(l1bData1%intField) ) then
      if ( any(l1bData1%intField /= l1bData2%intField) ) then
        call output( 'About to call dump with l1bData%intField', advance='yes' )
        call dump ( l1bData1%intField - l1bData2%intField, &
          & 'l1bData%intField (diff)', options=options )
        myNumDiffs = myNumDiffs + count(l1bData1%intField /= l1bData2%intField)
      else if ( .not. mySilent ) then
        call output('(integer arrays are exactly equal)', advance='yes')
      end if
    else
      if ( prntAssocStatus ) &
        & call output('(intField arrays not associated)', advance='yes')
    end if

    if ( present(l1bValues1) .and. present(l1bValues2) ) then
      l1b1NotFinite = .not. any(ieee_is_finite(l1bValues1))
      l1b2NotFinite = .not. any(ieee_is_finite(l1bValues2))
      if ( l1b1NotFinite .and. l1b2NotFinite ) then
        call output('both dpField arrays all NaNs', advance='yes')
      else if ( l1b1NotFinite ) then
        call output('l1bValues1 array all NaNs', advance='yes')
      else if ( l1b2NotFinite ) then
        call output('l1bValues2 array all NaNs', advance='yes')
      else
        if ( DEBUG ) call output( 'Calling direct diff with l1bvalues1 and 2', advance='yes' )
        myNumDiffs = myNumDiffs + count(l1bValues1 /= l1bValues2)
        if ( .not. myPeriodic ) then
          call DIFF ( &
            & l1bValues1, ' ', &
            & l1bValues2, ' ', &
            & options=options )
        else
          call dump ( &
            & diff_fun( l1bValues1, l1bValues1, &
            & auxvalue=myPeriod , options=options), &
            & options=options )
        end if
      end if
    end if

    if ( associated(l1bData1%dpField) &
      & .and. associated(l1bData1%dpField) ) then
      ! if ( any(l1bData1%dpField /= l1bData2%dpField) ) &
      l1b1NotFinite = .not. any(ieee_is_finite(l1bData1%dpField))
      l1b2NotFinite = .not. any(ieee_is_finite(l1bData2%dpField))
      if ( DEBUG ) then
        call outputNamedValue( 'l1b1NotFinite', l1b1NotFinite )
        call outputNamedValue( 'l1b2NotFinite', l1b2NotFinite )
      end if
      if ( l1b1NotFinite .and. l1b2NotFinite ) then
        call output('both dpField arrays all NaNs', advance='yes')
      else if ( l1b1NotFinite ) then
        call output('l1bData1%dpField array all NaNs', advance='yes')
      else if ( l1b2NotFinite ) then
        call output('l1bData2%dpField array all NaNs', advance='yes')
      else if ( myPeriodic ) then
        myNumDiffs = myNumDiffs + count(l1bData1%dpField /= l1bData2%dpField)
        call dump ( &
          & diff_fun( l1bData1%dpField, l1bData2%dpField, &
          & auxvalue=myPeriod , options=options), &
          & options=options )
      else if ( myDirect ) then
        if ( DEBUG ) call output( 'Calling direct diff', advance='yes' )
        myNumDiffs = myNumDiffs + count(l1bData1%dpField /= l1bData2%dpField)
        ! l1bData1%dpField = l1bData1%dpField - l1bData2%dpField
        ! call dump ( &
        !  & l1bData1%dpField, '(1)-(2)', &
        !  & options=options )
        call DIFF ( &
          & l1bData1%dpField, ' ', &
          & l1bData2%dpField, ' ', &
          & options=options )
      else if ( .not. EssentiallyEqual(l1bData1%dpField, l1bData2%dpField, &
        & FillValue=REAL(undefinedValue, R8)) ) then
        if ( DEBUG ) call output( 'Calling diff', advance='yes' )
        call diff ( &
        & l1bData1%dpField(:,:,mafStart1:mafEnd1), ' ', &
        & l1bData2%dpField(:,:,mafStart2:mafEnd2), ' ', &
        & FillValue=REAL(undefinedValue, R8), &
        & options=options )
        myNumDiffs = myNumDiffs + count(l1bData1%dpField /= l1bData2%dpField)
      end if
    else
      if ( prntAssocStatus ) &
        & call output('(dpField arrays not associated)', advance='yes')
    end if
    call doneHere

  contains
    subroutine doneHere
      ! Housekeeping
      ! call resumeOutput
      if ( present(numDiffs) ) numDiffs = myNumDiffs
      nameOnEachLine = ' '
    end subroutine doneHere
  end subroutine DiffL1BData

  ! ------------------------------------------------  DumpL1BData  -----
  subroutine DumpL1BData ( l1bData, details, options )
    ! Disclose pertinent, perhaps damning facts about an l1brad quantity
    type( L1BData_T ), intent(inout) :: L1bData
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even counterMAF
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1
    character(len=*), intent(in), optional :: options ! Passed dumping arrays

    ! Local variables
    integer :: MYDETAILS

    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details

    if ( trim(L1bData%NameInst) /= ' ' ) &
      & call output(trim(L1bData%NameInst), advance='yes')
    call output('L1B rad/oa quantity Name = ', advance='no')
    call output(trim(L1bData%L1BName), advance='yes')
    if ( myDetails < -1 ) return
    call output('  First major frame read = ', advance='no')
    call output(L1bData%FirstMAF, advance='yes')
    call output('  Num of MAFs read = ', advance='no')
    call output(L1bData%NoMAFs, advance='yes')
    call output('  Max # of MIFs/MAF in SD array = ', advance='no')
    call output(L1bData%MaxMIFs, advance='yes')
    call output('  Num of auxilliary indices = ', advance='no')
    call output(L1bData%NoAuxInds, advance='yes')
    call output('  First major frame counter = ', advance='no')
    call output(L1bData%FirstMAFCtr, advance='yes')
    call output('  Last major frame counter = ', advance='no')
    call output(L1bData%LastMAFCtr, advance='yes')

    if ( myDetails < 0 ) return
    if ( associated(l1bData%counterMAF) ) then
      call dump ( l1bData%counterMAF, 'l1bData%counterMAF', format='(i8)' )
    else
      call output('(CounterMAF array not associated)', advance='yes')
    end if

    if ( myDetails < 1 ) return
    if ( associated(l1bData%charField) ) then
      call dump ( l1bData%CharField, 'l1bData%CharField' )
    else
      call output('(CharField array not associated)', advance='yes')
    end if

    if ( associated(l1bData%intField) ) then
      call dump ( l1bData%intField, 'l1bData%intField', &
        & fillValue = int(undefinedValue), options=options )
    else
      call output('(intField array not associated)', advance='yes')
    end if

    if ( associated(l1bData%dpField) ) then
      call dump ( l1bData%dpField, 'l1bData%dpField', &
        & fillValue=undefinedValue*1.d0, options=options )
    else
      call output('(dpField array not associated)', advance='yes')
    end if
  end subroutine DumpL1BData

  ! -------------------------------------------------  GetL1BFile  -----
  function GetL1BFile ( filedatabase, fieldName, options, object ) &
    & result(item)
  ! Which MLSFile contains a given sd name or attribute name
  !
  ! options, if present, may influence how it works
  ! By convention, options is preceded by a "-"; it is ignored, however
  !
  ! option      effect
  ! ------      ------
  !  a       search for sttribute, not dataset
  !  d       search for dataset, not attribute
  !  /       search for attribute under FileID%grp_id; i.e., '/'
  !  f       search for attribute under FileID%f_id
  !  s       search for attribute under FileID%sd_id
  !  w       Wildcard * which allows 'a*' to match 'abcd'
  !  c       case insensitive which allows 'ABCD' to match 'abcd'

  ! Tries to adapt to where attributes are stored in file
  ! This means, e.g. if the attribute is attached to a ds, opening
  ! that ds.
  ! Afterwards we atempt to cleanup, closing the ds.
    use MLSHDF5, only: IsHDF5ItemPresent
    use MLSFiles, only: Dump
    use HDF5, only: H5DOpen_F, H5GClose_F, H5GOpen_F

    type (MLSFile_T), dimension(:), pointer :: FILEDATABASE
    character (len=*), intent(in)           :: fieldName ! Name of field
    character (len=*), optional, intent(in) :: options ! E.g., -a
    character (len=*), optional, intent(in) :: object  ! If we need to open
    type(MLSFile_T), pointer                :: item
    logical, parameter :: DEEBUG = .false.

    ! Externals
    integer, external :: SFN2INDEX

    ! Local variables
    logical :: alreadyOpen
    integer :: i
    integer :: locID
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: myhdfVersion
    character(len=80) :: myObject
    character(len=8) :: myOptions
    character(len=1) :: openedCode ! What type of object did we open? s,g,f
    integer :: returnStatus

    ! Executable code
    ! DEEBUG = (fieldname=='DACsDeconvolved')
    call trace_begin ( me, 'GetL1BFile' , cond=.false. )
    nullify(item)
    myOptions = '-d' ! By default, search for sd names
    if (present(options)) myOptions = options
    myObject = '/'
    if (present(object)) myObject = object
    openedCode = ' '
    do i = 1, size(filedatabase)
      returnStatus = 0
      if ( streq(filedatabase(i)%content, 'l1b*', '-w') ) then
        if ( DEEBUG ) call dump(filedatabase(i), details=1)
        ! This is an l1b file
        alreadyOpen = filedatabase(i)%StillOpen
        if ( .not. alreadyOpen ) &
          & call mls_openFile(filedatabase(i), returnStatus)
        ! print *, 'Oops--must open file?'
        if ( returnStatus /= 0 ) &
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'unable to open L1BFile searching for ' // trim(fieldname), &
            & MLSFile=filedatabase(i) )
        myHDFVersion = filedatabase(i)%HDFVersion
        ! Are we looking for an attribute?
        ! (Then we may need an object name)
        if ( index(myoptions, 'a') > 0 ) then
          if ( myhdfVersion == HDFVERSION_4 ) then
            openedCode = 'n' ! Because we can't read attributes from HDF4 files
          ! What object would the attribute be attached to?
          else if ( index(myOptions, 'f') > 0 ) then
            openedCode = 'f'
          else if ( index(myOptions, '/') > 0 ) then
            call h5gopen_f ( filedatabase(i)%fileID%f_id, trim(myObject), &
              & filedatabase(i)%fileID%grp_id, returnStatus )
            openedCode = 'g'
          else if ( index(myOptions, 's') > 0 ) then
            call h5dopen_f ( filedatabase(i)%fileID%f_id, trim(myObject), &
              & filedatabase(i)%fileID%sd_id, returnStatus )
            openedCode = 's'
          else
            call h5dopen_f ( filedatabase(i)%fileID%f_id, trim(myObject), &
              & filedatabase(i)%fileID%sd_id, returnStatus )
            openedCode = 's'
          end if
        end if

        if ( myhdfVersion == HDFVERSION_4 ) then
          if ( index(myOptions, 'a') > 0 ) then
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Cannot search for attribute names in hdf4 files' )
            cycle
          end if
          if ( sfn2index(filedatabase(i)%FileID%f_id, trim(fieldName)) /= -1 ) then
            item => filedatabase(i)
            if ( .not. alreadyOpen ) &
              & call mls_closeFile(filedatabase(i), returnStatus)
            call trace_end ( 'GetL1BFile' , cond=.false. )
            return
          end if
        else
          ! Are we looking for an attribute?
          if ( index(myoptions, 'a') > 0 ) then
            ! What object would the attribute be attached to?
            if ( index(myOptions, 'f') > 0 ) then
              locID = filedatabase(i)%FileID%f_id
            else if ( index(myOptions, '/') > 0 ) then
              locID = filedatabase(i)%FileID%grp_id
            else if ( index(myOptions, 's') > 0 ) then
              locID = filedatabase(i)%FileID%sd_id
            else
              locID = filedatabase(i)%FileID%f_id
            end if
          else
            locID = filedatabase(i)%FileID%f_id
          end if
          if ( DEEBUG ) call outputnamedValue('locID', locID)
          if ( IsHDF5ItemPresent( &
            & locID, trim(fieldName), myOptions ) &
            & ) then
            item => filedatabase(i)
            if ( .not. alreadyOpen ) call closeEverything(filedatabase(i))
            call trace_end ( 'GetL1BFile' , cond=.false. )
            return
          end if
        end if
        if ( .not. alreadyOpen ) call closeEverything(filedatabase(i))
      end if
    end do
    call trace_end ( 'GetL1BFile' , cond=.false. )
    ! if ( DEEBUG ) stop

  contains

    subroutine closeEverything(L1BFile)
      ! Args:
      type(MLSFile_T) :: L1BFile
      integer :: returnStatus
      ! Executable
       select case (openedCode)
       case ('s')
         call h5gclose_f ( L1BFile%fileID%sd_id, returnStatus )
       case ('g')
         call h5gclose_f ( L1BFile%fileID%grp_id, returnStatus )
       case default
       end select
       call mls_closeFile(L1BFile, returnStatus)
    end subroutine closeEverything
  end function GetL1BFile

  ! ------------------------------------------------  FindL1BData  -----
  integer function FindL1BData ( filedatabase, fieldName, hdfVersion )

  use MLSHDF5, only: IsHDF5DSPresent

    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    ! integer, dimension(:), intent(in) :: files ! File handles
    character (len=*), intent(in) :: fieldName ! Name of field
    integer, optional, intent(in) :: hdfVersion
    logical, parameter :: TRUSTDATABASE = .true.

    ! Externals
    integer, external :: SFN2INDEX

    ! Local variables
    integer :: i
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: myhdfVersion

    ! Executable code
    call trace_begin ( me, 'FindL1BData' , cond=.false. )
    if ( present(hdfVersion) ) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L1BDEFAULT_HDFVERSION
    end if

    findL1BData=0
    ! print *, 'Looking for ', trim(fieldname)
    do i = 1, size(filedatabase)
      if ( TRUSTDATABASE ) myHDFVersion = filedatabase(i)%HDFVersion
      if ( .not. streq(filedatabase(i)%content, 'l1b*', '-w') ) then
        ! This is not an l1b file
      else if ( myhdfVersion == HDFVERSION_4 ) then
        if ( sfn2index(filedatabase(i)%FileID%f_id, trim(fieldName)) /= -1 ) then
          findL1BData = filedatabase(i)%FileID%f_id
          return
        end if
      else
        ! print *, 'Looking in ', trim(filedatabase(i)%name)
        if ( IsHDF5DSPresent(filedatabase(i)%FileID%f_id, trim(fieldName)) ) then
          findL1BData = filedatabase(i)%FileID%f_id
          ! print *, 'Eureka '
          return
        end if
      end if
    end do
    call trace_end ( 'FindL1BData' , cond=.false. )

  end function FindL1BData

  ! -------------------------------------------------  FindMaxMAF_arr  -----
  integer function FindMaxMAF_arr ( files, minMAF )
  ! Find maximum MAF among files (using counterMAF arrays)

    type(MLSFile_T), dimension(:), intent(in), target    :: files ! File handles
    integer, optional, intent(out)                       :: minMAF

    ! Local variables
    integer :: i
    type(MLSFile_T), pointer :: L1BFile
    integer :: myMinMAF
    integer :: theMinMAF
    logical :: haveCtrMAF
    logical, parameter :: DEEBug = .false.
    ! Executable code

    FindMaxMAF_arr = 0
    myMinMAF = BIGGESTMAFCTR
    haveCtrMAF = .false.
    do i = 1, size(files)
      L1BFile => files(i)
      if ( streq( L1BFile%content, 'l1b*', '-w') ) then
        ! This is an l1b file
        FindMaxMAF_arr = max( FindMaxMAF_sca( L1BFile, theMinMAF ), FindMaxMAF_arr )
        myMinMAF = min( myMinMAF, theMinMAF )
      end if
    end do
    if ( DEEBug ) print *, 'FindMaxMAF_arr ', FindMaxMAF_arr
    if ( DEEBug ) print *, 'myMinMAF ', myMinMAF
    if ( present(minMAF) ) minMAF = myMinMAF
  end function FindMaxMAF_arr

  ! -------------------------------------------------  FindMaxMAF_sca  -----
  integer function FindMaxMAF_sca ( L1BFile, minMAF )
  ! Find maximum MAF among files (using counterMAF arrays)
  ! If there is no counterMAF, just return the isze of the MAFStartTimeTAI

  use MLSHDF5, only: IsHDF5DSPresent

    type(MLSFile_T), pointer               :: L1BFile ! File handles
    integer, optional, intent(out)         :: minMAF

    ! Externals
    integer, external :: SFN2INDEX

    ! Local variables
    logical :: alreadyOpen
    character(len=*), parameter :: fieldname = 'counterMAF'
    type(L1BData_T) :: l1bData
    integer :: myhdfVersion
    integer :: myMinMAF
    integer :: noMAFs
    integer :: status
    logical :: haveCounterMAF
    logical :: haveMAFStartTimeTAI
    logical, parameter :: DEEBug = .false.
    character(len=1024) :: DSNames
    ! Executable code
    alreadyOpen = L1BFile%StillOpen
    myMinMAF = 0
    if ( .not. alreadyOpen ) &
      & call mls_openFile( L1BFile, status )
    myHDFVersion = L1BFile%HDFVersion
    if ( myhdfVersion == HDFVERSION_4 ) then
      haveCounterMAF = sfn2index(L1BFile%FileID%f_id,trim(fieldName)) /= -1
      if ( haveCounterMAF ) then
        call ReadL1BData ( L1BFile, fieldName, L1bData, noMAFs, status, &
          & dontPad=.true. )
      end if
    else
      haveCounterMAF = IsHDF5DSPresent( L1BFile, trim(fieldName) )
      haveMAFStartTimeTAI = IsHDF5DSPresent( L1BFile, 'MAFStartTimeTAI' )
      if ( DEEBug ) call outputNamedValue( 'haveCounterMAF', haveCounterMAF )
      if ( DEEBug )  call outputNamedValue( 'haveMAFStartTimeTAI', haveMAFStartTimeTAI )
      if ( haveCounterMAF ) then
        call ReadL1BData ( L1BFile, fieldName, L1bData, noMAFs, status, &
          & dontPad=.true.)
      elseif ( haveMAFStartTimeTAI ) then
        call ReadL1BData ( L1BFile, 'MAFStartTimeTAI', L1bData, noMAFs, status, &
          & dontPad=.true.)
      end if
    end if
    if ( haveCounterMAF ) then
      FindMaxMAF_sca =  maxval(l1bData%counterMAF)
      myMinMAF =  myminval(l1bData%counterMAF)
      if ( DEEBug ) print *, 'counterMAF ', l1bData%counterMAF
      ! call output('Shape L1b counterMAF ')
      ! call output(shape(l1bData%intField), advance='yes')
      call deallocatel1bdata(L1bData)
    elseif ( haveMAFStartTimeTAI ) then
      FindMaxMAF_sca =  size(l1bData%dpField)
      myMinMAF = 1
      if ( DEEBug ) print *, 'MAFStartTimeTAI ', l1bData%dpField
      ! call output('Shape L1b counterMAF ')
      ! call output(shape(l1bData%intField), advance='yes')
      call deallocatel1bdata(L1bData)
    else
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Failed to find '//trim(fieldName)// &
      & ' in l1b files while FindMaxMAF_sca', MLSFile=L1BFile )
      call GetAllHDF5DSNames( L1BFile, DSNames )
      call Dump ( DSNames, 'DSNames' )
    end if
    if ( .not. alreadyOpen ) &
      & call mls_closeFile( L1BFile, status )
    if ( present(minMAF) ) minMAF = myMinMAF
  end function FindMaxMAF_sca

  ! -------------------------------------------------  IsL1BGappy  -----
  logical function IsL1BGappy ( l1bData, ignoreGlobalAttrs )
    ! Look for gaps in l1bData, returning true if any found
  use PCFHdr, only: GlobalAttributes

    ! Dummy arguments
    type (L1BData_T), intent(in) :: l1bData
    logical, optional, intent(in) :: ignoreGlobalAttrs  ! defaults to checking
    ! Internal variables
    logical :: skipGACHecks
    integer :: i
    integer :: nextMAF
    integer :: maxCtrMAF
    integer :: minCtrMAF
    ! Executable
    skipGACHecks = .false.
    if ( present(ignoreGlobalAttrs) ) skipGACHecks = ignoreGlobalAttrs
    IsL1BGappy = .false.
    if ( .not. skipGACHecks ) then
      maxCtrMAF = maxval(l1BData%counterMAF)
      minCtrMAF = myminval(l1BData%counterMAF)
      IsL1BGappy = (maxCtrMAF < GlobalAttributes%LastMAFCtr &
        & .or. &
        & minCtrMAF > GlobalAttributes%FirstMAFCtr)
    end if
    if ( IsL1BGappy ) return
    ! Search for interior gaps
    IsL1BGappy = .true. ! If exit from next loop prematurely, must be gappy
    nextMAF = minCtrMAF
    do i=1, size(l1BData%counterMAF)
      if ( nextMAF /= l1BData%counterMAF(i)) return
    end do
    IsL1BGappy = .false.
  end function IsL1BGappy

  ! -------------------------------------------------  showL1BGaps  -----
  subroutine showL1BGaps ( l1bData, ignoreGlobalAttrs )
    ! Look for gaps in l1bData, returning true if any found
  use PCFHdr, only: GlobalAttributes

    ! Dummy arguments
    type (L1BData_T), intent(in) :: l1bData
    logical, optional, intent(in) :: ignoreGlobalAttrs  ! defaults to checking
    ! Internal variables
    logical :: skipGACHecks
    integer :: i
    integer :: nextMAF
    integer :: maxCtrMAF
    integer :: minCtrMAF
    ! Executable
    skipGACHecks = .false.
    if ( present(ignoreGlobalAttrs) ) skipGACHecks = ignoreGlobalAttrs
    if ( .not. skipGACHecks ) then
      maxCtrMAF = maxval(l1BData%counterMAF)
      minCtrMAF = myminval(l1BData%counterMAF)
      if ( maxCtrMAF < GlobalAttributes%LastMAFCtr ) &
        & call outputnamedValue( 'maxCtrMAF, LastMAFCtr', &
        &                      (/ maxCtrMAF, GlobalAttributes%LastMAFCtr /) )
      if ( minCtrMAF > GlobalAttributes%FirstMAFCtr ) &
        & call outputnamedValue( 'minCtrMAF, FirstMAFCtr', &
        &                      (/ minCtrMAF, GlobalAttributes%FirstMAFCtr /) )
    end if
    ! Search for interior gaps
    nextMAF = minCtrMAF
    do i=1, size(l1BData%counterMAF)
      if ( nextMAF /= l1BData%counterMAF(i)) exit
      nextMAF = nextMAF + 1
    end do
    if ( i > size(l1BData%counterMAF) ) then
      call output( 'No interior gaps found in counterMAF array', advance='yes' )
    else
      call outputNamedValue( 'First counterMAF gap found at MAF ', i )
    endif
  end subroutine showL1BGaps

  !--------------------------------------------------  L1BOAHDF4DSName  -----
  function L1BOAHDF4DSName ( HDF5Name ) result( HDF4Name )
    ! Convert an l1boa ds name from hdf5 to hdf4

    ! Dummy arguments
    character(len=*), intent(in)      :: HDF5Name
    character(len=128)                :: HDF4Name
    character(len=128)                :: Name
    ! Executable
    Name = HDF5Name
    call ConvertL1BOADSNames ( HDF4Name, Name, '5to4' )
  end function L1BOAHDF4DSName

  !--------------------------------------------------  L1BOAHDF5DSName  -----
  function L1BOAHDF5DSName ( HDF4Name ) result( HDF5Name )
    ! Convert an l1boa ds name from hdf4 to hdf5

    ! Dummy arguments
    character(len=*), intent(in)      :: HDF4Name
    character(len=128)                :: HDF5Name
    character(len=128)                :: Name
    ! Executable
    Name = HDF4Name
    call ConvertL1BOADSNames ( Name, HDF5Name, '4to5' )
  end function L1BOAHDF5DSName

  !--------------------------------------------------  L1BOASetup  -----
  subroutine L1boaSetup ( root, filedatabase, F_FILE, hdfVersion )
    ! Take file name from l2cf, open, and add new file to filedatabase

    ! Dummy arguments
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    integer, intent(in) :: ROOT         ! of the l1brad file specification.
    integer, intent(in) :: F_FILE       ! From init_tables_module
    integer, optional, intent(inout) :: hdfVersion

    ! Local variables

    character(len=FileNameLen) :: FileName

    integer :: I                        ! Loop inductor, subscript
    integer :: numFiles
    integer :: returnStatus
    integer :: SON                      ! Some subtree of root.
    type(MLSFile_T) :: L1BFile

    ! Executable code
    error = 0
    if ( present(hdfVersion) ) hdfVersion = FILENOTFOUND
    ! Collect data from the fields. (only one legal field: file='...')
    do i = 2, nsons(root)
      son = subtree(i,root)
      if ( get_field_id(son) == f_file ) then
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        returnStatus = InitializeMLSFile(L1BFile, content = 'l1boa', &
          & type=l_hdf, access=DFACC_RDONLY, name=fileName)
        call mls_openFile(L1BFile, returnStatus)
        if ( returnStatus == 0 ) then
          numFiles = addFileToDatabase(filedatabase, L1BFile)
          if ( present(hdfVersion) ) hdfVersion = L1BFile%HDFVersion
        end if
      else
        call announce_error ( son, &
          & 'Unknown field specified in read l1boa' )
      end if
    end do
  end subroutine L1boaSetup

  ! ------------------------------------------------- L1BRadSetup  -----
  subroutine L1bradSetup ( Root, filedatabase, F_File, &
    & hdfVersion )
    ! Take file name from l2cf, open, and add new file to filedatabase
    ! Dummy arguments
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    integer, intent(in) :: ROOT         ! of the l1brad file specification.
    integer, intent(in) :: F_FILE
    integer, optional, intent(inout) :: hdfVersion

    ! Local variables
    character(len=FileNameLen) :: FILENAME

    integer :: I                        ! Loop inductor, subscript
    integer :: numFiles
    integer :: returnStatus
    integer :: SON                      ! Some subtree of root.

    type(MLSFile_T) :: L1BFile

    ! Executable code
    error = 0
    if ( present(hdfVersion) ) hdfVersion = FILENOTFOUND

    ! Collect data from the fields. (only one legal field: file='...')
    do i = 2, nsons(root)
      son = subtree(i,root)
      if ( get_field_id(son) == f_file ) then
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        returnStatus = InitializeMLSFile(L1BFile, content = 'l1brad', &
          & type=l_hdf, access=DFACC_RDONLY, name=fileName)
        call mls_openFile(L1BFile, returnStatus)
        if ( returnStatus == 0 ) then
          numFiles = addFileToDatabase(filedatabase, L1BFile)
          if ( present(hdfVersion) ) hdfVersion = L1BFile%HDFVersion
        end if
      else
        call announce_error ( son, &
          & 'Unknown field specified in read l1brad' )
      end if
    end do

  end subroutine L1bradSetup

  ! ----------------------------------  ReadL1BAttribute_intarr1l  -----
  subroutine ReadL1BAttribute_intarr1 ( L1BFile, value, AttrName, Flag )

    use MLSHDF5, only: IsHDF5AttributePresent, GetHDF5Attribute
    use HDF5, only: H5GClosE_F, H5GOpen_F

    ! Dummy arguments
    type(MLSFile_T), intent(in)            :: L1BFile
    integer, intent(out) :: value(:) ! Result
    character(len=*), intent(in) :: AttrName ! attribute name to retrieve
    integer, intent(out) :: Flag        ! Error flag

    ! Local variables
    integer :: myhdfVersion
    integer :: aID, status
    logical :: alreadyOpen
    ! Executable code
    alreadyOpen = L1BFile%StillOpen
    if ( .not. alreadyOpen ) &
      & call mls_openFile(L1BFile, status)
        myHDFVersion = L1BFile%HDFVersion

    Flag = 0

    if ( myhdfVersion == HDFVERSION_4 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Not implemented in hdf4 l1boa file', MLSFile=L1bFile)
        Flag = -1
    else
      call h5gOpen_f (L1BFile%FileID%f_id,'/', aID, status)
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to open group attribute in l1boa file', MLSFile=L1BFile )
        Flag = -1
      end if
      if ( .not. IsHDF5AttributePresent(aID, AttrName) ) then
        Flag = -1
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Failed to find attribute in l1boa file'//AttrName, MLSFile=L1bFile)
      else
        if ( DEEBUG ) then
          call output ('get attribute ', advance='no')
          call output (AttrName, advance='yes')
        end if
        call GetHDF5Attribute(aID, AttrName, value)
      end if
      call h5gClose_f (aID, status)
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to close group attribute in l1boa file', MLSFile=L1bFile )
      end if
    end if
    if ( .not. alreadyOpen ) &
      & call mls_closeFile(L1BFile, status)
  end subroutine ReadL1BAttribute_intarr1

  ! -----------------------------------  ReadL1BAttribute_dblarr1  -----
  subroutine ReadL1BAttribute_dblarr1 ( L1BFile, value, AttrName, Flag )

    use MLSHDF5, only: IsHDF5AttributePresent, GetHDF5Attribute
    use HDF5, only: H5GClose_F, H5GOpen_F

    ! Dummy arguments
    type(MLSFile_T), intent(in)            :: L1BFile
    real(r8), intent(out) :: value(:) ! Result
    character(len=*), intent(in) :: AttrName ! attribute name to retrieve
    integer, intent(out) :: Flag        ! Error flag

    ! Local variables
    integer :: myhdfVersion
    integer :: aID, status
    logical :: alreadyOpen
    ! Executable code
    alreadyOpen = L1BFile%StillOpen
    if ( .not. alreadyOpen ) &
      & call mls_openFile(L1BFile, status)
        myHDFVersion = L1BFile%HDFVersion

    Flag = 0

    if ( myhdfVersion == HDFVERSION_4 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Not implemented in hdf4 l1boa file', MLSFile=L1bFile)
      Flag = -1
    else
      call h5gOpen_f (L1BFile%FileID%f_id,'/', aID, status)
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to open group attribute in l1boa file', MLSFile=L1BFile )
        Flag = -1
      end if
      if ( .not. IsHDF5AttributePresent(aID, AttrName) ) then
        Flag = -1
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Failed to find attribute in l1boa file'//AttrName, MLSFile=L1bFile)
      else
        if ( DEEBUG ) then
          call output ('get attribute ', advance='no')
          call output (AttrName, advance='yes')
        endif
        call GetHDF5Attribute(aID, AttrName, value)
      end if
      call h5gClose_f (aID, status)
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to close group attribute in l1boa file', MLSFile=L1bFile )
      end if
    end if
    if ( .not. alreadyOpen ) &
      & call mls_closeFile(L1BFile, status)
  end subroutine ReadL1BAttribute_dblarr1

  ! -------------------------------------------------  PadL1BData  -----
  subroutine PadL1BData ( L1BDataIn, L1BDataOut, FirstMAF, LastMAF, NoMAFs )
    ! Pad l1bdataIn to fit NoMafs (assuming noMafs >= l1bdatain%NoMAFs)
    ! beginning with 1st MAF number of day
    ! (May need to pad at start, at end, or at both ends)
    ! Return newly created padded object as l1bdataOut
    ! If optionally supplied PrecisionIn it will return padded PrecisionOut
    ! where the padded values will all be cleverly set negative
    ! (meaning do not use)
    ! Just as important, fill any gaps
    ! where a gap is defined as a break in sequence in counterMAF array
    ! .., p-2, p-1, [p], [p+1], .., p+gap-1, ..
    ! such that the bracketed numbers above are absent from the counterMAF
    ! Dummy arguments
    type(L1BData_T), intent(in)  :: L1BDataIn
    type(L1BData_T), intent(out) :: L1BDataOut
    integer, intent(in)          :: FirstMAF ! 1st MAF
    integer, intent(in)          :: LastMAF  ! last MAF number of day
    integer, intent(out)         :: NoMAFs   ! number
    ! Internal variables
    integer, dimension(3) :: dims
    integer               :: maf
    integer               :: oldSize
    logical, parameter    :: DEEBug = .false.
    ! Executable
    if ( associated(l1bDataIn%charField)) then
      dims = shape(l1bDataIn%charField)
    else if ( associated(l1bDataIn%intField)) then
      dims = shape(l1bDataIn%intField)
    else if ( associated(l1bDataIn%dpField)) then
      dims = shape(l1bDataIn%dpField)
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Padl1bData was passed l1bDataIn w/o allocating it' )
    end if
    oldSize = dims(3)
    NoMAFs = max(dims(3), lastMAF-firstMAF+1)
    dims(3) = NoMAFs

    if ( DEEBug ) print *, 'Preparing to pad ', trim(l1bDataIn%L1BName)
    if ( DEEBug ) print *, 'NoMAFs ', NoMAFs
    if ( DEEBug ) print *, 'l1b(in)NoMAFs ', l1bDataIn%NoMAFs
    if ( DEEBug ) print *, 'dims ', dims
    if ( DEEBug ) print *, 'min(counterMAF) ', myminval(l1bDataIn%counterMAF)
    if ( DEEBug ) print *, 'min(dpField) ',minval(l1bDataIn%dpField(1,1,:))
    if ( DEEBug ) print *, 'max(dpField) ',maxval(l1bDataIn%dpField(1,1,:))
    call allocateL1BData ( l1bDataOut, dims, L1bDataSibling=l1bDataIn )
    do maf=1, NoMAFs
      l1bDataOut%counterMAF(maf) = FirstMAF - 1 + maf
    end do
    call zeroField(l1bdataOut, 1, undefinedValue, m=NoMAFs)
    call cpField(l1bdataIn, 1, l1bdataOut, 1, m=oldSize)
  end subroutine PadL1BData

  ! ----------------------------------------------  BadPadL1BData  -----
  subroutine BadPadL1BData ( L1BDataIn, L1BDataOut, FirstMAFCtr, NoMAFs, &
    & PrecisionIn, PrecisionOut, force )
    ! Pad l1bdataIn to fit NoMafs (assuming noMafs >= l1bdatain%NoMAFs)
    ! beginning with 1st MAF number of day
    ! (May need to pad at start, at end, or at both ends)
    ! Return newly created padded object as l1bdataOut
    ! If optionally supplied PrecisionIn it will return padded PrecisionOut
    ! where the padded values will all be cleverly set negative
    ! (meaning do not use)
    ! Just as important, fill any gaps
    ! where a gap is defined as a break in sequence in counterMAF array
    ! .., p-2, p-1, [p], [p+1], .., p+gap-1, ..
    ! such that the bracketed numbers above are absent from the counterMAF
    ! Dummy arguments
    type(L1BData_T), intent(in)  :: L1BDataIn
    type(L1BData_T), intent(out) :: L1BDataOut
    integer, intent(in)          :: FirstMAFCtr ! 1st MAF number of day
    integer, intent(in)          :: NoMAFs
    type(L1BData_T), optional, intent(in)  :: PrecisionIn
    type(L1BData_T), optional, intent(out) :: PrecisionOut
    logical, optional, intent(in) :: force   ! Relax assumption
    ! Internal variables
    logical :: MyForce
    integer :: Dims(3)
    integer :: Gap
    integer :: I
    integer :: IndexOut
    integer :: Cmindex
    integer :: Current
    integer :: Maf
    integer :: Rank
    logical, parameter :: DEEBug = .false.
    ! Executable
    myForce = .false.
    if ( present(force) ) myForce = force
    ! Check assumption that noMafs >= l1bdatain%NoMAFs
    if ( .not. myForce .and. noMAFs < l1bdatain%NoMAFs ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'noMAFs requested smaller than number stored in input l1bData' )
    else if ( present(PrecisionOut) .and. .not. present(PrecisionIn) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Must have PrecisionIn to calculate PrecisionOut' )
    end if
    if ( associated(l1bDataIn%charField)) then
      dims = shape(l1bDataIn%charField)
    else if ( associated(l1bDataIn%intField)) then
      dims = shape(l1bDataIn%intField)
    else if ( associated(l1bDataIn%dpField)) then
      dims = shape(l1bDataIn%dpField)
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'BadPadL1BData was passed l1bDataIn w/o allocating it' )
    end if
    dims(3) = noMAFs
    rank = l1bDataIn%trueRank
    if ( DEEBug ) print *, 'Preparing to pad ', trim(l1bDataIn%L1BName)
    if ( DEEBug ) print *, 'rank ', rank
    if ( DEEBug ) print *, 'NoMAFs ', NoMAFs
    if ( DEEBug ) print *, 'l1b(in)NoMAFs ', l1bDataIn%NoMAFs
    if ( DEEBug ) print *, 'dims ', dims
    if ( DEEBug ) print *, 'min(counterMAF) ', myminval(l1bDataIn%counterMAF)
    if ( DEEBug ) print *, 'min(dpField) ',minval(l1bDataIn%dpField(1,1,:))
    if ( DEEBug ) print *, 'max(dpField) ',maxval(l1bDataIn%dpField(1,1,:))
    call allocateL1BData ( l1bDataOut, dims, L1bDataSibling=l1bDataIn )
    do maf=1, NoMAFs
      l1bDataOut%counterMAF(maf) = FirstMAFCtr - 1 + maf
      if ( present(PrecisionOut) ) &
        & PrecisionOut%counterMAF(maf) = FirstMAFCtr - 1 + maf
    end do
    ! Now hunt for missing MAFs
    maf = FirstMAFCtr
    indexOut = 0
    if ( DEEBug ) print *, 'size(l1bDataOut%counterMAF) ', size(l1bDataOut%counterMAF)
    do cmindex = 1, l1bDataIn%noMAFs
      if ( DEEBug ) print *, 'l1bDataIn%counterMAF(cmindex) ', l1bDataIn%counterMAF(cmindex)
      if ( DEEBug ) print *, 'maf ', maf
      if ( DEEBug ) print *, 'indexOut ', indexOut
      if ( l1bDataIn%counterMAF(cmindex) < maf ) then
        ! counterMAF too small: must ignore (perhaps past end of array?)
      else if ( l1bDataIn%counterMAF(cmindex) == maf ) then
        ! Found expected maf: no gap from last maf
        indexOut = indexOut + 1
        if ( DEEBug ) print *, 'Found expected maf: no gap from last maf'
        l1bDataOut%counterMAF(indexOut) = maf
        call cpField(l1bDataIn, cmindex, l1bdataOut, indexOut)
        if ( present(PrecisionOut) ) &
          &  call cpField(PrecisionIn, cmindex, PrecisionOut, indexOut)
        maf = maf+1
      else
        current = l1bDataIn%counterMAF(cmindex)
        if ( DEEBug ) print *, 'Found a greater maf than expected: gap from maf to ', current-1
        ! Found a greater maf than expected: gap from maf to current-1
        gap = current-maf
        ! fill the gap with undefined
        do i=1, gap
          l1bDataOut%counterMAF(indexOut+i) = maf-1+i
          call zeroField(l1bdataOut, indexOut+i, undefinedValue)
          if ( present(PrecisionOut) ) &
            & call zeroField(PrecisionOut, indexOut+i, undefinedValue)
        end do
        if ( DEEBug ) print *, 'Now treat current normally ', current
        ! Now treat current normally
        indexOut = indexOut + gap + 1
        l1bDataOut%counterMAF(indexOut) = current
        call cpField(l1bDataIn, cmindex, l1bdataOut, indexOut)
        if ( present(PrecisionOut) ) &
          &  call cpField(PrecisionIn, cmindex, PrecisionOut, indexOut)
        maf = current+1
      end if
    end do
    if ( maf > FirstMAFCtr + noMAFs - 1 ) return
    gap = FirstMAFCtr + noMAFs - maf
    ! Apparently, we end too soon, so must pad
    if ( DEEBug ) print *, 'Apparently, we end too soon, so must pad ', gap, indexOut
    do i=1, gap
      l1bDataOut%counterMAF(indexOut+i) = maf-1+i
      call zeroField(l1bdataOut, indexOut+i, undefinedValue)
      if ( present(PrecisionOut) ) &
        & call zeroField(PrecisionOut, indexOut+i, undefinedValue)
    end do
  end subroutine BadPadL1BData

  ! -------------------------------------  ReadL1BData_fileHandle  -----
  ! In time we will do away with most file-handle based interfaces
  subroutine ReadL1BData_fileHandle ( L1FileHandle, QuantityName, L1bData, NoMAFs, Flag, &
    & FirstMAF, LastMAF, NeverFail, hdfVersion, dontPad, L2AUX )
    ! Dummy arguments
    character(len=*), intent(in)   :: QUANTITYNAME ! Name of SD to read
    integer, intent(in)            :: L1FILEHANDLE ! From HDF
    integer, intent(in), optional  :: FIRSTMAF ! First to read (default 0)
    integer, intent(in), optional  :: LASTMAF ! Last to read (default last/file)
    logical, intent(in), optional  :: NEVERFAIL ! Don't quit if TRUE
    logical, intent(in), optional  :: L2AUX     ! Don't even warn if TRUE
    type(l1bdata_t), intent(inout) :: L1BDATA ! Result
    integer, intent(out) :: FLAG        ! Error flag
    integer, intent(out) :: NOMAFS      ! Number actually read
    integer, optional, intent(in) :: HDFVERSION
    logical, intent(in), optional  :: DONTPAD ! Don't try to pad even if gappy

    ! Local variables
    type(MLSFile_T), target :: f
    type(MLSFile_T), pointer :: fp
    integer :: Me = -1                  ! String index for trace
    integer :: myhdfVersion
    ! logical, parameter :: DEEBug = .false.
    integer :: status

    ! Executable code
    call trace_begin ( me, "ReadL1BData_fileHandle", &
      & cond=toggle(gen) .and. levels(gen) > 1 )

    if ( present(hdfVersion) ) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L1BDEFAULT_HDFVERSION
    end if

    ! Set up an MLSFile instance
    status = InitializeMLSFile( f, content='l1b', &
      & shortname='l1b', type=l_hdf, access=DFACC_RDONLY )
    f%fileID%f_id = L1FileHandle
    f%stillOpen = .true.
    f%hdfVersion = myhdfVersion
    fp => f
    call ReadL1BData_MLSFile ( fp, QuantityName, L1bData, NoMAFs, Flag, &
      & FirstMAF, LastMAF, NeverFail, dontPad, L2AUX )
    call trace_end ( "ReadL1BData_fileHandle", &
      & cond=toggle(gen) .and. levels(gen) > 1 )
  end subroutine ReadL1BData_fileHandle

  ! ----------------------------------------  ReadL1BData_MLSFile  -----
  subroutine ReadL1BData_MLSFile ( L1BFile, QuantityName, L1bData, NoMAFs, Flag, &
    & FirstMAF, LastMAF, NeverFail, dontPad, L2AUX )
    use MLSFiles, only: Dump
    use MLSFillValues, only: IsFillValue
    use Optional_m, only: Default
    use PCFHdr, only: GlobalAttributes
    ! Dummy arguments
    character(len=*), intent(in)   :: QUANTITYNAME ! Name of SD to read
    type(MLSFile_T), pointer       :: L1BFile
    integer, intent(out)           :: FLAG        ! Error flag
    integer, intent(out)           :: NOMAFS      ! Number actually read
    integer, intent(in), optional  :: FIRSTMAF ! First to read (default 0)
    integer, intent(in), optional  :: LASTMAF ! Last to read (default last/file)
    logical, intent(in), optional  :: NEVERFAIL ! Don't quit if TRUE
    logical, intent(in), optional  :: L2AUX     ! Don't even warn if TRUE
    type(l1bdata_t), intent(inout) :: L1BDATA ! Result
    logical, intent(in), optional  :: DONTPAD ! Don't try to pad even if gappy

    ! Local variables
    logical :: alreadyOpen
    logical :: DEEBug
    character(len=MaxNDSNames*128) :: DSNames
    logical :: isScalar
    type(l1bdata_t) :: L1BDATATMP
    integer :: Me = -1                  ! String index for trace
    integer :: myhdfVersion
    logical :: myDontPad
    integer :: returnStatus
    logical, parameter :: ShowName = .false.
    ! Executable code
    DEEBug = .false. ! ( index(QuantityName, '/GHz/GeodAngle') > 0 )
    if ( .not. associated(L1BFile) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'null pointer passed to readl1bdata--probably cant find ' // &
        & trim(QuantityName) )
      return
    endif
    if ( showName ) call output ( '&&&&&&&&&&&&&&&&& ' // trim(Quantityname), advance='yes' )
    call trace_begin ( me, "ReadL1BData_MLSFile", &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    alreadyOpen = L1BFile%StillOpen
    if ( .not. alreadyOpen .and. DEEBug ) then
      print *, 'Oops--need to open l1b file before reading'
      call Dump( L1BFile, details=2 )
      call GetAllHDF5DSNames( L1BFile, DSNames )
      call Dump ( DSNames, 'DSNames' )
      call mls_openFile( L1BFile, returnStatus )
      call Dump( L1BFile, details=2 )
      call GetAllHDF5DSNames( L1BFile, DSNames )
      if ( returnStatus /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l1b file', MLSFile=L1BFile)
    end if
    myhdfVersion = L1BFile%HDFVersion
    ! print * 'hdfVersion: ', hdfVersion
    myDontPad = .false.
    if ( present(firstMAF) ) myDontPad = .true.
    if ( present(dontPad) ) myDontPad = dontPad
    if ( switchDetail(switches, 'l1bread') > -1 ) &
      & call output( 'About to read L1B dataset %% ' // &
      & trim(QuantityName) // &
      & ' %%', advance='yes' )

    if ( myhdfVersion == HDFVERSION_4 ) then
      call ReadL1BData_MF_hdf4 ( L1BFile, trim(QuantityName), L1bData, &
      & NoMAFs, Flag, FirstMAF, LastMAF, NeverFail, L2AUX )
      if ( flag /= 0 ) noMAFs = -1
    else
      call ReadL1BData_MF_hdf5 ( L1BFile, trim(QuantityName), L1bData, &
      & NoMAFs, Flag, FirstMAF, LastMAF, NeverFail, L2AUX )
      !Unfortunately, hdf5-formatted l1b data have different shapes from hdf4
      ! E.g., for MAFStartTimeTAI we obtain the following
      !  hdfVERSION      shape
      !     4           1   1   5
      !     5           5   1   1
      if ( Flag == 0 ) then
        call Reshape_for_hdf4(L1bData)
      else
        call MLSMessage ( MLSMSG_Warning, ModuleName // '/ReadL1BData_MLSFile', &
          & 'Failed to find '//trim(QuantityName)//' in l1b files')
        call OutputNamedValue ( 'read status flag', flag )
        call OutputNamedValue ( 'data type', trim(L1bData%data_type) )
        call Dump ( L1BFile, details=2 )
        call GetAllHDF5DSNames( L1BFile, DSNames )
        call Dump ( DSNames, 'DSNames' )
        if ( .not. Default( neverFail, .false. ) ) call crash_burn
        NOMAFS = -1
        go to 9
      end if
    end if
    isScalar = ( l1bData%noMAFs < 2 ) ! There may be a better way to decide
    if ( index(QuantityName, 'Angle') > 0 .and. DEEBUG ) then
      call outputnamedValue( 'isScalar', isScalar )
      call outputnamedValue( 'myDontPad', myDontPad )
      call outputnamedValue( 'IsL1BGappy(l1bData)', IsL1BGappy(l1bData) )
      if ( IsL1BGappy(l1bData) ) call showL1BGaps( l1bData )
    endif
    if ( myDontPad .or. flag /= 0 .or. isScalar .or. &
      & (GlobalAttributes%FirstMAFCtr > GlobalAttributes%LastMAFCtr) ) go to 9
    if ( .not. IsL1BGappy(l1bData) ) go to 9
    ! Must pad l1bdata; so first copy to temp l1bData
    if ( associated(l1bdata%dpField) .and. DEEBug ) then
      ! print *, 'max(l1bdata) ', maxval(l1bdata%dpField(1,1,:))
      ! print *, 'min(l1bdata) ', minval(l1bdata%dpField(1,1,:))
    end if
    call dupL1BData(l1bdata, l1bdataTmp)
    if ( associated(l1bdatatmp%dpField) .and. DEEBug ) then
      ! print *, 'max(l1bdatatmp) ', maxval(l1bdatatmp%dpField(1,1,:))
      ! print *, 'min(l1bdatatmp) ', minval(l1bdatatmp%dpField(1,1,:))
    end if
    call deallocatel1bData(l1bdata)
    if ( GlobalAttributes%LastMAFCtr > 0 ) &
      & NoMAFs = GlobalAttributes%LastMAFCtr - GlobalAttributes%FirstMAFCtr + 1
    if ( DEEBug ) print *, 'preparing to pad'
    if ( DEEBug ) print *, 'GlobalAttributes%FirstMAFCtr ', GlobalAttributes%FirstMAFCtr
    if ( DEEBug ) print *, 'GlobalAttributes%LastMAFCtr ', GlobalAttributes%LastMAFCtr
    if ( DEEBug ) print *, 'NoMAFs ', NoMAFs
    if ( present(firstMAF) .and. present(lastMAF)) then
      call PadL1BData(l1bdataTmp, l1bData, FirstMAF, LastMAF, NoMAFs)
      call deallocatel1bData(l1bdataTmp)
      if ( .not. alreadyOpen )  call mls_closeFile(L1BFile, returnStatus)
      go to 9
    end if
    call BadPadL1BData(l1bdataTmp, l1bData, GlobalAttributes%FirstMAFCtr, NoMAFs)
    if ( index(QuantityName, 'Angle') > 0 .and. DEEBUG ) then
      call output( 'Gone through padding', advance='yes' )
      call outputnamedValue( 'Fills in new l1bdata?', any(isFillValue(l1bdata%dpField)) )
      call outputnamedValue( 'Fills in l1bdataTmp?', any(isFillValue(l1bdataTmp%dpField)) )
    endif
    call deallocatel1bData(l1bdataTmp)
    if ( .not. present(firstMAF) .and. .not. present(lastMAF) ) go to 9
    ! Need to contract padded l1bdata to just those MAFs requested
    call dupL1BData(l1bdata, l1bdataTmp)
    call deallocatel1bData(l1bdata)
    if ( DEEBug ) print *, 'preparing to contract'
    call ContractL1BData(l1bdataTmp, l1bData, noMAFs, firstMAF, lastMAF)
    call deallocatel1bData(l1bdataTmp)
    if ( .not. alreadyOpen )  call mls_closeFile(L1BFile, returnStatus)
9   continue
    call trace_end ( "ReadL1BData_MLSFile", &
      & cond=toggle(gen) .and. levels(gen) > 1 )
  end subroutine ReadL1BData_MLSFile

  ! ----------------------------------------  ReadL1BData_MF_hdf4  -----
  subroutine ReadL1BData_MF_hdf4 ( L1BFile, QuantityName, L1bData, &
    & NoMAFs, Flag, FirstMAF, LastMAF, NeverFail, L2AUX )

    ! Dummy arguments
    character(len=*), intent(in)   :: QUANTITYNAME ! Name of SD to read
    type(MLSFile_T), pointer       :: L1BFile
    integer, intent(in), optional  :: FIRSTMAF ! First to read (default 0)
    integer, intent(in), optional  :: LASTMAF ! Last to read (default last in file)
    logical, intent(in), optional  :: NEVERFAIL ! Don't quit if TRUE
    logical, intent(in), optional  :: L2AUX     ! Don't even warn if TRUE
    type(l1bdata_t), intent(inout) :: L1BDATA ! Result
    integer, intent(out) :: FLAG        ! Error flag
    integer, intent(out) :: NOMAFS      ! Number actually read

    ! Local Parameters
    character (len=*), parameter :: INPUT_ERR = 'Error in input argument '
    integer, parameter :: MAX_VAR_DIMS = 32
    integer, parameter :: SD_NO_COUNTERMAF = -2

    ! Local Variables

    character (len=128) :: DUMMY        ! Dummy quantity name

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: ALLOC_ERR
    integer :: DATA_TYPE
    integer :: DIM_SIZES(MAX_VAR_DIMS)
    integer :: I
    integer :: L1FileHandle
    integer :: Me = -1                  ! String index for trace
    integer :: MLSMsgLevel              ! Eigher MLSMSG_Error or MLSMSG_Warning
    logical :: MyNeverFail
    integer :: N_ATTRS
    integer :: NUMMAFS
    integer :: RANK
    integer :: SDS1_ID
    integer :: SDS2_ID
    integer :: SDS_INDEX
    integer :: STATUS

    integer, dimension(:), allocatable :: EDGE
    integer, dimension(:), allocatable :: START
    integer, dimension(:), allocatable :: STRIDE

    real(r4), pointer, dimension(:,:,:) :: tmpR4Field

    ! Executable code
    call trace_begin ( me, "ReadL1BData_MF_hdf4", &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    L1FileHandle = L1BFile%FileID%f_id
    call deallocateL1BData ( l1bData ) ! Avoid memory leaks

    flag = 0
    MyNeverFail = .false.
    if ( present(NeverFail) ) MyNeverFail = NeverFail
    MLSMsgLevel = merge(MLSMSG_Warning, MLSMSG_Error, myNeverFail)

    ! Find data sets for counterMAF & quantity by name

    sds_index = sfn2index(L1FileHandle, 'counterMAF')
    if ( sds_index == -1 ) then
      if ( .not. JUSTLIKEL2AUX ) then
        flag = NOCOUNTERMAFINDX
        call MLSMessage ( MLSMsgLevel, ModuleName, &
        & 'Failed to find index of counterMAF data set.', MLSFile=L1BFile)
        go to 9
      end if
      sds1_id = SD_NO_COUNTERMAF
    else

      sds1_id = sfselect(L1FileHandle, sds_index)
      if ( sds1_id == -1 ) then
        flag = NOCOUNTERMAFID
        call MLSMessage ( MLSMsgLevel, ModuleName, &
        & 'Failed to find identifier of counterMAF data set.', MLSFile=L1BFile)
        go to 9
      end if
    end if

    sds_index = sfn2index(L1FileHandle, quantityName)
    if ( sds_index == -1 ) then
      flag = NOQUANTITYINDEX
      dummy = 'Failed to find index of quantity "' // trim(quantityName) // &
        & '" data set.'
      call MLSMessage ( MLSMsgLevel, ModuleName, dummy, MLSFile=L1BFile )
      go to 9
    end if

    sds2_id = sfselect(L1FileHandle, sds_index)
    if ( sds2_id == -1 ) then
      flag = NODATASETID
      call MLSMessage ( MLSMsgLevel, ModuleName, &
      & 'Failed to find identifier of data set matching the index.', MLSFile=L1BFile)
      go to 9
    end if

    ! Find rank (# of dimensions), dimension sizes of quantity data set
    status = sfginfo ( sds2_id, dummy, rank, dim_sizes, data_type, &
      n_attrs )

    if ( status == -1 ) then
      flag = NODATASETRANK
      call MLSMessage ( MLSMsgLevel, ModuleName,&
      & 'Failed to find rank of data set.', MLSFile=L1BFile)
      go to 9
    end if

    ! allocate, based on above SD, dim info; don't track allocatable sizes
    allocate ( edge(rank), start(rank), stride(rank), stat=status)
    call test_allocate ( status, moduleName, 'edge, start, stride' )

    ! Set "slab" dimensions
    edge = dim_sizes(1:rank)
    start = 0
    stride = 1

    ! Fill in "indexing" values of l1b object

    l1bData%L1BName = quantityName
    l1bData%TrueRank = rank

    l1bData%noAuxInds = product(dim_sizes(1:rank-2))
    l1bData%maxMIFs = 1
    if ( rank > 1 ) l1bData%maxMIFs = dim_sizes(rank-1)

    ! Check input arguments, set noMAFs

    numMAFs = dim_sizes(rank)

    if ( present ( firstMAF ) ) then
      if ( (firstMAF >= numMAFs) .or. (firstMAF < 0) ) then
        flag = FIRSTMAFNOTFOUND
        call MLSMEssage ( MLSMsgLevel, ModuleName, &
        & input_err // 'firstMAF (bad chunkNo?)', MLSFile=L1BFile )
        go to 9
      end if
      l1bData%firstMAF = firstMAF
    else
      l1bData%firstMAF = 0
    end if

    if ( present (lastMAF) ) then
      if ( lastMAF < l1bData%firstMAF ) then
        flag = LASTMAFNOTFOUND
        call MLSMEssage ( MLSMsgLevel, ModuleName, &
        & input_err // 'last' , MLSFile=L1BFile)
        go to 9
      end if
      if ( lastMAF >= numMAFs ) then
        l1bData%noMAFs = numMAFs - l1bData%firstMAF
      else
        l1bData%noMAFs = lastMAF - l1bData%firstMAF + 1
      end if
    else
      l1bData%noMAFs = numMAFs - l1bData%firstMAF
    end if

    noMAFs = l1bData%noMAFs
    edge(rank) = l1bData%noMAFs
    start(rank) = l1bData%firstMAF

    ! allocate, read counterMAF
    call Allocate_test ( l1bData%counterMaf, l1bData%noMAFs, &
      & 'l1bData%counterMAF', ModuleName )
    if ( sds1_id /= SD_NO_COUNTERMAF ) then
      status = sfrdata_f90(sds1_id,  (/ l1bData%firstMAF /) , (/1/), &
        & (/l1bData%noMAFs/), l1bData%counterMAF )
      if ( status == -1 ) then
        flag = CANTREADCOUNTERMAF
        call MLSMessage ( MLSMsgLevel, ModuleName, &
        & MLSMSG_L1BRead // 'counterMAF.', MLSFile=L1BFile )
        go to 9
      end if
    else
      ! Since we aren't reading these, just make them internally consistent
      do i = 1, l1bData%noMAFs
        l1bData%counterMAF(i) = l1bData%firstMAF + i - 1
      end do
    end if

    l1bData%FirstMAFCtr = myminval(l1bData%counterMAF)
    l1bData%LastMAFCtr = maxval(l1bData%counterMAF)
    ! allocate, read according to field type; nullify unused pointers
    select case ( data_type )

    case ( DFNT_CHAR8 ) ! ----------------------- character
      allocate ( l1bData%charField(l1bData%noAuxInds,l1bData%maxMIFs, &
        & l1bData%noMAFs), STAT=alloc_err )
      if ( alloc_err /= 0 ) then
        flag = CANTALLOCATECHARS
        ! Announce the error
        call test_allocate ( alloc_err, ModuleName, "l1bData%charField" )
        go to 9
      end if
      ! Account for the allocation size
      addr = 0
      if ( size(l1bData%charField) > 0 ) &
        & addr = transfer(c_loc(l1bData%charField(1,1,1)), addr)
      call test_allocate ( alloc_err, ModuleName, "l1bData%charField", &
        & uBounds = [l1bData%noAuxInds,l1bData%maxMIFs], &
        & elementSize = storage_size(l1bData%charField) / 8, address=addr )
      nullify (l1bData%intField, l1bData%dpField)
      l1bdata%data_type = 'character'
      status = sfrcdata ( sds2_id, start, stride, edge, &
        & l1bData%charField )

    case ( DFNT_INT32 ) ! ------------------------------ integer
      call allocate_test ( l1bData%intField, &
        & l1bData%noAuxInds, l1bData%maxMIFs, l1bData%noMAFs, &
        & 'l1bData%intField', ModuleName )
      nullify ( l1bData%charField, l1bData%dpField )
      l1bdata%data_type = 'integer'
      status = sfrdata_f90 ( sds2_id, start, stride, edge, &
          l1bData%intField )

    case ( DFNT_FLOAT64 ) ! ------------------------------ real (r8)
      call allocate_test ( l1bData%dpField, &
        & l1bData%noAuxInds, l1bData%maxMIFs, l1bData%noMAFs, &
        & 'l1bData%dpField', ModuleName )
      nullify ( l1bData%charField, l1bData%intField )
      l1bdata%data_type = 'double'
      status = sfrdata_f90 ( sds2_id, start, stride, edge, &
          l1bData%dpField )

    case ( DFNT_FLOAT32 ) ! ------------------------------ real (r4)
      call allocate_test ( l1bData%dpField, &
        & l1bData%noAuxInds, l1bData%maxMIFs, l1bData%noMAFs, &
        & 'l1bData%dpField', ModuleName )
      nullify ( l1bData%charField, l1bData%intField, tmpr4Field )
      l1bdata%data_type = 'double'

      call allocate_test ( tmpr4Field, &
        & l1bData%noAuxInds, l1bData%maxMIFs, l1bData%noMAFs, &
        & 'tmpR4Field', ModuleName )
      status = sfrdata_f90 ( sds2_id, start, stride, edge, &
          tmpR4Field )
      l1bData%dpField = tmpR4Field

      call deallocate_test ( tmpr4Field, 'tmpr4Field', ModuleName )

    case default
      status = -2
    end select ! ----------------------------------------------------

    if ( status == -1 ) then
      flag = CANTREAD3DFIELD
      call MLSMessage ( MLSMsgLevel, ModuleName, &
      & MLSMSG_L1BRead // quantityName, MLSFile=L1BFile )
      go to 9
    else if ( status == -2 ) then
      flag = UNKNOWNDATATYPE
      call MLSMessage ( MLSMsgLevel, ModuleName, &
      & 'Unknown data type in readl1bData', MLSFile=L1BFile   )
      go to 9
    end if

    ! Terminate access to the data sets

    if ( sds1_id /= SD_NO_COUNTERMAF ) then
      status = sfendacc(sds1_id)
      if ( status == -1 ) then
        flag = CANTENDCOUNTERMAF
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Failed to terminate access to data sets.', MLSFile=L1BFile )
        flag = -1
        if ( MyNeverFail ) go to 9
      end if
    end if

    status = sfendacc(sds2_id)
    if ( status == -1 ) then
      flag = CANTENDQUANTITY
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to terminate access to data sets.', MLSFile=L1BFile )
      flag = -1
    end if

9   continue

    call trace_end ( "ReadL1BData_MF_hdf4", &
      & cond=toggle(gen) .and. levels(gen) > 1 )

  end subroutine ReadL1BData_MF_hdf4

  ! ----------------------------------------  ReadL1BData_MF_hdf5  -----
  subroutine ReadL1BData_MF_hdf5 ( L1BFile, QuantityName, L1bData, NoMAFs, &
    & Flag, FirstMAF, LastMAF, NeverFail, L2AUX )
    use MLSFillValues, only: IsFillValue
    use HDF5, only: HSize_T
    use MLSHDF5, only: IsHDF5DSPresent, LoadFromHDF5DS, &
      & GetHDF5DSRank, GetHDF5DSDims, GetHDF5DSQType

    ! Dummy arguments
    character(len=*), intent(in)   :: QUANTITYNAME ! Name of SD to read
    type(MLSFile_T), pointer       :: L1BFile
    integer, intent(in), optional  :: FIRSTMAF ! First to read (default 0)
    integer, intent(in), optional  :: LASTMAF ! Last to read (default last/file)
    logical, intent(in), optional  :: NEVERFAIL ! Don't quit if TRUE
    logical, intent(in), optional  :: L2AUX     ! Don't even warn if TRUE
    type(l1bdata_t), intent(inout) :: L1BDATA ! Result
    integer, intent(out) :: FLAG        ! Error flag
    integer, intent(out) :: NOMAFS      ! Number actually read

    ! Local Parameters
    character (len=*), parameter :: INPUT_ERR = 'Error in input argument '

    ! Local Variables
    character (len=128) :: DUMMY        ! Dummy quantity name
    character(len=1) :: Char_rank
    real(r8), dimension(:,:,:,:), pointer :: DP4Buf => null()
    character(len=1024) :: DSNames
    integer :: I
    integer :: L1FileHandle
    integer :: MAFoffset
    integer :: Me = -1                  ! String index for trace
    integer :: MLSMsgLevel              ! Eigher MLSMSG_Error or MLSMSG_Warning
    logical :: MyNeverFail
    integer :: NUMMAFS
    integer :: RANK
    integer :: CMRANK
    character(len=16) :: QTYPE
    integer :: STATUS

    integer(kind=hSize_t), dimension(:), allocatable :: HDIMS
    integer,               dimension(:), allocatable :: DIMS
    integer(kind=hSize_t), dimension(:), allocatable :: MAXDIMS
    integer(kind=hSize_t), dimension(:), allocatable :: CMDIMS
    integer(kind=hSize_t), dimension(:), allocatable :: CMMAXDIMS

    logical, parameter           :: DEEBUG = .false.
    logical :: isL2AUX

    ! Executable code
    call trace_begin ( me, "ReadL1BData_MF_hdf5", &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    L1FileHandle = L1BFile%FileID%f_id
    if ( L1FileHandle < 1 .or. .not. L1BFile%stillOpen ) then
      call MLS_OpenFile ( L1BFile, status )
      if ( status /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open L1B File', MLSFile=L1BFile )
      L1FileHandle = L1BFile%FileID%f_id
    endif
    if ( present(FirstMAF) ) then
      MAFoffset = max(0, FirstMAF)   ! Never let this be < 0
    end if
    call deallocateL1BData ( l1bData ) ! Avoid memory leaks

    flag = 0
    MyNeverFail = .false.
    if ( present(NeverFail) ) MyNeverFail = NeverFail
    MLSMsgLevel = merge(MLSMSG_Warning, MLSMSG_Error, myNeverFail)
    isL2AUX = .false.
    if ( present(l2AUX) ) isL2AUX = L2AUX

    if ( .not. IsHDF5DSPresent(L1FileHandle, QuantityName) ) then
      flag = NOQUANTITYINDEX
      print *, 'Oops--' // trim(QuantityName) // ' not here'
      dummy = 'Failed to find index of quantity "' // trim(quantityName) // &
        & '" data set.'
      call MLSMessage ( MLSMsgLevel, ModuleName, dummy, MLSFile=L1BFile )
      go to 9
    end if

    ! Find Qtype, rank and dimensions of QuantityName
    ! print*, ' Find Qtype, rank and dimensions of QuantityName ', trim(QuantityName)
    call GetHDF5DSRank(L1FileHandle, QuantityName, rank)
    l1bData%TrueRank = rank
    allocate ( dims(rank), Hdims(rank), maxDims(rank), stat=status )
    if ( status /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Could not allocate dims in ReadL1BData_MF_hdf5' )
    call GetHDF5DSDims(L1FileHandle, QuantityName, Hdims, maxDims)
    dims = Hdims
    call GetHDF5DSQType ( L1FileHandle, QuantityName, Qtype )
    if ( DEEBug ) print *, 'L1FileHandle ', L1FileHandle
    if ( DEEBug ) print *, 'QuantityName ', trim(QuantityName)
    if ( DEEBug ) print *, 'rank ', rank
    if ( DEEBug ) print *, 'maxDims ', maxDims
    if ( DEEBug ) print *, 'dims ', dims
    if ( DEEBug ) print *, 'Qtype ', Qtype
    if ( present(firstMAF) .and. DEEBug ) print *, 'firstMAF ', firstMAF
    if ( present(LastMAF) .and. DEEBug ) print *, 'LastMAF ', LastMAF

    l1bData%noAuxInds = product(dims(1:rank-2))
    l1bData%maxMIFs = 1
    if ( rank > 1 ) l1bData%maxMIFs = dims(rank-1)

    ! Check input arguments, set noMAFs

    numMAFs = dims(rank)

    if ( present ( firstMAF ) ) then
      if ( (firstMAF >= numMAFs) .or. (firstMAF < 0) ) then
        flag = FIRSTMAFNOTFOUND
        if ( MyNeverFail ) go to 9
        call MLSMEssage ( MLSMsgLevel, ModuleName, &
        & input_err // 'firstMAF (bad chunkNo?)', MLSFile=L1BFile )
        go to 9
      end if
      l1bData%firstMAF = firstMAF
    else
      l1bData%firstMAF = 0
    end if

    if ( present (lastMAF) ) then
      if ( lastMAF < l1bData%firstMAF ) then
        flag = LASTMAFNOTFOUND
        call dump( (/lastMAF, l1bData%firstMAF/), 'lastMAF, l1bData%firstMAF')
        call MLSMEssage ( MLSMsgLevel, ModuleName, &
        & input_err // 'last', MLSFile=L1BFile )
        go to 9
      end if
      if ( lastMAF >= numMAFs ) then
        l1bData%noMAFs = numMAFs - l1bData%firstMAF
      else
        l1bData%noMAFs = lastMAF - l1bData%firstMAF + 1
      end if
    else
      l1bData%noMAFs = numMAFs - l1bData%firstMAF
    end if

    if ( DEEBug ) print *, 'l1bData%noMAFs ', l1bData%noMAFs
    l1bData%L1BName = quantityName

    noMAFs = l1bData%noMAFs
    if ( DEEBug ) print *, 'noMAFs ', noMAFs

    call Allocate_test ( l1bData%counterMaf, l1bData%noMAFs, &
      & 'l1bData%counterMAF', ModuleName )
    ! Find data sets for counterMAF & quantity by name

    if ( .not. IsHDF5DSPresent(L1FileHandle, '/counterMAF') ) then
      if ( DEEBUG ) print *, 'no counterMAF array in file'
      if ( .not. JUSTLIKEL2AUX ) then
        flag = NOCOUNTERMAFINDX
        call MLSMessage ( MLSMsgLevel, ModuleName, &
        & 'Failed to find index of counterMAF data set.', MLSFile=L1BFile)
        go to 9
      else
        if ( .not. isL2AUX ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to find index of counterMAF data set.')
        ! Since we aren't reading these, just make them internally consistent
        do i = 1, l1bData%noMAFs
          l1bData%counterMAF(i) = l1bData%firstMAF + i - 1
        end do
        ! Have we damaged L1BFile?
        if ( trim(Quantityname) == '/MAFStartTimeTAI' .and. .false. ) then
          call GetAllHDF5DSNames( L1BFile, DSNames )
          call Dump ( DSNames, 'DSNames' )
        endif
      end if
    else
      if ( DEEBUG) print *, 'Getting counterMAF rank'
      call GetHDF5DSRank(L1FileHandle, '/counterMAF', cmrank)
      if ( DEEBUG) print *, cmrank
      allocate ( cmdims(cmrank), cmmaxDims(cmrank) )
      if ( DEEBUG) print *, 'getting counterMAF dims'
      call GetHDF5DSDims(L1FileHandle, '/counterMAF', cmdims, cmmaxDims)
      if ( DEEBUG) print *, cmdims, cmmaxDims
      if ( DEEBUG) print *, 'l1bData%noMAFs', l1bData%noMAFs
      ! allocate(countermaf_ptr(MAX_NOMAFS))
      ! countermaf_ptr = 0
      if ( present(FirstMAF) ) then
        if ( DEEBUG) print *, 'reading counterMAF ', MAFoffset, l1bData%noMAFs
        call LoadFromHDF5DS(L1FileHandle, '/counterMAF', l1bData%counterMaf, &
          & (/MAFoffset/), (/l1bData%noMAFs/) )
      else if ( cmdims(1) /= l1bData%noMAFs ) then
        flag = CANTREADCOUNTERMAF
        if ( MyNeverFail ) then
          ! deallocate(dims, maxDims, cmdims, cmmaxdims)
          ! go to 9
          ! Since we aren't reading these, just make them internally consistent
          do i = 1, l1bData%noMAFs
            l1bData%counterMAF(i) = l1bData%firstMAF + i - 1
          end do
          flag = 0
        end if
        if ( DEEBug ) then
          call output('Quantity name: ', advance = 'no')
          call output(Quantityname, advance = 'yes')
          call output('Quantity rank: ', advance = 'no')
          call output(rank, advance = 'yes')
          call output('Quantity dims: ', advance = 'no')
          call output(int(dims), advance = 'yes')
          call output('Quantity noMAFs: ', advance = 'no')
          call output(l1bData%noMAFs, advance = 'yes')
          call output('counterMAF dims: ', advance = 'no')
          call output(int(cmdims), advance = 'yes')
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Sorry--ReadL1BData_hdf5 says counterMaf sized differently from ' &
            & // trim(QuantityName), MLSFile=L1BFile )
        endif
      else
        if ( DEEBUG) print *, 'reading counterMAF '
        ! This intermediate cm_array shouldn't be necessary
        ! unfortunately, w/o it sometimes hdf5 bombs
        ! allocate(cm_array(cmdims(1)), stat=status)
        ! allocate(cm_array(MAX_NOMAFS), stat=status)
        ! call LoadFromHDF5DS(L1FileHandle, '/counterMAF', cm_array)
        ! l1bData%counterMaf = cm_array(1:cmdims(1))
        ! deallocate(cm_array)
        call LoadFromHDF5DS(L1FileHandle, '/counterMAF', l1bData%counterMaf)
      end if
      if ( DEEBUG) print *, 'deallocating counterMAF dims'
      deallocate(cmdims, cmmaxdims)
    end if

    if ( DEEBUG) print *, 'computing FirstMAFCtr'
    l1bData%FirstMAFCtr = myminval(l1bData%counterMAF)
    if ( DEEBUG) print *, 'computing LastMAFCtr'
    l1bData%LastMAFCtr = maxval(l1bData%counterMAF)
    if ( DEEBUG) print *, 'l1bData%FirstMAFCtr ', l1bData%FirstMAFCtr
    if ( DEEBUG) print *, 'l1bData%LastMAFCtr ', l1bData%LastMAFCtr
    if ( DEEBUG) print *, 'About to use LoadFromHDF5DS to read ', trim(QuantityName)
    if ( DEEBUG) print *, 'dims ', dims
    if ( DEEBUG) print *, 'noMAFs ', l1bData%noMAFs
    write(Char_rank, '(i1)') rank
    if ( DEEBUG) print *, trim(Qtype) // Char_rank
    if ( DEEBUG .and. present(FirstMAF) ) then
      print *, 'FirstMAF-1 ', FirstMAF-1
      print *, 'MAFoffset ', MAFoffset
      print *, 'l1bData%noMAFs ', l1bData%noMAFs
    end if

    ! The following is a crude and despicable hack
    ! Please repair GetHDF5DSQType and delete these lines
    if ( Qtype == 'unknown' ) Qtype = 'character'

    select case (trim(Qtype) // Char_rank)
    case ('real1')
      call allocate_test ( l1bData%DpField, l1bData%noMAFs, 1, 1, "l1bData%DpField", &
        & moduleName )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1), &
          & (/MAFoffset/), (/l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1))
      end if
      l1bdata%data_type = 'double'
    case ('real2')
      call allocate_test ( l1bData%DpField, dims(1), l1bData%noMAFs, 1, "l1bData%DpField", &
        & moduleName )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1), &
          & (/0,MAFoffset/), (/dims(1),l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1))
      end if
      l1bdata%data_type = 'double'
      if ( index(QuantityName, 'Angle') < 1 .or. .true. ) then
      elseif ( any(isFillValue(l1bData%dpField) ) ) then
        call output( 'Fill values among 2d real l1bdata', advance='yes' )
        call outputNamedValue( 'shape(l1bData%dpField)', shape(l1bData%dpField) )
        call dump( l1bData%dpField, 'l1bData%dpField' )
      else
        call output( '2d real l1bdata ' // trim(QuantityName) // 'are clean', advance='yes' )
      endif
    case ('real3')
      call allocate_test ( l1bData%DpField, dims(1), dims(2), l1bData%noMAFs, &
        & "l1bData%DpField", moduleName )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField, &
          & (/0,0,MAFoffset/), (/dims(1),dims(2),l1bData%noMAFs/) )
      else
         call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField)
      end if
      if ( index(QuantityName, 'Angle') < 1 .or. .true. ) then
      elseif ( any(isFillValue(l1bData%dpField) ) ) then
        call output( 'Fill values among 3d real l1bdata', advance='yes' )
        call outputNamedValue( 'shape(l1bData%dpField)', shape(l1bData%dpField) )
        call dump( l1bData%dpField, 'l1bData%dpField' )
      else
        call output( '3d real l1bdata are clean', advance='yes' )
      endif
      l1bdata%data_type = 'double'
    case ('real4')
      call allocate_test ( l1bData%DpField, dims(1)*dims(2), dims(3), l1bData%noMAFs, &
        & "l1bData%DpField", moduleName )
      call allocate_test( DP4Buf, dims(1), dims(2), dims(3), l1bData%noMAFs, &
        & 'DP4Buf', ModuleName )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS( L1FileHandle, QuantityName, DP4Buf, &
          & (/0,0,0,MAFoffset/), &
          & (/dims(1),dims(2),dims(3),l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS( L1FileHandle, QuantityName, DP4Buf)
      end if
      l1bData%DpField = reshape( DP4Buf, &
        & (/dims(1)*dims(2), dims(3), l1bData%noMAFs /) )
      call deallocate_test( DP4Buf, 'DP4Buf', ModuleName )
      l1bdata%data_type = 'double'
    case ('double1')
      call allocate_test ( l1bData%DpField, l1bData%noMAFs, 1, 1, "l1bData%DpField", &
        & moduleName )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1), &
          & (/MAFoffset/), (/l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1))
      end if
      l1bdata%data_type = 'double'
    case ('double2')
      call allocate_test ( l1bData%DpField, dims(1), l1bData%noMAFs, 1, "l1bData%DpField", &
        & moduleName )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1), &
          & (/0,MAFoffset/), (/dims(1),l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1))
      end if
      l1bdata%data_type = 'double'
    case ('double3')
      call allocate_test ( l1bData%DpField, dims(1), dims(2), l1bData%noMAFs, &
        & "l1bData%DpField", moduleName )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField, &
          & (/0,0,MAFoffset/), (/ dims(1),dims(2),l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField)
      end if
      l1bdata%data_type = 'double'
    case ('double4')
      call allocate_test ( l1bData%DpField, dims(1)*dims(2), dims(3), l1bData%noMAFs, &
        & "l1bData%DpField", moduleName )
      call allocate_test( DP4Buf, dims(1), dims(2), dims(3), l1bData%noMAFs, &
        & 'DP4Buf', ModuleName )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS( L1FileHandle, QuantityName, DP4Buf, &
          & (/0,0,0,MAFoffset/), &
          & (/dims(1),dims(2),dims(3),l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS( L1FileHandle, QuantityName, DP4Buf)
      end if
      l1bData%DpField = reshape( DP4Buf, &
        & (/ dims(1)*dims(2), dims(3), l1bData%noMAFs /) )
      call deallocate_test( DP4Buf, 'DP4Buf', ModuleName )
      l1bdata%data_type = 'double'
    case ('integer1')
      allocate( l1bData%intField(l1bData%noMAFs, 1, 1),stat=status)
      call test_allocate ( status, moduleName, "l1bData%intField", &
        & uBounds = [l1bData%noMAFs, 1, 1], &
        & elementSize = storage_size(l1bData%intField) / 8 )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%intField(:,1,1), &
          & (/MAFoffset/), (/l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%intField(:,1,1))
      end if
      l1bdata%data_type = 'integer'
    case ('integer2')
      allocate( l1bData%IntField(dims(1),l1bData%noMAFs, 1),stat=status)
      call test_allocate ( status, moduleName, "l1bData%intField", &
        & uBounds = [int(dims(1)),l1bData%noMAFs, 1], &
        & elementSize = storage_size(l1bData%intField) / 8 )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%IntField(:,:,1), &
          & (/0,MAFoffset/), (/int(dims(1)),l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%IntField(:,:,1))
      end if
      l1bdata%data_type = 'integer'
    case ('integer3')
      allocate( l1bData%IntField(dims(1),dims(2),l1bData%noMAFs),stat=status)
      call test_allocate ( status, moduleName, "l1bData%intField", &
        & uBounds = [int(dims(1:2)),l1bData%noMAFs], &
        & elementSize = storage_size(l1bData%intField) / 8 )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%IntField, &
          & (/0,0,MAFoffset/), (/int(dims(1)),int(dims(2)),l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%IntField)
      end if
!     call MLSMessage ( MLSMSG_Error, ModuleName, &
!     & 'Sorry--LoadFromHDF5DS not yet written for type integer(:,:,:).', MLSFile=L1BFile)
      l1bdata%data_type = 'integer'
    case ('character1')
      allocate( l1bData%charField(l1bData%noMAFs, 1, 1),stat=status)
      call test_allocate ( status, moduleName, "l1bData%charField", &
        & uBounds = [l1bData%noMAFs, 1, 1], &
        & elementSize = storage_size(l1bData%charField) / 8 )
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%charField(:,1,1), &
          & (/MAFoffset/), (/l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%charField(:,1,1))
      end if
      ! call MLSMessage ( MLSMSG_Error, ModuleName, &
      ! & 'Sorry--LoadFromHDF5DS not yet written for type character(:).')
      l1bdata%data_type = 'character'
    case ('character3')
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--LoadFromHDF5DS not yet rewritten for type char(:,:,:).', MLSFile=L1BFile)
      l1bdata%data_type = 'integer'
    case default
      l1bdata%data_type = trim(Qtype) // Char_rank ! 'unknown'
      flag = UNKNOWNDATATYPE
      deallocate(dims, maxDims)
      call MLSMessage ( MLSMsgLevel, ModuleName, &
        & 'Sorry--ReadL1BData_hdf5 has encountered an unknown data type: ' &
        & // trim(Qtype) // Char_rank, MLSFile=L1BFile)
      go to 9
    end select
    if ( DEEBUG ) call dump(l1bData, 0)

9    continue
    ! No need to deallocate cmdims, cmmaxdims, dims and maxdims
    ! to avoid memory leaks because they're allocatable
    call trace_end ( "ReadL1BData_MF_hdf5", &
      & cond=toggle(gen) .and. levels(gen) > 1 )
  end subroutine ReadL1BData_MF_hdf5

  ! -------------------------------------------  Reshape_For_HDF4  -----
  subroutine Reshape_For_HDF4 ( L1bData )

    ! Here's what we'll assume:
    ! If the true rank of l1bdata is 3 or 0, no need to reshape
    ! If it's 1 or 2 then one of the following apply
    !  hdfVERSION         rank 1       rank 2(a)     rank 2(b)
    !      4               1 1 5         10 1 5       1 10 5
    !      5               5 1 1         10 5 1       10 5 1
    !   reorder            2 3 1         1  3 2       3  1 2
    ! (depending on whether it's (a) or (b) methods that apply;
    !   not sure which so we'll code for either; uh-oh, looking like it's (b))
    ! Arguments
    use MLSFillValues, only: IsFillValue
    type(l1bdata_t), intent(inout) :: L1BDATA ! Result
    ! Local variables
    character(len=MaxCharFieldLen), &
      &        dimension(:,:,:), pointer :: CharField => NULL()
    real(r8),  dimension(:,:,:), pointer :: DpField => NULL()
    integer,   dimension(:,:,:), pointer :: IntField => NULL()
    integer, dimension(3)                :: old_shape
    integer, dimension(3)                :: new_shape
    integer, dimension(3)                :: new_order
    integer, dimension(3, 3), parameter  :: reorder = &
      & reshape( &
      & source = (/ 2, 3, 1, 1, 3, 2, 1, 2, 3 /), &
      & shape = (/ 3, 3 /) &
      & )
!     & source = (/ 2, 3, 1, 1, 3, 2, 3, 1, 2 /), &  ! This transposed desired
    integer                              :: status   ! result (don't know why)
    logical, parameter           :: DEEBUG = .FALSE.
    logical                      :: FillInNew
    logical                      :: FillInOld
    character(len=1), parameter  :: method = 'b'
    integer                      :: mord
    integer                      :: S ! Size in bytes of object to deallocate

    ! Executable
    if ( DEEBUG ) then
      print *, 'l1b data_type: ', trim(l1bData%data_type)
      print *, 'rank: ', l1bData%TrueRank
    end if
    if ( l1bData%TrueRank > 2 .or. l1bData%TrueRank < 1 ) return
    new_shape = 1
    if ( method == 'a' .or. l1bData%TrueRank == 1 ) then
      new_order = reorder(:,l1bData%TrueRank)
      mord = 1
    else
      new_order = reorder(:,l1bData%TrueRank+1)
      mord = 2
    end if
    if ( DEEBUG ) then
      print *, 'new_order: ', new_order
      print *, 'mord: ', mord
      print *, '     o l d   f o r m'
      call DumpL1BData(l1bData)
    end if
    select case (l1bData%data_type)
    case('character')
      if ( .not. associated(l1bData%CharField) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 found ' // trim(l1bData%data_type) // &
        & 'unassociated')
      old_shape = shape(l1bData%CharField)
      new_shape(3) = old_shape(l1bData%TrueRank)
      if ( l1bData%TrueRank == 2 ) new_shape(mord) = old_shape(1)
      allocate(CharField(new_shape(1), new_shape(2), new_shape(3)), stat=status)
      call test_allocate ( status, moduleName, trim(l1bData%data_type), &
        & uBounds = new_shape(1:3), elementSize = storage_size(charField) / 8 )
      CharField = reshape(source=l1bData%CharField, &
        & shape=new_shape, order=new_order)
      s = size(l1bData%CharField) * storage_size(l1bData%CharField) / 8
      deallocate(l1bData%CharField, stat=status)
      call test_deallocate ( status, moduleName, trim(l1bData%data_type), s )
      l1bData%CharField => CharField
    case('double')
      if ( .not. associated(l1bData%DpField) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 found ' // trim(l1bData%data_type) // &
        & 'unassociated')
      old_shape = shape(l1bData%DpField)
      new_shape(3) = old_shape(l1bData%TrueRank)
      FillInOld = any(isFillValue(l1bData%DpField))
      if ( l1bData%TrueRank == 2 ) new_shape(mord) = old_shape(1)
      allocate(DpField(new_shape(1), new_shape(2), new_shape(3)), stat=status)
      call test_allocate ( status, moduleName, trim(l1bData%data_type), &
        & uBounds = new_shape(1:3), elementSize = storage_size(dpField) / 8 )
      DpField = reshape(source=l1bData%DpField, &
        & shape=new_shape, order=new_order)
      s = size(l1bData%DpField) * storage_size(l1bData%DpField) / 8
      deallocate(l1bData%DpField, stat=status)
      call test_deallocate ( status, moduleName, trim(l1bData%data_type), s )
      l1bData%DpField => DpField
      FillInNew = any(isFillValue(l1bData%DpField))
      if ( (FillInOld .neqv. FillInNew) .and. DEEBUG ) then
        call outputNamedValue( 'Fill values before reshaping', FillInOld )
        call outputNamedValue( 'Fill values after reshaping', FillInNew )
        call outputNamedValue( 'old shape', old_shape )
        call outputNamedValue( 'new shape', new_shape )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 has iuntroduced or removed Fill values')
      endif
      if ( DEEBUG ) then
        print *, 'old_shape: ', old_shape
        print *, 'new_shape: ', new_shape
      end if
    case('integer')
      if ( .not. associated(l1bData%IntField) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 found ' // trim(l1bData%data_type) // &
        & 'unassociated')
      old_shape = shape(l1bData%IntField)
      new_shape(3) = old_shape(l1bData%TrueRank)
      if ( l1bData%TrueRank == 2 ) new_shape(mord) = old_shape(1)
      allocate(IntField(new_shape(1), new_shape(2), new_shape(3)), stat=status)
      call test_allocate ( status, moduleName, trim(l1bData%data_type), &
        & uBounds = new_shape(1:3), elementSize = storage_size(intField) / 8 )
      IntField = reshape(source=l1bData%IntField, &
        & shape=new_shape, order=new_order)
      s = size(l1bData%IntField) * storage_size(l1bData%IntField) / 8
      deallocate(l1bData%IntField, stat=status)
      call test_deallocate ( status, moduleName, trim(l1bData%data_type), s )
      l1bData%IntField => IntField
    case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 has encountered an unknown data type: ' &
        & // trim(l1bData%data_type))
    end select
    if ( DEEBUG ) then
      print *, '     n e w   f o r m'
      call DumpL1BData(l1bData)
    end if
  end subroutine Reshape_For_HDF4

  ! ----------------------------------------------------  cpField  -----
  subroutine cpField(l1bDataIn, nIn, l1bdataOut, nOut, m)
    ! copy data from l1bDataIn to l1bdataOut starting from index n(In, Out)
    ! where the starting point may be different for the two
    ! Optionally continue copying for m indices
    ! (thus default is m==1)
    ! Arguments
    type (L1BData_T), intent(in)    :: l1bDataIn
    type (L1BData_T), intent(inout) :: l1bDataOut
    integer, intent(in)             :: nIn
    integer, intent(in)             :: nOut
    integer, optional, intent(in)   :: m
    ! Internal variables
    integer :: nn  ! number of indices to copy
    integer :: c1, c2, m1, m2 ! starting, ending channel nums and mif nums
    ! Executable
    nn = 0
    if ( present(m) ) nn=m-1
    ! print *, 'Copying from to ', nIn, nOut
    if ( associated(l1BDataIn%charField) .and. &
      & associated(l1BDataOut%charField) ) then
        c1 = 1
        c2 = min( size(l1BDataIn%charField, 1), size(l1BDataOut%charField, 1) )
        m1 = 1
        m2 = min( size(l1BDataIn%charField, 2), size(l1BDataOut%charField, 2) )
        l1BDataOut%charField(c1:c2,m1:m2,nOut:nOut+nn) = &
        & l1BDataIn%charField(c1:c2,m1:m2,nIn:nIn+nn)
    end if
    if ( associated(l1BDataIn%intField) .and. &
      & associated(l1BDataOut%intField) ) then
        c1 = 1
        c2 = min( size(l1BDataIn%intField, 1), size(l1BDataOut%intField, 1) )
        m1 = 1
        m2 = min( size(l1BDataIn%intField, 2), size(l1BDataOut%intField, 2) )
        l1BDataOut%intField(c1:c2,m1:m2,nOut:nOut+nn) = &
        & l1BDataIn%intField(c1:c2,m1:m2,nIn:nIn+nn)
    end if
    if ( associated(l1BDataIn%dpField) .and. &
      & associated(l1BDataOut%dpField) ) then
        if ( nOut+nn > size(l1BDataOut%dpField, 3) ) then
          call output('size(l1BDataOut%dpField, 3) ', advance='no')
          call output(size(l1BDataOut%dpField, 3), advance='yes')
          call output('nOut+nn ', advance='no')
          call output(nOut+nn, advance='yes')
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Tried to copy past size of output field in cpField')
        else if ( nIn+nn > size(l1BDataIn%dpField, 3) ) then
          call output('size(l1BDataIn%dpField, 3) ', advance='no')
          call output(size(l1BDataIn%dpField, 3), advance='yes')
          call output('nIn+nn ', advance='no')
          call output(nIn+nn, advance='yes')
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Tried to copy past size of input field in cpField')
        end if
        c1 = 1
        c2 = min( size(l1BDataIn%dpField, 1), size(l1BDataOut%dpField, 1) )
        m1 = 1
        m2 = min( size(l1BDataIn%dpField, 2), size(l1BDataOut%dpField, 2) )
        l1BDataOut%dpField(c1:c2,m1:m2,nOut:nOut+nn) = &
        & l1BDataIn%dpField(c1:c2,m1:m2,nIn:nIn+nn)
    ! print *, l1BDataIn%dpField(1,1,nIn:nIn+nn)
    end if
  end subroutine cpField

  ! --------------------------------------------------  ZeroField  -----
  subroutine ZeroField(l1bData, n, value, m, chValue)
    ! Set l1bData to zero starting from index n
    ! Optionally continue zeroing for m indices
    ! (thus default is m==1)
    ! Arguments
    type (L1BData_T), intent(inout)         :: l1bData
    integer, intent(in)                     :: n
    real, optional, intent(in)              :: value   ! default is 0.
    integer, optional, intent(in)           :: m
    character(len=1), optional, intent(in)  :: chValue ! default is ' '
    ! Internal variables
    real             :: zero
    character(len=1) :: blank
    integer          :: nn  ! number of indices to copy
    ! Executable
    nn = 0
    if ( present(m) ) nn=m-1
    zero = 0.
    if ( present(value) ) zero = value
    blank = ' '
    if ( present(chValue) ) blank = chValue
    if ( associated(l1BData%charField)  ) then
        l1BData%charField(:,:,n:n+nn) = blank
    end if
    if ( associated(l1BData%intField)  ) then
        l1BData%intField(:,:,n:n+nn) = int(zero)
    end if
    if ( associated(l1BData%dpField)  ) then
        l1BData%dpField(:,:,n:n+nn) = zero
    end if
  end subroutine ZeroField

  ! ---------------------------------------------------  MyMinval  -----
  function MyMinval ( ints ) result ( mymin )
  ! find minimum value of integer array ints
  ! subject to constraint that none is < 0
  ! (so we avoid picking any undefined, i.e. -999)
  ! Args
    integer, dimension(:), intent(in)  :: ints
    integer                            :: mymin
    ! Internal variables
    integer :: i
    ! Executable
    mymin= maxval(ints)
    if ( size(ints) < 2 ) return
    do i=1, size(ints)
      if ( ints(i) > 0 ) mymin=min(mymin, ints(i))
    end do
  end function MyMinval

  ! ---------------------------------------------  Announce_Error  -----
  subroutine Announce_Error ( lcf_where, full_message, use_toolkit, &
    & error_number )
    ! Arguments

    integer, intent(in)    :: lcf_where
    character(LEN=*), intent(in)    :: full_message
    logical, intent(in), optional :: use_toolkit
    integer, intent(in), optional    :: error_number
    ! Local
    !    character (len=80) :: msg, mnemonic
    !    integer :: status
    logical :: just_print_it
    logical, parameter :: default_output_by_toolkit = .true.

    if ( present(use_toolkit) ) then
      just_print_it = .not. use_toolkit
    else if ( default_output_by_toolkit ) then
      just_print_it = .false.
    else
      just_print_it = .true.
    end if

    if ( .not. just_print_it ) then
      error = max(error,1)
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
        call print_source ( where(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( ': ' )
      call output ( "The " );
      if ( lcf_where > 0 ) then
        call dump_tree_node ( lcf_where, 0 )
      else
        call output ( '(no lcf tree available)' )
      end if

      call output ( " Caused the following error: ", advance='yes', &
        & from_where=ModuleName )
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName )
      if ( present(error_number) ) then
        call output ( 'error number ', advance='no' )
        call output ( error_number, places=9, advance='yes' )
      end if
    else
      call output ( '***Error in module ' )
      call output ( ModuleName, advance='yes' )
      call output ( trim(full_message), advance='yes' )
      if ( present(error_number) ) then
        call output ( 'Error number ' )
        call output ( error_number, advance='yes' )
      end if
    end if
  end subroutine Announce_Error

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------
end module L1BData

! $Log$
! Revision 2.128  2020/01/27 18:28:28  pwagner
! Simplified AssembleL1BQtyName
!
! Revision 2.127  2020/01/16 00:18:45  pwagner
! Can now convert l1boa ds names between file hdf versions
!
! Revision 2.126  2020/01/09 22:24:05  pwagner
! Extra steps so AssembleL1BQtyName wont munge sids-related DS names
!
! Revision 2.125  2019/05/13 20:55:49  pwagner
! Prints less unless debugging
!
! Revision 2.124  2018/08/03 23:23:08  vsnyder
! Announce which file is being checked in CheckForCorruptFileDatabase.
! Instead of ignoring errors in ReadL1BData_MF_... if NeverFail is present
! and true, reduce message level from MLSMSG_Error to MLSMSG_Warning.
!
! Revision 2.123  2018/04/19 23:41:09  pwagner
! Dont crash in ReadL1BData if NeverFail set
!
! Revision 2.122  2018/04/19 02:00:36  vsnyder
! Compute address for allocate/deallocate tracking.  Remove USE statements for
! unused names.
!
! Revision 2.121  2018/03/05 19:25:38  pwagner
! Reduce non-debug printing
!
! Revision 2.120  2017/11/03 19:59:37  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.119  2017/10/18 22:46:52  pwagner
! Dont crash if a DS has more/fewer MAFs than counterMAF
!
! Revision 2.118  2017/02/09 23:46:46  pwagner
! Made more uses CamelCase
!
! Revision 2.117  2016/10/19 00:09:15  pwagner
! Added ConvertL1BData
!
! Revision 2.116  2016/10/11 23:27:18  pwagner
! Commented-out the crash_burn
!
! Revision 2.115  2016/10/05 20:13:58  pwagner
! Implemented Au (Gold) option
!
! Revision 2.114  2016/08/12 00:35:18  pwagner
! Seems to restore tthe gold brick
!
! Revision 2.113  2016/08/09 18:15:29  pwagner
! Made FindMaxMAF generic; survives encounter with non-satellite data files
!
! Revision 2.112  2016/07/28 19:23:10  pwagner
! Fixed error in GetAllHDF5GroupNames
!
! Revision 2.111  2016/07/28 01:42:27  vsnyder
! Refactoring dump and diff
!
! Revision 2.110  2016/07/22 20:03:52  pwagner
! Fixed errors in AssembleL1BQtyName, ReadL1BData_MLSFile
!
! Revision 2.109  2016/07/22 00:22:58  pwagner
! Improved pattern recognition
!
! Revision 2.108  2016/07/21 20:27:06  pwagner
! Can now handle ASMLS data better
!
! Revision 2.107  2016/04/20 00:05:15  pwagner
! cpL1BData now copies l1b datasets between files
!
! Revision 2.106  2016/03/23 00:19:31  pwagner
! DiffL1BData now able to print name on each line
!
! Revision 2.105  2015/07/31 20:40:34  pwagner
! Fixed error added with last commit
!
! Revision 2.104  2015/07/14 23:19:32  pwagner
! May debug case where gap in counterMAF causes Fill values in dp Field
!
! Revision 2.103  2015/03/28 01:07:53  vsnyder
! Changed the kind of dims to allow using ALlocate_Test.  Some spiffing.
! Added stuff to trace allocate/deallocate addresses -- mostly commented out
! because NAG build 1017 doesn't yet allow arrays as arguments to C_LOC.
!
! Revision 2.102  2014/09/04 23:46:21  vsnyder
! More complete and accurate allocate/deallocate size tracking.  Plug a
! memory leak using a final subroutine for L1BData_T.  Add some tracing.
! Convert some local pointers to allocatables.
!
! Revision 2.101  2014/07/23 21:57:59  pwagner
! Attempted to match names passed to allocate/deallocate
!
! Revision 2.100  2014/03/07 19:12:49  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 2.99  2014/01/09 00:25:06  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.98  2013/09/24 23:27:14  vsnyder
! Use Get_Where or Print_Source to start error messages
!
! Revision 2.97  2013/08/31 01:24:53  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end
!
! Revision 2.96  2013/06/29 00:16:37  pwagner
! Added -SL1bread switch
!
! Revision 2.95  2013/05/31 00:40:31  vsnyder
! Make test for SC case insensitive
!
! Revision 2.94  2012/09/11 18:53:40  pwagner
! Anonymized diffed arrays
!
! Revision 2.93  2012/09/05 21:41:20  pwagner
! Removed calls to suspend, resumeOutput when Silent
!
! Revision 2.92  2012/06/07 23:57:17  pwagner
! Limit printing while reading attributes to DEEBUGging
!
! Revision 2.91  2011/07/15 23:31:06  pwagner
! Can now read 4d l1b data types; removed redundant code
!
! Revision 2.90  2011/07/07 00:30:03  pwagner
! Treats diffs of l1bdata types with periods
!
! Revision 2.89  2011/02/05 01:37:05  pwagner
! Passes options to dump routines
!
! Revision 2.88  2010/03/12 21:12:07  pwagner
! Note when diffing if integer arrays are equal
!
! Revision 2.87  2010/01/11 18:33:19  pwagner
! Changed a comment
!
! Revision 2.86  2009/10/30 23:03:50  pwagner
! Commented-out more debugging lines
!
! Revision 2.85  2009/08/24 18:34:45  pwagner
! Repaired syntax errors only NAG complained about
!
! Revision 2.84  2009/08/17 16:53:40  pwagner
! May contract L1BData along any of its dimensions
!
! Revision 2.83  2009/08/04 20:43:18  pwagner
! Replaced silent optional arg
!
! Revision 2.82  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.81  2009/06/16 17:15:55  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.80  2008/06/06 22:52:21  pwagner
! EssentiallyEqual moved to MLSFillValues
!
! Revision 2.79  2008/02/08 00:00:04  pwagner
! Coded around another Intel compiler bug
!
! Revision 2.78  2008/02/07 18:48:16  pwagner
! Prevent appearance of non-ascii when diffing char valued l1bdata
!
! Revision 2.77  2008/01/07 21:37:05  pwagner
! Replace DEFAULTUNDEFINEDVALUE with user-settable undefinedValue
!
! Revision 2.76  2007/08/13 17:36:13  pwagner
! Push some procedures onto new MLSCallStack
!
! Revision 2.75  2007/07/17 00:25:50  pwagner
! Deal more gracefully with attempt to read rank 0 datasets
!
! Revision 2.74  2007/06/21 00:49:51  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.73  2007/02/07 20:56:38  pwagner
! Avoid opening hdf5-type groups in hdf4 files
!
! Revision 2.72  2007/01/26 23:58:02  pwagner
! Fixed bug affecting GetL1BFile
!
! Revision 2.71  2006/11/22 18:12:37  pwagner
! New optional args to diff l1b files with different number MAFs
!
! Revision 2.70  2006/04/11 23:13:29  pwagner
! More room needed in dumping counterMAF array
!
! Revision 2.69  2006/02/06 22:54:51  pwagner
! Should print warnings, not bomb if l1b files not found
!
! Revision 2.68  2006/01/26 00:33:09  pwagner
! demoted more use statements from module level to speed Lahey compiles
!
! Revision 2.67  2006/01/04 20:31:18  pwagner
! Diff procedures may keep silent, returning num of diffs only
!
! Revision 2.66  2005/11/18 01:25:16  pwagner
! L1BData%charField no longer of unit length (hope this breaks nothing else)
!
! Revision 2.65  2005/11/17 20:10:44  pwagner
! Can now read 2d and 3d integer-valued l1bdata; charaacter still fails
!
! Revision 2.64  2005/10/22 00:50:30  pwagner
! Warns but continues if attributes sought in hdf4 file
!
! Revision 2.63  2005/10/19 20:47:29  pwagner
! Fixed bug in GetL1BFile (do others remain?)
!
! Revision 2.62  2005/10/18 23:05:11  pwagner
! Added GetL1BFile to get MLSFile with DS, attribute, or group name
!
! Revision 2.61  2005/09/21 23:12:28  pwagner
! Unnecessary changes
!
! Revision 2.60  2005/08/05 20:36:29  pwagner
! L1BFile arg to ReadL1BData now a pointer
!
! Revision 2.59  2005/07/21 23:36:28  pwagner
! Simplified L1bradSetup
!
! Revision 2.58  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.57  2005/06/14 20:37:54  pwagner
! Many changes to accommodate the new fields in MLSFile_T
!
! Revision 2.56  2005/06/01 17:29:06  pwagner
! Disabled some superfluous printing
!
! Revision 2.55  2005/05/31 17:50:20  pwagner
! Began switch from passing file handles to passing MLSFiles
!
! Revision 2.54  2005/05/12 20:44:52  pwagner
! Uses dump0_m/diff to do diff
!
! Revision 2.53  2005/01/12 23:59:45  pwagner
! diff accepts options to print only stats, rms
!
! Revision 2.52  2004/12/21 22:05:14  pwagner
! Removed clunky temparray introduced to bypass memory issue
!
! Revision 2.51  2004/12/14 21:37:11  pwagner
! Some unnecessary debug printing and a temp array
!
! Revision 2.50  2004/08/26 22:34:41  pwagner
! No time to repair overly ambitious PadL1BData; quick patch instituted instead
!
! Revision 2.49  2004/08/19 00:10:43  pwagner
! L2AUX option to ReadL1BData stops warning msgs about counterMAF
!
! Revision 2.48  2004/08/17 23:47:48  pwagner
! Another guard against inappripriate padding
!
! Revision 2.47  2004/08/16 17:04:23  pwagner
! Pads L1bData after reading gappy data
!
! Revision 2.46  2004/08/04 23:19:01  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.45  2004/08/03 17:59:35  pwagner
! Gets DEFAULTUNDEFINEDVALUE from MLSCommon
!
! Revision 2.44  2004/03/12 00:36:33  pwagner
! Added diff subroutine; hdf version default increased to 5
!
! Revision 2.43  2003/09/12 16:37:10  cvuu
! add subroutine to read L1BOA attributes
!
! Revision 2.42  2003/05/07 01:05:57  vsnyder
! Remove three duplicated lines
!
! Revision 2.41  2003/05/05 23:00:04  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.40  2003/04/02 23:52:34  pwagner
! Checks for FILENOTFOUND
!
! Revision 2.39  2003/04/02 00:38:09  pwagner
! Cater for L1B quantities of rank > 3
!
! Revision 2.38  2003/03/10 17:27:10  pwagner
! Fixed same bug aslast; is it really fixed this time?
!
! Revision 2.37  2003/03/07 00:35:15  pwagner
! Fixed bad bug in reading hdf5 files with FirstMAF > 1
!
! Revision 2.36.2.1  2003/03/27 23:21:31  vsnyder
! Cater for L1B quantities of rank > 3, cosmetic changes
!
! Revision 2.36  2003/02/12 21:49:06  pwagner
! Some small fixes; seems to work with latest hdf5-formats
!
! Revision 2.35  2002/12/11 22:22:15  pwagner
! broadened error check on sd_id to any value lt 1
!
! Revision 2.34  2002/12/10 00:42:15  pwagner
! Stopped printing debugging stuff
!
! Revision 2.33  2002/12/06 01:10:08  pwagner
! No longer transposes hdf5-formatted 2d arrays
!
! Revision 2.32  2002/12/05 19:43:47  pwagner
! Unsuccessful at reading hdf5 so far; speeds up compiling tree-walker, however
!
! Revision 2.31  2002/12/04 01:15:11  pwagner
! Many changes; closer to be able to read hdf5 successfully
!
! Revision 2.30  2002/12/02 23:38:32  pwagner
! Should catch unknown data_type errors
!
! Revision 2.29  2002/11/22 21:48:54  pwagner
! Commented out USE of MLSAUXData
!
! Revision 2.28  2002/11/13 01:02:16  pwagner
! Actually reads hdf5 radiances
!
! Revision 2.27  2002/10/29 18:50:03  pwagner
! Modified AssembleL1BQtyName to convert complete hdf4-style names to hdf5
!
! Revision 2.26  2002/10/10 23:50:57  pwagner
! Passed 1st tests to read l1bdata from hdf5
!
! Revision 2.25  2002/10/09 00:05:02  pwagner
! Added trim function to fieldName
!
! Revision 2.24  2002/10/07 23:21:20  pwagner
! replaced undefined dim_sizes with dims
!
! Revision 2.23  2002/10/05 00:09:16  pwagner
! Finished hdf5 version of readl1bdata; untested however ..
!
! Revision 2.22  2002/10/03 23:04:11  pwagner
! hdfVersion now inout instead of intent(in)
!
! Revision 2.21  2002/09/27 23:37:54  pwagner
! More progress toward hdf5-capable l1b files
!
! Revision 2.20  2002/09/27 00:00:39  pwagner
! Began addings hdf5 functionality; incomplete
!
! Revision 2.19  2002/07/23 23:16:06  pwagner
! Added suggested cause of firstMAF error being bad chunk num
!
! Revision 2.18  2002/07/01 23:48:00  vsnyder
! Plug memory leaks, cosmetic changes
!
! Revision 2.17  2002/05/28 22:34:47  livesey
! Bug fix, wasn't properly allocating counterMAF in some circumstances.
!
! Revision 2.16  2002/01/09 23:42:23  pwagner
! Replaced discouraged print statements with favored calls to output
!
! Revision 2.15  2001/11/01 21:02:31  pwagner
! Willing to read l2aux files as if l1brad (untested); alphabetized procedures; added toc
!
! Revision 2.14  2001/10/26 23:11:10  pwagner
! Provides a single dump module interface
!
! Revision 2.13  2001/10/25 23:31:28  pwagner
! Fixed dump; readl1bData takes neverfail option
!
! Revision 2.12  2001/10/23 22:44:15  pwagner
! Added DumpL1BData
!
! Revision 2.11  2001/10/23 17:09:04  pwagner
! Added PRECISIONSUFFIX as parameter
!
! Revision 2.10  2001/10/03 22:50:03  vsnyder
! Add unfound-quantity name to an error message
!
! Revision 2.9  2001/06/01 02:06:45  livesey
! Bug fix with counterMAF
!
! Revision 2.8  2001/05/30 23:51:48  livesey
! New version, uses new HDF, also cleaner
!
! Revision 2.7  2001/05/06 20:53:47  pwagner
! Allocates l1binfo%filenames along with ids
!
! Revision 2.6  2001/05/04 22:51:06  pwagner
! Now sets L1B..FileName components in ..Setup
!
! Revision 2.5  2001/05/03 23:58:38  vsnyder
! Remove conflicts from merge
!
! Revision 2.4  2001/05/03 23:16:57  livesey
! Added use statement for lexer_core
!
! Revision 2.3  2001/05/03 22:32:25  pwagner
! Added L1B..Setup for Rad and OA
!
! Revision 2.2  2001/03/03 00:06:23  livesey
! Added FindL1BData
!
! Revision 2.1  2000/10/04 01:27:34  vsnyder
! Put ONLY clauses into the USE statements
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.11  2000/06/06 21:11:55  vsnyder
! Declare EXTERNAL attribute for HDF's sf* functions
!
