! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L1BData

  ! Reading and interacting with Level 1B data (HDF4)

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use Dump_0, only: DUMP
  use Hdf, only: DFACC_READ, SFSTART, SFGINFO, SFN2INDEX, SFSELECT, &
    & SFRDATA_f90, &
    & SFRCDATA, SFENDACC, DFNT_CHAR8, DFNT_INT32, DFNT_FLOAT64, &
    & DFNT_FLOAT32
  use Lexer_Core, only: PRINT_SOURCE
  use MLSCommon, only: R4, R8, L1BINFO_T, DEFAULTUNDEFINEDVALUE, FILENAMELEN
  use MLSFiles, only: FILENOTFOUND, HDFVERSION_4, HDFVERSION_5, &
    & MLS_HDF_VERSION, MLS_IO_GEN_OPENF
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR, &
    & MLSMSG_WARNING, MLSMSG_L1BREAD, MLSMSG_WARNING, MLSMSG_DEALLOCATE
  use MLSStrings, only: CompressString
  use MLSStringLists, only: NumStringElements
  use MoreTree, only: Get_Field_ID
  use Output_M, only: Output
  use PCFHdr, only: GlobalAttributes
  use SDPToolkit, only: max_orbits
  use String_Table, only: Get_String
  use TREE, only: NSONS, SUB_ROSA, SUBTREE, DUMP_TREE_NODE, SOURCE_REF

  implicit NONE

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (data types and parameters)
! L1BData_T                       Quantities from an L1B data file
! NAME_LEN                        Max length of l1b sds array name
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
! DeallocateL1BData               Called when an l1bData is finished with
! Dump                            Print facts about l1brad quantity
! findl1bdata                     Which file handle contains a given sd name
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
! char*32 AssembleL1BQtyName (char name, int hdfVersion, log isTngtQty,
!                         [char InstrumentName] ) 
! DeallocateL1BData (l1bData_T l1bData) 
! Dump (l1bData_T l1bData, int details)
! int FindL1BData (int files(:), char filedName, [int hdfVersion]) 
! L1boaSetup (int root, L1BInfo_T L1BInfo, int f_file, [int hdfVersion])
! L1bradSetup (int root, L1BInfo_T L1BInfo, int f_file, 
!               int MaxNumL1BRadIDs, int illegalL1BRadID, [int hdfVersion])
! ReadL1BData (int L1FileHandle, char QuantityName, l1bData_T l1bData,
!               int NoMAFs, int Flag, [int FirstMAF], [int LastMAF],
!               [log NeverFail], [int hdfVersion]) 
! ReadL1BAttribute (int L1FileHandle, value(:), nchar AttributeName, 
!               int Flag, [int hdfVersion]) 
! === (end of api) ===

  private

  public :: L1BData_T, NAME_LEN, PRECISIONSUFFIX, &
    & allocateL1BData, AssembleL1BQtyName, ContractL1BData, cpL1bData, &
    & DeallocateL1BData, DIFF, DUMP, &
    & FINDL1BDATA, FindMaxMAF, &
    & L1BRadSetup, L1BOASetup, PadL1BData, &
    & ReadL1BAttribute, ReadL1BData

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private not_used_here
  !---------------------------------------------------------------------------

  interface DIFF
    module procedure DiffL1BData
  end interface

  interface DUMP
    module procedure DumpL1BData
  end interface

  interface ReadL1BAttribute
    module procedure ReadL1BAttribute_intarr1, &
	& ReadL1BAttribute_dblarr1
  end interface

  ! Parameters
  integer, parameter :: NAME_LEN = 64  ! Max len of SDS array name
  ! suffix of sd precision; check against 'grep -i precision l1/OutputL1B.f90'
  character  (len=*), parameter :: PRECISIONSUFFIX = ' precision'
  ! The following should be same as l2/Fill
  real, parameter ::    UNDEFINED_VALUE = DEFAULTUNDEFINEDVALUE ! -999.99 
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

  type L1BData_T
    character (len=name_len) :: L1BName ! Name of field in file
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
    character, dimension(:,:,:), pointer :: CharField => NULL()
    real(r8),  dimension(:,:,:), pointer :: DpField => NULL()
    integer,   dimension(:,:,:), pointer :: IntField => NULL()
    ! all the above dimensioned (noAuxInds,maxMIFs,noMAFs)
    ! logical :: mustPad                  ! Gaps in counterMAF
  end type L1BData_T
  
  integer :: Error            ! Error level -- 0 = OK
  
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
  integer, public, parameter :: CANTALLOCATECHARS =  CANTREADCOUNTERMAF + 1
  integer, public, parameter :: CANTREAD3DFIELD =    CANTALLOCATECHARS + 1
  integer, public, parameter :: UNKNOWNDATATYPE =    CANTREAD3DFIELD + 1
  integer, public, parameter :: CANTENDCOUNTERMAF =  UNKNOWNDATATYPE + 1
  integer, public, parameter :: CANTENDQUANTITY =    CANTENDCOUNTERMAF + 1

contains ! ============================ MODULE PROCEDURES =======================

  !-------------------------------------------  allocateL1BData  -----
  subroutine allocateL1BData ( l1bData, &
    & indims, trueRank, datatype, &
    & L1bDataSibling )
    ! Allocate arrays in l1bData type
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
    elseif( present(L1bDataSibling) ) then
      myDataType = L1bDataSibling%data_type
      rank = L1bDataSibling%trueRank
      noMAFs = L1bDataSibling%noMAFs
      if ( associated(L1bDataSibling%charField)) then
        dims = shape(L1bDataSibling%charField)
      elseif ( associated(L1bDataSibling%intField)) then
        dims = shape(L1bDataSibling%intField)
      elseif ( associated(L1bDataSibling%dpField)) then
        dims = shape(L1bDataSibling%dpField)
      else
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'allocatel1bData was passed L1bDataSibling w/o allocating it' )
      endif
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'allocatel1bData must be passed either datatype or L1bDataSibling' )
    endif
    if ( present(indims) ) then
      dims = indims
      noMAFs = dims(3)
    endif
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
    endif
    call allocate_test ( l1bData%counterMAF, noMAFs, 'l1bData%counterMAF', &
      & ModuleName )
    l1bData%counterMAF = int(DEFAULTUNDEFINEDVALUE)
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
      l1bData%intField = int(DEFAULTUNDEFINEDVALUE)
    case ( 'd', 'r' ) ! double
      call allocate_test ( l1bData%dpField, dims(1), dims(2), dims(3), &
        & 'l1bData%dpField', ModuleName )
      l1bData%dpField = DEFAULTUNDEFINEDVALUE
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
  end subroutine allocateL1BData

  !-------------------------------------------  AssembleL1BQtyName  -----
  function AssembleL1BQtyName ( name, hdfVersion, isTngtQty, &
    & InstrumentName, dont_compress_name) &
    & result(QtyName)
    ! Returns a QtyName to be found in the L1b file
    ! If given InstrumentName, name should be a fragment:
    ! e.g., name='VelECI' and InstrumentName='sc'
    ! Without InstrumentName, name should be complete (in hdf4-style)
    ! e.g., name='scVelECI' or name='R1A:118.B1F:PT.S0.FB25-1'
    ! and isTngtQty is simply ignored
    character(len=*), intent(in) :: name        ! bare name; e.g. SolarZenith
    integer, intent(in)          :: hdfVersion  ! which QtyName must conform to
    logical, intent(in)          :: isTngtQty   ! T or F
    character(len=*), intent(in), optional :: InstrumentName ! e.g. THz
    logical, intent(in), optional :: dont_compress_name
    character(len=NAME_LEN)      :: QtyName
    logical, parameter           :: DEEBUG = .FALSE.

    ! Private
    character(len=1) :: head
    character(len=1) :: instr_tail
    character(len=1) :: tp_tail
    character(len=4) :: my_instrument
    character(len=NAME_LEN) :: the_rest
    logical          :: is_a_signal
    logical          :: compress
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
    if ( hdfVersion == HDFVERSION_5 ) then
      head = '/'
      instr_tail = '/'
      tp_tail = '/'
    else
      head = ''
      instr_tail = '.'
      if ( present(InstrumentName) ) then
        if ( trim(InstrumentName) == 'sc' ) instr_tail = ''
      end if
      tp_tail = ''
    end if
    if ( present(InstrumentName) ) then
      QtyName = head // trim(InstrumentName) // instr_tail
    else if ( is_a_signal .or. hdfVersion /= HDFVERSION_5 ) then
      QtyName = head
    else
      ! Need only to convert complete hdf4-name to hdf5-name
      ! This means we must parse hdf4-name fully, however
      QtyName = head
      ! 1st--is there an instrument prefixed to name?
      if ( name(1:2) == 'sc' ) then
        my_instrument = 'sc'
        the_rest = name(3:)
      else if ( name(1:4) == 'GHz.' ) then
        my_instrument = 'GHz'
        the_rest = name(5:)
      else if ( name(1:4) == 'THz.' ) then
        my_instrument = 'THz'
        the_rest = name(5:)
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
      else if ( DROPTPSUBGROUP ) then
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
      return
    end if
    if ( isTngtQty ) then
      if ( DROPTPSUBGROUP .and. hdfVersion == HDFVERSION_5 ) then
        QtyName = trim(QtyName)
      else
        QtyName = trim(QtyName) // 'tp' // tp_tail
      end if
    end if
    QtyName = trim(QtyName) // trim(name)
    if ( compress ) QtyName = CompressString(QtyName)
    if ( DEEBUG ) then                            
      print *, 'more converted name: ', trim(QtyName)  
    end if                                         
  end function AssembleL1BQtyName

  !-------------------------------------------------  ContractL1BData  -----
  subroutine ContractL1BData ( L1BDataIn, L1BDataOut, noMAFs, firstMAF, lastMAF )
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
    if ( present(firstMAF) ) myFirstMAF = firstMAF
    myLastMAF = l1bDataIn%NoMAFs - l1bDataIn%firstMAF - 1
    if ( present(lastMAF) ) myLastMAF = lastMAF
    if ( associated(l1bDataIn%charField)) then
      dims = shape(l1bDataIn%charField)
    elseif ( associated(l1bDataIn%intField)) then
      dims = shape(l1bDataIn%intField)
    elseif ( associated(l1bDataIn%dpField)) then
      dims = shape(l1bDataIn%dpField)
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Contractl1bData was passed l1bDataIn w/o allocating it' )
    endif
    noMAFs = myLastMAF - myFirstMAF + 1
    dims(3) = noMAFs
    rank = l1bDataIn%trueRank
    l1bDataOut%noMAFs = noMAFs
    l1bDataOut%firstMAF = myFirstMAF
    if ( DEEBug ) print *, 'Preparing to contract ', trim(l1bDataIn%L1BName)
    if ( DEEBug ) print *, 'rank ', rank
    if ( DEEBug ) print *, 'NoMAFs ', NoMAFs
    if ( DEEBug ) print *, 'l1b(in)NoMAFs ', l1bDataIn%NoMAFs
    if ( DEEBug ) print *, 'dims ', dims
    if ( DEEBug ) print *, 'min(counterMAF) ', myminval(l1bDataIn%counterMAF)
    if ( DEEBug ) print *, 'min(dpField) ',minval(l1bDataIn%dpField(1,1,:))
    if ( DEEBug ) print *, 'max(dpField) ',maxval(l1bDataIn%dpField(1,1,:))
    call allocateL1BData ( l1bDataOut, dims, L1bDataSibling=l1bDataIn )
    mafOffSet = firstMAF
    if ( mafOffSet+NoMAFs > size(l1bDataIn%counterMAF) ) &
      & mafOffSet = size(l1bDataIn%counterMAF) - NoMAFs
    do maf=1, NoMAFs
      l1bDataOut%counterMAF(maf) = l1bDataIn%counterMAF(mafOffSet+maf)
      call cpField(l1bDataIn, mafOffSet+maf, l1bdataOut, maf)
    enddo
    l1bDataOut%firstMAFCtr = myminval(l1bDataOut%counterMAF)
    l1bDataOut%lastMAFCtr = maxval(l1bDataOut%counterMAF)
  end subroutine ContractL1BData

  !-------------------------------------------  cpL1BData  -----
  subroutine cpL1BData ( l1bData1, l1bData2, offsetMAF )
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
    endif
    ! print *, 'l1bData2%counterMAF(1): ', l1bData2%counterMAF
    if ( associated(l1bdata1%charField) ) l1bdata2%charField=l1bdata1%charField
    if ( associated(l1bdata1%intField) ) l1bdata2%intField=l1bdata1%intField
    if ( associated(l1bdata1%dpField) ) l1bdata2%dpField=l1bdata1%dpField
  end subroutine cpL1BData

  !-------------------------------------------  DeallocateL1BData  -----
  subroutine DeallocateL1BData ( l1bData )
    ! This should be called when an l1bData is finished with
    type( L1BData_T ), intent(inout) :: L1bData

    ! Executable code
    call deallocate_test ( l1bData%counterMAF, 'l1bData%counterMAF', ModuleName )
    call deallocate_test ( l1bData%charField, 'l1bData%charField', ModuleName )
    call deallocate_test ( l1bData%intField, 'l1bData%intField', ModuleName )
    call deallocate_test ( l1bData%dpField, 'l1bData%dpField', ModuleName )
  end subroutine DeallocateL1BData

  !-------------------------------------------------  DiffL1BData  -----
  subroutine DiffL1BData ( l1bData1, l1bData2, details )
    ! Diiff two l1brad quantities
    type( L1BData_T ), intent(inout) :: L1bData1
    type( L1BData_T ), intent(inout) :: L1bData2
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even counterMAF
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1

    ! Local variables
    integer :: MYDETAILS

    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details
    if ( trim(L1bData1%NameInst) /= trim(L1bData2%NameInst) ) then
      call output(trim(L1bData1%NameInst), advance='yes')
      call output(trim(L1bData2%NameInst), advance='yes')
    endif
    if ( trim(L1bData1%L1BName) /= trim(L1bData2%L1BName) ) then
      call output('L1B rad quantity (1) Name = ', advance='no')
      call output(trim(L1bData1%L1BName), advance='yes')
      call output('L1B rad quantity (2) Name = ', advance='no')
      call output(trim(L1bData2%L1BName), advance='yes')
    endif
    if ( myDetails < -1 ) return
    if ( L1bData1%FirstMAF /= L1bData2%FirstMAF ) then
      call output(' (1) First major frame read = ', advance='no')
      call output(L1bData1%FirstMAF, advance='yes')
      call output(' (2) First major frame read = ', advance='no')
      call output(L1bData2%FirstMAF, advance='yes')
    endif
    if ( L1bData1%NoMAFs /= L1bData2%NoMAFs ) then
      call output(' (1) Num of MAFs read = ', advance='no')
      call output(L1bData1%NoMAFs, advance='yes')
      call output(' (2) Num of MAFs read = ', advance='no')
      call output(L1bData2%NoMAFs, advance='yes')
    endif
    if ( L1bData1%MaxMIFs /= L1bData2%MaxMIFs ) then
      call output(' (1) Max # of MIFs/MAF in SD array = ', advance='no')
      call output(L1bData1%MaxMIFs, advance='yes')
      call output(' (2) Max # of MIFs/MAF in SD array = ', advance='no')
      call output(L1bData2%MaxMIFs, advance='yes')
    endif
    if ( L1bData1%NoAuxInds /= L1bData2%NoAuxInds ) then
      call output(' (1) Num of auxilliary indices = ', advance='no')
      call output(L1bData1%NoAuxInds, advance='yes')
      call output(' (2) Num of auxilliary indices = ', advance='no')
      call output(L1bData2%NoAuxInds, advance='yes')
    endif
    if ( L1bData1%FirstMAFCtr /= L1bData2%FirstMAFCtr ) then
      call output(' (1) First major frame counter = ', advance='no')
      call output(L1bData1%FirstMAFCtr, advance='yes')
      call output(' (2) First major frame counter = ', advance='no')
      call output(L1bData2%FirstMAFCtr, advance='yes')
    endif
    if ( L1bData1%LastMAFCtr /= L1bData2%LastMAFCtr ) then
      call output(' (1) Last major frame counter = ', advance='no')
      call output(L1bData1%FirstMAFCtr, advance='yes')
      call output(' (2) Last major frame counter = ', advance='no')
      call output(L1bData2%LastMAFCtr, advance='yes')
    endif

    if ( myDetails < 0 ) return
    if ( associated(l1bData1%counterMAF) .and. &
      & associated(l1bData2%counterMAF) ) then
      if ( any(l1bData1%counterMAF /= l1bData2%counterMAF)) then
        call dump ( l1bData1%counterMAF - l1bData2%counterMAF, &
          & 'l1bData%counterMAF (diff)' )
       endif
    else
      call output('(CounterMAF arrays not associated)', advance='yes')
    end if

    if ( myDetails < 1 ) return
    if ( associated(l1bData1%charField) .and. &
      & associated(l1bData2%charField)) then
      if ( any(l1bData1%charField /= l1bData2%charField) ) then
        call dump ( l1bData1%CharField, 'l1bData1%CharField' )
        call dump ( l1bData2%CharField, 'l1bData2%CharField' )
      end if
    else
      call output('(CharField arrays not associated)', advance='yes')
    end if

    if ( associated(l1bData1%intField) &
      & .and. associated(l1bData1%intField) ) then
      if ( any(l1bData1%intField /= l1bData2%intField) ) then
        call dump ( l1bData1%intField - l1bData2%intField, &
          & 'l1bData%intField (diff)' )
      endif
    else
      call output('(intField arrays not associated)', advance='yes')
    end if

    if ( associated(l1bData1%dpField) &
      & .and. associated(l1bData1%dpField) ) then
      if ( any(l1bData1%dpField /= l1bData2%dpField) ) &
        & call dump ( l1bData1%dpField-l1bData2%dpField, &
        & 'l1bData%dpField (diff)' )
    else
      call output('(dpField arrays not associated)', advance='yes')
    end if

  end subroutine DiffL1BData

  !-------------------------------------------------  DumpL1BData  -----
  subroutine DumpL1BData ( l1bData, details )
    ! Disclose pertinent, perhaps damning facts about an l1brad quantity
    type( L1BData_T ), intent(inout) :: L1bData
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even counterMAF
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1

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
      call dump ( l1bData%counterMAF, 'l1bData%counterMAF' )
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
        & fillValue = int(DEFAULTUNDEFINEDVALUE) )
    else
      call output('(intField array not associated)', advance='yes')
    end if

    if ( associated(l1bData%dpField) ) then
      call dump ( l1bData%dpField, 'l1bData%dpField', &
        & fillValue=DEFAULTUNDEFINEDVALUE*1.d0 )
    else
      call output('(dpField array not associated)', advance='yes')
    end if
  end subroutine DumpL1BData

  ! ------------------------------------------------  FindL1BData  -----
  integer function FindL1BData ( files, fieldName, hdfVersion )

  use MLSHDF5, only: IsHDF5DSPresent

    integer, dimension(:), intent(in) :: files ! File handles
    character (len=*), intent(in) :: fieldName ! Name of field
    integer, optional, intent(in) :: hdfVersion

    ! Externals
    integer, external :: SFN2INDEX

    ! Local variables
    integer :: i
    integer :: myhdfVersion

    ! Executable code
    if ( present(hdfVersion) ) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L1BDEFAULT_HDFVERSION
    end if

    findL1BData=0
    do i = 1, size(files)
      if ( myhdfVersion == HDFVERSION_4 ) then
        if ( sfn2index(files(i),trim(fieldName)) /= -1 ) then
          findL1BData = files(i)
          return
        end if
      else
        if ( IsHDF5DSPresent(files(i),trim(fieldName)) ) then
          findL1BData = files(i)
          return
        end if
      end if
    end do
  end function FindL1BData

  ! ------------------------------------------------  FindMaxMAF  -----
  integer function FindMaxMAF ( files, hdfVersion, minMAF )
  ! Find maximum MAF among files (using counterMAF arrays)

  use MLSHDF5, only: IsHDF5DSPresent

    integer, dimension(:), intent(in) :: files ! File handles
    integer, optional, intent(in) :: hdfVersion
    integer, optional, intent(out) :: minMAF

    ! Externals
    integer, external :: SFN2INDEX

    ! Local variables
    integer :: i
    integer :: myhdfVersion
    character(len=*), parameter :: fieldname = 'counterMAF'
    type(L1BData_T) :: l1bData
    integer :: noMAFs
    integer :: myMinMAF
    integer :: status
    logical :: haveCtrMAF
    logical, parameter :: DEEBug = .false.
    ! Executable code
    if ( present(hdfVersion) ) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L1BDEFAULT_HDFVERSION
    end if

    FindMaxMAF=0
    myMinMAF = BIGGESTMAFCTR
    haveCtrMAF = .false.
    do i = 1, size(files)
      if ( myhdfVersion == HDFVERSION_4 ) then
        haveCtrMAF = sfn2index(files(i),trim(fieldName)) /= -1
        if ( haveCtrMAF ) then
          call ReadL1BData ( files(i), fieldName, L1bData, noMAFs, status, &
            & hdfVersion=HDFVERSION_4, dontPad=.true.)
        end if
      else
        haveCtrMAF = IsHDF5DSPresent(files(i),trim(fieldName))
        if ( haveCtrMAF ) then
          call ReadL1BData ( files(i), fieldName, L1bData, noMAFs, status, &
            & hdfVersion=HDFVERSION_5, dontPad=.true.)
        end if
      end if
      if ( haveCtrMAF ) then
        FindMaxMAF = max(FindMaxMAF, maxval(l1bData%counterMAF))
        myMinMAF = min(myMinMAF, myminval(l1bData%counterMAF))
        if ( DEEBug ) print *, 'counterMAF ', l1bData%counterMAF
        call deallocatel1bdata(L1bData)
      else
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to find '//trim(fieldName)//' in l1b files')
      endif
    end do
    if ( DEEBug ) print *, 'FindMaxMAF ', FindMaxMAF
    if ( DEEBug ) print *, 'myMinMAF ', myMinMAF
    if ( present(minMAF) ) minMAF = myMinMAF
  end function FindMaxMAF

  !--------------------------------------------------  IsL1BGappy  -----
  logical function IsL1BGappy ( l1bData, ignoreGlobalAttrs )
    ! Look for gaps in l1bData, returning true if any found

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
    endif
    if ( IsL1BGappy ) return
    ! Search for interior gaps
    IsL1BGappy = .true. ! If exit from next loop prematurely, must be gappy
    nextMAF = minCtrMAF
    do i=1, size(l1BData%counterMAF)
      if ( nextMAF /= l1BData%counterMAF(i)) return
    enddo
    IsL1BGappy = .false.
  end function IsL1BGappy

  !--------------------------------------------------  L1BOASetup  -----
  subroutine L1boaSetup ( root, l1bInfo, F_FILE, hdfVersion )
    ! Take file name from l2cf, open, and store unit no. in l1bInfo

    ! Dummy arguments
    type (L1BInfo_T) :: L1BINFO         ! File handles etc. for L1B dataset
    integer, intent(in) :: ROOT         ! of the l1brad file specification.
    integer, intent(in) :: F_FILE       ! From init_tables_module
    integer, optional, intent(inout) :: hdfVersion

    ! Local variables

    character(len=FileNameLen) :: FileName

    integer :: I                        ! Loop inductor, subscript
    integer :: record_length
    integer :: SON                      ! Some subtree of root.
    integer :: SD_ID                    ! From HDF
    integer :: the_hdf_version

    ! Executable code
    error = 0
    ! Collect data from the fields. (only one legal field: file='...')
    do i = 2, nsons(root)
      son = subtree(i,root)
      if ( get_field_id(son) == f_file ) then
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        the_hdf_version = mls_hdf_version(FileName)
        if ( the_hdf_Version == FILENOTFOUND ) &
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'File not found; make sure the name and path are correct' &
            & // trim(fileName) )
        if ( present(hdfVersion) ) then
          sd_id = mls_io_gen_openF('hg', .true., error, &
            & record_length, DFACC_READ, &
            & FileName, hdfVersion=the_hdf_version, debugOption=.false.)
        else
          sd_id = mls_io_gen_openF('hg', .true., error, &
            & record_length, DFACC_READ, &
            & FileName)
        end if
        if ( sd_id <= 0 ) then
          call announce_error ( son, &
            & 'Error opening L1BOA file: ' //Filename)
        else
          l1bInfo%L1BOAID = sd_id
          l1bInfo%L1BOAFileName = Filename
          if ( present(hdfVersion) ) hdfVersion = the_hdf_version
        end if
      else
        call announce_error ( son, &
          & 'Unknown field specified in read l1boa' )
      end if
    end do
  end subroutine L1boaSetup

  ! ------------------------------------------------- L1BRadSetup  -----
  subroutine L1bradSetup ( Root, L1bInfo, F_File, MaxNumL1BRadIDs, &
    & illegalL1BRadID, hdfVersion )
    ! Take file name from l2cf, open, and store unit no. in l1bInfo
    ! Dummy arguments
    type (L1BInfo_T) :: L1BINFO         ! File handles etc. for L1B dataset
    integer, intent(in) :: ROOT         ! of the l1brad file specification.
    integer, intent(in) :: F_FILE
    integer, intent(in) :: MAXNUML1BRADIDS
    integer, intent(in) :: ILLEGALL1BRADID
    integer, optional, intent(inout) :: hdfVersion

    ! Local variables
    character(len=FileNameLen) :: FILENAME

    integer :: I                        ! Loop inductor, subscript
    integer :: SON                      ! Some subtree of root.
    integer :: STATUS                   ! Flag
    integer :: SD_ID                    ! ID from HDF

    integer, save :: IFL1 = 0           ! num. of L1brad files opened so far
    integer :: record_length, the_hdf_version

    ! Executable code
    error = 0

    ! Collect data from the fields. (only one legal field: file='...')
    do i = 2, nsons(root)
      son = subtree(i,root)
      if ( get_field_id(son) == f_file ) then
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        if ( .NOT. associated(l1bInfo%L1BRADIDs) ) then
          allocate ( l1bInfo%L1BRADIDs(MAXNUML1BRADIDS), stat=status )
          allocate ( l1bInfo%L1BRADFileNames(MAXNUML1BRADIDS), stat=status )
          l1bInfo%L1BRADIDs = ILLEGALL1BRADID
          if ( status /= 0 ) &
            & call announce_error ( son, 'Allocation failed for l1bInfo' )
        endif
        the_hdf_version = mls_hdf_version(FileName)
        if ( the_hdf_Version == FILENOTFOUND ) &
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'File not found; make sure the name and path are correct' &
            & // trim(fileName) )

        if ( present(hdfVersion) ) then
          sd_id = mls_io_gen_openF('hg', .true., error, &
            & record_length, DFACC_READ, &
            & FileName, hdfVersion=the_hdf_version, debugOption=.false.)
        else
          sd_id = mls_io_gen_openF('hg', .true., error, &
            & record_length, DFACC_READ, &
            & FileName)
        end if
        if ( sd_id <= 0 ) then
          call announce_error ( son, &
            & 'Error opening L1BRAD file: ' //Filename)
        else if ( ifl1 == MAXNUML1BRADIDS ) then
          call announce_error ( son, "Cannot open any more L1BRAD files" )
          exit
        else
          ifl1 = ifl1 + 1
          l1bInfo%L1BRADIDs(ifl1) = sd_id
          l1bInfo%L1BRADFileNames(ifl1) = Filename
          if ( present(hdfVersion) ) hdfVersion = the_hdf_version
        end if
      else
        call announce_error ( son, &
          & 'Unknown field specified in read l1brad' )
      end if
    end do

  end subroutine L1bradSetup

  !-------------------------------------ReadL1BAttribute_intarr1l---------
  subroutine ReadL1BAttribute_intarr1 ( L1FileHandle, value, AttrName, Flag, &
     & hdfVersion )
    
    use MLSHDF5, only: IsHDF5AttributePresent, GetHDF5Attribute
    use HDF5, only: H5GCLOSE_F, H5GOPEN_F

    ! Dummy arguments
    integer, intent(in)            :: L1FileHandle ! From HDF
    integer, intent(out) :: value(:) ! Result
    character(len=*), intent(in) :: AttrName ! attribute name to retrieve
    integer, intent(out) :: Flag        ! Error flag
    integer, optional, intent(in) :: hdfVersion

    ! Local variables
    integer :: myhdfVersion
    integer :: aID, status

    ! Executable code
    if ( present(hdfVersion) ) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L1BDEFAULT_HDFVERSION
    end if

    Flag = 0

    if ( myhdfVersion == HDFVERSION_4 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Not implement in l1boa file')
	Flag = -1
    else 
	call h5gOpen_f (L1FileHandle,'/', aID, status)
        if ( status /= 0 ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to open group attribute in l1boa file' )
	   Flag = -1
	end if	
     	if ( .not. IsHDF5AttributePresent(aID, AttrName) ) then
     	   Flag = -1
           call MLSMessage ( MLSMSG_Warning, ModuleName, &
		& 'Failed to find attribute in l1boa file'//AttrName)
	else 
           call output ('get attribute', advance='no')
           call output (AttrName, advance='yes')
           call GetHDF5Attribute(aID, AttrName, value)
    	end if
	call h5gClose_f (aID, status)
        if ( status /= 0 ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to close group attribute in l1boa file' )
	end if	
    end if
  end subroutine ReadL1BAttribute_intarr1

  !-------------------------------------ReadL1BAttribute_dbtarr1l---------
  subroutine ReadL1BAttribute_dblarr1 ( L1FileHandle, value, AttrName, Flag, &
     & hdfVersion )
    
    use MLSHDF5, only: IsHDF5AttributePresent, GetHDF5Attribute
    use HDF5, only: H5GCLOSE_F, H5GOPEN_F

    ! Dummy arguments
    integer, intent(in)            :: L1FileHandle ! From HDF
    real(r8), intent(out) :: value(:) ! Result
    character(len=*), intent(in) :: AttrName   ! attribute name to retrieve
    integer, intent(out) :: Flag        ! Error flag
    integer, optional, intent(in) :: hdfVersion

    ! Local variables
    integer :: myhdfVersion
    integer :: aID, status

    ! Executable code
    if ( present(hdfVersion) ) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L1BDEFAULT_HDFVERSION
    end if

    Flag = 0

    if ( myhdfVersion == HDFVERSION_4 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Not implement in l1boa file')
	Flag = -1
    else 
	call h5gOpen_f (L1FileHandle,'/', aID, status)
        if ( status /= 0 ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to open group attribute in l1boa file' )
	   Flag = -1
	end if	
     	if ( .not. IsHDF5AttributePresent(aID, AttrName) ) then
     	   Flag = -1
           call MLSMessage ( MLSMSG_Warning, ModuleName, &
		& 'Failed to find attribute in l1boa file'//AttrName)
	else 
           call output ('get attribute', advance='no')
           call output (AttrName, advance='yes')
           call GetHDF5Attribute(aID, AttrName, value)
    	end if
	call h5gClose_f (aID, status)
        if ( status /= 0 ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to close group attribute in l1boa file' )
	end if	
    end if
  end subroutine ReadL1BAttribute_dblarr1

  !-------------------------------------------------  PadL1BData  -----
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
    elseif ( associated(l1bDataIn%intField)) then
      dims = shape(l1bDataIn%intField)
    elseif ( associated(l1bDataIn%dpField)) then
      dims = shape(l1bDataIn%dpField)
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Padl1bData was passed l1bDataIn w/o allocating it' )
    endif
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
    enddo
    call zeroField(l1bdataOut, 1, DEFAULTUNDEFINEDVALUE, m=NoMAFs)
    call cpField(l1bdataIn, 1, l1bdataOut, 1, m=oldSize)
  end subroutine PadL1BData

  !-------------------------------------------------  BadPadL1BData  -----
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
    logical :: myForce
    integer, dimension(3) :: dims
    integer :: gap
    integer :: i
    integer :: indexOut
    integer :: lastMAF
    integer :: cmindex
    integer :: current
    integer :: maf
    integer :: maxMAF
    integer :: rank
    logical, parameter :: DEEBug = .false.
    ! Executable
    myForce = .false.
    if ( present(force) ) myForce = force
    ! Check assumption that noMafs >= l1bdatain%NoMAFs
    if ( .not. myForce .and. noMAFs < l1bdatain%NoMAFs ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'noMAFs requested smaller than number stored in input l1bData' )
    elseif ( present(PrecisionOut) .and. .not. present(PrecisionIn) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Must have PrecisionIn to calculate PrecisionOut' )
    endif
    if ( associated(l1bDataIn%charField)) then
      dims = shape(l1bDataIn%charField)
    elseif ( associated(l1bDataIn%intField)) then
      dims = shape(l1bDataIn%intField)
    elseif ( associated(l1bDataIn%dpField)) then
      dims = shape(l1bDataIn%dpField)
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'BadPadL1BData was passed l1bDataIn w/o allocating it' )
    endif
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
    enddo
    ! Now hunt for missing MAFs
    maf = FirstMAFCtr
    lastmaf = maf - 1
    indexOut = 0
    if ( DEEBug ) print *, 'size(l1bDataOut%counterMAF) ', size(l1bDataOut%counterMAF)
    do cmindex = 1, l1bDataIn%noMAFs
      if ( DEEBug ) print *, 'l1bDataIn%counterMAF(cmindex) ', l1bDataIn%counterMAF(cmindex)
      if ( DEEBug ) print *, 'maf ', maf
      if ( DEEBug ) print *, 'indexOut ', indexOut
      if ( l1bDataIn%counterMAF(cmindex) < maf ) then
        ! counterMAF too small: must ignore (perhaps past end of array?)
      elseif ( l1bDataIn%counterMAF(cmindex) == maf ) then
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
          call zeroField(l1bdataOut, indexOut+i, DEFAULTUNDEFINEDVALUE)
          if ( present(PrecisionOut) ) &
            & call zeroField(PrecisionOut, indexOut+i, DEFAULTUNDEFINEDVALUE)
        enddo
        if ( DEEBug ) print *, 'Now treat current normally ', current
        ! Now treat current normally
        indexOut = indexOut + gap + 1
        l1bDataOut%counterMAF(indexOut) = current
        call cpField(l1bDataIn, cmindex, l1bdataOut, indexOut)
        if ( present(PrecisionOut) ) &
          &  call cpField(PrecisionIn, cmindex, PrecisionOut, indexOut)
        maf = current+1
      endif
    enddo
    if ( maf > FirstMAFCtr + noMAFs - 1 ) return
    gap = FirstMAFCtr + noMAFs - maf
    ! Apparently, we end too soon, so must pad
    if ( DEEBug ) print *, 'Apparently, we end too soon, so must pad ', gap, indexOut
    do i=1, gap
      l1bDataOut%counterMAF(indexOut+i) = maf-1+i
      call zeroField(l1bdataOut, indexOut+i, DEFAULTUNDEFINEDVALUE)
      if ( present(PrecisionOut) ) &
        & call zeroField(PrecisionOut, indexOut+i, DEFAULTUNDEFINEDVALUE)
    enddo
  end subroutine BadPadL1BData

  !-------------------------------------------------  ReadL1BData  -----
  subroutine ReadL1BData ( L1FileHandle, QuantityName, L1bData, NoMAFs, Flag, &
    & FirstMAF, LastMAF, NEVERFAIL, hdfVersion, dontPad, L2AUX )
    
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
    integer :: myhdfVersion
    logical :: myDontPad
    type(l1bdata_t) :: L1BDATATMP
    logical :: isScalar
    logical, parameter :: DEEBug = .false.

    ! Executable code
    if ( present(hdfVersion) ) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L1BDEFAULT_HDFVERSION
    end if
    ! print * 'hdfVersion: ', hdfVersion
    myDontPad = .false.
    if ( present(firstMAF) ) myDontPad = .true.
    if ( present(dontPad) ) myDontPad = dontPad

    if ( myhdfVersion == HDFVERSION_4 ) then
      call ReadL1BData_hdf4 ( L1FileHandle, trim(QuantityName), L1bData, &
      & NoMAFs, Flag, FirstMAF, LastMAF, NEVERFAIL, L2AUX )
    else
      call ReadL1BData_hdf5 ( L1FileHandle, trim(QuantityName), L1bData, &
      & NoMAFs, Flag, FirstMAF, LastMAF, NEVERFAIL, L2AUX )
      !Unfortunately, hdf5-formatted l1b data have different shapes from hdf4
      ! E.g., for MAFStartTimeTAI we obtain the following
      !  hdfVERSION      shape
      !     4           1   1   5
      !     5           5   1   1
      if ( Flag == 0 ) then
        call Reshape_for_hdf4(L1bData)
      else
        ! print *, 'Warning: ', trim(QuantityName) // ' not found in l1b files'
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to find '//trim(QuantityName)//' in l1b files')
      end if
    end if
    isScalar = ( l1bData%noMAFs < 2 ) ! There may be a better way to decide
    if ( myDontPad .or. flag /= 0 .or. isScalar .or. &
      & (GlobalAttributes%FirstMAFCtr > GlobalAttributes%LastMAFCtr) ) return
    if ( .not. IsL1BGappy(l1bData) ) return
    ! Must pad l1bdata; so first copy to temp l1bData
    if ( associated(l1bdata%dpField) .and. DEEBug ) then
      print *, 'max(l1bdata) ', maxval(l1bdata%dpField(1,1,:))
      print *, 'min(l1bdata) ', minval(l1bdata%dpField(1,1,:))
    endif
    call cpL1bData(l1bdata, l1bdataTmp)
    if ( associated(l1bdatatmp%dpField) .and. DEEBug ) then
      print *, 'max(l1bdatatmp) ', maxval(l1bdatatmp%dpField(1,1,:))
      print *, 'min(l1bdatatmp) ', minval(l1bdatatmp%dpField(1,1,:))
    endif
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
      return
    endif
    call BadPadL1BData(l1bdataTmp, l1bData, GlobalAttributes%FirstMAFCtr, NoMAFs)
    call deallocatel1bData(l1bdataTmp)
    if ( .not. present(firstMAF) .and. .not. present(lastMAF) ) return
    ! Need to contract padded l1bdata to just those MAFs requested
    call cpL1bData(l1bdata, l1bdataTmp)
    call deallocatel1bData(l1bdata)
    if ( DEEBug ) print *, 'preparing to contract'
    call ContractL1BData(l1bdataTmp, l1bData, noMAFs, firstMAF, lastMAF)
    call deallocatel1bData(l1bdataTmp)
  end subroutine ReadL1BData

  !--------------------------------------------  ReadL1BData_hdf4  -----
  subroutine ReadL1BData_hdf4 ( L1FileHandle, QuantityName, L1bData, &
    & NoMAFs, Flag, FirstMAF, LastMAF, NEVERFAIL, L2AUX )
    
    ! Dummy arguments
    character(len=*), intent(in)   :: QUANTITYNAME ! Name of SD to read
    integer, intent(in)            :: L1FILEHANDLE ! From HDF
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

    integer :: ALLOC_ERR
    integer :: DATA_TYPE
    integer :: DIM_SIZES(MAX_VAR_DIMS)
    integer :: I
    logical :: MyNeverFail
    integer :: N_ATTRS
    integer :: NUMMAFS
    integer :: RANK
    integer :: SDS1_ID
    integer :: SDS2_ID
    integer :: SDS_INDEX
    integer :: STATUS

    integer, dimension(:), pointer :: EDGE
    integer, dimension(:), pointer :: START
    integer, dimension(:), pointer :: STRIDE

    real(r4), pointer, dimension(:,:,:) :: tmpR4Field
    logical :: isL2AUX

    ! Executable code
    call deallocateL1BData ( l1bData ) ! Avoid memory leaks

    nullify ( edge, start, stride, tmpR4Field )
    flag = 0
    MyNeverFail = .false.
    if ( present(NeverFail) ) MyNeverFail = NeverFail
    isL2AUX = .false.
    if ( present(l2AUX) ) isL2AUX = L2AUX

    ! Find data sets for counterMAF & quantity by name

    sds_index = sfn2index(L1FileHandle, 'counterMAF')
    if ( sds_index == -1 ) then
      if ( .not. JUSTLIKEL2AUX ) then
        flag = NOCOUNTERMAFINDX
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Failed to find index of counterMAF data set.')
      else
        sds1_id = SD_NO_COUNTERMAF
      end if
    else

      sds1_id = sfselect(L1FileHandle, sds_index)
      if ( sds1_id == -1 ) then
        flag = NOCOUNTERMAFID
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Failed to find identifier of counterMAF data set.')
      end if
    end if

    sds_index = sfn2index(L1FileHandle, quantityName)
    if ( sds_index == -1 ) then
      flag = NOQUANTITYINDEX
      if ( MyNeverFail ) return
      dummy = 'Failed to find index of quantity "' // trim(quantityName) // &
        & '" data set.'
      call MLSMessage ( MLSMSG_Error, ModuleName, dummy )
    end if

    sds2_id = sfselect(L1FileHandle, sds_index)
    if ( sds2_id == -1 ) then
      flag = NODATASETID
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to find identifier of data set matching the index.')
    end if

    ! Find rank (# of dimensions), dimension sizes of quantity data set
    status = sfginfo ( sds2_id, dummy, rank, dim_sizes, data_type, &
      n_attrs )

    if ( status == -1 ) then
      flag = NODATASETRANK
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName,&
      & 'Failed to find rank of data set.')
    end if

    ! allocate, based on above SD, dim info
    call allocate_test ( edge,   rank, 'edge',   ModuleName )
    call allocate_test ( start,  rank, 'start',  ModuleName )
    call allocate_test ( stride, rank, 'stride', ModuleName )

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
        if ( MyNeverFail ) return
        call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & input_err // 'firstMAF (bad chunkNo?)' )
      end if
      l1bData%firstMAF = firstMAF
    else
      l1bData%firstMAF = 0
    end if

    if ( present (lastMAF) ) then
      if ( lastMAF < l1bData%firstMAF ) then
        flag = LASTMAFNOTFOUND
        if ( MyNeverFail ) return
        call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & input_err // 'last' )
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
      & 'counterMAF', ModuleName )
    if ( sds1_id /= SD_NO_COUNTERMAF ) then
      status = sfrdata_f90(sds1_id,  (/ l1bData%firstMAF /) , (/1/), &
        & (/l1bData%noMAFs/), l1bData%counterMAF )
      if ( status == -1 ) then
        flag = CANTREADCOUNTERMAF
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead // 'counterMAF.' )
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
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_allocate // ' charField pointer.' )
      end if
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
      nullify ( l1bData%charField, l1bData%intField )
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
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_L1BRead // quantityName )
    else if ( status == -2 ) then
      flag = UNKNOWNDATATYPE
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unknown data type in readl1bData'   )
    end if

    ! Terminate access to the data sets
    
    if ( sds1_id /= SD_NO_COUNTERMAF ) then
      status = sfendacc(sds1_id)
      if ( status == -1 ) then
        flag = CANTENDCOUNTERMAF
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Failed to terminate access to data sets.' )
        flag = -1
      end if
    end if

    status = sfendacc(sds2_id)
    if ( status == -1 ) then
      flag = CANTENDQUANTITY
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to terminate access to data sets.' )
      flag = -1
    end if

    ! deallocate local variables
    call Deallocate_test ( edge,   'edge',   ModuleName )
    call Deallocate_test ( start,  'start',  ModuleName )
    call Deallocate_test ( stride, 'stride', ModuleName )

  end subroutine ReadL1BData_hdf4

  !--------------------------------------------  ReadL1BData_hdf5  -----
  subroutine ReadL1BData_hdf5 ( L1FileHandle, QuantityName, L1bData, NoMAFs, &
    & Flag, FirstMAF, LastMAF, NEVERFAIL, L2AUX )
    use HDF5, only: HSIZE_T
    use MLSHDF5, only: IsHDF5DSPresent, LoadFromHDF5DS, &
      & GetHDF5DSRank, GetHDF5DSDims, GetHDF5DSQType
! use MLSAuxData, only: MLSAuxData_T, Read_MLSAuxData, Deallocate_MLSAuxData
    
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

    ! Local Parameters
    character (len=*), parameter :: INPUT_ERR = 'Error in input argument '
    ! integer, parameter :: MAX_NOMAFS = 7000     ! Expect ~3500 in one day
    integer, parameter :: SD_NO_COUNTERMAF = -2
    integer, dimension(:), pointer :: cm_array

    ! Local Variables

    character (len=128) :: DUMMY        ! Dummy quantity name

    character(len=1) :: Char_rank
    integer :: I
    integer :: MAFoffset
    logical :: MyNeverFail
    integer :: NUMMAFS
    integer :: RANK
    integer :: CMRANK
    character(len=16) :: QTYPE
    integer :: SDS1_ID
    integer :: STATUS

    integer(kind=hSize_t), dimension(:), pointer :: DIMS
    integer(kind=hSize_t), dimension(:), pointer :: MAXDIMS
    integer(kind=hSize_t), dimension(:), pointer :: CMDIMS
    integer(kind=hSize_t), dimension(:), pointer :: CMMAXDIMS

    logical, parameter           :: DEEBUG = .false.
    logical :: isL2AUX

    ! Executable code
    if ( present(FirstMAF) ) then
      MAFoffset = max(0, FirstMAF)   ! Never let this be < 0
    endif
    call deallocateL1BData ( l1bData ) ! Avoid memory leaks
    sds1_id = 0

    flag = 0
    MyNeverFail = .false.
    if ( present(NeverFail) ) MyNeverFail = NeverFail
    isL2AUX = .false.
    if ( present(l2AUX) ) isL2AUX = L2AUX

    if ( .not. IsHDF5DSPresent(L1FileHandle, QuantityName) ) then
      flag = NOQUANTITYINDEX
      ! print *, 'Oops--' // trim(QuantityName) // ' not here'
      if ( MyNeverFail ) return
      dummy = 'Failed to find index of quantity "' // trim(quantityName) // &
        & '" data set.'
      call MLSMessage ( MLSMSG_Error, ModuleName, dummy )
    end if
    
    ! Find Qtype, rank and dimensions of QuantityName
    ! print*, ' Find Qtype, rank and dimensions of QuantityName ', trim(QuantityName)
    call GetHDF5DSRank(L1FileHandle, QuantityName, rank)
    l1bData%TrueRank = rank
    allocate ( dims(rank), maxDims(rank) )
    call GetHDF5DSDims(L1FileHandle, QuantityName, dims, maxDims)
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
        if ( MyNeverFail ) return
        call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & input_err // 'firstMAF (bad chunkNo?)' )
      end if
      l1bData%firstMAF = firstMAF
    else
      l1bData%firstMAF = 0
    end if

    if ( present (lastMAF) ) then
      if ( lastMAF < l1bData%firstMAF ) then
        flag = LASTMAFNOTFOUND
        if ( MyNeverFail ) return
        call dump( (/lastMAF, l1bData%firstMAF/), 'lastMAF, l1bData%firstMAF')
        call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & input_err // 'last' )
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
      & 'counterMAF', ModuleName )
    ! Find data sets for counterMAF & quantity by name

    if ( .not. IsHDF5DSPresent(L1FileHandle, '/counterMAF') ) then
      if ( DEEBUG) print *, 'no counterMAF array in file'
      if ( .not. JUSTLIKEL2AUX ) then
        flag = NOCOUNTERMAFINDX
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Failed to find index of counterMAF data set.')
      else
        if ( .not. isL2AUX ) call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Failed to find index of counterMAF data set.')
        sds1_id = SD_NO_COUNTERMAF
        ! Since we aren't reading these, just make them internally consistent
        do i = 1, l1bData%noMAFs
          l1bData%counterMAF(i) = l1bData%firstMAF + i - 1
        end do
      end if
    else
      if ( DEEBUG) print *, 'Getting counterMAF rank'
      call GetHDF5DSRank(L1FileHandle, '/counterMAF', cmrank)
      if ( DEEBUG) print *, cmrank
      allocate ( cmdims(cmrank), cmmaxDims(cmrank) )
      if ( DEEBUG) print *, 'getting counterMAF dims'
      call GetHDF5DSDims(L1FileHandle, '/counterMAF', cmdims, cmmaxDims)
      if ( DEEBUG) print *, cmdims, cmmaxDims
      ! allocate(countermaf_ptr(MAX_NOMAFS))
      ! countermaf_ptr = 0
      if ( present(FirstMAF) ) then
        if ( DEEBUG) print *, 'reading counterMAF ', MAFoffset, l1bData%noMAFs
        call LoadFromHDF5DS(L1FileHandle, '/counterMAF', l1bData%counterMaf, &
          & (/MAFoffset/), (/l1bData%noMAFs/) )
      elseif ( cmdims(1) /= l1bData%noMAFs ) then
        flag = CANTREADCOUNTERMAF
        if ( MyNeverFail ) then
          deallocate(dims, maxDims, cmdims, cmmaxdims)
          return
        endif
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
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Sorry--ReadL1BData_hdf5 says counterMaf sized differently from ' &
          & // trim(QuantityName) )
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
      ! if ( sds1_id == -1 ) then
      !   flag = NOCOUNTERMAFID
      !   deallocate(dims, maxdims)
      !   if ( MyNeverFail ) return
      !   call MLSMessage ( MLSMSG_Error, ModuleName, &
      !   & 'Failed to find identifier of counterMAF data set.')
      ! end if
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
    endif
    select case (trim(Qtype) // Char_rank)
    case ('real1')
      allocate( l1bData%DpField(l1bData%noMAFs, 1, 1),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1), &
          & (/MAFoffset/), (/l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1))  
      end if                                                                  
      l1bdata%data_type = 'double'
    case ('real2')
      allocate( l1bData%DpField(dims(1),l1bData%noMAFs, 1),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1), &
          & (/0,MAFoffset/), (/int(dims(1)),l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1))  
      end if                                                                  
      l1bdata%data_type = 'double'
    case ('real3')
      allocate( l1bData%DpField(dims(1),dims(2),l1bData%noMAFs),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField, &
          & (/0,0,MAFoffset/), (/int(dims(1)),int(dims(2)),l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField)  
      endif                                                                  
      l1bdata%data_type = 'double'
    case ('double1')
      allocate( l1bData%DpField(l1bData%noMAFs, 1, 1),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1), &
          & (/MAFoffset/), (/l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1))  
      end if                                                                  
      l1bdata%data_type = 'double'
    case ('double2')
      allocate( l1bData%DpField(dims(1),l1bData%noMAFs, 1),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1), &
          & (/0,MAFoffset/), (/int(dims(1)),l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1))  
      end if                                                                  
      l1bdata%data_type = 'double'
    case ('double3')    
      allocate( l1bData%DpField(dims(1),dims(2),l1bData%noMAFs),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField, &
          & (/0,0,MAFoffset/), (/int(dims(1)),int(dims(2)),l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField)  
      endif                                                                  
      l1bdata%data_type = 'double'
    case ('integer1')  
      allocate( l1bData%intField(l1bData%noMAFs, 1, 1),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%intField(:,1,1), &
          & (/MAFoffset/), (/l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%intField(:,1,1))  
      end if                                                                  
      ! call MLSMessage ( MLSMSG_Error, ModuleName, &
      ! & 'Sorry--LoadFromHDF5DS not yet written for type integer(:).')
      l1bdata%data_type = 'integer'
    case ('integer3')  
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--LoadFromHDF5DS not yet written for type integer(:,:,:).')
      l1bdata%data_type = 'integer'
    case ('character3') 
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--LoadFromHDF5DS not yet rewritten for type char(:,:,:).')
      l1bdata%data_type = 'integer'
    case default 
      l1bdata%data_type = 'unknown'
      flag = UNKNOWNDATATYPE
      deallocate(dims, maxDims)
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--ReadL1BData_hdf5 has encountered an unknown data type: ' &
        & // trim(Qtype) // Char_rank)
    end select
    if ( DEEBUG ) call dump(l1bData, 0)
    
    ! Avoid memory leaks
    deallocate(dims, maxDims)
  end subroutine ReadL1BData_hdf5

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
    type(l1bdata_t), intent(inout) :: L1BDATA ! Result
    ! Local variables
    character, dimension(:,:,:), pointer :: CharField => NULL()
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
    character(len=1), parameter  :: method = 'b'
    integer                      :: mord
    
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
      if ( status /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 could not allocate ' // &
        & trim(l1bData%data_type))
      CharField = reshape(source=l1bData%CharField, &
        & shape=new_shape, order=new_order)
      deallocate(l1bData%CharField, stat=status)
      if ( status /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 could not deallocate ' // &
        & trim(l1bData%data_type))
      l1bData%CharField => CharField
    case('double')
      if ( .not. associated(l1bData%DpField) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 found ' // trim(l1bData%data_type) // &
        & 'unassociated')
      old_shape = shape(l1bData%DpField)
      new_shape(3) = old_shape(l1bData%TrueRank)
      if ( l1bData%TrueRank == 2 ) new_shape(mord) = old_shape(1)
      allocate(DpField(new_shape(1), new_shape(2), new_shape(3)), stat=status)
      if ( status /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 could not allocate ' // &
        & trim(l1bData%data_type))
      DpField = reshape(source=l1bData%DpField, &
        & shape=new_shape, order=new_order)
      deallocate(l1bData%DpField, stat=status)
      if ( status /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 could not deallocate ' // &
        & trim(l1bData%data_type))
      l1bData%DpField => DpField
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
      if ( status /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 could not allocate ' // &
        & trim(l1bData%data_type))
      IntField = reshape(source=l1bData%IntField, &
        & shape=new_shape, order=new_order)
      deallocate(l1bData%IntField, stat=status)
      if ( status /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--reshape_for_hdf4 could not deallocate ' // &
        & trim(l1bData%data_type))
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

  ! ---------------------------------------------  cpField  -----
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
    ! Executable
    nn = 0
    if ( present(m) ) nn=m-1
    ! print *, 'Copying from to ', nIn, nOut
    if ( associated(l1BDataIn%charField) .and. &
      & associated(l1BDataOut%charField) ) then
        l1BDataOut%charField(:,:,nOut:nOut+nn) = &
        & l1BDataIn%charField(:,:,nIn:nIn+nn)
    endif
    if ( associated(l1BDataIn%intField) .and. &
      & associated(l1BDataOut%intField) ) then
        l1BDataOut%intField(:,:,nOut:nOut+nn) = &
        & l1BDataIn%intField(:,:,nIn:nIn+nn)
    endif
    if ( associated(l1BDataIn%dpField) .and. &
      & associated(l1BDataOut%dpField) ) then
        if ( nOut+nn > size(l1BDataOut%dpField, 3) ) then
          call output('size(l1BDataOut%dpField, 3) ', advance='no')
          call output(size(l1BDataOut%dpField, 3), advance='yes')
          call output('nOut+nn ', advance='no')
          call output(nOut+nn, advance='yes')
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Tried to copy past size of output field in cpField')
        elseif ( nIn+nn > size(l1BDataIn%dpField, 3) ) then
          call output('size(l1BDataIn%dpField, 3) ', advance='no')
          call output(size(l1BDataIn%dpField, 3), advance='yes')
          call output('nIn+nn ', advance='no')
          call output(nIn+nn, advance='yes')
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Tried to copy past size of input field in cpField')
        endif
        l1BDataOut%dpField(:,:,nOut:nOut+nn) = &
        & l1BDataIn%dpField(:,:,nIn:nIn+nn)
    ! print *, l1BDataIn%dpField(1,1,nIn:nIn+nn)
    endif
  end subroutine cpField

  ! ---------------------------------------------  zeroField  -----
  subroutine zeroField(l1bData, n, value, m, chValue)
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
    endif
    if ( associated(l1BData%intField)  ) then
        l1BData%intField(:,:,n:n+nn) = int(zero)
    endif
    if ( associated(l1BData%dpField)  ) then
        l1BData%dpField(:,:,n:n+nn) = zero
    endif
  end subroutine zeroField

  ! ---------------------------------------------  myminval  -----
  function myminval(ints) result(mymin)
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
  enddo
  end function myminval

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
        call print_source ( source_ref(lcf_where) )
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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
end module L1BData

! $Log$
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
