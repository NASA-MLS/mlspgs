! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L1BData

  ! Reading and interacting with Level 1B data (HDF4)

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use Dump_0, only: DUMP
  use Hdf, only: DFACC_READ, SFSTART, SFGINFO, SFN2INDEX, SFSELECT, &
    & SFRDATA_f90, &
    & SFRCDATA, SFENDACC, DFNT_CHAR8, DFNT_INT32, DFNT_FLOAT64, &
    & DFNT_FLOAT32
  ! use HDF5, only: HSIZE_T
  use Lexer_Core, only: PRINT_SOURCE
  use MLSCommon, only: R4, R8, L1BINFO_T, FILENAMELEN
  use MLSFiles, only: HDFVERSION_4, HDFVERSION_5, &
    & MLS_HDF_VERSION, MLS_IO_GEN_OPENF
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR, &
    & MLSMSG_WARNING, MLSMSG_L1BREAD, MLSMSG_WARNING, MLSMSG_DEALLOCATE
  use MLSStrings, only: CompressString, NumStringElements
  use MoreTree, only: Get_Field_ID
  use Output_M, only: Output
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
! === (end of api) ===

  private

  public :: Dump, L1BData_T, L1BRadSetup, L1BOASetup, DeallocateL1BData, &
    & FINDL1BDATA, NAME_LEN, PRECISIONSUFFIX, ReadL1BData, &
    & AssembleL1BQtyName

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private not_used_here
  !---------------------------------------------------------------------------

  interface DUMP
    module procedure DumpL1BData
  end interface

  ! Parameters
  integer, parameter :: NAME_LEN = 64  ! Max len of SDS array name
  ! suffix of sd precision; check against 'grep -i precision l1/OutputL1B.f90'
  character  (len=*), parameter :: PRECISIONSUFFIX = ' precision'
  real, parameter ::    UNDEFINED_VALUE = -999.99 ! Same as l2/templates, Fill
  ! If TRUE, treats l1brad and l2aux files the same
  ! Among other things, this means forgiving absence of counterMAF
  ! so we may read l2aux and l1brad files alike
  logical, parameter            :: JUSTLIKEL2AUX = .true.

  ! Assume l1b files w/o explicit hdfVersion field are this
  ! 4 corresponds to hdf4, 5 to hdf5 in L2GP, L2AUX, etc. 
  integer, parameter :: L1BDEFAULT_HDFVERSION = HDFVERSION_4

  ! Unless the following is true, a subgroup named 'tp' will be
  ! created under each instrument Module for hdf5 versions
  ! (An early idea we later repented of)
  logical, parameter :: DROPTPSUBGROUP = .true.

  ! This data type is used to store quantities from an L1B data file.

  type L1BData_T
    character (len=name_len) :: L1BName ! Name of field in file
    character (len=16) :: data_type     ! 'character', 'double', or 'integer'
    integer :: FirstMAF                 ! First major frame read
    integer :: NoMAFs                   ! # of MAFs read
    integer :: MaxMIFs                  ! Max # of MIFs/MAF in SD array
    integer :: NoAuxInds                ! # of auxilliary indices
    integer :: TrueRank                 ! # necessary indices (e.g., 1 for MAFs)

    integer, dimension(:), pointer :: CounterMAF => NULL() ! dimensioned (noMAFs)

    character, dimension(:,:,:), pointer :: CharField => NULL()
    real(r8),  dimension(:,:,:), pointer :: DpField => NULL()
    integer,   dimension(:,:,:), pointer :: IntField => NULL()
    ! all the above dimensioned (noAuxInds,maxMIFs,noMAFs)
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


  !-------------------------------------------  DeallocateL1BData  -----
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
    endif
    if ( hdfVersion == HDFVERSION_5 ) then
      head = '/'
      instr_tail = '/'
      tp_tail = '/'
    else
      head = ''
      instr_tail = '.'
      if ( present(InstrumentName) ) then
        if ( trim(InstrumentName) == 'sc' ) instr_tail = ''
      endif
      tp_tail = ''
    endif
    if ( present(InstrumentName) ) then
      QtyName = head // trim(InstrumentName) // instr_tail
    elseif ( is_a_signal .or. hdfVersion /= HDFVERSION_5 ) then
      QtyName = head
    else
      ! Need only to convert complete hdf4-name to hdf5-name
      ! This means we must parse hdf4-name fully, however
      QtyName = head
      ! 1st--is there an instrument prefixed to name?
      if ( name(1:2) == 'sc' ) then
        my_instrument = 'sc'
        the_rest = name(3:)
      elseif ( name(1:4) == 'GHz.' ) then
        my_instrument = 'GHz'
        the_rest = name(5:)
      elseif ( name(1:4) == 'THz.' ) then
        my_instrument = 'THz'
        the_rest = name(5:)
      else
        my_instrument = ' '
        the_rest = name
      endif
      ! 2nd--Does the_rest start with 'tp'?
      if ( my_instrument == ' ' ) then
        QtyName = '/' // trim(name)
      elseif ( the_rest(1:2) /= 'tp' ) then
        QtyName = '/' // trim(my_instrument) // &
          & '/' // trim(the_rest)
      elseif (DROPTPSUBGROUP) then
        QtyName = '/' // trim(my_instrument) // &
          & '/' // trim(the_rest(3:))
      else
        QtyName = '/' // trim(my_instrument) // &
          & '/' // 'tp' // &
          & '/' // trim(the_rest(3:))
      endif
      ! We're done
      if ( DEEBUG ) then
        print *, 'converted name: ', trim(QtyName)
      endif
      return
    endif
    if ( isTngtQty ) then
      if ( DROPTPSUBGROUP .and. hdfVersion == HDFVERSION_5) then
        QtyName = trim(QtyName)
      else
        QtyName = trim(QtyName) // 'tp' // tp_tail
      endif
    endif
    QtyName = trim(QtyName) // trim(name)
    if ( compress ) QtyName = CompressString(QtyName)
    if ( DEEBUG ) then                            
      print *, 'more converted name: ', trim(QtyName)  
    endif                                         
  end function AssembleL1BQtyName

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
    
    if ( myDetails < -1 ) return
    call output('L1B rad quantity Name = ', advance='no')
    call output(L1bData%L1BName, advance='yes')
    call output('  First major frame read = ', advance='no')
    call output(L1bData%FirstMAF, advance='yes')
    call output('  Num of MAFs read = ', advance='no')
    call output(L1bData%NoMAFs, advance='yes')
    call output('  Max # of MIFs/MAF in SD array = ', advance='no')
    call output(L1bData%MaxMIFs, advance='yes')
    call output('  Num of auxilliary indices = ', advance='no')
    call output(L1bData%NoAuxInds, advance='yes')

    if ( myDetails < 0 ) return
    if ( associated(l1bData%counterMAF) ) then
      call dump ( l1bData%counterMAF,&
      & 'l1bData%counterMAF' )
    else
      call output('(CounterMAF array unassociated)', advance='yes')
    endif

    if ( myDetails < 1 ) return
    if ( associated(l1bData%charField) ) then
      call dump ( l1bData%CharField,&
      & 'l1bData%CharField' )
    else
      call output('(CharField array unassociated)', advance='yes')
    endif

    if ( associated(l1bData%intField) ) then
      call dump ( l1bData%intField,&
      & 'l1bData%intField' )
    else
      call output('(intField array unassociated)', advance='yes')
    endif

    if ( associated(l1bData%dpField) ) then
      call dump ( l1bData%dpField,&
      & 'l1bData%dpField' )
    else
      call output('(dpField array unassociated)', advance='yes')
    endif
  end subroutine DumpL1BData

  ! -------------------------------------------------  FindL1BData  ----
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
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L1BDEFAULT_HDFVERSION
    endif

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
      if(get_field_id(son) == f_file) then
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        ! sd_id = sfstart(Filename, DFACC_READ)
        if ( present(hdfVersion) ) then
          the_hdf_version = mls_hdf_version(FileName)
          sd_id = mls_io_gen_openF('hg', .true., error, &
            & record_length, DFACC_READ, &
            & FileName, hdfVersion=the_hdf_version, debugOption=.false.)
        else
          sd_id = mls_io_gen_openF('hg', .true., error, &
            & record_length, DFACC_READ, &
            & FileName)
        endif
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
      endif
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
      if(get_field_id(son) == f_file) then
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        if(.NOT. associated(l1bInfo%L1BRADIDs)) then
          allocate ( l1bInfo%L1BRADIDs(MAXNUML1BRADIDS), stat=status )
          allocate ( l1bInfo%L1BRADFileNames(MAXNUML1BRADIDS), stat=status )
          l1bInfo%L1BRADIDs = ILLEGALL1BRADID
          if ( status /= 0 ) &
            & call announce_error ( son, 'Allocation failed for l1bInfo' )
        endif
        ! sd_id = sfstart(Filename, DFACC_READ)
        if ( present(hdfVersion) ) then
          the_hdf_version = mls_hdf_version(FileName)
          sd_id = mls_io_gen_openF('hg', .true., error, &
            & record_length, DFACC_READ, &
            & FileName, hdfVersion=the_hdf_version, debugOption=.false.)
        else
          sd_id = mls_io_gen_openF('hg', .true., error, &
            & record_length, DFACC_READ, &
            & FileName)
        endif
        if ( sd_id <= 0 ) then
          call announce_error ( son, &
            & 'Error opening L1BRAD file: ' //Filename)
        elseif(ifl1 == MAXNUML1BRADIDS) then
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
      endif
    end do

  end subroutine L1bradSetup

  !-------------------------------------------------  ReadL1BData  -----
  subroutine ReadL1BData ( L1FileHandle, QuantityName, L1bData, NoMAFs, Flag, &
    & FirstMAF, LastMAF, NEVERFAIL, hdfVersion )
    
    ! Dummy arguments
    character(len=*), intent(in)   :: QUANTITYNAME ! Name of SD to read
    integer, intent(in)            :: L1FILEHANDLE ! From HDF
    integer, intent(in), optional  :: FIRSTMAF ! First to read (default 0)
    integer, intent(in), optional  :: LASTMAF ! Last to read (default last/file)
    logical, intent(in), optional  :: NEVERFAIL ! Don't call MLSMessage if TRUE
    type(l1bdata_t), intent(inout) :: L1BDATA ! Result
    integer, intent(out) :: FLAG        ! Error flag
    integer, intent(out) :: NOMAFS      ! Number actually read
    integer, optional, intent(in) :: hdfVersion

    ! Local variables
    integer :: myhdfVersion

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L1BDEFAULT_HDFVERSION
    endif
    ! print * 'hdfVersion: ', hdfVersion

    if ( myhdfVersion == HDFVERSION_4 ) then
      call ReadL1BData_hdf4 ( L1FileHandle, trim(QuantityName), L1bData, &
      & NoMAFs, Flag, FirstMAF, LastMAF, NEVERFAIL )
    else
      call ReadL1BData_hdf5 ( L1FileHandle, trim(QuantityName), L1bData, &
      & NoMAFs, Flag, FirstMAF, LastMAF, NEVERFAIL )
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
      endif
    endif
  end subroutine ReadL1BData

  !-------------------------------------------------  ReadL1BData_hdf4  -----
  subroutine ReadL1BData_hdf4 ( L1FileHandle, QuantityName, L1bData, &
    & NoMAFs, Flag, FirstMAF, LastMAF, NEVERFAIL )
    
    ! Dummy arguments
    character(len=*), intent(in)   :: QUANTITYNAME ! Name of SD to read
    integer, intent(in)            :: L1FILEHANDLE ! From HDF
    integer, intent(in), optional  :: FIRSTMAF ! First to read (default 0)
    integer, intent(in), optional  :: LASTMAF ! Last to read (default last in file)
    logical, intent(in), optional  :: NEVERFAIL ! Don't call MLSMessage if TRUE
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

    ! Executable code
    call deallocateL1BData ( l1bData ) ! Avoid memory leaks

    nullify ( edge, start, stride, tmpR4Field )
    flag = 0
    MyNeverFail = .false.
    if ( present(NeverFail) ) MyNeverFail = NeverFail

    ! Find data sets for counterMAF & quantity by name

    sds_index = sfn2index(L1FileHandle, 'counterMAF')
    if ( sds_index == -1) then
      if ( .not. JUSTLIKEL2AUX ) then
        flag = NOCOUNTERMAFINDX
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Failed to find index of counterMAF data set.')
      else
        sds1_id = SD_NO_COUNTERMAF
      endif
    else

      sds1_id = sfselect(L1FileHandle, sds_index)
      if ( sds1_id == -1) then
        flag = NOCOUNTERMAFID
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Failed to find identifier of counterMAF data set.')
      endif
    endif

    sds_index = sfn2index(L1FileHandle, quantityName)
    if ( sds_index == -1) then
      flag = NOQUANTITYINDEX
      if ( MyNeverFail ) return
      dummy = 'Failed to find index of quantity "' // trim(quantityName) // &
        & '" data set.'
      call MLSMessage ( MLSMSG_Error, ModuleName, dummy )
    end if

    sds2_id = sfselect(L1FileHandle, sds_index)
    if ( sds2_id == -1) then
      flag = NODATASETID
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to find identifier of data set matching the index.')
    endif

    ! Find rank (# of dimensions), dimension sizes of quantity data set
    status = sfginfo ( sds2_id, dummy, rank, dim_sizes, data_type, &
      n_attrs )

    if ( status == -1) then
      flag = NODATASETRANK
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName,&
      & 'Failed to find rank of data set.')
    endif

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

    if ( rank < 2 ) then
      l1bData%maxMIFs = 1
    else if ( rank > 2 ) then
      l1bData%maxMIFs = dim_sizes(2)
    else
      l1bData%maxMIFs = dim_sizes(1)
    end if

    if ( rank > 2 ) then
      l1bData%noAuxInds = dim_sizes(1)
    else
      l1bData%noAuxInds = 1
    end if

    ! Check input arguments, set noMAFs

    numMAFs = dim_sizes(rank)

    if ( present ( firstMAF ) ) then
      if ( (firstMAF >= numMAFs) .or. (firstMAF < 0) ) then
        flag = FIRSTMAFNOTFOUND
        if ( MyNeverFail ) return
        call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & input_err // 'firstMAF (bad chunkNo?)' )
      endif
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
      endif
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
      endif
    else
      ! Since we aren't reading these, just make them internally consistent
      do i = 1, l1bData%noMAFs
        l1bData%counterMAF(i) = l1bData%firstMAF + i - 1
      enddo
    endif

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
      endif
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
    elseif ( status == -2 ) then
      flag = UNKNOWNDATATYPE
      if ( MyNeverFail ) return
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unknown data type in readl1bData'   )
    endif

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

  !-------------------------------------------------  ReadL1BData_hdf5  -----
  subroutine ReadL1BData_hdf5 ( L1FileHandle, QuantityName, L1bData, NoMAFs, &
    & Flag, FirstMAF, LastMAF, NEVERFAIL )
  use HDF5, only: HSIZE_T
  use MLSHDF5, only: IsHDF5DSPresent, LoadFromHDF5DS, &
    & GetHDF5DSRank, GetHDF5DSDims, GetHDF5DSQType
! use MLSAuxData, only: MLSAuxData_T, Read_MLSAuxData, Deallocate_MLSAuxData
    
    ! Dummy arguments
    character(len=*), intent(in)   :: QUANTITYNAME ! Name of SD to read
    integer, intent(in)            :: L1FILEHANDLE ! From HDF
    integer, intent(in), optional  :: FIRSTMAF ! First to read (default 0)
    integer, intent(in), optional  :: LASTMAF ! Last to read (default last/file)
    logical, intent(in), optional  :: NEVERFAIL ! Don't call MLSMessage if TRUE
    type(l1bdata_t), intent(inout) :: L1BDATA ! Result
    integer, intent(out) :: FLAG        ! Error flag
    integer, intent(out) :: NOMAFS      ! Number actually read

    ! Local Parameters
    character (len=*), parameter :: INPUT_ERR = 'Error in input argument '
    integer, parameter :: MAX_VAR_DIMS = 32
    integer, parameter :: MAX_NOMAFS = 7000     ! Expect ~3500 in one day
    integer, parameter :: SD_NO_COUNTERMAF = -2

    ! Local Variables

    character (len=128) :: DUMMY        ! Dummy quantity name

    integer :: ALLOC_ERR
    character(len=1) :: Char_rank
    integer :: DATA_TYPE
    integer :: I
! > >      type(MLSAuxData_T) :: MLSAuxData
    logical :: MyNeverFail
    integer :: N_ATTRS
    integer :: NUMMAFS
    integer :: RANK
    character(len=16) :: QTYPE
    integer :: SDS1_ID
    integer :: SDS2_ID
    integer :: SDS_INDEX
    integer :: STATUS

    ! integer, dimension(:), pointer :: COUNTERMAF_PTR
    integer(kind=hSize_t), dimension(:), pointer :: DIMS
    integer(kind=hSize_t), dimension(:), pointer :: MAXDIMS

    real(r4), pointer, dimension(:,:,:) :: tmpR4Field
    logical, parameter           :: DEEBUG = .FALSE.

    ! Executable code
    call deallocateL1BData ( l1bData ) ! Avoid memory leaks

    flag = 0
    MyNeverFail = .false.
    if ( present(NeverFail) ) MyNeverFail = NeverFail

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
    allocate(dims(rank), maxDims(rank))
    call GetHDF5DSDims(L1FileHandle, QuantityName, dims, maxDims)
    call GetHDF5DSQType ( L1FileHandle, QuantityName, Qtype )
    if ( DEEBUG) print *, 'rank ', rank
    if ( DEEBUG) print *, 'maxDims ', maxDims
    if ( DEEBUG) print *, 'dims ', dims
    if ( DEEBUG) print *, 'Qtype ', Qtype
! > >     nullify(MLSAuxData%RealField, MLSAuxData%DpField, &
! > >       & MLSAuxData%IntField, MLSAuxData%CharField)
! > >     select case (trim(Qtype))
! > >     case ('real')
! > >       allocate( MLSAuxData%RealField(dims(1),dims(2),dims(3)),stat=status)
! > >       MLSAuxData%RealField = UNDEFINED_VALUE
! > >     case ('double')    
! > >       allocate( MLSAuxData%DpField(dims(1),dims(2),dims(3)),stat=status)
! > >       MLSAuxData%DpField = UNDEFINED_VALUE
! > >     case ('integer')  
! > >       allocate( MLSAuxData%IntField(dims(1),dims(2),dims(3)),stat=status)
! > >       MLSAuxData%IntField = int(UNDEFINED_VALUE)
! > >     case ('character') 
! > >       allocate( MLSAuxData%CharField(dims(1),dims(2),dims(3)),stat=status)
! > >       MLSAuxData%CharField = '(undefined)'
! > >     end select

    ! Check input arguments, set noMAFs

    numMAFs = dims(rank)

    if ( present ( firstMAF ) ) then
      if ( (firstMAF >= numMAFs) .or. (firstMAF < 0) ) then
        flag = FIRSTMAFNOTFOUND
        if ( MyNeverFail ) return
        call MLSMEssage ( MLSMSG_Error, ModuleName, &
        & input_err // 'firstMAF (bad chunkNo?)' )
      endif
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
      endif
      if ( lastMAF >= numMAFs ) then
        l1bData%noMAFs = numMAFs - l1bData%firstMAF
      else
        l1bData%noMAFs = lastMAF - l1bData%firstMAF + 1
      end if
    else
      l1bData%noMAFs = numMAFs - l1bData%firstMAF
    end if

    if ( DEEBUG) print *, 'l1bData%noMAFs ', l1bData%noMAFs
    l1bData%L1BName = quantityName

    noMAFs = l1bData%noMAFs
    if ( DEEBUG)  print *, 'noMAFs ', noMAFs

    if ( rank < 2 ) then
      l1bData%maxMIFs = 1
    else if ( rank > 2 ) then
      l1bData%maxMIFs = dims(2)
    else
      l1bData%maxMIFs = dims(1)
    end if
    if ( DEEBUG) print *, 'l1bData%maxMIFs ', l1bData%maxMIFs

    if ( rank > 2 ) then
      l1bData%noAuxInds = dims(1)
    else
      l1bData%noAuxInds = 1
    end if


    call Allocate_test ( l1bData%counterMaf, l1bData%noMAFs, &
      & 'counterMAF', ModuleName )
    ! Find data sets for counterMAF & quantity by name

    if ( .not. IsHDF5DSPresent(L1FileHandle, 'CounterMAF') ) then
      if ( .not. JUSTLIKEL2AUX ) then
        flag = NOCOUNTERMAFINDX
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Failed to find index of counterMAF data set.')
      else
        sds1_id = SD_NO_COUNTERMAF
        ! Since we aren't reading these, just make them internally consistent
        do i = 1, l1bData%noMAFs
          l1bData%counterMAF(i) = l1bData%firstMAF + i - 1
        enddo
      endif
    else

      ! allocate(countermaf_ptr(MAX_NOMAFS))
      ! countermaf_ptr = 0
      if ( present(FirstMAF) ) then
        call LoadFromHDF5DS(L1FileHandle, 'CounterMAF', l1bData%counterMaf, &
          & (/FirstMAF-1/), (/l1bData%noMAFs/) )
          ! & (/0/), (/l1bData%noMAFs/) )
      else
        call LoadFromHDF5DS(L1FileHandle, 'CounterMAF', l1bData%counterMaf)
      endif
      if ( sds1_id == -1) then
        flag = NOCOUNTERMAFID
        if ( MyNeverFail ) return
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Failed to find identifier of counterMAF data set.')
      endif
    endif

!    print *, 'About to use Read_MLSAuxData to read ', trim(QuantityName)
!    call Read_MLSAuxData(L1FileHandle, QuantityName, 'unknown', & 
!      & MLSAuxData, error, FirstMAF, LastMAF, read_attributes=.false.)
    ! print *, 'About to use LoadFromHDF5DS to read ', trim(QuantityName)
    ! print *, 'dims ', dims
    ! print * 'noMAFs ', l1bData%noMAFs
    write(Char_rank, '(i1)') rank
    select case (trim(Qtype) // Char_rank)
    case ('real1')
      allocate( l1bData%DpField(l1bData%noMAFs, 1, 1),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1), &
          & (/FirstMAF-1/), (/l1bData%noMAFs/) )                                      
          ! & (/0/), (/l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1))  
      endif                                                                  
      l1bdata%data_type = 'double'
    case ('real2')
      allocate( l1bData%DpField(dims(1),l1bData%noMAFs, 1),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1), &
          & (/0,FirstMAF-1/), (/int(dims(1)),l1bData%noMAFs/) )                                      
          ! & (/0,0/), (/int(dims(1)),l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1))  
      endif                                                                  
      l1bdata%data_type = 'double'
    case ('real3')
      allocate( l1bData%DpField(dims(1),dims(2),l1bData%noMAFs),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField, &
          & (/0,0,FirstMAF-1/), (/int(dims(1)),int(dims(2)),l1bData%noMAFs/) )                                      
          ! & (/0,0,0/), (/int(dims(1)),int(dims(2)),l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField)  
      endif                                                                  
!      l1bData%DpField = MLSAuxData%RealField
      l1bdata%data_type = 'double'
    case ('double1')
      allocate( l1bData%DpField(l1bData%noMAFs, 1, 1),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1), &
          & (/FirstMAF-1/), (/l1bData%noMAFs/) )                                      
          ! & (/0/), (/l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,1,1))  
      endif                                                                  
      l1bdata%data_type = 'double'
    case ('double2')
      allocate( l1bData%DpField(dims(1),l1bData%noMAFs, 1),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1), &
          & (/0,FirstMAF-1/), (/int(dims(1)),l1bData%noMAFs/) )                                      
          ! & (/0,0/), (/int(dims(1)),l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField(:,:,1))  
      endif                                                                  
      l1bdata%data_type = 'double'
    case ('double3')    
      allocate( l1bData%DpField(dims(1),dims(2),l1bData%noMAFs),stat=status)
      if ( present(FirstMAF) ) then                                          
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField, &
          & (/0,0,FirstMAF-1/), (/int(dims(1)),int(dims(2)),l1bData%noMAFs/) )                                      
          ! & (/0,0,0/), (/int(dims(1)),int(dims(2)),l1bData%noMAFs/) )                                      
      else                                                                   
        call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%DpField)  
      endif                                                                  
!      l1bData%DpField = MLSAuxData%DpField
      l1bdata%data_type = 'double'
    case ('integer3')  
! > >       allocate( l1bdata%IntField(dims(1),dims(2),l1bData%noMAFs),stat=status)
! > >       if ( present(FirstMAF) ) then                                          
! > >         call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%IntField, &
! > >           & (/0,0,0/), (/int(dims(1)),int(dims(2)),l1bData%noMAFs/) )                                      
! > >       else                                                                   
! > >         call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%IntField)  
! > >       endif                                                                  
! > !      l1bdata%IntField = MLSAuxData%IntField
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--LoadFromHDF5DS not yet written for type integer(:,:,:).')
      l1bdata%data_type = 'integer'
    case ('character3') 
! > >       allocate ( l1bData%charField(dims(1),dims(2),dims(3)), STAT=alloc_err )
! > >       if ( present(FirstMAF) ) then                                          
! > >         call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%CharField, &
! > >           & (/0,0,0/), (/int(dims(1)),int(dims(2)),l1bData%noMAFs/) )                                      
! > >       else                                                                   
! > >         call LoadFromHDF5DS(L1FileHandle, QuantityName, l1bData%CharField)  
! > >       endif                                                                  
!      l1bData%charField = MLSAuxData%CharField
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--LoadFromHDF5DS not yet rewritten for type char(:,:,:).')
      l1bdata%data_type = 'integer'
    case default 
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Sorry--ReadL1BData_hdf5 has encountered an unknown data type: ' &
        & // trim(Qtype) // Char_rank)
      l1bdata%data_type = 'unknown'
    end select
    ! print * 'datatype ', l1bdata%data_type
    if ( DEEBUG ) call dump(l1bData, 0)
    
    ! Avoid memory leaks
    deallocate(dims, maxDims)
! > >      call Deallocate_MLSAuxData(MLSAuxData)
  end subroutine ReadL1BData_hdf5

  ! ---------------------------------------------  reshape_for_hdf4  -----
  subroutine reshape_for_hdf4(L1bData)

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
    endif
    if ( l1bData%TrueRank > 2 .or. l1bData%TrueRank < 1 ) return
    new_shape = 1
    if ( method == 'a' .or. l1bData%TrueRank == 1 ) then
      new_order = reorder(:,l1bData%TrueRank)
      mord = 1
    else
      new_order = reorder(:,l1bData%TrueRank+1)
      mord = 2
    endif
    if ( DEEBUG ) then
      print *, 'new_order: ', new_order
      print *, 'mord: ', mord
      print *, '     o l d   f o r m'
      call DumpL1BData(l1bData)
    endif
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
      endif
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
    endif
  end subroutine reshape_for_hdf4

  ! ---------------------------------------------  announce_error  -----
  subroutine announce_error ( lcf_where, full_message, use_toolkit, &
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
  end subroutine announce_error

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
end module L1BData

! $Log$
! Revision 2.37  2003/03/07 00:35:15  pwagner
! Fixed bad bug in reading hdf5 files with FirstMAF > 1
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
