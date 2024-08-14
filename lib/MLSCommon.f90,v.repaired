! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module MLSCommon                ! Common definitions for the MLS software
!=============================================================================

  use IEEE_Arithmetic, only: IEEE_Is_Finite, IEEE_Is_Nan
  !   This doesn't result in a circular dependence.
  use Lexer_Types, only: Where_T ! Where is something in the L2CF
  use MLSKinds ! everything
  use MLSStrings_0,  only: Lowercase

  implicit none
  private

! Keep common parameters, flags, datatypes, and procedures here so
! that they may be shared by many higher-level modules
! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (data types and parameters)

! i1, i2, i4       integer types
! r4, r8           floating point types
! ip, rp           integer, floating point types used in forward model
! rv               floating point type used in vector quantity values
! rm               floating point type used in matrix values
! rt               floating point type used in topo-set values
! LineLen          character-length of most input
! FileNameLen      character-length of path/filenames
! BareFNLen        character-length of filenames
! namelen          Max length of hdf sds array name
! fill_signal      signal to is_what_ieee to check for undefined (fill) values
! finite_signal    signal to is_what_ieee to check for finite
! inf_signal       signal to is_what_ieee to check for inf
! nan_signal       signal to is_what_ieee to check for NaN
! ShortNameLen     character-length of short names (e.g., 'H2O')
! DefaultUndefinedValue 
!                  default fill values, e.g. when creating hdf arrays
! UndefinedValue    
!                  value to check for by is_what_ieee
! HDF_Acc_Create   signal to create the hdf file
! HDF_Acc_Rdonly   signal to open the hdf file for reading only
! HDF_Acc_RdWr     signal to open the hdf file for both reading and writing
!            Bits of MASK field as used in VectorValue_T
! M_Cloud = 2**4
! M_Fill = 2**2
! M_FullDerivatives = 2**1
! M_Ignore = 2**5
! M_LinAlg = 2**0      ! Don't use in linear algebra
! M_Spare = 2**6
! M_Tikhonov = 2**3    ! Where to do Tikhonov regularization
!
!     (severity levels)
! MLS_S_Success            if status not this, then something went wrong
! MLSMSG_Success           status returned when all went well
! MLSMSG_Pause             pause execution waiting for user input
! MLSMSG_Debug             should print only if debugging turned on
! MLSMSG_Info              fyi only
! MLSMSG_Testwarning       test to see if we would print this warning
! MLSMSG_Warning           not fatal, but deserving of attention
! MLSMSG_Error             quits after printing
! MLSMSG_Crash             should give traceback before quitting
!            Derived Types
! FileIDs_T        id numbers for file, group, swath or dataset
! L1BInfo_T        L1B data file names, etc. 
!                   (Should we replace these with FileIDs?)
! L2Metadata_T     Coords of (lon,lat) box to write as metadata
! MLSChunk_T       Chunk of level 1 data for level 2 to process independently
! MLSFile_T        File name, type, id, etc.
! MLSFill_T        Fill value data type
! MLSFills         database of MLSFill values
! Range_T          Bottom, Top of PCFID range
! Interval_T       Bottom, Top of a real interval (open)
! TAI93_Range_T    start, end times in TAI93 formatted r8
!       State-defining flags
! MLSVerbose       Should we print extra stuff?
! MLSDebug         Should we print even more stuff?
! MLSVerboseSticky
!                  Retain value of MLSVerbose thoughout phase
! MLSDebugSticky
!                  Retain value of MLSDebug thoughout phase

!     (subroutines and functions)
! accessType       Converts DFACC_* <-> {l_create, l_rdwr, l_rdonly}
! DontCrashHere    May we skip otherwise obligatory crash in named modules?
! inRange          Does an argument lie within a specified range or interval
! is_what_ieee     Is an argument a specified ieee type
! split_name_extension  
!                  Splits the input name.ext into name and exxt
! split_path_name  Splits the input path/name into path and name
! === (end of toc) ===                                                   
! === (start of api) ===
! char* accessType ( int dfacc )
! int accessType ( char* str )
! log DontCrashHere( char* arge )
! log inRange( int arg, Range_T range )
! log inRange( real arg, Interval_T range )
! log is_what_ieee( type_signal what, num arg )
! split_name_extension ( char* full_file_name, char* name, char* extension,
!  [char dot] )
! split_path_name ( char* full_file_name, char* path, char* name, [char slash] )
! === (end of api) ===

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains simple definitions that are common to all the MLS PGS
  ! f90 software.
  
  ! We also put in these functions.
  ! Should they be relocated to other modules? Where?

  public :: AccessType
  public :: DontCrashHere
  public :: InRange
  public :: Is_what_ieee
  public :: Split_name_extension
  public :: Split_path_name
  
  ! User-defined datatypes
  public :: FileIDs_T
  public :: MLSChunk_T
  public :: MLSFile_T
  public :: MLSFill_T
  public :: L1BInfo_T
  public :: Range_T
  public :: TAI93_Range_T
  public :: L2Metadata_T
  public :: Interval_T

  ! Make parameters gotten from MLSKinds public

  public :: i1
  public :: i2
  public :: i4
  public :: r4
  public :: r8
  public :: rm
  public :: rp
  public :: ip
  public :: rt
  public :: rv
  
  ! The following are integer flags 'signalling' what kind of ieee number
  ! we would want to test for, e.g. in a multi-purpose procedure
  ! Their principal use is as args to is_what_ieee
  ! Note that the same test value can be finite and a Fill value
  ! simultaneously
  integer, public, parameter :: Finite_Signal = 0
  integer, public, parameter :: Inf_Signal    = Finite_Signal + 1
  integer, public, parameter :: Nan_Signal    = Inf_Signal + 1
  integer, public, parameter :: Fill_Signal   = Nan_Signal + 1
  integer, public, parameter :: Oldfill_Signal= Fill_Signal + 1

  ! The following are integer flags that should be consistent
  ! with hdf settings; e.g. /software/toolkit/ifc17/hdf/include/hdf.inc
  ! You may overwrite them if you wish (and have a good reason to do so)
  integer, public, parameter :: HDF_Acc_Create = 4
  integer, public, parameter :: HDF_Acc_Rdonly = 1
  integer, public, parameter :: HDF_Acc_RdWr   = 3

  ! Bit of MASK field of VectorValue_T
  integer, public, parameter :: M_Cloud = 2**4
  integer, public, parameter :: M_Fill = 2**2
  integer, public, parameter :: M_FullDerivatives = 2**1
  integer, public, parameter :: M_Ignore = 2**5
  integer, public, parameter :: M_LinAlg = 2**0      ! Don't use in linear algebra
  integer, public, parameter :: M_Spare = 2**6
  integer, public, parameter :: M_Tikhonov = 2**3    ! Where to do Tikhonov regularization

  ! Define some low level parameters.  These are used by the calling code to
  ! indicate the severity or otherwise of the messages.
  ! Normally, we treat any severity of Error or worse as reason to stop.
  ! Any Warning is worth recording, and suppressed when too numerous.
  ! Info may be customized to show, phase name, chunk number, etc.
  ! Be advised, Crash may not properly close files opened by your run.
  ! Use it only for specific debugging where you need a walkback.
  ! See also MLSMessageConfig%crashOnAnyError

  integer, public, parameter :: MLS_S_Success = 0
  integer, public, parameter :: MLSMSG_Success     = MLS_S_Success ! == 0
  integer, public, parameter :: MLSMSG_Pause       = MLSMSG_Success + 1
  integer, public, parameter :: MLSMSG_Debug       = MLSMSG_Pause + 1
  integer, public, parameter :: MLSMSG_Info        = MLSMSG_Debug + 1
  integer, public, parameter :: MLSMSG_TestWarning = MLSMSG_Info + 1
  ! The next 3 should always be the highest, i.e. most sevre
  integer, public, parameter :: MLSMSG_Warning     = MLSMSG_TestWarning + 1
  integer, public, parameter :: MLSMSG_Error       = MLSMSG_Warning + 1
  integer, public, parameter :: MLSMSG_Crash       = MLSMSG_Error + 1

  ! Unless you fill the string table with l_ quantities from intrinsic
  ! make sure the next entry is .false. 
  ! (or else it will segment fault in get_string)
  ! MLSL2 does fill the string table, but many tools do not
  logical, save, public :: FILESTRINGTABLE = .FALSE.

  interface accessType
    module procedure accessDFACCToStr
    module procedure accessStrToDFACC
  end interface

  ! Not just "is" but "=", because "is" may satisfy looser conditions
  interface equalsFillValue
    module procedure equalsFillValue_r4, equalsFillValue_r8, equalsFillValue_int
  end interface

  interface inRange
    module procedure inRange_r4, inRange_r8, inRange_int, inRange_range
  end interface

  interface is_a_fill_value
    module procedure is_a_fill_value_r4, is_a_fill_value_r8, is_a_fill_value_int
  end interface

  interface is_what_ieee
    module procedure is_what_ieee_r4, is_what_ieee_r8, is_what_ieee_integer
    module procedure is_what_ieee_character
  end interface
  
  !-----------------------------------------------------------
  !   P a r a m e t e r s   a n d    D a t a t y p e s
  !-----------------------------------------------------------
  logical, public :: MLSVerbose = .false. ! Should we print extra stuff?
  logical, public :: MLSDebug =   .false. ! Should we print even more stuff?
  logical, public :: MLSVerboseSticky = .false. ! Always, not just by module
  logical, public :: MLSDebugSticky =   .false. ! Always, not just by module
  integer, parameter :: DEBUGNAMESLEN = 256
  character(len=DEBUGNAMESLEN), public, save :: MLSNamesAreDebug   = ' '
  character(len=DEBUGNAMESLEN), public, save :: MLSNamesAreVerbose = ' '
  character(len=DEBUGNAMESLEN), public, save :: MLSNamesDontCrash = ' '
  ! Because we'd like not to always drag the SDPToolkit with us
  ! everywhere we go
  integer, private, parameter :: PGSd_PC_FILE_PATH_MAX = 1024

  ! Now we have the lengths for various strings

  integer, public, parameter :: ShortNameLen = 32
  integer, public, parameter :: NameLen      = 64  ! Max len of SDS array name
  integer, public, parameter :: LineLen      = 132
  integer, public, parameter :: FileNameLen  = max(PGSd_PC_FILE_PATH_MAX, 132) ! was 132
  integer, public, parameter :: BareFNLen    = 64  ! Bare file name length (w/o path)

  ! The next are used for tracking allocated memory
  ! (The 1st is public to enable reporting finer or coarser grains)
  !  real, save, public  :: MEMORY_UNITS = 1024. ! Report nothing smaller than KB
  integer, save, public          :: NoBlocksAllocated = 0 ! Num allocate calls
  integer, save, public          :: NoBlocksDeAllocated = 0 ! Num deallocate calls
  double precision, save, public :: NoBytesAllocated = 0.0d0 ! Net MEMORY_UNITS allocated.
  double precision, save, public :: TotalAllocated = 0.0d0 ! Total allocated.
  double precision, save, public :: TotalDeAllocated = 0.0d0 ! Total deallocated.

  !----------------------------------------------------------------------
  !         Undefined value
  ! This can be set to a different value if that would be more convenient
  ! Undefined value has two possible uses:
  ! (1) on output, it may be pre-assigned to every numerically-valued
  ! entry. If any persist after the algorithm completes, that would
  ! be a sign the algorithm failed to calculate a valid result for that
  ! entry; e.g., it may have encountered input data outside its range of validity

  ! (2) on input, entries may be checked against it for a signal that the data
  ! should be ignored or that it is a Fill value
  
  ! On the other hand, we may phase (2) out in favor of the MLSFills
  ! mechanism; see below
  real(r4), public, parameter :: DefaultUndefinedValue = -999.99 ! Try to use in lib, l2
  integer, public, parameter  :: UndefinedIntegerValue = -999 ! Use for int fields
  real(r4), public, parameter :: UndefinedTolerance    = 0.2 ! Poss. could make it 1

  real(r4), public, save      :: UndefinedValue        = DefaultUnDefinedValue
  ! --------------------------------------------------------------------------

  !----------------------------------------------------------------------
  !         What do we mean by the start array?
  ! The hdf library, and any libraries based on hdf, were written with
  ! the c language in mind. So they treat the first element of an array
  ! as index number "0" instead of "1". Therefore, start[k]=0
  ! means the first element. To a Fortran mentality, start functions not like
  ! an array of first indexes but like an array of offsets. We have then two
  ! choices of how to treat the start array when not calling an hdf procedure:
  ! (1) Fort_Indx: treat start as the first element
  ! (2) Offset:    adopt the hdf convention, and treat it as an offset
  ! The choice we make is set by the MLS_HyperStart parameter below.
  integer, parameter             :: Fort_Indx            = 1
  integer, parameter             :: Offset               = 0
  integer, parameter, public     :: MLS_HyperStart       = Offset
  ! When used in a procedure that interfaces directly with hdf, nothing
  ! changes. Just remember that start[k]=0 means the first element.
  !
  ! When used with an mls procedure out of Hyperslabs.f90, or when used with
  ! L2GPData/ExtractL2GPRecord, the following hold
  ! MLS_HyperStart                 start[k]=this 
  !                                means the first element
  ! ----------                 ---------------------------
  ! Fort_Indx                             1
  ! Offset                                0
  !
  ! Up to the time of this writing (2021-04-02) we have allowed
  ! MLS_HyperStart == Fort_Index
  ! While intuitive for a Fortran user, it is inconsistent with the usage
  ! in the hdf library. Therefore, we may someday switch to 
  ! MLS_HyperStart == Offset
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! A type to hold the hdf file ids
  ! (Should we make it recursive in case dataset path something like
  ! "/grp_1/grp_2/../grp_n/sd"?)

  type FileIds_T
    integer :: f_id     = 0 ! File id, handle, or io unit
    integer :: grp_id   = 0 ! group id
    integer :: sd_id    = 0 ! sd or swath id
  end type Fileids_T

  ! A PCFid range
  ! Also could be a ranges of MAFs
  ! Should we put one inside the MLSChunk_T?
  ! Or should we create a type-limited function inside the MLSChunk_T?
  type Range_T
    integer :: Bottom   = 0
    integer :: Top      = 0
  end type Range_T

  ! An open real interval (But see inRange below)
  type Interval_T
    real(rt) :: Bottom   = 0._rt
    real(rt) :: Top      = 0._rt
  end type Interval_T

  ! --------------------------------------------------------------------------

  ! This datatype defines the `chunks' into which the input dataset is split

  ! Moved here from Chunks_m module
  type MLSChunk_T
    logical :: abandoned = .false. ! Did we abandon this chunk's retrieval?
    integer :: firstMAFIndex = -1  ! Index of first MAF in the chunk
    integer :: lastMAFIndex = -1   ! Index of last MAF in the chunk
    integer :: noMAFsLowerOverlap = 0 ! Number of MAFs in the lower overlap region
    integer :: noMAFsUpperOverlap = 0 ! Number of MAFs in the upper overlap region
    integer :: chunkNumber        = -1             ! Index of this chunk
    integer, dimension(:), pointer :: HGridOffsets => NULL()
    ! This for each chunk is the index of the first non-overlapped profile in 
    ! each hGrid into the relevant output (l2gp?) file.
    integer, dimension(:), pointer :: HGridTotals => NULL()
    ! This is somewhat repetitive.  It's the total number of profiles in
    ! the output hGrid.  It's only really used in parallel runs.
    real(rp) :: phiStart = 0. ! for use by regular HGrid
    real(rp) :: phiEnd   = 0.
    double precision :: StartTime = 0. ! for use by readGriddedData
    double precision :: EndTime   = 0.
  end type MLSChunk_T

  ! Information describing the files used by the mls software
  ! Instead of passing file handles or names back & forth between routines
  ! -- pass one of these instead
  type MLSFile_T
    character (len=16) :: content=""  ! e.g., 'l1brad', 'l2gp', 'l2aux', ..
    character (len=8) :: lastOperation=""  ! 'open','close','read','write'
    character (len=FileNameLen) :: Name=""  ! its name (usu. w/path)
    character (len=ShortNameLen) :: ShortName=""  ! its short name; e.g. 'H2O'
    character (len=8) :: typeStr=""  ! one of {'swath', 'hdf', ..}
    integer :: type=0  ! one of {l_swath, l_hdf, ..}
    integer :: access=0  ! one of {DFACC_RDONLY, DFACC_CREATE, ..}
    integer :: HDFVersion=0  ! its hdf version if hdf(eos)
    integer :: PCFId=0      ! its PCF ID (ref), if any
    integer :: recordLength=0! its max record_length, if any
    integer :: errorCode=0  ! non-zero usu. means trouble
    logical :: StillOpen=.false.
    type(Range_T) :: PCFidRange = Range_T()
    type(Fileids_T) :: FileID = FileIDs_T()
    type(where_t) :: Where  ! in l2cf, using Where function in Tree module
    ! Until MLSFile_T is moved out of MLSCommon, get these from a tree node
    ! using Source_Ref and The_File from the Tree module.
    integer :: L2CF = 0       ! String index of L2CF
    integer :: Source = 0     ! 256*line + column, in L2CF
  end type MLSFile_T

  ! This datatype describes the information on the L1B data files in use
  type L1BInfo_T
    integer :: L1BOAId=0     ! The HDF ID (handle) for the L1BOA file
    ! Id(s) for the L1BRAD file(s)
    integer, dimension(:), pointer :: L1BRADIds=>NULL()
    character (len=FileNameLen) :: L1BOAFileName=""  ! L1BOA file name
    character (len=FileNameLen), dimension(:), pointer :: &
         & L1BRADFileNames=>NULL()
  end type L1BInfo_T

  ! --------------------------------------------------------------------------

  ! The TAI93 time range
  ! Must be r8 because otherwise seconds since 1993 will be truncated
  type TAI93_Range_T
    real(r8) :: startTime ! TAI93 format
    real(r8) :: endTime   ! TAI93 format
  end type TAI93_Range_T
  ! --------------------------------------------------------------------------

  ! Datatype for the metadata files supplemental data
  type L2Metadata_T
    real :: minLat = -90.0       ! min Latitude
    real :: maxLat = 90.0      ! max Latitude
    real :: minLon = -180.0      ! min Longitude
    real :: maxLon = 180.0     ! max Longitude
    character(len=NameLen) :: doiIdentifier = ' '
  end type L2Metadata_T

  ! --------------------------------------------------------------------------
  ! Datatype for the Fill or undefined value
  ! (its value(s), whether Fills
  ! mean a value greater than or equal to, etc.)
  ! We will be migrating older is-Fill tests to use this instead
  ! to handle multiple Fill values and idea of range
  
  ! The idea is that a tested value is a fill value if it satisfies any
  ! of a set of conditions with respect to reference values
  ! This is an obvious generalization of the older test where we
  ! tested only whether it equaled one specified reference value
  
  ! Why do this? We use -999.99, both in mls level 1 and level 2
  ! Gloria uses 10^12
  ! GMAO uses 10^15

  type MLSFill_T
    real(r8)             :: value
    ! The condition for each value; possibilities are (testing a value "it")
    ! id   condition(i)      meaning (i.e., it is a Fill value if ..
    !  0    =              it = values(i)
    !  1    >              it > values(i)
    !  2    <              it < values(i)
    !  3    >=             it > or = values(i)
    
    !  4    <=             it < or = values(i)
    ! In addition, if '||' are found in the condition(i), we
    ! take its absolute value before applying the test
    ! (We also add 10 to its id--YACH "yet another crude hack")
    real(r8)             :: tol ! its tolerance, if any
    character(len=4) :: condition
  end type MLSFill_T
  
  type(MLSFill_T), public, save, dimension(:), pointer :: MLSFills => null()

contains
  ! Now any public procedures and functions
  ! Should we keep procedures and ffunctions here, or
  ! keep break them out to a separate module?
  !--------------------------------------------  accessType  -----
  function accessDFACCToStr ( dfacc ) result(str)

  ! This routine converts an hdf access type
  ! like DFACC_RDONLY into a string like 'rdonly'
  ! If access type is unrecognized, returns 'unknown'
  ! Args
  integer, intent(in)           :: dfacc
  character(len=8)              :: str
  ! Executable
  select case (dfacc)
  case (HDF_Acc_Create)
    str = 'create'
  case (HDF_Acc_Rdonly)
    str = 'rdonly'
  case (HDF_Acc_RdWr)
    str = 'rdwrite'
  case default
    str = 'unknown' ! Why not ' '? Or '?'
  end select
  end function accessDFACCToStr

  function accessStrToDFACC ( str ) result(dfacc)

  ! This routine converts a string like 'rdonly' into an hdf access type
  ! like DFACC_RDONLY
  ! If string is unrecognized, returns -999
  ! Args
  character(len=*), intent(in) :: str
  integer                      :: dfacc
  ! Executable
  select case (lowercase(str(1:4)))
  case ('crea')
    dfacc = HDF_Acc_Create
  case ('rdon')
    dfacc = HDF_Acc_Rdonly
  case ('rdwr')
    dfacc = HDF_Acc_RdWr
  case default
    dfacc = -999 ! why not DFACC_RDWR?
  end select
  end function accessStrToDFACC

  !--------------------------------------------  DontCrashHere  -----
  ! May we skip otherwise obligatory crashes in named modules?
  ! We decide on the basis of whether the module
  logical function DontCrashHere( arg )
    character(len=8), intent(in) :: arg
    ! Executable
    DontCrashHere = .false.
    if ( len_trim(MLSNamesDontCrash ) < 1 ) return
    if ( len_trim(arg) < 1 ) return
    DontCrashHere = index( lowercase(MLSNamesDontCrash), lowercase(trim(arg)) ) > 0
  end function DontCrashHere

  !--------------------------------------------  InRange  -----
  ! This family of functions returns TRUE for arg(s) within the given range
  ! Like FindInRange, the range includes its endpoints;
  ! e.g. if arg = range%Bottom = range%Top, returns TRUE
  ! Note:
  ! This runs counter to the spirit that the interval_t is an open interval
  elemental function inRange_int(arg, range) result(relation)
    ! Is arg in range?
    integer, intent(in)       :: arg
    type(Range_T), intent(in) :: range
    logical                   :: relation
    relation = (arg < (range%top + 1)) .and. (arg > (range%bottom - 1))
  end function inRange_int

  elemental function inRange_range(arg, range, butNotEqual) result(relation)
    ! Is arg wholly contained in range?
    type(Range_T), intent(in)       :: arg
    type(Range_T), intent(in)       :: range
    logical, optional, intent(in)   :: butNotEqual
    logical                         :: relation
    logical                         :: myButNot
    myButNot = .false.
    if ( present(butNotEqual) ) myButNot = butNotEqual
    relation = (arg%top < (range%top + 1)) .and. (arg%bottom > (range%bottom - 1))
    if ( myButNot .and. relation ) &
      & relation = .not. ( arg%top == range%top .and. arg%Bottom == range%Bottom )
  end function inRange_range

  elemental function inRange_r4(arg, range) result(relation)
    ! Is arg in range?
    real(r4), intent(in)       :: arg
    type(Interval_T), intent(in) :: range
    logical                   :: relation
    relation = (arg <= range%top) .and. (arg >= range%bottom)
  end function inRange_r4

  elemental function inRange_r8(arg, range) result(relation)
    ! Is arg in range?
    real(r8), intent(in)       :: arg
    type(Interval_T), intent(in) :: range
    logical                   :: relation
    relation = (arg <= range%top) .and. (arg >= range%bottom)
  end function inRange_r8

  !--------------------------------------------  Is_What_IEEE  -----
  elemental function is_what_ieee_r8( what, arg ) result( itIs )
    ! Args
    integer, intent(in)                       :: what ! a signal flag
    real(r8), intent(in)                      :: arg
    logical                                   :: itIs
    ! Executable
    select case (what)
    case (finite_signal)
      itIs = ieee_is_finite(arg)
    case (inf_signal)
      itIs = .not. ( ieee_is_finite(arg) .or. ieee_is_nan(arg) )
    case (nan_signal)
      itIs = ieee_is_nan(arg)
    case (oldfill_signal)
      itIs = abs(arg - UndefinedValue) < &
        & max( Real(undefinedTolerance, r8), abs(UndefinedValue/100000._r8) )
    case (fill_signal)
      itIs = is_a_fill_value( arg )
    case default
      itIs = .false. ! what signal flag did you mean? not recognized
    end select
  end function is_what_ieee_r8

  elemental function is_what_ieee_r4( what, arg ) result( itIs )
    ! Args
    integer, intent(in)                       :: what ! a signal flag
    real(r4), intent(in)                      :: arg
    logical                                   :: itIs
    ! Executable
    select case (what)
    case (finite_signal)
      itIs = ieee_is_finite(arg)
    case (inf_signal)
      itIs = .not. ( ieee_is_finite(arg) .or. ieee_is_nan(arg) )
    case (nan_signal)
      itIs = ieee_is_nan(arg)
    case (oldfill_signal)
      itIs = ( abs(arg-int(UndefinedValue)) < &
        & max(undefinedTolerance, abs(UndefinedValue/100000)) )
    case (fill_signal)
      itIs = is_a_fill_value( arg )
    case default
      itIs = .false. ! what signal flag did you mean? not recognized
    end select
  end function is_what_ieee_r4

  elemental function is_what_ieee_integer( what, arg ) result( itIs )
    ! Args
    integer, intent(in)                       :: what ! a signal flag
    integer, intent(in)                       :: arg
    logical                                   :: itIs
    ! Executable
    select case (what)
    case (finite_signal)
      itIs = .true. ! ieee_is_finite(arg)
    case (inf_signal)
      itIs = .false. ! .not. ( ieee_is_finite(arg) .or. ieee_is_nan(arg) )
    case (nan_signal)
      itIs = .false. ! ieee_is_nan(arg)
    case (oldfill_signal)
      itIs = ( abs(arg-int(UndefinedValue)) < 1 )
    case (fill_signal)
      itIs = is_a_fill_value( arg )
    case default
      itIs = .false. ! what signal flag did you mean? not recognized
    end select
  end function is_what_ieee_integer

  elemental function is_what_ieee_character( what, arg ) result( itIs )
    ! Args
    integer, intent(in)                       :: what ! a signal flag
    character(len=*), intent(in)              :: arg
    logical                                   :: itIs
    ! Executable
    select case (what)
    case (finite_signal)
      itIs = .false. ! ieee_is_finite(arg)
    case (inf_signal)
      itIs = .false. ! .not. ( ieee_is_finite(arg) .or. ieee_is_nan(arg) )
    case (nan_signal)
      itIs = .false. ! ieee_is_nan(arg)
    case (fill_signal)
      itIs = .false.
    case default
      itIs = .false. ! what signal flag did you mean? not recognized
    end select
  end function is_what_ieee_character
  
  !--------------------------------------------  Is_A_Fill_Value  -----
  elemental function is_a_fill_value_int( arg ) result ( itIs )
    integer, intent(in) :: arg
    logical :: itIs
    integer :: k
    integer :: x
    ! We can't use the standard real-valued undefined values because
    ! we may wish them to have values (like -1.e15) so big
    ! that the int function will fail
    ! include 'isafillvalue.f9h'
    itIs = ( abs(arg-UndefinedIntegerValue) < 1 )
  end function is_a_fill_value_int
  
  elemental function is_a_fill_value_r4( arg ) result ( itIs )
    real(r4), intent(in) :: arg
    logical :: itIs
    integer :: k
    real(r4) :: x
    include 'isafillvalue.f9h'
  end function is_a_fill_value_r4
  
  elemental function is_a_fill_value_r8( arg ) result ( itIs )
    real(r8), intent(in) :: arg
    logical :: itIs
    integer :: k
    real(r8) :: x
    include 'isafillvalue.f9h'
  end function is_a_fill_value_r8
  
  !--------------------------------------------  WhatCondition  -----
  elemental function whatCondition ( c ) result ( what )
    ! Returns what condition id
    character(len=*), intent(in) :: c
    integer :: what
    if ( index(c, '<=') > 0 ) then
      what = 4
    elseif ( index(c, '>=') > 0 ) then
      what = 3
    elseif ( index(c, '<') > 0 ) then
      what = 2
    elseif ( index(c, '>') > 0 ) then
      what = 1
    else
      what = 0
    endif
    if ( index(c, '|') > 0 ) what = what + 10
  end function whatCondition
  
  !--------------------------------------------  equalsFillValue  -----
  elemental function equalsFillValue_int( arg, theFillValue ) result( itsafill )
    integer, intent(in) :: arg
    real(r8), intent(in) :: theFillValue
    logical :: itsafill
    itsafill = ( abs(arg-int(theFillValue)) < 1 )
  end function equalsFillValue_int

  elemental function equalsFillValue_r4( arg, theFillValue ) result( itsafill )
    real(r4), intent(in) :: arg
    real(r8), intent(in) :: theFillValue
    logical :: itsafill
    itsafill = ( abs(arg-theFillValue) < &
    & max(real(undefinedTolerance, r8), abs(theFillValue*1.d-6)) )
  end function equalsFillValue_r4

  elemental function equalsFillValue_r8( arg, theFillValue ) result( itsafill )
    real(r8), intent(in) :: arg
    real(r8), intent(in) :: theFillValue
    logical :: itsafill
    itsafill = ( abs(arg-theFillValue) < &
    & max(real(undefinedTolerance, r8), abs(theFillValue*1.d-6)) )
  end function equalsFillValue_r8

  ! --------------------------------------------  split_name_extension  -----

  ! This routine splits the input full_file_name
  ! E.g. 'name.txt' -> 'name' + 'txt'

  ! Optionally you may supply the '.'
  ! which must be a single character

  elemental subroutine split_name_extension ( full_file_name, name, extension, &
    & dot )

    ! Arguments

    character (len=*), intent(in) :: full_file_name
    character (len=*), intent(out) :: name
    character (len=*), intent(out) :: extension
    character (len=1), optional, intent(in) :: dot

    ! Local

    character (len=1) :: myDot
    integer :: loc, n
!   logical, parameter :: DEBUG = .false.

    ! Begin
    myDot = '.'
    if ( present(dot) ) myDot = dot
    call split_path_name ( full_file_name, name, extension, slash=myDot )
    ! Now 
    ! (1) Typically name will end with '.', so snip that off
    ! (2) If the file name had no extension, its name has been returned
    !     in the extension arg, so that must be dealt with
    ! (3) If the file name was e.g., '.ext', the name has been returned as '.'
    if ( len_trim(name) < 1 ) then
      ! case (2)
      name = extension
      extension = ' '
    elseif ( name == myDot ) then
      ! case (3)
      name = ' '
    else
      ! Now check for case (1)
      loc = len_trim(name)
      if ( name(loc:loc) == myDot ) name(loc:loc) = ' '
    endif
      
  end subroutine split_name_extension

  ! --------------------------------------------  split_path_name  -----

  ! This routine splits the input full_file_name
  ! into its components path and name
  ! where path may include one or more "/" or slash elements
  ! (but one must be the terminating one; e.g., 'System/')
  ! while name must have none (actually "/" comes from Machine%filsep).
  ! special cases by example: full_file_name -> (path, name)
  ! look.ma.no.slash -> (' ', 'look.ma.no.slash')
  ! Luke/I/am/your/father/ -> ('Luke/I/am/your/father/', ' ')

  ! optionally you may supply the slash divider
  ! which must be a single character

  elemental subroutine split_path_name ( full_file_name, path, name, slash )

    ! Arguments

    use Machine, only: Filsep ! / or :\

    character (len=*), intent(in) :: full_file_name
    character (len=*), intent(out) :: path
    character (len=*), intent(out) :: name
    character (len=1), optional, intent(in) :: slash

    ! Local

    character (len=1) :: mySlash
    integer :: loc, n
!   logical, parameter :: DEBUG = .false.

    ! Begin

    n = len_trim(full_file_name)

    if ( n <= 0 ) then
      path = ' '
      name = ' '
      return
    end if

    if ( present(slash) ) then
      mySlash = slash
    else
      mySlash = filsep
    end if

    loc = scan(full_file_name(:n), mySlash, back=.true.)

    if ( loc <= 0 ) then
      path = ' '
      name = adjustl(full_file_name)
    else if ( loc == n ) then
      path = adjustl(full_file_name)
      name = ' '
    else
      path = adjustl(full_file_name(:loc))
      name = adjustl(full_file_name(loc+1:))
    end if

  end subroutine split_path_name

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSCommon
!=============================================================================

!
! $Log$
! Revision 2.63  2024/08/14 22:51:54  pwagner
! Added optional arg to inRange_range
!
! Revision 2.62  2024/08/08 20:37:33  pwagner
! inRange now works with arg that is a range,  too
!
! Revision 2.61  2023/09/28 20:57:28  pwagner
! Added split_name_extension
!
! Revision 2.60  2021/04/29 22:51:07  pwagner
! Makes MLS_HyperStart conform with hdf
!
! Revision 2.59  2021/04/15 22:42:28  pwagner
! Added MLS_HyperStart; should it be here?
!
! Revision 2.58  2020/07/09 23:51:12  pwagner
! Added Start,EndTime components to MLSChunk_T
!
! Revision 2.57  2020/04/27 21:28:55  pwagner
! Separately track allocates/deallocates
!
! Revision 2.56  2019/08/19 21:56:39  pwagner
! Moved NoBytesAllocated to MLSCommon from Allocate_Deallocate
!
! Revision 2.55  2019/05/13 23:31:15  pwagner
! Introduced UndefinedIntegerValue
!
! Revision 2.54  2019/04/09 20:30:59  pwagner
! Moved some procedures from MLSStrings to new MLSStrings_0
!
! Revision 2.53  2018/12/11 16:46:30  pwagner
! Changed parameter name to MLS_S_Success to avoid conflict in level 1
!
! Revision 2.52  2018/12/11 01:19:01  pwagner
! moved MLSMSG_sevrity parameters here
!
! Revision 2.51  2018/08/17 23:51:41  pwagner
! May use Where component in MLSFile_T
!
! Revision 2.50  2018/08/13 22:31:28  pwagner
! Add MLSChunk_T to toc; correct spelling error in comments
!
! Revision 2.49  2018/08/03 23:18:56  vsnyder
! Add L2CF and Source components to MLSFile_t
!
! Revision 2.48  2018/05/11 21:23:27  pwagner
! Stop Use-ing HDF; Moved M_ mask bit fields here
!
! Revision 2.47  2018/02/08 23:18:00  pwagner
! moved accessType and split_path_name here from MLSFiles
!
! Revision 2.46  2015/06/30 18:39:10  pwagner
! May keep list of module names to exempt from annoying crashes
!
! Revision 2.45  2015/06/19 00:32:38  pwagner
! Moved MLSChunk_T here
!
! Revision 2.44  2014/12/10 19:11:11  pwagner
! f.p. inRange now includes its endpoints like integer interface
!
! Revision 2.43  2014/09/05 00:00:27  vsnyder
! Remove ProcessID
!
! Revision 2.42  2014/08/05 00:17:07  pwagner
! Add --pId and --uId to set id_strings for slave task
!
! Revision 2.41  2014/03/26 17:42:52  pwagner
! Added ProductionLocation, identifier_product_DOI to metadata
!
! Revision 2.40  2014/03/07 19:11:33  pwagner
! Distinguish between ShortNameLen and NameLen
!
! Revision 2.39  2013/11/18 21:41:51  pwagner
! Sticky versions of verbose, debug available
!
! Revision 2.38  2013/08/30 23:12:21  pwagner
! Default values for all MLSFile_T fields
!
! Revision 2.37  2011/12/13 01:06:03  pwagner
! Generalized Fill by MLSFill type; added MLSVerbose and MLSDebug
!
! Revision 2.36  2011/08/26 00:27:48  pwagner
! Added Interval_T
!
! Revision 2.35  2010/11/30 00:33:27  pwagner
! Added instance for character arg to is_what_ieee
!
! Revision 2.34  2010/01/11 18:32:44  pwagner
! Completed toc; generic is_what_ieee_ now use r4, r8
!
! Revision 2.33  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.32  2009/06/02 17:47:26  cvuu
! Add L2Metadata structure
!
! Revision 2.31  2008/01/09 20:50:01  pwagner
! is_what_ieee supports integer args
!
! Revision 2.30  2008/01/07 21:33:19  pwagner
! Added ieee signal defs and id-ing functions
!
! Revision 2.29  2007/01/12 00:24:38  pwagner
! Tore ourselves loose from SDPToolkit, hoping we left no shreds of flesh
!
! Revision 2.28  2006/11/01 20:31:38  pwagner
! House-cleaning
!
! Revision 2.27  2005/12/16 00:02:05  pwagner
! FillValue-related stuff moved to new MLSFillValues module
!
! Revision 2.26  2005/10/19 22:53:01  vsnyder
! Move kinds to MLSKinds
!
! Revision 2.25  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.24  2005/06/14 20:31:10  pwagner
! Added and changed some fields of MLSFile_T
!
! Revision 2.23  2005/05/31 17:49:15  pwagner
! Added new fields to MLSFile_T
!
! Revision 2.22  2005/05/12 20:46:46  pwagner
! Added filterValues and isFinite procedures (Should they be elsewhere?)
!
! Revision 2.21  2004/08/03 17:58:25  pwagner
! Now holds DEFAULTUNDEFINEDVALUE to be used elsewhere
!
! Revision 2.20  2004/06/10 01:00:50  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.19  2004/05/19 19:16:40  vsnyder
! Move MLSChunks_t to Chunks_m
!
! Revision 2.18  2004/01/09 00:38:04  pwagner
! Added FindNext function
!
! Revision 2.17  2003/06/20 19:31:39  pwagner
! Changes to allow direct writing of products
!
! Revision 2.16  2003/02/17 03:52:49  livesey
! Bit the bullet and changed rm to r4.
!
! Revision 2.15  2002/12/05 19:44:24  pwagner
! Moved MLSFile_T from MLSFiles to MLSCommon
!
! Revision 2.14  2002/11/06 00:16:48  pwagner
! Added toc/api blocks
!
! Revision 2.13  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.12  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.11  2002/08/28 22:16:18  pwagner
! Added rm, rv types
!
! Revision 2.10  2002/02/19 23:10:54  pwagner
! Added BareFNLen
!
! Revision 2.9  2002/01/09 23:51:27  pwagner
! Connected FileNameLen with PGSd_PC_FILE_PATH_MAX
!
! Revision 2.8  2001/11/14 18:03:32  livesey
! Changed FindFirst to return 0 not -1 if not found
!
! Revision 2.7  2001/09/09 02:47:58  livesey
! Moved FindFirst into MLSCommon
!
! Revision 2.6.2.3  2001/09/09 01:53:27  livesey
! Bug fix
!
! Revision 2.6.2.2  2001/09/09 01:35:46  livesey
! Moved FindFirst in from MLSL2Common
!
! Revision 2.6.2.1  2001/09/08 22:32:24  livesey
! Added RP and IP
!
! Revision 2.6  2001/04/20 23:10:53  livesey
! Initialised parameters in L1BINFO
!
! Revision 2.5  2001/03/10 18:48:17  livesey
! Really nullified the pointer!
!
! Revision 2.4  2001/03/10 07:06:46  livesey
! Nullified L1BRadfileNames in L1BInfo
!
! Revision 2.3  2001/02/09 00:38:55  livesey
! Various changes
!
! Revision 2.2  2001/01/26 23:46:35  pwagner
! Restored L1BInfo from l1/MLSL1Common back to lib/MLSCommon
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.14  2000/09/02 01:58:30  vsnyder
! Cosmetic changes
!
