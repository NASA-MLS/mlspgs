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

  use IEEE_ARITHMETIC, only: IEEE_IS_FINITE, IEEE_IS_NAN
  use MLSKINDS ! EVERYTHING

  implicit none
  private

! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (data types and parameters)

! i1, i2, i4    integer types
! r4, r8        floating point types
! ip, rp        integer, floating point types used in forward model
! rv            floating point type used in vector quantity values
! rm            floating point type used in matrix values
! rt            floating point type used in topo-set values
! LineLen       character-length of most input
! FileNameLen   character-length of path/filenames
! BareFNLen     character-length of filenames
! namelen       Max length of hdf sds array name
! fill_signal   signal to is_what_ieee to check for undefined (fill) values
! finite_signal signal to is_what_ieee to check for finite
! inf_signal    signal to is_what_ieee to check for inf
! nan_signal    signal to is_what_ieee to check for NaN
! ShortNameLen  character-length of short names (e.g., 'H2O')
! DEFAULTUNDEFINEDVALUE 
!               default fill values, e.g. when creating hdf arrays
! UNDEFINEDVALUE 
!               value to check for by is_what_ieee
! FileIDs_T     id numbers for file, group, swath or dataset
! L1BInfo_T     L1B data file names, etc. 
!                (Should we replace these with FileIDs?)
! L2Metadata_T  Coords of (lon,lat) box to write as metadata
! MLSFile_T     File name, type, id, etc.
! MLSFill_T     Fill value data type
! MLSFills      database of MLSFill values
! Range_T       Bottom, Top of PCFID range
! Interval_T    Bottom, Top of a real interval (open)
! TAI93_Range_T start, end times in TAI93 formatted r8
! MLSVerbose    Should we print extra stuff?
! MLSDebug      Should we print even more stuff?
! MLSVerboseSticky
!               Retain value of MLSVerbose thoughout pahse
! MLSDebugSticky
!               Retain value of MLSDebug thoughout pahse

!     (subroutines and functions)
! inRange       does an argument lie with a specified range or interval
! is_what_ieee  is an argument a specified ieee type
! === (end of toc) ===                                                   
! === (start of api) ===
! log inRange( int arg, Range_T range )
! log inRange( real arg, Interval_T range )
! log is_what_ieee( type_signal what, num arg )
! === (end of api) ===

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains simple definitions that are common to all the MLS PGS
  ! f90 software.

  public :: INRANGE
  public :: IS_WHAT_IEEE

  public :: FILEIDS_T
  public :: MLSFILE_T
  public :: MLSFILL_T
  public :: L1BINFO_T
  public :: RANGE_T
  public :: TAI93_RANGE_T
  public :: L2METADATA_T
  public :: INTERVAL_T

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
  integer, public, parameter :: FINITE_SIGNAL = 0
  integer, public, parameter :: INF_SIGNAL    = FINITE_SIGNAL + 1
  integer, public, parameter :: NAN_SIGNAL    = INF_SIGNAL + 1
  integer, public, parameter :: FILL_SIGNAL   = NAN_SIGNAL + 1
  integer, public, parameter :: OLDFILL_SIGNAL= FILL_SIGNAL + 1

  ! Not just "is" but "=", because "is" may accept under broader conditions
  interface equalsFillValue
    module procedure equalsFillValue_r4, equalsFillValue_r8, equalsFillValue_int
  end interface

  interface inRange
    module procedure inRange_r4, inRange_r8, inRange_int
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
  ! Because we'd like not to always drag the SDPToolkit with us
  ! everywhere we go
  integer, private, parameter :: PGSd_PC_FILE_PATH_MAX = 1024

  ! Now we have the lengths for various strings

  integer, public, parameter :: ShortNameLen = 32
  integer, public, parameter :: namelen      = 64  ! Max len of SDS array name
  integer, public, parameter :: LineLen      = 132
  integer, public, parameter :: FileNameLen  = max(PGSd_PC_FILE_PATH_MAX, 132) ! was 132
  integer, public, parameter :: BareFNLen    = 64  ! Bare file name length (w/o path)

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
  real(r4), public, parameter :: DEFAULTUNDEFINEDVALUE = -999.99 ! Try to use in lib, l2
  real(r4), public, parameter :: UNDEFINEDTOLERANCE = 0.2 ! Poss. could make it 1

  real(r4), public, save      ::    UNDEFINEDVALUE = DEFAULTUNDEFINEDVALUE
  ! --------------------------------------------------------------------------
  
  ! A type to hold the hdf file ids
  ! (Should make it recursive in case dataset path something like
  ! "/grp_1/grp_2/../grp_n/sd")

  type FileIds_T
    integer :: f_id     = 0 ! File id, handle, or io unit
    integer :: grp_id   = 0 ! group id
    integer :: sd_id    = 0 ! sd or swath id
  end type Fileids_T

  ! A PCFid range
  type Range_T
    integer :: Bottom   = 0
    integer :: Top      = 0
  end type Range_T

  ! An open real interval
  type Interval_T
    real(rt) :: Bottom   = 0._rt
    real(rt) :: Top      = 0._rt
  end type Interval_T

  ! --------------------------------------------------------------------------

  ! Moved here from MLSFiles module
  ! Information describing the files used by the mls software
  ! Stop passing file handles or names back & forth between routines
  ! -- pass one of these instead
  type MLSFile_T
    character (LEN=16) :: content=""  ! e.g., 'l1brad', 'l2gp', 'l2aux', ..
    character (LEN=8) :: lastOperation=""  ! 'open','close','read','write'
    character (LEN=FileNameLen) :: Name=""  ! its name (usu. w/path)
    character (LEN=ShortNameLen) :: ShortName=""  ! its short name; e.g. 'H2O'
    character (LEN=8) :: typeStr=""  ! one of {'swath', 'hdf', ..}
    integer :: type=0  ! one of {l_swath, l_hdf, ..}
    integer :: access=0  ! one of {DFACC_RDONLY, DFACC_CREATE, ..}
    integer :: HDFVersion=0  ! its hdf version if hdf(eos)
    integer :: PCFId=0      ! its PCF ID (ref), if any
    integer :: recordLength=0! its max record_length, if any
    integer :: errorCode=0  ! non-zero usu. means trouble
    logical :: StillOpen=.false.
    type(Range_T) :: PCFidRange = Range_T()
    type(Fileids_T) :: FileID = FileIDs_T()
  end type MLSFile_T

  ! The next datatype describes the information on the L1B data files in use

  type L1BInfo_T
    integer :: L1BOAId=0     ! The HDF ID (handle) for the L1BOA file
    ! Id(s) for the L1BRAD file(s)
    integer, dimension(:), pointer :: L1BRADIds=>NULL()
    character (LEN=FileNameLen) :: L1BOAFileName=""  ! L1BOA file name
    character (LEN=FileNameLen), dimension(:), pointer :: &
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
  ! Datatype for the Fill 
  ! (its value(s), whether Fills
  ! mean a value greater than or equal to, etc.)
  ! We will be migrating older is-Fill tests to use this instead
  ! to handle multiple Fill values and idea of range
  
  ! The idea is that a tested value is a fill value if it satisfies any
  ! of a set of conditions with respect to reference values
  ! This is an obvious generalization of the older test where we
  ! tested only whether it equaled one specified reference value
  
  ! Why do this? We use -999.99, both in level 1 and level 2
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

  elemental function inRange_int(arg, range) result(relation)
    ! Is arg in range?
    integer, intent(in)       :: arg
    type(Range_T), intent(in) :: range
    logical                   :: relation
    relation = (arg < (range%top + 1)) .and. (arg > (range%bottom - 1))
  end function inRange_int

  elemental function inRange_r4(arg, range) result(relation)
    ! Is arg in range?
    real(r4), intent(in)       :: arg
    type(Interval_T), intent(in) :: range
    logical                   :: relation
    relation = (arg < range%top) .and. (arg > range%bottom)
  end function inRange_r4

  elemental function inRange_r8(arg, range) result(relation)
    ! Is arg in range?
    real(r8), intent(in)       :: arg
    type(Interval_T), intent(in) :: range
    logical                   :: relation
    relation = (arg < range%top) .and. (arg > range%bottom)
  end function inRange_r8

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
      itIs = abs(arg - undefinedValue) < &
        & max( Real(undefinedTolerance, r8), abs(undefinedValue/100000._r8) )
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
      itIs = ( abs(arg-int(undefinedvalue)) < &
        & max(undefinedTolerance, abs(undefinedvalue/100000)) )
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
      itIs = ( abs(arg-int(undefinedvalue)) < 1 )
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
  
  elemental function is_a_fill_value_int( arg ) result ( itIs )
    integer, intent(in) :: arg
    logical :: itIs
    integer :: k
    integer :: x
    include 'isafillvalue.f9h'
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
