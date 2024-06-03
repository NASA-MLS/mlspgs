! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MLSHDF5

  ! This module contains MLS specific routines to do lowish level common HDF5
  ! tasks.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dump_Options, only: Dopt_Laconic, Dopt_RMS, Dopt_Stats, Dopt_Verbose, &
    & NameOnEachLine
  use Dump_0, only: Dump
  use Dump_1, only: DumpNamedValues
  use Hdf, only: DFACC_RDOnly
  use HighOutput, only: OutputNamedValue
  use Intrinsic, only: L_Hdf
  use Machine, only: Crash_Burn
  use MLSCommon, only: MLSFile_T
  use MLSDataInfo, only: MLSDataInfo_T, Query_MLSData
  use MLSFiles, only: HDFVersion_5, &
    & Dump, InitializeMLSFile, MLS_CloseFile, MLS_OpenFile
  use MLSKinds, only: R8
  use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Info, MLSMSG_Warning, &
    & MLSMessage
  use MLSFinds, only: FindFirst
  use MLSStringLists, only: CatLists, IsInList, &
    & GetStringElement, Intersection, NumStringElements, StringElement
  use MLSStrings, only: Indexes, Lowercase, Replace, Trim_Safe
  use Optional_M, only: Print_Default
  use Output_M, only: Newline, Output
  use Trace_M, only: Trace_Begin, Trace_End
  ! Let's Break Down Our Use, Parameters First
  use HDF5, only: H5f_Acc_Rdonly_F, H5f_Acc_Rdwr_F, &
    & H5G_Group_F, H5p_Dataset_Create_F, &
    & H5sis_Simple_F, & ! H5soffset_Simple_F, &
    & H5s_Scalar_F, H5s_Select_Set_F, H5s_Unlimited_F, &
    & H5t_Ieee_F32le, H5t_Ieee_F64le, H5t_Ieee_F64le, &
    & H5t_Native_Double, H5t_Native_Real, H5t_Std_I32le, &
    & H5t_Native_Character, H5t_Native_Integer, H5t_String, H5t_String_F, &
    & Hid_T, Hsize_T, Size_T
  ! Now Routines
  use HDF5, only: H5aclose_F, H5acreate_F, &
    & H5aget_Name_F, H5aget_Num_Attrs_F, &
    & H5aget_Space_F, H5aget_Type_F, H5aopen_Idx_F, H5aopen_Name_F, &
    & H5aread_F, H5awrite_F, H5adelete_F, &
    & H5dcreate_F, H5dextend_F, H5dget_Space_F, H5dget_Type_F, H5dopen_F, &
    & H5dread_F, H5dwrite_F, H5dclose_F, H5dget_Create_Plist_F, &
    & H5eSet_Auto_F, &
    & H5fopen_F, H5fclose_F, &
    & H5gopen_F, H5gclose_F, h5gCreate_f, &
    & H5gn_Members_F, H5gget_Obj_Info_Idx_F, &
    & H5pcreate_F, H5pset_Chunk_F, H5pset_Fill_Value_F, &
    & H5pget_Chunk_F, H5pget_Fill_Value_F, &
    & H5sclose_F, &
    & H5screate_F, H5screate_Simple_F, H5sget_Simple_Extent_Ndims_F, &
    & H5sget_Simple_Extent_Dims_F, H5sselect_Hyperslab_F, &
    & H5tclose_F, H5tcopy_F, H5tequal_F, H5tget_Class_F, H5tget_Size_F, &
    & H5tset_Size_F

  implicit none
  private

  public :: CpHDF5Attribute, CpHDF5GlAttribute, &
    & DumpHDF5Attributes, DumpHDF5DS, &
    & GetAllHDF5AttrNames, GetAllHDF5DSNames, GetAllHDF5GroupNames, &
    & GetHDF5Attribute, GetHDF5AttributePtr, GetHDF5AttrDims, &
    & GetHDF5DSRank, GetHDF5DSDims, GetHDF5DSQType, &
    & IsHDF5AttributeInFile, IsHDF5AttributePresent, IsHDF5DSInFile, &
    & IsHDF5DSPresent, IsHDF5GroupPresent, IsHDF5ItemPresent, &
    & LoadAllocFromHDF5DS, LoadFromHDF5DS, LoadPtrFromHDF5DS, &
    & MakeHDF5Attribute, MakeNestedGroups, MatchHDF5Attributes, &
    & MLS_H5Open, MLS_H5Close, ReadHDF5Attribute, &
    & ReadLitIndexFromHDF5Attr, ReadStringIndexFromHDF5Attr, SaveAsHDF5DS, &
    & WriteLitIndexAsHDF5Attribute, WriteStringIndexAsHDF5Attribute

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! CpHDF5Attribute      Copies an attribute
! CpHDF5GlAttribute    Copies a global attribute
! DumpHDF5Attributes   Dumps attributes
! DumpHDF5DS           Dumps datasets
! GetAllHDF5AttrNames  Retrieves names of all attributes under itemID
! GetAllHDF5DSNames    Retrieves names of all DS in file (under group name)
! GetAllHDF5GroupNames Retrieves names of all groups in file (under group name)
! GetHDF5Attribute     Retrieves an attribute
! GetHDF5AttributePtr  Allocates an array for an attribute and retrieves it
! GetHDF5DSRank        How many dimensions in dataset
! GetHDF5DSDims        Size of the dimensions in dataset
! GetHDF5DSQType       What datatype is dataset?
! IsHDF5...Present     Is the (attribute, DS) in the locid?
! IsHDF5...InFile      Is the (attribute, DS) in the named file?
! LoadFromHDF5DS       Retrieves a dataset
! LoadPtrFromHDF5DS    Allocates an array and retrieves a dataset
! MakeHDF5Attribute    Turns an arg into an attribute
! MakeNestedGroups     Create a nested sequence of groups; e.g. a/b/c/../z
! MatchHDF5Attributes  Finds a dataset with matching attributes
! mls_h5close          Closes interface to hdf5; call once at end of run
! mls_h5open           Opens interface to hdf5; call once at start of run
! SaveAsHDF5DS         Turns an array into a dataset
! === (end of toc) ===

! === (start of api) ===
! CpHDF5Attribute (int fromitemID, int toitemID, char name,
!    [log skip_if_already_there])
! CpHDF5GlAttribute (char fromFile, char toFile, char name,
!    [log skip_if_already_there])
! DumpHDF5Attributes (int locID, [char names], [char groupName], [char DSName],
!    [char options])
! DumpHDF5DS (int locID, char groupame, [char names],
!    [real fillValue], [char options])
! GetAllHDF5AttrNames (int itemID, char DSNames)
! GetAllHDF5DSNames (file, char gname, char DSNames)
!     file can be one of:
!    {char* filename, int fileID}
! GetHDF5Attribute (int itemID, char name, value)
! GetHDF5AttributePtr (int itemID, char name, *value)
! GetHDF5DSDims (int FileID, char name, hsize_t dims(:), [hsize_t maxdims(:)])
! GetHDF5DSRank (int FileID, char name, int rank)
! GetHDF5DSQType (int FileID, char name, char QType)
! log IsHDF5AttributeInFile (char filename, char DSname, char name)
! log IsHDF5AttributeInFile (char filename, char name)
! log IsHDF5DSInFile (char filename, char name)
! log IsHDF5AttributePresent (int setid, char name)
! log IsHDF5AttributePresent (int fileid, char DSname, char name)
! log IsHDF5DSPresent (int locID, char name, [options])
! log IsHDF5GroupPresent (int locID, char name)
! log IsHDF5ItemPresent (int locID, char name, [options])
! LoadAllocFromHDF5DS (int locID, char name, *value [, lowBound] )
!       lowBound only for rank-1 arrays
! LoadFromHDF5DS (int locID, char name, value,
!       [int start(:), int count(:), [int stride(:), int block(:)] ] )
! LoadPtrFromHDF5DS (int locID, char name, *value [, lowBound] )
!       lowBound only for rank-1 arrays
! MakeHDF5Attribute ( int itemID, char name, value,
!       [log skip_if_already_there] )
! MakeNestedGroups ( int itemID, char groupnames(:) )
! MatchHDF5Attributes ( MLSFile_T MLSFile, char attrnames, char attrvalues,
!       char name )
! SaveAsHDF5DS ( int locID, char name, value,
!       [int start(:), int count(:), [int stride(:), int block(:)] ] )
!     value can be one of:
!    {char* value, int value, real value, double precision value,
!     int value(:),
!     real value(:), real value(:,:), real value(:,:,:), real value(:,:,:,:),
!     double precision value(:), double precision value(:,:),
!     double precision value(:,:,:)}

! Some of these may also be called with the int itemID or locID arg replaced
! by a MLSFile_T; e.g.
! GetHDF5Attribute ( MLSFile_T MLSFile, char name, value )
!        The attribute will be assumed to be an attribute of
!         MLSFile%FileID%sd_id
!         unless that component is zero, when it will try ..%grp_id
! LoadFromHDF5DS ( MLSFile_T MLSFile, char name, value,
!       [int start(:), int count(:), [int stride(:), int block(:)] ] )

! One standard is the character flag "options" which affects how loosely
! string matches may be interpreted
! it may include any of the following (poss. in combination, e.g. "-wc")
! w    Wildcard * which allows 'a*' to equal 'abcd'
! c    case insensitive which allows 'ABCD' to equal 'abcd'
! f    flush left which allows 'abcd' to equal '  abcd'
! In functions with "ITEM" in the title; for these expanded use allows
! a    Search for matches among attributes
! d    Search for matches among datasets
! g    Search for matches among groups
! (Don't confuse this use of options with the one below: here options
! affects how loosely strings may match; the one below affects how much
! and what details appear in dumps)

! The meaning of options has replaced the older logical arguments in Dumps
! if the options is present and contains the following characters:
!   character         meaning
!      ---            -------
!       c              Clean
!       L              Laconic
!       r              RMS       
!       s              Stats     
!       u              Unique    
!       w              WholeArray
! etc. (see dump_0 module)

! Note that when an argument is of type "char" and holds multiple items,
! e.g., DSnames in  GetAllHDF5DSNames or attrNames in MatchHDF5Attributes,
! it will be a comma-separated string list, as described in the module
! MLSStringLists
!
! Bugs and Gotchas:
! (1) How do we know where the wildcard "*" can appear in a field 
!    and where it is forbidden?
! (2) For which procedures can the arg be an MLSFile_T? Shouldn't you make it
!    all of them?
! (3) Beware that start uses the c position that array indexes start at 0, not 1.
!    It might be better named "offset". As it is, start is an unfortunate
!    misnomer that will probably continue to confuse.
!                      
! === (end of api) ===
  interface CpHDF5Attribute
    ! module procedure CpHDF5Attribute_int, CpHDF5Attribute_logical
    module procedure CpHDF5Attribute_string
    ! module procedure CpHDF5Attribute_sngl, CpHDF5Attribute_dbl
    ! module procedure CpHDF5Attribute_snglarr1, CpHDF5Attribute_intarr1
    ! module procedure CpHDF5Attribut _dblarr1, CpHDF5Attribute_stringarr1
  end interface

  interface CpHDF5GlAttribute
    module procedure CpHDF5GlAttribute_string
  end interface

  interface GetAllHDF5DSNames
    module procedure GetAllHDF5DSNames_fileID, GetAllHDF5DSNames_filename
    module procedure GetAllHDF5DSNames_MLSFile
  end interface

  interface GetHDF5Attribute
    module procedure GetHDF5Attribute_int, GetHDF5Attribute_logical, &
      & GetHDF5Attribute_logicalarr1, &
      & GetHDF5Attribute_string, GetHDF5Attribute_sngl, GetHDF5Attribute_dbl, &
      & GetHDF5Attribute_snglarr1, GetHDF5Attribute_intarr1, &
      & GetHDF5Attribute_dblarr1, GetHDF5Attribute_stringarr1
    module procedure GetHDF5Attr_ID_int, GetHDF5Attr_ID_logical, &
      & GetHDF5Attr_ID_logicalarr1, &
      & GetHDF5Attr_ID_string, GetHDF5Attr_ID_sngl, GetHDF5Attr_ID_dbl, &
      & GetHDF5Attr_ID_snglarr1, GetHDF5Attr_ID_intarr1, &
      & GetHDF5Attr_ID_dblarr1, GetHDF5Attr_ID_stringarr1
  end interface

  interface GetHDF5AttributePtr
      module procedure GetHDF5AttributePtr_snglarr1, GetHDF5AttributePtr_intarr1, &
          & GetHDF5AttributePtr_dblarr1, GetHDF5AttributePtr_stringarr1, &
          & GetHDF5AttributePtr_logicalarr1
  end interface

  interface ReadHDF5Attribute
      module procedure ReadHDF5Attr_FID_int, ReadHDF5Attr_FID_string
  end interface

  interface IsHDF5AttributeInFile
      module procedure IsHDF5AttributeInFile_DS, IsHDF5AttributeInFile_Grp
  end interface

  interface GetHDF5DSDims
    module procedure GetHDF5DSDims_ID, GetHDF5DSDims_MLSFile
  end interface

  interface GetHDF5DSRank
    module procedure GetHDF5DSRank_ID, GetHDF5DSRank_MLSFile
  end interface

  interface IsHDF5DSPresent
    module procedure IsHDF5DSPresent_in_fID, IsHDF5DSPresent_in_MLSFile
  end interface

  interface IsHDF5AttributePresent
    module procedure IsHDF5AttributePresent_in_fID, &
      & IsHDF5AttributePresent_in_MLSFile, &
      & IsHDF5AttributePresent_in_DSID, IsHDF5AttributePresent_in_grp
  end interface

  interface LoadAllocFromHDF5DS
    module procedure LoadAllocFromHDF5DS_chararr1, LoadAllocFromHDF5DS_chararr2, &
      & LoadAllocFromHDF5DS_intarr1, LoadAllocFromHDF5DS_intarr2, &
      & LoadAllocFromHDF5DS_intarr3, LoadAllocFromHDF5DS_intarr4, &
      & LoadAllocFromHDF5DS_logarr1, &
      & LoadAllocFromHDF5DS_dblarr1, LoadAllocFromHDF5DS_dblarr2, &
      & LoadAllocFromHDF5DS_dblarr3, LoadAllocFromHDF5DS_dblarr4, &
      & LoadAllocFromHDF5DS_snglarr1, LoadAllocFromHDF5DS_snglarr2, &
      & LoadAllocFromHDF5DS_snglarr3, LoadAllocFromHDF5DS_snglarr4
  end interface

  interface LoadFromHDF5DS
    module procedure LdFrmHDF5DS_ID_intarr1, LdFrmHDF5DS_ID_intarr2, &
      & LdFrmHDF5DS_ID_intarr3, LdFrmHDF5DS_ID_intarr4, &
      & LdFrmHDF5DS_ID_logarr1, &
      & LdFrmHDF5DS_ID_dblarr1, LdFrmHDF5DS_ID_dblarr2, &
      & LdFrmHDF5DS_ID_dblarr3, LdFrmHDF5DS_ID_dblarr4, &
      & LdFrmHDF5DS_ID_snglarr1, LdFrmHDF5DS_ID_snglarr2, &
      & LdFrmHDF5DS_ID_snglarr3, LdFrmHDF5DS_ID_snglarr4, &
      & LdFrmHDF5DS_ID_chararr1, LdFrmHDF5DS_ID_chararr2, &
      & LdFrmHDF5DS_ID_chararr3, LdFrmHDF5DS_ID_charscalar
    module procedure LoadFromHDF5DS_intarr1, LoadFromHDF5DS_intarr2, &
      & LoadFromHDF5DS_intarr3, LoadFromHDF5DS_intarr4, &
      & LoadFromHDF5DS_logarr1, &
      & LoadFromHDF5DS_dblarr1, LoadFromHDF5DS_dblarr2, &
      & LoadFromHDF5DS_dblarr3, LoadFromHDF5DS_dblarr4, &
      & LoadFromHDF5DS_snglarr1, LoadFromHDF5DS_snglarr2, &
      & LoadFromHDF5DS_snglarr3, LoadFromHDF5DS_snglarr4, &
      & LoadFromHDF5DS_chararr1, LoadFromHDF5DS_chararr2, &
      & LoadFromHDF5DS_chararr3, LoadFromHDF5DS_charscalar
  end interface

  interface LoadPtrFromHDF5DS
    module procedure LoadPtrFromHDF5DS_chararr1, LoadPtrFromHDF5DS_chararr2, &
      & LoadPtrFromHDF5DS_intarr1, LoadPtrFromHDF5DS_intarr2, &
      & LoadPtrFromHDF5DS_intarr3, LoadPtrFromHDF5DS_intarr4, &
      & LoadPtrFromHDF5DS_logarr1, &
      & LoadPtrFromHDF5DS_dblarr1, LoadPtrFromHDF5DS_dblarr2, &
      & LoadPtrFromHDF5DS_dblarr3, LoadPtrFromHDF5DS_dblarr4, &
      & LoadPtrFromHDF5DS_snglarr1, LoadPtrFromHDF5DS_snglarr2, &
      & LoadPtrFromHDF5DS_snglarr3, LoadPtrFromHDF5DS_snglarr4
  end interface

  interface MakeHDF5Attribute
    module procedure MakeHDF5Attribute_dbl, MakeHDF5Attribute_sngl, &
      & MakeHDF5Attribute_int, MakeHDF5Attribute_logical, &
      & MakeHDF5Attribute_logicalarr1, &
      & MakeHDF5Attribute_string, MakeHDF5Attribute_snglarr1, &
      & MakeHDF5Attribute_dblarr1, MakeHDF5Attribute_stringarr1, &
      & MakeHDF5Attribute_intarr1, MakeHDF5AttributeDSN_int, &
      & MakeHDF5AttributeDSN_string, MakeHDF5AttributeDSN_snglarr1, &
      & MakeHDF5AttributeDSN_st_arr1, MakeHDF5AttributeDSN_dblarr1, &
      & MakeHDF5AttributeDSN_single, MakeHDF5AttributeDSN_double, &
      & MakeHDF5Attribute_textFile
  end interface

  interface SaveAsHDF5DS
    module procedure &
      & SaveAsHDF5DS_intarr1, SaveAsHDF5DS_intarr2, SaveAsHDF5DS_intarr3, &
      & SaveAsHDF5DS_intarr4, &
      & SaveAsHDF5DS_logarr1, &
      & SaveAsHDF5DS_dblarr1, SaveAsHDF5DS_dblarr2, SaveAsHDF5DS_dblarr3, &
      & SaveAsHDF5DS_dblarr4, &
      & SaveAsHDF5DS_snglarr1, SaveAsHDF5DS_snglarr2, SaveAsHDF5DS_snglarr3, &
      & SaveAsHDF5DS_snglarr4, &
      & SaveAsHDF5DS_charsclr, SaveAsHDF5DS_chararr1, SaveAsHDF5DS_chararr2, &
      & SaveAsHDF5DS_textfile
  end interface

  integer, public, parameter :: MAXCHFIELDLENGTH = 2000000 ! max number of chars in l2cf
  integer, public, parameter :: MAXCHATTRLENGTH  =   40000 ! max number in attr
  integer, public, parameter :: MAXNDSNAMES = 20000   ! max number of DS names in a file

  ! Local parameters
  character, parameter :: Digits(7) = (/ '1', '2', '3', '4', '5', '6', '7' /)
  integer(hsize_t), dimension(7) :: ones = (/1,1,1,1,1,1,1/)
  logical, parameter :: countEmpty = .true.
  logical, parameter :: DEEBUG = .false.
  integer, save      :: cantGetDataspaceDims = 0
  integer, parameter :: DFLTMAXLINELENGTH = 1024
  integer, parameter :: MAXNUMWARNS = 40
  integer, parameter :: MAXATTRIBUTESIZE =   40000
  integer, parameter :: MAXTEXTSIZE      = 2000000
  integer, parameter :: MAXNAMELEN = 64
  logical, public    :: hdfVerbose = .false.
  ! logical, public    :: MayCrash = .false.
  character(len=*), dimension(2), parameter :: DONTDUMPTHESEDSNAMES = (/ &
    & 'wtfcoremetadata', 'wtfxmlmetadata ' /)
  ! Local variables
  integer(hid_t) :: cparms

contains ! ======================= Public Procedures =========================

  ! ------------------------------------------------  MLS_h5close  -----
  subroutine MLS_h5close ( error )
    ! Arguments
  ! To switch to/from hdfeos5.1.6(+) uncomment next line
    use H5Lib, only: H5Close_f
    integer, intent(out) :: error          ! Trouble if /= 0
    error = 0
    call h5close_f ( error )
  end subroutine MLS_h5close

  ! -------------------------------------------------  MLS_h5open  -----
  subroutine MLS_h5open ( error )
    ! Arguments
  ! To switch to/from hdfeos5.1.6(+) uncomment next line
    use H5Lib, only: H5Open_f
    integer, intent(out) :: error          ! Trouble if /= 0
    error = 0
    call h5open_f ( error )
  end subroutine MLS_h5open

  ! -------------------------------------  CpHDF5Attribute_string  -----
  subroutine CpHDF5Attribute_string ( fromitemID, toitemID, name, &
   & skip_if_already_there )
    integer, intent(in)           :: FROMITEMID ! Group etc. to Cp attr from
    integer, intent(in)           :: TOITEMID   ! Group etc. to Cp attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(in), optional :: skip_if_already_there ! Or if not in fromItemID

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_skip
    logical :: is_present
    character (len=4096) :: value1

    ! Executable code
    call trace_begin ( me, 'CpHDF5Attribute_string', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip = skip_if_already_there
    is_present = IsHDF5AttributePresent_in_DSID( fromitemID, name )
    if ( my_skip .and. .not. is_present ) go to 9
    if ( .not. is_present ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to cp: attribute not found ' // trim(name) )
    is_present = IsHDF5AttributePresent_in_DSID(toitemID, name)
    if ( my_skip .and. is_present ) go to 9
    call GetHDF5Attribute ( fromitemID, name, value1 )
    call MakeHDF5Attribute ( toitemID, name, trim(value1), &
      & skip_if_already_there )
  9 call trace_end ( cond=.false. )
  end subroutine CpHDF5Attribute_string

  ! -----------------------------------  CpHDF5GlAttribute_string  -----
  subroutine CpHDF5GlAttribute_string ( fromFileName, toFileName, name, &
   & skip_if_already_there )
    character (len=*), intent(in) :: FROMFILENAME     ! file name
    character (len=*), intent(in) :: TOFILENAME       ! file name
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: fromfileID
    integer :: fromgrpID
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag
    integer :: tofileID
    integer :: togrpID

    ! Executable code
    call trace_begin ( me, 'CpHDFGlAttribute_string', cond=.false. )
    call h5fopen_f ( trim(fromFilename), H5F_ACC_RDONLY_F, fromFileID, status )
    call h5fopen_f ( trim(toFilename), H5F_ACC_RDWR_F, toFileID, status )
    call h5gopen_f ( fromFileID, '/', fromgrpid, status )
    call h5gopen_f ( toFileID, '/', togrpid, status )
    call CpHDF5Attribute_string ( fromgrpID, togrpID, name, &
      & skip_if_already_there )
    call h5gclose_f ( fromgrpid, status )
    call h5gclose_f ( togrpid, status )
    call h5fclose_f ( fromFileID, status )
    call h5fclose_f ( toFileID, status )
    call trace_end ( 'CpHDFGlAttribute_string', cond=.false. )
  end subroutine CpHDF5GlAttribute_string

  ! -------------------------------  DumpHDF5Attributes  -----
  subroutine DumpHDF5Attributes ( locID, names, groupName, DSName, options )
    ! Dump attributes of locID (possibly followed by a group or DS name)
    ! All of them or only those in names string list
    integer, intent(in)                     :: locID ! attributes of what
    character (len=*), intent(in), optional :: NAMES  ! only these names
    character(len=*), intent(in), optional  :: groupName
    character(len=*), intent(in), optional  :: DSName
    character(len=*), intent(in), optional :: options

    ! Local variables
    integer :: attrID
    integer :: classID
    character(len=MAXCHFIELDLENGTH) :: chValue ! 1024; len may become MAXCHFIELDLENGTH
    ! logical, parameter :: DEEBUG = .true.
    integer, dimension(7) :: dims
    double precision, dimension(1024) :: dValue
    integer(kind=hSize_t), dimension(7) :: hdims
    integer :: i
    integer :: itemID
    integer, dimension(1024) :: iValue
    integer(kind=hSize_t), dimension(:), allocatable :: maxdims_ptr
    integer :: Me = -1                  ! String index for trace cacheing
    character(len=MAXNDSNAMES*MAXNAMELEN) :: myNames
    character(len=128) :: name
    integer :: numAttrs
    character(len=16) :: Qtype
    integer :: rank
    integer :: spaceID
    integer :: status
    integer :: type_id
    integer(kind=Size_t) :: type_size
    ! Executable
    call trace_begin ( me, 'DumpHDF5Attributes', cond=.false. )
    myNames = '*' ! Wildcard means 'all'
    if ( present(names) ) myNames = names
    if ( present(groupName) ) then
      ! call h5gOpen_f ( locID, '/', itemID, status )
      call h5gOpen_f ( locID, trim(groupName), itemID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open group' // trim(groupName) // &
        & ' while dumping its attributes' )
      if ( DEEBUG ) call outputNamedValue ( 'groupName', groupName )
    elseif(present(DSName) ) then
      call h5dOpen_f ( locID, trim(DSname), itemID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open dataset' // trim(DSName) // &
        & ' while dumping its attributes' )
      if ( DEEBUG ) call outputNamedValue ( 'DSName', DSName )
    else
      itemID = locID
      if ( DEEBUG ) call outputNamedValue ( 'locID', locID )
    endif
    if ( myNames == '*' ) then
      if ( present(groupName) ) then
        call GetAllHDF5AttrNames ( locID, myNames, groupname=groupName )
      else
        call GetAllHDF5AttrNames ( itemID, myNames )
      endif
    endif
    numAttrs = NumStringElements ( myNames, countEmpty )
    if ( DEEBUG ) call outputNamedValue ( 'myNames', myNames )
    if ( DEEBUG ) call outputNamedValue ( 'numAttrs', numAttrs )
    do i = 1, numAttrs
      if ( DEEBUG ) call outputNamedValue ( 'i', i )
      ! name = StringElement(myNames, i, countEmpty)
      call GetStringElement( myNames, name, i, countEmpty )
      if ( DEEBUG ) call outputNamedValue ( 'name', name )
      call h5aopen_name_f ( itemID, trim(name), attrID, status )
      if ( status /= 0 ) then
        call output ( trim(name), advance='no' )
        call output ( '  (not found)', advance='yes' )
        cycle
      endif
      call h5aGet_type_f ( attrID, type_id, status )
      call h5tGet_size_f ( type_id, type_size, status )
      call h5tget_class_f ( type_id, classID, status )
      call h5aget_space_f ( attrID, spaceID, status )
      call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
      allocate ( maxdims_ptr(rank) )
      call h5sget_simple_extent_dims_f ( spaceID, hdims(1:rank), &
           maxdims_ptr, status )
      if ( status /= rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dims for attr ' // trim(name) )
      hdims(1:rank) = maxdims_ptr
      deallocate ( maxdims_ptr )
      if ( DEEBUG ) then
        call outputNamedValue ( 'name', name )
        call outputNamedValue ( 'attrID', attrID )
        call outputNamedValue ( 'type_id', type_id )
        call outputNamedValue ( 'type_size', int(type_size) )
        call outputNamedValue ( 'classID', classID )
        call outputNamedValue ( 'H5T_STRING_F', H5T_STRING_F )
        call outputNamedValue ( 'spaceID', spaceID )
        call outputNamedValue ( 'rank', rank )
        call outputNamedValue ( 'hdims', int(hdims, kind(status)) )
      end if
      call h5aClose_f ( attrID, status )
      call GetHDF5AttrDims ( itemID, name, hdims )
      dims = hdims
      Qtype = WhatTypeAmI ( type_id )
      if ( Qtype == 'unknown' .and. classID == H5T_STRING_F ) QType='character'
      if ( DEEBUG ) then
        call outputNamedValue ( 'dims', dims )
        call outputNamedValue ( 'Qtype', Qtype )
      end if
      dims(1) = max(dims(1), 1)
      select case ( QType )
      case ( 'integer' )
        call GetHDF5Attribute ( itemID, name, iValue(1:dims(1)) )
        call dump ( iValue(1:dims(1)), trim(name), options=options )
      case ( 'double', 'real' )
        call GetHDF5Attribute ( itemID, name, dValue(1:dims(1)) )
        call dump ( dValue(1:dims(1)), trim(name), options=options )
      case ( 'character' )
        call GetHDF5Attribute ( itemID, name, chValue )
        call output( 'Dump of ' // trim(name), advance='yes' )
        call output( trim(chValue) )
        call newLine
      case default
        call output ( trim(name), advance='no' )
        call output ( '  (unrecognized type)', advance='yes' )
      end select
      
    enddo
    if ( present(groupName) ) then
      call h5gClose_f ( itemID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close group' // trim(groupName) // &
        & ' while dumping its attributes' )
    elseif(present(DSName) ) then
      call h5dClose_f ( itemID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open dataset' // trim(DSName) // &
        & ' while dumping its attributes' )
    endif
    call trace_end ( 'DumpHDF5Attributes', cond=.false. )
  end subroutine DumpHDF5Attributes
  
  ! -------------------------------------------------  DumpHDF5DS  -----
  subroutine DumpHDF5DS ( locID, groupName, &
    & names, fillValue, options )
    ! Dump datasets in groupID
    ! All of them or only those in names string list
    integer, intent(in)                     :: locID ! file or groupID
    character(len=*), intent(in)            :: groupName ! datasets in group
    character (len=*), intent(in), optional :: NAMES   ! only these names
    real, intent(in), optional              :: fillValue ! Show % = fill
    character(len=*), intent(in), optional  :: options

    ! Local variables
    logical :: addnameOnEachLine
    integer :: ch
    character(len=MAXCHFIELDLENGTH), dimension(1) :: chValue
    character(len=MAXCHFIELDLENGTH), dimension(:), pointer :: ch1dArray => null()
    character(len=1), dimension(:,:,:), pointer :: chArray => null()
    character(len=1), dimension(:), pointer :: charsArray  => null()
    character(len=80), dimension(:), pointer :: chLongArray => null()
    integer :: classID
    ! logical, parameter :: DEEBUG = .false.
    integer, parameter :: LONGCHARARRAYTHRSHLD = 10 ! was 999
    integer, dimension(7) :: dims
    logical :: dontPrintName
    double precision, dimension(:,:,:,:), pointer :: d4Value => null()
    double precision, dimension(:,:,:), pointer :: dValue => null()
    integer :: groupID
    integer(kind=hSize_t), dimension(7) :: hdims
    integer :: i
    integer :: ItemID
    integer, dimension(:,:,:), pointer :: iValue => null()
    integer :: k
    integer :: m
    integer :: Me = -1                  ! String index for trace cacheing
    character(len=MAXNDSNAMES*MAXNAMELEN) :: myNames
    character(len=128) :: name
    character(len=128) :: namePrinted
    integer :: Np1
    integer :: numDS
    character(len=16) :: Qtype
    integer :: rank
    character(len=MAXNDSNAMES*MAXNAMELEN) :: sdNames
    logical :: skipCharValues
    integer :: spaceID
    integer :: status
    integer, dimension(size(DONTDUMPTHESEDSNAMES) ) :: theIndexes
    integer :: type_id
    integer(kind=Size_t) :: type_size
    ! logical, parameter :: DEEBUG = .true.

    ! Executable
    call trace_begin ( me, 'DumpHDF5DS', cond=.false. )
    if ( present(options) .and. DEEBUG ) call outputNamedValue( 'options', options )
    skipCharValues = .false.
    if ( present(options) ) &
      & skipCharValues = any( indexes(options, (/dopt_rms, dopt_stats/)) > 0 )
    myNames = '*' ! Wildcard means 'all'
    if ( present(names) ) myNames = names
    call GetAllHDF5DSNames ( locID, groupName, sdNames )
    ! Did we ask for any datasets by name? With wildcard?
    ! print *, 'sdNames: ', trim(sdNames)
    ! print *, 'Names:   ', trim(Names)
    if ( myNames == '*' ) then
      myNames = sdNames
    elseif( index( myNames, '*') > 0 ) then
      myNames = Intersection( sdNames, myNames, options='-w' )
    endif
    ! print *, 'Intersection:   ', trim(myNames)
    dontPrintName = .false.
    if ( present(options) ) dontPrintName = index(options, dopt_laconic) > 0
    addnameOnEachLine = .false.
    if ( present(options) ) addnameOnEachLine = index(options, dopt_verbose ) > 0
    numDS = NumStringElements ( myNames, countEmpty )
    call h5gOpen_f ( locid, trim(groupName), groupID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open group' // trim(groupName) // &
      & ' while dumping its datasets' )
    do i = 1, numDS
      name = StringElement(myNames, i, countEmpty)
      if ( name == ',' .or. len_trim(name) < 1 ) cycle
      theIndexes = indexes( name, DONTDUMPTHESEDSNAMES )
      if ( any( theIndexes > 0 ) ) cycle
      namePrinted = name
      if ( dontPrintName ) namePrinted = ' '
      if ( addnameOnEachLine ) nameOnEachLine = namePrinted
      call h5dopen_f ( groupID, trim(name), ItemID, status )
      if ( status /= 0 ) then
        call output ( trim(name), advance='no' )
        call output ( '  (not found)', advance='yes' )
        cycle
      endif
      call h5dGet_type_f ( ItemID, type_id, status )
      call h5tGet_size_f ( type_id, type_size, status )
      call h5tget_class_f ( type_id, classID, status )
      call h5dget_space_f ( itemID, spaceID, status )
      call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
      if ( DEEBUG ) then
        call outputNamedValue ( 'name', name )
        call outputNamedValue ( 'ItemID', ItemID )
        call outputNamedValue ( 'type_id', type_id )
        call outputNamedValue ( 'type_size', int(type_size) )
        call outputNamedValue ( 'classID', classID )
        call outputNamedValue ( 'H5T_STRING_F', H5T_STRING_F )
        call outputNamedValue ( 'spaceID', spaceID )
        call outputNamedValue ( 'rank', rank )
      endif
      call h5dClose_f ( ItemID, status )
      call GetHDF5DSDims ( groupID, name, hdims )
      dims = hdims
      dims(1) = max(dims(1), 1)
      Qtype = WhatTypeAmI ( type_id )
      if ( Qtype == 'unknown' .and. classID == H5T_STRING_F ) QType='character'
      if ( DEEBUG ) then
        call outputNamedValue ( 'dims', dims )
        call outputNamedValue ( 'Qtype', Qtype )
      endif
      if ( QType == 'character' .and. skipCharValues ) cycle
      select case ( QType )
      case ( 'integer' )
        select case ( rank )
        case ( 3 )
          call allocate_test( iValue, dims(1), dims(2), dims(3), name, ModuleName )
          call LoadFromHDF5DS ( groupID, name, iValue )
          if ( present(fillvalue) ) then
            call dump ( iValue, trim(namePrinted), fillValue=int(fillvalue), &
              & options=options )
          else
            call dump ( iValue, trim(namePrinted), options=options )
          endif
          call deallocate_test( iValue, name, ModuleName )
        case ( 2 )
          call allocate_test( iValue, dims(1), dims(2), 1, name, ModuleName )
          call LoadFromHDF5DS ( groupID, name, iValue(:,:,1) )
          if ( present(fillvalue) ) then
            call dump ( iValue(:,:,1), trim(namePrinted), fillValue=int(fillvalue), &
              & options=options )
          else
            call dump ( iValue(:,:,1), trim(namePrinted), options=options )
          endif
          call deallocate_test( iValue, name, ModuleName )
        case default
          call allocate_test( iValue, dims(1), 1, 1, name, ModuleName )
          call LoadFromHDF5DS ( groupID, name, iValue(:,1,1) )
          if ( present(fillvalue) ) then
            call dump ( iValue(:,1,1), trim(namePrinted), fillValue=int(fillvalue), &
              & options=options )
          else
            call dump ( iValue(:,1,1), trim(namePrinted), options=options )
          endif
          call deallocate_test( iValue, name, ModuleName )
         end select
      case ( 'double', 'real' )
        select case ( rank )
        case ( 4 )
          call outputNamedValue( 'dims', dims )
          call output( '(Presently we only dump a 3-d slice of 4-d arrays)', &
            & advance='yes' )
          call allocate_test( d4Value, dims(1), dims(2), dims(3), dims(4), &
            & name, ModuleName )
          call LoadFromHDF5DS ( groupID, name, d4Value )
          if ( present(fillvalue) ) then
            if ( present(options) .and. DEEBUG ) &
              & call outputNamedValue( 'options', options )
            call dump ( d4Value(:,:,:,:), trim(namePrinted), fillValue=real(fillvalue, r8), &
              & options=options )
          else
            if ( present(options) .and. DEEBUG ) &
              & call outputNamedValue( 'options', options )
            call dump ( d4Value(:,:,:,:), trim(namePrinted), options=options )
          endif
          call deallocate_test( d4Value, name, ModuleName )
        case ( 3 )
          call allocate_test( dValue, dims(1), dims(2), dims(3), name, ModuleName )
          call LoadFromHDF5DS ( groupID, name, dValue )
          if ( present(fillvalue) ) then
            if ( present(options) .and. DEEBUG ) &
              & call outputNamedValue( 'options', options )
            call dump ( dValue, trim(namePrinted), fillValue=real(fillvalue, r8), &
              & options=options )
          else
            if ( present(options) .and. DEEBUG ) &
              & call outputNamedValue( 'options', options )
            call dump ( dValue, trim(namePrinted), options=options )
          endif
          call deallocate_test( dValue, name, ModuleName )
        case ( 2 )
          call allocate_test( dValue, dims(1), dims(2), 1, name, ModuleName )
          call LoadFromHDF5DS ( groupID, name, dValue(:,:,1) )
          if ( present(fillvalue) ) then
            call dump ( dValue(:,:,1), trim(namePrinted), &
              & fillValue=real(fillvalue, r8), options=options )
          else
            call dump ( dValue(:,:,1), trim(namePrinted), options=options )
          endif
          call deallocate_test( dValue, name, ModuleName )
        case default
          call allocate_test( dValue, dims(1), 1, 1, name, ModuleName )
          call LoadFromHDF5DS ( groupID, name, dValue(:,1,1) )
          if ( present(fillvalue) ) then
            call dump ( dValue(:,1,1), trim(namePrinted), fillValue=real(fillvalue, r8), &
              & options=options )
          else
            call dump ( dValue(:,1,1), trim(namePrinted), options=options )
          endif
          call deallocate_test( dValue, name, ModuleName )
         end select
      case ( 'character' )
        ! call outputNamedValue( 'type', QType )
        ! call outputNamedValue( 'rank', rank )
        select case ( rank )
        case ( 3 )
          call allocate_test( chArray, dims(1), dims(2), dims(3), 'chArray', ModuleName )
          call LoadFromHDF5DS ( groupID, name, chArray )
          call dump ( chArray, trim(namePrinted) )
          call deallocate_test( chArray, 'chArray', ModuleName )
        case ( 2 )
          call allocate_test( chArray, dims(1), dims(2), 1, 'chArray', ModuleName )
          call LoadFromHDF5DS ( groupID, name, chArray(:,:,1) )
          call dump ( chArray(:,:,1), trim(namePrinted) )
          call deallocate_test( chArray, 'chArray', ModuleName )
        case ( 0 )
          chvalue = ' ' ! 
          call LoadFromHDF5DS ( groupID, Name, chvalue(1) )
          if ( len_trim(namePrinted) > 0 ) call output( 'name: ' // trim(name), advance='yes' )
          call output( trim(chValue(1)), advance='yes' )
        case default
          if ( dims(1) < 2 ) then
            ! In case we have a very long dataset, e.g. the l2cf
            chvalue = ' ' ! 
            call LoadFromHDF5DS ( groupID, Name, chvalue )
            ! call outputNamedValue ( 'len before replacing nulls', len_trim(chvalue(1)) )
            ! Unfortunately, a lot of null characters sneak into this
            ! chValue(1) = Replace( chValue(1), char(0), char(32) ) ! Replace null with space
            ! call outputNamedValue ( 'len after replacing nulls', len_trim(chvalue(1)) )
            ! call dump ( trim(chValue(1)), trim(namePrinted) )
            if ( len_trim(namePrinted) > 0 ) call output( 'name: ' // trim(name), advance='yes' )
            call output( trim(chValue(1)), advance='yes' )
          elseif ( dims(1) > LONGCHARARRAYTHRSHLD ) then
            ! In case we have a very long array, e.g. dates
            if ( DEEBUG ) call outputNamedValue( 'dims', dims )
            call allocate_test( charsArray, 80*dims(1), 'charsArray', ModuleName )
            call allocate_test( chLongArray, dims(1), 'chLongArray', ModuleName )
            chLongArray = repeat( ACHAR(0), 80 )
            call LoadFromHDF5DS ( groupID, name, chLongArray )
            charsArray = ACHAR(0)
            m = 80
            do k=1, dims(1)
              do ch=1, m
                charsArray(ch + m*(k-1)) = chLongArray(k)(ch:ch)
              enddo
            enddo
            Np1 = FindFirst( charsArray, ACHAR(0) )
            m = Np1 / dims(1)
            chLongArray = ' '
            do k=1, dims(1)
              do ch=1, m
                chLongArray(k)(ch:ch) = charsArray(ch + m*(k-1))
              enddo
            enddo
            if ( DEEBUG ) call outputNamedValue( 'Np1', Np1 )
            if ( DEEBUG ) call outputNamedValue( 'm', m )
            call dump ( chLongArray, trim(namePrinted) )
            call deallocate_test( chLongArray, 'chLongArray', ModuleName )
            call deallocate_test( charsArray, 'charsArray', ModuleName )
          else
            call allocate_test( ch1dArray, dims(1), 'ch1dArray', ModuleName )
            call LoadFromHDF5DS ( groupID, name, ch1dArray )
            call dump ( ch1dArray, trim(namePrinted) )
            call deallocate_test( ch1dArray, 'ch1dArray', ModuleName )
          endif
         end select
      case default
        call output ( trim(name), advance='no' )
        call output ( '  (unrecognized type)', advance='yes' )
      end select
      
    enddo
    call h5gClose_f ( groupID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close group' // trim(groupName) // &
      & ' while dumping its datasets' )
    call trace_end ( 'DumpHDF5DS', cond=.false. )
  end subroutine DumpHDF5DS

  ! -------------------------------  GetAllHDF5AttrNames  -----
  subroutine GetAllHDF5AttrNames ( locID, names, groupName, DSName )
    ! Get all attribute names attributes of locID 
    ! (possibly followed by a group or DS name)
    ! returns string list names
    integer, intent(in)            :: locID ! attributes of what
    character (len=*), intent(out) :: NAMES ! Names of attribute
    character(len=*), intent(in), optional  :: groupName
    character(len=*), intent(in), optional  :: DSName

    ! Local variables
    integer :: attr_id
    integer :: i
    integer :: itemID
    integer :: Me = -1                  ! String index for trace cacheing
    character(len=128) :: name
    integer(kind=Size_t) :: namelength
    integer :: num
    integer :: status
    character(len=len(names)) :: tempnames

    ! Executable code
    call trace_begin ( me, 'GetAllHDF5AttrNames', cond=.false. )
    namelength = len(name)
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for all attr names' )
    if ( present(groupName) ) then
      call h5gOpen_f ( locID, trim(groupName), itemID, status )
      ! call outputNamedValue( 'groupName ', groupName )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open group' // trim(groupName) // &
        & ' while looking for all its attr names' )
    elseif(present(DSName) ) then
      call h5dOpen_f ( locID, trim(DSname), itemID, status )
      ! call outputNamedValue( 'DSname ', DSname )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open DS ' // trim(DSname) // &
        & ' while looking for all its attr names' )
    else
      itemID = locID
    endif
    Names = ' '
    ! call output( 'About to read num of attributes', advance='yes' )
    call h5aget_num_attrs_f( itemID, num, status )
    if ( status /= 0 ) go to 9
    ! call outputNamedValue( 'num ', num )
    do i=1, num
      call h5aopen_idx_f( itemid, i-1, attr_id, status )
      ! print *, 'attr_id, status ', attr_id, status
      if ( status == -1 .or. attr_id < 1 ) cycle
      call h5aget_name_f( attr_id, namelength, name, status )
      ! print *, 'name ', trim(name)
      if ( status == -1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to determine attr name' )
      ! print *, 'len(name) ', len_trim(name), status
      tempnames = catLists( names, name )
      names = tempnames
      call h5aclose_f( attr_id, status )
    enddo
    if ( present(groupName) ) then
      call h5gClose_f ( itemID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close group' // trim(groupName) // &
        & ' while looking for all its attr names' )
    elseif(present(DSName) ) then
      call h5dClose_f ( itemID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close DS' // trim(DSname) // &
        & ' while looking for all its attr names' )
    endif
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for all attr names' )
  9 call trace_end ( 'GetAllHDF5AttrNames', cond=.false. )
  end subroutine GetAllHDF5AttrNames

  ! -----------------------------------  GetAllHDF5DSNames_fileID  -----
  subroutine GetAllHDF5DSNames_fileID ( FileID, gname, DSNames, andSlash )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    integer, intent(in) :: FILEID          ! fileID
    character (len=*), intent(in) :: GNAME ! Name of group; e.g. '/'
    character (len=*), intent(out) :: DSNames ! Names of DS in file (,-separated)
    logical, optional, intent(in) :: andSlash ! Keep leading '/' if TRUE

    ! Local variables
    integer :: i                        ! loop counter
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: omitSlash
    integer :: STATUS                   ! Flag
    type(MLSDataInfo_T) :: dataset_info

    ! Executable code
    call trace_begin ( me, 'GetAllHDF5DSNames_fileID', cond=.false. )
    omitSlash = .true.
    if ( present(andSlash) ) omitSlash = .not. andSlash
    ! Initializing values returned if there was trouble
    DSNames = ''
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting DSNames' )
! The structure, dataset_info, is initialized below.
!
    nullify(dataset_info%name)
    call allocate_test ( dataset_info%name, MAXNDSNAMES, 'dataset_info%name', &
      & moduleName )
    dataset_info%name = ''
    dataset_info%number_of_entries = 0
    call Query_MLSData ( fileid, trim(gname), dataset_info )
    if ( dataset_info%number_of_entries < 0 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Unable to get DSNames' )
    else if ( dataset_info%number_of_entries == 1 ) then
      if ( omitSlash .and. dataset_info%name(1)(1:1) == '/' ) &
        & dataset_info%name(1) = dataset_info%name(1)(2:)
      DSNames = dataset_info%name(1)
    else if ( dataset_info%number_of_entries > 1 ) then
      ! DSNames = dataset_info%name(1)
      do i = 1, dataset_info%number_of_entries
        if ( omitSlash .and. dataset_info%name(i)(1:1) == '/' ) &
          & dataset_info%name(i) = dataset_info%name(i)(2:)
        DSNames = catLists(trim(DSNames), dataset_info%name(i))
      end do
    end if
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting DSNames' )
    call deallocate_test ( dataset_info%name, 'dataset_info%name', moduleName )
    call trace_end ( cond=.false. )
  end subroutine GetAllHDF5DSNames_fileID

  ! ---------------------------------  GetAllHDF5DSNames_filename  -----
  subroutine GetAllHDF5DSNames_filename ( FileName, gname, DSNames, andSlash )
    character (len=*), intent(in) :: FILENAME       ! file name
    character (len=*), intent(in) :: GNAME ! Name of group; e.g. '/'
    character (len=*), intent(out) :: DSNames ! Names of DS in file (,-separated)
    logical, optional, intent(in) :: andSlash ! Keep leading '/' if TRUE

    ! Local variables
    integer :: fileID
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'GetAllHDF5DSNames_filename', cond=.false. )
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting DSNames ' // trim(filename) )
    call h5fopen_f ( trim(filename), H5F_ACC_RDONLY_F, fileID, status )
    if ( status == 0 ) then
      call mls_eSet_auto ( 1, status )
      call GetAllHDF5DSNames_fileID ( fileID, gname, DSNames, andSlash )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open file for getting DSNames ' // trim(filename) )
    end if
    call mls_eSet_auto ( 0, status )
    call h5fclose_f ( fileID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close file after getting DSNames ' // trim(filename) )
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting DSNames ' // trim(filename) )
    call trace_end ( cond=.false. )
  end subroutine GetAllHDF5DSNames_filename

  ! ---------------------------------  GetAllHDF5DSNames_MLSFile  -----
  subroutine GetAllHDF5DSNames_MLSFile ( MLSFile, DSNames, andSlash )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(out) :: DSNames ! Names of DS in file (,-separated)
    logical, optional, intent(in) :: andSlash ! Keep leading '/' if TRUE

    ! Local variables
    character, parameter :: gname = '/'
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'GetAllHDF5DSNames_MLSFile', cond=.false. )
    if ( .not. MLSFile%stillOpen ) then
      call GetAllHDF5DSNames( MLSFile%Name, gname, DSNames, andSlash )
    elseif( MLSFile%FileId%grp_id > 0 ) then
      call GetAllHDF5DSNames( MLSFile%FileId%grp_id, gname, DSNames, andSlash )
    elseif( MLSFile%FileId%f_id > 0 ) then
      call GetAllHDF5DSNames( MLSFile%FileId%f_id, gname, DSNames, andSlash )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'all fields of MLSFile%Fileid are 0', MLSFile=MLSFile )
    endif
    call trace_end ( cond=.false. )
  end subroutine GetAllHDF5DSNames_MLSFile

  ! --------------------------------------  MakeHDF5Attribute_dbl  -----
  subroutine MakeHDF5Attribute_dbl ( itemID, name, value, &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID         ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    double precision, intent(in) :: VALUE ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: DSID                     ! ID for dataspace
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_dbl', cond=.false. )
    ! (Maybe) create the attribute
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_DOUBLE, dsID, attrID &
      &, skip_if_already_there ) ) then
      ! Write it
      call h5aWrite_f ( attrID, H5T_NATIVE_DOUBLE, value, ones, status )
      call finishMakeAttrib ( name, status, attrID, dsID )
    end if
    call trace_end ( cond=.false. )
  end subroutine MakeHDF5Attribute_dbl

  ! -------------------------------------  MakeHDF5Attribute_sngl  -----
  subroutine MakeHDF5Attribute_sngl ( itemID, name, value, &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID         ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    real, intent(in) :: VALUE             ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: DSID                     ! ID for dataspace
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_sngl', cond=.false. )
    ! (Maybe) create the attribute
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_REAL, dsID, attrID, &
      & skip_if_already_there ) ) then
      ! Write
      call h5aWrite_f ( attrID, H5T_NATIVE_REAL, value, ones, status )
      call finishMakeAttrib ( name, status, attrID, dsID )
    end if
    call trace_end ( cond=.false. )
  end subroutine MakeHDF5Attribute_sngl

  ! --------------------------------------  MakeHDF5Attribute_int  -----
  subroutine MakeHDF5Attribute_int ( itemID, name, value, &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID         ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(in) :: VALUE          ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: DSID                     ! ID for dataspace
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_int', cond=.false. )
    ! (Maybe) create the attribute
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_INTEGER, dsID, attrID, &
      & skip_if_already_there ) ) then
      ! Write
      call h5aWrite_f ( attrID, H5T_NATIVE_INTEGER, value, ones, status )
      call finishMakeAttrib ( name, status, attrID, dsID )
    end if
    call trace_end ( cond=.false. )
  end subroutine MakeHDF5Attribute_int

  ! ----------------------------------  MakeHDF5Attribute_logical  -----
  subroutine MakeHDF5Attribute_logical ( itemID, name, value, &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(in) :: VALUE        ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: IVALUE                   ! Value as integer
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_logical', cond=.false. )
    iValue = 0
    if ( value ) iValue = 1
    call MakeHDF5Attribute ( itemID, name, iValue, skip_if_already_there )
    call trace_end ( cond=.false. )
  end subroutine MakeHDF5Attribute_logical

  ! ------------------------------  MakeHDF5Attribute_logicalarr1  -----
  subroutine MakeHDF5Attribute_logicalarr1 ( itemID, name, value , &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID         ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(in) :: VALUE(:)       ! The attribute array itself
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: IValue(size(value))      ! 1 for true, 0 for false
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_logicalarr1', cond=.false. )
    where ( value )
      iValue = 1
    elsewhere
      iValue = 0
    end where
    call MakeHDF5Attribute ( itemID, name, iValue, skip_if_already_there )
    call trace_end ( cond=.false. )
  end subroutine MakeHDF5Attribute_logicalarr1

  ! -----------------------------------  MakeHDF5Attribute_string  -----
  subroutine MakeHDF5Attribute_string ( itemID, name, value , &
   & skip_if_already_there, DONT_TRIM )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    character (len=*), intent(in) :: VALUE ! Value of attribute
    logical, intent(in), optional :: SKIP_IF_ALREADY_THERE
    logical, intent(in), optional :: DONT_TRIM

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: DSID                     ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGTYPE               ! Type for string
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_dont_trim
    logical :: my_skip
    logical :: is_present
    logical, parameter :: NEVERDELETE = .false. ! .true.

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_string', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    my_dont_trim = .false.
    if ( present(DONT_TRIM) ) my_dont_trim=DONT_TRIM
    is_present = IsHDF5AttributePresent_in_DSID(itemID, name)
    if ( my_skip .and. is_present ) go to 9
    ! Setup
    ! Create a data type for this string
    call h5tcopy_f ( H5T_NATIVE_CHARACTER, stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype ' // trim(name) )
    if ( my_dont_trim) then
      call h5tset_size_f(stringtype, int(max(len(value), 1), size_t), status )
    else
      call h5tset_size_f(stringtype, int(max(len_trim(value), 1), size_t), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype ' // trim(name) )
    ! Create dataspace and attribute
    !call h5sCreate_F ( h5s_scalar_f, dsID, status )
    !if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
    !  & 'Unable to create dataspace for attribute ' // trim(name) )
    if ( is_present .and. .not. NEVERDELETE ) then
      ! print *, 'Deleting ' // trim(name)
      call h5adelete_f(itemID, trim(name), status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to delete ' )
    end if
    call h5sCreate_F ( h5s_scalar_f, dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for attribute ' // trim(name) )
    ! print *, 'itemID: ', itemID
    ! print *, 'stringtype: ', stringtype
    ! print *, 'dsID: ', dsID
    ! print *, 'name: ', trim(name)
    ! print *, 'was there: ', is_present
    ! print *, 'value: ', trim(value)
    if ( .not. ( is_present .and. NEVERDELETE ) ) then
      call h5aCreate_f ( itemID, trim(name), stringtype, dsID, attrID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create attribute ' // trim(name) )
    else
      call h5aopen_name_f ( itemID, trim(name), attrID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open attribute ' // trim(name) )
    endif
    ! print *, 'attrID: ', attrID
    ! print *, 'status: ', status
    ! Write
    if ( my_dont_trim) then
      call h5aWrite_f ( attrID, stringtype, value, ones, status )
    else
      call h5aWrite_f ( attrID, stringtype, trim_safe(value), ones, status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write attribute ' // trim(name) )
    ! Finish off
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
    if ( .not. is_present .or. .true. ) then
      call h5sClose_f ( dsID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close attribute dataspace ' // trim(name) )
    end if
    call h5tClose_f ( stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close stringtype ' // trim(name) )
  9 call trace_end ( 'MakeHDF5Attribute_string', cond=.false. )
  end subroutine MakeHDF5Attribute_string

  ! -----------------------------------  MakeHDF5Attribute_textFile  -----
  subroutine MakeHDF5Attribute_textFile ( textFile, itemID, name , &
   & skip_if_already_there, maxLineLen )
    use IO_Stuff, only: Read_TextFile
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    character (len=*), intent(in) :: TEXTFILE ! name of textfile
    logical, intent(in), optional :: SKIP_IF_ALREADY_THERE
    integer, optional, intent(in) :: maxLineLen

    ! Local variables
    logical, parameter :: NEVERDELETE = .true.
    integer :: ATTRID                   ! ID for attribute
    character(len=3) :: blockChar
    integer :: DSID                     ! ID for dataspace
    integer :: firstChar, lastChar
    integer :: iblock, nblocks
    logical :: is_present
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: myMaxLineLen
    character(len=len(name)+3) :: newname
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGTYPE               ! Type for string
    logical :: my_skip
    character(LEN=MAXTEXTSIZE)              :: value

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_textFile', cond=.false. )
    myMaxLineLen = DFLTMAXLINELENGTH
    if ( present(maxLineLen) ) myMaxLineLen = maxLineLen
    ! Try to read the textfile
    value = ' '
    call read_textFile( trim(textFile), value, myMaxLineLen )
    status = 0
    if ( status /= 0 ) then
      call MLSMessage(MLSMSG_Warning, ModuleName, &
        & 'Unable to write attribute--failed to read textfile' )
      goto 9
    endif
    ! Unfortunately, a lot of null characters sneak into this
    value = Replace( value, char(0), char(32) ) ! Replace null with space
    
    ! Be careful lest the attribute is too large
    if ( len_trim(value) > MAXATTRIBUTESIZE ) then
      nblocks = 1 + ( len_trim(value) - 1 ) / MAXATTRIBUTESIZE
      lastChar = 0
      do iblock=1, nblocks
        firstChar = lastChar + 1
        lastChar  = min( lastChar + MAXATTRIBUTESIZE, len_trim(value) )
        write(blockChar, '(i1)' ) iblock
        if ( iblock > 9 ) write(blockChar, '(i2)' ) iblock
        newName = trim(name) // adjustl(blockChar)
        call MakeHDF5Attribute( itemID, newName, value(firstChar:lastChar) )
      enddo
      goto 9
    endif

    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    is_present = IsHDF5AttributePresent_in_DSID(itemID, name)
    if ( my_skip .and. is_present ) go to 9
    ! Setup
    ! Create a data type for this string
    call h5tcopy_f ( H5T_NATIVE_CHARACTER, stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype ' // trim(name) )
    call h5tset_size_f(stringtype, int(max(len(value), 1), size_t), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype ' // trim(name) )
    ! Create dataspace and attribute
    if ( is_present .and. .not. NEVERDELETE ) then
      print *, 'Deleting ' // trim(name)
      call h5adelete_f(itemID, trim(name), status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to delete ' )
    end if
    call h5sCreate_F ( h5s_scalar_f, dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for attribute ' // trim(name) )
    ! print *, 'itemID: ', itemID
    ! print *, 'stringtype: ', stringtype
    ! print *, 'dsID: ', dsID
    ! print *, 'name: ', trim(name)
    ! print *, 'was there: ', is_present
    ! print *, 'value: ', trim(value)
    if ( .not. ( is_present .and. NEVERDELETE ) ) then
      call h5aCreate_f ( itemID, trim(name), stringtype, dsID, attrID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create attribute ' // trim(name) )
    else
      call h5aopen_name_f ( itemID, trim(name), attrID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open attribute ' // trim(name) )
    endif
    ! print *, 'attrID: ', attrID
    ! print *, 'status: ', status
    ! Write
    call h5aWrite_f ( attrID, stringtype, trim(value), ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write attribute ' // trim(name) )
    ! Finish off
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
    if ( .not. is_present .or. .true. ) then
      call h5sClose_f ( dsID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close attribute dataspace ' // trim(name) )
    end if
    call h5tClose_f ( stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close stringtype ' // trim(name) )
  9 call trace_end ( 'MakeHDF5Attribute_textFile', cond=.false. )
  end subroutine MakeHDF5Attribute_textFile

  ! -------------------------------  MakeHDF5Attribute_stringarr1  -----
  subroutine MakeHDF5Attribute_stringarr1 ( itemID, name, value , &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID         ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    character (len=*), dimension(:), intent(in) :: VALUE ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: DSID                     ! ID for dataspace
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGTYPE               ! Type for string
    logical :: my_skip
    integer(kind=hsize_t), dimension(1) :: SHP        ! Shape

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_stringarr1', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip = skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_DSID(itemID, name) ) go to 9
    end if
    ! Setup
    shp = shape(value)
    ! Create a data type for this string
    call h5tcopy_f ( H5T_NATIVE_CHARACTER, stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype for array' // trim(name) )
    call h5tset_size_f ( stringtype, int(len(value(1)), size_t), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype ' // trim(name) )
    ! Create dataspace and attribute
    ! call h5sCreate_F ( h5s_simple_f, dsID, status )
    call h5sCreate_simple_f ( 1, int(shp, hSize_T), dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for attribute ' // trim(name) )
    call h5aCreate_f ( itemID, trim(name), stringtype, dsID, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create attribute ' // trim(name) )
    ! Write
    call h5aWrite_f ( attrID, stringtype, value, &
      & ones, status )
    call finishMakeAttrib ( name, status, attrID, dsID, stringType )
  9 call trace_end ( cond=.false. )
  end subroutine MakeHDF5Attribute_stringarr1

  ! ---------------------------------  MakeHDF5Attribute_snglarr1  -----
  subroutine MakeHDF5Attribute_snglarr1 ( itemID, name, value , &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID         ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    real, intent(in) :: VALUE(:)          ! The attribute array itself
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: spaceID                  ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=hsize_t), dimension(1) :: SHP        ! Shape

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_snglarr1', cond=.false. )
    ! (Maybe) create the attribute
    shp = shape(value)
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_REAL, spaceID, &
      & attrID, skip_if_already_there, shp ) ) then
      ! Write
      call h5aWrite_f ( attrID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:6) /), hsize_t ), status )
      call finishMakeAttrib ( name, status, attrID, spaceID )
    end if
    call trace_end ( cond=.false. )
  end subroutine MakeHDF5Attribute_snglarr1

  ! ----------------------------------  MakeHDF5Attribute_dblarr1  -----
  subroutine MakeHDF5Attribute_dblarr1 ( itemID, name, value , &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID            ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME    ! Name of attribute
    double precision, intent(in) :: VALUE(:) ! The attribute array itself
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: spaceID                  ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=hsize_t), dimension(1) :: SHP        ! Shape

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_dblarr1', cond=.false. )
    ! (Maybe) create the attribute
    shp = shape(value)
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_DOUBLE, spaceID, &
      & attrID, skip_if_already_there, shp ) ) then
      ! Write
      call h5aWrite_f ( attrID, H5T_NATIVE_DOUBLE, value, &
        & int ( (/ shp, ones(1:6) /), hsize_t ), status )
      call finishMakeAttrib ( name, status, attrID, spaceID )
    end if
    call trace_end ( cond=.false. )
  end subroutine MakeHDF5Attribute_dblarr1

  ! ----------------------------------  MakeHDF5Attribute_intarr1  -----
  subroutine MakeHDF5Attribute_intarr1 ( itemID, name, value , &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID         ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(in) :: VALUE(:)       ! The attribute array itself
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: spaceID                  ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=hsize_t), dimension(1) :: SHP        ! Shape

    ! Executable code
    call trace_begin ( me, 'MakeHDF5Attribute_intarr1', cond=.false. )
    ! (Maybe) create the attribute
    shp = shape(value)
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_INTEGER, spaceID, &
      & attrID, skip_if_already_there, shp ) ) then
      ! Write
      call h5aWrite_f ( attrID, H5T_NATIVE_INTEGER, value, &
        & int ( (/ shp, ones(1:6) /), hsize_t ), status )
      call finishMakeAttrib ( name, status, attrID, spaceID )
    end if
    call trace_end ( cond=.false. )
  end subroutine MakeHDF5Attribute_intarr1

  ! -----------------------------------  MakeHDF5AttributeDSN_int  -----
  subroutine MakeHDF5AttributeDSN_int ( fileID, dataName, attrName, value, &
   & skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    integer, intent(in) :: VALUE        ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: dataID                   ! ID for dataspace
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag from HDF5
    logical :: my_skip

    ! Executable code
    call trace_begin ( me, 'MakeHDF5AttributeDSN_int', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) go to 9
    end if
    dataID = name_to_dataID( fileID, dataName )
    call MakeHDF5Attribute_int ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )

  9 call trace_end ( cond=.false. )
  end subroutine MakeHDF5AttributeDSN_int

  ! -------------------------------  MakeHDF5AttributeDSN_logical  -----
  ! Sorry--this could not be made part of MakeHDF5Attribute
  ! it conflicts with MakeHDF5Attribute_string generic
  ! all because of the skip_if_already_there optional argument
  subroutine MakeHDF5AttributeDSN_logical ( fileID, dataName, attrName, value, &
   & skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    logical, intent(in) :: VALUE        ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: dataID                   ! ID for data
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag from HDF5
    logical :: my_skip

    ! Executable code
    call trace_begin ( me, 'MakeHDF5AttributeDSN_logical', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) go to 9
    end if
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute_logical ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )

  9 call trace_end ( cond=.false. )
  end subroutine MakeHDF5AttributeDSN_logical

  ! --------------------------------  MakeHDF5AttributeDSN_string  -----
  subroutine MakeHDF5AttributeDSN_string ( fileID, dataName, attrName, value, &
   & skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    character (len=*), intent(in) :: VALUE ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: dataID                   ! ID for data
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_skip
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'MakeHDF5AttributeDSN_string', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) go to 9
    end if
    dataID = name_to_dataID( fileID, dataName)
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_DSID( dataID, attrName ) ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Inconsistency between file/dataName for attribute ' // trim(attrName) )
      end if
    end if
    call MakeHDF5Attribute_string ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )

  9 call trace_end ( cond=.false. )
  end subroutine MakeHDF5AttributeDSN_string

  ! -------------------------------  MakeHDF5AttributeDSN_st_arr1  -----
  subroutine MakeHDF5AttributeDSN_st_arr1 ( fileID, dataName, attrName, value, &
   & skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    character (len=*), dimension(:), intent(in) :: VALUE ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: dataID                   ! ID for data
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_skip
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'MakeHDF5AttributeDSN_st_arr1', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) go to 9
    end if
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute ( dataID, attrName, value )
    call h5dclose_f(dataID, status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )

  9 call trace_end ( cond=.false. )
  end subroutine MakeHDF5AttributeDSN_st_arr1

  ! ------------------------------  MakeHDF5AttributeDSN_single  -----
  subroutine MakeHDF5AttributeDSN_single ( fileID, dataName, attrName, value, &
   & skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    real, intent(in) :: VALUE        ! The attribute value itself
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: dataID                   ! ID for data
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_skip
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'MakeHDF5AttributeDSN_single', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) go to 9
    end if
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )
 9  call trace_end ( cond=.false. )

  end subroutine MakeHDF5AttributeDSN_single

  ! ------------------------------  MakeHDF5AttributeDSN_snglarr1  -----
  subroutine MakeHDF5AttributeDSN_snglarr1 ( fileID, dataName, attrName, value, &
   & skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    real, intent(in) :: VALUE(:)        ! The attribute array itself
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: dataID                   ! ID for data
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_skip
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'MakeHDF5AttributeDSN_snglarr1', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) go to 9
    end if
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )
 9  call trace_end ( cond=.false. )

  end subroutine MakeHDF5AttributeDSN_snglarr1

  ! --------------------------------  MakeHDF5AttributeDSN_double  -----
  subroutine MakeHDF5AttributeDSN_double ( fileID, dataName, attrName, value, &
   & skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    double precision, intent(in) :: VALUE  ! The attribute value itself
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: dataID                   ! ID for data
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_skip
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'MakeHDF5AttributeDSN_double', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) go to 9
    end if
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )
  9 call trace_end ( cond=.false. )

  end subroutine MakeHDF5AttributeDSN_double

  ! -------------------------------  MakeHDF5AttributeDSN_dblarr1  -----
  subroutine MakeHDF5AttributeDSN_dblarr1 ( fileID, dataName, attrName, value, &
   & skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    double precision, intent(in) :: VALUE(:)  ! The attribute array itself
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: dataID                   ! ID for data
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_skip
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'MakeHDF5AttributeDSN_dblarr1', cond=.false. )
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) go to 9
    end if
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )
  9 call trace_end ( cond=.false. )

  end subroutine MakeHDF5AttributeDSN_dblarr1

  ! --------------------------------------  MakeNestedGroups  -----
  ! Make a nested set of groups, beginning with groupNames(1)
  ! and ending with groupNames(n)
  ! If innermostID present, it will contain grpID(n)
  ! Otherwise, we'll close all the grpIDs
  subroutine MakeNestedGroups ( locID, groupNames, innermostID )
    integer, intent(in)                         :: locID  ! Where to base groups at
    character (len=*), intent(in), dimension(:) :: groupNames
    integer, intent(out), optional              :: innermostID ! grpID(n)

    ! Local variables
    integer :: containerID              ! ID for container
    integer :: grpID                    ! ID for new group
    integer :: i
    integer, dimension(size(groupNames)) :: IDs
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'MakeNestedGroups', cond=.false. )
    containerID = locID
    ! call dump( groupNames, 'groupNames' )
    ! call outputNamedValue( 'size(groupNames)', size(groupNames) )
    do i=1, size(groupNames)
      ! call outputNamedValue( 'groupName(i)', trim(groupnames(i)) )
      if ( len_trim(groupnames(i)) < 1 ) then
        iDs(i) = -1
        cycle
      elseif ( .not. IsHDF5GroupPresent ( containerID, trim(groupnames(i))) ) then
        !  Must create this group
        call h5gCreate_f ( containerID, trim(groupnames(i)), grpID, status )
      else
        ! Must open this group
        call h5gopen_f ( containerID, trim(groupnames(i)), grpID, status )
      endif
      containerID = grpID
      iDs(i) = grpID
    enddo
    if ( present(innermostID) ) then
      innermostID = containerID
    else
      do i=1, size(groupNames)
        if ( IDs(i) > 0 ) call h5gclose_f ( IDs(i), status )
      enddo
    endif
    call trace_end ( cond=.false. )
  end subroutine MakeNestedGroups

  ! ---------------------------------------  GetHDF5Attribute_int  -----
  subroutine GetHDF5Attribute_int ( MLSFile, name, value )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(out) :: VALUE         ! Result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attribute_int', cond=.false. )
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%sd_ID, name, attrID, status )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%grp_ID, name, attrID, status )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'no valid sd_id or grp_id to get attribute', MLSFile=MLSFile)
      go to 9
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name), &
      & MLSFile=MLSFile )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read an array into our one value.
    call h5aread_f ( attrID, H5T_NATIVE_INTEGER, value, ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name), &
      & MLSFile=MLSFile )
  9 call trace_end ( cond=.false. )
  end subroutine GetHDF5Attribute_int

  ! -----------------------------------  GetHDF5Attribute_intarr1  -----
  subroutine GetHDF5Attribute_intarr1 ( MLSFile, name, value )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(out) :: VALUE(:)      ! Result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: Me = -1                  ! String index for trace cacheing
    integer(kind=hsize_t), dimension(1) :: SHP        ! Shape
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attribute_intarr1', cond=.false. )
    shp = shape(value)
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%sd_ID, name, attrID, status )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%grp_ID, name, attrID, status )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'no valid sd_id or grp_id to get attribute', MLSFile=MLSFile)
      go to 9
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name), &
      & MLSFile=MLSFile )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read an array into our one value.
    call h5aread_f ( attrID, H5T_NATIVE_INTEGER, value, &
      & int ( (/ shp, ones(1:6) /), hsize_t ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name), &
      & MLSFile=MLSFile )
  9 call trace_end ( cond=.false. )
  end subroutine GetHDF5Attribute_intarr1

  ! ------------------------------------  GetHDF5Attribute_string  -----
  subroutine GetHDF5Attribute_string ( MLSFile, name, value )
    use HDF5, only: H5TIs_Variable_Str_F, H5SGet_Simple_Extent_NPoints_F
    use, intrinsic :: Iso_C_Binding
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME   ! Name of attribute
    character (len=*), intent(out) :: VALUE ! Result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    ! logical, parameter :: DEEBUG = .true.
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: VarLenStr                ! Is the string length variable?
    integer :: dspace_id
    integer(kind=Size_t) :: n_pts
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=Size_t) :: STRINGSIZE               ! String size
    integer :: STRINGTYPE               ! String type

    ! Executable code
    value = ''
    call trace_begin ( me, 'GetHDF5Attribute_string', cond=.false. )
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%sd_ID, name, attrID, status )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%grp_ID, name, attrID, status )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'no valid sd_id or grp_id to get attribute', MLSFile=MLSFile)
      go to 9
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5aGet_type_f ( attrID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5tis_variable_str_f ( stringType, VarLenStr, Status )
    if ( DEEBug ) then
      call outputnamedvalue ( 'string length variable?', VarLenStr )
      call h5aget_space_f ( attrID, dspace_id, status )
      call outputnamedvalue ( 'dspace_id', dspace_id )
      call h5sget_simple_extent_npoints_f ( dspace_id, n_pts, status )
      call outputnamedvalue ( 'n_pts', int(n_pts) )
    endif
    if ( stringSize > len(value) ) then
      call outputnamedValue( 'stringSize', int(stringSize) )
      call outputnamedValue( 'len(value)', len(value) )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Value too long to fit in space given for attribute ' // trim(name), &
        & MLSFile=MLSFile )
    elseif( DEEBUG ) then
      call outputnamedValue( 'stringSize', int(stringSize) )
      call outputnamedValue( 'len(value)', len(value) )
    endif

    ! Now actually read the data!
    if ( VarLenStr ) then
      call VariableStringLength
    else
      call FixedStringLength
    endif
    if( DEEBUG ) print *, trim(value)
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5tClose_f ( stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close string type for attribute ' // trim(name), &
      & MLSFile=MLSFile )
  9 call trace_end ( cond=.false. )
  contains
    subroutine FixedStringLength
      call h5aread_f ( attrID, stringType, value, ones, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to read attribute ' // trim(name), &
        & MLSFile=MLSFile )
      ! In case there are any nulls, we will replace them with blanks
      status = FindFirst( value, achar(0) )
      if ( status > 0 ) value( status: ) = ' '
    end subroutine FixedStringLength
    subroutine VariableStringLength
      integer, parameter :: max_len = 132
      character(len=max_len, kind=c_char), pointer :: cdata ! A pointer to a Fortran string
      integer(hsize_t), dimension(1:1) :: dims = (/1/)
      type(c_ptr) :: f_ptr
      integer :: length
      integer :: i
      integer(hsize_t), dimension(1:1) :: maxdims
      integer :: mem_space_id
      type(c_ptr), dimension(:), allocatable, target :: rdata ! Read buffer
      call h5aget_space_f ( attrID, mem_space_id, status )
      call h5sget_simple_extent_dims_f ( mem_space_id, dims, maxdims, status )

      allocate(rdata(1:dims(1)))
      !
      ! Read the data.
      !
      f_ptr = C_LOC( rdata(1) )
      CALL H5aread_f( attrID, H5T_STRING, f_ptr, status )
      !
      ! Output the variable-length data to the screen.
      !
      do i = 1, dims(1)
         call c_f_pointer( rdata(i), cdata )
         length = 0
         do while( cdata(length+1:length+1) .ne. c_null_char )
            length = length + 1
            value(length:length) = cdata(length:length)
         enddo
         if( DEEBUG ) write( *,'(a)' ) cdata( 1:length )
      end do
      !
      ! Close and release resources.
      !
      call h5sclose_f( mem_space_id, status )

    end subroutine VariableStringLength
  end subroutine GetHDF5Attribute_string

  ! --------------------------------  GetHDF5Attribute_stringarr1  -----
  subroutine GetHDF5Attribute_stringarr1 ( MLSFile, name, value )
    use HDF5, only: H5TIs_Variable_Str_F
    use, intrinsic :: Iso_C_Binding
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of attribute
    character (len=*), intent(out) :: VALUE(:) ! Result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    logical :: VarLenStr                ! Is the string length variable?
    integer :: Me = -1                  ! String index for trace cacheing
    integer(kind=hsize_t), dimension(1) :: SHP        ! Shape
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=Size_t) :: STRINGSIZE               ! String size
    integer :: STRINGTYPE               ! String type

    ! Executable code
    value = ''
    call trace_begin ( me, 'GetHDF5Attribute_stringarr1', cond=.false. )
    shp = shape(value)
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%sd_ID, name, attrID, status )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%grp_ID, name, attrID, status )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'no valid sd_id or grp_id to get attribute', MLSFile=MLSFile)
      go to 9
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5aGet_type_f ( attrID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for 1-d string attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for 1-d string attribute ' // trim(name), &
      & MLSFile=MLSFile )
    if ( stringSize > len(value) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for attribute ' // trim(name), &
      & MLSFile=MLSFile )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read an array into our one value.
    call h5tis_variable_str_f ( stringType, VarLenStr, Status )

    ! Now actually read the data!
    if ( VarLenStr ) then
      call VariableStringLength
    else
      call FixedStringLength
    endif
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5tClose_f ( stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close string type for attribute ' // trim(name), &
      & MLSFile=MLSFile )
  9 call trace_end ( cond=.false. )
  contains
    subroutine FixedStringLength
      call h5aread_f ( attrID, stringType, value, &
        & int ( (/ shp, ones(1:6) /), hsize_t ), status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to read attribute ' // trim(name), &
        & MLSFile=MLSFile )
    end subroutine FixedStringLength
    subroutine VariableStringLength
      integer, parameter :: max_len = 132
      character(len=max_len, kind=c_char), pointer :: cdata ! A pointer to a Fortran string
      integer(hsize_t), dimension(1:1) :: dims = (/1/)
      type(c_ptr) :: f_ptr
      integer :: length
      integer :: i
      integer(hsize_t), dimension(1:1) :: maxdims
      integer :: mem_space_id
      type(c_ptr), dimension(:), allocatable, target :: rdata ! Read buffer
      call h5aget_space_f ( attrID, mem_space_id, status )
      call h5sget_simple_extent_dims_f ( mem_space_id, dims, maxdims, status )

      allocate(rdata(1:dims(1)))
      !
      ! Read the data.
      !
      f_ptr = C_LOC( rdata(1) )
      CALL H5aread_f( attrID, H5T_STRING, f_ptr, status )
      !
      ! Output the variable-length data to the screen.
      !
      do i = 1, dims(1)
         call c_f_pointer( rdata(i), cdata )
         length = 0
         do while( cdata(length+1:length+1) .ne. c_null_char )
            length = length + 1
            value(length) = cdata(length:length)
         enddo
         if( DEEBUG ) write( *,'(a)' ) cdata( 1:length )
      end do
      !
      ! Close and release resources.
      !
      call h5sclose_f( mem_space_id, status )

    end subroutine VariableStringLength
  end subroutine GetHDF5Attribute_stringarr1

  ! -----------------------------------  GetHDF5Attribute_logical  -----
  subroutine GetHDF5Attribute_logical ( MLSFile, name, value )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(out) :: VALUE         ! Value of attribute

    ! Local variables
    integer :: IVALUE                     ! Value as integer
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attribute_logical', cond=.false. )
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call GetHDF5Attribute ( MLSFile%fileID%sd_ID, name, iValue )
      if ( DEEBUG ) call outputNamedValue( 'MLSFile%fileID%sd_ID', MLSFile%fileID%sd_ID )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call GetHDF5Attribute ( MLSFile%fileID%grp_ID, name, iValue )
      if ( DEEBUG ) call outputNamedValue( 'MLSFile%fileID%grp_ID', MLSFile%fileID%sd_ID )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'no valid sd_id or grp_id to get attribute', MLSFile=MLSFile)
      go to 9
    endif

    if ( DEEBUG ) call outputNamedValue( 'ivalue', ivalue )
    value = ( iValue == 1 )
  9 call trace_end ( cond=.false. )
  end subroutine GetHDF5Attribute_logical

  ! -------------------------------  GetHDF5Attribute_logicalarr1  -----
  subroutine GetHDF5Attribute_logicalarr1 ( MLSFile, name, value )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(out) :: VALUE(:)      ! Value of attribute

    ! Local variables
    integer :: IVALUE(size(VALUE,1))      ! Value as integer
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attribute_logicalarr1', cond=.false. )
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call GetHDF5Attribute ( MLSFile%fileID%sd_ID, name, iValue )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call GetHDF5Attribute ( MLSFile%fileID%grp_ID, name, iValue )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'no valid sd_id or grp_id to get attribute', MLSFile=MLSFile)
      go to 9
    endif
    value = ( iValue == 1 )
  9 call trace_end ( cond=.false. )
  end subroutine GetHDF5Attribute_logicalarr1

  ! ----------------------------------  GetHDF5Attribute_snglarr1  -----
  subroutine GetHDF5Attribute_snglarr1 ( MLSFile, name, value )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of attribute
    real, intent(out) :: VALUE(:)       ! The attribute array result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: Me = -1                  ! String index for trace cacheing
    integer(kind=hsize_t), dimension(1) :: SHP        ! Shape
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attribute_snglarr1', cond=.false. )
    shp = shape(value)
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%sd_ID, name, attrID, status )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%grp_ID, name, attrID, status )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'no valid sd_id or grp_id to get attribute', MLSFile=MLSFile)
      go to 9
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name), &
      & MLSFile=MLSFile )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_REAL, value, &
      & int ( (/ shp, ones(1:6) /), hsize_t ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read 1d attribute array ' // trim(name), &
      & MLSFile=MLSFile )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close 1d attribute array  ' // trim(name), &
      & MLSFile=MLSFile )
  9 call trace_end ( cond=.false. )
  end subroutine GetHDF5Attribute_snglarr1

  ! --------------------------------------  GetHDF5Attribute_sngl  -----
  subroutine GetHDF5Attribute_sngl ( MLSFile, name, value )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of attribute
    real, intent(out) :: VALUE          ! The attribute result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attribute_sngl', cond=.false. )
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%sd_ID, name, attrID, status )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%grp_ID, name, attrID, status )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'no valid sd_id or grp_id to get attribute', MLSFile=MLSFile)
      go to 9
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name), &
      & MLSFile=MLSFile )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_REAL, value, &
      & ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read sngl attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close sngl attribute  ' // trim(name) )
  9 call trace_end ( cond=.false. )
  end subroutine GetHDF5Attribute_sngl

  ! ---------------------------------------  GetHDF5Attribute_dbl  -----
  subroutine GetHDF5Attribute_dbl ( MLSFile, name, value )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME  ! Name of attribute
    double precision, intent(out) :: VALUE ! The attribute result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attribute_dbl', cond=.false. )
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%sd_ID, name, attrID, status )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%grp_ID, name, attrID, status )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'no valid sd_id or grp_id to get attribute', MLSFile=MLSFile)
      go to 9
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name), &
      & MLSFile=MLSFile )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_DOUBLE, value, ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dble attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dble attribute  ' // trim(name), &
      & MLSFile=MLSFile )
  9 call trace_end ( cond=.false. )
  end subroutine GetHDF5Attribute_dbl

  ! -----------------------------------  GetHDF5Attribute_dblarr1  -----
  subroutine GetHDF5Attribute_dblarr1 ( MLSFile, name, value )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME     ! Name of attribute
    double precision, intent(out) :: VALUE(:) ! The attribute result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: Me = -1                  ! String index for trace cacheing
    integer(kind=hsize_t), dimension(1) :: SHP        ! Shape
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attribute_dblarr1', cond=.false. )
    shp = shape(value)
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%sd_ID, name, attrID, status )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%grp_ID, name, attrID, status )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'no valid sd_id or grp_id to get attribute', MLSFile=MLSFile)
      go to 9
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name), &
      & MLSFile=MLSFile )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_DOUBLE, value, &
      & int ( (/ shp, ones(1:6) /), hsize_t ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dblarr1 attribute ' // trim(name), &
      & MLSFile=MLSFile )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dblarr1 attribute  ' // trim(name), &
      & MLSFile=MLSFile )
  9 call trace_end ( cond=.false. )
  end subroutine GetHDF5Attribute_dblarr1

  ! ---------------------------------------  GetHDF5Attr_ID_int  -----
  subroutine GetHDF5Attr_ID_int ( itemID, name, value )
    integer, intent(in) :: ITEMID         ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(out) :: VALUE         ! Result

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attr_ID_int', cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = itemID
    MLSFile%stillOpen = .true.
    call GetHDF5Attribute ( MLSFile, name, value )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5Attr_ID_int

  ! -----------------------------------  GetHDF5Attr_ID_intarr1  -----
  subroutine GetHDF5Attr_ID_intarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID         ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(out) :: VALUE(:)      ! Result

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attr_ID_intarr1', cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = itemID
    MLSFile%stillOpen = .true.
    call GetHDF5Attribute ( MLSFile, name, value )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5Attr_ID_intarr1

  ! ------------------------------------  GetHDF5Attr_ID_string  -----
  subroutine GetHDF5Attr_ID_string ( itemID, name, value )
    integer, intent(in) :: ITEMID           ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME   ! Name of attribute
    character (len=*), intent(out) :: VALUE ! Result

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attr_ID_string', cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = itemID
    MLSFile%stillOpen = .true.
    call GetHDF5Attribute ( MLSFile, name, value )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5Attr_ID_string

  ! --------------------------------  GetHDF5Attr_ID_stringarr1  -----
  subroutine GetHDF5Attr_ID_stringarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    character (len=*), intent(out) :: VALUE(:) ! Result

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attr_ID_stringarr1', cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = itemID
    MLSFile%stillOpen = .true.
    call GetHDF5Attribute ( MLSFile, name, value )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5Attr_ID_stringarr1

  ! -----------------------------------  GetHDF5Attr_ID_logical  -----
  subroutine GetHDF5Attr_ID_logical ( itemID, name, value )
    integer, intent(in) :: ITEMID         ! Group etc. to get attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(out) :: VALUE         ! Value of attribute

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attr_ID_logical', cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = itemID
    MLSFile%stillOpen = .true.
    call GetHDF5Attribute ( MLSFile, name, value )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5Attr_ID_logical

  ! -------------------------------  GetHDF5Attr_ID_logicalarr1  -----
  subroutine GetHDF5Attr_ID_logicalarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID         ! Group etc. to get attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(out) :: VALUE(:)      ! Value of attribute

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attr_ID_logicalarr1', cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = itemID
    MLSFile%stillOpen = .true.
    call GetHDF5Attribute ( MLSFile, name, value )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5Attr_ID_logicalarr1

  ! ----------------------------------  GetHDF5Attr_ID_snglarr1  -----
  subroutine GetHDF5Attr_ID_snglarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    real, intent(out) :: VALUE(:)       ! The attribute array result

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attr_ID_snglarr1', cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = itemID
    MLSFile%stillOpen = .true.
    call GetHDF5Attribute ( MLSFile, name, value )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5Attr_ID_snglarr1

  ! --------------------------------------  GetHDF5Attr_ID_sngl  -----
  subroutine GetHDF5Attr_ID_sngl ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    real, intent(out) :: VALUE          ! The attribute result

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attr_ID_sngl', cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = itemID
    MLSFile%stillOpen = .true.
    call GetHDF5Attribute ( MLSFile, name, value )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5Attr_ID_sngl

  ! ---------------------------------------  GetHDF5Attr_ID_dbl  -----
  subroutine GetHDF5Attr_ID_dbl ( itemID, name, value )
    integer, intent(in) :: ITEMID          ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME  ! Name of attribute
    double precision, intent(out) :: VALUE ! The attribute result

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attr_ID_dbl', cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = itemID
    MLSFile%stillOpen = .true.
    call GetHDF5Attribute ( MLSFile, name, value )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5Attr_ID_dbl

  ! -------------------------------------  GetHDF5Attr_ID_dblarr1  -----
  subroutine GetHDF5Attr_ID_dblarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME     ! Name of attribute
    double precision, intent(out) :: VALUE(:) ! The attribute result

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5Attr_ID_dblarr1', cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = itemID
    MLSFile%stillOpen = .true.
    call GetHDF5Attribute ( MLSFile, name, value )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5Attr_ID_dblarr1

  ! --------------------------------  GetHDF5AttributePtr_intarr1  -----
  subroutine GetHDF5AttributePtr_intarr1 ( itemID, name, value, LowBound )
    ! Allocate Value and read attribute Name into it
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: ITEMID         ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, pointer :: VALUE(:)          ! Result
    integer, intent(in), optional :: LowBound ! of allocated value, else 1

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: LB                       ! 1, else LowBound
    integer :: Me = -1                  ! String index for trace cacheing
    integer(kind=hSize_t), dimension(1) :: MaxShp, SHP  ! Shape
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5AttributePtr_intarr1', cond=.false. )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    call GetHDF5AttrDims ( itemID, trim(name), shp, maxShp )
    call allocate_test ( value, int(lb-1+maxShp(1)), name, moduleName )
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read an array into our one value.
    call h5aread_f ( attrID, H5T_NATIVE_INTEGER, value, &
      & (/ shp, ones(1:6) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5AttributePtr_intarr1

  ! -----------------------------  GetHDF5AttributePtr_stringarr1  -----
  subroutine GetHDF5AttributePtr_stringarr1 ( itemID, name, value, LowBound )
    ! Allocate Value and read attribute Name into it
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: ITEMID          ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME  ! Name of attribute
    character (len=*), pointer :: VALUE(:) ! Result
    integer, intent(in), optional :: LowBound ! of allocated value, else 1

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: LB                       ! 1, else LowBound
    integer :: Me = -1                  ! String index for trace cacheing
    integer(kind=hSize_t), dimension(1) :: MaxShp, SHP ! Shape
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=Size_t) :: STRINGSIZE               ! String size
    integer :: STRINGTYPE               ! String type

    ! Executable code
    call trace_begin ( me, 'GetHDF5AttributePtr_stringarr1', cond=.false. )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    call GetHDF5AttrDims ( itemID, trim(name), shp, maxShp )
    call allocate_test ( value, int(lb-1+maxShp(1)), name, moduleName )
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    call h5aGet_type_f ( attrID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for 1-d string attribute ' // trim(name) )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for 1-d string attribute ' // trim(name) )
    if ( stringSize > len(value) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read an array into our one value.
    call h5aread_f ( attrID, stringType, value, &
      & (/ shp, ones(1:6) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
    call h5tClose_f ( stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close string type for attribute ' // trim(name) )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5AttributePtr_stringarr1

  ! ----------------------------  GetHDF5AttributePtr_logicalarr1  -----
  subroutine GetHDF5AttributePtr_logicalarr1 ( itemID, name, value, LowBound )
    use Allocate_Deallocate, only: Deallocate_Test
    integer, intent(in) :: ITEMID         ! Group etc. to get attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, pointer :: VALUE(:)          ! Value of attribute
    integer, intent(in), optional :: LowBound ! of allocated value, else 1

    ! Local variables
    integer, pointer :: IVALUE(:)         ! Value as integer
    integer :: Me = -1                    ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'GetHDF5AttributePtr_logicalarr1', cond=.false. )
    nullify ( ivalue )
    call GetHDF5AttributePtr ( itemID, name, iValue, lowBound )
    value = ( iValue == 1 )
    call deallocate_test ( ivalue, name, moduleName )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5AttributePtr_logicalarr1

  ! -------------------------------  GetHDF5AttributePtr_snglarr1  -----
  subroutine GetHDF5AttributePtr_snglarr1 ( itemID, name, value, LowBound )
    ! Allocate Value and read attribute Name into it
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: ITEMID         ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    real, pointer :: VALUE(:)             ! The attribute array result
    integer, intent(in), optional :: LowBound ! of allocated value, else 1

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: LB                       ! 1, else LowBound
    integer :: Me = -1                  ! String index for trace cacheing
    integer(kind=hSize_t), dimension(1) :: MaxShp, SHP  ! Shape
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5AttributePtr_snglarr1', cond=.false. )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    call GetHDF5AttrDims ( itemID, trim(name), shp, maxShp )
    call allocate_test ( value, int(lb-1+maxShp(1)), name, moduleName )
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_REAL, value, &
      & (/ shp, int(ones(1:6), hsize_t) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read 1d attribute array ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close 1d attribute array  ' // trim(name) )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5AttributePtr_snglarr1

  ! ---------------------------------  GetHDF5AttributePtr_dblarr1  -----
  subroutine GetHDF5AttributePtr_dblarr1 ( itemID, name, value, LowBound )
    ! Allocate Value and read attribute Name into it
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: ITEMID         ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    double precision, pointer :: VALUE(:) ! The attribute result
    integer, intent(in), optional :: LowBound ! of allocated value, else 1

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: LB                       ! 1, else LowBound
    integer :: Me = -1                  ! String index for trace cacheing
    integer(kind=hSize_t), dimension(1) :: MaxShp, SHP  ! Shape
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'GetHDF5AttributePtr_dblarr1', cond=.false. )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    call GetHDF5AttrDims ( itemID, trim(name), shp, maxShp )
    call allocate_test ( value, int(lb-1+maxShp(1)), name, moduleName )
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_DOUBLE, value, &
      & (/ shp, int(ones(1:6), hsize_t) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dblarr1 attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dblarr1 attribute  ' // trim(name) )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5AttributePtr_dblarr1

  ! -----------------------------------  GetAllHDF5GroupNames  -----
  subroutine GetAllHDF5GroupNames ( FileID, GroupNames )
    integer, intent(in) :: FILEID          ! fileID
    character (len=*), intent(out) :: GroupNames ! Names of Group in file (,-separated)

    ! Local variables
    integer :: h5error
    integer :: i                        ! loop counter
    character(len=32) :: objName
    integer :: objType
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag
    integer :: nsubmembers
    logical, parameter :: DeeBug = .false.
    ! Executable code
    call trace_begin ( me, 'GetAllHDF5GroupNames', cond=.false. )
    GroupNames = ' '
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for all grp names' )
    call h5gn_members_f( FileID, '/', nsubmembers, h5error )
    if ( DeeBug ) call outputNamedValue ( 'nsubmembers', nsubmembers )
    ! This number is actually greater than it should be--
    ! another nasty hdf5 gotcha
    do i = 0, nsubmembers-1
      call h5gget_obj_info_idx_f( FileID, '/', i, &
        & objName, objType, h5error )
      if ( DeeBug ) call outputNamedValue ( 'objName', trim(objName) )
      if ( DeeBug ) call outputNamedValue ( 'objType', objType )
      if ( h5error /= 0 ) exit
      if ( objType /= H5G_Group_F ) cycle
      if ( len_trim( GroupNames ) < 1 ) then
        GroupNames = objName
      else
        GroupNames = trim(GroupNames) // ',' // objName
      endif
    enddo
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for all grp names' )
    call trace_end ( cond=.false. )
  end subroutine GetAllHDF5GroupNames

  ! ---------------------------------  MatchHDF5Attributes  -----
  subroutine MatchHDF5Attributes ( MLSFile, attrNames, attrValues, name )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: attrNames ! Names of attrs matched (,-separated)
    character (len=*), intent(in) :: attrValues ! Corresponding values (,-separated)
    character (len=*), intent(out) :: name ! DS name with matching attributes

    ! Local variables
    integer, parameter :: MAXATTRNAMELEN = 128
    logical, parameter :: DEEBUG = .false.
    integer :: attr
    character(len=MAXATTRNAMELEN) :: AttrName
    integer :: ds
    character(len=MAXATTRNAMELEN) :: DSName
    character(len=MAXNAMELEN*MAXNDSNAMES) :: DSNames
    integer :: error
    integer :: numattrs
    character(len=MAXATTRNAMELEN) :: status
    character(len=MAXATTRNAMELEN) :: value

    ! Executable code
    name = ' '
    numAttrs = NumStringElements ( attrNames, countEmpty )
    if ( DEEBUG ) call outputNamedValue( 'numAttrs', numAttrs )
    ! Do we have any FileID fields?
    If ( .not. MLSFile%stillOpen .or. all( &
      & (/ MLSFile%FileId%f_id, MLSFile%FileId%grp_id, MLSFile%FileId%sd_id /) &
      & == 0 ) ) then
      call mls_OpenFile( MLSFile, error )
      if ( error /= 0 ) then
        if ( DEEBUG ) call outputNamedValue( 'error in MatchHDF5Attributes', error )
      return
      endif
    endif
    call GetAllHDF5DSNames( MLSFile, DSNames )
    if ( NumStringElements( DSNames, countEmpty ) < 1 ) return
    allnames: do ds=1, NumStringElements( DSNames, countEmpty )
      DSName = StringElement(DSNames, ds, countEmpty)
      if ( DEEBUG ) call outputNamedValue( 'DS name', trim(DSName) )
      do attr=1, numAttrs
        attrname = StringElement(attrNames, attr, countEmpty)
        if( .not. ReadHDF5Attribute( MLSFile%fileID%f_id, &
          & trim(DSName), trim(attrName), &
          & value, status ) )  then
          if ( DEEBUG ) call outputNamedValue( 'status', status )
          cycle allnames
        endif
        if ( DEEBUG ) then
          call outputNamedValue( 'attr name', StringElement(attrNames, attr, countEmpty) )
          call outputNamedValue( 'its attr value', value )
          call outputNamedValue( 'compared with', StringElement(attrValues, attr, countEmpty) )
        endif
        if ( lowercase(value) /= &
          & lowercase(StringElement(attrValues, attr, countEmpty) ) &
          & ) cycle allnames
      end do
      ! If we survived this long, we must have a match
      name = StringElement(DSNames, ds, countEmpty)
      return
    end do allnames
  end subroutine MatchHDF5Attributes

  !---------------------- ReadHDF5Attr_FID_int ---------------------
  function ReadHDF5Attr_FID_int (fileid, dsName, attrName, value, error) &
    & result(success)
      integer, intent(in) :: fileid
      character(len=*), intent(in) :: dsName
      character(len=*), intent(in) :: attrName
      integer, intent(out) :: value
      character(len=*), optional, intent(out) :: error
      logical :: success

      integer :: setid, attrid, stat, minlen
      character(len=256) :: myerror ! our error message should not exceed 256 characters

      myerror = ' '
      if ( present(error) ) minlen = min(len(error), len(myerror))
      call h5dopen_f (fileid, dsName, setid, stat)
      if (stat /= 0) then
          success = .false.
          if  (present(error)) then
              myerror = 'Cannot open dataset ' // dsName
              error = myerror(1:minlen)
              return
          endif
      endif

      call h5aopen_name_f(setid, attrName, attrid, stat)
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Cannot open attribute ' // attrName // ' of dataset ' // dsName
              error = myerror(1:minlen)
              return
          endif
      endif

      call h5aread_f (attrid, H5T_NATIVE_INTEGER, value, ones, stat)
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Cannot read attribute ' // attrName // ' of dataset ' // dsName
              error = myerror(1:minlen)
              return
          endif
      endif

      call h5aClose_f ( attrID, stat )
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Cannot close attribute ' // attrName // ' of dataset ' // dsName
              error = myerror(1:minlen)
              return
          endif
      endif

      call h5dclose_f (setid, stat)
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Cannot close dataset ' // dsName
              error = myerror(1:minlen)
              return
          endif
      endif
      success = .true.
  end function

  !--------------------- ReadHDF5Attr_FID_string --------------------
  function ReadHDF5Attr_FID_string (fileid, dsName, attrName, value, error) &
    & result(success)
      integer, intent(in) :: fileid
      character(len=*), intent(in) :: dsName
      character(len=*), intent(in) :: attrName
      character(len=*), intent(out) :: value
      character(len=*), optional, intent(out) :: error
      logical :: success

      integer :: setid, attrid, stat, stringtype, minlen
      integer(kind=size_t) :: stringsize
      character(len=256) :: myerror ! error message shouldn't exceed 256 characters

      myerror = ' '
      value = ' '
      if ( present(error) ) minlen = min(len(error), len(myerror))
      call h5dopen_f (fileid, dsName, setid, stat)
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Cannot open dataset ' // dsName
              error = myerror (1:minlen)
              return
          endif
      endif

      call h5aopen_name_f(setid, attrName, attrid, stat)
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Cannot open attribute ' // attrName // ' of dataset ' &
                & // dsName
              error = myerror(1:minlen)
              return
          endif
      endif

      call h5aGet_type_f ( attrID, stringType, stat )
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Unable to find type info for ' // attrName // &
                & ' of dataset ' // dsname
              error = myerror(1:minlen)
              return
          endif
      endif

      call h5tGet_size_f ( stringType, stringSize, stat )
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Unable to find size info for ' // attrName // &
                & ' of dataset ' // dsname
              error = myerror(1:minlen)
              return
          endif
      endif

      if (stringsize > len(value)) then
          success = .false.
          if (present(error)) then
              myerror = 'String is not long enough to store attribute ' // &
                & attrName // ' of dataset ' // dsname
              error = myerror(1:minlen)
              return
          endif
      endif

      call h5aread_f ( attrID, stringType, value(1:stringSize), ones, stat)
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Cannot read attribute ' // attrName // ' of dataset ' &
                & // dsname
              error = myerror(1:minlen)
              return
          endif
      endif

      call h5tClose_f(stringtype, stat)
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Cannot close type for attribute ' // attrName // ' of dataset ' // dsname
              error = myerror(1:minlen)
              return
          endif
      endif

      call h5aClose_f ( attrID, stat )
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Cannot close attribute ' // attrName // ' of dataset ' // dsName
              error = myerror(1:minlen)
              return
          endif
      endif

      call h5dclose_f (setid, stat)
      if (stat /= 0) then
          success = .false.
          if (present(error)) then
              myerror = 'Cannot close dataset ' // dsName
              error = myerror(1:minlen)
              return
          endif
      endif

      success = .true.
  end function

  ! --------------------------------------------  GetHDF5AttrDims  -----
  subroutine GetHDF5AttrDims ( ItemID, name, DIMS, maxDims )
    integer, intent(in) :: ItemID         ! ItemID
    character (len=*), intent(in) :: NAME ! Name of DS
    integer(kind=hSize_t), dimension(:), intent(out) :: DIMS ! Values of dimensions
    integer(kind=hSize_t), dimension(:), optional, intent(out) :: MAXDIMS ! max Values

    ! Local variables
    ! logical, parameter :: DEEBUG = .true.
    integer :: AttrID                   ! ID for Attr
    integer :: classID
    integer :: dspace_id                ! spaceID for Attr
    integer(kind=hSize_t), dimension(:), allocatable :: maxdims_ptr
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: my_rank, rank
    integer :: STATUS                   ! Flag
    integer :: type_id
    integer(kind=Size_t) :: type_size

    ! Executable code
    call trace_begin ( me, 'GetHDF5AttrDims', cond=.false. )
    ! Initializing values returned if there was trouble
    dims = -1
    if ( present(maxDims)) maxDims = -1
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting dims of ' // trim(name) )
    call h5aOpen_name_f ( itemID, trim(name), attrID, status )
    call h5aGet_type_f ( attrID, type_id, status )
    call h5tGet_size_f ( type_id, type_size, status )
    call h5tget_class_f ( type_id, classID, status )
    call h5aget_space_f ( attrID, dspace_id, status )
    if ( status /= 0 ) call outputNamedValue ( 'h5dget_space_f status', status )
    call h5sget_simple_extent_ndims_f ( dspace_id, rank, status )
    if ( status /= 0 ) call outputNamedValue ( 'h5sget_simple_extent_ndims_f status', status )
    my_rank = min(rank, size(dims))
    allocate ( maxdims_ptr(my_rank) )
    call h5sget_simple_extent_dims_f ( dspace_id, dims(1:my_rank), &
         maxdims_ptr, status )
    if ( status /= my_rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dims for attr ' // trim(name) )
    if ( DEEBUG ) then
      call outputNamedValue ( 'name', name )
      call outputNamedValue ( 'attrID', attrID )
      call outputNamedValue ( 'type_id', type_id )
      call outputNamedValue ( 'type_size', int(type_size) )
      call outputNamedValue ( 'classID', classID )
      call outputNamedValue ( 'dspaceID', dspace_ID )
      call outputNamedValue ( 'rank', rank )
      call outputNamedValue ( 'my_rank', my_rank )
      call outputNamedValue ( 'dims', int(dims, kind(status)) )
      call outputNamedValue ( 'maxdims_ptr', int(maxdims_ptr, kind(status)) )
    endif
    call h5aClose_f ( attrID, status )
    if ( present(maxDims) ) maxdims = maxdims_ptr(1:my_rank)
    deallocate ( maxdims_ptr )
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting dims of ' // trim(name) )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5AttrDims

  ! ----------------------------------------------  GetHDF5DSDims  -----
  subroutine GetHDF5DSDims_ID ( FileID, name, DIMS, maxDims )
    integer, intent(in) :: FILEID         ! fileID
    character (len=*), intent(in) :: NAME ! Name of DS
    integer(kind=hSize_t), dimension(:), intent(out) :: DIMS ! Values of dimensions
    integer(kind=hSize_t), dimension(:), optional, intent(out) :: MAXDIMS ! max Values

    ! Local variables
    integer :: dspace_id                ! spaceID for DS
    integer(kind=hSize_t), dimension(:), pointer :: maxdims_ptr
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: my_rank, rank
    integer :: SETID                    ! ID for DS
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'GetHDF5DSDims', cond=.false. )
    ! Initializing values returned if there was trouble
    dims = -1
    if ( present(maxDims)) maxDims = -1
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting dims ' // trim(name) )
    call h5dOpen_f ( FileID, trim(name), setID, status )
    call h5dget_space_f ( setID, dspace_id, status )
    call h5sget_simple_extent_ndims_f ( dspace_id, rank, status )
    my_rank = min(rank, size(dims))
    allocate ( maxdims_ptr(my_rank) )
    call h5sget_simple_extent_dims_f ( dspace_id, dims(1:my_rank), &
         maxdims_ptr, status )
    if ( status /= my_rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dims for datset ' // trim(name) )
    call h5dClose_f ( setID, status )
    if ( present(maxDims) ) maxdims = maxdims_ptr
    deallocate ( maxdims_ptr )
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting dims ' // trim(name) )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5DSDims_ID

  subroutine GetHDF5DSDims_MLSFile ( MLSFile, name, DIMS, maxDims )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of DS
    integer(kind=hSize_t), dimension(:), intent(out) :: DIMS ! Values of dimensions
    integer(kind=hSize_t), dimension(:), optional, intent(out) :: MAXDIMS ! max Values

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'GetHDF5DSDims_MLSFile', cond=.false. )
    ! Do we have any FileID fields?
    If ( .not. MLSFile%stillOpen .or. all( &
      & (/ MLSFile%FileId%f_id, MLSFile%FileId%grp_id, MLSFile%FileId%sd_id /) &
      & == 0 ) ) then
      call h5fopen_f ( trim(MLSFile%name), H5F_ACC_RDONLY_F, MLSFile%fileID%f_id, status )
      call GetHDF5DSDims( MLSFile%fileID%f_id, name, DIMS, maxDims )
      call h5fclose_f ( MLSFile%fileID%f_id, status )
      MLSFile%fileID%f_id = 0
    elseif( MLSFile%FileId%grp_id > 0 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'How can we read Dimensions from a group?', MLSFile=MLSFile )
      ! call GetHDF5DSDims( MLSFile%fileID%grp_id, name, DIMS, maxDims )
    elseif( MLSFile%FileId%f_id > 0 ) then
      call GetHDF5DSDims( MLSFile%fileID%f_id, name, DIMS, maxDims )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'all fields of MLSFile%Fileid are 0', MLSFile=MLSFile )
    endif
    call trace_end ( cond=.false. )
  end subroutine GetHDF5DSDims_MLSFile

  ! ----------------------------------------------  GetHDF5DSRank  -----
  subroutine GetHDF5DSRank_MLSFile ( MLSFile, name, rank )
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name of DS
    integer, intent(out) :: rank        ! How many dimensions

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'GetHDF5DSRank_MLSFile', cond=.false. )
    ! Do we have any FileID fields?
    If ( .not. MLSFile%stillOpen .or. all( &
      & (/ MLSFile%FileId%f_id, MLSFile%FileId%grp_id, MLSFile%FileId%sd_id /) &
      & == 0 ) ) then
      call h5fopen_f ( trim(MLSFile%name), H5F_ACC_RDONLY_F, MLSFile%fileID%f_id, status )
      call GetHDF5DSRank_ID( MLSFile%fileID%f_id, name, rank )
      call h5fclose_f ( MLSFile%fileID%f_id, status )
      MLSFile%fileID%f_id = 0
    elseif( MLSFile%FileId%grp_id > 0 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'How can we read rank from a group?', MLSFile=MLSFile )
      ! call GetHDF5DSRank_ID( MLSFile%fileID%grp_id, name, rank )
    elseif( MLSFile%FileId%f_id > 0 ) then
      call GetHDF5DSRank_ID( MLSFile%fileID%f_id, name, rank )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'all fields of MLSFile%Fileid are 0', MLSFile=MLSFile )
    endif
    call trace_end ( cond=.false. )
  end subroutine GetHDF5DSRank_MLSFile

  subroutine GetHDF5DSRank_ID ( FileID, name, rank )
    integer, intent(in) :: FILEID       ! fileID
    character (len=*), intent(in) :: NAME ! Name of DS
    integer, intent(out) :: rank        ! How many dimensions

    ! Local variables
    integer :: dspace_id                ! spaceID for DS
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID for DS
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'GetHDF5DSRank', cond=.false. )
    rank = -1                           ! means trouble
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting rank ' // trim(name) )
    call h5dOpen_f ( FileID, trim(name), setID, status )
    call h5dget_space_f(setID,dspace_id,status)
    call h5sget_simple_extent_ndims_f(dspace_id,rank,status)
    call h5dClose_f ( setID, status )
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting rank ' // trim(name) )
    call trace_end ( cond=.false. )
  end subroutine GetHDF5DSRank_ID

  ! ---------------------------------------------  GetHDF5DSQType  -----
  subroutine GetHDF5DSQType ( FileID, name, Qtype )
    integer, intent(in) :: FILEID       ! fileID
    character (len=*), intent(in) :: NAME ! Name of DS
    character (len=*), intent(out) :: Qtype    ! 'real' or 'integer' or ..

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID for DS
    integer :: STATUS                   ! Flag
    integer :: type_id                  ! typeID for DS

    ! Executable code
    call trace_begin ( me, 'GetHDF5DSQType', cond=.false. )
    Qtype = 'unknown'                   ! means trouble
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting datatype of ' // trim(name) )
    call h5dOpen_f ( FileID, trim(name), setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open before getting datatype of ' // trim(name) )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get datatype ' // trim(name) )
    call h5dget_type_f( setID, type_id, status )
    Qtype = WhatTypeAmI ( type_id )
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close after datatype ' // trim(name) )
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting rank ' // trim(name))
    call trace_end ( cond=.false. )
  end subroutine GetHDF5DSQType

  ! --------------------------------------  IsHDF5AttributeInFile_DS  -----
  function IsHDF5AttributeInFile_DS ( filename, DSname, name ) &
    & result(SooDesu)
    ! This routine returns true if the given HDF5 DS is present
    ! Possible refiements:
    ! Add interface for MLSFile type
    character (len=*), intent(in) :: FILENAME ! Where to look
    character (len=*), intent(in) :: DSNAME ! Name for the dataset
    character (len=*), intent(in) :: NAME ! Name for the attribute
    logical                       :: SOODESU

    ! Local variables
    integer :: ATTRID                   ! ID for attribute if present
    integer :: fileID                   ! Where to look
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID for DS if present
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'IsHDF5AttributeInFile_DS', cond=.false. )
    SooDesu = .false.
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for attribute ' // trim(name) )
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, fileID, status)
    if ( status == 0 ) then
      call h5dOpen_f ( fileid, trim(DSname), setID, status )
      if ( status == 0 ) then
        call h5aopen_name_f(setID, trim(name), attrid, status)
        if ( status == 0 ) then
          SooDesu = .true.
          call h5aclose_f(attrid, status)
        end if
        call h5dClose_f ( setID, status )
      end if
      call h5fclose_f(fileID, status)
    end if
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute ' // trim(name) )
    call trace_end ( cond=.false. )
  end function IsHDF5AttributeInFile_DS

  ! --------------------------------------  IsHDF5AttributeInFile_Grp  -----
  function IsHDF5AttributeInFile_Grp ( filename, name ) &
    & result(SooDesu)
    ! This routine returns true if the given HDF5 DS is present in group '/'
    ! Possible refiements:
    ! Check for attributes to other groups beside root
    ! Add interface for MLSFile type
    character (len=*), intent(in) :: FILENAME ! Where to look
    character (len=*), intent(in) :: NAME ! Name for the attribute
    logical                       :: SOODESU
    ! Local variables
    integer :: ATTRID                   ! ID for attribute if present
    integer :: fileID                   ! Where to look
    integer :: GRPID                    ! ID for group if present
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'IsHDF5AttributeInFile_Grp', cond=.false. )
    SooDesu = .false.
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for attribute ' // trim(name) )
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, fileID, status)
    if ( status == 0 ) then
      call h5gOpen_f ( fileid, '/', grpID, status )
      if ( status == 0 ) then
        call h5aopen_name_f(grpID, trim(name), attrid, status)
        if ( status == 0 ) then
          SooDesu = .true.
          call h5aclose_f(attrid, status)
        end if
        call h5gClose_f ( grpID, status )
      end if
      call h5fclose_f(fileID, status)
    end if
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute ' // trim(name) )
    call trace_end ( cond=.false. )
  end function IsHDF5AttributeInFile_Grp

  ! ---------------------------------------------  IsHDF5DSInFile  -----
  logical function IsHDF5DSInFile ( filename, name )
    ! This routine returns true if the given HDF5 DS is present
    character (len=*), intent(in) :: FILENAME ! Where to look
    character (len=*), intent(in) :: NAME ! Name for the dataset

    ! Local variables
    integer :: FILEID                   ! Where to look
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID for DS if present
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'IsHDF5DSInFile', cond=.false. )
    IsHDF5DSInFile = .false.
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for DS ' // trim(name) )
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, fileID, status)
    if ( status == 0 ) then
      call h5dOpen_f ( fileid, trim(name), setID, status )
      if ( status == 0 ) then
        IsHDF5DSInFile = .true.
        call h5dClose_f ( setID, status )
      end if
      call h5fclose_f(fileID, status)
    end if
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for DS ' // trim(name) )
    call trace_end ( cond=.false. )
  end function IsHDF5DSInFile

  ! -----------------------------  IsHDF5AttributePresent_in_DSID  -----
  logical function IsHDF5AttributePresent_in_DSID ( SETID, name )
    ! This routine returns true if the given HDF5 attribute is present
    integer, intent(in) :: SETID        ! Dataset ID--Where to look
    character (len=*), intent(in) :: NAME ! Name for the attribute

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: ATTRID                   ! ID for attribute if present
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'IsHDF5AttributePresent_in_DSID', cond=.false. )
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for attribute ' // trim(name) )
    call h5aOpen_name_f ( SETID, name, attrID, status )
    if ( status /= 0 ) then
      IsHDF5AttributePresent_in_DSID = .false.
    else
      IsHDF5AttributePresent_in_DSID = .true.
      call h5aClose_f ( attrID, status )
    end if
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute ' // trim(name) )
    call trace_end ( cond=.false. )
  end function IsHDF5AttributePresent_in_DSID

  ! ------------------------------  IsHDF5AttributePresent_in_fID  -----
  logical function IsHDF5AttributePresent_in_fID ( fileID, DSname, name, &
    & is_grpattr )
    ! This routine returns true if the given HDF5 attribute is present
    ! If present it would be as an attribute either of the dataset DSname
    ! or, if is_grpattr is TRUE, of the group DSname
    integer, intent(in) :: fileID        ! file ID--Where to look
    character (len=*), intent(in) :: DSNAME ! Name for the dataset
    character (len=*), intent(in) :: NAME ! Name for the attribute
    logical, optional, intent(in) :: is_grpattr ! DSNAME is a group name

    ! Local variables
    integer :: ATTRID                   ! ID for attribute if present
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_grpattr
    integer :: SETID                    ! ID for dataset if present
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'IsHDF5AttributePresent_in_fID', cond=.false. )
    my_grpattr = .false.
    if ( present(is_grpattr) ) my_grpattr = is_grpattr
    if ( my_grpattr ) then
      IsHDF5AttributePresent_in_fID = IsHDF5AttributePresent_in_grp ( &
        & DSname, fileID, name)
      go to 9
    end if
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for attribute ' // trim(name) )
    call h5dOpen_f ( fileid, trim(DSname), setID, status )
    if ( status /= 0 ) then
      IsHDF5AttributePresent_in_fID = .false.
    else
      call h5aOpen_name_f ( SETID, name, attrID, status )
      if ( status /= 0 ) then
        IsHDF5AttributePresent_in_fID = .false.
      else
        IsHDF5AttributePresent_in_fID = .true.
        call h5aClose_f ( attrID, status )
      end if
      call h5dClose_f ( setID, status )
    end if
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute ' // trim(name) )
  9 call trace_end ( cond=.false. )
  end function IsHDF5AttributePresent_in_fID

  ! ------------------------------  IsHDF5AttributePresent_in_grp  -----
  logical function IsHDF5AttributePresent_in_grp ( grpName, fileID, name )
    ! This routine returns true if the given HDF5 attribute is present
    ! (Note the unusual ordering of args--to make unambiguous)
    integer, intent(in) :: fileID        ! file ID--Where to look
    character (len=*), intent(in) :: grpName ! Name for the group
    character (len=*), intent(in) :: NAME ! Name for the attribute

    ! Local variables
    integer :: ATTRID                   ! ID for attribute if present
    integer :: GRPID                    ! ID for group if present
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'IsHDF5AttributePresent_in_grp', cond=.false. )
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for attribute ' // trim(name) )
    call h5gopen_f(fileID, trim(grpName), grpid, status)
    if ( status /= 0 ) then
      IsHDF5AttributePresent_in_grp = .false.
    else
      call h5aOpen_name_f ( GRPID, name, attrID, status )
      if ( status /= 0 ) then
        IsHDF5AttributePresent_in_grp = .false.
      else
        IsHDF5AttributePresent_in_grp = .true.
        call h5aClose_f ( attrID, status )
      end if
      call h5gclose_f(grpid, status)
    end if
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute ' // trim(name) )
    call trace_end ( cond=.false. )
  end function IsHDF5AttributePresent_in_grp

  ! ------------------------------  IsHDF5AttributePresent_in_MLSFile  -----
  logical function IsHDF5AttributePresent_in_MLSFile ( MLSFile, name )
    ! This routine returns true if the given HDF5 attribute is present
    ! somewhere under the MLSFile tree
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for the attribute

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'IsHDF5AttributePresent_in_MLSFile', cond=.false. )
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for attribute ' // trim(name) )
    if ( MLSFile%fileID%sd_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%sd_ID, name, attrID, status )
    elseif ( MLSFile%fileID%grp_ID > 0 ) then
      call h5aOpen_name_f ( MLSFile%fileID%grp_ID, name, attrID, status )
    else
      status = -1
    endif
    IsHDF5AttributePresent_in_MLSFile = ( status == 0 )
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute ' // trim(name) )
    call trace_end ( cond=.false. )
  end function IsHDF5AttributePresent_in_MLSFile

  ! --------------------------------------------  IsHDF5DSPresent_in_MLSFile  -----
  logical function IsHDF5DSPresent_in_MLSFile ( MLSFile, name, options )
    ! This routine returns true if the given HDF5 DS is present
    ! options such as "-cw" are explained in notes above
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for the dataset
    character (len=*), optional, intent(in) :: options ! E.g., -c

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: STATUS                   ! Flag
    ! logical, parameter :: DEEBUG = .true.

    ! Executable code
    IsHDF5DSPresent_in_MLSFile = .false.
    call trace_begin ( me, 'IsHDF5DSPresent_in_MLSFile', cond=.false. )
    ! if ( .not. MLSFile%stillOpen .or. all( &
    !   & (/ MLSFile%FileId%f_id, MLSFile%FileId%grp_id, MLSFile%FileId%sd_id /) &
    !   & == 0 ) ) then
    if ( .not. MLSFile%stillOpen ) then
      call MLS_Openfile ( MLSFile, status )
      if ( status /= 0 ) then
        if ( DEEBUG ) &
          & call outputNamedValue( 'error in IsHDF5DSPresent_in_MLSFile', status )
      return
      endif
    endif
    if ( MLSFile%fileID%grp_ID > 0 ) then
      IsHDF5DSPresent_in_MLSFile = &
        & IsHDF5DSPresent_in_fID ( MLSFile%fileID%grp_ID, name, options )
    elseif ( MLSFile%FileID%f_id > 0 ) then
      IsHDF5DSPresent_in_MLSFile = &
        & IsHDF5DSPresent_in_fID ( MLSFile%FileID%f_id, name, options )
    else
      status = -1
    endif
    if ( DEEBUG .and.  .not. IsHDF5DSPresent_in_MLSFile ) &
      & call Dump( MLSFile, details=2 )
    call trace_end ( cond=.false. )
  end function IsHDF5DSPresent_in_MLSFile

  logical function IsHDF5DSPresent_in_fID ( locID, name, options )
    ! This routine returns true if the given HDF5 DS is present
    ! options such as "-cw" are explained in notes above
    integer, intent(in) :: LOCID        ! Where to look
    character (len=*), intent(in) :: NAME ! Name for the dataset
    character (len=*), optional, intent(in) :: options ! E.g., -c
    ! Local variables
    character (len=MAXNDSNAMES*MAXNAMELEN) :: DSNames ! Names of DS in file (,-separated)
    character(len=8) :: myOptions

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID for DS if present
    integer :: STATUS                   ! Flag
    logical, parameter :: DEEBUG = .false.

    ! Executable code
    IsHDF5DSPresent_in_fID = .false.
    call trace_begin ( me, 'IsHDF5DSPresent_in_fID', cond=.false. )
    myOptions = ' ' ! By default, match name exactly
    if (present(options)) myOptions = options
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for DS ' // trim(name) )
    if ( myOptions == ' ' ) then
      call h5dOpen_f ( locID, name, setID, status )
      IsHDF5DSPresent_in_fID = ( status == 0 )
      if ( IsHDF5DSPresent_in_fID ) then
        call h5dClose_f ( setID, status )
      elseif ( DEEBUG ) then
        call MLSMessage ( MLSMSG_Info, ModuleName, &
          & 'Unable to open DS ' // trim(name) )
      endif
    elseif ( index(myOptions, 'a') < 1 ) then
      call GetAllHDF5DSNames_fileID ( locID, '/', DSNames )
      IsHDF5DSPresent_in_fID = IsInList( DSNames, name, options )
    else
      call GetAllHDF5AttrNames( locID, DSNames )
      IsHDF5DSPresent_in_fID = IsInList( DSNames, name, options )
    endif
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for DS ' // trim(name) )
    call trace_end ( cond=.false. )
  end function IsHDF5DSPresent_in_fID

  ! -----------------------------------------  IsHDF5GroupPresent  -----
  logical function IsHDF5GroupPresent ( locID, name )
    ! This routine returns true if the given HDF5 Group is present
    integer, intent(in) :: LOCID        ! Where to look
    character (len=*), intent(in) :: NAME ! Name for the dataset

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID for DS if present
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'IsHDF5GroupPresent', cond=.false. )
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for Group ' // trim(name) )
    call h5GOpen_f ( locID, name, setID, status )
    IsHDF5groupPresent = status == 0
    if ( IsHDF5groupPresent ) call h5GClose_f ( setID, status )
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for Group ' // trim(name) )
    call trace_end ( cond=.false. )
  end function IsHDF5GroupPresent

  ! --------------------------------------------  IsHDF5ItemPresent  -----
  logical function IsHDF5ItemPresent ( locID, name, options )
    ! This routine returns true if the given HDF5 itemis present
    ! options such as "-cw" are explained in notes above
    ! with additional options specifying the item type as follows
    ! d => DS, g => group, a => attribute
    integer, intent(in) :: LOCID        ! Where to look
    character (len=*), intent(in) :: NAME ! Name for the dataset
    character (len=*), optional, intent(in) :: options ! E.g., -c

    ! Local variables
    character (len=MAXNDSNAMES*MAXNAMELEN) :: DSNames ! Names of DS in file (,-separated)
    integer :: Me = -1                  ! String index for trace cacheing
    character(len=8) :: myOptions
    integer :: SETID                    ! ID for DS if present
    integer :: STATUS                   ! Flag

    ! Executable code
    call trace_begin ( me, 'IsHDF5ItemPresent', cond=.false. )
    myOptions = ' ' ! By default, match name exactly, search for sd names
    if (present(options)) myOptions = options
    call mls_eSet_auto ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for item ' // trim(name) )
    if ( myOptions == ' ' .or. myOptions == '-d' ) then
      call h5dOpen_f ( locID, name, setID, status )
      IsHDF5ItemPresent = ( status == 0 )
      if ( IsHDF5ItemPresent ) call h5dClose_f ( setID, status )
    elseif ( index(myOptions, 'g') > 0 ) then
      IsHDF5ItemPresent = IsHDF5GroupPresent (locID, name)
    elseif ( index(myOptions, 'a') < 1 ) then
      call GetAllHDF5DSNames_fileID ( locID, '/', DSNames )
      IsHDF5ItemPresent = IsInList( DSNames, name, options )
      ! print *, 'DSNames ', trim(DSNames)
    else
      call GetAllHDF5AttrNames( locID, DSNames )
      IsHDF5ItemPresent = IsInList( DSNames, name, options )
      ! print *, 'attr Names ', trim(DSNames)
    endif
    call mls_eSet_auto ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for item ' // trim(name) )
    call trace_end ( cond=.false. )
  end function IsHDF5ItemPresent

  ! --------------------------------------  SaveAsHDF5DS_charsclr  -----
  subroutine SaveAsHDF5DS_charsclr ( locID, name, value, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID           ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME  ! Name for this dataset
    character (len=*), intent(in) :: VALUE ! The scalar char string
    logical, optional, intent(in)     :: adding_to

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myAddingTo
    integer :: nsubmembers
    integer (HID_T) :: setID            ! ID for dataset
    integer(kind=hsize_t), dimension(1) :: SHP        ! Shape
    integer :: spaceID                  ! ID for dataspace
    integer :: status                   ! Flag from HDF5
    integer(hid_t) :: s_type_id
    integer(hid_t) :: type_id

    ! Executable code
    call trace_begin ( me, 'SaveAsHDF5DS_charsclr', cond=.false. )
    ! Create the dataspace
    myAddingTo = .false.
    if ( present(adding_to) ) myAddingTo = adding_to
    shp = 1
    call h5sCreate_simple_f ( 1, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for scalar character ' // trim(name) )
    type_id = H5T_NATIVE_CHARACTER
    call h5tcopy_f ( type_id, s_type_id, status )
    call h5tset_size_f ( s_type_id, int(len(value), size_t), status )
    ! Create the dataset
    if ( myAddingTo ) then
      call h5dOpen_f ( locID, trim(name), setID, status )
    else
      call h5dCreate_f ( locID, trim(name), s_type_id, spaceID, setID, &
        & status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for scalar character ' // trim(name) )
    endif
    ! Write the data
    call h5dWrite_f ( setID, s_type_id, value, &
      & int ( (/ shp, ones(1:6) /), hsize_t ), status )
    call finishSaveDS ( Name, status, setID, spaceID )
    call trace_end ( cond=.false. )
  end subroutine SaveAsHDF5DS_charsclr

  ! --------------------------------------  SaveAsHDF5DS_chararr1  -----
  subroutine SaveAsHDF5DS_chararr1 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(in) :: VALUE(:)     ! The array itself
    character (len=*), optional, intent(in) :: FILLVALUE
    integer, dimension(1), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    include 'SaveAsHDF5DS_chararr.f9h'

  end subroutine SaveAsHDF5DS_chararr1

  ! --------------------------------------  SaveAsHDF5DS_chararr2  -----
  subroutine SaveAsHDF5DS_chararr2 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(in) :: VALUE(:,:)     ! The array itself
    character (len=*), optional, intent(in) :: FILLVALUE
    integer, dimension(2), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    include 'SaveAsHDF5DS_chararr.f9h'

  end subroutine SaveAsHDF5DS_chararr2

  ! --------------------------------------  SaveAsHDF5DS_textFile  -----
  subroutine SaveAsHDF5DS_textFile ( textFile, locID, name, &
    & maxLineLen, fromNull, adding_to, cr2lf )
    use IO_Stuff, only: Read_TextFile
    ! This routine writes the contents of a textfile as a char-valued dataset
    integer, intent(in) :: LOCID           ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME  ! Name for this dataset
    character (len=*), intent(in) :: textFile ! Name of the textfile
    integer, optional, intent(in) :: maxLineLen
    character(len=1), optional, intent(in) :: fromNull
    logical, optional, intent(in)     :: adding_to ! name may already exits
    logical, optional, intent(in)     :: cr2lf ! convert <cr> to <lf>

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myAddingTo
    logical :: myConvertCR
    integer :: myMaxLineLen
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer(kind=hsize_t), dimension(1) :: SHP        ! Shape
    integer(hid_t) :: s_type_id
    integer(hid_t) :: type_id
    character(LEN=MAXTEXTSIZE)              :: value
    character(len=1) :: wasNull

    ! Executable code
    call trace_begin ( me, 'SaveAsHDF5DS_textFile', cond=.false. )
    myMaxLineLen = DFLTMAXLINELENGTH
    if ( present(maxLineLen) ) myMaxLineLen = maxLineLen
    wasNull = char(32) ! blank space
    if ( present(fromNull) ) wasNull = fromNull
    myAddingTo = .false.
    if ( present(adding_to) ) myAddingTo = adding_to
    myConvertCR = .false.
    if ( present(cr2lf) ) myConvertCR = cr2lf
    ! Try to read the textfile
    value = ' '
    call read_textFile( trim(textFile), value, myMaxLineLen )
    ! Unfortunately, a lot of null characters sneak into this
    value = Replace( value, char(0), wasNull ) ! Replace null with wasNull
    if ( myConvertCR ) &
      & value = Replace( value, char(13), char(10) ) ! Replace null with wasNull
    
    ! Create the dataspace
    shp = 1
    call h5sCreate_simple_f ( 1, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for scalar character ' // trim(name) )
    type_id = H5T_NATIVE_CHARACTER
    call h5tcopy_f ( type_id, s_type_id, status )
    call h5tset_size_f ( s_type_id, int(len_trim(value), size_t), status )

    if ( myAddingTo ) then
      call h5dOpen_f ( locID, trim(name), setID, status )
    else
      ! Create the dataset
      call h5dCreate_f ( locID, trim(name), s_type_id, spaceID, setID, &
        & status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for scalar character ' // trim(name) )
    endif
    ! Write the data
    call h5dWrite_f ( setID, s_type_id, trim(value), &
      & int ( (/ shp, ones(1:6) /), hsize_t ), status )
    call finishSaveDS ( name, status, setID, spaceID )
    call trace_end ( cond=.false. )
  end subroutine SaveAsHDF5DS_textFile

  ! ---------------------------------------  SaveAsHDF5DS_intarr1  -----
  subroutine SaveAsHDF5DS_intarr1 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(in) :: VALUE(:)     ! The array itself
    integer, optional, intent(in) :: FILLVALUE
    integer, dimension(1), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    include 'SaveAsHDF5DS_intarr.f9h'

  end subroutine SaveAsHDF5DS_intarr1

  ! ---------------------------------------  SaveAsHDF5DS_intarr2  -----
  subroutine SaveAsHDF5DS_intarr2 ( locID, name, value, &
    & fillValue, finalShape, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(in) :: VALUE(:,:)     ! The array itself
    integer, optional, intent(in) :: FILLVALUE
    integer, dimension(2), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    include 'SaveAsHDF5DS_intarr.f9h'

  end subroutine SaveAsHDF5DS_intarr2

  ! ---------------------------------------  SaveAsHDF5DS_intarr3  -----
  subroutine SaveAsHDF5DS_intarr3 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(in) :: VALUE(:,:,:)     ! The array itself
    integer, optional, intent(in) :: FILLVALUE
    integer, dimension(3), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    include 'SaveAsHDF5DS_intarr.f9h'

  end subroutine SaveAsHDF5DS_intarr3

  ! ---------------------------------------  SaveAsHDF5DS_intarr4  -----
  subroutine SaveAsHDF5DS_intarr4 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(in) :: VALUE(:,:,:,:) ! The array itself
    integer, optional, intent(in) :: FILLVALUE
    integer, dimension(4), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    include 'SaveAsHDF5DS_intarr.f9h'

  end subroutine SaveAsHDF5DS_intarr4

  ! ---------------------------------------  SaveAsHDF5DS_logarr1  -----
  subroutine SaveAsHDF5DS_logarr1 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    logical, intent(in) :: VALUE(:)     ! The array itself
    logical, optional, intent(in) :: FILLVALUE
    integer, dimension(1), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    character :: MyFillValue, MyValue(size(value))   ! T = true, F = false

    ! Executable code
    call trace_begin ( me, 'SaveAsHDF5DS_logarr1', cond=.false. )
    ! Turn fillValue into a character
    myFillValue = 'F'
    if ( present(fillValue) ) then
      if ( fillValue ) myFillValue = 'T'
    end if
    ! Turn value into a character
    where ( value )
      myValue = 'T'
    elsewhere
      myValue = 'F'
    endwhere
    ! write ( myValue, '(L1)' ) value
    call saveAsHDF5DS ( locID, name, myValue, finalShape, myFillValue, adding_to )
    call trace_end ( cond=.false. )
  end subroutine SaveAsHDF5DS_logarr1

  ! ---------------------------------------  SaveAsHDF5DS_dblarr1  -----
  subroutine SaveAsHDF5DS_dblarr1 ( locID, name, value, &
    & start, count, stride, block, may_add_to, adding_to, fillValue )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(in) :: VALUE(:)     ! The array itself
    double precision, optional, intent(in) :: FILLVALUE

    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, optional, intent(in)               :: may_add_to
                                 ! May call again to add to same dataset
    logical, optional, intent(in)               :: adding_to
                                 ! Calling again to add to same dataset
    include 'SaveAsHDF5DS_dblarr.f9h'

  end subroutine SaveAsHDF5DS_dblarr1

  ! ---------------------------------------  SaveAsHDF5DS_dblarr2  -----
  subroutine SaveAsHDF5DS_dblarr2 ( locID, name, value, &
    & start, count, stride, block, may_add_to, adding_to, fillValue )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(in) :: VALUE(:,:)  ! The array itself
    double precision, optional, intent(in) :: FILLVALUE

    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, optional, intent(in)               :: may_add_to
                                 ! May call again to add to same dataset
    logical, optional, intent(in)               :: adding_to
                                 ! Calling again to add to same dataset
    include 'SaveAsHDF5DS_dblarr.f9h'

  end subroutine SaveAsHDF5DS_dblarr2

  ! ---------------------------------------  SaveAsHDF5DS_dblarr3  -----
  subroutine SaveAsHDF5DS_dblarr3 ( locID, name, value, &
    & start, count, stride, block, may_add_to, adding_to, fillValue )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(in) :: VALUE(:,:,:)  ! The array itself
    double precision, optional, intent(in) :: FILLVALUE

    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, optional, intent(in)               :: may_add_to
                                 ! May call again to add to same dataset
    logical, optional, intent(in)               :: adding_to
                                 ! Calling again to add to same dataset
    include 'SaveAsHDF5DS_dblarr.f9h'

  end subroutine SaveAsHDF5DS_dblarr3

  ! ---------------------------------------  SaveAsHDF5DS_dblarr4  -----
  subroutine SaveAsHDF5DS_dblarr4 ( locID, name, value, &
    & start, count, stride, block, may_add_to, adding_to, fillValue )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(in) :: VALUE(:,:,:,:)  ! The array itself
    double precision, optional, intent(in) :: FILLVALUE

    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, optional, intent(in)               :: may_add_to
                                 ! May call again to add to same dataset
    logical, optional, intent(in)               :: adding_to
                                 ! Calling again to add to same dataset
    include 'SaveAsHDF5DS_dblarr.f9h'

  end subroutine SaveAsHDF5DS_dblarr4

  ! --------------------------------------  SaveAsHDF5DS_snglarr1  -----
  subroutine SaveAsHDF5DS_snglarr1 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(in) :: VALUE(:)        ! The array itself
    real, optional, intent(in) :: FILLVALUE
    integer, dimension(1), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer (HID_T) :: setID            ! ID for dataset
    integer(kind=hsize_t), dimension(1) :: SHP, MAXDIMS        ! Shape
    integer :: spaceID                  ! ID for dataspace
    integer :: status                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'SaveAsHDF5DS_snglarr1', cond=.false. )
    shp = shape(value)
    maxdims=shp
    if ( present(finalShape) ) maxdims = finalShape
    ! Create the dataspace
    call createSpaceSet ( locID, name, maxdims, H5T_NATIVE_REAL, &
      & spaceID, setID, status, adding_to, rFill=fillValue )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 1D real array ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
      & int ( (/ shp, ones(1:6) /), hsize_t ), status )
    call finishSaveDS ( name, status, setID, spaceID )
    call trace_end ( cond=.false. )
  end subroutine SaveAsHDF5DS_snglarr1

  ! --------------------------------------  SaveAsHDF5DS_snglarr2  -----
  subroutine SaveAsHDF5DS_snglarr2 ( locID, name, value, &
    & start, count, stride, block, may_add_to, adding_to, fillValue )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(in) :: VALUE(:,:)      ! The array itself
    real, optional, intent(in) :: FILLVALUE

    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, optional, intent(in)               :: may_add_to
                                 ! May call again to add to same dataset
    logical, optional, intent(in)               :: adding_to
                                 ! Calling again to add to same dataset

    include 'SaveAsHDF5DS_snglarr.f9h'

  end subroutine SaveAsHDF5DS_snglarr2

  ! --------------------------------------  SaveAsHDF5DS_snglarr3  -----
  subroutine SaveAsHDF5DS_snglarr3 ( locID, name, value, &
    & start, count, stride, block, may_add_to, adding_to, fillValue )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(in) :: VALUE(:,:,:)    ! The array itself
    real, optional, intent(in) :: FILLVALUE
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, optional, intent(in)               :: may_add_to
                                 ! May call again to add to same dataset
    logical, optional, intent(in)               :: adding_to
                                 ! Calling again to add to same dataset

    include 'SaveAsHDF5DS_snglarr.f9h'

  end subroutine SaveAsHDF5DS_snglarr3

  ! --------------------------------------  SaveAsHDF5DS_snglarr4  -----
  subroutine SaveAsHDF5DS_snglarr4 ( locID, name, value, &
    & start, count, stride, block, may_add_to, adding_to, fillValue )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(in) :: VALUE(:,:,:,:)  ! The array itself
    real, optional, intent(in) :: FILLVALUE
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, optional, intent(in)               :: may_add_to
                                 ! May call again to add to same dataset
    logical, optional, intent(in)               :: adding_to
                                 ! Calling again to add to same dataset

    include 'SaveAsHDF5DS_snglarr.f9h'

  end subroutine SaveAsHDF5DS_snglarr4

  ! ------------------------------------  LdFrmHDF5DS_ID_charscalar  -----
  subroutine LdFrmHDF5DS_ID_charscalar ( locID, name, value )
    ! This routine loads a scalar with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(out) :: VALUE    ! The scalar itself

    character(len=*), parameter :: Sfx = 'charscalar'

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    type (MLSFile_T)   :: MLSFile
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'LdFromHDF5DS_ID_'//sfx , cond=.false. )
    status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_RDONLY, &
      & name='unknown', shortName='unknown', HDFVersion=HDFVERSION_5 )
    MLSFile%fileID%sd_id = locID
    MLSFile%stillOpen = .true.
    call LoadFromHDF5DS ( MLSFile, name, value )
    call trace_end ( cond=.false. )

  end subroutine LdFrmHDF5DS_ID_charscalar

  ! ------------------------------------  LdFrmHDF5DS_ID_chararr1  -----
  subroutine LdFrmHDF5DS_ID_chararr1 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(out) :: VALUE(:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinates of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'chararr1'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_chararr1

  ! ------------------------------------  LdFrmHDF5DS_ID_chararr2  -----
  subroutine LdFrmHDF5DS_ID_chararr2 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(out) :: VALUE(:,:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'chararr2'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_chararr2

  ! ------------------------------------  LdFrmHDF5DS_ID_chararr3  -----
  subroutine LdFrmHDF5DS_ID_chararr3 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(out) :: VALUE(:,:,:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'chararr3'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_chararr3

  ! -------------------------------------  LdFrmHDF5DS_ID_intarr1  -----
  subroutine LdFrmHDF5DS_ID_intarr1 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(out) :: VALUE(:)      ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'intarr1'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_intarr1

  ! -------------------------------------  LdFrmHDF5DS_ID_intarr2  -----
  subroutine LdFrmHDF5DS_ID_intarr2 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(out) :: VALUE(:,:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'intarr2'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_intarr2

  ! -------------------------------------  LdFrmHDF5DS_ID_intarr3  -----
  subroutine LdFrmHDF5DS_ID_intarr3 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(out) :: VALUE(:,:,:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'intarr3'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_intarr3

  ! -------------------------------------  LdFrmHDF5DS_ID_intarr4  -----
  subroutine LdFrmHDF5DS_ID_intarr4 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID           ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME  ! Name for this dataset
    integer, intent(out) :: VALUE(:,:,:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'intarr4'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_intarr4

  ! -------------------------------------  LdFrmHDF5DS_ID_logarr1  -----
  subroutine LdFrmHDF5DS_ID_logarr1 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    logical, intent(out) :: VALUE(:)      ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    character :: MyValue(size(value))   ! 'F' = false, 'T' = true

    ! Executable code
    call trace_begin ( me, 'LdFromHDF5DS_ID_logarr1', cond=.false. )
    call LoadFromHDF5DS ( locID, name, myValue, start, count, stride, block )
    value = myValue == 'T'
    call trace_end ( cond=.false. )
  end subroutine LdFrmHDF5DS_ID_logarr1

  ! -------------------------------------  LdFrmHDF5DS_ID_dblarr1  -----
  subroutine LdFrmHDF5DS_ID_dblarr1 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME     ! Name for this dataset
    double precision, intent(out) :: VALUE(:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'dblarr1'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_dblarr1

  ! -------------------------------------  LdFrmHDF5DS_ID_dblarr2  -----
  subroutine LdFrmHDF5DS_ID_dblarr2 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME       ! Name for this dataset
    double precision, intent(out) :: VALUE(:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'dblarr2'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_dblarr2

  ! -------------------------------------  LdFrmHDF5DS_ID_dblarr3  -----
  subroutine LdFrmHDF5DS_ID_dblarr3 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(out) :: VALUE(:,:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'dblarr3'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_dblarr3

  ! -------------------------------------  LdFrmHDF5DS_ID_dblarr4  -----
  subroutine LdFrmHDF5DS_ID_dblarr4 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(out) :: VALUE(:,:,:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'dblarr4'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_dblarr4

  ! ------------------------------------  LdFrmHDF5DS_ID_snglarr1  -----
  subroutine LdFrmHDF5DS_ID_snglarr1 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(out) :: VALUE(:)         ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'snglarr1'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_snglarr1

  ! ------------------------------------  LdFrmHDF5DS_ID_snglarr2  -----
  subroutine LdFrmHDF5DS_ID_snglarr2 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(out) :: VALUE(:,:)       ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'snglarr2'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_snglarr2

  ! ------------------------------------  LdFrmHDF5DS_ID_snglarr3  -----
  subroutine LdFrmHDF5DS_ID_snglarr3 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(out) :: VALUE(:,:,:)     ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'snglarr3'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_snglarr3

  ! ------------------------------------  LdFrmHDF5DS_ID_snglarr4  -----
  subroutine LdFrmHDF5DS_ID_snglarr4 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(out) :: VALUE(:,:,:,:)   ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    character(len=*), parameter :: Sfx = 'snglarr4'

    include 'LdFrmHDF5DS_ID.f9h'

  end subroutine LdFrmHDF5DS_ID_snglarr4

  ! ------------------------------------  LoadFromHDF5DS_charscalar  -----
  subroutine LoadFromHDF5DS_charscalar ( MLSFile, name, value )
    ! This routine loads a scalar with values from a DS
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(out) :: VALUE    ! The scalar itself

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer :: SPACEID                  ! ID of dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=Size_t) :: STRINGSIZE               ! String size
    integer :: STRINGTYPE               ! String type

    value = ' '
    call trace_begin ( me, 'LoadFromHDF5DS_charscalar', cond=.false. )
    call h5dOpen_f ( MLSFile%fileID%sd_id, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name), &
      & MLSFile=MLSFile )
    call h5dGet_type_f ( setID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for scalar ' // trim(name), &
      & MLSFile=MLSFile )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for scalar ' // trim(name), &
      & MLSFile=MLSFile )
    if ( stringSize > len(value) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for scalar ' // trim(name), &
      & MLSFile=MLSFile )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name), &
      & MLSFile=MLSFile )
    call h5dread_f ( setID, stringtype, value, ones(1:7), status )
    ! In case there are any nulls, we will replace them with blanks
    status = FindFirst( value, achar(0) )
    if ( status > 0 ) value( status: ) = ' '
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset ' // trim(name), &
      & MLSFile=MLSFile )
    call trace_end ( cond=.false. )

  end subroutine LoadFromHDF5DS_charscalar

  ! ------------------------------------  LoadFromHDF5DS_chararr1  -----
  subroutine LoadFromHDF5DS_chararr1 ( MLSFile, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(out) :: VALUE(:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinates of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    include 'LoadFromHDF5DS_chararr.f9h'

  end subroutine LoadFromHDF5DS_chararr1

  ! ------------------------------------  LoadFromHDF5DS_chararr2  -----
  subroutine LoadFromHDF5DS_chararr2 ( MLSFile, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(out) :: VALUE(:,:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    include 'LoadFromHDF5DS_chararr.f9h'

  end subroutine LoadFromHDF5DS_chararr2

  ! ------------------------------------  LoadFromHDF5DS_chararr3  -----
  subroutine LoadFromHDF5DS_chararr3 ( MLSFile, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(out) :: VALUE(:,:,:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    include 'LoadFromHDF5DS_chararr.f9h'

  end subroutine LoadFromHDF5DS_chararr3

  ! -------------------------------------  LoadFromHDF5DS_intarr1  -----
  subroutine LoadFromHDF5DS_intarr1 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Integer ! HDF type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(out) :: VALUE(:)      ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'intarr1' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_intarr1

  ! -------------------------------------  LoadFromHDF5DS_intarr2  -----
  subroutine LoadFromHDF5DS_intarr2 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Integer ! HDF type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(out) :: VALUE(:,:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'intarr2' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_intarr2

  ! -------------------------------------  LoadFromHDF5DS_intarr3  -----
  subroutine LoadFromHDF5DS_intarr3 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Integer ! HDF type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(out) :: VALUE(:,:,:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'intarr3' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_intarr3

  ! -------------------------------------  LoadFromHDF5DS_intarr4  -----
  subroutine LoadFromHDF5DS_intarr4 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Integer ! HDF type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME  ! Name for this dataset
    integer, intent(out) :: VALUE(:,:,:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'intarr4' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_intarr4

  ! -------------------------------------  LoadFromHDF5DS_logarr1  -----
  subroutine LoadFromHDF5DS_logarr1 ( MLSFile, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    logical, intent(out) :: VALUE(:)      ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    character :: MyValue(size(value))   ! 'F' = false, 'T' = true

    ! Executable code
    call trace_begin ( me, 'LoadFromHDF5DS_logarr1', cond=.false. )
    call LoadFromHDF5DS ( MLSFile, name, myValue, start, count, stride, block )
    value = myValue == 'T'
    call trace_end ( cond=.false. )
  end subroutine LoadFromHDF5DS_logarr1

  ! -------------------------------------  LoadFromHDF5DS_dblarr1  -----
  subroutine LoadFromHDF5DS_dblarr1 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Double ! hdf type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME     ! Name for this dataset
    double precision, intent(out) :: VALUE(:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'dblarr1' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_dblarr1

  ! -------------------------------------  LoadFromHDF5DS_dblarr2  -----
  subroutine LoadFromHDF5DS_dblarr2 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Double ! hdf type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME       ! Name for this dataset
    double precision, intent(out) :: VALUE(:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'dblarr2' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_dblarr2

  ! -------------------------------------  LoadFromHDF5DS_dblarr3  -----
  subroutine LoadFromHDF5DS_dblarr3 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Double ! hdf type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(out) :: VALUE(:,:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'dblarr3' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_dblarr3

  ! -------------------------------------  LoadFromHDF5DS_dblarr4  -----
  subroutine LoadFromHDF5DS_dblarr4 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Double ! hdf type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(out) :: VALUE(:,:,:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'dblarr4' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_dblarr4

  ! ------------------------------------  LoadFromHDF5DS_snglarr1  -----
  subroutine LoadFromHDF5DS_snglarr1 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Real ! hdf type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(out) :: VALUE(:)         ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'snglarr1' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_snglarr1

  ! ------------------------------------  LoadFromHDF5DS_snglarr2  -----
  subroutine LoadFromHDF5DS_snglarr2 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Real ! hdf type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(out) :: VALUE(:,:)       ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'snglarr2' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_snglarr2

  ! ------------------------------------  LoadFromHDF5DS_snglarr3  -----
  subroutine LoadFromHDF5DS_snglarr3 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Real ! hdf type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(out) :: VALUE(:,:,:)     ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'snglarr3' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_snglarr3

  ! ------------------------------------  LoadFromHDF5DS_snglarr4  -----
  subroutine LoadFromHDF5DS_snglarr4 ( MLSFile, name, value, &
    & start, count, stride, block, MissingOK )
    ! This routine loads a predefined array with values from a DS
    use H5Global, only: H5T =>  H5T_Native_Real ! hdf type
    type (MLSFile_T)   :: MLSFile
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(out) :: VALUE(:,:,:,:)   ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    logical, intent(in), optional :: MissingOK

    character(len=*), parameter :: Sfx = 'snglarr4' ! Type and rank in subr name

    include 'LoadFromHDF5DS.f9h'

  end subroutine LoadFromHDF5DS_snglarr4

  ! ---------------------------------  LoadAllocFromHDF5DS_chararr1  -----
  subroutine LoadAllocFromHDF5DS_chararr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), allocatable :: VALUE(:) ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: LB, UB
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: STATUS                   ! Flag from HDF5
    integer :: SPACEID                  ! ID of dataspace
    integer(kind=Size_t) :: STRINGSIZE               ! String size
    integer :: STRINGTYPE               ! String type

    ! Executable code
    call trace_begin ( me, 'LoadAllocFromHDF5DS_chararr1', cond=.false. )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dGet_type_f ( setID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for 1D char array ' // trim(name) )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get string size for 1D char array ' // trim(name) )
    if ( stringSize > len(value) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for 1D char array ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_DS_shape ( spaceID, shp, name )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + shp(1)
    call allocate_test ( value, ub, name, moduleName, lowBound=lb )
    call h5dread_f ( setID, stringtype, value, (/ shp(1), ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID, stringType=stringType )
    call trace_end ( cond=.false. )
  end subroutine LoadAllocFromHDF5DS_chararr1

  ! ---------------------------------  LoadAllocFromHDF5DS_chararr2  -----
  subroutine LoadAllocFromHDF5DS_chararr2 ( locID, name, value )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), allocatable :: VALUE(:,:) ! The array itself

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer(hsize_t) :: SHP(2)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=Size_t) :: STRINGSIZE               ! String size
    integer :: STRINGTYPE               ! String type

    ! Executable code
    call trace_begin ( me, 'LoadAllocFromHDF5DS_chararr2', cond=.false. )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dGet_type_f ( setID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for 1D char array ' // trim(name) )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for 1D char array ' // trim(name) )
    if ( stringSize > len(value) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for 1D char array ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_DS_shape ( spaceID, shp, name )
    call allocate_test ( value, int(shp(1)), int(shp(2)), name, moduleName )
    call h5dread_f ( setID, stringtype, value, (/ shp, ones(1:5) /), status )
    call finishLoad ( name, status, spaceID, setID, stringType=stringType )
    call trace_end ( cond=.false. )
  end subroutine LoadAllocFromHDF5DS_chararr2

  ! ----------------------------------  LoadAllocFromHDF5DS_intarr1  -----
  subroutine LoadAllocFromHDF5DS_intarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, allocatable :: VALUE(:)          ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: LB, UB
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'LoadAllocFromHDF5DS_intarr1', cond=.false. )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + shp(1)
    call allocate_test ( value, ub, name, moduleName, lowBound=lb )
    call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
      & (/ shp, ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID )
    call trace_end ( cond=.false. )
  end subroutine LoadAllocFromHDF5DS_intarr1

  ! ----------------------------------  LoadAllocFromHDF5DS_intarr2  -----
  subroutine LoadAllocFromHDF5DS_intarr2 ( locID, name, value )
    ! This routine allocates a allocatable array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Integer ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, allocatable :: VALUE(:,:)        ! The array itself
    character(len=*), parameter :: Sfx = 'intarr2'

    include 'LoadAllocFromHDF5DS.f9h'

  end subroutine LoadAllocFromHDF5DS_intarr2

  ! ----------------------------------  LoadAllocFromHDF5DS_intarr3  -----
  subroutine LoadAllocFromHDF5DS_intarr3 ( locID, name, value )
    ! This routine allocates a allocatable array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Integer ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, allocatable :: VALUE(:,:,:)      ! The array itself
    character(len=*), parameter :: Sfx = 'intarr3'

    include 'LoadAllocFromHDF5DS.f9h'

  end subroutine LoadAllocFromHDF5DS_intarr3

  ! ----------------------------------  LoadAllocFromHDF5DS_intarr4  -----
  subroutine LoadAllocFromHDF5DS_intarr4 ( locID, name, value )
    ! This routine allocates a allocatable array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Integer ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, allocatable :: VALUE(:,:,:,:)    ! The array itself
    character(len=*), parameter :: Sfx = 'intarr4'

    include 'LoadAllocFromHDF5DS.f9h'

  end subroutine LoadAllocFromHDF5DS_intarr4

  ! ----------------------------------  LoadAllocFromHDF5DS_logarr1  -----
  subroutine LoadAllocFromHDF5DS_logarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    logical, allocatable :: VALUE(:)          ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: LB, UB
    integer :: Me = -1                  ! String index for trace cacheing
    character, allocatable :: MyValue(:)      ! 'F' = false, 'T' = true

    ! Executable code
    call trace_begin ( me, 'LoadAllocFromHDF5DS_logarr1', cond=.false. )
    call LoadAllocFromHDF5DS ( locID, name, myValue )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + size(myValue)
    call allocate_test ( value, ub, name, moduleName, lowBound=lb )
    value = myValue == 'T'
    call deallocate_test ( myValue, name, moduleName )
    call trace_end ( cond=.false. )
  end subroutine LoadAllocFromHDF5DS_logarr1

  ! ----------------------------------  LoadAllocFromHDF5DS_dblarr1  -----
  subroutine LoadAllocFromHDF5DS_dblarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, allocatable :: VALUE(:) ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: LB, UB
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'LoadAllocFromHDF5DS_dblarr1', cond=.false. )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + shp(1)
    call allocate_test ( value, ub, name, moduleName, lowBound=lb )
    call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & (/ shp, ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID )
    call trace_end ( cond=.false. )
  end subroutine LoadAllocFromHDF5DS_dblarr1

  ! ----------------------------------  LoadAllocFromHDF5DS_dblarr2  -----
  subroutine LoadAllocFromHDF5DS_dblarr2 ( locID, name, value )
    ! This routine allocates a allocatable array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Double ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, allocatable :: VALUE(:,:) ! The array itself
    character(len=*), parameter :: Sfx = 'dblarr2'

    include 'LoadAllocFromHDF5DS.f9h'

  end subroutine LoadAllocFromHDF5DS_dblarr2

  ! ----------------------------------  LoadAllocFromHDF5DS_dblarr3  -----
  subroutine LoadAllocFromHDF5DS_dblarr3 ( locID, name, value )
    ! This routine allocates a allocatable array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Double ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, allocatable :: VALUE(:,:,:) ! The array itself
    character(len=*), parameter :: Sfx = 'dblarr3'

    include 'LoadAllocFromHDF5DS.f9h'

  end subroutine LoadAllocFromHDF5DS_dblarr3

  ! ----------------------------------  LoadAllocFromHDF5DS_dblarr4  -----
  subroutine LoadAllocFromHDF5DS_dblarr4 ( locID, name, value )
    ! This routine allocates a allocatable array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Double ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, allocatable :: VALUE(:,:,:,:) ! The array itself
    character(len=*), parameter :: Sfx = 'dblarr4'

    include 'LoadAllocFromHDF5DS.f9h'

  end subroutine LoadAllocFromHDF5DS_dblarr4

  ! ---------------------------------  LoadAllocFromHDF5DS_snglarr1  -----
  subroutine LoadAllocFromHDF5DS_snglarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, allocatable :: VALUE(:)             ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: LB, UB
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'LoadAllocFromHDF5DS_snglarr1', cond=.false. )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + shp(1)
    call allocate_test ( value, ub, name, moduleName, lowBound=lb )
    call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
      & (/ shp, ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID )
    call trace_end ( cond=.false. )
  end subroutine LoadAllocFromHDF5DS_snglarr1

  ! ---------------------------------  LoadAllocFromHDF5DS_snglarr2  -----
  subroutine LoadAllocFromHDF5DS_snglarr2 ( locID, name, value )
    ! This routine allocates a allocatable array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Real ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, allocatable :: VALUE(:,:)           ! The array itself
    character(len=*), parameter :: Sfx = 'snglarr2'

    include 'LoadAllocFromHDF5DS.f9h'

  end subroutine LoadAllocFromHDF5DS_snglarr2

  ! ---------------------------------  LoadAllocFromHDF5DS_snglarr3  -----
  subroutine LoadAllocFromHDF5DS_snglarr3 ( locID, name, value )
    ! This routine allocates a allocatable array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Real ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, allocatable :: VALUE(:,:,:)         ! The array itself
    character(len=*), parameter :: Sfx = 'snglarr3'

    include 'LoadAllocFromHDF5DS.f9h'

  end subroutine LoadAllocFromHDF5DS_snglarr3

  ! ---------------------------------  LoadAllocFromHDF5DS_snglarr4  -----
  subroutine LoadAllocFromHDF5DS_snglarr4 ( locID, name, value )
    ! This routine allocates a allocatable array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Real ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, allocatable :: VALUE(:,:,:,:)       ! The array itself
    character(len=*), parameter :: Sfx = 'snglarr4'

    include 'LoadAllocFromHDF5DS.f9h'

  end subroutine LoadAllocFromHDF5DS_snglarr4

  ! ---------------------------------  LoadPtrFromHDF5DS_chararr1  -----
  subroutine LoadPtrFromHDF5DS_chararr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), pointer :: VALUE(:) ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: LB, UB
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: STATUS                   ! Flag from HDF5
    integer :: SPACEID                  ! ID of dataspace
    integer(kind=Size_t) :: STRINGSIZE               ! String size
    integer :: STRINGTYPE               ! String type

    ! Executable code
    call trace_begin ( me, 'LoadPtrFromHDF5DS_chararr1', cond=.false. )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dGet_type_f ( setID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for 1D char array ' // trim(name) )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get string size for 1D char array ' // trim(name) )
    if ( stringSize > len(value) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for 1D char array ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_DS_shape ( spaceID, shp, name )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + shp(1)
    call allocate_test ( value, ub, name, moduleName, lowBound=lb )
    call h5dread_f ( setID, stringtype, value, (/ shp(1), ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID, stringType=stringType )
    call trace_end ( cond=.false. )
  end subroutine LoadPtrFromHDF5DS_chararr1

  ! ---------------------------------  LoadPtrFromHDF5DS_chararr2  -----
  subroutine LoadPtrFromHDF5DS_chararr2 ( locID, name, value )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), pointer :: VALUE(:,:) ! The array itself

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer(hsize_t) :: SHP(2)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=Size_t) :: STRINGSIZE               ! String size
    integer :: STRINGTYPE               ! String type

    ! Executable code
    call trace_begin ( me, 'LoadPtrFromHDF5DS_chararr2', cond=.false. )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dGet_type_f ( setID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for 1D char array ' // trim(name) )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for 1D char array ' // trim(name) )
    if ( stringSize > len(value) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for 1D char array ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_DS_shape ( spaceID, shp, name )
    call allocate_test ( value, int(shp(1)), int(shp(2)), name, moduleName )
    call h5dread_f ( setID, stringtype, value, (/ shp, ones(1:5) /), status )
    call finishLoad ( name, status, spaceID, setID, stringType=stringType )
    call trace_end ( cond=.false. )
  end subroutine LoadPtrFromHDF5DS_chararr2

  ! ----------------------------------  LoadPtrFromHDF5DS_intarr1  -----
  subroutine LoadPtrFromHDF5DS_intarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, pointer :: VALUE(:)          ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: LB, UB
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'LoadPtrFromHDF5DS_intarr1', cond=.false. )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + shp(1)
    call allocate_test ( value, ub, name, moduleName, lowBound=lb )
    call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
      & (/ shp, ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID )
    call trace_end ( cond=.false. )
  end subroutine LoadPtrFromHDF5DS_intarr1

  ! ----------------------------------  LoadPtrFromHDF5DS_intarr2  -----
  subroutine LoadPtrFromHDF5DS_intarr2 ( locID, name, value )
    ! This routine allocates a pointer array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Integer ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, pointer :: VALUE(:,:)        ! The array itself
    character(len=*), parameter :: Sfx = 'intarr2'

    include 'LoadPtrFromHDF5DS.f9h'

  end subroutine LoadPtrFromHDF5DS_intarr2

  ! ----------------------------------  LoadPtrFromHDF5DS_intarr3  -----
  subroutine LoadPtrFromHDF5DS_intarr3 ( locID, name, value )
    ! This routine allocates a pointer array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Integer ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, pointer :: VALUE(:,:,:)      ! The array itself
    character(len=*), parameter :: Sfx = 'intarr3'

    include 'LoadPtrFromHDF5DS.f9h'

  end subroutine LoadPtrFromHDF5DS_intarr3

  ! ----------------------------------  LoadPtrFromHDF5DS_intarr4  -----
  subroutine LoadPtrFromHDF5DS_intarr4 ( locID, name, value )
    ! This routine allocates a pointer array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Integer ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, pointer :: VALUE(:,:,:,:)    ! The array itself
    character(len=*), parameter :: Sfx = 'intarr4'

    include 'LoadPtrFromHDF5DS.f9h'

  end subroutine LoadPtrFromHDF5DS_intarr4

  ! ----------------------------------  LoadPtrFromHDF5DS_logarr1  -----
  subroutine LoadPtrFromHDF5DS_logarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    logical, pointer :: VALUE(:)          ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: LB, UB
    integer :: Me = -1                  ! String index for trace cacheing
    character, pointer :: MyValue(:)      ! 'F' = false, 'T' = true

    ! Executable code
    call trace_begin ( me, 'LoadPtrFromHDF5DS_logarr1', cond=.false. )
    nullify ( myValue )
    call LoadPtrFromHDF5DS ( locID, name, myValue )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + size(myValue)
    call allocate_test ( value, ub, name, moduleName, lowBound=lb )
    value = myValue == 'T'
    call deallocate_test ( myValue, name, moduleName )
    call trace_end ( cond=.false. )
  end subroutine LoadPtrFromHDF5DS_logarr1

  ! ----------------------------------  LoadPtrFromHDF5DS_dblarr1  -----
  subroutine LoadPtrFromHDF5DS_dblarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, pointer :: VALUE(:) ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: LB, UB
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'LoadPtrFromHDF5DS_dblarr1', cond=.false. )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + shp(1)
    call allocate_test ( value, ub, name, moduleName, lowBound=lb )
    call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & (/ shp, ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID )
    call trace_end ( cond=.false. )
  end subroutine LoadPtrFromHDF5DS_dblarr1

  ! ----------------------------------  LoadPtrFromHDF5DS_dblarr2  -----
  subroutine LoadPtrFromHDF5DS_dblarr2 ( locID, name, value )
    ! This routine allocates a pointer array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Double ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, pointer :: VALUE(:,:) ! The array itself
    character(len=*), parameter :: Sfx = 'dblarr2'

    include 'LoadPtrFromHDF5DS.f9h'

  end subroutine LoadPtrFromHDF5DS_dblarr2

  ! ----------------------------------  LoadPtrFromHDF5DS_dblarr3  -----
  subroutine LoadPtrFromHDF5DS_dblarr3 ( locID, name, value )
    ! This routine allocates a pointer array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Double ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, pointer :: VALUE(:,:,:) ! The array itself
    character(len=*), parameter :: Sfx = 'dblarr3'

    include 'LoadPtrFromHDF5DS.f9h'

  end subroutine LoadPtrFromHDF5DS_dblarr3

  ! ----------------------------------  LoadPtrFromHDF5DS_dblarr4  -----
  subroutine LoadPtrFromHDF5DS_dblarr4 ( locID, name, value )
    ! This routine allocates a pointer array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Double ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, pointer :: VALUE(:,:,:,:) ! The array itself
    character(len=*), parameter :: Sfx = 'dblarr4'

    include 'LoadPtrFromHDF5DS.f9h'

  end subroutine LoadPtrFromHDF5DS_dblarr4

  ! ---------------------------------  LoadPtrFromHDF5DS_snglarr1  -----
  subroutine LoadPtrFromHDF5DS_snglarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, pointer :: VALUE(:)             ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: LB, UB
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: SETID                    ! ID of dataset
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call trace_begin ( me, 'LoadPtrFromHDF5DS_snglarr1', cond=.false. )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + shp(1)
    call allocate_test ( value, ub, name, moduleName, lowBound=lb )
    call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
      & (/ shp, ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID )
    call trace_end ( cond=.false. )
  end subroutine LoadPtrFromHDF5DS_snglarr1

  ! ---------------------------------  LoadPtrFromHDF5DS_snglarr2  -----
  subroutine LoadPtrFromHDF5DS_snglarr2 ( locID, name, value )
    ! This routine allocates a pointer array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Real ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, pointer :: VALUE(:,:)           ! The array itself
    character(len=*), parameter :: Sfx = 'snglarr2'

    include 'LoadPtrFromHDF5DS.f9h'

  end subroutine LoadPtrFromHDF5DS_snglarr2

  ! ---------------------------------  LoadPtrFromHDF5DS_snglarr3  -----
  subroutine LoadPtrFromHDF5DS_snglarr3 ( locID, name, value )
    ! This routine allocates a pointer array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Real ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, pointer :: VALUE(:,:,:)         ! The array itself
    character(len=*), parameter :: Sfx = 'snglarr3'

    include 'LoadPtrFromHDF5DS.f9h'

  end subroutine LoadPtrFromHDF5DS_snglarr3

  ! ---------------------------------  LoadPtrFromHDF5DS_snglarr4  -----
  subroutine LoadPtrFromHDF5DS_snglarr4 ( locID, name, value )
    ! This routine allocates a pointer array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    use H5Global, only: H5T =>  H5T_Native_Real ! hdf type
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, pointer :: VALUE(:,:,:,:)       ! The array itself
    character(len=*), parameter :: Sfx = 'snglarr4'

    include 'LoadPtrFromHDF5DS.f9h'

  end subroutine LoadPtrFromHDF5DS_snglarr4

  ! -----------------------------------  ReadLitIndexFromHDF5Attr  -----
  subroutine ReadLitIndexFromHDF5Attr ( itemID, name, index )
    use MoreTree, only: GetLitIndexFromString
    ! Dummy arguments
    integer, intent(in) :: ITEMID        ! Group etc. to make attr. for
    character(len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(out) :: INDEX        ! String index

    ! Local variables
    character (len=1024) :: LINE
    integer :: L
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'ReadLitIndexFromHDF5Attr', cond=.false. )
    call GetHDF5Attribute ( itemID, name, line )
    l = len_trim(line)
    if ( l > 0 ) then
      index = GetLitIndexFromString ( line(:l) )
    else
      index = 0
    end if
    call trace_end ( cond=.false. )
  end subroutine ReadLitIndexFromHDF5Attr

  ! --------------------------------  ReadStringIndexFromHDF5Attr  -----
  subroutine ReadStringIndexFromHDF5Attr ( itemID, name, index )
    use MoreTree, only: GetStringIndexFromString
    ! Dummy arguments
    integer, intent(in) :: ITEMID        ! Group etc. to make attr. for
    character(len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(out) :: INDEX        ! String index

    ! Local variables
    character (len=1024) :: LINE
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'ReadStringIndexFromHDF5Attr', cond=.false. )
    call GetHDF5Attribute ( itemID, name, line )
    if ( len_trim ( line ) > 0 ) then
      index = GetStringIndexFromString ( trim(line) )
    else
      index = 0
    end if
    call trace_end ( cond=.false. )
  end subroutine ReadStringIndexFromHDF5Attr

  ! -------------------------------  WriteLitIndexAsHDF5Attribute  -----
  subroutine WriteLitIndexAsHDF5Attribute ( itemID, name, index )
    use Intrinsic, only: Lit_Indices
    integer, intent(in) :: ITEMID        ! Group etc. to make attr. for
    character(len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(in) :: INDEX         ! String index

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'WriteLitIndexAsHDF5Attribute', cond=.false. )
    if ( index == 0 ) then
      call MakeHDF5Attribute ( itemID, name, '' )
    else
      call WriteStringIndexAsHDF5Attribute ( itemID, name, lit_indices ( index ) )
    end if
    call trace_end ( cond=.false. )
  end subroutine WriteLitIndexAsHDF5Attribute

  ! ----------------------------  WriteStringIndexAsHDF5Attribute  -----
  subroutine WriteStringIndexAsHDF5Attribute ( itemID, name, index )
    use String_table, only: Get_string
    integer, intent(in) :: ITEMID        ! Group etc. to make attr. for
    character(len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(in) :: INDEX         ! String index

    ! Local variables
    character(len=1024) :: LINE
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'WriteStringIndexAsHDF5Attribute', cond=.false. )
    if ( index == 0 ) then
      call MakeHDF5Attribute ( itemID, name, '' )
    else
      call get_string ( index, line, strip=.true., noError=.true. )
      call MakeHDF5Attribute ( itemID, name, trim(line) )
    end if
    call trace_end ( cond=.false. )
  end subroutine WriteStringIndexAsHDF5Attribute

! ======================= Private Procedures ===========================
! --------------------------------------------  AreThe2TypesEqual  -----
  logical function AreThe2TypesEqual ( type1, type2 )
    ! This routine returns true if the two datatypes are "the same"
    integer, intent(in) :: type1
    integer, intent(in) :: type2
    ! Local variables
    integer :: status
    ! Initialize in case something goes wrong with hdf5 routine
    AreThe2TypesEqual = .false.
    call h5tEqual_f ( type1, type2, AreThe2TypesEqual, status )
    if ( status /= 0 ) AreThe2TypesEqual = .false.
  end function AreThe2TypesEqual

! ------------------------------------------------  Check_for_fit  -----
  subroutine Check_for_fit ( spaceID, value_dims, name )
  ! Checks that dataspace will fit into values before LoadFromHDF5DS
    integer, intent(in)               :: spaceID
    integer(hsize_t), dimension(:), intent(in) :: value_dims
    character (len=*), intent(in)     :: NAME ! Name for this dataset
    integer(hsize_t), dimension(size(value_dims)) :: dims
    integer                           :: i
    call get_ds_shape ( spaceID, dims, name )
    if ( any ( dims > value_dims ) ) &
      & call my_message ( MLSMSG_Error, ModuleName // '%Check_for_fit', &
      & 'Dataspace too large for destination value of ' // trim(name) , &
      & 'dims(space), dims(value)', (/ (int(dims(i)), int(value_dims(i)), i=1, size(dims)) /), &
      & no_pairs=.true. )
  end subroutine Check_for_fit

! -----------------------------------------------  CreateSpaceSet  -----
  subroutine CreateSpaceSet ( locID, name, maxdims, datatype, &
    & spaceID, setID, status, adding_to, cFill, dFill, iFill, rFill, chunk_dims )
  ! Creates or optionally attaches dataspace and dataset
    integer, intent(in)               :: locID
    character (len=*), intent(in)     :: NAME ! Name for this dataset
    integer(hsize_t), dimension(:), intent(in) :: maxdims
    integer, intent(in)               :: datatype
    integer, intent(out)              :: spaceID
    integer, intent(out)              :: setID
    integer, intent(out)              :: status
    logical, optional, intent(in)     :: adding_to
    character(len=*), optional, intent(in) :: cFill
    double precision, optional, intent(in) :: dFill
    integer, optional, intent(in)     :: iFill
    real, optional, intent(in)        :: rFill
    integer, dimension(:), optional, intent(in) :: chunk_dims

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_adding_to
    logical :: my_fill
    integer(hsize_t), dimension(size(maxdims,1)) :: my_chunkdims, my_maxdims
    integer :: Rank
    character(len=1) :: cFilled
    double precision :: dFilled
    integer :: iFilled
    real :: rFilled

    ! Executable
    call trace_begin ( me, 'CreateSpaceSet', cond=.false. )
    my_adding_to = .false.
    if ( present(adding_to) ) my_adding_to = adding_to
    my_chunkdims = maxdims
    rank = size(maxdims,1)
    my_chunkdims(rank) = 1
    if ( present(chunk_dims) ) my_chunkdims = chunk_dims
    my_fill = present(cFill) .or. present(dFill) .or. present(iFill) .or. &
      & present(rFill) .or. present(chunk_dims)
    if ( my_adding_to ) then
      call h5dopen_f ( locID, trim(name), setID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open dataset for ' // trim(name) )
      call h5dget_space_f ( setID, spaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dataspace for ' // trim(name) )
      if ( my_fill ) then
        ! print *, 'dataspace before extending ', spaceID
        my_maxdims = maxdims
        call h5dextend_f ( setID, my_maxdims, status )
        call h5dget_space_f ( setID, spaceID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to get dataspace for ' // trim(name) )
        ! print *, 'dataspace after extending ', spaceID
        call h5dget_create_plist_f ( setID, cparms, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to get property list for ' // trim(name) )
        call h5pget_chunk_f ( cparms, rank, my_chunkdims, status )
        ! print *, '(read) rank, my_chunkdims, status'
        ! print *, rank, my_chunkdims, status
        if ( present(cFill) ) then
          call h5pget_fill_value_f ( cparms, H5T_NATIVE_CHARACTER, cFilled, status )
        else if ( present(dFill) ) then
          call h5pget_fill_value_f ( cparms, datatype, dFilled, status )
        else if ( present(iFill) ) then
          call h5pget_fill_value_f ( cparms, datatype, iFilled, status )
          ! print *, 'cparms, datatype, iFilled, status'
          ! print *, cparms, datatype, iFilled, status
        else if ( present(rFill) ) then
          call h5pget_fill_value_f ( cparms, datatype, rFilled, status )
        end if
      end if
    else
      ! Create the dataset
      if ( my_fill ) then
        my_maxdims = maxdims
        my_maxdims(rank) = H5S_UNLIMITED_F
        call h5sCreate_simple_f ( rank, int(maxdims,hSize_T), spaceID, status, &
          & my_maxdims )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to create extensible dataspace ' // trim(name) )
        call h5pcreate_f ( H5P_DATASET_CREATE_F, cparms, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to create property list for ' // trim(name) )
        call h5pset_chunk_f ( cparms, rank, my_chunkdims, status )
        ! print *, 'rank, my_chunkdims, status'
        ! print *, rank, my_chunkdims, status
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to set chunking for ' // trim(name) )
        if ( present(cFill) ) then
          call h5pset_fill_value_f ( cparms, H5T_NATIVE_CHARACTER, cFill, status )
        else if ( present(dFill) ) then
          call h5pset_fill_value_f ( cparms, datatype, dFill, status )
        else if ( present(iFill) ) then
          call h5pset_fill_value_f ( cparms, datatype, iFill, status )
          ! print *, 'cparms, datatype, iFill, status'
          ! print *, cparms, datatype, iFill, status
        else if ( present(rFill) ) then
          call h5pset_fill_value_f (cparms, datatype, rFill, status)
        end if
        if ( status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
          &  "Unable to set Fill value for " // trim (name) )
        call h5dCreate_f ( locID, trim(name), datatype, spaceID, setID, &
          & status, cparms )
        call h5dextend_f ( setID, int(maxdims,hSize_T), status )
        call h5dget_space_f( setID, spaceID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to get dataspace for ' // trim(name) )
      else
        call h5sCreate_simple_f ( rank, int(maxdims,hSize_T), spaceID, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to create dataspace ' // trim(name) )
        call h5dCreate_f ( locID, trim(name), datatype, spaceID, setID, &
          & status )
      end if
    end if
    call trace_end ( cond=.false. )
  end subroutine CreateSpaceSet

! ---------------------------------------------------  Dump_space  -----
  subroutine Dump_space ( spaceID )
  ! Dumps dataspace info
    integer, intent(in)               :: spaceID
    integer                           :: rank
    integer(hsize_t), dimension(7)    :: dims, maxdims
    integer                           :: status
    logical                           :: is_simple
    call h5sis_simple_f ( spaceID, is_simple, status )
    if ( DEEBUG ) print *, 'is_simple ', is_simple
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( DEEBUG ) print *, 'rank ', rank
    call h5sget_simple_extent_dims_f ( spaceID, dims(1:rank), maxdims(1:rank), &
      &  status )
    if ( status /= rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dims for dumping space' )
    if ( DEEBUG ) print *, 'dims ', dims(1:rank)
    if ( DEEBUG ) print *, 'maxdims ', maxdims(1:rank)
    ! call h5soffset_simple_f ( spaceID, offset(1:rank), &
    !  &  status )
    ! print *, 'offset ', offset(1:rank)
  end subroutine Dump_space

! ---------------------------------------------------  FinishLoad  -----
  subroutine FinishLoad ( Name, Status, SpaceID, SetID, &
    & MemspaceID, StringType, MLSFile )
    ! Checks status and closes stuff to finish Load...
    character(len=*), intent(in) :: Name    ! of the dataset
    integer, intent(inout) :: Status        ! from last read operation
    integer, intent(in) :: SpaceID          ! dataspace ID
    integer, intent(in) :: SetID            ! dataset ID
    integer, intent(in), optional :: MemspaceID ! dataspace ID
    integer, intent(in), optional :: StringType ! stringtype ID
    type (MLSFile_T), optional   :: MLSFile

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable
    call trace_begin ( me, 'FinishLoad', cond=.false. )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset ' // trim(name), &
      & MLSFile=MLSFile )
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset ' // trim(name), &
      & MLSFile=MLSFile )
    if ( present(memspaceID) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset ' // trim(name), &
        & MLSFile=MLSFile )
    end if
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset ' // trim(name), &
      & MLSFile=MLSFile )
    if ( present(stringType) ) then
      call h5tClose_f ( stringType, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close string type for ' // trim(name), &
        & MLSFile=MLSFile )
    end if
    call trace_end ( cond=.false. )
  end subroutine FinishLoad

! ---------------------------------------------  FinishMakeAttrib  -----
  subroutine FinishMakeAttrib ( Name, Status, AttrID, DSID, StringType )
  ! Check status and close stuff for MakeHDF5Attribute_....
    character(len=*), intent(in) :: Name        ! of the attrib
    integer, intent(inout) :: Status            ! from writing the attrib
    integer, intent(in) :: AttrID               ! attrib ID
    integer, intent(in) :: DSID                 ! dataset ID
    integer, intent(in), optional :: StringType ! stringType ID

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable
    call trace_begin ( me, 'FinishMakeHDF5Attrib', cond=.false. )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write attribute ' // trim(name) )
    ! Finish off
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
    call h5sClose_f ( dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute dataspace ' // trim(name) )
    if ( present(stringType) ) then
      call h5tClose_f ( stringType, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close stringType ' // trim(name) )
    end if
    call trace_end ( cond=.false. )
  end subroutine FinishMakeAttrib

! -------------------------------------------------  FinishSaveDS  -----
  subroutine FinishSaveDS ( Name, Status, SetID, SpaceID, StringType, &
    & MemspaceID )
  ! Check status and close stuff for SaveAsHDF5DS_...
    character(len=*), intent(in) :: Name        ! of the DS
    integer, intent(inout) :: Status            ! from writing the DS
    integer, intent(in) :: SetID                ! dataset ID
    integer, intent(in) :: SpaceID              ! dataspace ID
    integer, intent(in), optional :: StringType ! stringtype ID
    integer, intent(in), optional :: MemspaceID ! arrayspace ID

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable
    call trace_begin ( me, 'FinishSaveDS', cond=.false. )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write to dataset for ' // trim(name) )
    ! Close things
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for ' // trim(name) )
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for ' // trim(name) )
    if ( present(memspaceID) ) then
      call h5sClose_F ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close mem dataspace for ' // trim(name) )
    end if
    if ( present(stringType) ) then
      call h5tClose_f ( stringtype, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close stringtype for ' // trim(name) )
    end if
    call trace_end ( cond=.false. )
  end subroutine FinishSaveDS

! -------------------------------------------------  Get_DS_Shape  -----
  subroutine Get_DS_Shape ( spaceID, dims, name )
  ! Checks that dataspace has the same rank as size(dims), then gets
  ! its dimensions.
    integer, intent(in)               :: spaceID
    integer(hsize_t), dimension(:), intent(out) :: dims
    character (len=*), intent(in)     :: NAME ! Name for this dataset
    integer                           :: rank
    integer                           :: value_rank
    integer(hsize_t)                  :: maxdims(size(dims))
    integer                           :: status

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable
    call trace_begin ( me, 'Get_DS_Shape', cond=.false. )
    value_rank = size(dims)
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( status /= 0 )  call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get rank for dataset ' // trim(name) )
    if ( rank /= value_rank ) call my_message ( MLSMSG_Error, ModuleName, &
      & 'Inconsistent rank for dataset ' // trim(name) , &
      & 'rank(space), rank(values)', (/rank, value_rank/) )
    call h5sget_simple_extent_dims_f ( spaceID, dims, maxdims, status )
    if ( status /= rank ) call my_message ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dimension information for dataset ' // trim(name) , &
      & 'rank(space), h5s status', (/rank, status/) )
    call trace_end ( cond=.false. )
  end subroutine Get_DS_Shape

! ---------------------------------------------------  MLS_eSet_auto  -----
  subroutine MLS_eSet_auto ( value, status )
  ! Calls h5eSet_auto_f to turn on (if value==1) or off (if value==0)
  ! error messages printed by the hdf5 library
  ! Unless hdfVerbose, in which case always print
    integer, intent(in)                 :: value
    integer, intent(out)                :: status
    !
    integer                             :: myValue
    ! Executable
    myValue = value
    if ( hdfVerbose ) myValue = 1
    call h5eSet_auto_f ( myvalue, status )
  end subroutine MLS_eSet_auto

! ---------------------------------------------------  MLS_extend  -----
  subroutine MLS_extend ( setID, newCount, start, dataSpaceID )
  ! Checks whether we need to extend setID to accommodate newDims
  ! plus any offsets in start array
  ! Does the extending if necessary
    integer, intent(in)                :: setID
    integer(hsize_t), dimension(:), intent(in)  :: newCount
    integer, dimension(:), optional, intent(in) :: start
    integer, optional, intent(inout) :: dataSpaceID
                                 ! Starting coordinatess of hyperslab
  ! Local variables
    integer(hsize_t), dimension(7)    :: dims, maxdims, my_start
    integer                           :: i
    logical                           :: is_simple
    logical                           :: itFits
    integer :: Me = -1                ! String index for trace cacheing
    integer                           :: rank
    integer                           :: spaceID
    integer                           :: status

  ! Executable code
    call trace_begin ( me, 'MLS_extend', cond=.false. )
    my_start = 0
    dims = 0
    if ( present(start) ) my_start(1:size(start)) = start
    if ( present(dataSpaceID) ) then
      spaceID = dataSpaceID
    else
      call h5dget_space_f(setID, spaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get data space ID to extend data set' )
    end if
    if ( DEEBUG ) print *, 'spaceID ', spaceID
    call h5sis_simple_f ( spaceID, is_simple, status )
    if ( DEEBUG ) print *, 'is simple? ', is_simple
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get data space rank to extend data set' )
    if ( DEEBUG ) print *, 'rank ', rank
    call h5sget_simple_extent_dims_f ( spaceID, dims, maxdims, status )
    if ( status /= rank .and. cantGetDataspaceDims <= MAXNUMWARNS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Unable to get data space dims to extend data set' )
      cantGetDataspaceDims = cantGetDataspaceDims + 1
      if ( cantGetDataspaceDims > MAXNUMWARNS ) &
        & call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Max no. of warnings reached--suppressing further ones')
    end if
    if ( DEEBUG ) print *, 'dims ', dims(1:rank)
    if ( DEEBUG ) print *, 'maxdims ', maxdims(1:rank)
    if ( DEEBUG ) print *, 'status ', status
    itFits = .true.
    do i = 1, min(rank, size(newCount))
      if ( my_start(i) + newCount(i) > dims(i) ) then
        itFits = .false.
        dims(i) = my_start(i) + newCount(i)
      end if
    enddo
    if ( itFits ) go to 9
    if ( DEEBUG ) print *, 'Need to extend dataset'
    if ( DEEBUG ) print *, '(New dims) ', dims
    call h5dextend_f ( setID, dims, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to extend data set in mls_extend' )
    if ( present(dataspaceID) ) then
      call h5dget_space_f ( setID, dataspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get return dataspaceID in mls_extend' )
    end if
  9 call trace_end ( cond=.false. )
  end subroutine MLS_extend

! ------------------------------------------------  MLS_hyperslab  -----
  subroutine MLS_hyperslab ( spaceID, value_dims, name, memspaceID, &
    & start, count, stride, block )
  ! Sits between LoadFromHDF5DS and h5sselect_hyperslab_f
  ! Restriction:
  ! The optional parameters must be present in 1 of the two patterns below
  ! pattern(1): start, count)
  ! pattern(2): start, count, stride, block)
    integer, intent(in)                :: spaceID
    integer, intent(out)               :: memspaceID
    integer(hsize_t), dimension(:), intent(in)  :: value_dims
    character (len=*), intent(in)     :: NAME ! Name for this dataset
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    integer                           :: value_rank, rank
    integer(hsize_t), dimension(7)    :: dims, maxdims
    integer                           :: status
                                 ! Starting coordinatess of hyperslab
    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    ! logical, parameter :: DeeBug = .true.

    ! Begin execution
    call trace_begin ( me, 'MLS_hyperslab', cond=.false. )
    ! Check that pattern 1 or pattern 2 is satisfied
    status = 0
    if ( present(start) ) status = status + 1
    if ( present(count) ) status = status + 2
    if ( present(stride) ) status = status + 4
    if ( present(block) ) status = status + 8
    if ( status /= 3 .and. status /= 15 ) call my_message ( MLSMSG_Error, ModuleName, &
      & 'Impossible optional parameters pattern for dataset ' // trim(name), &
      & 'status', (/status/) )
    value_rank = size(value_dims)
    if ( DEEBUG ) print *, 'value_rank: ', value_rank
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( DEEBUG ) print *, 'dataspace rank: ', rank
    if ( status /= 0 )  call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get rank for dataset ' // trim(name) )
    call h5sget_simple_extent_dims_f ( spaceID, dims(1:rank), maxdims(1:rank), &
      &  status )
    if ( DEEBUG ) print *, 'dataspace dims: ', dims(1:rank)
    if ( status /= rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dimension information for dataset ' // trim(name) )
    call h5screate_simple_f(value_rank, &
      & int(value_dims(1:value_rank), hsize_t), memspaceID, &
      & status)
    if ( DEEBUG ) print *, 'memspace id: ', memspaceID
    if ( status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create memspace for dataset ' // trim(name) )
    if ( present(stride) ) then
     if ( DEEBUG ) print *, 'trying to select hyperslab: ', &
       & int(start(1:rank), hsize_t), int(count(1:rank), hsize_t), &
       & int(stride(1:rank), hsize_t), int(block(1:rank), hsize_t)
      if ( rank > size(start) ) then
        print *, 'Error in selecting hyperslab--rank too big for size(start)'
        print *, 'rank   ', rank
        print *, 'start  ', start 
        print *, 'count  ', count 
        print *, 'stride ', stride
        print *, 'block  ', block 
      endif
      call h5sselect_hyperslab_f ( spaceID, H5S_SELECT_SET_F, &
        & int(start(1:rank), hsize_t), int(count(1:rank), hsize_t), status, &
        & int(stride(1:rank), hsize_t), int(block(1:rank), hsize_t) )
    else
      if ( DEEBUG ) print *, 'trying to select hyperslab: ', &
        & int(start(1:rank), hsize_t), int(count(1:rank), hsize_t)
      call h5sselect_hyperslab_f ( spaceID, H5S_SELECT_SET_F, &
        & int(start(1:rank), hsize_t), int(count(1:rank), hsize_t), status )
    end if
    if ( status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set hyperslab for dataset ' // trim(name) )
    if ( DEEBUG ) print *, 'Returning memspaceID ', memspaceID
    call trace_end ( cond=.false. )
  end subroutine MLS_hyperslab

! -------------------------------------------  MLS_hyperslab_save  -----
  subroutine MLS_hyperslab_save ( spaceID, start, count, stride, block )
  ! Sits between SaveAsHDF5DS and h5sselect_hyperslab_f
  ! Restriction:
  ! The optional parameters must be present in 1 of the two patterns below
  ! pattern(1): start, count)
  ! pattern(2): start, count, stride, block)
    integer, intent(in)                :: spaceID
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    integer                           :: rank
    integer(hsize_t), dimension(7)    :: dims, maxdims
    integer                           :: status
    character(len=*), parameter :: name = 'mls_hyperslab_save'
                                 ! Starting coordinatess of hyperslab
    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Begin execution
    call trace_begin ( me, 'MLS_hyperslab_save', cond=.false. )
    ! Check that pattern 1 or pattern 2 is satisfied
    status = 0
    if ( present(start) ) status = status + 1
    if ( present(count) ) status = status + 2
    if ( present(stride) ) status = status + 4
    if ( present(block) ) status = status + 8
    if ( status /= 3 .and. status /= 15 ) call my_message ( MLSMSG_Error, ModuleName, &
      & 'Impossible optional parameters pattern for dataset ' // trim(name), &
      & 'status', (/status/) )
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( DEEBUG ) print *, 'dataspace rank: ', rank
    if ( status /= 0 )  call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get rank for dataset ' // trim(name) )
    call h5sget_simple_extent_dims_f ( spaceID, dims(1:rank), maxdims(1:rank), &
      &  status )
    if ( DEEBUG ) print *, 'dataspace dims: ', dims(1:rank)
    if ( status /= rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dimension information for dataset ' // trim(name) )
    if ( present(stride) ) then
     if ( DEEBUG ) print *, 'trying to select hyperslab: ', &
       & int(start(1:rank), hsize_t), int(count(1:rank), hsize_t), &
       & int(stride(1:rank), hsize_t), int(block(1:rank), hsize_t)
      call h5sselect_hyperslab_f ( spaceID, H5S_SELECT_SET_F, &
        & int(start(1:rank), hsize_t), int(count(1:rank), hsize_t), status, &
        & int(stride(1:rank), hsize_t), int(block(1:rank), hsize_t) )
    else
      if ( DEEBUG ) print *, 'trying to select hyperslab: ', &
        & int(start(1:rank), hsize_t), int(count(1:rank), hsize_t)
      call h5sselect_hyperslab_f ( spaceID, H5S_SELECT_SET_F, &
        & int(start(1:rank), hsize_t), int(count(1:rank), hsize_t), status )
    end if
    if ( status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set hyperslab for dataset ' // trim(name) )
    call trace_end ( cond=.false. )
  end subroutine MLS_hyperslab_save

! ---------------------------------------------------  My_message  -----
  subroutine My_message ( severity, ModuleNameIn, Message, &
    & names, ints, reals, doubles, no_pairs )
    ! Take opportunity to dump a diagnostic table of values before stopping.
    ! Dummy arguments:
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character (len=*), intent(in) :: names   ! comma-separated list of names
    integer, dimension(:), optional, intent(in)          :: ints
    real, dimension(:), optional, intent(in)             :: reals
    double precision, dimension(:), optional, intent(in) :: doubles
    logical, optional, intent(in)                        :: no_pairs
    ! Local variables
    character(len=*), parameter :: int_format = '(i12)'
    character(len=*), parameter :: dbl_format = '(1pd12.2)'
    character(len=*), parameter :: real_format = '(1pe12.2)'
    integer, parameter          :: width = 2
    logical                     :: my_no_pairs
    my_no_pairs = .false.
    if ( present(no_pairs) ) my_no_pairs = no_pairs
    if ( my_no_pairs ) then
      if ( present(ints) ) then
        call dump ( ints, names, &
          & format=int_format, width=width )
      else if ( present(reals) ) then
        call dump ( reals, names )
      else if ( present(doubles) ) then
        call dump ( doubles, names )
      end if
    else
      if ( present(ints) ) then
        call dumpNamedValues ( ints, names, &
          & format=int_format, width=width )
      else if ( present(reals) ) then
        call dumpNamedValues ( reals, names, &
          & format=real_format, width=width )
      else if ( present(doubles) ) then
        call dumpNamedValues ( doubles, names, &
          & format=dbl_format, width=width )
      end if
    end if
    call MLSMessage ( severity, ModuleNameIn, message )
  end subroutine My_message

! -----------------------------------------------  Name_to_attrID  -----
  function Name_to_attrID ( fileID, dataName, attrName, dont_close ) &
   & result(attrID)
  ! Given a file, dataname, attribute name, return attribute ID of that name
  ! unless not there, when return -1
  ! Note that by default both attrName and dataname are left open
  ! (otherwise the returned number becomes invalid for closed names)
  ! (But how will we know what object number to close afterward?
  !  This is a bug.)
    integer, intent(in)               :: fileID
    character (len=*), intent(in)     :: dataNAME ! Name for this dataset
    character (len=*), intent(in)     :: attrNAME ! Name for this attribute
    logical, optional, intent(in)     :: dont_close ! default is TRUE
    integer                           :: attrID
    ! Local variables
    integer :: STATUS                   ! Flag
    integer                           :: dataID
    logical :: my_dont_close
    ! Executable code
    my_dont_close = .TRUE.
    if ( present(dont_close) ) my_dont_close = dont_close
    dataID = -1
    call h5dOpen_f ( fileID, dataName, dataID, status )
    if ( status /= 0 ) then
      attrID = -1
    else
      call h5aOpen_name_f ( dataID, attrName, attrID, status )
      if ( status /= 0 ) then
        attrID = -1
      else
        if ( .not. my_dont_close ) call h5aClose_f ( attrID, status )
      end if
      if ( .not. my_dont_close ) call h5dClose_f ( dataID, status )
    end if
  end function Name_to_attrID

! -----------------------------------------------  Name_to_dataID  -----
  function Name_to_dataID ( fileID, name, dont_close ) result(dataID)
  ! Given a file, dataname, return dataID of that name
  ! unless not there, when return -1
  ! Note that by default the dataname is left open
  ! (otherwise the returned number becomes invalid for a closed dataname)
    integer, intent(in)               :: fileID
    character (len=*), intent(in)     :: NAME ! Name for this dataset
    integer                           :: dataID
    logical, optional, intent(in)     :: dont_close ! default is TRUE
    ! Local variables
    integer :: STATUS                   ! Flag
    logical :: my_dont_close
    ! Executable code
    my_dont_close = .TRUE.
    if ( present(dont_close) ) my_dont_close = dont_close
    call h5dOpen_f ( fileID, name, dataID, status )
    if ( status /= 0 ) then
      dataID = -1
    else
      if ( .not. my_dont_close ) call h5dClose_f ( dataID, status )
    end if
  end function name_to_dataID

! ----------------------------------------------  StartMakeAttrib  -----
  logical function StartMakeAttrib ( ItemID, Name, Type, SpaceID, AttrID, &
    & Skip_if_already_there, Shp )
    ! Start up an attribute.  Return false if present(skip_if_already_there)
    ! and skip_if_already_there and the attribute is already there.
    integer, intent(in) :: ItemID         ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer(hid_t), intent(in) :: Type    ! Data type of attrib
    integer, intent(out) :: SpaceID       ! Dataspace ID
    integer, intent(out) :: AttrId        ! Attrib ID
    logical, intent(in), optional :: Skip_if_already_there ! Duh...
    integer(kind=hsize_t), intent(in), optional :: Shp(:) ! Shape of array attrib, else scalar

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: my_skip
    integer :: STATUS                   ! Flag from HDF5
    ! Executable
    call trace_begin ( me, 'StartMakeAttrib', cond=.false. )
    startMakeAttrib = .false.
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip = skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_dsID(itemID, name) ) go to 9
    end if
    startMakeAttrib = .true.
    if ( present(shp) ) then
      call h5sCreate_simple_f ( size(shp,1), int(shp,hSize_T), spaceID, status )
    else
      call h5sCreate_F ( h5s_scalar_f, spaceID, status )
    end if
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for attribute ' // trim(name) )
    ! Now create the attribute
    call h5aCreate_f ( itemID, trim(name), type, spaceID, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create attribute ' // trim(name) )
  9 call trace_end ( cond=.false. )
  end function StartMakeAttrib

! --------------------------------------------  WhatTypeAmI  -----
  function WhatTypeAmI ( me ) result ( Qtype )
    ! This routine returns the type of a type_id
    ! as an easily-recognized string
    ! E.g., 'character', 'real', 'double', or 'integer'
    integer, intent(in) :: me
    character(len=16)   :: Qtype
    ! Initialize in case something goes wrong with hdf5 routines
    Qtype = 'unknown'
    if ( AreThe2TypesEqual(me, H5T_NATIVE_INTEGER) .or. &
      &  AreThe2TypesEqual(me, H5T_STD_I32LE) ) then
      Qtype = 'integer'
    else if ( AreThe2TypesEqual(me, H5T_NATIVE_CHARACTER) ) then
      Qtype = 'character'
    else if ( AreThe2TypesEqual(me, H5T_STRING) ) then
      Qtype = 'character'
    else if ( AreThe2TypesEqual(me, H5T_NATIVE_REAL) .or. &
      &      AreThe2TypesEqual(me, H5T_IEEE_F32LE) ) then
      Qtype = 'real'
    else if ( AreThe2TypesEqual(me, H5T_NATIVE_DOUBLE) .or. &
      &      AreThe2TypesEqual(me, H5T_IEEE_F64LE) ) then
      Qtype = 'double'
    end if
  end function WhatTypeAmI

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSHDF5

! $Log$
! Revision 2.151  2024/06/03 20:39:45  pwagner
! Dont reference optional arg error
!
! Revision 2.150  2021/04/01 23:50:02  pwagner
! Strong words now warn of mis-named 'start' array
!
! Revision 2.149  2021/02/05 05:13:15  pwagner
! Prints more if rank too large
!
! Revision 2.148  2019/06/26 22:56:10  pwagner
! Needed H5SGet_Simple_Extent_NPoints_F to fix sometime link error
!
! Revision 2.147  2019/06/06 23:50:31  pwagner
! Can now read string-valued attributes with variable lengths
!
! Revision 2.146  2018/11/12 23:12:27  pwagner
! Added arg to convert carriagereturns to linefeeds
!
! Revision 2.145  2018/05/15 03:16:59  vsnyder
! Add LoadAllocFromHDF5DS
!
! Revision 2.144  2018/05/11 21:29:52  pwagner
! May be adding_to during SaveAsHDF5DS_textFile; if hdfVerbose print all warnings and quibbles
!
! Revision 2.143  2018/04/13 00:15:41  pwagner
! Added MLSFile_T api for GetHDF5DSRank; improved comments
!
! Revision 2.142  2017/08/10 22:42:23  pwagner
! Make apis for d.p. SaveAsHDF5DS like s.p.
!
! Revision 2.141  2017/08/04 19:42:14  pwagner
! Fixed case when groupName starts with /
!
! Revision 2.140  2017/07/25 22:30:34  pwagner
! Added MakeNestedGroups
!
! Revision 2.139  2017/02/09 23:45:48  pwagner
! Made more uses CamelCase
!
! Revision 2.138  2016/08/09 18:12:16  pwagner
! Made IsHDF5DSPresent generic
!
! Revision 2.137  2016/07/28 19:24:37  pwagner
! Fixed error in GetAllHDF5GroupNames
!
! Revision 2.136  2016/07/28 01:42:27  vsnyder
! Refactoring dump and diff
!
! Revision 2.135  2016/07/27 22:14:06  pwagner
! Added GetAllHDF5GroupNames
!
! Revision 2.134  2016/04/05 23:53:51  pwagner
! If -v option present, DumpHDF5DS will print DSName on each line
!
! Revision 2.133  2016/01/20 00:23:09  pwagner
! Added wildcard matches to DumpHDF5DS
!
! Revision 2.132  2015/07/14 23:17:17  pwagner
! May specify how to replace nulls in text file saved as DS
!
! Revision 2.131  2014/12/09 01:25:30  pwagner
! Comment-out a debugging print
!
! Revision 2.130  2014/10/27 23:04:06  pwagner
! needed to extend existing hdf attributes
!
! Revision 2.129  2014/09/05 00:05:14  vsnyder
! Convert some local pointer temps to allocatable
!
! Revision 2.128  2014/07/18 21:59:17  pwagner
! Pass item name to allocate_test
!
! Revision 2.127  2014/04/29 17:10:12  pwagner
! UFixed bugs regarding my_dont_trim in MakeHDF5Attribute_string
!
! Revision 2.126  2014/03/07 19:13:47  pwagner
! Increased MAXNDSNAMES; should there even be a limit?
!
! Revision 2.125  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.124  2013/11/01 00:05:13  pwagner
! Match trace_begin, _end in MakeHDF5Attribute_textFile
!
! Revision 2.123  2013/09/14 00:37:20  pwagner
! Corrected erroneous double call to trace_begin
!
! Revision 2.122  2013/08/31 01:24:53  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end
!
! Revision 2.121  2013/08/20 00:30:38  pwagner
! May also skip if attribute not where we cp from
!
! Revision 2.120  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.119  2012/05/17 20:19:07  pwagner
! Fixed bug preventing cleanly dumping MAFStartTimeUTC
!
! Revision 2.118  2012/05/08 01:15:42  vsnyder
! Option allows numeric DS silently to be absent
!
! Revision 2.117  2012/01/25 01:12:07  pwagner
! Fixed most of the bugs in MatchHDF5Attributes
!
! Revision 2.116  2012/01/13 01:10:14  pwagner
! Added MatchHDF5Attributes; untested yet
!
! Revision 2.115  2011/11/04 23:41:16  pwagner
! Removed unused items
!
! Revision 2.114  2011/11/03 23:48:10  pwagner
! Added MLSFile api for GetHDF5DSDims
!
! Revision 2.113  2011/11/01 21:01:34  pwagner
! GetAllHDF5DSNames now accepts MLSFile as arg
!
! Revision 2.112  2011/10/25 17:56:09  pwagner
! May make scalar real and d.p.-valued attributes
!
! Revision 2.111  2011/08/04 17:38:18  pwagner
! Corrected choice of kind for stringsize
!
! Revision 2.110  2011/08/04 16:57:35  honghanh
! Revert the previous bug fix
! because the issue is with toolkit
!
! Revision 2.109  2011/08/04 16:45:06  honghanh
! Fix a bug in the new read attribute function
!
! Revision 2.108  2011/08/02 16:53:00  honghanh
! Add ReadHDF5Attr_FID_string and ReadHDF5Attr_FID_int functions
!
! Revision 2.107  2011/07/15 23:33:11  pwagner
! Can now dump rank 4 datasets
!
! Revision 2.106  2011/05/06 00:35:33  pwagner
! Fixed bug in IsHDF5AttributePresent
!
! Revision 2.105  2011/05/05 15:14:45  pwagner
! Added MLSFile api to IsHDF5AttributePresent
!
! Revision 2.104  2011/03/16 18:19:47  pwagner
! Workaround for bug in Intel v12.0.0-2011.1.107
!
! Revision 2.103  2011/02/05 01:35:34  pwagner
! Consistent with new dopt_ dump options
!
! Revision 2.102  2011/01/12 18:08:33  pwagner
! Wont truncate dumped long string attributes
!
! Revision 2.101  2010/11/23 01:11:07  pwagner
! Improved dumps of character arrays with many elements
!
! Revision 2.100  2010/08/27 20:59:06  pwagner
! Reduced character length array name chLongArray for long arrays so we dont exhaust memory
!
! Revision 2.99  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.98  2010/01/11 18:33:57  pwagner
! Fixed bug in dumping array-valued attributes
!
! Revision 2.97  2009/11/16 21:55:55  pwagner
! orked around lack of generic outputnamedValue for longint
!
! Revision 2.96  2009/11/10 00:31:43  pwagner
! Raised character string size in CpHDF5Attribute_string consistent with PCFHdr%MiscNotesLENGTH
!
! Revision 2.95  2009/10/05 23:38:59  pwagner
! Moved use h5lib statements from module scope to speedup Lahey; this is the last time we do that
!
! Revision 2.94  2009/09/29 23:31:24  pwagner
! Changes needed by 64-bit build
!
! Revision 2.93  2009/08/04 20:44:08  pwagner
! Now able to dump character scalar ds
!
! Revision 2.92  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.91  2009/06/16 17:14:24  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.90  2009/05/08 00:42:04  pwagner
! New optional arg laconic to prevent DumpHDF5DS from printing names
!
! Revision 2.89  2009/04/01 23:27:27  pwagner
! Worked around bug in saving logical array
!
! Revision 2.88  2008/06/18 20:56:55  pwagner
! New optional arg 'unique' dumps print unique elements, counts only
!
! Revision 2.87  2008/05/24 00:54:46  vsnyder
! Remove unused declarations
!
! Revision 2.86  2008/05/24 00:34:05  pwagner
! Removed scraps from older reading of text files
!
! Revision 2.85  2008/05/22 01:06:42  vsnyder
! Use new allocate_test to simplify LoadPtrFromDS_*
!
! Revision 2.84  2008/05/20 01:59:28  vsnyder
! Add 4d integer, real, double
!
! Revision 2.83  2008/05/10 01:14:04  vsnyder
! Use includes to generate multiple-rank routines
!
! Revision 2.82  2008/05/02 00:05:19  pwagner
! Reads textFile using io stuff
!
! Revision 2.81  2008/04/25 22:51:23  pwagner
! Passes optional arg maxLineLen through to textFile_to_Chars
!
! Revision 2.80  2008/04/22 17:17:38  pwagner
! Should not print extra debugging stuff
!
! Revision 2.79  2008/04/18 16:29:06  pwagner
! Now works properly with NAG, Lahey, and Intel
!
! Revision 2.78  2008/03/07 01:35:47  pwagner
! Will subdivide attributes into blocks if too large
!
! Revision 2.77  2008/02/22 21:29:13  pwagner
! Can now save entire textfile as attribute or dataset
!
! Revision 2.76  2008/01/04 01:00:13  pwagner
! More successful at dumping l1b files
!
! Revision 2.75  2007/10/04 20:43:36  vsnyder
! Remove unused symbols
!
! Revision 2.74  2007/10/03 23:51:53  vsnyder
! Don't overflow MLSCallStack
!
! Revision 2.73  2007/08/17 00:27:48  pwagner
! push more procedures onto MLSCallStack
!
! Revision 2.72  2007/01/31 00:07:23  pwagner
! Compatible with renamed dumpNamedValues
!
! Revision 2.71  2007/01/26 23:59:13  pwagner
! Added more unused debugging
!
! Revision 2.70  2007/01/12 00:29:28  pwagner
! Renamed routine outputNamedValue
!
! Revision 2.69  2006/08/23 18:04:23  pwagner
! NAG hates splitting /) operator
!
! Revision 2.68  2006/08/22 20:41:23  pwagner
! Fixed a bug in DumpHDF5Attributes
!
! Revision 2.67  2006/07/11 00:24:36  pwagner
! use fillValue properly when computing rms etc.
!
! Revision 2.66  2006/06/29 20:38:20  pwagner
! Added a few extra error checks
!
! Revision 2.65  2006/06/27 23:59:05  pwagner
! May dump attributes, datasets
!
! Revision 2.64  2006/04/12 20:51:24  pwagner
! Attempts to work around hdf5-1.6.5 bugs rewriting string attributes
!
! Revision 2.63  2006/01/25 00:57:39  pwagner
! Removed some troublesome, superfluous calls from IsHDF5ItemPresent
!
! Revision 2.62  2005/11/17 20:09:24  pwagner
! LoadFromHDF5DS can now read 3d integer arrays
!
! Revision 2.61  2005/11/11 21:36:00  pwagner
! Added new interface for IsHDF5AttributeInFile when attached to grp, not dataset
!
! Revision 2.60  2005/10/28 23:13:42  pwagner
! Prevent ref to undefined num in GetAllHDF5AttrNames
!
! Revision 2.59  2005/10/19 20:46:05  pwagner
! GetAllHDF5DSNames will omit leading slashes by default
!
! Revision 2.58  2005/10/18 23:01:43  pwagner
! Added GetAllHDF5AttrNames, IsHDF5ItemPresent; introduced options
!
! Revision 2.57  2005/10/11 17:35:00  pwagner
! MLSFile interface to GetHDF5Attribute can get attribute from grp_id component if sd_id is 0
!
! Revision 2.56  2005/07/12 17:12:50  pwagner
! New hdf5 library will drop integer dimension interfaces
!
! Revision 2.55  2005/06/29 00:40:48  pwagner
! New interfaces for GetHDF5Attribute and LoadFromHDF5DS accept MLSFiles
!
! Revision 2.54  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.53  2005/04/29 21:55:17  pwagner
! Nullified name component of data_set info before allocating
!
! Revision 2.52  2005/01/12 03:04:15  vsnyder
! Use correct data type in MakeHDF5Attribute_sngl, some cannonball polishing
!
! Revision 2.51  2005/01/07 01:57:53  vsnyder
! Don't try to get DS dims from an Attrib
!
! Revision 2.50  2005/01/07 01:14:26  vsnyder
! Allocate the correct upper bound in LoadPtrFromHDF5DS_logarr1
!
! Revision 2.49  2005/01/07 01:04:37  vsnyder
! Use generics for some internal references
!
! Revision 2.48  2005/01/07 00:38:17  vsnyder
! Add IsHDF5GroupPresent, simplify some stuff, delete unused stuff
!
! Revision 2.47  2004/12/31 02:38:53  vsnyder
! Added LoadPtrFromHDF5DS, simplified a lot of stuff
!
! Revision 2.46  2004/12/13 20:26:53  vsnyder
! Added MakeHDF5Attribute_sngl.  Change specifics for generics to use Real
! and Double Precision instead of r4 (which isn't guaranteed to be default
! real) and r8 (which isn't guaranteed to be different from default real).
! Changed a few references to specifics to refer to generics.
!
! Revision 2.45  2004/09/23 23:00:10  pwagner
! Added CpHDF5GlAttribute, CpHDF5Attribute
!
! Revision 2.44  2004/08/04 23:19:01  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.43  2004/08/03 18:00:28  pwagner
! Optionally sets fill value for some datasets
!
! Revision 2.42  2004/07/22 20:37:11  cvuu
! Change in MakeHDF5Attribute_string to allow re-write the char string in the file attribute
!
! Revision 2.41  2004/06/29 00:07:00  pwagner
! Exploit catlist function
!
! Revision 2.40  2004/06/08 18:59:59  pwagner
! Another break with old toolkit; this time may be correct
!
! Revision 2.39  2004/05/25 20:35:53  pwagner
! Reverted to older hdf5; temporarily we hope
!
! Revision 2.38  2004/05/19 19:08:05  pwagner
! After hdf5-1.6.2, h5open_f, h5close_f become module procedures
!
! Revision 2.37  2004/03/25 18:38:13  pwagner
! May save 3d integer and double dsets
!
! Revision 2.36  2004/03/24 23:50:27  pwagner
! Added mls_h5open/close
!
! Revision 2.35  2004/02/26 21:59:09  pwagner
! Added GetAllHDF5DSNames
!
! Revision 2.34  2003/10/02 23:09:47  pwagner
! Some small fixes; can get string array attribute
!
! Revision 2.33  2003/09/30 18:29:41  perun
! Change len_trim to len in MakeHDF5Attribute_string_arr1.
!
! Revision 2.32  2003/09/12 16:40:27  cvuu
! Add subroutines to get L1BOA attributes
!
! Revision 2.31  2003/08/07 15:44:19  perun
! Add MakeHDF5Attribute_intarr1
!
! Revision 2.30  2003/07/24 22:10:45  pwagner
! Fixed another bug preventing multiple chunks from writing to same dataset
!
! Revision 2.29  2003/07/21 23:30:23  pwagner
! Check on returnstatus from h5sget_simple_extent_dims_f being rank (marking success)
!
! Revision 2.28  2003/07/18 16:04:19  pwagner
! Fixed some bugs in DirectWriting 3-d datasets
!
! Revision 2.27  2003/07/15 23:37:57  pwagner
! No changes I can see, but cvs says so, so ..
!
! Revision 2.26  2003/05/19 22:06:31  pwagner
! Shortened names to Read..IndexFromHDF5Attr to comply with namelength standard
!
! Revision 2.25  2003/05/13 04:46:42  livesey
! Bug fix, added more specifics to generic
!
! Revision 2.24  2003/05/12 18:08:10  pwagner
! Added 1d, 2d char array dsets, 2d int dset, gets of sngl, dbl scalar attrs
!
! Revision 2.23  2003/05/12 16:51:20  livesey
! Added the lit and string index stuff
!
! Revision 2.22  2003/04/28 23:08:20  pwagner
! Added ability to Load, Save rank4 s.p. arrays
!
! Revision 2.21  2003/03/20 19:21:11  pwagner
! Fixed simple bug in saving 3d arrays
!
! Revision 2.20  2003/02/21 23:41:02  pwagner
! Additional MakeHDF5Attribute interface
!
! Revision 2.19  2003/02/12 21:38:14  pwagner
! May make dbl scalar and array attributes; find if name is an attribute of a group
!
! Revision 2.18  2003/01/30 00:56:01  pwagner
! Added string arrays as possible attributes (untested)
!
! Revision 2.17  2003/01/27 21:38:44  pwagner
! May make 1d s.p. array attributes; may make attributes with sdname attrname call
!
! Revision 2.16  2003/01/23 23:30:49  pwagner
! May add to same 2d, 3d single-precision datasets
!
! Revision 2.15  2002/12/07 00:24:40  pwagner
! Added SaveAsHDF5DS_snglarr3
!
! Revision 2.14  2002/12/02 23:35:57  pwagner
! Should provide more info when something goes awry
!
! Revision 2.13  2002/10/29 01:01:05  pwagner
! Can save a char scalar as a DS
!
! Revision 2.12  2002/10/11 23:42:04  pwagner
! Remembered to close hyperslab memspaceID if created one
!
! Revision 2.11  2002/10/10 23:51:57  pwagner
! Optional hyperslab args to LoadFromHDF5DS
!
! Revision 2.10  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.9  2002/10/04 22:22:52  pwagner
! Fixed bug in GetHDF5DSQType; can retrieve rank3 datasets
!
! Revision 2.8  2002/10/02 23:20:07  livesey
! Bug fix in single precision stuff
!
! Revision 2.7  2002/09/27 23:39:26  pwagner
! Added GetHDF5DSQType
!
! Revision 2.6  2002/09/26 23:56:15  pwagner
! Added some things for MLSAux and l1bdata
!
! Revision 2.5  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.4  2002/08/26 16:42:09  livesey
! Bug fix with error messages in IsHDF5DSPresent
!
! Revision 2.3  2002/08/23 01:23:21  livesey
! Added IsHDF5DSPresent and IsHDF5AttributePresent
!
! Revision 2.2  2002/07/17 06:00:21  livesey
! Got hdf5 l2pc reading stuff working
!
! Revision 2.1  2002/07/11 22:18:26  pwagner
! First commit in this directory--welcome old friendshe5*.f90
!
! Revision 1.1  2002/06/18 21:57:09  livesey
! First version
!
