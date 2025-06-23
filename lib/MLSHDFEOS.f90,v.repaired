! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MLSHDFEOS

  ! This module contains MLS specific routines to do common HDFEOS
  ! tasks not specifically found in he5_swapi or he5_gdapi.
  
  ! We are in mid-transition between the older FileID/FileName interfaces
  ! and the newer MLSFile defined-type interfaces

  use HDF, only: Dfnt_Char8, Dfnt_Float32, Dfnt_Float64, &
    & Dfnt_Int8, Dfnt_Int16, Dfnt_Int32, Dfnt_Int64
  use HDFEOS, only: Gdattach, Gdcreate, &
    & Swattach, Swcreate, Swdefdfld, Swdefgfld, Swdefdim, Swdetach, &
    & Swdiminfo, Swinqdflds, Swinqswath
  use HDFEOS5, only: HE5t_Native_Float, HE5t_Native_Double, &
    & HE5t_Native_Int, HE5t_Native_Int8, HE5t_Native_Int16, HE5t_Native_Int64
  use HDFEOS5, only: HE5_Gdattach, HE5_Gdcreate, &
    & HE5_Swattach, HE5_Swcreate, HE5_Swdefchunk, HE5_Swdetach, HE5_Swdefdfld, &
    & HE5_Swdefgfld, HE5_Swdefdim, HE5_Swdiminfo, HE5_Swinqdflds, &
    & HE5_Swinqswath, MLS_Chartype
  use HE5_Swapi, only: HE5_Swsetfill, HE5_Swwrattr, HE5_Swwrlattr
  use HE5_Swapi_Character_Array, only: HE5_Ehwrglatt_Character_Array, &
    & HE5_Ehrdglatt_Character_Array
  use HE5_Swapi_Character_Scalar, only: HE5_Ehwrglatt_Character_Scalar, &
    & HE5_Ehrdglatt_Character_Scalar
  use HE5_Swapi_Double, only: HE5_Ehwrglatt_Double, HE5_Ehrdglatt_Double, &
    & HE5_Swrdfld_Double, HE5_Swrdfld_Double_2d, HE5_Swrdfld_Double_3d, &
    & HE5_Swwrfld_Double, HE5_Swwrfld_Double_2d, HE5_Swwrfld_Double_3d
  use HE5_Swapi_Integer, only: HE5_Ehwrglatt_Integer, &
    & HE5_Swrdfld_Integer, HE5_Swwrfld_Integer, HE5_Ehrdglatt
  use HE5_Swapi_Real, only: HE5_Ehwrglatt_Real, HE5_Ehrdglatt_Real, &
    & HE5_Swrdfld_Real, HE5_Swrdfld_Real_2d, HE5_Swrdfld_Real_3d, &
    & HE5_Swwrfld_Real, HE5_Swwrfld_Real_2d, HE5_Swwrfld_Real_3d
  use HighOutput, only: OutputnamedValue
  use MLSCommon, only: MLSFile_T
  use MLSFiles, only: HDFversion_4, HDFversion_5, WildcardHDFversion, &
    & MLS_HDF_Version
  use MLSMessagemodule, only: MLSMSG_Error, MLSMSG_Warning, MLSMessage
  use MLSStringlists, only: StringElementnum
  use MLSStrings, only: Replace
  use Swapi_Double, only: Swrdfld_Double, Swrdfld_Double_2d, Swrdfld_Double_3d, &
    & Swwrfld_Double, Swwrfld_Double_2d, Swwrfld_Double_3d
  use Swapi_Integer, only: Swrdfld_Integer, &
    & Swwrfld_Integer
  use Swapi_Real, only: Swrdfld_Real, Swrdfld_Real_2d, Swrdfld_Real_3d, &
    & Swwrfld_Real, Swwrfld_Real_2d, Swwrfld_Real_3d
  use Trace_M, only: Trace_Begin, Trace_End

  implicit none
  private

  public :: He5_Ehwrglatt, He5_Ehrdglatt, Hsize, Hsizes, &
    & MLS_Ehwrglatt, MLS_Isglatt, &
    & MLS_Dfldsetup, MLS_Gfldsetup, &
    & MLS_Swdefdim, MLS_Swdiminfo, MLS_Swrdfld, MLS_Swwrfld, &
    & MLS_Swattach, MLS_Swcreate, MLS_Swdetach, MLS_Gdcreate, MLS_Gdwrattr, &
    & MLS_Swwrattr, MLS_Swwrlattr, &
    & MLS_Swath_In_File

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! HE5_EHRDGLATT     Reads the global attributes at file level
! HE5_EHWRGLATT     Writes the global attributes at file level
! HSIZE             Return arg with same integer value
!                     but with kind value size_t
! HSIZES            Return array as HSIZE
! MLS_DFLDSETUP     Sets up a data field in a swath
! MLS_GDCREATE      Creates a grid (perhaps prior to writing)
! MLS_GDWRATTR      Writes a grid-level attribute (e.g., pressures)
! MLS_GFLDSETUP     Sets up a geolocation field in a swath
! MLS_ISGLATT       Is the named global attribute present in the file?
! MLS_SWATH_IN_FILE Is named swath in file?
! MLS_SWATTACH      Attaches a swath (perhaps prior to reading)
! MLS_SWCREATE      Creates a swath (perhaps prior to writing)
! MLS_SWDEFDIM      Sets up a dimension for a swath
! MLS_SWDETACH      Detaches from a swath (perhaps prior to closing the file)
! MLS_SWDIMINFO     Gets info on a dimension for a swath
! MLS_SWRDFLD       Reads a field from a swath, data or geolocation
! MLS_SWWRFLD       Writes a field to a swath, data or geolocation
! MLS_SWWRATTR      Writes a swath-level attribute (e.g., pressures)
! MLS_SWWRLATTR     Writes a datafield-level attribute (e.g., units)
! === (end of toc) ===

! === (start of api) ===
! int HE5_EHRDGLATT (int fileID, char* attrName, value)
!     value can be one of:
!    {char* value, char* value(:), int value(:), r4 value(:), r8 value(:)}
! int HE5_EHWRGLATT (int fileID, char* attrName, int datatype, int count, value)
!     value can be one of:
!    {char* value, char* value(:), int value(:), r4 value(:), r8 value(:)}
! size_t HSIZE (int arg)
! size_t(:) HSIZES (int ints(:))
! int MLS_SWSETFILL (int swathID, char* names, int datatype, value) 
!     value can be one of:
!    {char* value, int value, r4 value, r8 value}
! int MLS_DFLDSETUP (int swathID, char* fieldname, char* dimname, 
!     char* MAXDIMList, int datatype, int merge, int chunk_rank, 
!     int chunk_dims(:), [char* filename], [int hdfVersion], [log dontFail],
!     [int iFill], [r4 rFill], [r8 dFill]) 
! int MLS_GDCREATE (int fileID, char* gridname, 
!     int xdimsize, int ydimsize, r8 upleft(2), r8 lowright(2), 
!     [char* filename], [int hdfVersion]) 
! int MLS_GDWRATTR (int gridID, char* attrName, int datatype, int count, 
!     char* buffer)
! int MLS_GFLDSETUP (int swathID, char* fieldname, char* dimname, 
!     char* MAXDIMList, int datatype, int merge, int chunk_rank, 
!     int chunk_dims(:), [char* filename], [int hdfVersion], [log dontFail],
!     [int iFill], [r4 rFill], [r8 dFill]) 
! log MLS_ISGLATT (file, char* attrname)
!     file can be one of: int fileID, or char* filename
! log MLS_SWATH_IN_FILE (file, char* attrname)
! int MLS_SWATTACH (file, char* swathname, [char* filename], [int hdfVersion], 
!     [log dontFail]) 
!     file can be one of: int fileID, or type(MLSFile_T) MLSFile
! int MLS_SWCREATE (file, char* swathname, [char* filename], [int hdfVersion]) 
!     file can be one of: int fileID, or type(MLSFile_T) MLSFile
! int MLS_SWdefdim (int swathID, char* dimname, 
!     int dimsize, [char* filename], [int hdfVersion], [log dontFail]) 
! int MLS_SWDETACH (int swathID, [char* filename], [int hdfVersion])
! int MLS_SWdiminfo (int swathID, char* dimname, 
!     [char* filename], [int hdfVersion], [log dontFail]) 
! int MLS_SWRDFLD (int swathID, char* fieldname, 
!     int start(rank), int stride(rank), int edge(rank), values,
!     [char* filename], [int hdfVersion], [log dontFail]) 
!     values can be one of:
!    {char* value, char* value(:), int value(:), r4 value(shp), r8 value(shp)}
!     and shp is an array of ints of size rank
! int MLS_SWWRFLD (int swathID, char* fieldname, 
!     int start(rank), int stride(rank), int edge(rank), values,
!     [char* filename], [int hdfVersion], [log dontFail]) 
!     values can be one of types in MLS_SWRDFLD
! int MLS_SWWRATTR (int swathID, char* attrName, int datatype, int count, 
!     char* buffer)
! int MLS_SWWRLATTR (int swathID, char* fieldName, char* attrName, 
!     int datatype, int count, char* buffer)
! === (end of api) ===
  interface HE5_EHWRGLATT   ! From its name, might better be in he5_ehapi.f90
    module procedure HE5_EHWRGLATT_CHARACTER_SCALAR, HE5_EHWRGLATT_DOUBLE, &
    HE5_EHWRGLATT_INTEGER, HE5_EHWRGLATT_REAL, HE5_EHWRGLATT_CHARACTER_ARRAY
  end interface

  interface HE5_EHRDGLATT   ! From its name, might better be in he5_ehapi.f90
    module procedure HE5_EHRDGLATT_CHARACTER_SCALAR, &
      & HE5_EHRDGLATT_DOUBLE, &
      & HE5_EHRDGLATT_REAL, HE5_EHRDGLATT_CHARACTER_ARRAY
  end interface

  interface MLS_EHWRGLATT
    module procedure MLS_EHWRGLATT_char, MLS_EHWRGLATT_textfile
  end interface

  interface MLS_ISGLATT
    module procedure MLS_ISGLATT_FID, &
      & MLS_ISGLATT_FN
  end interface

  interface MLS_SWATH_IN_FILE
    module procedure MLS_SWATH_IN_FILE_ARR, &
      & MLS_SWATH_IN_FILE_SCA
  end interface

  interface MLS_SWATTACH
    module procedure MLS_SWATTACH_ID, &
      & MLS_SWATTACH_MF
  end interface

  interface MLS_SWCREATE
    module procedure MLS_SWCREATE_ID, &
      & MLS_SWCREATE_MF
  end interface

  interface MLS_SWRDFLD
    module procedure &
      & MLS_SWRDFLD_CHAR_1D, &
      & MLS_SWRDFLD_DOUBLE_1D, MLS_SWRDFLD_DOUBLE_2D, MLS_SWRDFLD_DOUBLE_3D, &
      & MLS_SWRDFLD_INTEGER, &
      & MLS_SWRDFLD_REAL_1D, MLS_SWRDFLD_REAL_2D, MLS_SWRDFLD_REAL_3D
  end interface

  interface MLS_SWWRFLD
    module procedure &
      & MLS_SWWRFLD_CHAR_1D, &
      & MLS_SWWRFLD_DOUBLE_1D, MLS_SWWRFLD_DOUBLE_2D, MLS_SWWRFLD_DOUBLE_3D, &
      & MLS_SWWRFLD_INTEGER, &
      & MLS_SWWRFLD_REAL_1D, MLS_SWWRFLD_REAL_2D, MLS_SWWRFLD_REAL_3D
  end interface

  integer, public, parameter :: MAXNODFIELDS = 1000
  integer, public, parameter :: MAXDLISTLENGTH = 400000
  integer, parameter         :: DFLTMAXLINELENGTH = 1024

  ! Print debugging stuff?
  logical, parameter :: DEEBUG = .false.  
  character(len=1), parameter :: BLANK = ' '

contains ! ======================= Public Procedures =========================

  ! ----------------------------------------------  MLS_EHWRGLATT  -----
  integer function MLS_EHWRGLATT_char ( FILEID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_CHARTYPE
    integer, intent(in) :: COUNT   ! How many
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'MLS_EHWRGLATT_char' , cond=.false. )
    if ( len_trim(buffer) > 0 ) then
      MLS_EHWRGLATT_char = he5_ehwrglatt_character_scalar( FILEID, &
      & ATTRNAME, DATATYPE, hsize(max(COUNT, len_trim(BUFFER))), BUFFER )
    else
      MLS_EHWRGLATT_char = he5_ehwrglatt_character_scalar( FILEID, &
      & ATTRNAME, DATATYPE, hsize(1), BLANK )
    endif
    call trace_end ( cond=.false. )

  end function MLS_EHWRGLATT_char

  integer function MLS_EHWRGLATT_textfile ( textFile, FILEID, &
    & ATTRNAME, maxLineLen )
    ! Write contents of text file as global attribute
    ! E.g., may wish to include leap sec file or utc pole
    use IO_STUFF, only: read_textFile
    use MLSHDF5, only: MAXCHATTRLENGTH, MAXCHFIELDLENGTH

    character (len=*), intent(in) :: TEXTFILE ! name of textfile
    integer, intent(in) :: FILEID      ! File ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, optional, intent(in) :: maxLineLen

    ! Internal variables
    integer, parameter :: DATATYPE = MLS_CHARTYPE
    character(len=3) :: blockChar
    character(len=MAXCHFIELDLENGTH) :: BUFFER  ! Buffer to hold contents
    integer :: firstChar, lastChar
    integer :: iblock, nblocks
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: myMaxLineLen
    character(len=len(attrname)+3) :: newname
    integer :: status

    ! Executable code
    call trace_begin ( me, 'MLS_EHWRGLATT_textfile' , cond=.false. )
    myMaxLineLen = DFLTMAXLINELENGTH
    if ( present(maxLineLen) ) myMaxLineLen = maxLineLen
    ! Try to read the textfile
    BUFFER = ' '
    call read_textFile( trim(textFile), BUFFER, myMaxLineLen )
    status = 0
    if ( status /= 0 ) then
      call MLSMessage(MLSMSG_Warning, ModuleName, &
        & 'Unable to write attribute--failed to read textfile' )
      call trace_end ( cond=.false. )
      return
    endif
    ! Unfortunately, a lot of null characters sneak into this
    BUFFER = Replace( BUFFER, char(0), char(32) ) ! Replace null with space
    ! Be careful lest the attribute is too large
    if ( len_trim(buffer) > MAXCHATTRLENGTH ) then
      nblocks = 1 + ( len_trim(buffer) - 1 ) / MAXCHATTRLENGTH
      lastChar = 0
      do iblock=1, nblocks
        firstChar = lastChar + 1
        lastChar  = min( lastChar + MAXCHATTRLENGTH, len_trim(buffer) )
        write(blockChar, '(i1)' ) iblock
        if ( iblock > 9 ) write(blockChar, '(i2)' ) iblock
        newName = trim(attrname) // adjustl(blockChar)
        MLS_EHWRGLATT_textfile = &
          & MLS_EHWRGLATT_char( fileID, newName, &
          & MLS_CHARTYPE, lastChar-firstChar+1, buffer(firstChar:lastChar) )
      enddo
      call trace_end ( cond=.false. )
      return
    elseif ( len_trim(buffer) > 0 ) then
      MLS_EHWRGLATT_textfile = he5_ehwrglatt_character_scalar( FILEID, &
      & ATTRNAME, DATATYPE, hsize(max(1, len_trim(BUFFER))), BUFFER )
    else
      MLS_EHWRGLATT_textfile = he5_ehwrglatt_character_scalar( FILEID, &
      & ATTRNAME, DATATYPE, hsize(1), BLANK )
    endif
    call trace_end ( cond=.false. )

  end function MLS_EHWRGLATT_textfile

  ! ---------------------------------------------  MLS_ISGLATT_FN  -----
  function MLS_ISGLATT_FN ( FILENAME, ATTRNAME ) result(isThere)
    ! Is the named attribute a global attribute of the named file?
    use HDFEOS5, only: he5_swclose, he5_swopen, HE5F_ACC_RDONLY
    ! Args
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    logical                      :: isThere

    ! Local variables
    integer :: fileID
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: status

    ! Executable code
    call trace_begin ( me, 'MLS_ISGLATT_FN' , cond=.false. )
    isThere = .false.
    fileID = he5_swopen(trim(fileName), HE5F_ACC_RDONLY)
    if ( fileID > 0 ) then
      isThere = MLS_ISGLATT ( fileID, ATTRNAME )
      status = he5_swclose( fileID )
    endif
    call trace_end ( cond=.false. )
  end function MLS_ISGLATT_FN

  ! --------------------------------------------  MLS_ISGLATT_FID  -----
  function MLS_IsGlatt_fid ( Fileid, Attrname ) result(isThere)
    ! Is the named attribute a global attribute of the file?
    use HDFEOS5, only: HE5_EHinqglatts
    use MLSStringLists, only: StringElementNum
    ! Args
    integer, intent(in)           :: FileID
    character(len=*), intent(in)  :: Attrname     ! Attribute name
    logical                       :: IsThere
    ! Local variables
    character(len=MAXDLISTLENGTH) :: AttrList
    logical, parameter            :: Countempty = .true.
    integer                       :: ListSize
    integer                       :: Status
    logical, parameter            :: DEEBUG = .false.
    ! Executable code
    isThere = .false.
    attrList = ' ' ! Why is this necssary?
    listSize = 0
    status = he5_EHinqglatts(fileID, attrList, listSize)
    ! if ( status /= 0 ) return
    if ( len_trim(attrList) > 0 ) then
      if ( DEEBUG ) print *, 'glattr list: ', trim(attrList)
      status = StringElementNum( trim(attrList), trim(attrName), countempty)
      isThere = ( status > 0 )
    endif
  end function MLS_ISGLATT_FID

  ! -----------------------------------------------  MLS_GDCREATE  -----
  integer function MLS_GDCREATE ( FILEID, GRIDNAME, &
   &  xdimsize, ydimsize, upleft, lowright, FileName, hdfVersion )
    integer, intent(in) :: FILEID      ! ID returned by MLS_swopen
    character(len=*), intent(in) :: GRIDNAME       ! Swath name
    integer, intent(in) :: xdimsize
    integer, intent(in) :: ydimsize
    double precision, dimension(2), intent(in) :: upleft
    double precision, dimension(2), intent(in) :: lowright
    character(len=*), optional, intent(in) :: FILENAME  ! File name
    integer, optional, intent(in) :: hdfVersion

    ! Internal variables
    logical :: alreadyThere
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: myHdfVersion
    logical :: needsFileName
    logical, parameter :: MUSTCREATE = .true.

    ! Executable code
    call trace_begin ( me, 'MLS_GDCREATE' , cond=.false. )
    MLS_GDCREATE = 0
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_GDCREATE' )
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      if ( MUSTCREATE ) then
        alreadyThere = .false.
      else
        alreadyThere = (gdattach(FileID, trim(GRIDNAME)) >= 0)
      endif
      if ( alreadyThere ) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'GRIDNAME call to MLS_GDCREATE already there: ' // trim(GRIDNAME))
      endif
      MLS_GDCREATE = gdcreate(Fileid, trim(GRIDNAME), &
        & xdimsize, ydimsize, upleft, lowright)
    case (HDFVERSION_5)
      if ( MUSTCREATE ) then
        alreadyThere = .false.
      else
        alreadyThere = (he5_gdattach(FileID, trim(GRIDNAME)) >= 0)
      endif
      if ( alreadyThere ) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'GRIDNAME call to MLS_GDCREATE already there: ' // trim(GRIDNAME))
      endif
      MLS_GDCREATE = he5_gdcreate(Fileid, trim(GRIDNAME), &
        & hsize(xdimsize), hsize(ydimsize), upleft, lowright)
    case default
      MLS_GDCREATE = -1
    end select
    if ( MLS_GDCREATE == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to create grid name ' // trim(GRIDNAME) )

    call trace_end ( cond=.false. )
  end function MLS_GDCREATE

  ! -----------------------------------------------  MLS_gdwrattr  -----
  integer function MLS_gdwrattr ( GRIDID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: GRIDID      ! Grid ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_CHARTYPE
    integer, intent(in) :: COUNT   ! How many
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    integer, external ::   he5_GDwrattr
    call trace_begin ( me, 'MLS_gdwrattr' , cond=.false. )
    if ( len_trim(buffer) > 0 ) then
      MLS_gdwrattr = HE5_GDWRATTR( GRIDID, &
      & ATTRNAME, DATATYPE, hsize(max(COUNT, len_trim(BUFFER))), BUFFER )
    else
      MLS_gdwrattr = HE5_GDWRATTR( GRIDID, &
      & ATTRNAME, DATATYPE, hsize(1), BLANK )
    endif

    call trace_end ( cond=.false. )
  end function MLS_gdwrattr

  ! --------------------------------------------  MLS_SWATTACH_ID  -----
  function MLS_SWATTACH_ID ( FILEID, SWATHNAME, FileName, &
    &  hdfVersion, DONTFAIL ) result(MLS_SWATTACH)
    integer, intent(in) :: FILEID      ! ID returned by MLS_swopen
    character(len=*), intent(in) :: SWATHNAME       ! Swath name
    character(len=*), optional, intent(in) :: FILENAME  ! File name
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    integer :: MLS_SWATTACH

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_swattach_id' , cond=.false. )
    MLS_SWATTACH = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWATTACH = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWATTACH' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
      if ( myHdfVersion < 0 ) print *, 'uh-oh, MLS_hdf_version: ', myhdfVersion
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWATTACH = swattach(Fileid, trim(swathName))
    case (HDFVERSION_5)
      MLS_SWATTACH = he5_swattach(Fileid, trim(swathName))
    case default
      MLS_SWATTACH = -1
    end select
    if (Deebug) print *, ' (swath id is ', MLS_SWATTACH, ')'
    if ( MLS_SWATTACH /= -1 ) then
      call trace_end ( cond=.false. )
      return
    elseif ( myDontFail ) then
      CALL MLSMessage ( MLSMSG_Warning, moduleName,  &
          & 'Failed to attach swath name ' // trim(swathname) )
    else
      CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to attach swath name ' // trim(swathname) )
    endif
    call trace_end ( cond=.false. )
  end function MLS_SWATTACH_ID

  ! --------------------------------------------  MLS_SWATTACH_MF  -----
  function MLS_SWATTACH_MF ( MLSFile, SWATHNAME, DONTFAIL ) &
    & result(MLS_SWATTACH)
    type (MLSFile_T)   :: MLSFile
    character(len=*), intent(in) :: SWATHNAME       ! Swath name
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: MLS_SWATTACH
    logical :: myDontFail

    ! Executable code
    call trace_begin ( me, 'MLS_swattach_mf' , cond=.false. )
    MLS_SWATTACH = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    select case (MLSFile%HdfVersion)
    case (HDFVERSION_4)
      MLS_SWATTACH = swattach(MLSFile%FileID%f_id, trim(swathName))
    case (HDFVERSION_5)
      MLS_SWATTACH = he5_swattach(MLSFile%FileID%f_id, trim(swathName))
    case default
      MLS_SWATTACH = -1
    end select
    if ( MLS_SWATTACH /= -1 ) then
      call trace_end ( cond=.false. )
      return
    elseif ( myDontFail ) then
      CALL MLSMessage ( MLSMSG_Warning, moduleName,  &
          & 'Failed to attach swath name ' // trim(swathname), &
          & MLSFile=MLSFile )
    else
      CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to attach swath name ' // trim(swathname), &
          & MLSFile=MLSFile )
    endif
    call trace_end ( cond=.false. )
  end function MLS_SWATTACH_MF

  ! --------------------------------------------  MLS_SWCREATE_MF  -----
  function MLS_SWCREATE_MF ( MLSFile, SWATHNAME ) &
    & result(MLS_SWCREATE)
    use hdf5, only: h5eSet_auto_f
    type (MLSFile_T)   :: MLSFile
    character(len=*), intent(in) :: SWATHNAME       ! Swath name
    integer :: MLS_SWCREATE

    ! Internal variables
    logical :: alreadyThere
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: status

    logical, parameter :: ALWAYSTRYSWATTACH = .true.
    logical, parameter :: MUSTCREATE = .true.

    ! Executable code
    call trace_begin ( me, 'MLS_swcreate_mf' , cond=.false. )
    MLS_SWCREATE = 0
    ! All necessary input supplied?
    select case (MLSFile%HdfVersion)
    case (HDFVERSION_4)
      if ( MUSTCREATE ) then
        alreadyThere = .false.
      elseif ( ALWAYSTRYSWATTACH ) then
        if(DEEBUG) print *, 'About to call swattach with FileID: ', &
          & MLSFile%FileID%f_id, ' swathname ', &
          & trim(swathName)
        alreadyThere = (swattach(MLSFile%FileId%f_id, trim(swathName)) >= 0)
      else
        if(DEEBUG) print *, 'About to call MLS_swath_in_file'
        alreadyThere = MLS_swath_in_file(MLSFile%Name, swathName, HDFVERSION_4)
      endif
      if ( alreadyThere ) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Swathname call to MLS_SWCREATE already there: ' // trim(swathName), &
          & MLSFile=MLSFile )
      endif
      if(DEEBUG) print *, 'About to call swcreate with FileID: ', &
        &  MLSFile%Fileid%f_id, ' swathname ', &
        & trim(swathName)
      MLS_SWCREATE = swcreate(MLSFile%FileID%f_id, trim(swathName))
    case (HDFVERSION_5)
      call h5eSet_auto_f ( 0, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to turn error messages off in MLS_SWCREATE', MLSFile=MLSFile )
      if ( MUSTCREATE ) then
        alreadyThere = .false.
      elseif ( ALWAYSTRYSWATTACH  ) then
        if(DEEBUG) print *, 'About to call he5_swattach with FileID: ', &
          & MLSFile%FileID%f_id, ' swathname ', &
          & trim(swathName)
        alreadyThere = (he5_swattach(MLSFile%FileID%f_id, trim(swathName)) >= 0)
      else
        if(DEEBUG) print *, 'About to call MLS_swath_in_file'
        alreadyThere = MLS_swath_in_file(MLSFile%Name, swathName, HDFVERSION_5)
      endif
      call h5eSet_auto_f ( 1, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to turn error messages back on in MLS_SWCREATE')
      if ( alreadyThere ) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Swathname call to MLS_SWCREATE already there: ' // trim(swathName))
      endif
      MLS_SWCREATE = he5_swcreate(MLSFile%FileID%f_id, trim(swathName))
    case default
      MLS_SWCREATE = -1
    end select
    if ( MLS_SWCREATE == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to create swath name ' // trim(swathname), &
          & MLSFile=MLSFile )
    call trace_end ( cond=.false. )

  end function MLS_SWCREATE_MF

  ! --------------------------------------------  MLS_SWCREATE_ID  -----
  function MLS_SWCREATE_ID ( FILEID, SWATHNAME, FileName, hdfVersion ) &
    & result(MLS_SWCREATE)
    use hdf5, only: h5eSet_auto_f
    integer, intent(in) :: FILEID      ! ID returned by MLS_swopen
    character(len=*), intent(in) :: SWATHNAME       ! Swath name
    character(len=*), optional, intent(in) :: FILENAME  ! File name
    integer, optional, intent(in) :: hdfVersion
    integer  :: MLS_SWCREATE

    ! Internal variables
    logical :: alreadyThere
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: myHdfVersion
    logical :: needsFileName
    integer :: status
    logical, parameter :: ALWAYSTRYSWATTACH = .true.
    logical, parameter :: MUSTCREATE = .true.

    ! Executable code
    call trace_begin ( me, 'MLS_swcreate_id' , cond=.false. )
    MLS_SWCREATE = 0
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FileName)) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWCREATE' )
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      if ( MUSTCREATE ) then
        alreadyThere = .false.
      elseif ( .not. present(filename) .or. ALWAYSTRYSWATTACH ) then
        if(DEEBUG) print *, 'About to call swattach with FileID: ', Fileid, ' swathname ', &
          & trim(swathName)
        alreadyThere = (swattach(FileID, trim(swathName)) >= 0)
      else
        if(DEEBUG) print *, 'About to call MLS_swath_in_file'
        alreadyThere = MLS_swath_in_file(FileName, swathName, HDFVERSION_4)
      endif
      if ( alreadyThere ) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Swathname call to MLS_SWCREATE already there: ' // trim(swathName))
      endif
      if(DEEBUG) print *, 'About to call swcreate with FileID: ', Fileid, ' swathname ', &
        & trim(swathName)
      MLS_SWCREATE = swcreate(Fileid, trim(swathName))
    case (HDFVERSION_5)
      call h5eSet_auto_f ( 0, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to turn error messages off in MLS_SWCREATE')
      if ( MUSTCREATE ) then
        alreadyThere = .false.
      elseif ( .not. present(filename) .or. ALWAYSTRYSWATTACH  ) then
        if(DEEBUG) print *, 'About to call he5_swattach with FileID: ', Fileid, ' swathname ', &
          & trim(swathName)
        alreadyThere = (he5_swattach(FileID, trim(swathName)) >= 0)
      else
        if(DEEBUG) print *, 'About to call MLS_swath_in_file'
        alreadyThere = MLS_swath_in_file(FileName, swathName, HDFVERSION_5)
      endif
      call h5eSet_auto_f ( 1, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to turn error messages back on in MLS_SWCREATE')
      if ( alreadyThere ) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Swathname call to MLS_SWCREATE already there: ' // trim(swathName))
      endif
      MLS_SWCREATE = he5_swcreate(Fileid, trim(swathName))
    case default
      MLS_SWCREATE = -1
    end select
    if ( MLS_SWCREATE == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to create swath name ' // trim(swathname) )
    call trace_end ( cond=.false. )

  end function MLS_SWCREATE_ID

  ! -----------------------------------------------  MLS_SWdefdim  -----
  integer function MLS_SWdefdim ( SWATHID, DIMNAME, DIMSIZE, FILENAME, &
    & hdfVersion, DONTFAIL )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    integer, intent(in)          :: DIMSIZE
    character(len=*), optional :: FileName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_swdefdim' , cond=.false. )
    MLS_SWdefdim = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWdefdim = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWdefdim' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWdefdim = swdefdim(swathid, DIMNAME, DIMSIZE)
    case (HDFVERSION_5)
        MLS_SWdefdim = HE5_SWdefdim(swathid, DIMNAME, hsize(DIMSIZE))
    case default
      MLS_SWdefdim = -1
    end select
    if ( .not. myDontFail .and. MLS_SWdefdim == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to define dimension ' // trim(dimname) )
    call trace_end ( cond=.false. )

  end function MLS_SWdefdim

  ! -----------------------------------------------  MLS_SWDETACH  -----
  integer function MLS_SWDETACH ( SWATHID, FileName, hdfVersion )
    integer, intent(in) :: SWATHID      ! ID returned by MLS_swattach
    character(len=*), optional, intent(in) :: FILENAME  ! File name
    integer, optional, intent(in) :: hdfVersion

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_swdetach' , cond=.false. )
    if (Deebug) print *, 'swdetaching ', swathid
    MLS_SWDETACH = 0
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWDETACH' )
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWDETACH = swdetach(swathid)
    case (HDFVERSION_5)
      MLS_SWDETACH = he5_swdetach(swathid)
    case default
      MLS_SWDETACH = -1
    end select
    if ( MLS_SWDETACH == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to detach swath id ' )
    call trace_end ( cond=.false. )

  end function MLS_SWDETACH

  ! ----------------------------------------------  MLS_SWdiminfo  -----
  integer function MLS_SWdiminfo ( SWATHID, DIMNAME, FILENAME, &
    & hdfVersion, DONTFAIL )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    character(len=*), optional :: FileName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_swdiminfo' , cond=.false. )
    MLS_SWdiminfo = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWdiminfo = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWdiminfo' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWdiminfo = SWdiminfo(swathid, DIMNAME)
    case (HDFVERSION_5)
        MLS_SWdiminfo = HE5_SWdiminfo(swathid, DIMNAME)
    case default
      MLS_SWdiminfo = -1
    end select
    if ( .not. myDontFail .and. MLS_SWdiminfo == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to get info on dimension ' // trim(dimname) )
    call trace_end ( cond=.false. )

  end function MLS_SWdiminfo

  ! ----------------------------------------------  MLS_DFLDSETUP  -----
  integer function MLS_DFLDSETUP ( SWATHID, FIELDNAME, DIMNAME, MAXDIMList, &
    & DATATYPE, MERGE, CHUNK_RANK, CHUNK_DIMS, &
    & FILENAME, hdfVersion, DONTFAIL, iFill, rFill, dFill )
    integer, parameter :: RANK = 7
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    character(len=*), intent(in) :: MAXDIMLIST
    integer, intent(in)          :: DATATYPE
    integer, intent(in)          :: MERGE
    integer, intent(in)          :: CHUNK_RANK
    integer, dimension(RANK), intent(in) :: CHUNK_DIMS
    character(len=*), optional :: FileName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    integer, optional, intent(in) :: iFill
    real, optional, intent(in) :: rFill
    double precision, optional, intent(in) :: dFill

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_dfldsetup' , cond=.false. )
    MLS_dfldsetup = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_dfldsetup = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_dfldsetup' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_dfldsetup = swdefdfld(swathid, FIELDName, DIMNAME, Datatype, &
         & MERGE)
    case (HDFVERSION_5)
      if ( chunk_rank /= 0 ) &
        & MLS_dfldsetup = HE5_SWdefchunk( swathid, chunk_rank, hsizes(chunk_dims) )
      if ( MLS_dfldsetup == 0 ) then
        if ( present(iFill) ) then
          MLS_dfldsetup = HE5_SWsetfill(swathid, trim(fieldname), &
            & he2he5_DataType(Datatype), iFill)
        elseif ( present(rFill) ) then
          MLS_dfldsetup = HE5_SWsetfill(swathid, trim(fieldname), &
            & he2he5_DataType(Datatype), rFill)
        elseif ( present(dFill) ) then
          MLS_dfldsetup = HE5_SWsetfill(swathid, trim(fieldname), &
            & he2he5_DataType(Datatype), dFill)
        endif
      endif
      if ( MLS_dfldsetup == 0 ) &
        & MLS_dfldsetup = HE5_SWdefdfld(swathid, FIELDName, DIMNAME, MAXDIMLIST, &
        & he2he5_DataType(Datatype), MERGE)
    case default
      MLS_dfldsetup = -1
    end select
    if ( DEEBUG ) then
      print *, 'MLS_dfldsetup: FIELDName, DIMNAME, MAXDIMLIST, TYPE, MERGE'
      print *, FIELDName, DIMNAME, trim(MAXDIMLIST), DATATYPE, MERGE
      print *, chunk_rank, chunk_dims
    endif
    if ( .not. myDontFail .and. MLS_dfldsetup == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to set up data field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_DFLDSETUP

  ! ----------------------------------------------  MLS_GFLDSETUP  -----
  integer function MLS_GFLDSETUP ( SWATHID, FIELDNAME, DIMNAME, MAXDIMList, &
    & DATATYPE, MERGE, CHUNK_RANK, CHUNK_DIMS, &
    & FILENAME, hdfVersion, DONTFAIL, iFill, rFill, dFill )
    integer, parameter :: RANK = 7
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    character(len=*), intent(in) :: MAXDIMLIST
    integer, intent(in)          :: DATATYPE
    integer, intent(in)          :: MERGE
    integer, intent(in)          :: CHUNK_RANK
    integer, dimension(RANK), intent(in) :: CHUNK_DIMS
    character(len=*), optional :: FileName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL
    integer, optional, intent(in) :: iFill
    real, optional, intent(in) :: rFill
    double precision, optional, intent(in) :: dFill

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_gfldsetup' , cond=.false. )
    MLS_gfldsetup = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_gfldsetup = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_gfldsetup' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_gfldsetup = swdefgfld(swathid, trim(FIELDName), trim(DIMNAME), &
        & Datatype, MERGE)
    case (HDFVERSION_5)
      if ( chunk_rank /= 0 ) &
        & MLS_gfldsetup = HE5_SWdefchunk( swathid, chunk_rank, hsizes(chunk_dims) )
      if ( MLS_gfldsetup == 0 ) then
        if ( present(iFill) ) then
          MLS_gfldsetup = HE5_SWsetfill(swathid, trim(fieldname), &
            & he2he5_DataType(Datatype), iFill)
        elseif ( present(rFill) ) then
          MLS_gfldsetup = HE5_SWsetfill(swathid, trim(fieldname), &
            & he2he5_DataType(Datatype), rFill)
        elseif ( present(dFill) ) then
          MLS_gfldsetup = HE5_SWsetfill(swathid, trim(fieldname), &
            & he2he5_DataType(Datatype), dFill)
        endif
      endif
      if ( MLS_gfldsetup == 0 ) &
        & MLS_gfldsetup = HE5_SWdefgfld(swathid, trim(FIELDName), &
        & trim(DIMNAME), trim(MAXDIMLIST), &
        & he2he5_DataType(Datatype), MERGE)
    case default
      MLS_gfldsetup = -1
    end select
    if ( DEEBUG ) then
      print *, 'MLS_gfldsetup: FIELDName, DIMNAME, MAXDIMLIST, TYPE, MERGE'
      print *, FIELDName, DIMNAME, trim(MAXDIMLIST), DATATYPE, MERGE
      print *, chunk_rank, chunk_dims
    endif
    if ( .not. myDontFail .and. MLS_gfldsetup == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to set up geoloc field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_GFLDSETUP

  ! ----------------------------------------  MLS_SWRDFLD_CHAR_1D  -----
  integer function MLS_SWRDFLD_CHAR_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    character(len=*), dimension(:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    ! Declare these as externals to try to fool hdfeos(5)
    integer, external :: swrdfld
    integer, external :: he5_swrdfld

    ! Executable code
    call trace_begin ( me, 'MLS_SWRDFLD_CHAR_1D' , cond=.false. )
    MLS_swrdfld_CHAR_1D = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_swrdfld_CHAR_1D = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_CHAR_1D' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    if (myDontFail) then
      if (.not. is_datafield_in_swath(swathid, trim(fieldname), myHdfVersion) ) then
        call trace_end ( cond=.false. )
        return
      end if
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      if (myDontFail) then
        if (.not. is_swath_datatype_right(swathid, trim(fieldname), &
          & DFNT_CHAR8, myHdfVersion) ) then
          call trace_end ( cond=.false. )
          return
        end if
      endif
      MLS_swrdfld_CHAR_1D = SWRDFLD(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_swrdfld_CHAR_1D = HE5_SWRDFLD(swathid, trim(fieldname), &
        & start, stride, edge, values)
      if (myDontFail) then
        if (.not. is_swath_datatype_right(swathid, trim(fieldname), &
          & MLS_CHARTYPE, myHdfVersion) ) then
          call trace_end ( cond=.false. )
          return
        end if
      endif
    case default
      MLS_swrdfld_CHAR_1D = -1
    end select
    if ( .not. myDontFail .and. MLS_swrdfld_CHAR_1D == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1-d char field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWRDFLD_CHAR_1D

  ! --------------------------------------  MLS_SWRDFLD_DOUBLE_1D  -----
  integer function MLS_SWRDFLD_DOUBLE_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWRDFLD_DOUBLE_1D' , cond=.false. )
    MLS_swrdfld_double_1d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_swrdfld_double_1d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_DOUBLE_1D' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_swrdfld_double_1d = SWRDFLD_DOUBLE(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_swrdfld_double_1d = HE5_SWRDFLD_DOUBLE(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_swrdfld_double_1d = -1
    end select
    if ( .not. myDontFail .and. MLS_swrdfld_double_1d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1d double field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWRDFLD_DOUBLE_1D

  ! --------------------------------------  MLS_SWRDFLD_DOUBLE_2d  -----
  integer function MLS_SWRDFLD_DOUBLE_2d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 2
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:,:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWRDFLD_DOUBLE_2D' , cond=.false. )
    MLS_swrdfld_double_2d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_swrdfld_double_2d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_DOUBLE_2d' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_swrdfld_double_2d = SWRDFLD_DOUBLE_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_swrdfld_double_2d = HE5_SWRDFLD_DOUBLE_2D(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_swrdfld_double_2d = -1
    end select
    if ( .not. myDontFail .and. MLS_swrdfld_double_2d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 2-d double field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWRDFLD_DOUBLE_2d

  ! --------------------------------------  MLS_SWRDFLD_DOUBLE_3d  -----
  integer function MLS_SWRDFLD_DOUBLE_3d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 3
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:,:,:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWRDFLD_DOUBLE_3D' , cond=.false. )
    MLS_swrdfld_double_3d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_swrdfld_double_3d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_DOUBLE_3d' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_swrdfld_double_3d = SWRDFLD_DOUBLE_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_swrdfld_double_3d = HE5_SWRDFLD_DOUBLE_3D(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_swrdfld_double_3d = -1
    end select
    if ( .not. myDontFail .and. MLS_swrdfld_double_3d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 3d double field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWRDFLD_DOUBLE_3d

  ! ----------------------------------------  MLS_SWRDFLD_INTEGER  -----
  integer function MLS_SWRDFLD_INTEGER ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    integer, dimension(:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWRDFLD_INTEGER' , cond=.false. )
    MLS_swrdfld_integer = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_swrdfld_integer = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_integer' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_swrdfld_integer = SWRDFLD_INTEGER(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_swrdfld_integer = HE5_SWRDFLD_INTEGER(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_swrdfld_integer = -1
    end select
    if ( .not. myDontFail .and. MLS_swrdfld_integer == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read integer field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWRDFLD_INTEGER

  ! ----------------------------------------  MLS_SWRDFLD_REAL_1D  -----
  integer function MLS_SWRDFLD_REAL_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWRDFLD_REAL_1D' , cond=.false. )
    MLS_swrdfld_REAL_1d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_swrdfld_REAL_1d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_REAL_1D' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_swrdfld_REAL_1d = SWRDFLD_REAL(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_swrdfld_REAL_1d = HE5_SWRDFLD_REAL(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_swrdfld_REAL_1d = -1
    end select
    if ( .not. myDontFail .and. MLS_swrdfld_REAL_1d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 1d real field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWRDFLD_REAL_1D

  ! ----------------------------------------  MLS_SWRDFLD_REAL_2d  -----
  integer function MLS_SWRDFLD_REAL_2d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
  ! use hdf5, only: hsize_t
    integer, parameter :: RANK = 2
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:,:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWRDFLD_REAL_2D' , cond=.false. )
    MLS_swrdfld_REAL_2d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_swrdfld_REAL_2d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_REAL_2d' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_swrdfld_REAL_2d = SWRDFLD_REAL_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      if ( DEEBUG ) then
        print *, 'About to read 2-d field: ', trim(fieldname)
        print *, 'start:  ', start
        print *, 'stride: ', stride
        print *, 'edge: ', edge
        print *, 'hsize_t(start):  ', hsizes(start)
        print *, 'hsize_t(stride): ', hsizes(stride)
        print *, 'hsize_t(edge): ', hsizes(edge)
      endif
      MLS_swrdfld_REAL_2d = HE5_SWRDFLD_REAL_2D(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_swrdfld_REAL_2d = -1
    end select
    if ( .not. myDontFail .and. MLS_swrdfld_REAL_2d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 2-d REAL field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWRDFLD_REAL_2d

  ! ----------------------------------------  MLS_SWRDFLD_REAL_3d  -----
  integer function MLS_SWRDFLD_REAL_3d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 3
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:,:,:), intent(out) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWRDFLD_REAL_3D' , cond=.false. )
    MLS_swrdfld_REAL_3d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_swrdfld_REAL_3d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWRDFLD_REAL_3d' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_swrdfld_REAL_3d = SWRDFLD_REAL_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_swrdfld_REAL_3d = HE5_SWRDFLD_REAL_3D(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_swrdfld_REAL_3d = -1
    end select
    if ( .not. myDontFail .and. MLS_swrdfld_REAL_3d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to read 3d real field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWRDFLD_REAL_3d

  ! -----------------------------------------------  MLS_swwrattr  -----
  integer function MLS_swwrattr ( SWATHID, &
    & ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_CHARTYPE
    integer, intent(in) :: COUNT   ! How many
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'MLS_swwrattr' , cond=.false. )
    if ( len_trim(buffer) > 0 ) then
      MLS_swwrattr = HE5_SWWRATTR( SWATHID, &
      & ATTRNAME, DATATYPE, hsize(max(COUNT, len_trim(BUFFER))), BUFFER )
    else
      MLS_swwrattr = HE5_SWWRATTR( SWATHID, &
      & ATTRNAME, DATATYPE, hsize(1), BLANK )
    endif
    call trace_end ( cond=.false. )

  end function MLS_swwrattr

  ! ----------------------------------------------  MLS_swwrlattr  -----
  integer function MLS_swwrlattr ( SWATHID, &
    & FIELDNAME, ATTRNAME, DATATYPE, COUNT, BUFFER )
    integer, intent(in) :: SWATHID      ! Swath ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_CHARTYPE
    integer, intent(in) :: COUNT   ! How many
    character(len=*), intent(in) :: BUFFER  ! Buffer for write

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'MLS_swwrlattr' , cond=.false. )
    if ( len_trim(buffer) > 0 ) then
      MLS_swwrlattr = HE5_SWWRLATTR( SWATHID, FIELDNAME, &
      & ATTRNAME, DATATYPE, hsize(max(COUNT, len_trim(BUFFER))), BUFFER )
    else
      MLS_swwrlattr = HE5_SWWRLATTR( SWATHID, FIELDNAME, &
      & ATTRNAME, DATATYPE, hsize(1), BLANK )
    endif
    call trace_end ( cond=.false. )

  end function MLS_swwrlattr

  ! ----------------------------------------  MLS_SWWRFLD_CHAR_1D  -----
  integer function MLS_SWWRFLD_CHAR_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    character(len=*), dimension(:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    ! Declare these as externals to try to fool hdfeos(5)
    integer, external :: swwrfld
    integer, external :: he5_swwrfld
    integer, dimension(12) :: dfrank
    integer, dimension(12) :: numbertype
    character(len=80)      :: fieldlist
    integer                :: nflds

    ! Executable code
    call trace_begin ( me, 'MLS_SWWRFLD_CHAR_1D' , cond=.false. )
    MLS_SWWRFLD_CHAR_1D = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWWRFLD_CHAR_1D = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_CHAR_1D' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWWRFLD_CHAR_1D = SWWRFLD(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_SWWRFLD_CHAR_1D = HE5_SWWRFLD(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_SWWRFLD_CHAR_1D = -1
    end select
    needsFileName = is_swath_datatype_right(swathid, fieldname, DFNT_CHAR8, &
      & myHdfVersion, rank_out=dfrank, numbertype_out=numbertype, &
      & fieldlist_out=fieldlist, nflds_out=nflds)
    if ( .not. myDontFail .and. MLS_SWWRFLD_CHAR_1D == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to write 1d char field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWWRFLD_CHAR_1D

  ! --------------------------------------  MLS_SWWRFLD_DOUBLE_1D  -----
  integer function MLS_SWWRFLD_DOUBLE_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWWRFLD_DOUBLE_1D' , cond=.false. )
    MLS_SWWRFLD_double_1d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWWRFLD_double_1d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_DOUBLE_1D' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWWRFLD_double_1d = SWWRFLD_DOUBLE(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_SWWRFLD_double_1d = HE5_SWWRFLD_DOUBLE(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_SWWRFLD_double_1d = -1
    end select
    if ( .not. myDontFail .and. MLS_SWWRFLD_double_1d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to write 1d double field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWWRFLD_DOUBLE_1D

  ! --------------------------------------  MLS_SWWRFLD_DOUBLE_2d  -----
  integer function MLS_SWWRFLD_DOUBLE_2d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 2
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:,:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWWRFLD_DOUBLE_2D' , cond=.false. )
    MLS_SWWRFLD_double_2d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWWRFLD_double_2d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_DOUBLE_2d' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWWRFLD_double_2d = SWWRFLD_DOUBLE_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_SWWRFLD_double_2d = HE5_SWWRFLD_DOUBLE_2D(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_SWWRFLD_double_2d = -1
    end select
    if ( .not. myDontFail .and. MLS_SWWRFLD_double_2d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to write 2-d double field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWWRFLD_DOUBLE_2d

  ! --------------------------------------  MLS_SWWRFLD_DOUBLE_3d  -----
  integer function MLS_SWWRFLD_DOUBLE_3d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 3
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    double precision, dimension(:,:,:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWWRFLD_DOUBLE_3D' , cond=.false. )
    MLS_SWWRFLD_double_3d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWWRFLD_double_3d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_DOUBLE_3d' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWWRFLD_double_3d = SWWRFLD_DOUBLE_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_SWWRFLD_double_3d = HE5_SWWRFLD_DOUBLE_3D(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_SWWRFLD_double_3d = -1
    end select
    if ( .not. myDontFail .and. MLS_SWWRFLD_double_3d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to write 3d double field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWWRFLD_DOUBLE_3d

  ! ----------------------------------------  MLS_SWWRFLD_INTEGER  -----
  integer function MLS_SWWRFLD_INTEGER ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    integer, dimension(:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWWRFLD_INTEGER' , cond=.false. )
    MLS_SWWRFLD_integer = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWWRFLD_integer = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_integer' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWWRFLD_integer = SWWRFLD_INTEGER(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_SWWRFLD_integer = HE5_SWWRFLD_INTEGER(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_SWWRFLD_integer = -1
    end select
    if ( .not. myDontFail .and. MLS_SWWRFLD_integer == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to write integer field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWWRFLD_INTEGER

  ! ----------------------------------------  MLS_SWWRFLD_REAL_1D  -----
  integer function MLS_SWWRFLD_REAL_1D ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 1
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    logical, parameter :: DEEBUG = .false.
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SWWRFLD_REAL_1D' , cond=.false. )
    MLS_SWWRFLD_REAL_1d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWWRFLD_REAL_1d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_REAL_1D' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    if(DEEBUG) print *, 'swathid: ', swathid
    if(DEEBUG) print *, 'myHdfVersion: ', myHdfVersion
    if(DEEBUG) print *, 'fieldname: ', trim(fieldname)
    if(DEEBUG) print *, 'start: ', hsizes(start)
    if(DEEBUG) print *, 'stride: ', hsizes(stride)
    if(DEEBUG) print *, 'edge: ', hsizes(edge)
    if(DEEBUG) print *, 'shape(values): ', shape(values)
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWWRFLD_REAL_1d = SWWRFLD_REAL(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_SWWRFLD_REAL_1d = HE5_SWWRFLD_REAL(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_SWWRFLD_REAL_1d = -1
    end select
    if ( .not. myDontFail .and. MLS_SWWRFLD_REAL_1d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to write 1d real field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWWRFLD_REAL_1D

  ! ----------------------------------------  MLS_SWWRFLD_REAL_2d  -----
  integer function MLS_SWWRFLD_REAL_2d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 2
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:,:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    logical, parameter :: DEEBUG = .false.
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable
    call trace_begin ( me, 'MLS_SWWRFLD_REAL_2D' , cond=.false. )
    MLS_SWWRFLD_REAL_2d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWWRFLD_REAL_2d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_REAL_2d' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    if(DEEBUG) print *, 'swathid: ', swathid
    if(DEEBUG) print *, 'myHdfVersion: ', myHdfVersion
    if(DEEBUG) print *, 'fieldname: ', trim(fieldname)
    if(DEEBUG) print *, 'start: ', hsizes(start)
    if(DEEBUG) print *, 'stride: ', hsizes(stride)
    if(DEEBUG) print *, 'edge: ', hsizes(edge)
    if(DEEBUG) print *, 'shape(values): ', shape(values)
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWWRFLD_REAL_2d = SWWRFLD_REAL_2D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_SWWRFLD_REAL_2d = HE5_SWWRFLD_REAL_2D(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_SWWRFLD_REAL_2d = -1
    end select
    if ( .not. myDontFail .and. MLS_SWWRFLD_REAL_2d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to write 2-d REAL field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWWRFLD_REAL_2d

  ! ----------------------------------------  MLS_SWWRFLD_REAL_3d  -----
  integer function MLS_SWWRFLD_REAL_3d ( SWATHID, FIELDNAME, &
    & START, STRIDE, EDGE, VALUES, FILENAME, hdfVersion, DONTFAIL )
    integer, parameter :: RANK = 3
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    integer, dimension(RANK), intent(in) :: START
    integer, dimension(RANK), intent(in) :: STRIDE
    integer, dimension(RANK), intent(in) :: EDGE
    real, dimension(:,:,:), intent(in) :: VALUES
    character(len=*), optional :: FIleName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: DONTFAIL

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName
    call trace_begin ( me, 'MLS_SWWRFLD_REAL_3D' , cond=.false. )
    MLS_SWWRFLD_REAL_3d = 0
    myDontFail = .false.
    if ( present(DontFail) ) myDontFail = DontFail
    ! All necessary input supplied?
    needsFileName = (.not. present(hdfVersion))
    if ( present(hdfVersion) ) &
      & needsFileName = (hdfVersion == WILDCARDHDFVERSION)
    if ( needsFileName .and. .not. present(FIleName)) then
      if ( myDontFail ) then
        MLS_SWWRFLD_REAL_3d = -1
      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Missing needed arg FILENAME from call to MLS_SWWRFLD_REAL_3d' )
      endif
      call trace_end ( cond=.false. )
      return
    endif
    if ( needsFileName ) then
      myHdfVersion = MLS_hdf_version ( trim(FileName) )
    else
      myHdfVersion = hdfVersion
    endif
    select case (myHdfVersion)
    case (HDFVERSION_4)
      MLS_SWWRFLD_REAL_3d = SWWRFLD_REAL_3D(swathid, trim(fieldname), &
        & start, stride, edge, values)
    case (HDFVERSION_5)
      MLS_SWWRFLD_REAL_3d = HE5_SWWRFLD_REAL_3D(swathid, trim(fieldname), &
        & hsizes(start), hsizes(stride), hsizes(edge), values)
    case default
      MLS_SWWRFLD_REAL_3d = -1
    end select
    if ( .not. myDontFail .and. MLS_SWWRFLD_REAL_3d == -1 ) &
      & CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'Failed to write 3d real field ' // trim(fieldname) )
    call trace_end ( cond=.false. )

  end function MLS_SWWRFLD_REAL_3d

  ! --------------------------------------  MLS_swath_in_file_sca  -----
  logical function MLS_swath_in_file_sca( filename, swath, HdfVersion, error )
    ! Returns .true. if swath found in file, .false. otherwise
    use HDF5, only: Size_t
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: swath
    integer, intent(in) :: HdfVersion
    integer, optional, intent(out) :: error

    ! Internal variables
    logical, parameter               :: Deebug = .false.
    character(len=MAXDLISTLENGTH)    :: Fieldlist
    integer(kind=size_t)             :: Hlistsize
    integer :: Me = -1               !   String index for trace cacheing
    integer                          :: Listsize
    integer                          :: Nswaths

    ! Begin execution
    ! print*,'Scalar version'
    call trace_begin ( me, 'MLS_swath_in_file_sca' , cond=.false. )
    nswaths = 0
    MLS_swath_in_file_sca = .false.
    fieldlist = ''
    listsize = -1 ! So we can check later whether it has been set
    hlistsize = -1
    select case (HdfVersion)
    case (HDFVERSION_4)
      nswaths = swinqswath(trim(filename), fieldlist, listsize)
    case (HDFVERSION_5)
      nswaths = HE5_swinqswath(trim(filename), fieldlist, hlistsize)
    end select
    if ( deebug ) then
      call outputnamedValue ( 'nswaths', nswaths )
      call outputnamedValue ( 'fieldlist', trim(fieldlist) )
    endif
    if ( present(error) ) error = min(0, nswaths)
    if ( nswaths < 1 ) then
      call trace_end ( cond=.false. )
      return
    endif
    if ( listsize < 0 ) listsize = hlistsize
    if ( listsize > MAXDLISTLENGTH ) then
       CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'list size too big in MLS_swath_in_file_sca ' // trim(filename) )
    elseif ( listsize < MAXDLISTLENGTH .and. listsize > 0 ) then
      fieldlist = fieldlist(1:listsize) // ' '
    endif
    MLS_swath_in_file_sca = &
      & ( StringElementNum(fieldlist, trim(swath), .true.) > 0 )
    call trace_end ( cond=.false. )
  end function MLS_swath_in_file_sca

  ! --------------------------------------  MLS_swath_in_file_arr  -----
  logical function MLS_swath_in_file_arr(filename, swaths, HdfVersion, &
    & which, error )
    use HDF5, only: Size_t
    ! Array version of the above
    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: swaths
    integer, intent(in) :: HdfVersion
    logical, dimension(:), intent(out) :: which
    integer, optional, intent(out) :: error

    ! Internal variables
    character(len=MAXDLISTLENGTH)    :: Fieldlist
    integer(kind=size_t)             :: Hlistsize
    integer                          :: I
    integer                          :: Listsize
    integer :: Me = -1               !   String index for trace cacheing
    integer                          :: Nswaths

    ! Begin execution
    call trace_begin ( me, 'MLS_swath_in_file_arr' , cond=.false. )
    nswaths = 0
    MLS_swath_in_file_arr = .false.
    which = .false.
    fieldlist = ' '
    if ( size(swaths) > size(which) ) &
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'array to small to hold values in MLS_swath_in_file_arr' )
    select case (HdfVersion)
    case (HDFVERSION_4)
      nswaths = swinqswath(trim(filename), fieldlist, listsize)
    case (HDFVERSION_5)
      nswaths = HE5_swinqswath(trim(filename), fieldlist, hlistsize)
      listsize = hlistsize
    end select
    if ( present(error) ) error = min(0, nswaths)
    if ( nswaths < 1 ) then
      call trace_end ( cond=.false. )
      return
    endif
    if (Deebug) print *, ' nswaths is ', nswaths
    if (Deebug) print *, ' listsize is ', listsize
    if (Deebug) print *, ' fieldlist is ', trim(fieldlist)
    if ( listsize > MAXDLISTLENGTH ) then
       CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'list size too big in MLS_swath_in_file_arr ' // trim(filename) )
    elseif ( listsize < MAXDLISTLENGTH .and. listsize > 0 ) then
      fieldlist = fieldlist(1:listsize) // ' '
    endif
    do i=1, size(swaths)
      which(i) = &
      & ( StringElementNum(fieldlist, trim(swaths(i)), .true.) > 0 )
    enddo
    MLS_swath_in_file_arr = any( which )
    call trace_end ( cond=.false. )
  end function MLS_swath_in_file_arr

  ! ---------------- hsize ------------
  function hsize ( arg ) result ( h )
    use hdf5, only: size_t
    ! Return arg with same integer value
    ! but with kind value size_t or, someday, hdfeos5size_t
    integer,  intent(in)               :: arg
    integer(kind=size_t)               :: h
    h = int(arg, size_t)
  end function hsize

  ! ---------------- hsizes ------------
  function hsizes ( ints ) result ( h )
    use hdf5, only: size_t
    ! Return array with same integer values as ints
    ! but with kind value size_t or, someday, hdfeos5size_t
    integer, dimension(:), intent(in)               :: ints
    integer(kind=size_t), dimension(size(ints))     :: h
    integer :: i
    h = 0
    do i=1, size(ints)
      h(i) = int(ints(i), size_t)
    enddo
  end function hsizes

! ======================= Private Procedures ===========================

  ! --------------------------------------  is_datafield_in_swath  -----
  logical function is_datafield_in_swath(swathid, field, HdfVersion)
    ! Returns .true. if datafield found in swath, .false. otherwise
    integer, intent(in) :: swathid
    character(len=*), intent(in) :: field
    integer, intent(in) :: HdfVersion
    ! Internal variables
    integer, dimension(MAXNODFIELDS) :: Rank
    integer, dimension(MAXNODFIELDS) :: Numbertype
    character(len=MAXDLISTLENGTH)    :: Fieldlist
    integer                          :: Nflds
    ! Begin execution
    nflds = 0
    is_datafield_in_swath = .false.
    fieldlist = ''
    select case (HdfVersion)
    case (HDFVERSION_4)
      nflds = swinqdflds(swathid, fieldlist, rank, numbertype)
    case (HDFVERSION_5)
      nflds = HE5_swinqdflds(swathid, fieldlist, rank, numbertype)
    end select
    if ( nflds == 0 ) return
    is_datafield_in_swath = &
      & ( StringElementNum(fieldlist, trim(field), .true.) > 0 )
  end function is_datafield_in_swath

  ! ------------------------------------  is_swath_datatype_right  -----
  logical function is_swath_datatype_right(swathid, field, datatype, &
    & HdfVersion, rank_out, numbertype_out, fieldlist_out, nflds_out)
    integer, intent(in)                          :: Swathid
    character(len=*), intent(in)                 :: Field
    integer, intent(in)                          :: Datatype
    integer, intent(in)                          :: HdfVersion
    integer, dimension(:), intent(out), optional :: Rank_out
    integer, dimension(:), intent(out), optional :: Numbertype_out
    character(len=*), intent(out), optional      :: Fieldlist_out
    integer, intent(out), optional               :: Nflds_out
    ! Internal variables
    integer, dimension(MAXNODFIELDS) :: rank
    integer, dimension(MAXNODFIELDS) :: numbertype
    character(len=MAXDLISTLENGTH)    :: fieldlist
    integer                         :: nflds
    ! Begin execution
    nflds = 0
    is_swath_datatype_right = .false.
    ! Set defaults for inputs to HDF routines so that
    ! we get appropriate answers at all times.
    rank = 0
    fieldlist = ''
    numbertype = 0
    select case (HdfVersion)
    case (HDFVERSION_4)
      nflds = swinqdflds(swathid, fieldlist, rank, numbertype)
    case (HDFVERSION_5)
      nflds = HE5_swinqdflds(swathid, fieldlist, rank, numbertype)
    end select
    ! Set defaults for ouptuts
    if ( present(nflds_out) ) nflds_out = nflds
    if ( present(rank_out) ) rank_out = 0
    if ( present(numbertype_out) ) numbertype_out = 0
    if ( present(fieldlist_out) ) fieldlist_out = ' '
    if ( nflds == 0 ) return
    ! Now update our outputs
    if ( present(rank_out) ) rank_out = rank(1:size(rank_out))
    if ( present(numbertype_out) ) &
      & numbertype_out = numbertype(1:size(numbertype_out))
    if ( present(fieldlist_out) ) fieldlist_out = fieldlist
    nflds = StringElementNum(fieldlist, trim(field), .true.)
    is_swath_datatype_right = ( nflds > 0 .and. nflds <= MAXNODFIELDS )
    if ( .not. is_swath_datatype_right ) return
    is_swath_datatype_right = &
      & ( numbertype(nflds) == datatype )
  end function is_swath_datatype_right

  ! --------------------------------------------  he2he5_DataType  -----

  ! This function converts hdfeos2 datatypes to
  ! corresponding hdfeos5 numbers
  ! (unless they are already HE5 types in which it returns them unchanged)

  function he2he5_DataType(dataType) result (HE5_dataType)

    ! Arguments
    integer, intent(IN)       :: dataType
    integer                   :: HE5_dataType
    ! begin
    HE5_dataType = dataType
    select case (dataType)
    case(DFNT_FLOAT32)
      HE5_dataType = HE5T_NATIVE_FLOAT
    case(DFNT_FLOAT64)
      HE5_dataType = HE5T_NATIVE_DOUBLE
    case(DFNT_CHAR8)
      HE5_dataType = MLS_CHARTYPE
    case(DFNT_INT8)
      HE5_dataType = HE5T_NATIVE_INT8
    case(DFNT_INT16)
      HE5_dataType = HE5T_NATIVE_INT16
    case(DFNT_INT32)
      HE5_dataType = HE5T_NATIVE_INT
    case(DFNT_INT64)
      HE5_dataType = HE5T_NATIVE_INT64
    !case default
      !HE5_dataType = HE5T_NATIVE_INT
    end select
  end function he2he5_DataType

  ! ----------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSHDFEOS

! $Log$
! Revision 2.48  2019/01/29 21:44:15  pwagner
! Initialized attrList in MLS_IsGlatt_fid to prevent bleed-thru from last call
!
! Revision 2.47  2018/03/22 16:56:35  pwagner
! Debug printing in MLS_swath_in_file_sca; CamelCase use statements
!
! Revision 2.46  2013/08/31 01:24:53  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end
!
! Revision 2.45  2013/06/12 02:11:27  vsnyder
! Cruft removal
!
! Revision 2.44  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.43  2010/01/29 01:16:06  pwagner
! Fixed undefined hlistsize bug Lahey complained about
!
! Revision 2.42  2010/01/14 23:29:35  pwagner
! Uses separate scalar and array generic forms for HE5_EHRDGLATT_INTEGER
!
! Revision 2.41  2010/01/11 18:34:50  pwagner
! Added more debug printing
!
! Revision 2.40  2009/10/05 23:38:21  pwagner
! Moved use hdf5 statements from module scope to speedup Lahey; this is the last time we do that
!
! Revision 2.39  2009/09/29 23:33:49  pwagner
! Changes needed by 64-bit build
!
! Revision 2.38  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.37  2008/05/24 00:54:46  vsnyder
! Remove unused declarations
!
! Revision 2.36  2008/05/02 00:05:36  pwagner
! Reads textFile using io stuff
!
! Revision 2.35  2008/04/25 22:51:44  pwagner
! Passes optional arg maxLineLen through to textFile_to_Chars
!
! Revision 2.34  2008/04/18 16:29:30  pwagner
! Now works properly with NAG, Lahey, and Intel
!
! Revision 2.33  2008/03/07 01:36:08  pwagner
! Will subdivide attributes into blocks if too large
!
! Revision 2.32  2008/02/22 21:30:18  pwagner
! Can now save entire textfile as global attribute
!
! Revision 2.31  2007/08/20 22:00:47  pwagner
! More procedures push their names onto MLSCallStack
!
! Revision 2.30  2007/08/17 00:27:48  pwagner
! push more procedures onto MLSCallStack
!
! Revision 2.29  2005/11/15 00:18:20  pwagner
! present a poor choice for argument name; esp when optional args exist
!
! Revision 2.28  2005/11/11 21:41:03  pwagner
! MLS_swath_in_file_sca optionally returns an error flag, too
!
! Revision 2.27  2005/10/11 17:29:15  pwagner
! Added MLS_ISGLATT function
!
! Revision 2.26  2005/06/29 00:39:30  pwagner
! New interfaces for MLS_SWCREATE MLS_SWATTACH accept MLSFiles
!
! Revision 2.25  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.24  2005/06/14 20:34:20  pwagner
! Small changes to guard against zero-length characters
!
! Revision 2.23  2004/09/23 22:58:59  pwagner
! Fix for when passed blank BUFFER (do we need to fix this elsewhere\?)
!
! Revision 2.22  2004/08/04 23:19:01  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.21  2004/07/22 17:07:15  pwagner
! Fixed set fill values
!
! Revision 2.20  2004/05/05 21:28:42  pwagner
! More debug printing
!
! Revision 2.19  2004/03/24 23:53:02  pwagner
! Switched from HE5T_NATIVE_SCHAR to MLS_CHARTYPE
!
! Revision 2.18  2004/02/13 00:16:39  pwagner
! New stuff for reading swath attributes
!
! Revision 2.17  2004/02/05 23:31:05  pwagner
! Extra debugging when appropriate
!
! Revision 2.16  2004/01/23 01:12:28  pwagner
! Some care taken in handling ...inq.. functions
!
! Revision 2.15  2003/10/30 00:01:57  pwagner
! Added MLS_GDWRATTR
!
! Revision 2.14  2003/10/28 00:28:53  pwagner
! Added MLS_EHWRGLATT,MLS_SWWRATTR,MLS_SWWRLATTR
!
! Revision 2.13  2003/07/15 23:36:28  pwagner
! Disabled most printing; trims args to (HE5_)SWdefgfld
!
! Revision 2.12  2003/07/11 21:50:31  livesey
! Minor bug fix, probably of little consquence
!
! Revision 2.11  2003/07/11 01:22:55  livesey
! Bug fixes and tidyups
!
! Revision 2.10  2003/07/09 21:49:07  pwagner
! Wont try swattaching just to see if swath already there
!
! Revision 2.9  2003/07/02 00:54:25  pwagner
! Failed attempt to turn off hdfeos5 warnings by turning off hdf5 ones
!
! Revision 2.8  2003/06/26 00:05:40  pwagner
! Added optional DONTFAIL arg to MLS_SWATTACH
!
! Revision 2.7  2003/06/20 19:31:39  pwagner
! Changes to allow direct writing of products
!
! Revision 2.6  2003/06/06 22:49:12  pwagner
! Added MLS_sw(gd)create
!
! Revision 2.5  2003/04/21 19:32:27  pwagner
! Can read/write 1-d char fields
!
! Revision 2.4  2003/04/17 23:04:51  pwagner
! Now sits between L2GPData and HDFEOS(5)
!
! Revision 2.3  2003/04/15 23:16:46  pwagner
! Removed prints; begun MLS_SWRDFLD stuff
!
! Revision 2.2  2003/04/15 21:58:54  pwagner
! Now sets _FillValue attribute because swsetfill seems broken
!
! Revision 2.1  2003/04/11 23:28:44  pwagner
! FIrst commit
!
