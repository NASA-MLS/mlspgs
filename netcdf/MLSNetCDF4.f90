! Copyright 2020, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MLSNetCDF4

  ! This module contains MLS specific routines to do common NetCDF
  ! tasks not specifically found in the NetCDF module.
  
  ! Many concern what is equivalent to an hdfeos swath
  ! but implemented with NetCDF groups, variables, etc.
  ! hdfeos                 NetCDF
  ! -----------            ----------------
  ! swath name             named group 
  ! Data Fields            'Data Fields' group
  ! L2gpdata               'L2gpdata' var
  ! L2gpprecision          'L2gpprecision' var
  !  .    .    .             .    .    .
  ! Geolocation Fields     'geolocation Fields' group
  ! ChunkNumber            'ChunkNumber' var
  ! Latitude               'Latitude' var
  !  .    .    .             .    .    .
  

  use NetCDF
  use HDF, only: Dfnt_Char8, Dfnt_Float32, Dfnt_Float64, &
    & Dfnt_Int8, Dfnt_Int16, Dfnt_Int32, Dfnt_Int64
  use HighOutput, only: OutputNamedValue
  use MLSCommon, only: MLSFile_T
  use MLSFiles, only: HDFversion_4, HDFversion_5, WildcardHDFversion, &
    & MLS_HDF_Version
  use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning, MLSMessage
  use MLSStringlists, only: StringElementnum
  use MLSStrings, only: Replace
  use Output_M, only: Output
  use Trace_M, only: Trace_Begin, Trace_End

  implicit none
  private

  ! These have been given the same names as their counterparts
  ! in MLSHDFEOS.
  ! While that's simple and clear, will it be a problem later on?
  public :: MLS_Ehwrglatt, MLS_Isglatt, &
    & MLS_Dfldsetup, MLS_Gfldsetup, &
    & MLS_Swdefdim, MLS_Swdiminfo, MLS_Swrdfld, MLS_Swwrfld, &
    & MLS_Swattach, MLS_Swcreate, MLS_Swdetach, &
    & MLS_Swrdattr, MLS_Swwrattr, MLS_Swrdlattr, MLS_Swwrlattr, &
    & MLS_Swath_In_File
  public :: mls_InqSwath

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! MLS_DfldSetup     Sets up a data field in a swath
! MLS_GfldSetup     Sets up a geolocation field in a swath
! MLS_InqSwath      Return the list of swaths present in the file
! MLS_IsglAtt       Is the named global attribute present in the file?
! MLS_Swath_In_File Is named swath in file?
! MLS_SwAttach      Attaches a swath (perhaps prior to reading)
! MLS_SwCreate      Creates a swath (perhaps prior to writing)
! MLS_SwDefdim      Sets up a dimension for a swath
! MLS_SwDetach      Detaches from a swath (perhaps prior to closing the file)
! MLS_SwDiminfo     Gets info on a dimension for a swath
! MLS_SwRdattr      Reads a swath-level attribute (e.g., pressures)
! MLS_Swrdfld       Reads a field from a swath, data or geolocation
! MLS_SwRdLattr     Reads a datafield-level attribute (e.g, units)
! MLS_SwWrfld       Writes a field to a swath, data or geolocation
! MLS_SwWrattr      Writes a swath-level attribute (e.g., pressures)
!                     if passed FileID instead of swathID: global attributes
! MLS_SwWrlattr     Writes a datafield-level attribute (e.g., units)
! === (end of toc) ===

! === (start of api) ===
! int MLS_DFLDSETUP (int swathID, char* fieldname, char* dimname, 
!     char* MAXDIMList, int datatype, int chunk_rank, 
!     int chunk_dims(:), [char* filename], [int hdfVersion], [log dontFail],
!     [int iFill], [r4 rFill], [r8 dFill]) 
! int MLS_GFLDSETUP (int swathID, char* fieldname, char* dimname, 
!     char* MAXDIMList, int datatype, int chunk_rank, 
!     int chunk_dims(:), [char* filename], [int hdfVersion], [log dontFail],
!     [int iFill], [r4 rFill], [r8 dFill]) 
! log MLS_ISGLATT (file, char* attrname)
!     file can be one of: int fileID, or char* filename
! log MLS_InqSwath (char* filename, char* swathlist, int listSize)
!     swathlist is a comma-separated strlist, listsize is its number of strings
! log MLS_Swath_In_File (file, char* attrname)
! int MLS_SwATTACH (file, char* swathname, [char* filename], [int hdfVersion], 
!     [log dontFail]) 
!     file can be one of: int fileID, or type(MLSFile_T) MLSFile
! int MLS_SwCREATE (file, char* swathname, [char* filename], [int hdfVersion]) 
!     file can be one of: int fileID, or type(MLSFile_T) MLSFile
! int MLS_Swdefdim (int swathID, char* dimname, 
!     int dimsize) 
! int MLS_SwDETACH (int swathID)
! int MLS_Swdiminfo (int swathID, char* dimname, 
!     [char* filename], [int hdfVersion], [log dontFail]) 
! int MLS_SwRdAttr (int swathID, char* attrName, int datatype, int count, 
!     char* buffer)
! int MLS_Swrdfld (int swathID, char* fieldname, 
!     int start(rank), int stride(rank), int edge(rank), values,
!     [char* filename], [int hdfVersion], [log dontFail]) 
!     values can be one of:
!    {char* value, char* value(:), int value(:), r4 value(shp), r8 value(shp)}
!     and shp is an array of ints of size rank
! int MLS_SwWRFLD (int swathID, char* fieldname, 
!     int start(rank), int stride(rank), int edge(rank), values,
!     [char* filename], [int hdfVersion], [log dontFail]) 
!     values can be one of types in MLS_Swrdfld
! int MLS_SwWRATTR (int swathID, char* attrName, int datatype, int count, 
!     char* buffer)
! int MLS_SwWRLATTR (int swathID, char* fieldName, char* attrName, 
!     int datatype, int count, char* buffer)
! === (end of api) ===
  interface MLS_ISGLATT
    module procedure MLS_ISGLATT_FID, &
      & MLS_ISGLATT_FN
  end interface

  interface MLS_Swath_In_File
    module procedure MLS_Swath_In_File_Arr, &
      & MLS_Swath_In_File_Sca, MLS_Swath_In_File_FID
  end interface

  interface MLS_SwATTACH
    module procedure MLS_SwATTACH_ID, &
      & MLS_SwATTACH_MF
  end interface

  interface MLS_SwCREATE
    module procedure MLS_SwCREATE_ID, &
      & MLS_SwCREATE_MF
  end interface

  interface MLS_Swrdfld
    module procedure &
      & MLS_Swrdfld_CHAR_1D, &
      & MLS_Swrdfld_DOUBLE_1D, MLS_Swrdfld_DOUBLE_2D, MLS_Swrdfld_DOUBLE_3D, &
      & MLS_Swrdfld_INTEGER, &
      & MLS_Swrdfld_REAL_1D, MLS_Swrdfld_REAL_2D, MLS_Swrdfld_REAL_3D
  end interface

  interface MLS_SwWRFLD
    module procedure &
      & MLS_SwWRFLD_CHAR_1D, &
      & MLS_SwWRFLD_DOUBLE_1D, MLS_SwWRFLD_DOUBLE_2D, MLS_SwWRFLD_DOUBLE_3D, &
      & MLS_SwWRFLD_INTEGER, &
      & MLS_SwWRFLD_REAL_1D, MLS_SwWRFLD_REAL_2D, MLS_SwWRFLD_REAL_3D
  end interface

  integer, public, save :: NCError
  character(len=*), dimension(2), parameter :: SWGroupNames = (/&
    & 'Data Fields          ', &
    & 'Geolocation Fields   ' &
    & /)
  ! Print debugging stuff?
  logical, parameter          :: DEEBUG = .false.  
  character(len=1), parameter :: BLANK = ' '
  integer, public, parameter  :: MAXNUMSWATHS = 300

contains ! ======================= Public Procedures =========================

  ! --------------------------------------------  MLS_EHwrglAtt  -----
  function MLS_EHwrglAtt ( Fileid, Attrname, value ) result( returnStatus )
    ! Writes the named attribute as global attribute of the file
    ! Args
    integer, intent(in)           :: FileID
    character(len=*), intent(in)  :: Attrname     ! Attribute name
    character(len=*), intent(in)  :: value        ! Attribute value
    integer                       :: returnStatus
    ! Local variables
    logical, parameter            :: Countempty = .true.
    integer                       :: ListSize
    integer                       :: Status
    logical, parameter            :: DEEBUG = .false.
    ! Executable code
    call check ( nf90_put_att( FileID, nf90_global, AttrName, value ), &
      & 'MLS_EHwrglAtt' // trim(Attrname), FailureOK=.true. )
     returnStatus = NCError
  end function MLS_EHwrglAtt

  ! --------------------------------------------  MLS_EHwrglAtt  -----
  function mls_InqSwath ( FileName, SwathList, listSize ) result( returnStatus )
    character(len=*), intent(in)     :: FileName     ! File name
    character(len=*), intent(out)    :: SwathList    ! swath names
    integer, intent(out)             :: listSize
    integer                          :: returnStatus
    integer                          :: FileID
    integer                          :: grpID
    integer                          :: i
    integer, dimension(MAXNUMSWATHS) :: ncids
    character(len=256)               :: grpName
    integer                          :: Me = -1 !String index for trace cacheing
    character, parameter             :: NULL = achar(0)
    !
    call trace_begin ( me, 'mls_InqSwath' , cond=.false. )
    SwathList = ""
    call check( nf90_open(Filename, NF90_NoWrite, FileId) , &
      & 'mls_InqSwath ' // trim(FileName) )
!     print *, 'FileID:   ', FileID
    !
    ! Despite the netcdf docs implying that we must open the "/" root group
    ! the following command returns an error message complaining about an
    ! illegal character. The character is "illegal" even if it's a
    ! BLANK instead of a NULL.
    ! Directly counting the datasets under FileID seems to work
    ! properly so that's how we're implementing it instead.
!     call check( nf90_inq_ncid( FileId, NULL, grpId) , &
!       & 'mls_InqSwath ' // "/", &
!       & FailureOK=.false., silent=.false. )
!     print *, 'grpID:    ', grpID
!     call check ( nf90_inq_grps( grpId, listSize, ncids ), &
!       & 'mls_InqSwath' // '/inq grps', FailureOK=.true. )
    call check ( nf90_inq_grps( FileId, listSize, ncids ), &
      & 'mls_InqSwath' // '/inq grps', FailureOK=.true. )
!     print *, 'listSize: ', listSize
    do i=1, listSize
      call check( nf90_inq_grpname( ncids(i), grpName) , &
        & 'mls_InqSwath ' // "/grpName", &
        & FailureOK=.true., silent=.true. )
      if ( len_trim(SwathList) < 1 ) then
        SwathList = grpName
      else
        SwathList = trim(Swathlist) // "," // grpName
      endif
!       print *, 'grpName: ', trim(grpName)
    enddo
    
    call check( nf90_close(FileId), &
      & 'mls_InqSwath close ', FailureOK=.true. )
    returnStatus = (NCError == nf90_noerr)
    call trace_end ( cond=.false. )
  end function mls_InqSwath

  ! ---------------------------------------------  MLS_IsGlatt_Fn  -----
  function MLS_IsGlatt_Fn ( Filename, Attrname ) result(isThere)
    ! Is the named attribute a global attribute of the named file?
    ! Args
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    logical                      :: isThere

    ! Local variables
    integer :: fileID
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: status

    ! Executable code
    call trace_begin ( me, 'MLS_IsGlatt_Fn' , cond=.false. )
    call check( nf90_open(FILENAME, NF90_NoWrite, FileId), &
      & 'MLS_IsGlatt_Fn open ' // trim(Attrname), FailureOK=.true. )
    isThere = MLS_IsGlatt_fid ( Fileid, Attrname )
    call check( nf90_close(FileId), &
      & 'MLS_IsGlatt_Fn close ' // trim(Attrname), FailureOK=.true. )
    call trace_end ( cond=.false. )
  end function MLS_IsGlatt_Fn

  ! --------------------------------------------  MLS_IsGlatt_fid  -----
  function MLS_IsGlatt_fid ( Fileid, Attrname ) result(isThere)
    ! Is the named attribute a global attribute of the file?
    ! Args
    integer, intent(in)           :: FileID
    character(len=*), intent(in)  :: Attrname     ! Attribute name
    logical                       :: IsThere
    ! Local variables
    logical, parameter            :: Countempty = .true.
    integer                       :: ListSize
    integer                       :: Status
    logical, parameter            :: DEEBUG = .false.
    ! Executable code
    print *, 'Checking for global attribute: ', trim(Attrname)
    call check ( nf90_inquire_attribute( FileID, nf90_global, AttrName), &
      & 'MLS_IsGlatt_fid' // trim(Attrname), FailureOK=.true., silent=.true. )
    isThere = (NCError == nf90_noerr)
  end function MLS_IsGlatt_fid

  ! --------------------------------------------  MLS_SwAttach_ID  -----
  function MLS_SwAttach_ID ( FILEID, SWATHNAME ) result(MLS_SwATTACH)
    integer, intent(in) :: FILEID      ! ID returned by MLS_Swopen
    character(len=*), intent(in) :: SWATHNAME       ! Swath name
    integer :: MLS_SwATTACH

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    logical :: myDontFail
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_SwAttach_ID' , cond=.false. )
    call check( nf90_inq_ncid( FileId, swathName, MLS_SwAttach ), &
      & 'MLS_SwAttach_ID' // trim(swathname) )
    call trace_end ( cond=.false. )
  end function MLS_SwAttach_ID

  ! --------------------------------------------  MLS_Swattach_mf  -----
  function MLS_Swattach_mf ( MLSFile, SWATHNAME ) result(MLS_SwATTACH)
    type (MLSFile_T)   :: MLSFile
    character(len=*), intent(in) :: SWATHNAME       ! Swath name

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: MLS_SwATTACH
    logical :: myDontFail

    ! Executable code
    call trace_begin ( me, 'MLS_Swattach_mf' , cond=.false. )
    MLS_SwAttach = MLS_Swattach_id( MLSFile%FileID%f_id, trim(swathName) )
    call trace_end ( cond=.false. )
  end function MLS_Swattach_mf

  ! --------------------------------------------  MLS_Swcreate_mf  -----
  function MLS_Swcreate_mf ( MLSFile, SWATHNAME ) &
    & result(MLS_SwCREATE)
    use hdf5, only: h5eSet_auto_f
    type (MLSFile_T)   :: MLSFile
    character(len=*), intent(in) :: SWATHNAME       ! Swath name
    integer :: MLS_SwCREATE

    ! Internal variables
    logical :: alreadyThere
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: status

    logical, parameter :: ALWAYSTRYSWATTACH = .true.
    logical, parameter :: MUSTCREATE = .true.

    ! Executable code
    call trace_begin ( me, 'MLS_Swcreate_mf' , cond=.false. )
    MLS_SwCREATE = MLS_Swcreate_id( MLSFile%FileID%f_id, trim(swathName) )
    call trace_end ( cond=.false. )

  end function MLS_Swcreate_mf

  ! --------------------------------------------  MLS_Swcreate_id  -----
  function MLS_Swcreate_id ( FILEID, SWATHNAME ) &
    & result(MLS_SwCREATE)
    use hdf5, only: h5eSet_auto_f
    integer, intent(in) :: FILEID      ! ID returned by MLS_Swopen
    character(len=*), intent(in) :: SWATHNAME       ! Swath name
    integer  :: MLS_SwCREATE

    ! Internal variables
    logical :: alreadyThere
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: dfid
    integer :: glid

    ! Executable code
    call trace_begin ( me, 'MLS_Swcreate_id' , cond=.false. )
    MLS_SwCREATE = 0
    call check( nf90_def_grp( FileId, swathName, MLS_SwCREATE ), &
      & 'MLS_Swcreate_id' // trim(swathname) )
    call check( nf90_def_grp( MLS_SwCREATE, 'Data Fields', dfid ), &
      & 'MLS_Swcreate_id Data Fields'  // trim(swathname) )
    call check( nf90_def_grp( MLS_SwCREATE, 'Geolocation Fields', glid ), &
      & 'MLS_Swcreate_id Geolocation Fields'  // trim(swathname) )
    call trace_end ( cond=.false. )

  end function MLS_Swcreate_id

  ! -----------------------------------------------  MLS_Swdefdim  -----
  integer function MLS_Swdefdim ( swathid, dimname, dimsize )
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    integer, intent(in)          :: DIMSIZE

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: status

    ! Executable code
    call trace_begin ( me, 'MLS_Swdefdim' , cond=.false. )
    MLS_Swdefdim = 0
    call check( nf90_def_dim( swathid, dimName, dimSize, MLS_Swdefdim ), &
      & 'MLS_Swdefdim'  // trim(dimname) )
    call trace_end ( cond=.false. )

  end function MLS_Swdefdim

  ! -----------------------------------------------  MLS_Swdetach  -----
  ! We don't really have a model for what it means to detach from a NetCDF
  ! swath, which is just a group id
  integer function MLS_Swdetach ( SWATHID )
    integer, intent(in) :: SWATHID      ! ID returned by MLS_Swattach

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: myHdfVersion
    logical :: needsFileName

    ! Executable code
    call trace_begin ( me, 'MLS_Swdetach' , cond=.false. )
    call trace_end ( cond=.false. )
    MLS_Swdetach = nf90_NoErr

  end function MLS_Swdetach

  ! ----------------------------------------------  MLS_Swdiminfo  -----
  integer function MLS_Swdiminfo ( FileID, Dimname )
    integer, intent(in) :: FileID      ! File ID (not swathID)
    character(len=*), intent(in) :: Dimname       ! Dimension name

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'MLS_Swdiminfo' , cond=.false. )
    call check( nf90_inq_dimid( FileId, DimName, MLS_Swdiminfo ), &
      & 'MLS_Swdiminfo' // trim(Dimname) )
    call trace_end ( cond=.false. )

  end function MLS_Swdiminfo

  ! ----------------------------------------------  MLS_dfldsetup  -----
  integer function MLS_dfldsetup ( swathid, fieldname, dimname, maxdimlist, &
    & datatype, chunk_rank, chunksizes, DimIDs, &
    & iFill, rFill, dFill )
    integer, parameter :: RANK = 7
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    character(len=*), intent(in) :: MAXDIMLIST
    integer, intent(in)          :: DATATYPE
    integer, intent(in)          :: CHUNK_RANK
    integer, dimension(:), intent(in) :: chunksizes
    integer, dimension(:), intent(in) :: DimIDs
    integer, optional, intent(in) :: iFill
    real, optional, intent(in) :: rFill
    double precision, optional, intent(in) :: dFill

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: grpId
    integer :: varId

    ! Executable code
    call trace_begin ( me, 'MLS_dfldsetup' , cond=.false. )
    ! Must get grpid for 'Data Fields'
    call check( nf90_inq_ncid( swathId, 'Data Fields', grpId ), &
      & 'MLS_dfldsetup Data Fields group' )
    call check( nf90_def_var(grpid, Fieldname, DATATYPE, dimids, varID), &
      & 'MLS_dfldsetup' // trim(Fieldname) )
    ! Chunked?
    call check( nf90_def_var_chunking(grpId, varid, NF90_CHUNKED, chunksizes), &
      & 'MLS_dfldsetup chunking' // trim(Fieldname) )
    ! Fill values?
    if ( present(IFill) ) then
      call check( nf90_def_var_fill(grpId, varid, 0, IFill), &
      & 'MLS_dfldsetup IFill' )
    elseif ( present(RFill) ) then
      call check( nf90_def_var_fill(grpId, varid, 0, RFill), &
      & 'MLS_dfldsetup RFill' )
    elseif ( present(DFill) ) then
      call check( nf90_def_var_fill(grpId, varid, 0, DFill), &
      & 'MLS_dfldsetup DFill' )
    endif
    MLS_dfldsetup = nf90_noErr
    call trace_end ( cond=.false. )

  end function MLS_dfldsetup

  ! ----------------------------------------------  MLS_gfldsetup  -----
  integer function MLS_gfldsetup ( swathid, fieldname, dimname, maxdimlist, &
    & datatype, chunk_rank, chunksizes, DimIDs, &
    & iFill, rFill, dFill )
    integer, parameter :: RANK = 7
    integer, intent(in) :: SWATHID      ! Swath structure ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: DIMNAME       ! Dimension name
    character(len=*), intent(in) :: MAXDIMLIST
    integer, intent(in)          :: DATATYPE
    integer, intent(in)          :: CHUNK_RANK
    integer, dimension(:), intent(in) :: chunksizes
    integer, dimension(:), intent(in) :: DimIDs
    integer, optional, intent(in) :: iFill
    real, optional, intent(in) :: rFill
    double precision, optional, intent(in) :: dFill

    ! Internal variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: grpId
    integer :: varId

    ! Executable code
    call trace_begin ( me, 'MLS_gfldsetup' , cond=.false. )
    ! Must get grpid for 'Geolocation Fields'
    call check( nf90_inq_ncid( swathId, 'Geolocation Fields', grpId ), &
      & 'MLS_gfldsetup Data Fields group' )
    call check( nf90_def_var(grpid, Fieldname, DATATYPE, dimids, varId), &
      & 'MLS_gfldsetup' // trim(Fieldname) )
    ! Chunked?
    call check( nf90_def_var_chunking(grpId, varid, NF90_CHUNKED, chunksizes), &
      & 'MLS_gfldsetup chunking' // trim(Fieldname) )
    ! Fill values?
    if ( present(IFill) ) then
      call check( nf90_def_var_fill(grpId, varid, 0, IFill), &
      & 'MLS_gfldsetup IFill' )
    elseif ( present(RFill) ) then
      call check( nf90_def_var_fill(grpId, varid, 0, RFill), &
      & 'MLS_gfldsetup RFill' )
    elseif ( present(DFill) ) then
      call check( nf90_def_var_fill(grpId, varid, 0, DFill), &
      & 'MLS_gfldsetup DFill' )
    endif
    MLS_gfldsetup = nf90_noErr
    call trace_end ( cond=.false. )

  end function MLS_gfldsetup

  ! ----------------------------------------  MLS_SwRdattr  -----  
  integer function MLS_SwRdattr ( swathid, &
    & attrname, datatype, count, CharBuffer, &
    & int, intbuffer, realbuffer, realsca, doublesca, doublebuffer )
    integer, intent(in) :: SWATHID      ! Swath ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_CHARTYPE
    integer, intent(in) :: COUNT   ! How many
    ! The following buffers are for reading chars or scalars or arrays
    ! of valid types
    character(len=*), intent(out), optional :: CharBuffer
    integer, intent(out), optional          :: int
    integer, intent(out), optional, dimension(:) :: intbuffer
    double precision, intent(out), optional, dimension(:)    :: doublebuffer
    real, intent(out), optional, dimension(:)    :: realbuffer
    real, intent(out), optional                  :: realsca
    double precision, intent(out), optional      :: doublesca

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'MLS_SwRdattr' , cond=.false. )
    if ( present(CharBuffer) ) then
      call check( nf90_get_att( swathid, nf90_Global, AttrName, CharBuffer ), &
        & 'MLS_SwRdattr chars(:)' // trim(AttrName) )
    elseif ( present(int) ) then
      call check( nf90_get_att( swathid, nf90_Global, AttrName, Int ), &
        & 'MLS_SwRdattr int ' // trim(AttrName) )
    elseif ( present(intbuffer) ) then
      call check( nf90_get_att( swathid, nf90_Global, AttrName, IntBuffer ), &
        & 'MLS_SwRdattr ints(:)' // trim(AttrName) )
    elseif ( present(realbuffer) ) then
      call check( nf90_get_att( swathid, nf90_Global, AttrName, realBuffer ), &
        & 'MLS_SwRdattr reals(:)' // trim(AttrName) )
    elseif ( present(doublebuffer) ) then
      call check( nf90_get_att( swathid, nf90_Global, AttrName, doubleBuffer ), &
        & 'MLS_SwRdattr doubles(:)' // trim(AttrName) )
    elseif ( present(realsca) ) then
      call check( nf90_get_att( swathid, nf90_Global, AttrName, realsca ), &
        & 'MLS_SwRdattr real' // trim(AttrName) )
    elseif ( present(doublesca) ) then
      call check( nf90_get_att( swathid, nf90_Global, AttrName, doublesca ), &
        & 'MLS_SwRdattr double' // trim(AttrName) )
    endif
    MLS_SwRdattr = nf90_NoErr
    call trace_end ( cond=.false. )

  end function MLS_SwRdattr

  ! ----------------------------------------  MLS_SwRdLattr  -----  
  integer function MLS_SwRdLattr ( swathid, &
    & FieldName, Attrname, datatype, CharBuffer, &
    & int, realsca, intbuffer, realbuffer, doublesca, doublebuffer )
    integer, intent(in) :: SWATHID      ! Swath ID
    character(len=*), intent(in) :: FieldName    ! Field name
    character(len=*), intent(in) :: AttrName     ! Attribute name
    integer, intent(in) :: Datatype    ! E.g., MLS_CHARTYPE
    ! integer, intent(in) :: Count   ! How many
    ! The following buffers are for reading chars or scalars or arrays
    ! of valid types
    character(len=*), intent(out), optional      :: CharBuffer
    integer, intent(out), optional               :: int
    real, intent(out), optional                  :: realsca
    integer, intent(out), optional, dimension(:) :: intbuffer
    real, intent(out), optional, dimension(:)    :: realbuffer
    double precision, intent(out), optional, dimension(:)    :: doublebuffer
    double precision, intent(out), optional      :: doublesca

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: GrpId
    integer :: VarId
    integer :: Which_Group

    ! Executable code
    call trace_begin ( me, 'MLS_SwRdattr' , cond=.false. )
    MLS_SwRdLattr = -1
    if ( Is_datafield_in_swath(swathId, fieldname, which_group) ) then
      ! Must get grpid for Groupname
      call check( &
        & nf90_inq_ncid( swathId, trim(SWGroupNames(which_group)), grpId ), &
        & 'MLS_SwRdLattr' )
      ! Try to find varId
      call check( nf90_inq_varid(grpid, Fieldname, varid), &
        & 'MLS_SwRdLattr:' // trim(Fieldname) )
      if ( present(CharBuffer) ) then
        call check( nf90_get_att( grpid, varId, AttrName, CharBuffer ), &
          & 'MLS_SwRdLattr' // trim(AttrName) )
      elseif ( present(IntBuffer) ) then
        call check( nf90_get_att( grpid, varId, AttrName, IntBuffer ), &
          & 'MLS_SwRdLattr' // trim(AttrName) )
      elseif ( present(Int) ) then
        call check( nf90_get_att( grpid, varId, AttrName, Int ), &
          & 'MLS_SwRdLattr' // trim(AttrName) )
      elseif ( present(Realsca) ) then
        call check( nf90_get_att( grpid, varId, AttrName, Realsca ), &
          & 'MLS_SwRdLattr' // trim(AttrName) )
      elseif ( present(RealBuffer) ) then
        call check( nf90_get_att( grpid, varId, AttrName, RealBuffer ), &
          & 'MLS_SwRdLattr' // trim(AttrName) )
      elseif ( present(doubleBuffer) ) then
        call check( nf90_get_att( grpid, varId, AttrName, doubleBuffer ), &
          & 'MLS_SwRdLattr' // trim(AttrName) )
      elseif ( present(Doublesca) ) then
        call check( nf90_get_att( grpid, varId, AttrName, Doublesca ), &
          & 'MLS_SwRdLattr' // trim(AttrName) )
      endif
      MLS_SwRdLattr = nf90_NoErr
    endif
    call trace_end ( cond=.false. )

  end function MLS_SwRdLattr

  function MLS_Swrdfld_Char_1d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 1
    character(len=*), parameter :: ProcName = 'MLS_Swrdfld_Char_1d'
    character(len=*), dimension(:), intent(out) :: values
    include 'MLS_Swrdfld.f9h'

  end function MLS_Swrdfld_Char_1d

  ! ----------------------------------------  MLS_Swrdfld_Double_1d  -----
  function MLS_Swrdfld_Double_1d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 1
    character(len=*), parameter :: ProcName = 'MLS_Swrdfld_Double_1d'
    double precision, dimension(:), intent(out) :: values
    include 'MLS_Swrdfld.f9h'

  end function MLS_Swrdfld_Double_1d

  ! ----------------------------------------  MLS_Swrdfld_Double_2d  -----
  function MLS_Swrdfld_Double_2d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 2
    character(len=*), parameter :: ProcName = 'MLS_Swrdfld_Double_2d'
    double precision, dimension(:,:), intent(out) :: values
    include 'MLS_Swrdfld.f9h'

  end function MLS_Swrdfld_Double_2d

  ! ----------------------------------------  MLS_Swrdfld_Double_3d  -----
  function MLS_Swrdfld_Double_3d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 3
    character(len=*), parameter :: ProcName = 'MLS_Swrdfld_Double_3d'
    double precision, dimension(:,:,:), intent(out) :: values
    include 'MLS_Swrdfld.f9h'

  end function MLS_Swrdfld_Double_3d

  ! ----------------------------------------  MLS_Swrdfld_Integer  -----
  function MLS_Swrdfld_Integer ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 1
    character(len=*), parameter :: ProcName = 'MLS_Swrdfld_Integer'
    integer, dimension(:), intent(out)   :: values
    include 'MLS_Swrdfld.f9h'

  end function MLS_Swrdfld_Integer

  ! ----------------------------------------  MLS_Swrdfld_Real_1d  -----
  function MLS_Swrdfld_Real_1d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 1
    character(len=*), parameter :: ProcName = 'MLS_Swrdfld_Real_1d'
    real, dimension(:), intent(out)      :: values
    include 'MLS_Swrdfld.f9h'

  end function MLS_Swrdfld_Real_1d

  ! ----------------------------------------  MLS_Swrdfld_Real_2d  -----
  function MLS_Swrdfld_Real_2d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 2
    character(len=*), parameter :: ProcName = 'MLS_Swrdfld_Real_2d'
    real, dimension(:,:), intent(out)    :: values
    include 'MLS_Swrdfld.f9h'

  end function MLS_Swrdfld_Real_2d

  ! ----------------------------------------  MLS_Swrdfld_Real_3d  -----
  function MLS_Swrdfld_Real_3d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 3
    character(len=*), parameter :: ProcName = 'MLS_Swrdfld_Real_3d'
    real, dimension(:,:,:), intent(out)  :: values
    include 'MLS_Swrdfld.f9h'

  end function MLS_Swrdfld_Real_3d

  ! -----------------------------------------------  MLS_Swwrattr  -----
  integer function MLS_Swwrattr ( swathid, &
    & attrname, datatype, count, CharBuffer, &
    & int, intbuffer, realsca, realbuffer, doublesca, doublebuffer )
    integer, intent(in) :: SWATHID      ! Swath ID
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_CHARTYPE
    integer, intent(in) :: COUNT   ! How many
    character(len=*), intent(in), optional :: CharBuffer
    integer, intent(in), optional          :: int
    real, intent(in), optional             :: realsca
    integer, intent(in), optional, dimension(:) :: intbuffer
    real, intent(in), optional, dimension(:)    :: realbuffer
    double precision, intent(in), optional, dimension(:)    :: doublebuffer
    double precision, intent(in), optional      :: doublesca

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'MLS_Swwrattr' , cond=.false. )
    if ( present(CharBuffer) ) then
      if ( len_trim(CharBuffer) < 1 ) then
        call check( nf90_put_att( swathid, nf90_Global, AttrName, Blank ), &
        & 'MLS_Swwrattr' // trim(AttrName) )
      else
        call check( nf90_put_att( swathid, nf90_Global, AttrName, CharBuffer ), &
        & 'MLS_Swwrattr' // trim(AttrName) )
      endif
    elseif ( present(IntBuffer) ) then
      call check( nf90_put_att( swathid, nf90_Global, AttrName, IntBuffer ), &
      & 'MLS_Swwrattr' // trim(AttrName) )
    elseif ( present(RealSca) ) then
      call check( nf90_put_att( swathid, nf90_Global, AttrName, RealSca ), &
      & 'MLS_Swwrattr' // trim(AttrName) )
    elseif ( present(RealBuffer) ) then
      call check( nf90_put_att( swathid, nf90_Global, AttrName, RealBuffer ), &
      & 'MLS_Swwrattr' // trim(AttrName) )
    elseif ( present(doubleBuffer) ) then
      call check( nf90_put_att( swathid, nf90_Global, AttrName, doubleBuffer ), &
      & 'MLS_Swwrattr' // trim(AttrName) )
    elseif ( present(Int) ) then
      call check( nf90_put_att( swathid, nf90_Global, AttrName, Int ), &
      & 'MLS_Swwrattr' // trim(AttrName) )
    elseif ( present(DoubleSca) ) then
      call check( nf90_put_att( swathid, nf90_Global, AttrName, DoubleSca ), &
      & 'MLS_Swwrattr' // trim(AttrName) )
    endif
    MLS_Swwrattr = nf90_NoErr
    call trace_end ( cond=.false. )

  end function MLS_Swwrattr

  ! ----------------------------------------------  MLS_Swwrlattr  -----
  ! Write an attrname's value
  integer function MLS_Swwrlattr ( swathid, &
    & fieldname, attrname, datatype, count, CharBuffer, &
    & int, intbuffer, realsca, realbuffer, doublesca, doublebuffer )
    integer, intent(in) :: SWATHID      ! Swath ID
    character(len=*), intent(in) :: FIELDNAME     ! Field name
    character(len=*), intent(in) :: ATTRNAME     ! Attribute name
    integer, intent(in) :: DATATYPE    ! E.g., MLS_CHARTYPE
    integer, intent(in) :: COUNT   ! How many
    character(len=*), intent(in), optional :: CharBuffer
    integer, intent(in), optional          :: int
    real, intent(in), optional             :: realsca
    integer, intent(in), optional, dimension(:) :: intbuffer
    real, intent(in), optional, dimension(:)    :: realbuffer
    double precision, intent(in), optional, dimension(:)    :: doublebuffer
    double precision, intent(in), optional      :: doublesca

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: grpId
    integer :: varId
    integer :: which_group

    ! Executable code
    call trace_begin ( me, 'MLS_Swwrlattr' , cond=.false. )
    MLS_Swwrlattr = -1
    if ( Is_datafield_in_swath(swathId, fieldname, which_group) ) then
      ! Must get grpid for Groupname
      call check( &
        & nf90_inq_ncid( swathId, trim(SWGroupNames(which_group)), grpId ), &
        & 'MLS_dfldsetup Data Fields group' )
      ! Try to find varId
      call check( nf90_inq_varid( grpid, Fieldname, varid), &
        & 'MLS_Swwrlattr:' // trim(Fieldname) )
      if ( present(CharBuffer) ) then
        if ( len_trim(CharBuffer) > 0 ) then
          call check( nf90_put_att( grpid, varid, AttrName, CharBuffer), &
          & 'MLS_Swwrlattr:' // trim(Fieldname) )
        else
          call check( nf90_put_att( grpid, varid, AttrName, Blank), &
          & 'MLS_Swwrlattr:' // trim(Fieldname) )
        endif
      elseif ( present(IntBuffer) ) then
        call check( nf90_put_att( grpid, varid, AttrName, IntBuffer ), &
        & 'MLS_Swwrlattr' // trim(FieldName) )
      elseif ( present(RealSca) ) then
        call check( nf90_put_att( grpid, varid, AttrName, RealSca ), &
        & 'MLS_Swwrlattr' // trim(FieldName) )
      elseif ( present(RealBuffer) ) then
        call check( nf90_put_att( grpid, varid, AttrName, RealBuffer ), &
        & 'MLS_Swwrlattr' // trim(FieldName) )
      elseif ( present(DoubleBuffer) ) then
        call check( nf90_put_att( grpid, varid, AttrName, DoubleBuffer ), &
        & 'MLS_Swwrlattr' // trim(FieldName) )
      elseif ( present(DoubleSca) ) then
        call check( nf90_put_att( grpid, varid, AttrName, DoubleSca ), &
        & 'MLS_Swwrlattr' // trim(FieldName) )
      endif
      MLS_Swwrlattr = nf90_NoErr
    endif
    call trace_end ( cond=.false. )

  end function MLS_Swwrlattr

  ! ----------------------------------------  MLS_Swwrfld_Char_1D  -----
  function MLS_Swwrfld_Char_1D ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 1
    character(len=*), parameter :: ProcName = 'MLS_Swwrfld_Char_1d'
    character(len=*), dimension(:), intent(in) :: values
    include 'MLS_Swwrfld.f9h'

  end function MLS_Swwrfld_Char_1D

  ! ----------------------------------------  MLS_Swwrfld_Double_1d  -----
  function MLS_Swwrfld_Double_1d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 1
    character(len=*), parameter :: ProcName = 'MLS_Swwrfld_Double_1d'
    double precision, dimension(:), intent(in) :: values
    include 'MLS_Swwrfld.f9h'

  end function MLS_Swwrfld_Double_1d

  ! ----------------------------------------  MLS_Swwrfld_Double_2d  -----
  function MLS_Swwrfld_Double_2d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 2
    character(len=*), parameter :: ProcName = 'MLS_Swwrfld_Double_2d'
    double precision, dimension(:,:), intent(in) :: values
    include 'MLS_Swwrfld.f9h'

  end function MLS_Swwrfld_Double_2d

  ! ----------------------------------------  MLS_Swwrfld_Double_3d  -----
  function MLS_Swwrfld_Double_3d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 3
    character(len=*), parameter :: ProcName = 'MLS_Swwrfld_Double_3d'
    double precision, dimension(:,:,:), intent(in) :: values
    include 'MLS_Swwrfld.f9h'

  end function MLS_Swwrfld_Double_3d

  ! ----------------------------------------  MLS_Swwrfld_Integer  -----
  function MLS_Swwrfld_Integer ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 1
    character(len=*), parameter :: ProcName = 'MLS_Swwrfld_Integer'
    integer, dimension(:), intent(in)   :: values
    include 'MLS_Swwrfld.f9h'

  end function MLS_Swwrfld_Integer

  ! ----------------------------------------  MLS_Swwrfld_Real_1d  -----
  function MLS_Swwrfld_Real_1d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 1
    character(len=*), parameter :: ProcName = 'MLS_Swwrfld_Real_1d'
    real, dimension(:), intent(in)      :: values
    include 'MLS_Swwrfld.f9h'

  end function MLS_Swwrfld_Real_1d

  ! ----------------------------------------  MLS_Swwrfld_Real_2d  -----
  function MLS_Swwrfld_Real_2d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 2
    character(len=*), parameter :: ProcName = 'MLS_Swwrfld_Real_2d'
    real, dimension(:,:), intent(in)    :: values
    include 'MLS_Swwrfld.f9h'

  end function MLS_Swwrfld_Real_2d

  ! ----------------------------------------  MLS_Swwrfld_Real_3d  -----
  function MLS_Swwrfld_Real_3d ( swathid, fieldname, &
    & start, stride, count, values ) result ( returnStatus )
    integer, parameter :: RANK = 3
    character(len=*), parameter :: ProcName = 'MLS_Swwrfld_Real_d'
    real, dimension(:,:,:), intent(in)  :: values
    include 'MLS_Swwrfld.f9h'

  end function MLS_Swwrfld_Real_3d

  ! --------------------------------------  MLS_Swath_In_File_sca  -----
  logical function MLS_Swath_In_File_sca( filename, swath, error )
    ! Returns .true. if swath found in file, .false. otherwise
    use HDF5, only: Size_t
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: swath
    integer, optional, intent(out) :: error

    ! Internal variables
    logical, parameter               :: Deebug = .false.
    integer :: Me = -1               !   String index for trace cacheing
    integer                          :: FileId
    integer                          :: grpId

    ! Begin execution
    ! print*,'Scalar version'
    call trace_begin ( me, 'MLS_Swath_In_File_sca' , cond=.false. )
    ! Now let's try to read the file we just created
    call check( nf90_open(Filename, NF90_NoWrite, FileId) , &
      & 'MLS_Swath_In_File_sca ' // trim(FileName) )
    ! Any group with that name?
    call check( nf90_inq_ncid( FileId, swath, grpId) , &
      & 'MLS_Swath_In_File_sca ' // trim(swath), &
      & FailureOK=.true., silent=.true. )
    MLS_Swath_In_File_sca = (NCError == nf90_noerr)
    call trace_end ( cond=.false. )
  end function MLS_Swath_In_File_sca

  ! --------------------------------------  MLS_Swath_In_File_fid  -----
  logical function MLS_Swath_In_File_fid( fileId, swath, error )
    ! Returns .true. if swath found in file, .false. otherwise
    use HDF5, only: Size_t
    integer, intent(in)            :: FileId
    character(len=*), intent(in)   :: swath
    integer, optional, intent(out) :: error

    ! Internal variables
    logical, parameter               :: Deebug = .false.
    integer :: Me = -1               !   String index for trace cacheing
    integer                          :: grpId

    ! Begin execution
    ! print*,'Scalar version'
    call trace_begin ( me, 'MLS_Swath_In_File_fid' , cond=.false. )
    ! Any group with that name?
    call check( nf90_inq_ncid( FileId, swath, grpId) , &
      & 'MLS_Swath_In_File_fid ' // trim(swath), &
      & FailureOK=.true., silent=.true. )
    MLS_Swath_In_File_fid = (NCError == nf90_noerr)
    call trace_end ( cond=.false. )
  end function MLS_Swath_In_File_fid

  ! --------------------------------------  MLS_Swath_In_File_arr  -----
  logical function MLS_Swath_In_File_arr(filename, swaths, &
    & which, error )
    use HDF5, only: Size_t
    ! Array version of the above
    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: swaths
    logical, dimension(:), intent(out) :: which
    integer, optional, intent(out) :: error

    ! Internal variables
    integer                          :: I
    integer :: Me = -1               !   String index for trace cacheing
    integer                          :: Nswaths

    ! Begin execution
    call trace_begin ( me, 'MLS_Swath_In_File_arr' , cond=.false. )
    nswaths = 0
    MLS_Swath_In_File_arr = .false.
    which = .false.
    if ( size(swaths) > size(which) ) &
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
          & 'array to small to hold values in MLS_Swath_In_File_arr' )
    do i=1, size(swaths)
      which(i) = &
      & MLS_Swath_In_File_sca( filename, swaths(i), error )
    enddo
    MLS_Swath_In_File_arr = any( which )
    call trace_end ( cond=.false. )
  end function MLS_Swath_In_File_arr

! ======================= Private Procedures ===========================

  ! --------------------------------------  check  -----
  ! Check that the NetCDF operation was successful; i.e., nf90_noerr
  ! If not, print mesg (unless silent)
  ! and halt with error (unless FailureOK)
  subroutine check( status, mesg, FailureOK, silent )
    integer, intent ( in )                          :: status ! how was it?
    character(len=*), intent ( in )                 :: mesg
    logical, optional, intent ( in )                :: FailureOK ! Don't halt
    logical, optional, intent ( in )                :: silent ! Don't print
    !
    logical                                         :: MyFailOK
    logical                                         :: MySilent
    ! Executable
    NCError = status ! So the caller knows in case we don't exit
    MyFailOK = .false.
    if ( present(FailureOK) ) MyFailOK = FailureOK
    MySilent = .false.
    if ( present(Silent) ) MySilent = Silent
    if( status /= nf90_noerr ) then 
      if ( .not. MySilent ) &
        & call output( trim(nf90_strerror(status)), advance='yes' )
      if ( MyFailOK ) return
      call MLSMessage ( MLSMSG_Error, moduleName,  &
          & trim(mesg) )
    end if
  end subroutine check  
  
  ! --------------------------------------  is_datafield_in_group  -----
  logical function is_datafield_in_group( swathid, field, groupname )
    ! Returns .true. if datafield found in group, .false. otherwise
    integer, intent(in) :: swathid
    character(len=*), intent(in) :: field
    character(len=*), intent(in) :: groupName
    ! Internal variables
    integer                      :: grpId
    integer                      :: varId
    ! Begin execution
    ! Must get grpid for Groupname
    call check( nf90_inq_ncid( swathId, trim(GroupName), grpId ), &
      & 'MLS_dfldsetup Data Fields group' )
    ! Try to find varId
    call check( nf90_inq_varid(grpid, Field, varid), &
      & 'is_datafield_in_group:' // trim(Field) // ':' // trim(Groupname), &
      & FailureOK=.true., silent=.true. )
    is_datafield_in_group = (NCError == nf90_noerr)
  end function is_datafield_in_group

  ! --------------------------------------  is_datafield_in_swath  -----
  logical function is_datafield_in_swath( swathid, field, which_group )
    ! Returns .true. if datafield found in swath, .false. otherwise
    ! Also say which group of vars the field is part of
    ! 'Data Fields' or 'Geolocation Fields'
    integer, intent(in)          :: swathid
    character(len=*), intent(in) :: field
    integer, intent(out)         :: which_group ! 1 if Data, 2 if Geolocations
    ! Internal variables
    ! Begin execution
    which_group = 1
    is_datafield_in_swath = is_datafield_in_group ( swathId, field, &
      & 'Data Fields' )
    if ( is_datafield_in_swath ) return
    which_group = 2
    is_datafield_in_swath = is_datafield_in_group ( swathId, field, &
      & 'Geolocation Fields' )
    if ( is_datafield_in_swath ) return
    which_group = 0
  end function is_datafield_in_swath

  ! --------------------------------------  is_swath_in_file  -----
  logical function is_swath_in_file( Fileid, swath, swid )
    ! Returns .true. if swath found in file, .false. otherwise
    ! Also say what the swath id is
    integer, intent(in)          :: Fileid
    character(len=*), intent(in) :: swath
    integer, intent(out)         :: swid
    ! Internal variables
    ! Begin execution
    swid = 0
    call check ( nf90_inq_ncid ( FileId, swath, swid), &
      & 'is_swath_in_file:' // trim(swath), &
      & FailureOK=.true., silent=.true. )
    is_swath_in_file = (NCError == nf90_noerr)
  end function is_swath_in_file

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

end module MLSNetCDF4

! $Log$
! Revision 1.2  2020/03/19 22:33:29  pwagner
! Repaired many bugs in writing global attrs
!
! Revision 1.1  2020/03/06 00:24:19  pwagner
! First commit
!
