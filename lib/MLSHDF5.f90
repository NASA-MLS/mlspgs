! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MLSHDF5

  ! This module contains MLS specific routines to do lowish level common HDF5
  ! tasks.

  ! Usage note: due to idiosyncrasy (euphemism for bug) in Lahey f95 compiler,
  ! don't put a "USE MLSHDF5" statement at module level in any modules that
  ! will use this. Instead bury each "USE MLSHDF5" down inside the
  ! procedures that need them. Otherwise, Lahey will take nearly forever
  ! to finish compiling the highest modules.

  use DUMP_0, only: DUMP, DUMP_NAME_V_PAIRS
  use MLSDataInfo, only: MLSDataInfo_T, Query_MLSData
  ! To switch to/from hdfeos5.1.6(+) uncomment next line
  use H5LIB, ONLY: h5open_f, h5close_f
  ! Lets break down our use, parameters first
  use HDF5, only: H5F_ACC_RDONLY_F, &
    & H5P_DATASET_CREATE_F, &
    & H5SIS_SIMPLE_F, & ! H5SOFFSET_SIMPLE_F, &
    & H5S_SCALAR_F, H5S_SELECT_SET_F, H5S_UNLIMITED_F, &
    & H5T_IEEE_F32LE, H5T_IEEE_F64LE, &
    & H5T_NATIVE_DOUBLE, H5T_NATIVE_REAL, H5T_STD_I32LE, &
    & H5T_NATIVE_CHARACTER, H5T_NATIVE_INTEGER, &
    & HID_T, HSIZE_T ! , HSSIZE_T
  ! Now routines
  use HDF5, only: H5ACLOSE_F, H5ACREATE_F, H5AGET_TYPE_F, H5AOPEN_NAME_F, &
    & H5AREAD_F, H5AWRITE_F, H5ADELETE_F, &
    & H5DCREATE_F, H5DEXTEND_F, H5DGET_SPACE_F, H5DGET_TYPE_F, H5DOPEN_F, &
    & H5DREAD_F, H5DWRITE_F, H5DCLOSE_F, H5DGET_CREATE_PLIST_F, &
    & H5ESET_AUTO_F, &
    & H5FOPEN_F, H5FCLOSE_F, &
    & H5GOPEN_F, H5GCLOSE_F, &
    & H5PCREATE_F, H5PSET_CHUNK_F, H5PSET_FILL_VALUE_F, &
    & H5PGET_CHUNK_F, H5PGET_FILL_VALUE_F, &
    & H5SCLOSE_F, &
    & H5SCREATE_F, H5SCREATE_SIMPLE_F, H5SGET_SIMPLE_EXTENT_NDIMS_F, &
    & H5SGET_SIMPLE_EXTENT_DIMS_F, H5SSELECT_HYPERSLAB_F, &
    & H5TCLOSE_F, H5TCOPY_F, H5TEQUAL_F, H5TGET_SIZE_F, H5TSET_SIZE_F
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_WARNING
  use MLSStringLists, only: catLists

  implicit NONE
  private

  public :: CpHDF5Attribute, CpHDF5GlAttribute, &
    & GetHDF5Attribute, GetHDF5AttributePtr, GetHDF5AttrDims, &
    & GetAllHDF5DSNames, GetHDF5DSRank, GetHDF5DSDims, GetHDF5DSQType, &
    & IsHDF5AttributeInFile, IsHDF5AttributePresent, IsHDF5DSInFile, &
    & IsHDF5DSPresent, IsHDF5GroupPresent, &
    & LoadFromHDF5DS, LoadPtrFromHDF5DS, MakeHDF5Attribute, &
    & MLS_H5Open, MLS_H5Close, &
    & ReadLitIndexFromHDF5Attr, ReadStringIndexFromHDF5Attr, SaveAsHDF5DS, &
    & WriteLitIndexAsHDF5Attribute, WriteStringIndexAsHDF5Attribute

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
  !---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! CpHDF5Attribute      Copies an attribute
! CpHDF5GlAttribute    Copies a global attribute
! GetAllHDF5DSNames    Retrieves names of all DS in file (under group name)
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
! mls_h5close          Closes interface to hdf5; call once at end of run
! mls_h5open           Opens interface to hdf5; call once at start of run
! SaveAsHDF5DS         Turns an array into a dataset
! === (end of toc) ===

! === (start of api) ===
! CpHDF5Attribute (int fromitemID, int toitemID, char name,
!    [log skip_if_already_there])
! CpHDF5GlAttribute (char fromFile, char toFile, char name,
!    [log skip_if_already_there])
! GetAllHDF5DSNames (file, char gname, char DSNames)
!     file can be one of:
!    {char* filename, int fileID}
! GetHDF5Attribute (int itemID, char name, value)
! GetHDF5AttributePtr (int itemID, char name, *value)
! GetHDF5DSDims (int FileID, char name, hsize_t dims(:), [hsize_t maxdims(:)])
! GetHDF5DSRank (int FileID, char name, int rank)
! GetHDF5DSQType (int FileID, char name, char QType)
! log IsHDF5AttributeInFile (char filename, char DSname, char name)
! log IsHDF5DSInFile (char filename, char name)
! log IsHDF5AttributePresent (int setid, char name)
! log IsHDF5AttributePresent (int fileid, char DSname, char name)
! log IsHDF5DSPresent (int locID, char name)
! log IsHDF5GroupPresent (int locID, char name)
! LoadFromHDF5DS (int locID, char name, value,
!       [int start(:), int count(:), [int stride(:), int block(:)] ] )
! LoadPtrFromHDF5DS (int locID, char name, *value [, lowBound] )
!       lowBound only for rank-1 arrays
! MakeHDF5Attribute (int itemID, char name, value,
!       [log skip_if_already_there])
! SaveAsHDF5DS (int locID, char name, value)
!     value can be one of:
!    {char* value, int value, real value, double precision value,
!     int value(:),
!     real value(:), real value(:,:), real value(:,:,:), real value(:,:,:,:),
!     double precision value(:), double precision value(:,:),
!     double precision value(:,:,:)}
! === (end of api) ===
  interface GetAllHDF5DSNames
    module procedure GetAllHDF5DSNames_fileID, GetAllHDF5DSNames_filename
  end interface

  interface MakeHDF5Attribute
    module procedure MakeHDF5Attribute_dbl, MakeHDF5Attribute_sngl, &
      & MakeHDF5Attribute_int, MakeHDF5Attribute_logical, &
      & MakeHDF5Attribute_logicalarr1, &
      & MakeHDF5Attribute_string, MakeHDF5Attribute_snglarr1, &
      & MakeHDF5Attribute_dblarr1, MakeHDF5Attribute_stringarr1, &
      & MakeHDF5Attribute_intarr1, MakeHDF5AttributeDSN_int, &
      & MakeHDF5AttributeDSN_string, MakeHDF5AttributeDSN_snglarr1, &
      & MakeHDF5AttributeDSN_st_arr1, MakeHDF5AttributeDSN_dblarr1
  end interface

  interface GetHDF5Attribute
    module procedure GetHDF5Attribute_int, GetHDF5Attribute_logical, &
      & GetHDF5Attribute_logicalarr1, &
      & GetHDF5Attribute_string, GetHDF5Attribute_sngl, GetHDF5Attribute_dbl, &
      & GetHDF5Attribute_snglarr1, GetHDF5Attribute_intarr1, &
      & GetHDF5Attribute_dblarr1, GetHDF5Attribute_stringarr1
  end interface

  interface GetHDF5AttributePtr
    module procedure GetHDF5AttributePtr_snglarr1, GetHDF5AttributePtr_intarr1, &
      & GetHDF5AttributePtr_dblarr1, GetHDF5AttributePtr_stringarr1, &
      & GetHDF5AttributePtr_logicalarr1
  end interface

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

  interface IsHDF5AttributePresent
    module procedure IsHDF5AttributePresent_in_fID, &
      & IsHDF5AttributePresent_in_DSID, IsHDF5AttributePresent_in_grp
  end interface

  interface SaveAsHDF5DS
    module procedure &
      & SaveAsHDF5DS_intarr1, SaveAsHDF5DS_intarr2, SaveAsHDF5DS_intarr3, &
      & SaveAsHDF5DS_logarr1, &
      & SaveAsHDF5DS_dblarr1, SaveAsHDF5DS_dblarr2, SaveAsHDF5DS_dblarr3, &
      & SaveAsHDF5DS_snglarr1, SaveAsHDF5DS_snglarr2, SaveAsHDF5DS_snglarr3, &
      & SaveAsHDF5DS_snglarr4, &
      & SaveAsHDF5DS_charsclr, SaveAsHDF5DS_chararr1, SaveAsHDF5DS_chararr2
  end interface

  interface LoadFromHDF5DS
    module procedure LoadFromHDF5DS_intarr1, LoadFromHDF5DS_intarr2, &
      & LoadFromHDF5DS_logarr1, &
      & LoadFromHDF5DS_dblarr1, LoadFromHDF5DS_dblarr2, LoadFromHDF5DS_dblarr3, &
      & LoadFromHDF5DS_snglarr1, LoadFromHDF5DS_snglarr2, &
      & LoadFromHDF5DS_snglarr3, LoadFromHDF5DS_snglarr4, &
      & LoadFromHDF5DS_chararr1, LoadFromHDF5DS_chararr2
  end interface

  interface LoadPtrFromHDF5DS
    module procedure LoadPtrFromHDF5DS_chararr1, LoadPtrFromHDF5DS_chararr2, &
      & LoadPtrFromHDF5DS_intarr1, LoadPtrFromHDF5DS_intarr2, &
      & LoadPtrFromHDF5DS_logarr1, LoadPtrFromHDF5DS_dblarr1, &
      & LoadPtrFromHDF5DS_dblarr2, LoadPtrFromHDF5DS_dblarr3, &
      & LoadPtrFromHDF5DS_snglarr1, LoadPtrFromHDF5DS_snglarr2, &
      & LoadPtrFromHDF5DS_snglarr3, LoadPtrFromHDF5DS_snglarr4
  end interface

  ! Local parameters
  integer, dimension(7) :: ones = (/1,1,1,1,1,1,1/)
  logical, parameter    :: DEEBUG = .false.
  integer, save :: cantGetDataspaceDims = 0
  integer, parameter :: MAXNUMWARNS = 40
  integer, parameter :: MAXNDSNAMES = 1000   ! max number of DS names in a file
  ! Local variables
  integer(hid_t) :: cparms

contains ! ======================= Public Procedures =========================

  ! ------------------------------------------------  MLS_h5close  -----
  subroutine MLS_h5close ( error )
    ! Arguments
    integer, intent(out) :: error          ! Trouble if /= 0
    error = 0
    call h5close_f ( error )
  end subroutine MLS_h5close

  ! -------------------------------------------------  MLS_h5open  -----
  subroutine MLS_h5open ( error )
    ! Arguments
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
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    logical :: my_skip
    logical :: is_present
    character (len=2000) :: value1
    ! character (len=*), intent(in) :: VALUE ! Value of attribute

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip = skip_if_already_there
    is_present = IsHDF5AttributePresent_in_DSID(fromitemID, name)
    if ( .not. is_present ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to cp: attribute not found ' // trim(name) )
    is_present = IsHDF5AttributePresent_in_DSID(toitemID, name)
    if ( my_skip .and. is_present ) return
    call GetHDF5Attribute ( fromitemID, name, value1 )
    call MakeHDF5Attribute ( toitemID, name, trim(value1), &
      & skip_if_already_there )
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
    integer :: STATUS                   ! Flag
    integer :: tofileID
    integer :: togrpID

    ! Executable code
    call h5fopen_f ( trim(fromFilename), H5F_ACC_RDONLY_F, fromFileID, status )
    call h5fopen_f ( trim(toFilename), H5F_ACC_RDONLY_F, toFileID, status )
    call h5gopen_f ( fromFileID, '/', fromgrpid, status )
    call h5gopen_f ( toFileID, '/', togrpid, status )
    call CpHDF5Attribute_string ( fromgrpID, togrpID, name, &
      & skip_if_already_there )
    call h5gclose_f ( fromgrpid, status )
    call h5gclose_f ( togrpid, status )
    call h5fclose_f ( fromFileID, status )
    call h5fclose_f ( toFileID, status )
  end subroutine CpHDF5GlAttribute_string

  ! -----------------------------------  GetAllHDF5DSNames_fileID  -----
  subroutine GetAllHDF5DSNames_fileID ( FileID, gname, DSNames )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    integer, intent(in) :: FILEID          ! fileID
    character (len=*), intent(in) :: GNAME ! Name of group; e.g. '/'
    character (len=*), intent(out) :: DSNames ! Names of DS in file (,-separated)

    ! Local variables
    integer :: i                        ! loop counter
    integer :: STATUS                   ! Flag
    type(MLSDataInfo_T) :: dataset_info

    ! Executable code
    ! Initializing values returned if there was trouble
    DSNames = ''
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting DSNames' )
! The structure, dataset_info, is initialized below.
!
    call allocate_test ( dataset_info%name, MAXNDSNAMES, 'dataset_info%name', &
      & moduleName )
    dataset_info%name = ''
    dataset_info%number_of_entries = 0
    call Query_MLSData ( fileid, trim(gname), dataset_info )
    if ( dataset_info%number_of_entries < 0 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Unable to get DSNames' )
    else if ( dataset_info%number_of_entries == 1 ) then
      DSNames = dataset_info%name(1)
    else if ( dataset_info%number_of_entries > 1 ) then
      ! DSNames = dataset_info%name(1)
      do i = 1, dataset_info%number_of_entries
        DSNames = catLists(trim(DSNames), dataset_info%name(i))
      end do
    end if
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting DSNames' )
    call deallocate_test ( dataset_info%name, 'dataset_info%name', moduleName )
  end subroutine GetAllHDF5DSNames_fileID

  ! ---------------------------------  GetAllHDF5DSNames_filename  -----
  subroutine GetAllHDF5DSNames_filename ( FileName, gname, DSNames )
    character (len=*), intent(in) :: FILENAME       ! file name
    character (len=*), intent(in) :: GNAME ! Name of group; e.g. '/'
    character (len=*), intent(out) :: DSNames ! Names of DS in file (,-separated)

    ! Local variables
    integer :: fileID
    integer :: STATUS                   ! Flag

    ! Executable code
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting DSNames ' // trim(filename) )
    call h5fopen_f ( trim(filename), H5F_ACC_RDONLY_F, fileID, status )
    if ( status == 0 ) then
      call h5eSet_auto_f ( 1, status )
      call GetAllHDF5DSNames_fileID ( fileID, gname, DSNames )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open file for getting DSNames ' // trim(filename) )
    end if
    call h5eSet_auto_f ( 0, status )
    call h5fclose_f ( fileID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close file after getting DSNames ' // trim(filename) )
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting DSNames ' // trim(filename) )
  end subroutine GetAllHDF5DSNames_filename

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
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    ! (Maybe) create the attribute
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_DOUBLE, dsID, attrID &
      &, skip_if_already_there ) ) then
      ! Write it
      call h5aWrite_f ( attrID, H5T_NATIVE_DOUBLE, value, ones, status )
      call finishMakeAttrib ( name, status, attrID, dsID )
    end if
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
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    ! (Maybe) create the attribute
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_REAL, dsID, attrID, &
      & skip_if_already_there ) ) then
      ! Write
      call h5aWrite_f ( attrID, H5T_NATIVE_REAL, value, ones, status )
      call finishMakeAttrib ( name, status, attrID, dsID )
    end if
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
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    ! (Maybe) create the attribute
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_INTEGER, dsID, attrID, &
      & skip_if_already_there ) ) then
      ! Write
      call h5aWrite_f ( attrID, H5T_NATIVE_INTEGER, value, ones, status )
      call finishMakeAttrib ( name, status, attrID, dsID )
    end if
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

    ! Executable code
    iValue = 0
    if ( value ) iValue = 1
    call MakeHDF5Attribute ( itemID, name, iValue, skip_if_already_there )
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

    ! Executable code
    where ( value )
      iValue = 1
    elsewhere
      iValue = 0
    end where
    call MakeHDF5Attribute ( itemID, name, iValue, skip_if_already_there )
  end subroutine MakeHDF5Attribute_logicalarr1

  ! -----------------------------------  MakeHDF5Attribute_string  -----
  subroutine MakeHDF5Attribute_string ( itemID, name, value , &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    character (len=*), intent(in) :: VALUE ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: DSID                     ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGTYPE               ! Type for string
    logical :: my_skip
    logical :: is_present

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    is_present = IsHDF5AttributePresent_in_DSID(itemID, name)
    if ( my_skip .and. is_present ) return
    ! Setup
    ! Create a data type for this string
    call h5tcopy_f ( H5T_NATIVE_CHARACTER, stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype ' // trim(name) )
    call h5tset_size_f(stringtype, max(len_trim(value), 1), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype ' // trim(name) )
    ! Create dataspace and attribute
    !call h5sCreate_F ( h5s_scalar_f, dsID, status )
    !if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
    !  & 'Unable to create dataspace for attribute ' // trim(name) )
    ! print *, 'itemID: ', itemID
    ! print *, 'stringtype: ', stringtype
    ! print *, 'dsID: ', dsID
    ! print *, 'name: ', trim(name)
    if ( is_present ) then
      call h5adelete_f(itemID, trim(name), status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to delete ' )
    end if
    call h5sCreate_F ( h5s_scalar_f, dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for attribute ' // trim(name) )
    call h5aCreate_f ( itemID, trim(name), stringtype, dsID, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create attribute ' // trim(name) )
    ! print *, 'attrID: ', attrID
    ! print *, 'status: ', status
    ! Write
    call h5aWrite_f ( attrID, stringtype, value, ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write attribute ' // trim(name) )
    ! Finish off
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
    if ( .not. is_present ) then
      call h5sClose_f ( dsID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close attribute dataspace ' // trim(name) )
    end if
    call h5tClose_f ( stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close stringtype ' // trim(name) )
  end subroutine MakeHDF5Attribute_string

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
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGTYPE               ! Type for string
    logical :: my_skip
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip = skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_DSID(itemID, name) ) return
    end if
    ! Setup
    shp = shape(value)
    ! Create a data type for this string
    call h5tcopy_f ( H5T_NATIVE_CHARACTER, stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype for array' // trim(name) )
    call h5tset_size_f ( stringtype, len(value(1)), status )
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
    integer :: spaceID                  ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    ! (Maybe) create the attribute
    shp = shape(value)
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_REAL, spaceID, &
      & attrID, skip_if_already_there, shp ) ) then
      ! Write
      call h5aWrite_f ( attrID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:6) /), hID_T ), status )
      call finishMakeAttrib ( name, status, attrID, spaceID )
    end if
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
    integer :: spaceID                  ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    ! (Maybe) create the attribute
    shp = shape(value)
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_DOUBLE, spaceID, &
      & attrID, skip_if_already_there, shp ) ) then
      ! Write
      call h5aWrite_f ( attrID, H5T_NATIVE_DOUBLE, value, &
        & int ( (/ shp, ones(1:6) /), hID_T ), status )
      call finishMakeAttrib ( name, status, attrID, spaceID )
    end if
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
    integer :: spaceID                  ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    ! (Maybe) create the attribute
    shp = shape(value)
    if ( startMakeAttrib ( itemId, name, H5T_NATIVE_INTEGER, spaceID, &
      & attrID, skip_if_already_there, shp ) ) then
      ! Write
      call h5aWrite_f ( attrID, H5T_NATIVE_INTEGER, value, &
        & int ( (/ shp, ones(1:6) /), hID_T ), status )
      call finishMakeAttrib ( name, status, attrID, spaceID )
    end if
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
    integer :: STATUS                   ! Flag from HDF5
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) return
    end if
    dataID = name_to_dataID( fileID, dataName )
    call MakeHDF5Attribute_int ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )

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
    integer :: STATUS                   ! Flag from HDF5
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) return
    end if
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute_logical ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )

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
    integer :: STATUS                   ! Flag from HDF5
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) return
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
    integer :: STATUS                   ! Flag from HDF5
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) return
    end if
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute ( dataID, attrName, value )
    call h5dclose_f(dataID, status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )

  end subroutine MakeHDF5AttributeDSN_st_arr1

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
    integer :: STATUS                   ! Flag from HDF5
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) return
    end if
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )

  end subroutine MakeHDF5AttributeDSN_snglarr1

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
    integer :: STATUS                   ! Flag from HDF5
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_fID( fileID, dataName, attrName ) ) return
    end if
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute ( dataID, attrName, value )
    call h5dclose_f ( dataID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data ' // trim(dataName) )

  end subroutine MakeHDF5AttributeDSN_dblarr1

  ! ---------------------------------------  GetHDF5Attribute_int  -----
  subroutine GetHDF5Attribute_int ( itemID, name, value )
    integer, intent(in) :: ITEMID         ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(out) :: VALUE         ! Result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read an array into our one value.
    call h5aread_f ( attrID, H5T_NATIVE_INTEGER, value, ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
  end subroutine GetHDF5Attribute_int

  ! -----------------------------------  GetHDF5Attribute_intarr1  -----
  subroutine GetHDF5Attribute_intarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID         ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(out) :: VALUE(:)      ! Result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    shp = shape(value)
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read an array into our one value.
    call h5aread_f ( attrID, H5T_NATIVE_INTEGER, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
  end subroutine GetHDF5Attribute_intarr1

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
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=hSize_t), dimension(1) :: MaxShp, SHP  ! Shape

    ! Executable code
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    call GetHDF5AttrDims ( itemID, trim(name), shp, maxShp )
    call allocate_test ( value, int(lb-1+maxShp(1)), 'Value', moduleName )
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read an array into our one value.
    call h5aread_f ( attrID, H5T_NATIVE_INTEGER, value, &
      & (/ int(shp,hid_t), int(ones(1:6),hid_t) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
  end subroutine GetHDF5AttributePtr_intarr1

  ! ------------------------------------  GetHDF5Attribute_string  -----
  subroutine GetHDF5Attribute_string ( itemID, name, value )
    integer, intent(in) :: ITEMID           ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME   ! Name of attribute
    character (len=*), intent(out) :: VALUE ! Result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGTYPE               ! String type
    integer :: STRINGSIZE               ! String size

    ! Executable code
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    call h5aGet_type_f ( attrID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for attribute ' // trim(name) )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for attribute ' // trim(name) )
    if ( stringSize > len(value) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for attribute ' // trim(name) )
    ! Now actually read the data!
    value = ''
    call h5aread_f ( attrID, stringType, value(1:stringSize), ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
    call h5tClose_f ( stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close string type for attribute ' // trim(name) )
  end subroutine GetHDF5Attribute_string

  ! --------------------------------  GetHDF5Attribute_stringarr1  -----
  subroutine GetHDF5Attribute_stringarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    character (len=*), intent(out) :: VALUE(:) ! Result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGTYPE               ! String type
    integer :: STRINGSIZE               ! String size
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    shp = shape(value)
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
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
    call h5tClose_f ( stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close string type for attribute ' // trim(name) )
  end subroutine GetHDF5Attribute_stringarr1

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
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGTYPE               ! String type
    integer :: STRINGSIZE               ! String size
    integer(kind=hSize_t), dimension(1) :: MaxShp, SHP ! Shape

    ! Executable code
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    call GetHDF5AttrDims ( itemID, trim(name), shp, maxShp )
    call allocate_test ( value, int(lb-1+maxShp(1)), 'Value', moduleName )
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
      & (/ int(shp,hid_t), int(ones(1:6),hid_t) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute ' // trim(name) )
    call h5tClose_f ( stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close string type for attribute ' // trim(name) )
  end subroutine GetHDF5AttributePtr_stringarr1

  ! -----------------------------------  GetHDF5Attribute_logical  -----
  subroutine GetHDF5Attribute_logical ( itemID, name, value )
    integer, intent(in) :: ITEMID         ! Group etc. to get attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(out) :: VALUE         ! Value of attribute

    ! Local variables
    integer :: IVALUE                     ! Value as integer

    ! Executable code
    call GetHDF5Attribute ( itemID, name, iValue )
    value = ( iValue == 1 )
  end subroutine GetHDF5Attribute_logical

  ! -------------------------------  GetHDF5Attribute_logicalarr1  -----
  subroutine GetHDF5Attribute_logicalarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID         ! Group etc. to get attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(out) :: VALUE(:)      ! Value of attribute

    ! Local variables
    integer :: IVALUE(size(VALUE,1))      ! Value as integer

    ! Executable code
    call GetHDF5Attribute ( itemID, name, iValue )
    value = ( iValue == 1 )
  end subroutine GetHDF5Attribute_logicalarr1

  ! ----------------------------  GetHDF5AttributePtr_logicalarr1  -----
  subroutine GetHDF5AttributePtr_logicalarr1 ( itemID, name, value, LowBound )
    use Allocate_Deallocate, only: Deallocate_Test
    integer, intent(in) :: ITEMID         ! Group etc. to get attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, pointer :: VALUE(:)          ! Value of attribute
    integer, intent(in), optional :: LowBound ! of allocated value, else 1

    ! Local variables
    integer, pointer :: IVALUE(:)         ! Value as integer

    ! Executable code
    nullify ( ivalue )
    call GetHDF5AttributePtr ( itemID, name, iValue, lowBound )
    value = ( iValue == 1 )
    call deallocate_test ( ivalue, 'IValue', moduleName )
  end subroutine GetHDF5AttributePtr_logicalarr1

  ! ----------------------------------  GetHDF5Attribute_snglarr1  -----
  subroutine GetHDF5Attribute_snglarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    real, intent(out) :: VALUE(:)       ! The attribute array result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    shp = shape(value)
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_REAL, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read 1d attribute array ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close 1d attribute array  ' // trim(name) )
  end subroutine GetHDF5Attribute_snglarr1

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
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=hSize_t), dimension(1) :: MaxShp, SHP  ! Shape

    ! Executable code
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    call GetHDF5AttrDims ( itemID, trim(name), shp, maxShp )
    call allocate_test ( value, int(lb-1+maxShp(1)), 'Value', moduleName )
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_REAL, value, &
      & (/ int(shp,hid_t), int(ones(1:6),hid_t) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read 1d attribute array ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close 1d attribute array  ' // trim(name) )
  end subroutine GetHDF5AttributePtr_snglarr1

  ! --------------------------------------  GetHDF5Attribute_sngl  -----
  subroutine GetHDF5Attribute_sngl ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    real, intent(out) :: VALUE          ! The attribute result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_REAL, value, &
      & ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read sngl attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close sngl attribute  ' // trim(name) )
  end subroutine GetHDF5Attribute_sngl

  ! ---------------------------------------  GetHDF5Attribute_dbl  -----
  subroutine GetHDF5Attribute_dbl ( itemID, name, value )
    integer, intent(in) :: ITEMID          ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME  ! Name of attribute
    double precision, intent(out) :: VALUE ! The attribute result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_DOUBLE, value, ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dble attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dble attribute  ' // trim(name) )
  end subroutine GetHDF5Attribute_dbl

  ! -----------------------------------  GetHDF5Attribute_dblarr1  -----
  subroutine GetHDF5Attribute_dblarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME     ! Name of attribute
    double precision, intent(out) :: VALUE(:) ! The attribute result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    shp = shape(value)
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_DOUBLE, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dblarr1 attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dblarr1 attribute  ' // trim(name) )
  end subroutine GetHDF5Attribute_dblarr1

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
    integer :: STATUS                   ! Flag from HDF5
    integer(kind=hSize_t), dimension(1) :: MaxShp, SHP  ! Shape

    ! Executable code
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    call GetHDF5AttrDims ( itemID, trim(name), shp, maxShp )
    call allocate_test ( value, int(lb-1+maxShp(1)), 'Value', moduleName )
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute ' // trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_DOUBLE, value, &
      & (/ int(shp,hid_t), int(ones(1:6),hid_t) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dblarr1 attribute ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dblarr1 attribute  ' // trim(name) )
  end subroutine GetHDF5AttributePtr_dblarr1

  ! --------------------------------------------  GetHDF5AttrDims  -----
  subroutine GetHDF5AttrDims ( ItemID, name, DIMS, maxDims )
    integer, intent(in) :: ItemID         ! ItemID
    character (len=*), intent(in) :: NAME ! Name of DS
    integer(kind=hSize_t), dimension(:), intent(out) :: DIMS ! Values of dimensions
    integer(kind=hSize_t), dimension(:), optional, intent(out) :: MAXDIMS ! max Values

    ! Local variables
    integer :: dspace_id                ! spaceID for Attr
    integer(kind=hSize_t), dimension(:), pointer :: maxdims_ptr
    integer :: my_rank, rank
    integer :: AttrID                   ! ID for Attr
    integer :: STATUS                   ! Flag

    ! Executable code
    ! Initializing values returned if there was trouble
    dims = -1
    if ( present(maxDims)) maxDims = -1
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting dims of ' // trim(name) )
    call h5aOpen_name_f ( itemID, trim(name), attrID, status )
    call h5dget_space_f ( attrID, dspace_id, status )
    call h5sget_simple_extent_ndims_f ( dspace_id, rank, status )
    my_rank = min(rank, size(dims))
    allocate ( maxdims_ptr(my_rank) )
    call h5sget_simple_extent_dims_f ( dspace_id, dims(1:my_rank), &
         maxdims_ptr, status )
    if ( status /= my_rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dims for datset ' // trim(name) )
    call h5aClose_f ( attrID, status )
    if ( present(maxDims) ) maxdims = maxdims_ptr(1:my_rank)
    deallocate ( maxdims_ptr )
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting dims of ' // trim(name) )
  end subroutine GetHDF5AttrDims

  ! ----------------------------------------------  GetHDF5DSDims  -----
  subroutine GetHDF5DSDims ( FileID, name, DIMS, maxDims )
    integer, intent(in) :: FILEID         ! fileID
    character (len=*), intent(in) :: NAME ! Name of DS
    integer(kind=hSize_t), dimension(:), intent(out) :: DIMS ! Values of dimensions
    integer(kind=hSize_t), dimension(:), optional, intent(out) :: MAXDIMS ! max Values

    ! Local variables
    integer :: dspace_id                ! spaceID for DS
    integer(kind=hSize_t), dimension(:), pointer :: maxdims_ptr
    integer :: my_rank, rank
    integer :: SETID                    ! ID for DS
    integer :: STATUS                   ! Flag

    ! Executable code
    ! Initializing values returned if there was trouble
    dims = -1
    if ( present(maxDims)) maxDims = -1
    call h5eSet_auto_f ( 0, status )
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
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting dims ' // trim(name) )
  end subroutine GetHDF5DSDims

  ! ----------------------------------------------  GetHDF5DSRank  -----
  subroutine GetHDF5DSRank ( FileID, name, rank )
    integer, intent(in) :: FILEID       ! fileID
    character (len=*), intent(in) :: NAME ! Name of DS
    integer, intent(out) :: rank        ! How many dimensions

    ! Local variables
    integer :: dspace_id                ! spaceID for DS
    integer :: SETID                    ! ID for DS
    integer :: STATUS                   ! Flag

    ! Executable code
    rank = -1                           ! means trouble
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting rank ' // trim(name) )
    call h5dOpen_f ( FileID, trim(name), setID, status )
    call h5dget_space_f(setID,dspace_id,status)
    call h5sget_simple_extent_ndims_f(dspace_id,rank,status)
    call h5dClose_f ( setID, status )
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting rank ' // trim(name) )
  end subroutine GetHDF5DSRank

  ! ---------------------------------------------  GetHDF5DSQType  -----
  subroutine GetHDF5DSQType ( FileID, name, Qtype )
    integer, intent(in) :: FILEID       ! fileID
    character (len=*), intent(in) :: NAME ! Name of DS
    character (len=*), intent(out) :: Qtype    ! 'real' or 'integer' or ..

    ! Local variables
    integer :: type_id                  ! typeID for DS
    integer :: SETID                    ! ID for DS
    integer :: STATUS                   ! Flag

    ! Executable code
    Qtype = 'unknown'                   ! means trouble
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting rank ' // trim(name) )
    call h5dOpen_f ( FileID, trim(name), setID, status )
    call h5dget_type_f(setID,type_id,status)
    if ( AreThe2TypesEqual(type_id, H5T_STD_I32LE) .or. &
      &  AreThe2TypesEqual(type_id, H5T_STD_I32LE) ) then
      Qtype = 'integer'
    else if ( AreThe2TypesEqual(type_id, H5T_NATIVE_CHARACTER) ) then
      Qtype = 'character'
    else if ( AreThe2TypesEqual(type_id, H5T_NATIVE_REAL) .or. &
      &      AreThe2TypesEqual(type_id, H5T_IEEE_F32LE) ) then
      Qtype = 'real'
    else if ( AreThe2TypesEqual(type_id, H5T_NATIVE_DOUBLE) .or. &
      &      AreThe2TypesEqual(type_id, H5T_IEEE_F64LE) ) then
      Qtype = 'double'
    end if
    call h5dClose_f ( setID, status )
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting rank ' // trim(name))
  end subroutine GetHDF5DSQType

  ! --------------------------------------  IsHDF5AttributeInFile  -----
  logical function IsHDF5AttributeInFile ( filename, DSname, name )
    ! This routine returns true if the given HDF5 DS is present
    character (len=*), intent(in) :: FILENAME ! Where to look
    character (len=*), intent(in) :: DSNAME ! Name for the dataset
    character (len=*), intent(in) :: NAME ! Name for the attribute
    ! Local variables
    integer :: fileID                   ! Where to look
    integer :: ATTRID                   ! ID for attribute if present
    integer :: SETID                    ! ID for DS if present
    integer :: STATUS                   ! Flag

    ! Executable code
    IsHDF5AttributeInFile = .false.
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for attribute ' // trim(name) )
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, fileID, status)
    if ( status == 0 ) then
      call h5dOpen_f ( fileid, trim(DSname), setID, status )
      if ( status == 0 ) then
        call h5aopen_name_f(setID, trim(name), attrid, status)
        if ( status == 0 ) then
          IsHDF5AttributeInFile = .true.
          call h5aclose_f(attrid, status)
        end if
        call h5dClose_f ( setID, status )
      end if
      call h5fclose_f(fileID, status)
    end if
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute ' // trim(name) )
  end function IsHDF5AttributeInFile

  ! ---------------------------------------------  IsHDF5DSInFile  -----
  logical function IsHDF5DSInFile ( filename, name )
    ! This routine returns true if the given HDF5 DS is present
    character (len=*), intent(in) :: FILENAME ! Where to look
    character (len=*), intent(in) :: NAME ! Name for the dataset
    ! Local variables
    integer :: FILEID                   ! Where to look
    integer :: SETID                    ! ID for DS if present
    integer :: STATUS                   ! Flag

    ! Executable code
    IsHDF5DSInFile = .false.
    call h5eSet_auto_f ( 0, status )
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
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for DS ' // trim(name) )
  end function IsHDF5DSInFile

  ! -----------------------------  IsHDF5AttributePresent_in_DSID  -----
  logical function IsHDF5AttributePresent_in_DSID ( SETID, name )
    ! This routine returns true if the given HDF5 attribute is present
    integer, intent(in) :: SETID        ! Dataset ID--Where to look
    character (len=*), intent(in) :: NAME ! Name for the attribute
    ! Local variables
    integer :: ATTRID                   ! ID for attribute if present
    integer :: STATUS                   ! Flag

    ! Executable code
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for attribute ' // trim(name) )
    call h5aOpen_name_f ( SETID, name, attrID, status )
    if ( status /= 0 ) then
      IsHDF5AttributePresent_in_DSID = .false.
    else
      IsHDF5AttributePresent_in_DSID = .true.
      call h5aClose_f ( attrID, status )
    end if
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute ' // trim(name) )
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
    integer :: SETID                    ! ID for dataset if present
    integer :: ATTRID                   ! ID for attribute if present
    integer :: STATUS                   ! Flag
    logical :: my_grpattr

    ! Executable code
    my_grpattr = .false.
    if ( present(is_grpattr) ) my_grpattr = is_grpattr
    if ( my_grpattr ) then
      IsHDF5AttributePresent_in_fID = IsHDF5AttributePresent_in_grp ( &
        & DSname, fileID, name)
      return
    end if
    call h5eSet_auto_f ( 0, status )
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
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute ' // trim(name) )
  end function IsHDF5AttributePresent_in_fID

  ! ------------------------------  IsHDF5AttributePresent_in_grp  -----
  logical function IsHDF5AttributePresent_in_grp ( grpName, fileID, name )
    ! This routine returns true if the given HDF5 attribute is present
    ! (Note the unusual ordering of args--to make unambiguous)
    integer, intent(in) :: fileID        ! file ID--Where to look
    character (len=*), intent(in) :: grpName ! Name for the group
    character (len=*), intent(in) :: NAME ! Name for the attribute
    integer :: GRPID                    ! ID for group if present
    integer :: ATTRID                   ! ID for attribute if present
    integer :: STATUS                   ! Flag

    ! Executable code
    call h5eSet_auto_f ( 0, status )
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
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute ' // trim(name) )
  end function IsHDF5AttributePresent_in_grp

  ! --------------------------------------------  IsHDF5DSPresent  -----
  logical function IsHDF5DSPresent ( locID, name )
    ! This routine returns true if the given HDF5 DS is present
    integer, intent(in) :: LOCID        ! Where to look
    character (len=*), intent(in) :: NAME ! Name for the dataset
    ! Local variables
    integer :: SETID                    ! ID for DS if present
    integer :: STATUS                   ! Flag

    ! Executable code
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for DS ' // trim(name) )
    call h5dOpen_f ( locID, name, setID, status )
    IsHDF5DSPresent = status == 0
    if ( IsHDF5DSPresent ) call h5dClose_f ( setID, status )
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for DS ' // trim(name) )
  end function IsHDF5DSPresent

  ! -----------------------------------------  IsHDF5GroupPresent  -----
  logical function IsHDF5GroupPresent ( locID, name )
    ! This routine returns true if the given HDF5 DS is present
    integer, intent(in) :: LOCID        ! Where to look
    character (len=*), intent(in) :: NAME ! Name for the dataset
    ! Local variables
    integer :: SETID                    ! ID for DS if present
    integer :: STATUS                   ! Flag

    ! Executable code
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for Group ' // trim(name) )
    call h5GOpen_f ( locID, name, setID, status )
    IsHDF5groupPresent = status == 0
    if ( IsHDF5groupPresent ) call h5GClose_f ( setID, status )
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for Group ' // trim(name) )
  end function IsHDF5GroupPresent

  ! --------------------------------------  SaveAsHDF5DS_charsclr  -----
  subroutine SaveAsHDF5DS_charsclr ( locID, name, value )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID           ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME  ! Name for this dataset
    character (len=*), intent(in) :: VALUE ! The scalar char string

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape
    integer(hid_t) :: s_type_id
    integer(hid_t) :: type_id

    ! Executable code
    ! Create the dataspace
    shp = 1
    call h5sCreate_simple_f ( 1, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for scalar character ' // trim(name) )
    type_id = H5T_NATIVE_CHARACTER
    call h5tcopy_f ( type_id, s_type_id, status )
    call h5tset_size_f ( s_type_id, len(value), status )
    ! Create the dataset
    call h5dCreate_f ( locID, trim(name), s_type_id, spaceID, setID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for scalar character ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, s_type_id, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    call finishSaveDS ( name, status, setID, spaceID )
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

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(1) :: SHP, MAXDIMS        ! Shape
    integer :: STRINGTYPE               ! Type for string

    ! Executable code
    shp = shape(value)
    maxdims = shp
    if ( present(finalShape) ) maxdims = finalShape
    ! Create a data type for this string
    call h5tcopy_f ( H5T_NATIVE_CHARACTER, stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype for array' // trim(name) )
    call h5tset_size_f ( stringtype, max(len(value(1)), 1), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype ' // trim(name) )
    ! Create the dataspace
 !    call h5sCreate_simple_f ( 1, int(maxdims,hSize_T), spaceID, status )
 !    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
 !      & 'Unable to create dataspace for 1D char array ' // trim(name) )
 !      ! Create the dataset
 !    if ( present(fillValue) ) then
 !      call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, status)
 !      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
 !        & 'Unable to create property list for 1D char array ' // trim(name) )
 !      call h5pset_fill_value_f (cparms, H5T_NATIVE_CHARACTER, fillValue, status)
 !      if ( status /= 0) call MLSMessage (MLSMSG_Error, ModuleName, &
 !        &  "Unable to set Fill value for 1D char array " // trim (name))
 !      call h5dCreate_f ( locID, trim(name), stringtype, spaceID, setID, &
 !        & status, cparms )
 !    else
 !      call h5dCreate_f ( locID, trim(name), stringtype, spaceID, setID, &
 !        & status )
 !    end if
    call createSpaceSet ( locID, name, maxdims, stringtype, &
      & spaceID, setID, status, adding_to, cFill=FillValue)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 1D char array ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, stringtype, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    call finishSaveDS ( name, status, setID, spaceID, stringType )
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

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(2) :: SHP, MAXDIMS        ! Shape
    integer :: STRINGTYPE               ! Type for string

    ! Executable code
    shp = shape(value)
    maxdims = shp
    if ( present(finalShape) ) maxdims = finalShape
    ! Create a data type for this string
    call h5tcopy_f ( H5T_NATIVE_CHARACTER, stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype for array' // trim(name) )
    call h5tset_size_f ( stringtype, max(len(value(1,1)), 1), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype ' // trim(name) )
    ! Create the dataspace
 !    call h5sCreate_simple_f ( 2, int(maxdims,hSize_T), spaceID, status )
 !    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
 !      & 'Unable to create dataspace for 2D char array ' // trim(name) )
 !    ! Create the dataset
 !    if ( present(fillValue) ) then
 !      call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, status)
 !      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
 !        & 'Unable to create property list for 2D char array ' // trim(name) )
 !      call h5pset_fill_value_f (cparms, H5T_NATIVE_CHARACTER, fillValue, status)
 !      if ( status /= 0) call MLSMessage (MLSMSG_Error, ModuleName, &
 !        &  "Unable to set Fill value for 2D charr array " // trim (name))
 !      call h5dCreate_f ( locID, trim(name), stringtype, spaceID, setID, &
 !        & status, cparms )
 !    else
 !      call h5dCreate_f ( locID, trim(name), stringtype, spaceID, setID, &
 !        & status )
 !    end if
    call createSpaceSet ( locID, name, maxdims, stringtype, &
      & spaceID, setID, status, adding_to, cFill=fillValue)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 2D char array ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, stringtype, value, &
      & int ( (/ shp, ones(1:5) /), hID_T ), status )
    call finishSaveDS ( name, status, setID, spaceID, stringType )
  end subroutine SaveAsHDF5DS_chararr2

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

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(1) :: SHP, MAXDIMS        ! Shape

    ! Executable code
    shp = shape(value)
    maxdims = shp
    if ( present(finalShape) ) maxdims = finalShape
    ! Create the dataspace
    call createSpaceSet ( locID, name, maxdims, H5T_NATIVE_INTEGER, &
      & spaceID, setID, status, adding_to, iFill=fillValue)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 1D integer array ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_INTEGER, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    call finishSaveDS ( name, status, setID, spaceID )
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

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(2) :: SHP, MAXDIMS        ! Shape

    ! Executable code
    shp = shape(value)
    maxdims = shp
    if ( present(finalShape) ) maxdims = finalShape
    ! Create the dataspace
    call createSpaceSet ( locID, name, maxdims, H5T_NATIVE_INTEGER, &
      & spaceID, setID, status, adding_to, iFill=fillValue)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 2D integer array ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_INTEGER, value, &
      & int ( (/ shp, ones(1:5) /), hID_T ), status )
    call finishSaveDS ( name, status, setID, spaceID )
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

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(3) :: SHP, MAXDIMS        ! Shape

    ! Executable code
    shp = shape(value)
    maxdims = shp
    if ( present(finalShape) ) maxdims = finalShape
    ! Create the dataspace
    call createSpaceSet ( locID, name, maxdims, H5T_NATIVE_INTEGER, &
      & spaceID, setID, status, adding_to, iFill=fillValue )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 3D integer array ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_INTEGER, value, &
      & int ( (/ shp, ones(1:4) /), hID_T ), status )
    call finishSaveDS ( name, status, setID, spaceID )
  end subroutine SaveAsHDF5DS_intarr3

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
    integer :: I
    character :: MyFillValue, MyValue(size(value))   ! T = true, F = false

    ! Executable code
    ! Turn fillValue into a character
    myFillValue = 'F'
    if ( present(fillValue) ) then
      if ( fillValue ) myFillValue = 'T'
    end if
    ! Turn value into a character
    do i = 1, size(value)
      myValue(i) = 'F'
      if ( value(i) ) myValue(i) = 'T'
    end do
    call saveAsHDF5DS ( locID, name, myValue, finalShape, myFillValue, adding_to )
  end subroutine SaveAsHDF5DS_logarr1

  ! ---------------------------------------  SaveAsHDF5DS_dblarr1  -----
  subroutine SaveAsHDF5DS_dblarr1 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(in) :: VALUE(:)     ! The array itself
    double precision, optional, intent(in) :: FILLVALUE
    integer, dimension(1), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(1) :: SHP, MAXDIMS        ! Shape

    ! Executable code
    shp = shape(value)
    maxdims = shp
    if ( present(finalShape) ) maxdims = finalShape
    ! Create the dataspace
    call createSpaceSet ( locID, name, maxdims, H5T_NATIVE_DOUBLE, &
      & spaceID, setID, status, adding_to, dFill=fillValue )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 1D double array ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    call finishSaveDS ( name, status, setID, spaceID )
  end subroutine SaveAsHDF5DS_dblarr1

  ! ---------------------------------------  SaveAsHDF5DS_dblarr2  -----
  subroutine SaveAsHDF5DS_dblarr2 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(in) :: VALUE(:,:)  ! The array itself
    double precision, optional, intent(in) :: FILLVALUE
    integer, dimension(2), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(2) :: SHP, MAXDIMS        ! Shape

    ! Executable code
    shp = shape(value)
    maxdims = shp
    if ( present(finalShape) ) maxdims = finalShape
    ! Create the dataspace
    call createSpaceSet ( locID, name, maxdims, H5T_NATIVE_DOUBLE, &
      & spaceID, setID, status, adding_to, dFill=fillValue )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 2D double array ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & int ( (/ shp, ones(1:5) /), hID_T ), status )
    call finishSaveDS ( name, status, setID, spaceID )
  end subroutine SaveAsHDF5DS_dblarr2

  ! ---------------------------------------  SaveAsHDF5DS_dblarr3  -----
  subroutine SaveAsHDF5DS_dblarr3 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, intent(in) :: VALUE(:,:,:)  ! The array itself
    double precision, optional, intent(in) :: FILLVALUE
    integer, dimension(3), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(3) :: SHP, MAXDIMS        ! Shape

    ! Executable code
    shp = shape(value)
    maxdims = shp
    if ( present(finalShape) ) maxdims = finalShape
    ! Create the dataspace
    call createSpaceSet ( locID, name, maxdims, H5T_NATIVE_DOUBLE, &
      & spaceID, setID, status, adding_to, dFill=fillValue )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 3D double array ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & int ( (/ shp, ones(1:4) /), hID_T ), status )
    call finishSaveDS ( name, status, setID, spaceID )
  end subroutine SaveAsHDF5DS_dblarr3

  ! --------------------------------------  SaveAsHDF5DS_snglarr1  -----
  subroutine SaveAsHDF5DS_snglarr1 ( locID, name, value, &
    & finalShape, fillValue, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, intent(in) :: VALUE(:)     ! The array itself
    real, optional, intent(in) :: FILLVALUE
    integer, dimension(1), optional, intent(in) :: FINALSHAPE
    logical, optional, intent(in)     :: adding_to

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(1) :: SHP, MAXDIMS        ! Shape

    ! Executable code
    shp = shape(value)
    maxdims=shp
    if ( present(finalShape) ) maxdims = finalShape
    ! Create the dataspace
    call createSpaceSet ( locID, name, maxdims, H5T_NATIVE_REAL, &
      & spaceID, setID, status, adding_to, rFill=fillValue )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 1D double array ' // trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    call finishSaveDS ( name, status, setID, spaceID )
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
    ! Local variables
    integer(hsize_t), dimension(2) :: chunk_dims, dims, maxdims
    integer :: spaceID                  ! ID for filespace
    integer :: memSpaceID               ! ID for arrayspace
    integer :: filespaceID              ! ID for filespace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape
    integer :: style   ! fixed static (0), dynamic 1st (1), adding to (2)
    ! integer :: test_rank
    ! integer(hsize_t), dimension(7) :: test_dims, test_maxdims
    ! integer(hssize_t), dimension(7) :: test_offset
    ! logical :: test_issimple

    ! Executable code
    style = 0  ! The default
    if ( present(may_add_to) ) then
      if ( may_add_to ) style = 1
    end if
    if ( present(adding_to) ) then
      if ( adding_to ) style = 2
    end if
    if ( DEEBUG ) print *, 'style:   ', style

    ! Create the dataspace
    shp = shape(value)
    maxdims = shp
    dims = shp
    chunk_dims = shp
    ! Create the dataset
    if ( style == 0 ) then
     ! print *, 'Creating a static data space'
      call h5sCreate_simple_f ( 2, int(shp,hSize_T), spaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataspace for 2D real array ' // trim(name) )
      if ( present(fillValue) ) then
        call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, status)
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to create property list for 2D real array ' // trim(name) )
        call h5pset_fill_value_f (cparms, H5T_NATIVE_REAL, fillValue, status)
        if ( status /= 0) call MLSMessage (MLSMSG_Error, ModuleName, &
          &  "Unable to set Fill value for 2D real array " // trim (name))
        call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
          & status, cparms )
      else
        call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
          & status )
      end if
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 2D real array ' // trim(name) )
      memspaceID = spaceID
      filespaceID = spaceID
    else if ( style == 1 ) then
     ! print *, 'Creating a dynamic data space'
      maxdims(2) = H5S_UNLIMITED_F
      dims(2) = max(1, shp(2))
      chunk_dims(2) = 1
      if ( DEEBUG ) print *, 'shape ', shp
      if ( DEEBUG ) print *, 'maxdims ', maxdims
      if ( DEEBUG ) print *, 'dims ', dims
      if ( DEEBUG ) print *, 'chunk_dims ', chunk_dims
      call h5screate_simple_f ( 2, dims, memspaceID, status, maxdims )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create memspace for 2D real array ' // trim(name) )
      spaceID = memspaceID
      call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create property list for 2D real array ' // trim(name) )
      call h5pset_chunk_f(cparms, 2, chunk_dims, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to set chunking 2D real array ' // trim(name) )
      if ( present(fillValue) ) then
        call h5pset_fill_value_f (cparms, H5T_NATIVE_REAL, fillValue, status)
        if ( status /= 0) call MLSMessage (MLSMSG_Error, ModuleName, &
          &  "Unable to set Fill value for 2D real array " // trim (name))
      end if
      call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
        & status, cparms )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 2D real array ' // trim(name) )
      call h5dextend_f ( setID, dims, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to extend dataset creating 2D real array ' // trim(name) )
      call h5dget_space_f(setID, filespaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get filespaceID creating 2D real array ' // trim(name) )
      if ( DEEBUG ) print *, 'filespaceID ', filespaceID
      if ( DEEBUG ) print *, 'memspaceID ', memspaceID
      if ( DEEBUG ) print *, 'cparms ', cparms
      if ( DEEBUG ) print *, 'locID ', locID
      if ( DEEBUG ) print *, 'setID ', setID
    else
     ! print *, 'Reopening/extending a dynamic data space'
      call h5dopen_f ( locID, trim(name), setID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open dataset for 2D real array ' // trim(name) )
     ! print *, 'setID:   ', setID
      call h5dget_space_f ( setID, spaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dataspace for 2D real array ' // trim(name) )
     ! print *, 'spaceID:   ', spaceID
      ! call h5sget_simple_extent_ndims_f ( spaceID, test_rank, status )
      ! if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !  & 'Unable to get space rank for 2D real array ' // trim(name) )
     ! print *, 'rank:   ', test_rank
      ! call h5sget_simple_extent_dims_f ( spaceID, test_dims(1:test_rank), &
      ! &  test_maxdims(1:test_rank), status )
      ! if ( status /= test_rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      ! & 'Unable to get dims for datset ' // trim(name) )
     ! print *, 'status                  : ', status
     ! print *, 'dims (before extending) : ', test_dims(1:test_rank)
     ! print *, 'max_dims                : ', test_maxdims(1:test_rank)
      ! Can't test on status--it's set to the test_rank
      ! test_dims(2) = test_dims(2) + shp(2)
     ! print *, 'dims (sfter extending)  : ', test_dims(1:test_rank)
      ! call h5dextend_f( setID, test_dims, status)
      ! if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !  & 'Unable to extend dims for 2D real array ' // trim(name) )
      if ( .not. present(start) ) &
        & call mls_extend ( setID, shp )
      memspaceID = spaceID
      call h5dget_space_f( setID, filespaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dataspace for 2D real array ' // trim(name) )
     ! print *, 'new spaceID             : ', spaceID
    end if
    if ( present(start) ) then
     ! print *, 'start ', start
      call mls_extend ( setID, Count, start, filespaceID )
      call mls_hyperslab_save ( filespaceID, start, count, stride, block )
!     Check before actually writing
     ! print *, 'About to write'
      ! space defined for data
     ! print *, 'We received mem space (defined for data)', memspaceID
      if ( style == 2 ) &
        & call h5screate_simple_f(2, dims, memspaceID, status, maxdims)
     ! call h5sis_simple_f ( memspaceID, test_issimple, status )
     ! call h5sget_simple_extent_ndims_f ( memspaceID, test_rank, status )
     ! call h5sget_simple_extent_dims_f ( memspaceID, test_dims(1:test_rank), &
     !   &  test_maxdims(1:test_rank), status )
     ! if ( status /= test_rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
     !   & 'Unable to get dims for datset ' // trim(name) )
     ! call h5soffset_simple_f( memspaceID, test_offset(1:test_rank), status)
     ! print *, 'Is simple? ', test_issimple
     ! print *, 'rank    : ', test_rank
     ! print *, 'dims    : ', test_dims(1:test_rank)
     ! print *, 'max_dims: ', test_maxdims(1:test_rank)
     ! print *, 'offsets : ', test_offset(1:test_rank)
      ! space defined for file
     ! print *, 'We received file space', spaceID
     ! call h5sis_simple_f ( spaceID, test_issimple, status )
     ! call h5sget_simple_extent_ndims_f ( spaceID, test_rank, status )
     ! call h5sget_simple_extent_dims_f ( spaceID, test_dims(1:test_rank), &
     !   &  test_maxdims(1:test_rank), status )
     ! if ( status /= test_rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
     !   & 'Unable to get dims for datset ' // trim(name) )
     ! call h5soffset_simple_f( spaceID, test_offset(1:test_rank), status)
     ! print *, 'Is simple? ', test_issimple
     ! print *, 'rank    : ', test_rank
     ! print *, 'dims    : ', test_dims(1:test_rank)
     ! print *, 'max_dims: ', test_maxdims(1:test_rank)
     ! print *, 'offsets : ', test_offset(1:test_rank)
      call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:5) /), hID_T ), status, &
        & memspaceID, filespaceID )
     ! print *, 'After writing mem space ', memspaceID
      spaceID = filespaceID
    else
      ! Write the data
      call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:5) /), hID_T ), status )
    end if
    if ( present(start) ) then
      call finishSaveDS ( name, status, setID, spaceID, memspaceID=memspaceID )
    else
      call finishSaveDS ( name, status, setID, spaceID )
    end if
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

    ! Local variables
    integer(hsize_t), dimension(3) :: chunk_dims, dims, maxdims
    integer :: spaceID                  ! ID for filespace
    integer :: filespaceID              ! ID for filespace
    integer :: memSpaceID               ! ID for arrayspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(3) :: SHP        ! Shape
    integer :: style   ! fixed static (0), dynamic 1st (1), adding to (2)
    ! integer :: test_rank
    ! integer(hsize_t), dimension(7) :: test_dims, test_maxdims
    ! integer(hssize_t), dimension(7) :: test_offset
    ! logical :: test_issimple

    ! Executable code
    style = 0  ! The default
    if ( present(may_add_to) ) then
      if ( may_add_to ) style = 1
    end if
    if ( present(adding_to) ) then
      if ( adding_to ) style = 2
    end if

    if ( DEEBUG ) print *, 'style:   ', style
    ! Create the dataspace
    shp = shape(value)
    maxdims = shp
    dims = shp
    chunk_dims = shp
    ! Create the dataset
    if ( style == 0 ) then
      call h5sCreate_simple_f ( 3, int(shp,hSize_T), spaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataspace for 3D real array ' // trim(name) )
      ! Create the dataset
      if ( present(fillValue) ) then
        call h5pcreate_f ( H5P_DATASET_CREATE_F, cparms, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to create property list for 3D real array ' // trim(name) )
        call h5pset_fill_value_f ( cparms, H5T_NATIVE_REAL, fillValue, status)
        if ( status /= 0 ) call MLSMessage (MLSMSG_Error, ModuleName, &
          &  "Unable to set Fill value for 3D real array " // trim (name))
        call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
          & status, cparms )
      else
        call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
          & status )
      end if
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 3D real array ' // trim(name) )
      memspaceID = spaceID
      filespaceID = spaceID
    else if ( style == 1 ) then
      maxdims(3) = H5S_UNLIMITED_F
      dims(3) = max(1, shp(3))
      chunk_dims(3) = 1
      if ( DEEBUG ) print *, 'shape ', shp
      if ( DEEBUG ) print *, 'maxdims ', maxdims
      if ( DEEBUG ) print *, 'dims ', dims
      if ( DEEBUG ) print *, 'chunk_dims ', chunk_dims
      ! call h5screate_simple_f(3, dims, filespaceID, status, maxdims)
      ! if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !  & 'Unable to create filespace for 3D real array ' // trim(name) )
      call h5screate_simple_f ( 3, dims, memspaceID, status, maxdims )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create memspace for 3D real array ' // trim(name) )
      spaceID = memspaceID
      call h5pcreate_f ( H5P_DATASET_CREATE_F, cparms, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create property list for 3D real array ' // trim(name) )
      call h5pset_chunk_f(cparms, 3, chunk_dims, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to set chunking 3D real array ' // trim(name) )
      if ( present(fillValue) ) then
        call h5pset_fill_value_f ( cparms, H5T_NATIVE_REAL, fillValue, status )
        if ( status /= 0) call MLSMessage (MLSMSG_Error, ModuleName, &
          &  "Unable to set Fill value for 3D real array " // trim (name))
      end if
      call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, &
        & setID, status, cparms )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 3D real array ' // trim(name) )
      call h5dextend_f ( setID, dims, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to extend dataset creating 3D real array ' // trim(name) )
      call h5dget_space_f(setID, filespaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get filespaceID creating 3D real array ' // trim(name) )
      if ( DEEBUG ) print *, 'filespaceID ', filespaceID
      if ( DEEBUG ) print *, 'memspaceID ', memspaceID
      if ( DEEBUG ) print *, 'cparms ', cparms
      if ( DEEBUG ) print *, 'locID ', locID
      if ( DEEBUG ) print *, 'setID ', setID
    else
      call h5dopen_f ( locID, trim(name), setID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open dataset for 3D real array ' // trim(name) )
      call h5dget_space_f ( setID, spaceID, status)
      ! if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !   & 'Unable to get dataspace for 3D real array ' // trim(name) )
      ! call h5sget_simple_extent_ndims_f ( spaceID, test_rank, status )
      ! if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !   & 'Unable to get space rank for 3D real array ' // trim(name) )
      ! call h5sget_simple_extent_dims_f ( spaceID, test_dims(1:test_rank), &
      !  &  test_maxdims(1:test_rank), status )
      ! Can't test on status--it's set to the test_rank
      ! test_dims(test_rank) = test_dims(test_rank) + shp(test_rank)
      ! call h5dextend_f( setID, test_dims, status)
      ! if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !   & 'Unable to extend dims for 3D real array ' // trim(name) )
      if ( .not. present(start) ) &
        & call mls_extend ( setID, shp )
      memspaceID = spaceID
      call h5dget_space_f ( setID, filespaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dataspace for 3D real array ' // trim(name) )
    end if
    if ( present(start) ) then
      if ( DEEBUG ) print *, 'name ', name
      if ( DEEBUG ) print *, 'shape(value) ', shp
      if ( DEEBUG ) print *, 'start ', start
      if ( DEEBUG ) print *, 'count ', count
      call mls_extend ( setID, Count, start, filespaceID )
      call mls_hyperslab_save ( filespaceID, start, count, stride, block )
      if ( style == 2 ) &
        & call h5screate_simple_f(3, dims, memspaceID, status, maxdims)
      if ( DEEBUG ) print *, 'filespaceID ', filespaceID
      call dump_space(filespaceID)
      if ( DEEBUG ) print *, 'memspaceID ', memspaceID
      call dump_space(memspaceID)
      call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:4) /), hID_T ), status, &
        & memspaceID, filespaceID )
      spaceID = filespaceID
    else
      ! Write the data
      call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:4) /), hID_T ), status )
    end if
    if ( present(start) ) then
      call finishSaveDS ( name, status, setID, spaceID, memspaceID=memspaceID )
    else
      call finishSaveDS ( name, status, setID, spaceID )
    end if
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

    ! Local variables
    integer(hsize_t), dimension(4) :: chunk_dims, dims, maxdims
    integer :: spaceID                  ! ID for filespace
    integer :: filespaceID              ! ID for filespace
    integer :: memSpaceID               ! ID for arrayspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(4) :: SHP        ! Shape
    integer :: style   ! fixed static (0), dynamic 1st (1), adding to (2)

    ! Executable code
    style = 0  ! The default
    if ( present(may_add_to) ) then
      if ( may_add_to ) style = 1
    end if
    if ( present(adding_to) ) then
      if ( adding_to ) style = 2
    end if

    if ( DEEBUG ) print *, 'style:   ', style
    ! Create the dataspace
    shp = shape(value)
    maxdims = shp
    dims = shp
    chunk_dims = shp
    ! Create the dataset
    if ( style == 0 ) then
      call h5sCreate_simple_f ( 4, int(shp,hSize_T), spaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataspace for 4D real array ' // trim(name) )
      ! Create the dataset
      if ( present(fillValue) ) then
        call h5pcreate_f ( H5P_DATASET_CREATE_F, cparms, status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unable to create property list for 4D real array ' // trim(name) )
        call h5pset_fill_value_f ( cparms, H5T_NATIVE_REAL, fillValue, status )
        if ( status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
          &  "Unable to set Fill value for 4D real array " // trim (name) )
        call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
          & status, cparms )
      else
        call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
        & status )
      end if
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 4D real array ' // trim(name) )
      memspaceID = spaceID
      filespaceID = spaceID
    else if ( style == 1 ) then
      maxdims(4) = H5S_UNLIMITED_F
      dims(4) = max(1, shp(4))
      chunk_dims(4) = 1
      if ( DEEBUG ) print *, 'shape ', shp
      if ( DEEBUG ) print *, 'maxdims ', maxdims
      if ( DEEBUG ) print *, 'dims ', dims
      if ( DEEBUG ) print *, 'chunk_dims ', chunk_dims
      call h5screate_simple_f ( 4, dims, memspaceID, status, maxdims )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create memspace for 4D real array ' // trim(name) )
      spaceID = memspaceID
      call h5pcreate_f ( H5P_DATASET_CREATE_F, cparms, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create property list for 4D real array ' // trim(name) )
      call h5pset_chunk_f ( cparms, 4, chunk_dims, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to set chunking 4D real array ' // trim(name) )
      if ( present(fillValue) ) then
        call h5pset_fill_value_f ( cparms, H5T_NATIVE_REAL, fillValue, status )
        if ( status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
          &  "Unable to set Fill value for 4D real array " // trim (name) )
      end if
      call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, &
        & setID, status, cparms )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 4D real array ' // trim(name) )
      call h5dextend_f ( setID, dims, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to extend dataset creating 4D real array ' // trim(name) )
      call h5dget_space_f ( setID, filespaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get filespaceID creating 4D real array ' // trim(name) )
      if ( DEEBUG ) print *, 'filespaceID ', filespaceID
      if ( DEEBUG ) print *, 'memspaceID ', memspaceID
      if ( DEEBUG ) print *, 'cparms ', cparms
      if ( DEEBUG ) print *, 'locID ', locID
      if ( DEEBUG ) print *, 'setID ', setID
    else
      call h5dopen_f ( locID, trim(name), setID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open dataset for 4D real array ' // trim(name) )
      call h5dget_space_f ( setID, spaceID, status )
      if ( .not. present(start) ) call mls_extend ( setID, shp )
      memspaceID = spaceID
      call h5dget_space_f ( setID, filespaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dataspace for 4D real array ' // trim(name) )
    end if
    if ( present(start) ) then
      if ( DEEBUG ) print *, 'name ', name
      if ( DEEBUG ) print *, 'shape(value) ', shp
      if ( DEEBUG ) print *, 'start ', start
      if ( DEEBUG ) print *, 'count ', count
      call mls_extend ( setID, Count, start, filespaceID )
      call mls_hyperslab_save ( filespaceID, start, count, stride, block )
      if ( style == 2 ) &
        & call h5screate_simple_f(4, dims, memspaceID, status, maxdims)
      if ( DEEBUG ) print *, 'filespaceID ', filespaceID
      call dump_space ( filespaceID )
      if ( DEEBUG ) print *, 'memspaceID ', memspaceID
      call dump_space ( memspaceID )
      call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:3) /), hID_T ), status, &
        & memspaceID, filespaceID )
      spaceID = filespaceID
    else
      ! Write the data
      call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:3) /), hID_T ), status )
    end if
    if ( present(start) ) then
      call finishSaveDS ( name, status, setID, spaceID, memspaceID=memspaceID )
    else
      call finishSaveDS ( name, status, setID, spaceID )
    end if
  end subroutine SaveAsHDF5DS_snglarr4

  ! ------------------------------------  LoadFromHDF5DS_chararr1  -----
  subroutine LoadFromHDF5DS_chararr1 ( locID, name, value, &
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

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: STRINGTYPE               ! String type
    integer :: STRINGSIZE               ! String size

    ! Executable code
    shp = shape ( value )
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
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f ( setID, stringtype, value, &
        & (/ shp, ones(1:6) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID, stringType )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, stringtype, value, (/ shp, ones(1:6) /), status )
      call finishLoad ( name, status, spaceID, setID, stringType=stringType )
    end if
  end subroutine LoadFromHDF5DS_chararr1

  ! ---------------------------------  LoadPtrFromHDF5DS_chararr1  -----
  subroutine LoadPtrFromHDF5DS_chararr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), pointer :: VALUE(:) ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: STRINGTYPE               ! String type
    integer :: STRINGSIZE               ! String size
    integer :: LB, UB

    ! Executable code
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
    call allocate_test ( value, ub, 'Value', moduleName, lowBound=lb )
    call h5dread_f ( setID, stringtype, value, (/ int(shp(1)), ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID, stringType=stringType )
  end subroutine LoadPtrFromHDF5DS_chararr1

  ! ------------------------------------  LoadFromHDF5DS_chararr2  -----
  subroutine LoadFromHDF5DS_chararr2 ( locID, name, value, &
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

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: STRINGTYPE               ! String type
    integer :: STRINGSIZE               ! String size

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dGet_type_f ( setID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for 2D char array ' // trim(name) )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for 2D char array ' // trim(name) )
    if ( stringSize > len(value) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for 2D char array ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f ( setID, stringtype, value, &
        & (/ shp, ones(1:5) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID, stringType )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, stringtype, value, (/ shp, ones(1:5) /), status )
      call finishLoad ( name, status, spaceID, setID, stringType=stringType )
    end if
  end subroutine LoadFromHDF5DS_chararr2

  ! ---------------------------------  LoadPtrFromHDF5DS_chararr2  -----
  subroutine LoadPtrFromHDF5DS_chararr2 ( locID, name, value )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), pointer :: VALUE(:,:) ! The array itself

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(2)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: STRINGTYPE               ! String type
    integer :: STRINGSIZE               ! String size

    ! Executable code
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
   call allocate_test ( value, int(shp(1)), int(shp(2)), 'Value', moduleName )
   call h5dread_f ( setID, stringtype, value, (/ int(shp), ones(1:5) /), status )
   call finishLoad ( name, status, spaceID, setID, stringType=stringType )
  end subroutine LoadPtrFromHDF5DS_chararr2

  ! -------------------------------------  LoadFromHDF5DS_intarr1  -----
  subroutine LoadFromHDF5DS_intarr1 ( locID, name, value, &
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

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
        & (/ shp, ones(1:6) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
        & (/ shp, ones(1:6) /), status )
      call finishLoad ( name, status, spaceID, setID )
    end if
  end subroutine LoadFromHDF5DS_intarr1

  ! ----------------------------------  LoadPtrFromHDF5DS_intarr1  -----
  subroutine LoadPtrFromHDF5DS_intarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, pointer :: VALUE(:)          ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: LB, UB

    ! Executable code
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
    call allocate_test ( value, ub, 'Value', moduleName, lowBound=lb )
    call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
      & (/ int(shp), ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID )
  end subroutine LoadPtrFromHDF5DS_intarr1

  ! -------------------------------------  LoadFromHDF5DS_intarr2  -----
  subroutine LoadFromHDF5DS_intarr2 ( locID, name, value, &
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
    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
        & (/ shp, ones(1:5) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
        & (/ shp, ones(1:5) /), status )
      call finishLoad ( name, status, spaceID, setID )
    end if
  end subroutine LoadFromHDF5DS_intarr2

  ! ----------------------------------  LoadPtrFromHDF5DS_intarr2  -----
  subroutine LoadPtrFromHDF5DS_intarr2 ( locID, name, value )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, pointer :: VALUE(:,:)        ! The array itself

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(2)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    call allocate_test ( value, int(shp(1)), int(shp(2)), 'Value', moduleName )
    call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
      & (/ int(shp), ones(1:5) /), status )
    call finishLoad ( name, status, spaceID, setID )
  end subroutine LoadPtrFromHDF5DS_intarr2

  ! -------------------------------------  LoadFromHDF5DS_logarr1  -----
  subroutine LoadFromHDF5DS_logarr1 ( locID, name, value, &
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
    character :: MyValue(size(value))   ! 'F' = false, 'T' = true

    ! Executable code
    call LoadFromHDF5DS ( locID, name, myValue, start, count, stride, block )
    value = myValue == 'T'
  end subroutine LoadFromHDF5DS_logarr1

  ! ----------------------------------  LoadPtrFromHDF5DS_logarr1  -----
  subroutine LoadPtrFromHDF5DS_logarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    logical, pointer :: VALUE(:)          ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    character, pointer :: MyValue(:)      ! 'F' = false, 'T' = true
    integer :: LB, UB

    ! Executable code
    nullify ( myValue )
    call LoadPtrFromHDF5DS ( locID, name, myValue )
    lb = 1
    if ( present(lowBound) ) lb = lowBound
    ub = lb - 1 + size(myValue)
    call allocate_test ( value, ub, 'Value', moduleName, lowBound=lb )
    value = myValue == 'T'
    call deallocate_test ( myValue, 'myValue', moduleName )
  end subroutine LoadPtrFromHDF5DS_logarr1

  ! -------------------------------------  LoadFromHDF5DS_dblarr1  -----
  subroutine LoadFromHDF5DS_dblarr1 ( locID, name, value, &
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

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f(setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shp, ones(1:6) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shp, ones(1:6) /), status )
      call finishLoad ( name, status, spaceID, setID )
    end if
  end subroutine LoadFromHDF5DS_dblarr1

  ! ----------------------------------  LoadPtrFromHDF5DS_dblarr1  -----
  subroutine LoadPtrFromHDF5DS_dblarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, pointer :: VALUE(:) ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: LB, UB

    ! Executable code
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
    call allocate_test ( value, ub, 'Value', moduleName, lowBound=lb )
    call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & (/ int(shp), ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID )
  end subroutine LoadPtrFromHDF5DS_dblarr1

  ! -------------------------------------  LoadFromHDF5DS_dblarr2  -----
  subroutine LoadFromHDF5DS_dblarr2 ( locID, name, value, &
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

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shp, ones(1:5) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shp, ones(1:5) /), status )
      call finishLoad ( name, status, spaceID, setID )
    end if
  end subroutine LoadFromHDF5DS_dblarr2

  ! ----------------------------------  LoadPtrFromHDF5DS_dblarr2  -----
  subroutine LoadPtrFromHDF5DS_dblarr2 ( locID, name, value )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, pointer :: VALUE(:,:) ! The array itself

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(2)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    call allocate_test ( value, int(shp(1)), int(shp(2)), 'Value', moduleName )
    call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & (/ int(shp), ones(1:5) /), status )
    call finishLoad ( name, status, spaceID, setID )
  end subroutine LoadPtrFromHDF5DS_dblarr2

  ! -------------------------------------  LoadFromHDF5DS_dblarr3  -----
  subroutine LoadFromHDF5DS_dblarr3 ( locID, name, value, &
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

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(3) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shp, ones(1:4) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shp, ones(1:4) /), status )
      call finishLoad ( name, status, spaceID, setID )
    end if
  end subroutine LoadFromHDF5DS_dblarr3

  ! ----------------------------------  LoadPtrFromHDF5DS_dblarr3  -----
  subroutine LoadPtrFromHDF5DS_dblarr3 ( locID, name, value )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    double precision, pointer :: VALUE(:,:,:) ! The array itself

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(3)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    call allocate_test ( value, int(shp(1)), int(shp(2)), int(shp(3)), 'Value', &
      & moduleName )
    call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & (/ int(shp), ones(1:4) /), status )
    call finishLoad ( name, status, spaceID, setID )
  end subroutine LoadPtrFromHDF5DS_dblarr3

  ! ------------------------------------  LoadFromHDF5DS_snglarr1  -----
  subroutine LoadFromHDF5DS_snglarr1 ( locID, name, value, &
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

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shp, ones(1:6) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shp, ones(1:6) /), status )
      call finishLoad ( name, status, spaceID, setID )
    end if
  end subroutine LoadFromHDF5DS_snglarr1

  ! ---------------------------------  LoadPtrFromHDF5DS_snglarr1  -----
  subroutine LoadPtrFromHDF5DS_snglarr1 ( locID, name, value, lowBound )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, pointer :: VALUE(:)             ! The array itself
    integer, intent(in), optional :: LowBound

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(1)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: LB, UB

    ! Executable code
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
    call allocate_test ( value, ub, 'Value', moduleName, lowBound=lb )
    call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
      & (/ int(shp), ones(1:6) /), status )
    call finishLoad ( name, status, spaceID, setID )
  end subroutine LoadPtrFromHDF5DS_snglarr1

  ! ------------------------------------  LoadFromHDF5DS_snglarr2  -----
  subroutine LoadFromHDF5DS_snglarr2 ( locID, name, value, &
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

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shp, ones(1:5) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shp, ones(1:5) /), status )
      call finishLoad ( name, status, spaceID, setID )
    end if
  end subroutine LoadFromHDF5DS_snglarr2

  ! ---------------------------------  LoadPtrFromHDF5DS_snglarr2  -----
  subroutine LoadPtrFromHDF5DS_snglarr2 ( locID, name, value )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, pointer :: VALUE(:,:)           ! The array itself

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(2)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    call allocate_test ( value, int(shp(1)), int(shp(2)), 'Value', moduleName )
    call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
      & (/ int(shp), ones(1:5) /), status )
    call finishLoad ( name, status, spaceID, setID )
  end subroutine LoadPtrFromHDF5DS_snglarr2

  ! ------------------------------------  LoadFromHDF5DS_snglarr3  -----
  subroutine LoadFromHDF5DS_snglarr3 ( locID, name, value, &
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

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(3) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shp, ones(1:4) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shp, ones(1:4) /), status )
      call finishLoad ( name, status, spaceID, setID )
    end if
  end subroutine LoadFromHDF5DS_snglarr3

  ! ---------------------------------  LoadPtrFromHDF5DS_snglarr3  -----
  subroutine LoadPtrFromHDF5DS_snglarr3 ( locID, name, value )
    ! This routine allocates an array and loads it with values from a DS
    use Allocate_Deallocate, only: Allocate_Test
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, pointer :: VALUE(:,:,:)         ! The array itself

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(3)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    call allocate_test ( value, int(shp(1)), int(shp(2)), int(shp(3)), 'Value', &
      & moduleName )
    call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
      & (/ int(shp), ones(1:4) /), status )
    call finishLoad ( name, status, spaceID, setID )
  end subroutine LoadPtrFromHDF5DS_snglarr3

  ! ------------------------------------  LoadFromHDF5DS_snglarr4  -----
  subroutine LoadFromHDF5DS_snglarr4 ( locID, name, value, &
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

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(4) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID               ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    if ( present(start) ) then
      call mls_hyperslab ( spaceID, shp, name, memspaceID, &
        & start, count, stride, block )
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shp, ones(1:3) /), status, memspaceID, spaceID )
      call finishLoad ( name, status, spaceID, setID, memspaceID )
    else
      call check_for_fit ( spaceID, shp, name )
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shp, ones(1:3) /), status )
      call finishLoad ( name, status, spaceID, setID )
    end if
  end subroutine LoadFromHDF5DS_snglarr4

  ! ---------------------------------  LoadPtrFromHDF5DS_snglarr4  -----
  subroutine LoadPtrFromHDF5DS_snglarr4 ( locID, name, value )
    ! This routine allocates an array and loads it with values from a DS
    use MLSMessageModule, only: MLSMSG_Allocate
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real, pointer :: VALUE(:,:,:,:)       ! The array itself

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer(hsize_t) :: SHP(4)          ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset ' // trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset ' // trim(name) )
    call get_ds_shape ( spaceID, shp, name )
    allocate ( value(shp(1),shp(2),shp(3),shp(4)), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Allocate // 'Value' )
    call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
      & (/ int(shp), ones(1:3) /), status )
    call finishLoad ( name, status, spaceID, setID )
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
    ! Executable code
    call GetHDF5Attribute ( itemID, name, line )
    l = len_trim(line)
    if ( l > 0 ) then
      index = GetLitIndexFromString ( line(:l) )
    else
      index = 0
    end if
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
    ! Executable code
    call GetHDF5Attribute ( itemID, name, line )
    if ( len_trim ( line ) > 0 ) then
      index = GetStringIndexFromString ( trim(line) )
    else
      index = 0
    end if
  end subroutine ReadStringIndexFromHDF5Attr

  ! -------------------------------  WriteLitIndexAsHDF5Attribute  -----
  subroutine WriteLitIndexAsHDF5Attribute ( itemID, name, index )
    use Intrinsic, only: LIT_INDICES
    integer, intent(in) :: ITEMID        ! Group etc. to make attr. for
    character(len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(in) :: INDEX         ! String index
    ! Executable code
    if ( index == 0 ) then
      call MakeHDF5Attribute ( itemID, name, '' )
    else
      call WriteStringIndexAsHDF5Attribute ( itemID, name, lit_indices ( index ) )
    end if
  end subroutine WriteLitIndexAsHDF5Attribute

  ! ----------------------------  WriteStringIndexAsHDF5Attribute  -----
  subroutine WriteStringIndexAsHDF5Attribute ( itemID, name, index )
    use String_table, only: GET_STRING
    integer, intent(in) :: ITEMID        ! Group etc. to make attr. for
    character(len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(in) :: INDEX         ! String index
    ! Local variables
    character(len=1024) :: LINE
    ! Executable code
    if ( index == 0 ) then
      call MakeHDF5Attribute ( itemID, name, '' )
    else
      call get_string ( index, line, strip=.true., noError=.true. )
      call MakeHDF5Attribute ( itemID, name, trim(line) )
    end if
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
    integer, dimension(:), intent(in) :: value_dims
    character (len=*), intent(in)     :: NAME ! Name for this dataset
    integer(hsize_t), dimension(size(value_dims)) :: dims
    integer                           :: i
    call get_ds_shape ( spaceID, dims, name )
    if ( any ( dims > value_dims ) ) &
      & call my_message ( MLSMSG_Error, ModuleName, &
      & 'Dataspace too large for destination value of ' // trim(name) , &
      & 'dims(space), dims(value)', (/ (int(dims(i)), value_dims(i), i=1, size(dims)) /), &
      & no_pairs=.true. )
  end subroutine Check_for_fit

! -----------------------------------------------  CreateSpaceSet  -----
  subroutine CreateSpaceSet ( locID, name, maxdims, datatype, &
    & spaceID, setID, status, adding_to, cFill, dFill, iFill, rFill, chunk_dims )
  ! Creates or optionally attaches dataspace and dataset
    integer, intent(in)               :: locID
    character (len=*), intent(in)     :: NAME ! Name for this dataset
    integer, dimension(:), intent(in) :: maxdims
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
    logical :: my_adding_to
    logical :: my_fill
    integer(hsize_t), dimension(size(maxdims,1)) :: my_chunkdims, my_maxdims
    integer :: Rank
    character(len=1) :: cFilled
    double precision :: dFilled
    integer :: iFilled
    real :: rFilled
    ! Executable
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
  subroutine FinishLoad ( Name, Status, SpaceID, SetID, MemspaceID, StringType )
    ! Checks status and closes stuff to finish Load...
    character(len=*), intent(in) :: Name    ! of the dataset
    integer, intent(inout) :: Status        ! from last read operation
    integer, intent(in) :: SpaceID          ! dataspace ID
    integer, intent(in) :: SetID            ! dataset ID
    integer, intent(in), optional :: MemspaceID ! dataspace ID
    integer, intent(in), optional :: StringType ! stringtype ID
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset ' // trim(name) )
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset ' // trim(name) )
    if ( present(memspaceID) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset ' // trim(name) )
    end if
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset ' // trim(name) )
    if ( present(stringType) ) then
      call h5tClose_f ( stringType, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close string type for ' // trim(name) )
    end if
  end subroutine FinishLoad

! ---------------------------------------------  FinishMakeAttrib  -----
  subroutine FinishMakeAttrib ( Name, Status, AttrID, DSID, StringType )
  ! Check status and close stuff for MakeHDF5Attribute_....
    character(len=*), intent(in) :: Name        ! of the attrib
    integer, intent(inout) :: Status            ! from writing the attrib
    integer, intent(in) :: AttrID               ! attrib ID
    integer, intent(in) :: DSID                 ! dataset ID
    integer, intent(in), optional :: StringType ! stringType ID
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
  end subroutine Get_DS_Shape

! ---------------------------------------------------  MLS_extend  -----
  subroutine MLS_extend ( setID, newCount, start, dataSpaceID )
  ! Checks whether we need to extend setID to accommodate newDims
  ! plus any offsets in start array
  ! Does the extending if necessary
    integer, intent(in)                :: setID
    integer, dimension(:), intent(in)  :: newCount
    integer, dimension(:), optional, intent(in) :: start
    integer, optional, intent(inout) :: dataSpaceID
                                 ! Starting coordinatess of hyperslab
  ! Local variables
    integer                           :: spaceID
    integer                           :: rank
    integer(hsize_t), dimension(7)    :: dims, maxdims, my_start
    integer                           :: status
    integer                           :: i
    logical                           :: itFits
    logical                           :: is_simple
  ! Executable code
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
    if ( itFits ) return
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
    integer, dimension(:), intent(in)  :: value_dims
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

    ! Begin execution
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

    ! Begin execution
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
    logical, parameter          :: clean=.false.
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
          & clean=clean, format=int_format, width=width )
      else if ( present(reals) ) then
        call dump ( reals, names, clean=clean )
      else if ( present(doubles) ) then
        call dump ( doubles, names, clean=clean )
      end if
    else
      if ( present(ints) ) then
        call dump_name_v_pairs ( ints, names, &
          & clean=clean, format=int_format, width=width )
      else if ( present(reals) ) then
        call dump_name_v_pairs ( reals, names, &
          & clean=clean, format=real_format, width=width )
      else if ( present(doubles) ) then
        call dump_name_v_pairs ( doubles, names, &
          & clean=clean, format=dbl_format, width=width )
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
    integer, intent(in), optional :: Shp(:) ! Shape of array attrib, else scalar
    ! Local variables
    logical :: my_skip
    integer :: STATUS                   ! Flag from HDF5
    startMakeAttrib = .false.
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip = skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_dsID(itemID, name) ) return
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
  end function StartMakeAttrib

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSHDF5

! $Log$
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
