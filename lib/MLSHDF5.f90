! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
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
  ! Lets break down our use, parameters first
  use HDF5, only: H5F_ACC_RDONLY_F, &
    & H5P_DATASET_CREATE_F, &
    & H5SIS_SIMPLE_F, H5SOFFSET_SIMPLE_F, &
    & H5S_SCALAR_F, H5S_SELECT_SET_F, H5S_UNLIMITED_F, &
    & H5T_IEEE_F32LE, H5T_IEEE_F64LE, &
    & H5T_NATIVE_DOUBLE, H5T_NATIVE_REAL, H5T_STD_I32LE, &
    & H5T_NATIVE_CHARACTER, H5T_NATIVE_INTEGER, &
    & HID_T, HSIZE_T, HSSIZE_T
  ! Now routines
  use HDF5, only: H5ACLOSE_F, H5ACREATE_F, H5AGET_TYPE_F, H5AOPEN_NAME_F, &
    & H5AREAD_F, H5AWRITE_F, &
    & H5DCREATE_F, H5DEXTEND_F, H5DGET_SPACE_F, H5DGET_TYPE_F, H5DOPEN_F, &
    & H5DREAD_F, H5DWRITE_F, H5DCLOSE_F, &
    & H5ESET_AUTO_F, &
    & H5FOPEN_F, H5FCLOSE_F, &
    & H5GOPEN_F, H5GCLOSE_F, &
    & H5PCREATE_F, H5PSET_CHUNK_F, &
    & H5SCLOSE_F, &
    & H5SCREATE_F, H5SCREATE_SIMPLE_F, H5SGET_SIMPLE_EXTENT_NDIMS_F, &
    & H5SGET_SIMPLE_EXTENT_DIMS_F, H5SSELECT_HYPERSLAB_F, &
    & H5TCLOSE_F, H5TCOPY_F, H5TEQUAL_F, H5TGET_SIZE_F, H5TSET_SIZE_F
  use MLSCommon, only: r4, r8
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_WARNING

  implicit NONE
  private

  public :: MakeHDF5Attribute, SaveAsHDF5DS, IsHDF5AttributePresent, &
    & IsHDF5DSPresent, GetHDF5Attribute, LoadFromHDF5DS, &
    & IsHDF5DSInFile, IsHDF5AttributeInFile, &
    & GetHDF5DSRank, GetHDF5DSDims, GetHDF5DSQType, &
    & ReadLitIndexFromHDF5Attr, ReadStringIndexFromHDF5Attr, &
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

! GetHDF5Attribute     Retrieves an attribute
! GetHDF5DSRank        How many dimensions in dataset
! GetHDF5DSDims        Size of the dimensions in dataset
! GetHDF5DSQType       What datatype is dataset?
! IsHDF5...Present     Is the (attribute, DS) in the locid?
! IsHDF5...InFile      Is the (attribute, DS) in the named file?
! LoadFromHDF5DS       Retrieves a dataset
! MakeHDF5Attribute    Turns an arg into an attribute
! SaveAsHDF5DS         Turns an array into a dataset
! === (end of toc) ===

! === (start of api) ===
! GetHDF5Attribute (int itemID, char name, value) 
! GetHDF5DSDims (int FileID, char name, hsize_t dims(:), [hsize_t maxdims(:)]) 
! GetHDF5DSRank (int FileID, char name, int rank) 
! GetHDF5DSQType (int FileID, char name, char QType) 
! log IsHDF5AttributeInFile (char filename, char DSname, char name) 
! log IsHDF5DSInFile (char filename, char name) 
! log IsHDF5AttributePresent (int setid, char name) 
! log IsHDF5AttributePresent (int fileid, char DSname, char name) 
! log IsHDF5DSPresent (int locID, char name) 
! LoadFromHDF5DS (int locID, char name, value,  
!       [int start(:), int count(:), [int stride(:), int block(:)] ] ) 
! MakeHDF5Attribute (int itemID, char name, value) 
! SaveAsHDF5DS (int locID, char name, value) 
!     value can be one of:
!    {char* value, int value, r4 value, r8 value, 
!     int value(:),
!     r4 value(:), r4 value(:,:), r4 value(:,:,:),
!     r8 value(:), r8 value(:,:), r8 value(:,:,:)}
! === (end of api) ===
  interface MakeHDF5Attribute
    module procedure MakeHDF5Attribute_dbl, &
      & MakeHDF5Attribute_int, MakeHDF5Attribute_logical, &
      & MakeHDF5Attribute_string, MakeHDF5Attribute_snglarr1, &
      & MakeHDF5Attribute_dblarr1, MakeHDF5Attribute_string_arr1, &
      & MakeHDF5AttributeDSN_int, &
      & MakeHDF5AttributeDSN_string, MakeHDF5AttributeDSN_snglarr1, &
      & MakeHDF5AttributeDSN_st_arr1, MakeHDF5AttributeDSN_dblarr1
  end interface

  interface GetHDF5Attribute
    module procedure GetHDF5Attribute_int, GetHDF5Attribute_logical, &
      & GetHDF5Attribute_string, GetHDF5Attribute_sngl, GetHDF5Attribute_dbl, &
      & GetHDF5Attribute_snglarr1
  end interface

  interface IsHDF5AttributePresent
    module procedure IsHDF5AttributePresent_in_fID, &
      & IsHDF5AttributePresent_in_DSID, IsHDF5AttributePresent_in_grp
  end interface

  interface SaveAsHDF5DS
    module procedure SaveAsHDF5DS_intarr1, SaveAsHDF5DS_intarr2, &
      & SaveAsHDF5DS_dblarr1, SaveAsHDF5DS_dblarr2, &
      & SaveAsHDF5DS_snglarr1, SaveAsHDF5DS_snglarr2, SaveAsHDF5DS_snglarr3, &
      & SaveAsHDF5DS_snglarr4, &
      & SaveAsHDF5DS_charsclr, SaveAsHDF5DS_chararr1, SaveAsHDF5DS_chararr2
  end interface

  interface LoadFromHDF5DS
    module procedure LoadFromHDF5DS_intarr1, LoadFromHDF5DS_intarr2, &
      & LoadFromHDF5DS_dblarr1, LoadFromHDF5DS_dblarr2, LoadFromHDF5DS_dblarr3, &
      & LoadFromHDF5DS_snglarr1, LoadFromHDF5DS_snglarr2, &
      & LoadFromHDF5DS_snglarr3, LoadFromHDF5DS_snglarr4, &
      & LoadFromHDF5DS_chararr1, LoadFromHDF5DS_chararr2
  end interface

  ! Local parameters
  integer, dimension(7) :: ones = (/1,1,1,1,1,1,1/)
  logical, parameter    :: DEEBUG = .false.
  integer, save :: cantGetDataspaceDims = 0
  integer, parameter :: MAXNUMWARNS = 40

contains ! ======================= Public Procedures =========================

  ! ------------------------------------- MakeHDF5Attribute_dbl
  subroutine MakeHDF5Attribute_dbl ( itemID, name, value, &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    real(r8), intent(in) :: VALUE        ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: DSID                     ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_DSID(itemID, name) ) return
    endif
    call h5sCreate_F ( h5s_scalar_f, dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for attribute '//trim(name) )
    ! Now create the attribute
    call h5aCreate_f ( itemID, trim(name), H5T_NATIVE_DOUBLE, dsID, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create attribute '//trim(name) )
    ! Write
    call h5aWrite_f ( attrID, H5T_NATIVE_DOUBLE, value, ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write attribute '//trim(name) )
    ! Finish off
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute '//trim(name) )
    call h5sClose_f ( dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute dataspace '//trim(name) )
  end subroutine MakeHDF5Attribute_dbl

  ! ------------------------------------- MakeHDF5Attribute_int
  subroutine MakeHDF5Attribute_int ( itemID, name, value, &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(in) :: VALUE        ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: DSID                     ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_DSID(itemID, name) ) return
    endif
    call h5sCreate_F ( h5s_scalar_f, dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for attribute '//trim(name) )
    ! Now create the attribute
    call h5aCreate_f ( itemID, trim(name), H5T_NATIVE_INTEGER, dsID, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create attribute '//trim(name) )
    ! Write
    call h5aWrite_f ( attrID, H5T_NATIVE_INTEGER, value, ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write attribute '//trim(name) )
    ! Finish off
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute '//trim(name) )
    call h5sClose_f ( dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute dataspace '//trim(name) )
  end subroutine MakeHDF5Attribute_int

  ! ------------------------------------- MakeHDF5Attribute_logical
  subroutine MakeHDF5Attribute_logical ( itemID, name, value , &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(in) :: VALUE        ! Value of attribute
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: IVALUE                   ! Value as integer
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_DSID(itemID, name) ) return
    endif
    iValue = 0
    if ( value ) iValue = 1
    call MakeHDF5Attribute ( itemID, name, iValue )
  end subroutine MakeHDF5Attribute_logical

  ! ------------------------------------- MakeHDF5Attribute_string
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

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_DSID(itemID, name) ) return
    endif
    ! Setup
    ! Create a data type for this string
    call h5tcopy_f( H5T_NATIVE_CHARACTER, stringtype, status ) 
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype '//trim(name) )
    call h5tset_size_f(stringtype, max(len_trim(value), 1), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype '//trim(name) )
    ! Create dataspace and attribute
    call h5sCreate_F ( h5s_scalar_f, dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for attribute '//trim(name) )
    ! print *, 'itemID: ', itemID
    ! print *, 'stringtype: ', stringtype
    ! print *, 'dsID: ', dsID
    ! print *, 'name: ', trim(name)
    call h5aCreate_f ( itemID, trim(name), stringtype, dsID, attrID, status )
    ! print *, 'attrID: ', attrID
    ! print *, 'status: ', status
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create attribute '//trim(name) )
    ! Write
    call h5aWrite_f ( attrID, stringtype, value, &
      & ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write attribute '//trim(name) )
    ! Finish off
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute '//trim(name) )
    call h5sClose_f ( dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute dataspace '//trim(name) )
    call h5tClose_f ( stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close stringtype '//trim(name) )
  end subroutine MakeHDF5Attribute_string

  ! ------------------------------------- MakeHDF5Attribute_string_arr1
  subroutine MakeHDF5Attribute_string_arr1 ( itemID, name, value , &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
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
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_DSID(itemID, name) ) return
    endif
    ! Setup
    shp = shape(value)
    ! Create a data type for this string
    call h5tcopy_f( H5T_NATIVE_CHARACTER, stringtype, status ) 
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype for array'//trim(name) )
    call h5tset_size_f(stringtype, max(len_trim(value(1)), 1), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype '//trim(name) )
    ! Create dataspace and attribute
    ! call h5sCreate_F ( h5s_simple_f, dsID, status )
    call h5sCreate_simple_f ( 1, int(shp, hSize_T), dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for attribute '//trim(name) )
    call h5aCreate_f ( itemID, trim(name), stringtype, dsID, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create attribute '//trim(name) )
    ! Write
    call h5aWrite_f ( attrID, stringtype, value, &
      & ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write attribute '//trim(name) )
    ! Finish off
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute '//trim(name) )
    call h5sClose_f ( dsID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute dataspace '//trim(name) )
    call h5tClose_f ( stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close stringtype '//trim(name) )
  end subroutine MakeHDF5Attribute_string_arr1

  ! ------------------------------------- MakeHDF5Attribute_snglarr1
  subroutine MakeHDF5Attribute_snglarr1 ( itemID, name, value , &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    real(r4), intent(in) :: VALUE(:)     ! The attribute array itself
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: spaceID                  ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_DSID(itemID, name) ) return
    endif
    shp = shape(value)
    call h5sCreate_simple_f ( 1, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for 1D attribute array '//trim(name) )
    ! Now create the attribute
    call h5aCreate_f ( itemID, trim(name), H5T_NATIVE_REAL, spaceID, &
      & attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create 1d attribute array '//trim(name) )
    ! Write
    call h5aWrite_f ( attrID, H5T_NATIVE_REAL, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write 1d attribute array '//trim(name) )
    ! Finish off
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close 1d attribute array'//trim(name) )
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close 1d attribute array dataspace '//trim(name) )
  end subroutine MakeHDF5Attribute_snglarr1

  ! ------------------------------------- MakeHDF5Attribute_dblarr1
  subroutine MakeHDF5Attribute_dblarr1 ( itemID, name, value , &
   & skip_if_already_there )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    real(r8), intent(in) :: VALUE(:)     ! The attribute array itself
    logical, intent(in), optional :: skip_if_already_there

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: spaceID                  ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape
    logical :: my_skip

    ! Executable code
    my_skip = .false.
    if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
    if ( my_skip ) then
      if ( IsHDF5AttributePresent_in_DSID(itemID, name) ) return
    endif
    shp = shape(value)
    call h5sCreate_simple_f ( 1, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for 1D attribute array '//trim(name) )
    ! Now create the attribute
    call h5aCreate_f ( itemID, trim(name), H5T_NATIVE_DOUBLE, spaceID, &
      & attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create 1d attribute array '//trim(name) )
    ! Write
    call h5aWrite_f ( attrID, H5T_NATIVE_DOUBLE, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write 1d attribute array '//trim(name) )
    ! Finish off
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close 1d attribute array'//trim(name) )
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close 1d attribute array dataspace '//trim(name) )
  end subroutine MakeHDF5Attribute_dblarr1

  ! ------------------------------------- MakeHDF5AttributeDSN_int
  subroutine MakeHDF5AttributeDSN_int ( fileID, &
   & dataName, attrName, value, skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    integer, intent(in) :: VALUE        ! Value of attribut
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
    endif
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute_int( &
     & dataID, attrName, value )
    call h5dclose_f(dataID, status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data '//trim(dataName) )

  end subroutine MakeHDF5AttributeDSN_int

  ! ------------------------------------- MakeHDF5AttributeDSN_logical
  ! Sorry--this could not be made part of MakeHDF5Attribute
  ! it conflicts with MakeHDF5Attribute_string generic
  ! all becuase of the skip_if_already_there optional argument
  subroutine MakeHDF5AttributeDSN_logical ( fileID, &
   & dataName, attrName, value, skip_if_already_there )
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
    endif
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute_logical( &
     & dataID, attrName, value )
    call h5dclose_f(dataID, status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data '//trim(dataName) )

  end subroutine MakeHDF5AttributeDSN_logical

  ! ------------------------------------- MakeHDF5AttributeDSN_string
  subroutine MakeHDF5AttributeDSN_string ( fileID, &
   & dataName, attrName, value, skip_if_already_there )
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
    endif
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute_string( &
     & dataID, attrName, value )
    call h5dclose_f(dataID, status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data '//trim(dataName) )

  end subroutine MakeHDF5AttributeDSN_string

  ! ------------------------------------- MakeHDF5AttributeDSN_st_arr1
  subroutine MakeHDF5AttributeDSN_st_arr1 ( fileID, &
   & dataName, attrName, value, skip_if_already_there )
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
    endif
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute_string_arr1( &
     & dataID, attrName, value )
    call h5dclose_f(dataID, status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data '//trim(dataName) )

  end subroutine MakeHDF5AttributeDSN_st_arr1

  ! ------------------------------------- MakeHDF5AttributeDSN_snglarr1
  subroutine MakeHDF5AttributeDSN_snglarr1 ( fileID, &
   & dataName, attrName, value, skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    real(r4), intent(in) :: VALUE(:)     ! The attribute array itself
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
    endif
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute_snglarr1( &
     & dataID, attrName, value )
    call h5dclose_f(dataID, status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data '//trim(dataName) )

  end subroutine MakeHDF5AttributeDSN_snglarr1

  ! ------------------------------------- MakeHDF5AttributeDSN_dblarr1
  subroutine MakeHDF5AttributeDSN_dblarr1 ( fileID, &
   & dataName, attrName, value, skip_if_already_there )
    integer, intent(in) :: FILEID       ! FIle where to find them
    character (len=*), intent(in) :: DATANAME ! Name of data set
    character (len=*), intent(in) :: ATTRNAME ! Name of attribute
    real(r8), intent(in) :: VALUE(:)     ! The attribute array itself
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
    endif
    dataID = name_to_dataID( fileID, dataName)
    call MakeHDF5Attribute_dblarr1( &
     & dataID, attrName, value )
    call h5dclose_f(dataID, status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close data '//trim(dataName) )

  end subroutine MakeHDF5AttributeDSN_dblarr1

  ! ------------------------------------------- GetHDF5Attribute_int
  subroutine GetHDF5Attribute_int ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(out) :: VALUE       ! Result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute '//trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read an array into our one value.
    call h5aread_f ( attrID, H5T_NATIVE_INTEGER, value, ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute '//trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute '//trim(name) )
  end subroutine GetHDF5Attribute_int
   
  ! ------------------------------------------- GetHDF5Attribute_string
  subroutine GetHDF5Attribute_string ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    character (len=*), intent(out) :: VALUE ! Result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGTYPE               ! String type
    integer :: STRINGSIZE               ! String size

    ! Executable code
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute '//trim(name) )
    call h5aGet_type_f ( attrID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for attribute '//trim(name) )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for attribute '//trim(name) )
    if ( stringSize > len(value) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for attribute '//trim(name) )
    ! Now actually read the data!
    value = ''
    call h5aread_f ( attrID, stringType, value(1:stringSize), ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read attribute '//trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close attribute '//trim(name) )
    call h5tClose_f ( stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close string type for attribute '//trim(name) )
  end subroutine GetHDF5Attribute_string
   
  ! ------------------------------------- GetHDF5Attribute_logical
  subroutine GetHDF5Attribute_logical ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(out) :: VALUE        ! Value of attribute

    ! Local variables
    integer :: IVALUE                   ! Value as integer

    ! Executable code
    call GetHDF5Attribute ( itemID, name, iValue )
    value = ( iValue == 1 )
  end subroutine GetHDF5Attribute_logical

  ! ------------------------------------------- GetHDF5Attribute_snglarr1
  subroutine GetHDF5Attribute_snglarr1 ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    real(r4), intent(out) :: VALUE(:)     ! The attribute array result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    shp = shape(value)
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute '//trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_REAL, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read 1d attribute array '//trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close 1d attribute array  '//trim(name) )
  end subroutine GetHDF5Attribute_snglarr1
   
  ! ------------------------------------------- GetHDF5Attribute_sngl
  subroutine GetHDF5Attribute_sngl ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    real(r4), intent(out) :: VALUE     ! The attribute result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute '//trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_REAL, value, &
      & ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read sngl attribute '//trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close sngl attribute  '//trim(name) )
  end subroutine GetHDF5Attribute_sngl
   
  ! ------------------------------------------- GetHDF5Attribute_dbl
  subroutine GetHDF5Attribute_dbl ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to get attribute from
    character (len=*), intent(in) :: NAME ! Name of attribute
    real(r8), intent(out) :: VALUE     ! The attribute result

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
    call h5aOpen_name_f ( itemID, name, attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open attribute '//trim(name) )
    ! Note we're going to assume here that the attribute indeed represents the
    ! right type, and that we won't overflow memory etc. by accidentally trying
    ! to read too big array into ours.
    call h5aread_f ( attrID, H5T_NATIVE_DOUBLE, value, &
      & ones, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dble attribute '//trim(name) )
    call h5aClose_f ( attrID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dble attribute  '//trim(name) )
  end subroutine GetHDF5Attribute_dbl
   
  ! ------------------------------------- GetHDF5DSDims
  subroutine GetHDF5DSDims ( FileID, name, DIMS, maxDims )
    integer, intent(in) :: FILEID       ! fileID
    character (len=*), intent(in) :: NAME ! Name of DS
    integer(kind=hSize_t), dimension(:), intent(out) :: DIMS        ! Values of dimensions
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
    if (present(maxDims)) maxDims = -1
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting dims '//trim(name) )
    call h5dOpen_f ( FileID, trim(name), setID, status ) 
    call h5dget_space_f(setID,dspace_id,status)
    call h5sget_simple_extent_ndims_f(dspace_id,rank,status)
    my_rank = min(rank, size(dims))
    allocate(maxdims_ptr(my_rank))
    call h5sget_simple_extent_dims_f(dspace_id,dims(1:my_rank),&
         maxdims_ptr(1:my_rank),status)
    call h5dClose_f ( setID, status )
    if ( present(maxDims) ) then
      maxdims = maxdims_ptr(1:my_rank)
    endif
    deallocate(maxdims_ptr)
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting dims '//trim(name) )
  end subroutine GetHDF5DSDims

  ! ------------------------------------- GetHDF5DSRank
  subroutine GetHDF5DSRank ( FileID, name, rank )
    integer, intent(in) :: FILEID       ! fileID
    character (len=*), intent(in) :: NAME ! Name of DS
    integer, intent(out) :: rank        ! How many dimensions

    ! Local variables
    integer :: dspace_id                ! spaceID for DS
    integer :: SETID                    ! ID for DS
    integer :: STATUS                   ! Flag

    ! Executable code
    rank = -1                          ! means trouble
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting rank '//trim(name) )
    call h5dOpen_f ( FileID, trim(name), setID, status ) 
    call h5dget_space_f(setID,dspace_id,status)
    call h5sget_simple_extent_ndims_f(dspace_id,rank,status)
    call h5dClose_f ( setID, status )
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting rank '//trim(name) )
  end subroutine GetHDF5DSRank

  ! ------------------------------------- GetHDF5DSQType
  subroutine GetHDF5DSQType ( FileID, name, Qtype )
    integer, intent(in) :: FILEID       ! fileID
    character (len=*), intent(in) :: NAME ! Name of DS
    character (len=*), intent(out) :: Qtype    ! 'real' or 'integer' or ..

    ! Local variables
    integer :: type_id                 ! typeID for DS
    integer :: SETID                    ! ID for DS
    integer :: STATUS                   ! Flag

    ! Executable code
    Qtype = 'unknown'                          ! means trouble
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before getting rank '//trim(name) )
    call h5dOpen_f ( FileID, trim(name), setID, status ) 
    call h5dget_type_f(setID,type_id,status)
    if ( AreThe2TypesEqual(type_id, H5T_STD_I32LE) .or. &
      &  AreThe2TypesEqual(type_id, H5T_STD_I32LE) ) then
      Qtype = 'integer'                                        
    elseif ( AreThe2TypesEqual(type_id, H5T_NATIVE_CHARACTER) ) then
      Qtype = 'character'                                      
    elseif ( AreThe2TypesEqual(type_id, H5T_NATIVE_REAL) .or. &
      &      AreThe2TypesEqual(type_id, H5T_IEEE_F32LE) ) then
      Qtype = 'real'                                           
    elseif ( AreThe2TypesEqual(type_id, H5T_NATIVE_DOUBLE) .or. &
      &      AreThe2TypesEqual(type_id, H5T_IEEE_F64LE) ) then
      Qtype = 'double'                                         
    endif                                                               
    call h5dClose_f ( setID, status )
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after getting rank '//trim(name))
  end subroutine GetHDF5DSQType

  ! --------------------------------------------- IsHDF5AttributeInFile ---
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
      & 'Unable to turn error messages off before looking for attribute '//trim(name) )
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, fileID, status)
    if ( status == 0 ) then
      call h5dOpen_f ( fileid, trim(DSname), setID, status )
      if ( status == 0 ) then
        call h5aopen_name_f(setID, trim(name), attrid, status)
        if ( status == 0 ) then
          IsHDF5AttributeInFile = .true.
          call h5aclose_f(attrid, status)
        endif
        call h5dClose_f ( setID, status )
      end if
      call h5fclose_f(fileID, status)
    end if
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute '//trim(name) )
  end function IsHDF5AttributeInFile

  ! --------------------------------------------- IsHDF5DSInFile ---
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
      & 'Unable to turn error messages off before looking for DS '//trim(name) )
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
      & 'Unable to turn error messages back on after looking for DS '//trim(name) )
  end function IsHDF5DSInFile

  ! -------------------------------- IsHDF5AttributePresent_in_DSID ---
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
      & 'Unable to turn error messages off before looking for attribute '//trim(name) )
    call h5aOpen_name_f ( SETID, name, attrID, status ) 
    if ( status /= 0 ) then
      IsHDF5AttributePresent_in_DSID = .false.
    else
      IsHDF5AttributePresent_in_DSID = .true.
      call h5aClose_f ( attrID, status )
    end if
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute '//trim(name) )
  end function IsHDF5AttributePresent_in_DSID

  ! -------------------------------- IsHDF5AttributePresent_in_fID ---
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
    my_grpattr = .true.
    if ( present(is_grpattr) ) my_grpattr = is_grpattr
    if ( my_grpattr ) then
      IsHDF5AttributePresent_in_fID = IsHDF5AttributePresent_in_grp ( &
        & DSname, fileID, name)
      return
    endif
    call h5eSet_auto_f ( 0, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages off before looking for attribute '//trim(name) )
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
    endif
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute '//trim(name) )
  end function IsHDF5AttributePresent_in_fID

  ! -------------------------------- IsHDF5AttributePresent_in_grp ---
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
      & 'Unable to turn error messages off before looking for attribute '//trim(name) )
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
    endif
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for attribute '//trim(name) )
  end function IsHDF5AttributePresent_in_grp

  ! --------------------------------------------- IsHDF5DSPresent ---
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
      & 'Unable to turn error messages off before looking for DS '//trim(name) )
    call h5dOpen_f ( locID, name, setID, status ) 
    if ( status /= 0 ) then
      IsHDF5DSPresent = .false.
    else
      IsHDF5DSPresent = .true.
      call h5dClose_f ( setID, status )
    end if
    call h5eSet_auto_f ( 1, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to turn error messages back on after looking for DS '//trim(name) )
  end function IsHDF5DSPresent

  ! --------------------------------------------- SaveAsHDF5DS_charsclr
  subroutine SaveAsHDF5DS_charsclr ( locID, name, value )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(in) :: VALUE     ! The scalar char string

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
      & 'Unable to create dataspace for 1D integer array '//trim(name) )
    type_id = H5T_NATIVE_CHARACTER
    call h5tcopy_f(type_id, s_type_id, status)
    call h5tset_size_f(s_type_id, len(value), status)
    ! Create the dataset
    call h5dCreate_f ( locID, trim(name), s_type_id, spaceID, setID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 1D integer array '//trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, s_type_id, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write to dataset for 1D integer array '//trim(name) )
    ! Close things
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 1D integer array '//trim(name) )
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 1D integer array '//trim(name) )
  end subroutine SaveAsHDF5DS_charsclr

  ! --------------------------------------------- SaveAsHDF5DS_chararr1
  subroutine SaveAsHDF5DS_chararr1 ( locID, name, value )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(in) :: VALUE(:)     ! The array itself

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape
    integer :: STRINGTYPE               ! Type for string

    ! Executable code
    shp = shape(value)
    ! Create a data type for this string
    call h5tcopy_f( H5T_NATIVE_CHARACTER, stringtype, status ) 
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype for array'//trim(name) )
    call h5tset_size_f(stringtype, max(len_trim(value(1)), 1), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype '//trim(name) )
    ! Create the dataspace
    call h5sCreate_simple_f ( 1, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for 1D char array '//trim(name) )
    ! Create the dataset
    call h5dCreate_f ( locID, trim(name), stringtype, spaceID, setID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 1D char array '//trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, stringtype, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write to dataset for 1D char array '//trim(name) )
    ! Close things
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 1D char array '//trim(name) )
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 1D char array '//trim(name) )
    call h5tClose_f ( stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close stringtype '//trim(name) )
  end subroutine SaveAsHDF5DS_chararr1

  ! --------------------------------------------- SaveAsHDF5DS_chararr2
  subroutine SaveAsHDF5DS_chararr2 ( locID, name, value )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(in) :: VALUE(:,:)     ! The array itself

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape
    integer :: STRINGTYPE               ! Type for string

    ! Executable code
    shp = shape(value)
    ! Create a data type for this string
    call h5tcopy_f( H5T_NATIVE_CHARACTER, stringtype, status ) 
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype for array'//trim(name) )
    call h5tset_size_f(stringtype, max(len_trim(value(1,1)), 1), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype '//trim(name) )
    ! Create the dataspace
    call h5sCreate_simple_f ( 2, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for 2D char array '//trim(name) )
    ! Create the dataset
    call h5dCreate_f ( locID, trim(name), stringtype, spaceID, setID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 2D char array '//trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, stringtype, value, &
      & int ( (/ shp, ones(1:5) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write to dataset for 2D char array '//trim(name) )
    ! Close things
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 2D char array '//trim(name) )
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 2D char array '//trim(name) )
    call h5tClose_f ( stringtype, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close stringtype '//trim(name) )
  end subroutine SaveAsHDF5DS_chararr2

  ! --------------------------------------------- SaveAsHDF5DS_intarr1
  subroutine SaveAsHDF5DS_intarr1 ( locID, name, value )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(in) :: VALUE(:)     ! The array itself

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    ! Create the dataspace
    shp = shape(value)
    call h5sCreate_simple_f ( 1, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for 1D integer array '//trim(name) )
    ! Create the dataset
    call h5dCreate_f ( locID, trim(name), H5T_NATIVE_INTEGER, spaceID, setID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 1D integer array '//trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_INTEGER, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write to dataset for 1D integer array '//trim(name) )
    ! Close things
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 1D integer array '//trim(name) )
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 1D integer array '//trim(name) )
  end subroutine SaveAsHDF5DS_intarr1

  ! --------------------------------------------- SaveAsHDF5DS_intarr2
  subroutine SaveAsHDF5DS_intarr2 ( locID, name, value )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(in) :: VALUE(:,:)     ! The array itself

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape

    ! Executable code
    ! Create the dataspace
    shp = shape(value)
    call h5sCreate_simple_f ( 2, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for 2D integer array '//trim(name) )
    ! Create the dataset
    call h5dCreate_f ( locID, trim(name), H5T_NATIVE_INTEGER, spaceID, setID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 2D integer array '//trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_INTEGER, value, &
      & int ( (/ shp, ones(1:5) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write to dataset for 2D integer array '//trim(name) )
    ! Close things
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 2D integer array '//trim(name) )
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 2D integer array '//trim(name) )
  end subroutine SaveAsHDF5DS_intarr2

  ! --------------------------------------------- SaveAsHDF5DS_dblarr1
  subroutine SaveAsHDF5DS_dblarr1 ( locID, name, value )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r8), intent(in) :: VALUE(:)     ! The array itself

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    ! Create the dataspace
    shp = shape(value)
    call h5sCreate_simple_f ( 1, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for 1D double array '//trim(name) )
    ! Create the dataset
    call h5dCreate_f ( locID, trim(name), H5T_NATIVE_DOUBLE, spaceID, setID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 1D double array '//trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write to dataset for 1D double array '//trim(name) )
    ! Close things
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 1D double array '//trim(name) )
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 1D double array '//trim(name) )
  end subroutine SaveAsHDF5DS_dblarr1

  ! --------------------------------------------- SaveAsHDF5DS_dblarr2
  subroutine SaveAsHDF5DS_dblarr2 ( locID, name, value )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r8), intent(in) :: VALUE(:,:)  ! The array itself

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape

    ! Executable code
    ! Create the dataspace
    shp = shape(value)
    call h5sCreate_simple_f ( 2, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for 2D integer array '//trim(name) )
    ! Create the dataset
    call h5dCreate_f ( locID, trim(name), H5T_NATIVE_DOUBLE, spaceID, setID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 2D integer array '//trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & int ( (/ shp, ones(1:5) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write to dataset for 2D integer array '//trim(name) )
    ! Close things
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 2D integer array '//trim(name) )
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 2D integer array '//trim(name) )
  end subroutine SaveAsHDF5DS_dblarr2

  ! --------------------------------------------- SaveAsHDF5DS_snglarr1
  subroutine SaveAsHDF5DS_snglarr1 ( locID, name, value )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r4), intent(in) :: VALUE(:)     ! The array itself

    ! Local variables
    integer :: spaceID                  ! ID for dataspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape

    ! Executable code
    ! Create the dataspace
    shp = shape(value)
    call h5sCreate_simple_f ( 1, int(shp,hSize_T), spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataspace for 1D double array '//trim(name) )
    ! Create the dataset
    call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
      & status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create dataset for 1D double array '//trim(name) )
    ! Write the data
    call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
      & int ( (/ shp, ones(1:6) /), hID_T ), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to write to dataset for 1D double array '//trim(name) )
    ! Close things
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 1D double array '//trim(name) )
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 1D double array '//trim(name) )
  end subroutine SaveAsHDF5DS_snglarr1

  ! --------------------------------------------- SaveAsHDF5DS_snglarr2
  subroutine SaveAsHDF5DS_snglarr2 ( locID, name, value , &
    & start, count, stride, block, may_add_to, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r4), intent(in) :: VALUE(:,:)  ! The array itself

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
    integer(hid_t) :: cparms
    integer(hsize_t), dimension(2) :: chunk_dims, dims, maxdims
    integer :: spaceID                  ! ID for filespace
    integer :: memSpaceID               ! ID for arrayspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape
    integer :: style   ! fixed static (0), dynamic 1st (1), adding to (2)
    integer :: test_rank
    integer(hsize_t), dimension(7) :: test_dims, test_maxdims
    integer(hssize_t), dimension(7) :: test_offset
    logical :: test_issimple

    ! Executable code
    style = 0  ! The default
    if ( present(may_add_to) ) then
      if ( may_add_to ) style = 1
    endif
    if ( present(adding_to) ) then
      if ( adding_to ) style = 2
    endif
   ! print *, 'style:   ', style
      
    ! Create the dataspace
    shp = shape(value)
    maxdims = shp
    dims = shp
    chunk_dims = shp
    ! Create the dataset
    if (style == 0) then
     ! print *, 'Creating a static data space'
      call h5sCreate_simple_f ( 2, int(shp,hSize_T), spaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataspace for 2D real array '//trim(name) )
      call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
        & status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 2D real array '//trim(name) )
    elseif (style == 1) then
     ! print *, 'Creating a dynamic data space'
      maxdims(2) = H5S_UNLIMITED_F
      dims(2) = max(1, shp(2))
      chunk_dims(2) = 1
      call h5screate_simple_f(2, dims, spaceID, status, maxdims)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataspace for 2D real array '//trim(name) )
      call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create property list for 2D real array '//trim(name) )
      call h5pset_chunk_f(cparms, 2, chunk_dims, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to set chunking 2D real array '//trim(name) )
      call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
        & status, cparms )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 2D real array '//trim(name) )
    else
     ! print *, 'Reopening/extending a dynamic data space'
      call h5dopen_f ( locID, trim(name), setID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open dataset for 2D real array '//trim(name) )
     ! print *, 'setID:   ', setID
      call h5dget_space_f( setID, spaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dataspace for 2D real array '//trim(name) )
     ! print *, 'spaceID:   ', spaceID
      call h5sget_simple_extent_ndims_f ( spaceID, test_rank, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get space rank for 2D real array '//trim(name) )
     ! print *, 'rank:   ', test_rank
      call h5sget_simple_extent_dims_f ( spaceID, test_dims(1:test_rank), &
       &  test_maxdims(1:test_rank), status )
     ! print *, 'status                  : ', status
     ! print *, 'dims (before extending) : ', test_dims(1:test_rank)
     ! print *, 'max_dims                : ', test_maxdims(1:test_rank)
      ! Can't test on status--it's set to the test_rank
      test_dims(2) = test_dims(2) + shp(2)
     ! print *, 'dims (sfter extending)  : ', test_dims(1:test_rank)
      call h5dextend_f( setID, test_dims, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to extend dims for 2D real array '//trim(name) )
      call h5dget_space_f( setID, spaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dataspace for 2D real array '//trim(name) )
     ! print *, 'new spaceID             : ', spaceID
    endif
    if ( present(start) ) then
     ! print *, 'start ', start
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
!     Check before actually writing
     ! print *, 'About to write'
      ! space defined for data
     ! print *, 'We received mem space (defined for data)', memspaceID
      call h5sis_simple_f ( memspaceID, test_issimple, status )
      call h5sget_simple_extent_ndims_f ( memspaceID, test_rank, status )
      call h5sget_simple_extent_dims_f ( memspaceID, test_dims(1:test_rank), &
        &  test_maxdims(1:test_rank), status )
      call h5soffset_simple_f( memspaceID, test_offset(1:test_rank), status)
     ! print *, 'Is simple? ', test_issimple
     ! print *, 'rank    : ', test_rank
     ! print *, 'dims    : ', test_dims(1:test_rank)
     ! print *, 'max_dims: ', test_maxdims(1:test_rank)
     ! print *, 'offsets : ', test_offset(1:test_rank)
      ! space defined for file
     ! print *, 'We received file space', spaceID
      call h5sis_simple_f ( spaceID, test_issimple, status )
      call h5sget_simple_extent_ndims_f ( spaceID, test_rank, status )
      call h5sget_simple_extent_dims_f ( spaceID, test_dims(1:test_rank), &
        &  test_maxdims(1:test_rank), status )
      !if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !  & 'Unable to get space dims for 2D real array '//trim(name) )
      call h5soffset_simple_f( spaceID, test_offset(1:test_rank), status)
     ! print *, 'Is simple? ', test_issimple
     ! print *, 'rank    : ', test_rank
     ! print *, 'dims    : ', test_dims(1:test_rank)
     ! print *, 'max_dims: ', test_maxdims(1:test_rank)
     ! print *, 'offsets : ', test_offset(1:test_rank)
      call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:5) /), hID_T ), status, &
        & memspaceID, spaceID )
     ! print *, 'After writing mem space ', memspaceID
    else
      ! Write the data
      call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:5) /), hID_T ), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &   
      & 'Unable to write to dataset for 2D real array '//trim(name) )  
    ! Close things
   ! print *, 'Before closing other, mem space ', memspaceID
   ! print *, 'Attempting to close data space ', spaceID
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 2D real array '//trim(name) )
    if ( present(start) ) then
     ! print *, 'Attempting to close mem space ', memspaceID
      call h5sClose_F ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close mem dataspace for 2D real array '//trim(name) )
    endif
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 2D real array '//trim(name) )
  end subroutine SaveAsHDF5DS_snglarr2

  ! --------------------------------------------- SaveAsHDF5DS_snglarr3
  subroutine SaveAsHDF5DS_snglarr3 ( locID, name, value, &
    & start, count, stride, block, may_add_to, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r4), intent(in) :: VALUE(:,:,:)  ! The array itself
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
    integer(hid_t) :: cparms
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
    endif
    if ( present(adding_to) ) then
      if ( adding_to ) style = 2
    endif

    if ( DEEBUG ) print *, 'style:   ', style
    ! Create the dataspace
    shp = shape(value)
    maxdims = shp
    dims = shp
    chunk_dims = shp
    ! Create the dataset
    if (style == 0) then
      call h5sCreate_simple_f ( 3, int(shp,hSize_T), spaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataspace for 3D real array '//trim(name) )
      ! Create the dataset
      call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
        & status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 3D real array '//trim(name) )
      memspaceID = spaceID
      filespaceID = spaceID
    elseif (style == 1) then
      maxdims(3) = H5S_UNLIMITED_F
      dims(3) = max(1, shp(3))
      chunk_dims(3) = 1
      if ( DEEBUG ) print *, 'shape ', shp
      if ( DEEBUG ) print *, 'maxdims ', maxdims
      if ( DEEBUG ) print *, 'dims ', dims
      if ( DEEBUG ) print *, 'chunk_dims ', chunk_dims
      ! call h5screate_simple_f(3, dims, filespaceID, status, maxdims)
      ! if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !  & 'Unable to create filespace for 3D real array '//trim(name) )
      call h5screate_simple_f(3, dims, memspaceID, status, maxdims)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create memspace for 3D real array '//trim(name) )
      spaceID = memspaceID
      call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create property list for 3D real array '//trim(name) )
      call h5pset_chunk_f(cparms, 3, chunk_dims, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to set chunking 3D real array '//trim(name) )
      call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, &
        & setID, status, cparms )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 3D real array '//trim(name) )
      call h5dextend_f ( setID, dims, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to extend dataset creating 3D real array '//trim(name) )
      call h5dget_space_f(setID, filespaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get filespaceID creating 3D real array '//trim(name) )
      if ( DEEBUG ) print *, 'filespaceID ', filespaceID
      if ( DEEBUG ) print *, 'memspaceID ', memspaceID
      if ( DEEBUG ) print *, 'cparms ', cparms
      if ( DEEBUG ) print *, 'locID ', locID
      if ( DEEBUG ) print *, 'setID ', setID
    else
      call h5dopen_f ( locID, trim(name), setID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open dataset for 3D real array '//trim(name) )
      call h5dget_space_f( setID, spaceID, status)
      ! if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !   & 'Unable to get dataspace for 3D real array '//trim(name) )
      ! call h5sget_simple_extent_ndims_f ( spaceID, test_rank, status )
      ! if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !   & 'Unable to get space rank for 3D real array '//trim(name) )
      ! call h5sget_simple_extent_dims_f ( spaceID, test_dims(1:test_rank), &
      !  &  test_maxdims(1:test_rank), status )
      ! Can't test on status--it's set to the test_rank
      ! test_dims(test_rank) = test_dims(test_rank) + shp(test_rank)
      ! call h5dextend_f( setID, test_dims, status)
      ! if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      !   & 'Unable to extend dims for 3D real array '//trim(name) )
      if ( .not. present(start) ) &
        & call mls_extend ( setID, shp )
      memspaceID = spaceID
      call h5dget_space_f( setID, filespaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dataspace for 3D real array '//trim(name) )
    endif
    if ( present(start) ) then
      if ( DEEBUG ) print *, 'name ', name
      if ( DEEBUG ) print *, 'shape(value) ', shape(value)
      if ( DEEBUG ) print *, 'start ', start
      if ( DEEBUG ) print *, 'count ', count
      call mls_extend ( setID, Count, start, filespaceID )
      call mls_hyperslab_save(filespaceID, &
        & start, count, stride, block)
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
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &   
      & 'Unable to write to dataset for 3D real array '//trim(name) )  
    ! Close things
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 3D real array '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 3D real array '//trim(name) )
  end subroutine SaveAsHDF5DS_snglarr3

  ! --------------------------------------------- SaveAsHDF5DS_snglarr4
  subroutine SaveAsHDF5DS_snglarr4 ( locID, name, value, &
    & start, count, stride, block, may_add_to, adding_to )
    ! This routine does the initial work of creating a dataset
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r4), intent(in) :: VALUE(:,:,:,:)  ! The array itself
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
    integer(hid_t) :: cparms
    integer(hsize_t), dimension(4) :: chunk_dims, dims, maxdims
    integer :: spaceID                  ! ID for filespace
    integer :: memSpaceID               ! ID for arrayspace
    integer (HID_T) :: setID            ! ID for dataset
    integer :: status                   ! Flag from HDF5
    integer, dimension(4) :: SHP        ! Shape
    integer :: style   ! fixed static (0), dynamic 1st (1), adding to (2)
    integer :: test_rank
    integer(hsize_t), dimension(7) :: test_dims, test_maxdims

    ! Executable code
    style = 0  ! The default
    if ( present(may_add_to) ) then
      if ( may_add_to ) style = 1
    endif
    if ( present(adding_to) ) then
      if ( adding_to ) style = 2
    endif

    ! Create the dataspace
    shp = shape(value)
    maxdims = shp
    dims = shp
    chunk_dims = shp
    ! Create the dataset
    if (style == 0) then
      call h5sCreate_simple_f ( 4, int(shp,hSize_T), spaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataspace for 4D real array '//trim(name) )
      ! Create the dataset
      call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
        & status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 4D real array '//trim(name) )
    elseif (style == 1) then
      maxdims(4) = H5S_UNLIMITED_F
      dims(4) = max(1, shp(4))
      chunk_dims(4) = 1
      call h5screate_simple_f(4, dims, spaceID, status, maxdims)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataspace for 4D real array '//trim(name) )
      call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create property list for 4D real array '//trim(name) )
      call h5pset_chunk_f(cparms, 4, chunk_dims, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to set chunking 4D real array '//trim(name) )
      call h5dCreate_f ( locID, trim(name), H5T_NATIVE_REAL, spaceID, setID, &
        & status, cparms )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to create dataset for 4D real array '//trim(name) )
    else
      call h5dopen_f ( locID, trim(name), setID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open dataset for 4D real array '//trim(name) )
      call h5dget_space_f( setID, spaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dataspace for 4D real array '//trim(name) )
      call h5sget_simple_extent_ndims_f ( spaceID, test_rank, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get space rank for 4D real array '//trim(name) )
      call h5sget_simple_extent_dims_f ( spaceID, test_dims(1:test_rank), &
       &  test_maxdims(1:test_rank), status )
      ! Can't test on status--it's set to the test_rank
      test_dims(test_rank) = test_dims(test_rank) + shp(test_rank)
      call h5dextend_f( setID, test_dims, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to extend dims for 4D real array '//trim(name) )
      call h5dget_space_f( setID, spaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get dataspace for 4D real array '//trim(name) )
    endif
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:3) /), hID_T ), status, &
        & memspaceID, spaceID )
    else
      ! Write the data
      call h5dWrite_f ( setID, H5T_NATIVE_REAL, value, &
        & int ( (/ shp, ones(1:3) /), hID_T ), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &   
      & 'Unable to write to dataset for 4D real array '//trim(name) )  
    ! Close things
    call h5sClose_F ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for 4D real array '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_F ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset for 4D real array '//trim(name) )
  end subroutine SaveAsHDF5DS_snglarr4

  ! ----------------------------------- LoadFromHDF5DS_chararr1
  subroutine LoadFromHDF5DS_chararr1 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    character (len=*), intent(out) :: VALUE(:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: STRINGTYPE               ! String type
    integer :: STRINGSIZE               ! String size

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dGet_type_f ( setID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for 1D char array '//trim(name) )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for 1D char array '//trim(name) )
    if ( stringSize > len(value(1)) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for 1D car array '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, stringtype, value, &
        & (/ shape(value), ones(1:6) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, stringtype, value, &
        & (/ shape(value), ones(1:6) /), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
    call h5tClose_f ( stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close string type for 1D char array '//trim(name) )
  end subroutine LoadFromHDF5DS_chararr1

  ! ----------------------------------- LoadFromHDF5DS_chararr2
  subroutine LoadFromHDF5DS_chararr2 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
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
    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: STRINGTYPE               ! String type
    integer :: STRINGSIZE               ! String size

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dGet_type_f ( setID, stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get type for 2D char array '//trim(name) )
    call h5tGet_size_f ( stringType, stringSize, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get size for 2D char array '//trim(name) )
    if ( stringSize > len(value(1,1)) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Value too long to fit in space given for 2D char array '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, stringtype, value, &
        & (/ shape(value), ones(1:5) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, stringtype, value, &
        & (/ shape(value), ones(1:5) /), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
    call h5tClose_f ( stringType, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close string type for 2D char array '//trim(name) )
  end subroutine LoadFromHDF5DS_chararr2

  ! ----------------------------------- LoadFromHDF5DS_intarr1
  subroutine LoadFromHDF5DS_intarr1 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(out) :: VALUE(:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block
    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, H5T_NATIVE_INTEGER, value, &
        & (/ shape(value), ones(1:6) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
        & (/ shape(value), ones(1:6) /), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_intarr1

  ! ----------------------------------- LoadFromHDF5DS_intarr2
  subroutine LoadFromHDF5DS_intarr2 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
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
    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, H5T_NATIVE_INTEGER, value, &
        & (/ shape(value), ones(1:5) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
        & (/ shape(value), ones(1:5) /), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_intarr2

  ! ----------------------------------- LoadFromHDF5DS_dblarr1
  subroutine LoadFromHDF5DS_dblarr1 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r8), intent(out) :: VALUE(:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shape(value), ones(1:6) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shape(value), ones(1:6) /), status )
    endif
    ! Now, (at last!) read the data
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_dblarr1

  ! ----------------------------------- LoadFromHDF5DS_dblarr2
  subroutine LoadFromHDF5DS_dblarr2 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r8), intent(out) :: VALUE(:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shape(value), ones(1:5) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shape(value), ones(1:5) /), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_dblarr2

  ! ----------------------------------- LoadFromHDF5DS_dblarr3
  subroutine LoadFromHDF5DS_dblarr3 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r8), intent(out) :: VALUE(:,:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(3) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shape(value), ones(1:4) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
        & (/ shape(value), ones(1:4) /), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_dblarr3

  ! ----------------------------------- LoadFromHDF5DS_snglarr1
  subroutine LoadFromHDF5DS_snglarr1 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r4), intent(out) :: VALUE(:)    ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, H5T_NATIVE_REAL, value, &
        & (/ shape(value), ones(1:6) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shape(value), ones(1:6) /), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_snglarr1

  ! ----------------------------------- LoadFromHDF5DS_snglarr2
  subroutine LoadFromHDF5DS_snglarr2 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r4), intent(out) :: VALUE(:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, H5T_NATIVE_REAL, value, &
        & (/ shape(value), ones(1:5) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shape(value), ones(1:5) /), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_snglarr2

  ! ----------------------------------- LoadFromHDF5DS_snglarr3
  subroutine LoadFromHDF5DS_snglarr3 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r4), intent(out) :: VALUE(:,:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(3) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, H5T_NATIVE_REAL, value, &
        & (/ shape(value), ones(1:4) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shape(value), ones(1:4) /), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_snglarr3

  ! ----------------------------------- LoadFromHDF5DS_snglarr4
  subroutine LoadFromHDF5DS_snglarr4 ( locID, name, value, &
    & start, count, stride, block )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID          ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r4), intent(out) :: VALUE(:,:,:,:) ! The array itself
    integer, dimension(:), optional, intent(in) :: start
                                 ! Starting coordinatess of hyperslab
    integer, dimension(:), optional, intent(in) :: count
                                 ! Num of blocks to select from dataspace
    integer, dimension(:), optional, intent(in) :: stride
                                 ! How many elements to move in each direction
    integer, dimension(:), optional, intent(in) :: block
                                 ! Size of element block

    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(4) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: MEMSPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call mls_hyperslab(spaceID, shape(value), name, memspaceID, &
        & start, count, stride, block)
      call h5dread_f(setID, H5T_NATIVE_REAL, value, &
        & (/ shape(value), ones(1:3) /), status, memspaceID, spaceID )
    else
      call check_for_fit(spaceID, shape(value), name)
      call h5dread_f ( setID, H5T_NATIVE_REAL, value, &
        & (/ shape(value), ones(1:3) /), status )
    endif
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    if ( present(start) ) then
      call h5sClose_f ( memspaceID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to close hyperslab dataspace for dataset '//trim(name) )
    endif
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_snglarr4

  ! ------------------------------------------ ReadLitIndexFromHDF5Attr ---
  subroutine ReadLitIndexFromHDF5Attr ( itemID, name, index )
    use MoreTree, only: GetLitIndexFromString
    ! Dummy arguments
    integer, intent(in) :: ITEMID       ! Group etc. to make attr. for
    character(len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(out) :: INDEX        ! String index
    ! Local variables
    character (len=1024) :: LINE
    ! Executable code
    call GetHDF5Attribute ( itemID, name, line )
    if ( len_trim ( line ) > 0 ) then
      index = GetLitIndexFromString ( trim(line) )
    else
      index = 0
    end if
  end subroutine ReadLitIndexFromHDF5Attr

  ! ------------------------------------------ ReadStringIndexFromHDF5Attr ---
  subroutine ReadStringIndexFromHDF5Attr ( itemID, name, index )
    use MoreTree, only: GetStringIndexFromString
    ! Dummy arguments
    integer, intent(in) :: ITEMID       ! Group etc. to make attr. for
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

  ! --------------------------------------- WriteStringIndexAsHDF5Attribute
  subroutine WriteLitIndexAsHDF5Attribute ( itemID, name, index )
    use Intrinsic, only: LIT_INDICES
    integer, intent(in) :: ITEMID       ! Group etc. to make attr. for
    character(len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(in) :: INDEX        ! String index
    ! Executable code
    if ( index == 0 ) then
      call MakeHDF5Attribute ( itemID, name, '' )
    else
      call WriteStringIndexAsHDF5Attribute ( itemID, name, lit_indices ( index ) )
    end if
  end subroutine WriteLitIndexAsHDF5Attribute

  ! --------------------------------------- WriteStringIndexAsHDF5Attribute
  subroutine WriteStringIndexAsHDF5Attribute ( itemID, name, index )
    use String_table, only: GET_STRING
    integer, intent(in) :: ITEMID       ! Group etc. to make attr. for
    character(len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(in) :: INDEX        ! String index
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

! ======================= Private Procedures =========================  
! --------------------------------------------- AreThe2TypesEqual ---
  logical function AreThe2TypesEqual ( type1, type2 )
    ! This routine returns true if the two datatypes are "the same"
    integer, intent(in) :: type1
    integer, intent(in) :: type2
    ! Local variables
    integer :: status
    ! Initialize in case something goes wrong with hdf5 routine
    AreThe2TypesEqual = .false.
    call h5tEqual_f(type1, type2, AreThe2TypesEqual, status)
    if (status /= 0 ) AreThe2TypesEqual = .false.
  end function AreThe2TypesEqual

  subroutine check_for_fit(spaceID, value_dims, name)
  ! Checks that dataspace will fit into values before LoadFromHDF5DS
    integer, intent(in)               :: spaceID 
    integer, dimension(:), intent(in) :: value_dims
    character (len=*), intent(in)     :: NAME ! Name for this dataset
    integer                           :: rank
    integer                           :: value_rank
    integer(hsize_t), dimension(7)    :: dims, maxdims
    integer                           :: status  
    integer                           :: i  
    value_rank = size(value_dims)
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( status /= 0 )  call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get rank for dataset '//trim(name) )
    if ( rank /= value_rank ) call my_message ( MLSMSG_Error, ModuleName, &
      & 'Inconsistant rank for dataset '//trim(name) , &
      & 'rank(space),rank(values)', (/rank, value_rank/) )
    call h5sget_simple_extent_dims_f ( spaceID, dims(1:rank), maxdims(1:rank), &
      &  status )
    if ( status /= rank ) call my_message ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dimension information for dataset '//trim(name) , &
      & 'rank(space),h5s status', (/rank, status/) )
    if ( any ( dims(1:rank) > value_dims ) ) &
      & call my_message ( MLSMSG_Error, ModuleName, &
      & 'Dataspace too large for destination value of '//trim(name) , &
      & 'dims(space), dims(value)', (/ (int(dims(i)), value_dims(i), i=1, rank) /), &
      & no_pairs=.true. )
  end subroutine check_for_fit


  subroutine dump_space(spaceID)
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
    if ( DEEBUG ) print *, 'dims ', dims(1:rank)
    if ( DEEBUG ) print *, 'maxdims ', maxdims(1:rank)
    ! call h5soffset_simple_f ( spaceID, offset(1:rank), &
    !  &  status )
    ! print *, 'offset ', offset(1:rank)
  end subroutine dump_space

  subroutine mls_extend ( setID, newCount, start, dataSpaceID )
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
    endif
    if ( DEEBUG ) print *, 'spaceID ', spaceID
    call h5sis_simple_f ( spaceID, is_simple, status )
    if ( DEEBUG ) print *, 'is simple? ', is_simple
    call h5sget_simple_extent_ndims_f(spaceID, rank, status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get data space rank to extend data set' )
    if ( DEEBUG ) print *, 'rank ', rank
    call h5sget_simple_extent_dims_f(spaceID, dims,&
         maxdims, status)
    if ( status /= 0 .and. cantGetDataspaceDims <= MAXNUMWARNS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Unable to get data space dims to extend data set' )
      cantGetDataspaceDims = cantGetDataspaceDims + 1
      if ( cantGetDataspaceDims > MAXNUMWARNS ) &
        & call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Max no. of warnings reached--suppressing further ones')
    endif
    if ( DEEBUG ) print *, 'dims ', dims(1:rank)
    if ( DEEBUG ) print *, 'maxdims ', maxdims(1:rank)
    if ( DEEBUG ) print *, 'status ', status
    itFits = .true.
    do i = 1, min(rank, size(newCount))
      if ( my_start(i) + newCount(i) > dims(i) ) then
        itFits = .false.
        dims(i) = my_start(i) + newCount(i)
      endif
    enddo
    if ( itFits ) return
    if ( DEEBUG ) print *, 'Need to extend dataset'
    if ( DEEBUG ) print *, '(New dims) ', dims
    call h5dextend_f(setID, dims, status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to extend data set in mls_extend' )
    if ( present(dataspaceID) ) then
      call h5dget_space_f(setID, dataspaceID, status)
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to get return dataspaceID in mls_extend' )
    endif
  end subroutine mls_extend

  subroutine mls_hyperslab ( spaceID, value_dims, name, memspaceID, &
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
      & 'Impossible optional parameters pattern for dataset '//trim(name), &
      & 'status', (/status/) )
    value_rank = size(value_dims)
    if ( DEEBUG ) print *, 'value_rank: ', value_rank
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( DEEBUG ) print *, 'dataspace rank: ', rank
    if ( status /= 0 )  call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get rank for dataset '//trim(name) )
    call h5sget_simple_extent_dims_f ( spaceID, dims(1:rank), maxdims(1:rank), &
      &  status )
    if ( DEEBUG ) print *, 'dataspace dims: ', dims(1:rank)
    if ( status /= rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dimension information for dataset '//trim(name) )
    call h5screate_simple_f(value_rank, &
      & int(value_dims(1:value_rank), hsize_t), memspaceID, &
      & status)
    if ( DEEBUG ) print *, 'memspace id: ', memspaceID
    if (status /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
      & 'Unable to create memspace for dataset '//trim(name) )
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
    endif
    if (status /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
      & 'Unable to set hyperslab for dataset '//trim(name) )
    if ( DEEBUG ) print *, 'Returning memspaceID ', memspaceID
  end subroutine mls_hyperslab

  subroutine mls_hyperslab_save ( spaceID, &
    & start, count, stride, block )
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
      & 'Impossible optional parameters pattern for dataset '//trim(name), &
      & 'status', (/status/) )
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( DEEBUG ) print *, 'dataspace rank: ', rank
    if ( status /= 0 )  call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get rank for dataset '//trim(name) )
    call h5sget_simple_extent_dims_f ( spaceID, dims(1:rank), maxdims(1:rank), &
      &  status )
    if ( DEEBUG ) print *, 'dataspace dims: ', dims(1:rank)
    if ( status /= rank ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dimension information for dataset '//trim(name) )
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
    endif
    if (status /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
      & 'Unable to set hyperslab for dataset '//trim(name) )
  end subroutine mls_hyperslab_save

  subroutine my_message(severity, ModuleNameIn, Message, &
    & names, ints, reals, doubles, no_pairs)
    ! Take opportunity to dump a diagnostic table of values before stopping
    ! Dummy arguments
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text
    character (len=*), intent(in) :: names   ! comma-separated list of names
    integer, dimension(:), optional, intent(in)          :: ints
    real, dimension(:), optional, intent(in)             :: reals
    real(r8), dimension(:), optional, intent(in)         :: doubles
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
        call dump(ints, names, &
          & clean=clean, format=int_format, width=width)
      elseif ( present(reals) ) then
        call dump(reals, names, clean=clean)
      elseif ( present(doubles) ) then
        call dump(doubles, names, clean=clean)
      endif
    else
      if ( present(ints) ) then
        call dump_name_v_pairs(ints, names, &
          & clean=clean, format=int_format, width=width)
      elseif ( present(reals) ) then
        call dump_name_v_pairs(reals, names, &
          & clean=clean, format=real_format, width=width)
      elseif ( present(doubles) ) then
        call dump_name_v_pairs(doubles, names, &
          & clean=clean, format=dbl_format, width=width)
      endif
    endif
    call MLSMessage ( severity, ModuleNameIn, &
      & message )
  end subroutine my_message

  function name_to_attrID(fileID, dataName, attrName, dont_close) &
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
      dataID = -1
    else
      call h5aOpen_name_f ( dataID, attrName, attrID, status )
      if ( status /= 0 ) then
        attrID = -1
      else
        if (.not. my_dont_close) call h5aClose_f ( attrID, status )
      end if
      if (.not. my_dont_close) call h5dClose_f ( dataID, status )
    end if
  end function name_to_attrID

  function name_to_dataID(fileID, name, dont_close) result(dataID)
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
      if (.not. my_dont_close) call h5dClose_f ( dataID, status )
    end if
  end function name_to_dataID

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSHDF5

! $Log$
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
 
