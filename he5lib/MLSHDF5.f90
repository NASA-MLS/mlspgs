! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MLSHDF5

  ! This module contains MLS specific routines to do lowish level common HDF5
  ! tasks.  Initially mainly to do with simple attributes, but more will no
  ! doubt be added.

  ! Lets break down our use, parameters first
  use HDF5, only: H5S_SCALAR_F, H5T_NATIVE_INTEGER, H5T_NATIVE_CHARACTER, &
    & H5T_NATIVE_DOUBLE, HID_T, HSIZE_T
  ! Now routines
  use HDF5, only: H5ACREATE_F, H5AGET_TYPE_F, H5AOPEN_NAME_F, H5AREAD_F, &
    & H5AWRITE_F, H5ACLOSE_F, &
    & H5DCREATE_F, H5DGET_SPACE_F, H5DOPEN_F, &
    & H5DREAD_F, H5DWRITE_F, H5DCLOSE_F, &
    & H5SCREATE_F, H5SCREATE_SIMPLE_F, H5SCLOSE_F, &
    & H5SGET_SIMPLE_EXTENT_NDIMS_F, H5SGET_SIMPLE_EXTENT_DIMS_F, &
    & H5TCLOSE_F, H5TCOPY_F, H5TGET_SIZE_F, H5TSET_SIZE_F
  use MLSCommon, only: r4, r8
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR

  implicit NONE
  private

  public :: MakeHDF5Attribute, SaveAsHDF5DS, GetHDF5Attribute, LoadFromHDF5DS

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  !---------------------------------------------------------------------------

  interface MakeHDF5Attribute
    module procedure MakeHDF5Attribute_int, MakeHDF5Attribute_logical, &
      & MakeHDF5Attribute_string
  end interface

  interface GetHDF5Attribute
    module procedure GetHDF5Attribute_int, GetHDF5Attribute_logical, &
      & GetHDF5Attribute_string
  end interface

  interface SaveAsHDF5DS
    module procedure SaveAsHDF5DS_intarr1, &
      SaveAsHDF5DS_dblarr1, SaveAsHDF5DS_dblarr2
  end interface

  interface LoadFromHDF5DS
    module procedure LoadFromHDF5DS_intarr1, &
      LoadFromHDF5DS_dblarr1, LoadFromHDF5DS_dblarr2
  end interface

  ! Local parameters
  integer, dimension(7) :: ones = (/1,1,1,1,1,1,1/)

contains ! ======================= Public Procedures =========================

  ! ------------------------------------- MakeHDF5Attribute_int
  subroutine MakeHDF5Attribute_int ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    integer, intent(in) :: VALUE        ! Value of attribut

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: DSID                     ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5

    ! Executable code
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
  subroutine MakeHDF5Attribute_logical ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    logical, intent(in) :: VALUE        ! Value of attribute

    ! Local variables
    integer :: IVALUE                   ! Value as integer

    ! Executable code
    iValue = 0
    if ( value ) iValue = 1
    call MakeHDF5Attribute ( itemID, name, iValue )
  end subroutine MakeHDF5Attribute_logical

  ! ------------------------------------- MakeHDF5Attribute_string
  subroutine MakeHDF5Attribute_string ( itemID, name, value )
    integer, intent(in) :: ITEMID       ! Group etc. to make attr to.
    character (len=*), intent(in) :: NAME ! Name of attribute
    character (len=*), intent(in) :: VALUE ! Value of attribute

    ! Local variables
    integer :: ATTRID                   ! ID for attribute
    integer :: DSID                     ! ID for dataspace
    integer :: STATUS                   ! Flag from HDF5
    integer :: STRINGTYPE               ! Type for string

    ! Executable code
    ! Setup
    ! Create a data type for this string
    call h5tcopy_f( H5T_NATIVE_CHARACTER, stringtype, status ) 
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to create stringtype '//trim(name) )
    call h5tset_size_f(stringtype, len_trim(value), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to set size for stringtype '//trim(name) )
    ! Create dataspace and attribute
    call h5sCreate_F ( h5s_scalar_f, dsID, status )
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
  end subroutine MakeHDF5Attribute_string

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
    call h5aread_f ( attrID, H5T_NATIVE_INTEGER, value, ones, status )
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

  ! ----------------------------------- LoadFromHDF5DS_intarr1
  subroutine LoadFromHDF5DS_intarr1 ( locID, name, value )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    integer, intent(out) :: VALUE(:)    ! The array itself

    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: RANK                     ! Rank in file
    integer (kind=hSize_t) , dimension(maxDimensions) :: DIMS, MAXDIMS ! Dimensions in file

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    ! Check that dimensions are all ok.
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get rank for dataset '//trim(name) )
    if ( rank /= 1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Inconsistant rank for dataset '//trim(name) )
    call h5sget_simple_extent_dims_f ( spaceID, dims(1:rank), maxdims(1:rank), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dimension information for dataset '//trim(name) )
    if ( any ( dims(1:rank) > shape(value) ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Dataspace too large for destination value of '//trim(name) )
    ! Now, (at last!) read the data
    call h5dread_f ( setID, H5T_NATIVE_INTEGER, value, &
      & (/ shape(value), ones(1:6) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_intarr1

  ! ----------------------------------- LoadFromHDF5DS_dblarr1
  subroutine LoadFromHDF5DS_dblarr1 ( locID, name, value )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r8), intent(out) :: VALUE(:)    ! The array itself

    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(1) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: RANK                     ! Rank in file
    integer (kind=hSize_t) , dimension(maxDimensions) :: DIMS, MAXDIMS ! Dimensions in file

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    ! Check that dimensions are all ok.
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get rank for dataset '//trim(name) )
    if ( rank /= 1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Inconsistant rank for dataset '//trim(name) )
    call h5sget_simple_extent_dims_f ( spaceID, dims(1:rank), maxdims(1:rank), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dimension information for dataset '//trim(name) )
    if ( any ( dims(1:rank) > shape(value) ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Dataspace too large for destination value of '//trim(name) )
    ! Now, (at last!) read the data
    call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & (/ shape(value), ones(1:6) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_dblarr1

  ! ----------------------------------- LoadFromHDF5DS_dblarr2
  subroutine LoadFromHDF5DS_dblarr2 ( locID, name, value )
    ! This routine loads a predefined array with values from a DS
    integer, intent(in) :: LOCID        ! Where to place it (group/file)
    character (len=*), intent(in) :: NAME ! Name for this dataset
    real(r8), intent(out) :: VALUE(:,:) ! The array itself

    ! Local parameters
    integer, parameter :: MAXDIMENSIONS = 7

    ! Local variables
    integer :: STATUS                   ! Flag from HDF5
    integer, dimension(2) :: SHP        ! Shape of value
    integer :: SPACEID                  ! ID of dataspace
    integer :: SETID                    ! ID of dataset
    integer :: RANK                     ! Rank in file
    integer (kind=hSize_t) , dimension(maxDimensions) :: DIMS, MAXDIMS ! Dimensions in file

    ! Executable code
    shp = shape ( value )
    call h5dOpen_f ( locID, name, setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataset '//trim(name) )
    call h5dget_space_f ( setID, spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to open dataspace for dataset '//trim(name) )
    ! Check that dimensions are all ok.
    call h5sget_simple_extent_ndims_f ( spaceID, rank, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get rank for dataset '//trim(name) )
    if ( rank /= 2 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Inconsistant rank for dataset '//trim(name) )
    call h5sget_simple_extent_dims_f ( spaceID, dims(1:rank), maxdims(1:rank), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to get dimension information for dataset '//trim(name) )
    if ( any ( dims(1:rank) > shape(value) ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Dataspace too large for destination value of '//trim(name) )
    ! Now, (at last!) read the data
    call h5dread_f ( setID, H5T_NATIVE_DOUBLE, value, &
      & (/ shape(value), ones(1:5) /), status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to read dataset '//trim(name) )
    ! Close up
    call h5sClose_f ( spaceID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataspace for dataset '//trim(name) )
    call h5dClose_f ( setID, status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to close dataset '//trim(name) )
  end subroutine LoadFromHDF5DS_dblarr2


end module MLSHDF5

! $Log$
