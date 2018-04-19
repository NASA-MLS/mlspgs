! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MLSAuxData

  ! Reading and interacting with Level 1B data (HDF5)
  use MLS_DataProducts, only: DATAPRODUCTS_T
  use HDF5, only: HID_T, HSIZE_T, H5T_NATIVE_CHARACTER, H5T_NATIVE_DOUBLE, &
       H5T_NATIVE_INTEGER, H5T_IEEE_F32LE, H5S_UNLIMITED_F, H5T_STD_I32LE, & 
       H5T_NATIVE_REAL,H5T_IEEE_F64LE,H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, &
       H5DOPEN_F, H5DCLOSE_F, H5SCREATE_SIMPLE_F, H5PCREATE_F, H5DCREATE_F, &
       H5PSET_CHUNK_F, H5SCLOSE_F, H5DGET_SPACE_F, H5PCLOSE_F, & 
       H5SGET_SIMPLE_EXTENT_NDIMS_F, H5SGET_SIMPLE_EXTENT_DIMS_F, &
       H5DGET_CREATE_PLIST_F, H5PGET_CHUNK_F, H5DGET_TYPE_F, & 
       H5SSELECT_HYPERSLAB_F, H5DREAD_F, H5DWRITE_F, H5DEXTEND_F, &
       H5ACREATE_F, H5AWRITE_F, H5AREAD_F, H5ACLOSE_F, H5TCOPY_F, &
       H5TSET_SIZE_F, H5AOPEN_NAME_F, H5AGET_SPACE_F, &
       H5TEQUAL_F, H5ESET_AUTO_F, H5GCREATE_F, H5GCLOSE_F, &
       H5GOPEN_F, H5PSET_FILL_VALUE_F, SIZE_T
  use MLSCommon, only: NAMELEN
  use MLSKINDS, only: R4, R8
  use MLSStrings, only: LOWERCASE
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
  use output_m, only: switchOutput, revertOutput
  
  implicit NONE

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!     (data types and parameters)
! MLSAuxData_T                   Quantities from an L1B/L2 data file
!
!                 (subroutines and functions)
!
! AddMLSAuxDataToDatabase      Creates and updates an array of MLSAuxData_T.
!                              Returns size of the array.
!
! DestroyMLSAuxDataDatabase    Deallocates each element of an array of 
!                              MLSAuxData_T and then deallocates the array.
!
! Build_MLSAuxData             Incorporates the routines stated below to 
!                              create a dataset of any type and any dimension
!                              so that user does not need to worry about 
!                              MLSAuxData_T data structure.
!
! Recall_MLSAuxData            Incorporates the routines stated below to 
!                              read a dataset of any type and any dimension
!                              so that user does not need to worry about 
!                              MLSAuxData_T data structure.
!
! Allocate_MLSAuxData          Allocates a MLSAuxData data structure.
! 
! Create_MLSAuxData            Creates MLSAuxData in a file. 
!
! CreateGroup_MLSAuxData       Creates a group or subgroup in a file. 
!
! Read_MLSAuxAttributes        Reads an attribute of a MLSAuxData quantity
!                              from a file.
!
! Read_MLSAuxData              Reads all info concerning a MLSAuxData quantity
!                              from a file.
!
! Write_MLSAuxAttributes       Writes an attribute of a MLSAuxData quantity
!                              to a file.
!
! Write_MLSAuxData             Writes all info concerning a MLSAuxData quantity
!                              to a file.
!
! Deallocate_MLSAuxData        Deallocates memory of a MLSAuxData_T.
!
! === (end of toc) ===
  private
  public :: MLSAuxData_T, Create_MLSAuxData, Read_MLSAuxData, & 
   & Deallocate_MLSAuxData, Write_MLSAuxData, Read_MLSAuxAttributes, &
   & Write_MLSAuxAttributes, CreateGroup_MLSAuxData, Allocate_MLSAuxData, &
   & Build_MLSAuxData, AddMLSAuxDataToDatabase, DestroyMLSAuxDataDatabase, &
   MaxCharFieldLen
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  !
  ! Parameters
  INTEGER, PARAMETER :: MaxCharFieldLen = 100000  ! max char field length
  !
  type MLSAuxData_T
    character (len=namelen) :: name, type_name  
    character(len=20),  dimension(:), pointer :: Dimensions => NULL()
    !
    character(len=MaxCharFieldLen), dimension(:,:,:), pointer :: CharField => &
         NULL()
    real(r8),  dimension(:,:,:), pointer :: DpField          => NULL()   
    real(r4),  dimension(:,:,:), pointer :: RealField        => NULL()   
    integer,   dimension(:,:,:), pointer :: IntField         => NULL()
    ! all the above dimensioned (noAuxInds,maxMIFs,noMAFs)
    !
    real(r8), dimension(:), pointer :: FrequencyCoordinates  => NULL()
    real(r8), dimension(:), pointer :: VerticalCoordinates   => NULL()
    real(r8), dimension(:), pointer :: HorizontalCoordinates => NULL()
    !
    integer :: FrequencyDimension          ! Enumerated type
    integer :: VerticalDimension           ! Enumerated type  
    integer :: HorizontalDimension         ! Enumerated type
    !
    integer :: rank ! rank of data set
    !
  end type MLSAuxData_T
!------------------------------------------------------------------- Exceptions
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_OPEN   = 'HDF5 Error Opening Dataset '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_CREATE = 'HDF5 Error Creating Dataset '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_WRITE  = 'HDF5 Error Writing Dataset '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_READ  = 'HDF5 Error Reading Dataset '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_EXTEND  = 'HDF5 Error Extending Dataset '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_GETSPACE = 'HDF5 Error Retrieving Dataspace of Dataset '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_GET_TYPE = 'HDF5 Error Retrieving Data Type of Dataset '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_CLOSE  = 'HDF5 Error Closing Dataset '

  CHARACTER(len=*), PUBLIC, PARAMETER :: &
    H5_ERROR_GROUP_OPEN   = 'HDF5 Error Opening Group '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_GROUP_CREATE = 'HDF5 Error Creating Group '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_GROUP_CLOSE  = 'HDF5 Error Closing Group '

  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_ATT_WRITE  = 'HDF5 Error Writing Attribute '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_ATT_GETSPACE  = 'HDF5 Error Getting Space of Attribute '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_ATT_READ  = 'HDF5 Error Reading Attribute '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_ATT_CREATE = 'HDF5 Error Creating Attribute '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_ATT_CLOSE  = 'HDF5 Error Closing Attribute '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_ATT_OPEN  = 'HDF5 Error Opening Attribute '

  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_TYPE_COPY = 'HDF5 Error Copying Datatype '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_TYPE_SET  = 'HDF5 Error Setting Datatype '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_TYPE_EQUAL  = 'HDF5 Error Equating Datatype '

  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSPACE_OPEN   = 'HDF5 Error Opening Dataspace '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSPACE_CREATE = 'HDF5 Error Creating Dataspace '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSPACE_GET_DIMS = & 
    'HDF5 Error Obtaining Dimensions of Dataspace '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSPACE_GET_NDIMS = & 
    'HDF5 Error Obtaining Number of Dimensions of Dataspace '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSPACE_HYPERSLAB = 'HDF5 Error Selecting Hyperslab '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSPACE_CLOSE  = 'HDF5 Error Closing Dataspace '

  CHARACTER(len=*), PUBLIC, PARAMETER :: &
    H5_ERROR_FILE_OPEN   = 'HDF5 Error Opening File '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_FILE_CLOSE  = 'HDF5 Error Closing File '

  CHARACTER(len=*), PUBLIC, PARAMETER :: &
    H5_ERROR_OPEN   = 'HDF5 Error Opening FORTRAN File API '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_CLOSE  = 'HDF5 Error Closing FORTRAN File API '

  CHARACTER(len=*), PUBLIC, PARAMETER :: &
    H5_ERROR_PROPERTY_CREATE = 'HDF5 Error Creating Property List '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_PROPERTY_CHUNK_SET  = 'HDF5 Error Setting Chunk Property '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_PROPERTY_CHUNK_GET  = 'HDF5 Error Getting Chunk Property '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_PROPERTY_CLOSE  = 'HDF5 Error Closing Property List '

  logical, parameter :: DEBUG = .false.

  interface Build_MLSAuxData
    module procedure Build_MLSAuxData_Character
    module procedure Build_MLSAuxData_Double, &
         Build_MLSAuxData_Double_1d, Build_MLSAuxData_Double_2d
    module procedure Build_MLSAuxData_Double_3d
    module procedure Build_MLSAuxData_Real, Build_MLSAuxData_Real_1d, & 
         Build_MLSAuxData_Real_2d
    module procedure Build_MLSAuxData_Real_3d
    module procedure Build_MLSAuxData_Integer, Build_MLSAuxData_Integer_1d, & 
         Build_MLSAuxData_Integer_2d
    module procedure Build_MLSAuxData_Integer_3d
  end interface

  interface Recall_MLSAuxData
    module procedure Recall_MLSAuxData_Character
    module procedure Recall_MLSAuxData_Double, &
         Recall_MLSAuxData_Double_1d, Recall_MLSAuxData_Double_2d
    module procedure Recall_MLSAuxData_Double_3d
    module procedure Recall_MLSAuxData_Real, Recall_MLSAuxData_Real_1d, & 
         Recall_MLSAuxData_Real_2d
    module procedure Recall_MLSAuxData_Real_3d
    module procedure Recall_MLSAuxData_Integer, Recall_MLSAuxData_Integer_1d, & 
         Recall_MLSAuxData_Integer_2d
    module procedure Recall_MLSAuxData_Integer_3d
  end interface

contains ! ============================ MODULE PROCEDURES ====================
!-------------------------------------------------------AddMLSAuxDataToDatabase
  integer function AddMLSAuxDataToDatabase ( database, item )
  ! Add a MLSAuxData to a database, and/or create the database.
  ! Returns the size of the array.
    
    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    type (MLSAuxData_T), dimension(:), pointer :: database
    type (MLSAuxData_T), intent(in) :: item

    ! Local variables
    type (MLSAuxData_T), dimension(:), pointer :: tempDatabase => NULL()

    include "addItemToDatabase.f9h"

    AddMLSAuxDataToDatabase = newSize

  end function AddMLSAuxDataToDatabase
!-----------------------------------------------------DestroyMLSAuxDataDatabase
  subroutine DestroyMLSAuxDataDatabase ( database )
  ! Deallocates the elements of array of MLSAuxData_T and then the array.
    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type (MLSAuxData_T), dimension(:), pointer :: database
    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: index, s, status

    if ( associated(database) ) then
       do index = 1, SIZE(database)
          call  Deallocate_MLSAuxData( database(index) )
       end do
       s = size(database) * storage_size(database) / 8
       addr = 0
       if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
       deallocate ( database, stat=status )
       call test_deallocate ( status, ModuleName, "database", s, address=addr )
    end if

  end subroutine DestroyMLSAuxDataDatabase
!---------------------------------------------------------- Allocate_MLSAuxData
 subroutine Allocate_MLSAuxData( name, data_type, dims, MLSData  )
    ! This should be called when allocating a MLSAuxData structure.
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( MLSAuxData_T ), intent(inout) :: MLSData
    integer, dimension(3), intent(in) :: dims
    character (len=*), intent(in) :: name, data_type
    ! internal variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: status

    MLSData%name        = trim(name) 
    MLSData%type_name   = trim(data_type)

    if ( data_type == 'real') then
      allocate(MLSData%RealField(dims(1),dims(2),dims(3)), stat=status)
      addr = 0
      if ( status == 0 ) then
        addr = transfer(c_loc(MLSData%RealField(1,1,1)), addr)
      end if
      call test_allocate ( status, ModuleName, "MLSData%RealField", &
        & uBounds = dims(1:3), elementSize = storage_size(MLSData%RealField) / 8, &
        & address=addr )
    end if

    if ( data_type == 'double') then
      allocate(MLSData%DpField(dims(1),dims(2),dims(3)), stat=status)
      addr = 0
      if ( status == 0 ) then
        addr = transfer(c_loc(MLSData%DpField(1,1,1)), addr)
      end if
      call test_allocate ( status, ModuleName, "MLSData%DpField", &
        & uBounds = dims(1:3), elementSize = storage_size(MLSData%DpField) / 8, &
        & address=addr )
    end if

    if ( data_type == 'integer') then
      allocate(MLSData%IntField(dims(1),dims(2),dims(3)), stat=status)
      addr = 0
      if ( status == 0 ) then
        addr = transfer(c_loc(MLSData%IntField(1,1,1)), addr)
      end if
      call test_allocate ( status, ModuleName, "MLSData%IntField", &
        & uBounds = dims(1:3), elementSize = storage_size(MLSData%IntField) / 8, &
        & address=addr )
    end if

    if ( data_type == 'character') then
      allocate(MLSData%CharField(dims(1),dims(2),dims(3)), stat=status)
      addr = 0
      if ( status == 0 ) then
        addr = transfer(c_loc(MLSData%CharField(1,1,1)), addr)
      end if
      call test_allocate ( status, ModuleName, "MLSData%CharField", &
        & uBounds = dims(1:3), elementSize = storage_size(MLSData%CharField) / 8, &
        & address=addr )
    end if

 end subroutine Allocate_MLSAuxData
!---------------------------------------------------------Deallocate_MLSAuxData
 subroutine Deallocate_MLSAuxData( MLSAuxData )
    ! This should be called when deallocating a MLSAuxData structure.
    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( MLSAuxData_T ), intent(inout) :: MLSAuxData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: s, status

    if (associated(MLSAuxData%RealField)) then
      s = size(MLSAuxData%RealField) * storage_size(MLSAuxData%RealField) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MLSAuxData%RealField(1,1,1)), addr)
      deallocate(MLSAuxData%RealField, stat=status)
      call test_deallocate ( status, ModuleName, &
        & ' MLSAuxData%RealField in ' // trim(MLSAuxData%name), s, address=addr )
    end if

    if (associated(MLSAuxData%IntField)) then 
      s = size(MLSAuxData%IntField) * storage_size(MLSAuxData%IntField) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MLSAuxData%IntField(1,1,1)), addr)
      deallocate(MLSAuxData%IntField, stat=status)
      call test_deallocate ( status, ModuleName, &
        & ' MLSAuxData%IntField in ' // trim(MLSAuxData%name), s, address=addr )
    end if

    if (associated(MLSAuxData%CharField)) then 
      s = size(MLSAuxData%CharField) * storage_size(MLSAuxData%CharField) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MLSAuxData%CharField(1,1,1)), addr)
      deallocate(MLSAuxData%CharField, stat=status)
      call test_deallocate ( status, ModuleName, &
        & ' MLSAuxData%CharField in ' // trim(MLSAuxData%name), s, address=addr )
    end if

    if (associated(MLSAuxData%DpField)) then 
      s = size(MLSAuxData%DpField) * storage_size(MLSAuxData%DpField) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MLSAuxData%DpField(1,1,1)), addr)
      deallocate(MLSAuxData%DpField, stat=status)
      call test_deallocate ( status, ModuleName, &
        & ' MLSAuxData%DpField in ' // trim(MLSAuxData%name), s, address=addr )
    end if

    if (associated(MLSAuxData%FrequencyCoordinates)) then 
      s = size(MLSAuxData%FrequencyCoordinates) * &
        & storage_size(MLSAuxData%FrequencyCoordinates) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MLSAuxData%FrequencyCoordinates(1)), addr)
      deallocate(MLSAuxData%FrequencyCoordinates, stat=status)
      call test_deallocate ( status, ModuleName, &
        & ' MLSAuxData%FrequencyCoordinates in ' // trim(MLSAuxData%name), s, &
        & address=addr )
    end if

    if (associated(MLSAuxData%VerticalCoordinates)) then 
      s = size(MLSAuxData%VerticalCoordinates) * &
        & storage_size(MLSAuxData%VerticalCoordinates) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MLSAuxData%VerticalCoordinates(1)), addr)
      deallocate(MLSAuxData%VerticalCoordinates, stat=status)
      call test_deallocate ( status, ModuleName, &
        & ' MLSAuxData%VerticalCoordinates in ' // trim(MLSAuxData%name), s, &
        & address=addr )
    end if

    if (associated(MLSAuxData%HorizontalCoordinates)) then 
      s = size(MLSAuxData%HorizontalCoordinates) * &
        & storage_size(MLSAuxData%HorizontalCoordinates) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MLSAuxData%HorizontalCoordinates(1)), addr)
      deallocate(MLSAuxData%HorizontalCoordinates, stat=status)
      call test_deallocate ( status, ModuleName, &
        & ' MLSAuxData%HorizontalCoordinates in ' // trim(MLSAuxData%name), s, &
        & address=addr )
    end if

    if (associated(MLSAuxData%Dimensions)) then 
      s = size(MLSAuxData%Dimensions) * &
        & storage_size(MLSAuxData%Dimensions) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(MLSAuxData%Dimensions(1)), addr)
      deallocate(MLSAuxData%Dimensions, stat=status)
      call test_deallocate ( status, ModuleName, &
        & ' MLSAuxData%Dimensions in ' // trim(MLSAuxData%name), s, address=addr )
    end if

 end subroutine Deallocate_MLSAuxData
!--------------------------------------------------------------Build_MLSAuxData
 subroutine Build_MLSAuxData_Character( file_id, dataset, char_data, & 
      char_length, lastIndex, disable_attrib)
    type( DataProducts_T ), intent(in) :: dataset
    character (len=*), intent(in) :: char_data
    integer, intent(in) :: char_length
    integer, intent(in), optional ::lastIndex
    logical, intent(in), optional :: disable_attrib
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error
    logical :: attribenabled

    dims = 1

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),&
         trim(dataset%data_type),dims,MLSData)

    MLSData%CharField(1,1,1) = char_data
    MLSData%rank = 1

    call CopyFromDataProducts (dataset, MLSData)

    attribenabled = .false.
    if (present(lastIndex)) then 
       if (lastIndex == 1) then 
          if (present (disable_attrib)) then
             attribenabled = .false.
          else
             attribenabled = .true.
          end if
       end if
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=attribenabled, & 
            string_length=char_length, index=lastIndex)
    else
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=.true., & 
            string_length=char_length)
    end if

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Character
!--------------------------------------------------------------Build_MLSAuxData
 subroutine Build_MLSAuxData_Integer( file_id, dataset, int_data, & 
      lastIndex, fill_value, disable_attrib)
    type( DataProducts_T ), intent(in) :: dataset
    integer, intent(in) :: int_data
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id
    integer, intent(in), optional :: fill_value
    logical, intent(in), optional :: disable_attrib

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error
    logical :: attribenabled

    dims = 1

    call deallocate_mlsauxdata(MLSData)

    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dims,MLSData)

    MLSData%IntField(1,1,1) = int_data
    MLSData%rank = 1

    call CopyFromDataProducts (dataset, MLSData)

    attribenabled = .false.
    if ( present (lastIndex)) then
       if (lastIndex == 1) then 
          if (present (disable_attrib)) then
             attribenabled = .false.
          else
             attribenabled = .true.
          end if
       end if
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=attribenabled, index=lastIndex, &
               fill_value_i=fill_value)
    else

       call Write_MLSAuxData(file_id, MLSData, error, &
            write_attributes=.true., fill_value_i=fill_value)

    end if

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) )
    call deallocate_mlsauxdata(MLSData)

 end subroutine Build_MLSAuxData_Integer
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Real (file_id, dataset, real_data, lastIndex, &
      fill_value, disable_attrib)
    type( DataProducts_T ), intent(inout) :: dataset
    real, intent(in) :: real_data
    integer, intent(in), optional :: lastIndex
    real, intent(in), optional :: fill_value
    integer(hid_t), intent(in) :: file_id
    logical, intent(in), optional :: disable_attrib

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error
    logical :: attribenabled

    dims = 1

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dims,MLSData)

    MLSData%RealField(1,1,1) = real_data
    MLSData%rank = 1

    call CopyFromDataProducts (dataset, MLSData)

    attribenabled = .false.
    if ( present(lastIndex) ) then
       if (lastIndex == 1) then 
          if (present (disable_attrib)) then
             attribenabled = .false.
          else
             attribenabled = .true.
          end if
       end if
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=attribenabled, index=lastIndex, &
               fill_value_r=fill_value)
    else
       call Write_MLSAuxData (file_id, MLSData, error, fill_value_r=fill_value)
    end if

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Real
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Double( file_id, dataset, double_data, lastIndex, &
      fill_value, disable_attrib)
    type( DataProducts_T ), intent(in) :: dataset
    real(r8), intent(in) :: double_data
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id
    real(r8), intent(in), optional :: fill_value
    logical, intent(in), optional :: disable_attrib

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error
    logical :: attribenabled

    dims = 1

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dims,MLSData)
    MLSData%DpField(1,1,1) = double_data
    MLSData%rank = 1

    call CopyFromDataProducts (dataset, MLSData)

    attribenabled = .false.
    if ( present (lastIndex) ) then 
        if (lastIndex == 1) then 
          if (present (disable_attrib)) then
             attribenabled = .false.
          else
             attribenabled = .true.
          end if
       end if
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=attribenabled, index=lastIndex, &
               fill_value_d=fill_value)
     else
       call Write_MLSAuxData(file_id, MLSData, error,&           
            write_attributes=.true., fill_value_d=fill_value)
    end if

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Double
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Real_1d( file_id, dataset, real_data, lastIndex, &
      dims, fill_value, disable_attrib)
    type( DataProducts_T ), intent(in) :: dataset
    real, dimension(:), intent(in) :: real_data
    integer, dimension(:), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    real, intent(in), optional :: fill_value
    integer(hid_t), intent(in) :: file_id
    logical, intent(in), optional :: disable_attrib

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error
    logical :: attribenabled

    dim_array = 1
    if (present(dims) ) then
       dim_array(:size(dims)) = dims
    else
       do i=1,size(shape(real_data))
          dim_array(i) = size(real_data,i)
       end do
    end if

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    MLSData%RealField(1:dim_array(1),1,1) = real_data(1:dim_array(1)) 

    if ( present(lastIndex)) then 
       MLSData%rank = 2
    else
       MLSData%rank = 1
    end if

    call CopyFromDataProducts (dataset, MLSData)

    attribenabled = .false.
    if ( present(lastIndex)) then
       if (lastIndex == 1) then 
          if (present (disable_attrib)) then
             attribenabled = .false.
          else
             attribenabled = .true.
          end if
       end if
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=attribenabled, index=lastIndex, &
               fill_value_r=fill_value)
    else
       call Write_MLSAuxData(file_id, MLSData, error,&           
            write_attributes=.true., fill_value_r=fill_value) 
    end if

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) ) 

    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Real_1d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Double_1d( file_id, dataset, double_data, & 
      lastIndex, dims, fill_value, disable_attrib)
    type( DataProducts_T ), intent(in) :: dataset
    real(r8), dimension(:), intent(in) :: double_data
    integer, dimension(:), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id
    real(r8), intent(in), optional :: fill_value
    logical, intent(in), optional :: disable_attrib

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error
    logical :: attribenabled

    dim_array = 1
    if (present(dims) ) then
       dim_array(:size(dims)) = dims
    else
       do i=1,size(shape(double_data))
          dim_array(i) = size(double_data,i)
       end do
    end if

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    MLSData%DpField(1:dim_array(1),1,1) = double_data(1:dim_array(1)) 

    if ( present (lastIndex) ) then 
       MLSData%rank = 2
    else
       MLSData%rank = 1
    end if

    call CopyFromDataProducts (dataset, MLSData)

    attribenabled = .false.
    if ( present(lastIndex)) then
       if (lastIndex == 1) then 
          if (present (disable_attrib)) then
             attribenabled = .false.
          else
             attribenabled = .true.
          end if
       end if
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=attribenabled, index=lastIndex, &
               fill_value_d=fill_value)
    else
       call Write_MLSAuxData(file_id, MLSData, error, &
            write_attributes=.true., fill_value_d=fill_value)
    end if

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Double_1d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Integer_1d( file_id, dataset, integer_data, &
      lastIndex, dims, fill_value, disable_attrib)
    type( DataProducts_T ), intent(in) :: dataset
    integer, dimension(:), intent(in) :: integer_data
    integer, dimension(:), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id
    integer, intent(in), optional :: fill_value
    logical, intent(in), optional :: disable_attrib

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error
    logical :: attribenabled

    dim_array = 1
    if (present(dims) ) then
       dim_array(:size(dims)) = dims
    else
       do i=1,size(shape(integer_data))
          dim_array(i) = size(integer_data,i)
       end do
    end if

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    MLSData%IntField(1:dim_array(1),1,1) = integer_data(1:dim_array(1)) 

    if (present (lastIndex) ) then
       MLSData%rank = 2 
    else
       MLSData%rank = 1
    end if

    call CopyFromDataProducts (dataset, MLSData)

    attribenabled = .false.
    if ( present(lastIndex)) then
       if (lastIndex == 1) then 
          if (present (disable_attrib)) then
             attribenabled = .false.
          else
             attribenabled = .true.
          end if
       end if
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=attribenabled, index=lastIndex, &
               fill_value_i=fill_value)
    else
       call Write_MLSAuxData(file_id, MLSData, error, &
            write_attributes=.true., fill_value_i=fill_value)
    end if

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Integer_1d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Real_2d( file_id, dataset, real_data, lastIndex, &
      dims, fill_value, disable_attrib)
    type( DataProducts_T ), intent(in) :: dataset
    real, dimension(:,:), intent(in) :: real_data
    integer, dimension(:), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    real, intent(in), optional :: fill_value
    integer(hid_t), intent(in) :: file_id
    logical, intent(in), optional :: disable_attrib

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error
    logical :: attribenabled

    dim_array = 1
    if (present(dims) ) then
       dim_array(:size(dims)) = dims
    else
       do i=1,size(shape(real_data))
          dim_array(i) = size(real_data,i)
       end do
    end if

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    MLSData%RealField(1:dim_array(1),1:dim_array(2),1) = &
      & real_data(1:dim_array(1),1:dim_array(2))

    if (present (lastIndex) ) then 
       MLSData%rank = 3
    else
       MLSData%rank = 2
    end if

    call CopyFromDataProducts (dataset, MLSData)

    if ( MLSData%name == "R1A:118.B22D:PT.S0.DACS-4 precision" .and. DEBUG ) then
      print *, '2d real dim_array: ', dim_array
    end if
    attribenabled = .false.
    if (present (lastIndex) ) then 
       if (lastIndex == 1) then 
          if (present (disable_attrib)) then
             attribenabled = .false.
          else
             attribenabled = .true.
          end if
       end if
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=attribenabled, index=lastIndex, &
               fill_value_r=fill_value)
    else
       call Write_MLSAuxData(file_id, MLSData, error, &
            write_attributes=.true., fill_value_r=fill_value)
    end if

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Real_2d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Double_2d( file_id, dataset, double_data, & 
      lastIndex, dims, fill_value, disable_attrib)
    type( DataProducts_T ), intent(in) :: dataset
    real(r8), dimension(:,:), intent(in) :: double_data
    integer, dimension(:), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id
    real(r8), intent(in), optional :: fill_value
    logical, intent(in), optional :: disable_attrib

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error
    logical :: attribenabled

    dim_array = 1
    if (present(dims) ) then
       dim_array(:size(dims)) = dims
    else
       do i=1,size(shape(double_data))
          dim_array(i) = size(double_data,i)
       end do
    end if

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    MLSData%DpField(1:dim_array(1),1:dim_array(2),1) = &
      & double_data(1:dim_array(1),1:dim_array(2))

    if (present (lastIndex) ) then 
       MLSData%rank = 3
    else
       MLSData%rank = 2
    end if

    call CopyFromDataProducts (dataset, MLSData)

    if ( MLSData%name == "R1A:118.B22D:PT.S0.DACS-4 precision"  .and. DEBUG ) then
      print *, '2d double dim_array: ', dim_array
    end if
    attribenabled = .false.
    if (present (lastIndex) ) then
       if (lastIndex == 1) then 
          if (present (disable_attrib)) then
             attribenabled = .false.
          else
             attribenabled = .true.
          end if
       end if
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=attribenabled, index=lastIndex, &
               fill_value_d=fill_value)
    else
       call Write_MLSAuxData(file_id, MLSData, error, &
            write_attributes=.true., fill_value_d=fill_value)
    end if

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Double_2d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Integer_2d( file_id, dataset, integer_data, &
      lastIndex, dims, fill_value, disable_attrib)
    type( DataProducts_T ), intent(in) :: dataset
    integer, dimension(:,:), intent(in) :: integer_data
    integer, dimension(:), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id
    integer, intent(in), optional :: fill_value
    logical, intent(in), optional :: disable_attrib

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error
    logical :: attribenabled

    dim_array = 1
    if (present(dims) ) then
       dim_array(:size(dims)) = dims
    else
       do i=1,size(shape(integer_data))
          dim_array(i) = size(integer_data,i)
       end do
    end if

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    MLSData%IntField(1:dim_array(1),1:dim_array(2),1) = &
      & integer_data(1:dim_array(1),1:dim_array(2))

    if (present (lastIndex)) then 
       MLSData%rank = 3
    else
       MLSData%rank = 2
    end if

    call CopyFromDataProducts (dataset, MLSData)

    attribenabled = .false.
    if (present (lastIndex)) then
       if (lastIndex == 1) then 
          if (present (disable_attrib)) then
             attribenabled = .false.
          else
             attribenabled = .true.
          end if
       end if
       call Write_MLSAuxData(file_id, MLSData, error, & 
            write_attributes=attribenabled, index=lastIndex, &
               fill_value_i=fill_value)
    else
       call Write_MLSAuxData(file_id, MLSData, error, &
            write_attributes=.true., fill_value_i=fill_value)
    end if

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Integer_2d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Real_3d( file_id, dataset, real_data, &
      dims, fill_value)
    type( DataProducts_T ), intent(in) :: dataset
    real, dimension(:,:,:), intent(in) :: real_data
    integer, dimension(:), intent(in), optional :: dims
    real, intent(in), optional :: fill_value
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error

    dim_array = 1
    if (present(dims) ) then
       dim_array(:size(dims)) = dims
    else
       do i=1,size(shape(real_data))
          dim_array(i) = size(real_data,i)
       end do
    end if

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    MLSData%rank = 3

    MLSData%RealField(1:dim_array(1),1:dim_array(2),1:dim_array(3)) = &
      & real_data(1:dim_array(1),1:dim_array(2),1:dim_array(3))

    call CopyFromDataProducts (dataset, MLSData)

    if ( MLSData%name == "R1A:118.B22D:PT.S0.DACS-4 precision"  .and. DEBUG ) then
      print *, '3d real dim_array: ', dim_array
    end if

    call Write_MLSAuxData(file_id, MLSData, error,write_attributes=.true., &
         fill_value_r=fill_value)

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Real_3d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Double_3d( file_id, dataset, double_data, &
      dims, fill_value)
    type( DataProducts_T ), intent(in) :: dataset
    real(r8), dimension(:,:,:), intent(in) :: double_data
    integer, dimension(:), intent(in), optional :: dims
    integer(hid_t), intent(in) :: file_id
    real(r8), intent(in), optional :: fill_value

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error

    dim_array = 1
    if (present(dims) ) then
       dim_array(:size(dims)) = dims
    else
       do i=1,size(shape(double_data))
          dim_array(i) = size(double_data,i)
       end do
    end if

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    MLSData%rank = 3

    MLSData%DpField(1:dim_array(1),1:dim_array(2),1:dim_array(3)) = &
      & double_data(1:dim_array(1),1:dim_array(2),1:dim_array(3))

    call CopyFromDataProducts (dataset, MLSData)
    if ( MLSData%name == "R1A:118.B22D:PT.S0.DACS-4 precision"  .and. DEBUG ) then
      print *, '3d double dim_array: ', dim_array
    end if

    call Write_MLSAuxData(file_id, MLSData, error,write_attributes=.true., &
         fill_value_d=fill_value)
    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData.' ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Double_3d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Integer_3d (file_id,dataset, integer_data, dims, &
      fill_value)
    type( DataProducts_T ), intent(in) :: dataset
    integer, dimension(:,:,:), intent(in) :: integer_data
    integer, dimension(:), intent(in), optional :: dims
    integer(hid_t), intent(in) :: file_id
    integer, intent(in), optional :: fill_value

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error

    dim_array = 1
    if (present(dims) ) then
       dim_array(:size(dims)) = dims
    else
       do i=1,size(shape(integer_data))
          dim_array(i) = size(integer_data,i)
       end do
    end if

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    MLSData%rank = 3

    MLSData%IntField(1:dim_array(1),1:dim_array(2),1:dim_array(3)) = &
      & integer_data(1:dim_array(1),1:dim_array(2),1:dim_array(3))

    call CopyFromDataProducts (dataset, MLSData)

    call Write_MLSAuxData(file_id, MLSData, error, write_attributes=.true., &
         fill_value_i=fill_value)
    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         'Error Writing MLSAuxData.' ) 
    call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Integer_3d
  !------------------------------------------------------- Recall_MLSAuxData
 subroutine Recall_MLSAuxData_Character( file_id, dataset, char_data, &
      char_length)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    character (len=*), intent(inout) :: char_data
    integer, intent(in) :: char_length
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dims
    integer :: error, status

    dims = 1

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),&
         trim(dataset%data_type),dims,MLSData)

    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type), MLSData, error, read_attributes=.true.)

    if (error == 0) then 
       char_data = MLSData%CharField(1,1,1) 

       allocate(dataset%Dimensions(1), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = 1, elementSize = storage_size(dataset%Dimensions) / 8, &
         & address=addr )
       call CopyToDataProducts (MLSData, dataset)

    else 
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) )
    end if

    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Character
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Integer( file_id, dataset, int_data )
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    integer, intent(inout) :: int_data
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dims
    integer :: error, status

    dims = 1

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dims,MLSData)

    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type),MLSData, error, read_attributes=.true.)

    if (error == 0) then 
       int_data = MLSData%IntField(1,1,1)

       allocate(dataset%Dimensions(1), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = 1, elementSize = storage_size(dataset%Dimensions) / 8, &
         & address=addr )
       call CopyToDataProducts(MLSData, dataset)

    else 
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    end if

    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Integer
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Real(file_id,dataset,real_data)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    real, intent(inout) :: real_data
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dims
    integer :: error, status

    dims = 1

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dims,MLSData)

    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type),MLSData, error, read_attributes=.true.)
    if (error == 0) then 
       real_data = MLSData%RealField(1,1,1)

       allocate(dataset%Dimensions(1), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = 1, elementSize = storage_size(dataset%Dimensions) / 8, &
         & address=addr )
       call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) )
    end if
    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Real
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Double( file_id, dataset, double_data)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    real(r8), intent(inout) :: double_data
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dims
    integer :: error, status

    dims = 1

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dims,MLSData)
    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type),MLSData, error, read_attributes=.true.)

    if (error == 0) then 
       double_data = MLSData%DpField(1,1,1)
       allocate(dataset%Dimensions(1), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = 1, elementSize = storage_size(dataset%Dimensions) / 8, &
         & address=addr )
       call CopyToDataProducts(MLSData, dataset)
    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) )
    end if

    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Double
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Real_1d( file_id, dataset, real_data, & 
      firstIndex, lastIndex)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    real, dimension(:), intent(inout) :: real_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dim_array
    integer :: i,error, status, i_first, i_last

    dim_array = 1

    do i=1,size(shape(real_data))
      dim_array(i) = size(real_data,i)
    end do

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)
    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type),MLSData, error, read_attributes=.true.)

    if (error == 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          i_first = firstIndex
          i_last  = lastIndex
       else
          i_first = 1
          i_last  = dim_array(1)
       end if

       allocate(dataset%Dimensions(MLSData%rank), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = MLSData%rank, &
         & elementSize = storage_size(dataset%Dimensions) / 8, address=addr )
       call CopyToDataProducts(MLSData, dataset)

       real_data(i_first:i_last) = MLSData%RealField(i_first:i_last,1,1)  

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) )
    end if

    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Real_1d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Double_1d(file_id, dataset, double_data, & 
      firstIndex, lastIndex)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    real(r8), dimension(:), intent(inout) :: double_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dim_array
    integer :: i,error, status, i_first, i_last

    dim_array = 1

    do i=1,size(shape(double_data))
       dim_array(i) = size(double_data,i)
    end do

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type),MLSData, error, read_attributes=.true.)

    if (error == 0) then 

       if (present (firstIndex) .and. present (lastIndex) ) then
          i_first = firstIndex
          i_last  = lastIndex
       else
          i_first = 1
          i_last  = dim_array(1)
       end if

       allocate(dataset%Dimensions(MLSData%rank), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = MLSData%rank, &
         & elementSize = storage_size(dataset%Dimensions) / 8, address=addr )
       call CopyToDataProducts(MLSData, dataset)

       do i = i_first, i_last
          double_data(i) = MLSData%DpField(i,1,1) 
       end do


    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) )
    end if

    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Double_1d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Integer_1d(file_id,dataset,integer_data, &
      firstIndex, lastIndex)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    integer, dimension(:), intent(inout) :: integer_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dim_array
    integer :: i,error,status,i_first,i_last

    dim_array = 1

    do i=1,size(shape(integer_data))
       dim_array(i) = size(integer_data,i)
    end do

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type),MLSData, error, read_attributes=.true.)

    if (error == 0) then
       if (present (firstIndex) .and. present (lastIndex) ) then
          i_first = firstIndex
          i_last  = lastIndex
       else
          i_first = 1
          i_last  = dim_array(1)
       end if
       do i = i_first, i_last
          integer_data(i) = MLSData%IntField(i,1,1)  
       end do

       allocate(dataset%Dimensions(MLSData%rank), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = MLSData%rank, &
         & elementSize = storage_size(dataset%Dimensions) / 8, address=addr )
       call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) )
    end if

    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Integer_1d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Real_2d( file_id, dataset, real_data, & 
      firstIndex, lastIndex)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    real, dimension(:,:), intent(inout) :: real_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dim_array
    integer :: i,error,status,j_first,j_last

    dim_array = 1

    do i=1,size(shape(real_data))
       dim_array(i) = size(real_data,i)
    end do

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type),MLSData, error, read_attributes=.true.)
    if (error == 0) then 

       if (present (firstIndex) .and. present (lastIndex) ) then
          j_first = firstIndex
          j_last  = lastIndex
       else
          j_first = 1
          j_last  = dim_array(2)
       end if

       real_data(1:dim_array(1),j_first:j_last) = &
         & MLSData%RealField(1:dim_array(1),j_first:j_last,1) 

       allocate(dataset%Dimensions(MLSData%rank), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = MLSData%rank, &
         & elementSize = storage_size(dataset%Dimensions) / 8, address=addr )
       call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) )
    end if

    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Real_2d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Double_2d( file_id, dataset, double_data,& 
      firstIndex, lastIndex)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    real(r8), dimension(:,:), intent(inout) :: double_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dim_array
    integer :: i,error,status,j_first,j_last

    dim_array = 1

    do i=1,size(shape(double_data))
       dim_array(i) = size(double_data,i)
    end do

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type), MLSData, error, read_attributes=.true.)
    if (error == 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          j_first = firstIndex
          j_last  = lastIndex
       else
          j_first = 1
          j_last  = dim_array(2)
       end if

       double_data(1:dim_array(1),j_first:j_last) = &
         & MLSData%DpField(1:dim_array(1),j_first:j_last,1) 

       allocate(dataset%Dimensions(MLSData%rank), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = MLSData%rank, &
         & elementSize = storage_size(dataset%Dimensions) / 8, address=addr )
       call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) )
    end if
    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Double_2d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Integer_2d( file_id, dataset, integer_data, &
      firstIndex, lastIndex)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    integer, dimension(:,:), intent(inout) :: integer_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dim_array
    integer :: i,error, status,j_first, j_last

    dim_array = 1

    do i=1,size(shape(integer_data))
       dim_array(i) = size(integer_data,i)
    end do

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type), MLSData, error, read_attributes=.true.)
    if (error == 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          j_first = firstIndex
          j_last  = lastIndex
       else
          j_first = 1
          j_last  = dim_array(2)
       end if

       integer_data(1:dim_array(1),j_first:j_last) = &
         & MLSData%IntField(1:dim_array(1),j_first:j_last,1) 

       allocate(dataset%Dimensions(MLSData%rank), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = MLSData%rank, &
         & elementSize = storage_size(dataset%Dimensions) / 8, address=addr )
       call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) )
    end if
    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Integer_2d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Real_3d( file_id, dataset, real_data, & 
      firstIndex, lastIndex)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    real, dimension(:,:,:), intent(inout) :: real_data
    integer(hid_t), intent(in) :: file_id
    integer, intent(in), optional :: firstIndex, lastIndex

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: i,error, status,k_first, k_last

    dim_array = 1

    do i=1,size(shape(real_data))
       dim_array(i) = size(real_data,i)
    end do

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type), MLSData, error, read_attributes=.true.)

    if (error == 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          k_first = firstIndex
          k_last  = lastIndex
       else
          k_first = 1
          k_last  = dim_array(3)
       end if

       real_data(1:dim_array(1),1:dim_array(2),1:dim_array(3)) = &
         & MLSData%RealField(1:dim_array(1),1:dim_array(2),1:dim_array(3))

       allocate(dataset%Dimensions(3), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = 3, elementSize = storage_size(dataset%Dimensions) / 8, &
         & address=addr )
       call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
            'Error Writing MLSAuxData for '// trim(dataset%name) )
    end if
    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Real_3d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Double_3d( file_id, dataset, double_data, &
      firstIndex, lastIndex)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    real(r8), dimension(:,:,:), intent(inout) :: double_data
    integer(hid_t), intent(in) :: file_id
    integer, intent(in), optional :: firstIndex, lastIndex

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dim_array
    integer :: i,error, status,k_first, k_last

    dim_array = 1

    do i=1,size(shape(double_data))
       dim_array(i) = size(double_data,i)
    end do

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)
    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type), MLSData, error, read_attributes=.true.)
    if (error == 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          k_first = firstIndex
          k_last  = lastIndex
       else
          k_first = 1
          k_last  = dim_array(3)
       end if

       double_data(1:dim_array(1),1:dim_array(2),1:dim_array(3)) = &
         & MLSData%DpField(1:dim_array(1),1:dim_array(2),1:dim_array(3))

       allocate(dataset%Dimensions(3), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = 3, elementSize = storage_size(dataset%Dimensions) / 8, &
         & address=addr )
       call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error,ModuleName,'Error Writing MLSAuxData.')
    end if
    call deallocate_mlsauxdata(MLSData)
  end subroutine Recall_MLSAuxData_Double_3d
  !------------------------------------------------------------------------------
  subroutine Recall_MLSAuxData_Integer_3d(file_id,dataset,integer_data, &
       firstIndex, lastIndex)
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( DataProducts_T ), intent(inout) :: dataset
    integer, dimension(:,:,:), intent(inout) :: integer_data
    integer(hid_t), intent(in) :: file_id
    integer, intent(in), optional :: firstIndex, lastIndex

    type( MLSAuxData_T ) :: MLSData
    integer(c_intptr_t) :: Addr         ! For tracing
    integer, dimension(3) :: dim_array
    integer :: i,error,status,k_first, k_last

    dim_array = 1

    do i=1,size(shape(integer_data))
       dim_array(i) = size(integer_data,i)
    end do

    call deallocate_mlsauxdata(MLSData)
    call Allocate_MLSAuxData(trim(dataset%name),& 
         trim(dataset%data_type),dim_array,MLSData)

    call Read_MLSAuxData(file_id, trim(dataset%name), &
         trim(dataset%data_type), MLSData, error, read_attributes=.true.)
    if (error == 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          k_first = firstIndex
          k_last  = lastIndex
       else
          k_first = 1
          k_last  = dim_array(3)
       end if

       ! dump the data into the array.

       integer_data(1:dim_array(1),1:dim_array(2),1:dim_array(3)) = &
         & MLSData%IntField(1:dim_array(1),1:dim_array(2),1:dim_array(3))

       allocate(dataset%Dimensions(3), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(dataset%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "dataset%Dimensions", &
         & ubounds = 3, elementSize = storage_size(dataset%Dimensions) / 8, &
         & address=addr )
       call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error,ModuleName,'Error Writing MLSAuxData.')
    end if
    call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Integer_3d
!------------------------------------------------------------------------------
  subroutine CreateGroup_MLSAuxData(loc_id, group_name, subgroup_name)
!
! External variables
!
    character(len=*), intent(in) :: group_name
    character(len=*), intent(in), optional :: subgroup_name
    integer(hid_t), intent(in)   :: loc_id ! From HDF
!
! Internal variables
!
    integer(hid_t) :: group_id, subgroup_id
    integer :: h5error

    if ( .not. present(subgroup_name) ) then 

       call h5gcreate_f(loc_id, trim(group_name), group_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_GROUP_CREATE // trim(group_name))

       call h5gclose_f(group_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_GROUP_CLOSE // trim(group_name)) 

    else

       call h5gopen_f(loc_id, trim(group_name), group_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_GROUP_OPEN // trim(group_name))

       call h5gcreate_f(group_id, trim(subgroup_name), subgroup_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_GROUP_CREATE // trim(subgroup_name))

       call h5gclose_f(subgroup_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_GROUP_CLOSE // trim(subgroup_name)) 

       call h5gclose_f(group_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_GROUP_CLOSE // trim(group_name)) 

    end if

  end subroutine CreateGroup_MLSAuxData
! -------------------------------------------------  Create_MLSAuxData ----
  subroutine Create_MLSAuxData(file_id, MLSAuxData, string_length, & 
       write_attributes)

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
! This subroutine creates an entry in the HDF5 file.
!----------------------------------------------------------------------
! External variables

    type( MLSAuxData_T ), intent(in) :: MLSAuxData
    integer(hid_t), intent(in)       :: file_id ! From HDF
    integer, intent(in), optional    :: string_length
    logical, intent(in), optional    :: write_attributes
!-----------------------------------------------------------------------
! Internal variables
!
    character(len=namelen) :: aname
    real, dimension(:), allocatable :: attr_data
    integer(hsize_t), dimension(7) :: adims
    integer(hsize_t), dimension(3) :: chunk_dims, dims, maxdims
    integer(hsize_t), dimension(1) :: adims_create
    integer(hid_t) :: cparms,dspace_id,dset_id,type_id, &
         attr_id, atype_id, aspace_id, s_type_id
    integer :: i, rank, arank, h5error, status
!-----------------------------------------------------------------------

    test_type: select case (trim(MLSAuxData%type_name))
    case ('real')
       if (associated(MLSAuxData%RealField)) then 
          rank = size(shape(MLSAuxData%RealField))
          type_id = H5T_IEEE_F32LE
             do i=1,rank-1
                dims(i) = size(MLSAuxData%RealField, i)
                chunk_dims(i) = dims(i)
                maxdims(i) = chunk_dims(i)
             end do
             dims(rank) = 1
             chunk_dims(rank) = 1
             maxdims(rank) = H5S_UNLIMITED_F
          end if
       case ('double')
          if (associated(MLSAuxData%DpField)) then 
             rank = size(shape(MLSAuxData%DpField))
             type_id = H5T_NATIVE_DOUBLE
                do i=1,rank-1
                   dims(i) = size(MLSAuxData%DpField, i)
                   chunk_dims(i) = dims(i)
                   maxdims(i) = chunk_dims(i)
                end do
                dims(rank) = 1
                chunk_dims(rank) = 1
                maxdims(rank) = H5S_UNLIMITED_F
             end if
          case ('integer')
             if (associated(MLSAuxData%IntField)) then 
                rank = size(shape(MLSAuxData%IntField))
                type_id = H5T_NATIVE_INTEGER
                   do i=1,rank-1
                      dims(i) = size(MLSAuxData%IntField, i)
                      chunk_dims(i) = dims(i)
                      maxdims(i) = chunk_dims(i)
                   end do
                   dims(rank) = 1
                   chunk_dims(rank) = 1
                   maxdims(rank) = H5S_UNLIMITED_F
                end if
             case ('character')
                if (associated(MLSAuxData%CharField)) then 
                   rank = size(shape(MLSAuxData%CharField))
                   type_id = H5T_NATIVE_CHARACTER
                      do i=1,rank-1
                         dims(i) = size(MLSAuxData%CharField, i)
                         chunk_dims(i) = dims(i)
                         maxdims(i) = chunk_dims(i)
                      end do
                      dims(rank) = 1
                      chunk_dims(rank) = 1
                      maxdims(rank) = H5S_UNLIMITED_F
                   end if
                end select test_type
!--------------------------------------------------------------------
                call h5screate_simple_f(rank, dims(1:rank),dspace_id,h5error,& 
                     maxdims(1:rank))
                if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
                     H5_ERROR_DSPACE_CREATE // trim(MLSAuxData%name) )

                call h5pcreate_f(H5P_DATASET_CREATE_F,cparms,h5error)
                if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
                     H5_ERROR_PROPERTY_CREATE // trim(MLSAuxData%name) )

                call h5pset_chunk_f(cparms,rank,chunk_dims(1:rank),h5error)
                if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
                     H5_ERROR_PROPERTY_CHUNK_SET // trim(MLSAuxData%name) )

                if (trim(MLSAuxData%type_name)/='character') then

                   call h5dcreate_f(file_id,trim(MLSAuxData%name),type_id,&
                        dspace_id, dset_id,h5error,cparms)
                   if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                        ModuleName, H5_ERROR_DSET_CREATE//trim(MLSAuxData%name))

                else

                   call h5tcopy_f(type_id, s_type_id, h5error)
                   if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                        ModuleName, H5_ERROR_TYPE_COPY // trim(MLSAuxData%name))

                   call h5tset_size_f(s_type_id, int(string_length, size_t), h5error)
                   if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                        & ModuleName, H5_ERROR_TYPE_SET // trim(MLSAuxData%name))

                   call h5dcreate_f(file_id, trim(MLSAuxData%name), s_type_id, &
                        dspace_id, dset_id, h5error, cparms)
                   if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                        ModuleName, H5_ERROR_DSET_CREATE // &
                        trim(MLSAuxData%name))

                end if

!---------------- write the attributes -----------------------------

                if (present(write_attributes)) then 

                   if (write_attributes) then 

                      do i = 2, 7
                         adims(i) = 0
                      end do


                      if (associated(MLSAuxData%HorizontalCoordinates)) then 

                         arank = size(shape(MLSAuxData%HorizontalCoordinates))
                         atype_id = H5T_IEEE_F32LE
                         adims(1) = size(MLSAuxData%HorizontalCoordinates)
                         adims_create(1) = adims(1)

                         allocate(attr_data(adims(1)),STAT=status)
                         call test_allocate ( status, ModuleName, 'attr_data' )

                         do i=1, adims(1)
                            attr_data(i) = MLSAuxData%HorizontalCoordinates(i)
                         end do

                         aname = "HorizontalCoordinates"  

                         call h5screate_simple_f (arank, adims_create, &
                              aspace_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_DSPACE_CREATE//trim(aname)//&
                              ' in '//trim(MLSAuxData%name) )

                         call h5acreate_f (dset_id, trim(aname), atype_id, &
                              aspace_id, attr_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_ATT_CREATE //trim(aname)// &
                              ' in ' // trim(MLSAuxData%name) )

                         call h5awrite_f (attr_id, atype_id, attr_data, adims, &
                              h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_ATT_WRITE // &
                              trim(MLSAuxData%name))

                         call h5aclose_f (attr_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_ATT_CLOSE // &
                              trim(MLSAuxData%name))

                         call h5sclose_f (aspace_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_DSPACE_CLOSE // &
                              trim(MLSAuxData%name))

                         deallocate(attr_data,STAT=status)
                         call test_deallocate ( status, ModuleName,'attr_data' )

                      end if

                      if (associated(MLSAuxData%VerticalCoordinates) ) then 
                         arank = size(shape(MLSAuxData%VerticalCoordinates))
                         atype_id = H5T_IEEE_F32LE
                         adims(1) = size(MLSAuxData%VerticalCoordinates)
                         adims_create(1) = adims(1)

                         allocate(attr_data(adims(1)),STAT=status)
                         call test_allocate ( status, ModuleName, 'attr_data' )

                         do i=1,adims(1)
                            attr_data(i) = MLSAuxData%VerticalCoordinates(i)
                         end do

                         call h5screate_simple_f (arank, adims_create, &
                              aspace_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_DSPACE_CREATE // &
                              trim(MLSAuxData%name))

                         aname = "VerticalCoordinates"   
                         call h5acreate_f (dset_id, trim(aname), atype_id, &
                              aspace_id, attr_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_ATT_CREATE // &
                              trim(MLSAuxData%name))

                         call h5awrite_f (attr_id, atype_id, attr_data, adims, &
                              h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_ATT_WRITE // &
                              trim(MLSAuxData%name))

                         call h5aclose_f (attr_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_ATT_CLOSE // &
                              trim(MLSAuxData%name))

                         call h5sclose_f (aspace_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_DSPACE_CLOSE // &
                              trim(MLSAuxData%name))

                         deallocate(attr_data,STAT=status)
                         call test_deallocate ( status, ModuleName,'attr_data' )

                      end if

                      if (associated(MLSAuxData%FrequencyCoordinates)) then 
                         arank = size(shape(MLSAuxData%FrequencyCoordinates))
                         atype_id = H5T_IEEE_F32LE
                         adims(1) = size(MLSAuxData%FrequencyCoordinates)
                         adims_create(1) = adims(1)

                         allocate(attr_data(adims(1)),STAT=status)
                         call test_allocate ( status, ModuleName, 'attr_data' )

                         do i=1,adims(1)
                            attr_data(i) = MLSAuxData%FrequencyCoordinates(i)
                         end do

                         call h5screate_simple_f (arank, adims_create, &
                              aspace_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_DSPACE_CREATE // &
                              trim(MLSAuxData%name) )

                         aname = "FrequencyCoordinates"
                         call h5acreate_f (dset_id, trim(aname), atype_id, &
                              aspace_id, attr_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_ATT_CREATE // &
                              trim(MLSAuxData%name))

                         call h5awrite_f (attr_id, atype_id, attr_data, adims, &
                              h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_ATT_WRITE // &
                              trim(MLSAuxData%name) )

                         call h5aclose_f (attr_id, h5error)
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_ATT_CLOSE // &
                              trim(MLSAuxData%name) )

                         call h5sclose_f (aspace_id, h5error) 
                         if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
                              ModuleName, H5_ERROR_DSPACE_CLOSE // &
                              trim(MLSAuxData%name) )

                         deallocate(attr_data,STAT=status)
                         call test_deallocate ( status, ModuleName,'attr_data' )
                      end if
                   end if
                end if

!---------------- close all structures----------------------------------
!
                call h5sclose_f( dspace_id, h5error )
                if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
                     H5_ERROR_DSPACE_CLOSE // trim(MLSAuxData%name) )

                call h5dclose_f( dset_id,   h5error )
                if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
                     H5_ERROR_DSET_CLOSE // trim(MLSAuxData%name) )

                call h5pclose_f( cparms,    h5error )
                if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
                     H5_ERROR_PROPERTY_CLOSE // trim(MLSAuxData%name) )

  end subroutine Create_MLSAuxData
!------------Read_MLSAuxAttributes ----
  subroutine Read_MLSAuxAttributes(dset_id, QuantityName, AttributeName, & 
       error, MLSAuxData, AttributeData)
  ! This subroutine reads attributes from the HDF5 file;
  ! e.g., FrequencyCoordinates
  ! They are returned either as:
  ! components of MLSAuxData (if supplied); or
  ! in the array AttributeData
  ! It is an error unless exactly one of the two is supplied
    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    character(len=*), intent(in)                  :: QuantityName
    character(len=*), intent(in)                  :: AttributeName
    integer(hid_t), intent(in)                    :: dset_id ! From h5dopen_f
    integer, intent(out)                          :: error   ! 0 unless trouble
    type( MLSAuxData_T ), intent(inout), optional :: MLSAuxData
    real, dimension(:), intent(out), optional     :: AttributeData
  ! Private
    character, dimension(:), allocatable :: char_data
    real, dimension(:), allocatable :: attr_data
    integer(hsize_t), dimension(7) :: adims    
    integer(hid_t) :: attr_id, atype_id, aspace_id
    integer, parameter :: OPTARGISMAD = 1
    integer, parameter :: OPTARGISARR = 2
    integer :: which_opt_arg, i, arank, h5error, status
    logical :: is_char
!------------------------------------------------------------------------
  ! Executable
    error = 0
    is_char = .false.
    atype_id = H5T_IEEE_F32LE
    do i = 2, 7
       adims(i) = 0
    end do
    ! Check that exactly 1 of the optional args is supplied
    which_opt_arg = 0
    if ( present(MLSAuxData) ) which_opt_arg = which_opt_arg + OPTARGISMAD
    if ( present(AttributeData) ) which_opt_arg = which_opt_arg + OPTARGISARR
    select case (which_opt_arg)
    case (OPTARGISMAD)
       select case (trim(AttributeName))
       case ('FrequencyCoordinates')
          if ( associated ( MLSAuxData%FrequencyCoordinates ) ) then
             arank = 1
             adims(1) = size(MLSAuxData%FrequencyCoordinates)
          else
             error = 1
             return
          end if
       case ('HorizontalCoordinates')
          if ( associated ( MLSAuxData%HorizontalCoordinates ) ) then
             arank = 1
             adims(1) = size(MLSAuxData%HorizontalCoordinates)
          else
             error = 1
             return
          end if
       case ('VerticalCoordinates')
          if ( associated ( MLSAuxData%VerticalCoordinates ) ) then
             arank = 1
             adims(1) = size(MLSAuxData%VerticalCoordinates)
          else
             error = 1
             return
          end if
       case ('Dimensions')
          if ( associated (MLSAuxData%Dimensions) ) then
             arank = 1
             is_char  = .true.
             atype_id = H5T_NATIVE_CHARACTER
             do i = 1, arank
                adims(i) = size(MLSAuxData%Dimensions)
             end do
          else
             error = 1
             return
          end if
       case default
          error = 1
          return
       end select
    case (OPTARGISARR)
       adims(1) = size(AttributeData)
    case default
       error = 1
       return
    end select

    call h5aopen_name_f(dset_id, trim(AttributeName), attr_id, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_ATT_OPEN // trim(AttributeName) // ' in ' // &
         QuantityName )

    call h5aget_space_f(attr_id, aspace_id, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_ATT_GETSPACE // QuantityName )

    if (is_char) then 

       allocate(char_data(adims(1)),stat=status)
       call test_allocate ( status, ModuleName, 'char_data' )

       call h5aread_f(attr_id, atype_id, char_data, adims, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_ATT_READ // trim(AttributeName) // ' in ' // &
            QuantityName )

    else

       allocate(attr_data(adims(1)),stat=status)
       call test_allocate ( status, ModuleName, 'attr_data' )

       call h5aread_f(attr_id, atype_id, attr_data, adims, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_ATT_READ // trim(AttributeName) // ' in ' // &
            QuantityName )

    end if

    select case (which_opt_arg)
    case (OPTARGISMAD)
       select case (trim(AttributeName))
       case ('FrequencyCoordinates')
          MLSAuxData%FrequencyCoordinates = attr_data(:adims(1))
       case ('HorizontalCoordinates')
          MLSAuxData%HorizontalCoordinates = attr_data(:adims(1))
       case ('VerticalCoordinates')
          MLSAuxData%VerticalCoordinates = attr_data(:adims(1))
       case ('Dimensions')
          MLSAuxData%Dimensions = char_data(:adims(1))
       end select

    case (OPTARGISARR)
       AttributeData = attr_data(:adims(1))
    end select

    call h5aclose_f(attr_id, h5error)                                  
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_ATT_CLOSE // trim(AttributeName) // ' in ' // &
         QuantityName )

    call h5sclose_f(aspace_id, h5error)                                
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_CLOSE // trim(AttributeName) // ' in ' // &
         QuantityName )

    if (allocated(attr_data)) then 
       deallocate(attr_data,STAT=status)
       call test_deallocate ( status, ModuleName, 'attr_data' )
    end if

    if (allocated(char_data)) then      
       deallocate(char_data,STAT=status)
       call test_deallocate ( status, ModuleName, 'char_data' )
    end if

  end subroutine Read_MLSAuxAttributes
  !-------------------------------------------------Write_MLSAuxAttributes ----
  subroutine Write_MLSAuxAttributes(dset_id, QuantityName, AttributeName, & 
       error, MLSAuxData, AttributeData)

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    ! This subroutine writes attributes to the HDF5 file;
    ! e.g., FrequencyCoordinates
    ! They are supplied either as:
    ! components of MLSAuxData (if supplied); or
    ! in the array AttributeData
    ! It is an error unless exactly one of the two is supplied
    character(len=*), intent(in)                  :: QuantityName
    character(len=*), intent(in)                  :: AttributeName
    integer(hid_t), intent(in)                    :: dset_id ! From h5dcreate_f
    integer, intent(out)                          :: error
    type( MLSAuxData_T ), intent(in), optional :: MLSAuxData
    real, dimension(:), intent(in), optional      :: AttributeData
    ! Private
    character(len=20), dimension(:), allocatable :: char_data
    real, dimension(:), allocatable :: attr_data
    integer(hsize_t), dimension(7) :: adims
    integer(hsize_t), dimension(1) :: adims_create
    integer(hid_t) :: atype_id, aspace_id, attr_id
    integer, parameter :: OPTARGISMAD = 1
    integer, parameter :: OPTARGISARR = 2
    integer :: which_opt_arg, i, arank, attrlen, h5error, status  
    logical :: is_char
    !-----------------------------------------------------------------------
    ! Executable
    error = 0
    atype_id = H5T_IEEE_F32LE
    is_char = .false.
    adims = 1
    ! Check that exactly 1 of the optional args is supplied
    ! The following will produce 0 if neither, 3 if both; else 1 or 2
    which_opt_arg = 0
    if ( present(MLSAuxData) ) which_opt_arg = which_opt_arg + OPTARGISMAD
    if ( present(AttributeData) ) which_opt_arg = which_opt_arg + OPTARGISARR
    select case (which_opt_arg)
    case (OPTARGISMAD)
       select case (trim(AttributeName))
       case ('FrequencyCoordinates')
          arank = 1
          adims(1) = size(MLSAuxData%FrequencyCoordinates)
          allocate(attr_data(adims(1)),STAT=status)
          call test_allocate ( status, ModuleName, 'attr_data' )

          do i=1,adims(1)
             attr_data(i) = MLSAuxData%FrequencyCoordinates(i)
          end do
       case ('HorizontalCoordinates')
          arank = 1
          adims(1) = size(MLSAuxData%HorizontalCoordinates)
          allocate(attr_data(adims(1)),STAT=status)
          call test_allocate ( status, ModuleName, 'attr_data' )

          do i=1,adims(1)
             attr_data(i) = MLSAuxData%HorizontalCoordinates(i)
          end do
       case ('VerticalCoordinates')
          arank = 1
          adims(1) = size(MLSAuxData%VerticalCoordinates)
          allocate(attr_data(adims(1)),STAT=status)
          call test_allocate ( status, ModuleName, 'attr_data' )

          do i=1,adims(1)
             attr_data(i) = MLSAuxData%VerticalCoordinates(i)
          end do
       case ('Dimensions')
          arank = 1
          atype_id = H5T_NATIVE_CHARACTER
          is_char  = .true.
          adims(1) = size(MLSAuxData%Dimensions)
          allocate(char_data(adims(1)),STAT=status)
          call test_allocate ( status, ModuleName, 'char_data' )
          do i=1,adims(1)
             char_data(i) = trim(MLSAuxData%Dimensions(i))
          end do
       case default
          error = 1
          return
       end select
    case (OPTARGISARR)
       arank = 1
       adims(1) = size(AttributeData)
       allocate(attr_data(adims(1)),STAT=status)
       call test_allocate ( status, ModuleName, 'attr_data' )

       do i=1,adims(1)
          attr_data(i) = AttributeData(i)
       end do
    case default
       error = 1
       return
    end select

    adims_create(1) = adims(1)                                         
    call h5screate_simple_f(arank, adims_create(1:arank), aspace_id, h5error)  
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_CREATE // trim(AttributeName) // ' in ' // & 
         trim(MLSAuxData%name) )

    if (is_char) then 

       attrlen = 20 

       CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_TYPE_COPY // trim(AttributeName) // ' in ' // &
            trim(MLSAuxData%name))

       CALL h5tset_size_f(atype_id, int(attrlen, size_t), h5error)
       IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_TYPE_SET // trim(AttributeName) // ' in ' // &
            trim(MLSAuxData%name))

       call h5acreate_f(dset_id, trim(AttributeName), atype_id, &         
            & aspace_id, attr_id, h5error)                                   
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_ATT_CREATE // trim(AttributeName) // ' in ' // &
            trim(MLSAuxData%name) )

       call h5awrite_f(attr_id, atype_id, char_data, adims, h5error)      
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_ATT_WRITE // trim(AttributeName) // ' in ' // &
            trim(MLSAuxData%name) )

    else

       call h5acreate_f(dset_id, trim(AttributeName), atype_id, &         
            & aspace_id, attr_id, h5error)                                   
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_ATT_CREATE // trim(AttributeName) // ' in ' // &
            trim(MLSAuxData%name) )

       call h5awrite_f(attr_id, atype_id, attr_data, adims, h5error)      
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_ATT_WRITE // trim(AttributeName) // ' in ' // &
            trim(MLSAuxData%name) )

    end if

    call h5aclose_f(attr_id, h5error)                                  
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_ATT_CLOSE // trim(AttributeName) // ' in ' // &
         trim(MLSAuxData%name) )

    call h5sclose_f(aspace_id, h5error)                                
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_CLOSE // trim(AttributeName) // ' in ' // & 
         trim(MLSAuxData%name) )

    if (allocated(attr_data)) then 
       deallocate(attr_data,STAT=status)
       call test_deallocate ( status, moduleName, 'attr_data' )
    end if

    if (allocated(char_data)) then 
       deallocate(char_data,STAT=status)
       call test_deallocate ( status, moduleName, 'char_data' )
    end if

  end subroutine Write_MLSAuxAttributes
  !-------------------------------------------------Read_MLSAuxData ----
  subroutine Read_MLSAuxData(file_id, QuantityName, QuantityType, & 
       MLSAuxData, error, read_attributes, FirstIndex, LastIndex)

    ! This subroutine reads an entry from the HDF5 file.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    type( MLSAuxData_T ), intent(inout) :: MLSAuxData
    character(len=*), intent(in)   :: QuantityName
    character(len=*), intent(in)   :: QuantityType ! if 'unknown' will derive
    integer(hid_t), intent(in)     :: file_id ! From HDF
    integer, intent(out)           :: error   ! 0 unless trouble
    integer, intent(in), optional  :: FirstIndex, LastIndex
    logical, intent(in), optional  :: read_attributes  ! read freqCoords, etc 
!
    !-------------------------------------------------------------------------
!
    character, dimension(:,:,:), allocatable :: character_buffer
    character(len=16) :: myQuantityType
!
    real, dimension(:,:,:), allocatable :: real_buffer, double_buffer
!
    integer, dimension(:,:,:), allocatable :: integer_buffer
    integer(hsize_t), dimension(7) :: dims
    integer(hsize_t), dimension(3) :: chunk_dims, dims_create, maxdims, start
    integer(hid_t) :: cparms, dset_id, dspace_id, type_id, memspace
    integer        :: rank, h5error, i, status
!
    logical :: myRead_attributes, is_integer, is_real, is_double, & 
         is_character, is_int32, is_float32, is_float64
!
    error = 0
    myRead_attributes = .false.
    dims = 1
    if ( present(read_attributes) ) myRead_attributes=read_attributes
    !-- assign name and type for dataset

    MLSAuxData%name = trim(QuantityName)
    MLSAuxData%type_name = trim(QuantityType)

    call h5dopen_f(file_id,trim(QuantityName),dset_id,h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_OPEN // trim(MLSAuxData%name) )

    call h5dget_space_f(dset_id,dspace_id,h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_GETSPACE // trim(MLSAuxData%name) )

    call h5sget_simple_extent_ndims_f(dspace_id,rank,h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_GET_NDIMS // trim(MLSAuxData%name) )

    MLSAuxData%rank = rank

    call h5sget_simple_extent_dims_f(dspace_id,dims_create(1:rank),&
         maxdims(1:rank),h5error)
    if (h5error /= rank) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_GET_DIMS // trim(MLSAuxData%name) )

    call h5dget_create_plist_f(dset_id,cparms,h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_CREATE // trim(MLSAuxData%name) )

    call h5pget_chunk_f(cparms,rank,chunk_dims(1:rank),h5error)
    if (h5error /= rank) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_PROPERTY_CHUNK_GET // trim(MLSAuxData%name) )

    call h5screate_simple_f(rank, dims_create(1:rank), memspace, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_CREATE // trim(MLSAuxData%name) )

    call h5dget_type_f(dset_id, type_id, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_GET_TYPE // trim(MLSAuxData%name) )

    if ( lowercase(trim(QuantityType)) /= 'unknown' ) then
       myQuantityType = QuantityType
    else

       is_integer = .FALSE.
       is_real = .FALSE.
       is_double = .FALSE.
       is_int32 = .FALSE.
       is_float32 = .FALSE.
       is_float64 = .FALSE.
       is_character = .FALSE.

       call h5tequal_f(type_id, H5T_NATIVE_INTEGER, is_integer, h5error)
       call h5tequal_f(type_id, H5T_STD_I32LE, is_int32, h5error)
       call h5tequal_f(type_id, H5T_NATIVE_CHARACTER, is_character, h5error)
       call h5tequal_f(type_id, H5T_NATIVE_REAL, is_real, h5error)
       call h5tequal_f(type_id, H5T_IEEE_F32LE, is_float32, h5error)
       call h5tequal_f(type_id, H5T_IEEE_F64LE, is_float64, h5error)
       call h5tequal_f(type_id, H5T_NATIVE_DOUBLE, is_double, h5error)

       ! Please add any new or compound datatypes here.

       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
            H5_ERROR_TYPE_EQUAL // trim(MLSAuxData%name) )

       if ( is_integer .or. is_int32 ) then
          myQuantityType = 'integer'
       else if ( is_real .or. is_float32 ) then
          myQuantityType = 'real'
       else if ( is_double .or. is_float64 ) then
          myQuantityType = 'double'
       else if ( is_character ) then
          myQuantityType = 'character'
       else
          error = 1
          return
       end if

    end if

    do i = 1, 7
       if (i .le. rank) then 
          dims(i) = dims_create(i) 
       else
          dims(i) = 0
       end if
    end do

    do i = 1, rank
       start(i) = 0
    end do

    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
         start(1:rank), dims_create(1:rank), h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_HYPERSLAB // trim(MLSAuxData%name) )

    ! based on type chose which h5dread to call.

    if ( .NOT.(PRESENT(FirstIndex) ) ) then 

       test_type: select case (trim(myQuantityType))
       case ('real')
          call h5dread_f(dset_id,type_id,MLSAuxData%RealField,dims,& 
               h5error, memspace,dspace_id)
       case ('double')    
          call h5dread_f(dset_id,type_id,MLSAuxData%DpField,dims,&
               h5error,memspace,dspace_id)
       case ('integer')  
          call h5dread_f(dset_id,type_id,MLSAuxData%IntField,dims,& 
               h5error,memspace,dspace_id)
       case ('character') 
          call h5dread_f(dset_id,type_id,MLSAuxData%CharField,dims,&
               h5error,memspace,dspace_id)
       end select test_type

    else

       test_type_dims: select case (trim(myQuantityType))
       case ('real')
          allocate( real_buffer(dims(1),dims(2),dims(3)),stat=status)
          call test_allocate ( status, ModuleName, 'real_buffer' )

          call h5dread_f(dset_id,type_id,& 
               real_buffer,dims,h5error, memspace, dspace_id)

          MLSAuxData%RealField(1:dims(1),1:dims(2),firstIndex:lastIndex) = &
            & real_buffer(1:dims(1),1:dims(2),firstIndex:lastIndex)

          if (allocated(real_buffer)) then 
             deallocate(real_buffer, stat=status)
             call test_deallocate ( status, ModuleName, 'real_buffer' )
          end if

       case ('double')    

          allocate( double_buffer(dims(1),dims(2),dims(3)),stat=status)
          call test_allocate ( status, ModuleName, 'double_buffer' )

          call h5dread_f(dset_id,type_id,& 
               double_buffer,dims,h5error, memspace, dspace_id)

          MLSAuxData%DpField(1:dims(1),1:dims(2),firstIndex:lastIndex) = &
            & double_buffer(1:dims(1),1:dims(2),firstIndex:lastIndex)

          if (allocated(double_buffer)) then 
             deallocate(double_buffer, stat=status)
             call test_deallocate ( status, ModuleName, 'double_buffer' )
          end if

       case ('integer')  
          allocate( integer_buffer(dims(1),dims(2),dims(3)),stat=status)
          call test_allocate ( status, ModuleName, 'integer_buffer' )

          call h5dread_f(dset_id,type_id,& 
               integer_buffer,dims,h5error, memspace, dspace_id)

          MLSAuxData%IntField(1:dims(1),1:dims(2),firstIndex:lastIndex) = &
            & integer_buffer(1:dims(1),1:dims(2),firstIndex:lastIndex)

          if (allocated(integer_buffer)) then 
             deallocate(integer_buffer, stat=status)
             call test_deallocate ( status, ModuleName, 'integer_buffer' )
          end if

       case ('character') 
          allocate( character_buffer(dims(1),dims(2),dims(3)),stat=status)
          call test_allocate ( status, ModuleName, 'character_buffer' )

          call h5dread_f(dset_id,type_id,& 
               character_buffer,dims,h5error, memspace, dspace_id)

          MLSAuxData%CharField(1:dims(1),1:dims(2),firstIndex:lastIndex) = &
            & character_buffer(1:dims(1),1:dims(2),firstIndex:lastIndex)

          if (allocated(character_buffer)) then 
             deallocate(character_buffer, stat=status)
             call test_deallocate ( status, ModuleName, 'character_buffer' )
          end if

       end select test_type_dims

    end if

    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_READ // trim(MLSAuxData%name) )

    ! Do we also read the attributes
    if ( myRead_attributes ) then
       call Read_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
            & 'FrequencyCoordinates', h5error, MLSAuxData)
       call Read_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
            & 'HorizontalCoordinates', h5error, MLSAuxData)
       call Read_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
            & 'VerticalCoordinates', h5error, MLSAuxData)
       call Read_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
            & 'Dimensions', h5error, MLSAuxData)
    end if
!
    ! Close all identifiers.
!
    call h5sclose_f(memspace, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_CLOSE // trim(MLSAuxData%name) )

    call h5sclose_f(dspace_id, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_CLOSE // trim(MLSAuxData%name) )

    call h5dclose_f(dset_id, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_CLOSE // trim(MLSAuxData%name) )

    call h5pclose_f(cparms, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_PROPERTY_CLOSE // trim(MLSAuxData%name) )

  end subroutine Read_MLSAuxData
! -------------------------------------------------  Write_MLSAuxData ----
  subroutine Write_MLSAuxData (file_id, MLSAuxData, error, &
       & write_attributes, string_length, index, fill_value_r, fill_value_d, &
       fill_value_i)
!
! This routine assumes that the MLSAuxData dataset is created.
! This subroutine writes an entry to the HDF5 file.
!
    type( MLSAuxData_T ), intent(in) :: MLSAuxData
    integer(hid_t), intent(in)     :: file_id ! From HDF
    integer, intent(out)           :: error   ! 0 unless trouble
    integer, intent(in), optional  :: index ! number of major frames
    logical, intent(in), optional  :: write_attributes ! write freqCoords, etc 
    integer, intent(in), optional  :: string_length ! length of string  
    real, intent(in), optional     :: fill_value_r    ! value to use for filling
    real(r8), intent(in), optional :: fill_value_d    ! value to use for filling
    integer, intent(in), optional  :: fill_value_i    ! value to use for filling
!-------------------------------------------------------------------------
! Internal variables.
!
    character(len=1), dimension(:), pointer :: char_data => NULL()
    real, dimension(:), pointer :: attr_data => NULL()
    integer(hsize_t), dimension(7) :: dims
    integer(hsize_t), dimension(3) :: chunk_dims, dims_create, maxdims, start
    integer(hsize_t), dimension(3) :: dims_file
    integer(hid_t) :: cparms,dspace_id,dset_id,type_id, &
         filespace, memspace, s_type_id 
    integer :: i, rank, h5error
    logical :: myWrite_attributes
    !--------------------------------------------------------------------------
      call switchOutput( 'stdout' )
    error = 0
    myWrite_attributes = .false.
    if ( present(write_attributes) ) myWrite_attributes=write_attributes

    nullify (char_data, attr_data)
    dims = 1
    test_type: select case (trim(MLSAuxData%type_name))

    case ('real')
       if (associated(MLSAuxData%RealField)) then 
          rank = MLSAuxData%rank
          type_id = H5T_IEEE_F32LE
          if (present(index)) then 
             do i=1,rank-1
                dims(i) = size(MLSAuxData%RealField, i)
                chunk_dims(i) = dims(i)
                dims_create(i) = dims(i)
                maxdims(i) = dims(i)
             end do
             dims(rank) = 1
             chunk_dims(rank) = 1
             dims_create(rank) = 1
             maxdims(rank) = H5S_UNLIMITED_F
          else
             do i=1,rank
                dims(i) = size(MLSAuxData%RealField, i)
                chunk_dims(i) = dims(i)
                dims_create(i) = dims(i)
                maxdims(i) = dims(i)
             end do
          end if
       end if

    case ('double')
       if (associated(MLSAuxData%DpField)) then 
          rank = MLSAuxData%rank
          type_id = H5T_NATIVE_DOUBLE
          if (present(index)) then 
             do i=1,rank-1
                dims(i) = size(MLSAuxData%DpField, i)
                chunk_dims(i) = dims(i)
                dims_create(i) = dims(i)
                maxdims(i) = dims(i)
             end do
             dims(rank) = 1
             chunk_dims(rank) = 1
             dims_create(rank) = 1
             maxdims(rank) = H5S_UNLIMITED_F
          else
             do i=1,rank
                dims(i) = size(MLSAuxData%DpField, i)
                chunk_dims(i) = dims(i)
                dims_create(i) = dims(i)
                maxdims(i) = dims(i)
             end do
          end if
       end if
    case ('integer')

       if (associated(MLSAuxData%IntField)) then 
          rank = MLSAuxData%rank
          type_id = H5T_NATIVE_INTEGER
          if (present(index)) then 
             do i=1,rank-1
                dims(i) = size(MLSAuxData%IntField, i)
                chunk_dims(i) = dims(i)
                dims_create(i) = dims(i)
                maxdims(i) = dims(i)
             end do
             dims(rank) = 1
             chunk_dims(rank) = 1
             dims_create(rank) = 1
             maxdims(rank) = H5S_UNLIMITED_F
          else
             do i=1,rank
                dims(i) = size(MLSAuxData%IntField, i)
                chunk_dims(i) = dims(i)
                dims_create(i) = dims(i)
                maxdims(i) = dims(i)
             end do
          end if
       end if

    case ('character')
       if (associated(MLSAuxData%CharField)) then 
          rank = MLSAuxData%rank
          type_id = H5T_NATIVE_CHARACTER
          if (present (index) ) then 
             do i=1, rank-1
                dims(i) = 1
                chunk_dims(i) = dims(i)
                dims_create(i) = dims(i)
                maxdims(i) = dims(i)
             end do
             dims(rank) = 1
             chunk_dims(rank) = 1
             dims_create(rank) = 1
             maxdims(rank) = H5S_UNLIMITED_F
          else
             do i=1,rank
                dims(i) = size(MLSAuxData%CharField, i)
                chunk_dims(i) = dims(i)
                dims_create(i) = dims(i)
                maxdims(i) = dims(i)
             end do
          end if
       end if
    end select test_type

    call h5eset_auto_f(0, h5error)
    call h5dopen_f(file_id, trim(MLSAuxData%name), dset_id, h5error)

    if (h5error /= 0) then 

       if ( MLSAuxData%name == 'R1A:118.B22D:PT.S0.DACS-4 precision'  .and. DEBUG ) then
         print *, 'The dataset doesnt exist yet--must create it', h5error, dset_id
       end if
       call h5eset_auto_f(1, h5error)
!--------------------------------------------------------------------
       call h5screate_simple_f (rank, dims_create(1:rank), &
            dspace_id,h5error, maxdims(1:rank))
       if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
            ModuleName, H5_ERROR_DSPACE_CREATE // trim(MLSAuxData%name))

       call h5pcreate_f (H5P_DATASET_CREATE_F,cparms,h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, &
            H5_ERROR_PROPERTY_CREATE // trim(MLSAuxData%name))

       call h5pset_chunk_f (cparms,rank,chunk_dims(1:rank),h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, &
            H5_ERROR_PROPERTY_CHUNK_SET // trim(MLSAuxData%name) )

       IF (PRESENT (fill_value_r)) THEN
          CALL H5pSet_Fill_Value_f (cparms, type_id, fill_value_r, h5error)
          IF (h5error /= 0) CALL MLSMessage (MLSMSG_Error, &
               ModuleName, "Fill value " // TRIM (MLSAuxData%name))
       end if

       IF (PRESENT (fill_value_d)) THEN
          CALL H5pSet_Fill_Value_f (cparms, type_id, fill_value_d, h5error)
          IF (h5error /= 0) CALL MLSMessage (MLSMSG_Error, &
               ModuleName, "Fill value " // TRIM (MLSAuxData%name))
       end if

       IF (PRESENT (fill_value_i)) THEN
          CALL H5pSet_Fill_Value_f (cparms, type_id, fill_value_i, h5error)
          IF (h5error /= 0) CALL MLSMessage (MLSMSG_Error, &
               ModuleName, "Fill value " // TRIM (MLSAuxData%name))
       end if

       if (trim(MLSAuxData%type_name)/='character') then

          call h5dcreate_f (file_id, trim(MLSAuxData%name), &
               type_id, dspace_id, dset_id, h5error, cparms)
          if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
               ModuleName, H5_ERROR_DSET_CREATE // trim(MLSAuxData%name))

       else

          call h5tcopy_f (type_id, s_type_id, h5error)
          if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
               ModuleName, H5_ERROR_TYPE_COPY // trim(MLSAuxData%name))

          call h5tset_size_f (s_type_id, int(string_length, size_t), h5error)
          IF (h5error /= 0) CALL MLSMessage (MLSMSG_Error, &
               ModuleName, H5_ERROR_TYPE_SET//trim(MLSAuxData%name))

          CALL h5dcreate_f (file_id, trim(MLSAuxData%name), &
               s_type_id, dspace_id, dset_id, h5error, cparms)
          IF (h5error /= 0) CALL MLSMessage (MLSMSG_Error, &
               ModuleName, H5_ERROR_DSET_CREATE // trim(MLSAuxData%name))

       end if
!---------------- close all structures----------------------------------
       call h5sclose_f (dspace_id, h5error )
       if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
            ModuleName, H5_ERROR_DSPACE_CLOSE // trim(MLSAuxData%name))

       call h5pclose_f (cparms, h5error)
       if (h5error /= 0) call MLSMessage (MLSMSG_Error, &
            ModuleName, H5_ERROR_PROPERTY_CLOSE // trim(MLSAuxData%name))

    end if

    if (associated(MLSAuxData%RealField)) then

       rank = MLSAuxData%rank

       do i=1,rank
          dims_create(i) = size(MLSAuxData%RealField,i)
          dims(i) = dims_create(i)
          chunk_dims(i) = dims(i)
       end do
    end if

    if (associated(MLSAuxData%DpField)) then 

       rank = MLSAuxData%rank

       do i=1,rank
          dims_create(i) = size(MLSAuxData%DpField,i)
          dims(i) = dims_create(i)
          chunk_dims(i) = dims(i)
       end do
    end if

    if (associated(MLSAuxData%CharField)) then 

       rank = MLSAuxData%rank

       do i=1,rank
          dims_create(i) = 1
          dims(i) = 1
          chunk_dims(i) = 1
       end do

    end if

    if (associated(MLSAuxData%IntField)) then 

       rank = MLSAuxData%rank

       do i=1,rank
          dims_create(i) = size(MLSAuxData%IntField,i)
          dims(i) = dims_create(i)
          chunk_dims(i) = dims(i)
       end do
    end if

    do i = 1, rank
       start(i) = 0
    end do

    if (PRESENT(index)) then 
       start(rank) = index-1
       dims_create(rank) = index
       dims(rank) = index
       chunk_dims(rank) = index
    end if
    call h5dget_space_f(dset_id, filespace, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_GETSPACE // trim(MLSAuxData%name) )

  ! Wonder what the existing dims are?
    call h5sget_simple_extent_dims_f(filespace, dims_file(1:rank),&
         maxdims(1:rank),h5error)
    if ( MLSAuxData%name == 'R1A:118.B22D:PT.S0.DACS-4 precision' .and. DEBUG ) then
      print *, 'Before call to extend'
      write(*, '(a12, 8I8)') 'start      : ', start(1:rank)
      write(*, '(a12, 8I8)') 'dims       : ', dims(1:rank)
      write(*, '(a12, 8I8)') 'dims_create: ', dims_create(1:rank)
      write(*, '(a12, 8I8)') 'dims_file: ', dims_file(1:rank)
      write(*, '(a12, 8I8)') 'maxdims    : ', maxdims(1:rank)
    end if
    ! call h5dextend_f(dset_id, dims_create, h5error)
    do i=1, rank
      dims_file(i) = max( dims_file(i), dims_create(i) )
    end do
    ! hdf5 1.8 introduced h5dset_extent_f to rplace h5dextend_f
    ! which has been deprecated
    ! Therefore swap comments between the following two statements
    ! when necessary.
    ! call h5dset_extent_f(dset_id, dims_file, h5error)
    call h5dextend_f(dset_id, dims_file, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_EXTEND // trim(MLSAuxData%name) )

    call h5dget_space_f(dset_id, filespace, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_GETSPACE // trim(MLSAuxData%name) )

  ! Wonder what the existing dims are?
    call h5sget_simple_extent_dims_f(filespace, dims_file(1:rank),&
         maxdims(1:rank),h5error)
    if ( MLSAuxData%name == 'R1A:118.B22D:PT.S0.DACS-4 precision'  .and. DEBUG ) then
      print *, 'After dextend'
      write(*, '(a12, 8I8)') 'dims_file: ', dims_file(1:rank)
      write(*, '(a12, 8I8)') 'maxdims    : ', maxdims(1:rank)
    end if

    call h5dget_type_f(dset_id, type_id, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_GET_TYPE // trim(MLSAuxData%name) )

    dims_create(rank)=1

    call h5screate_simple_f (rank, dims_create(1:rank), memspace, h5error)
    if (h5error /= 0) call MLSMessage (MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_CREATE // trim(MLSAuxData%name))

  ! Wonder what the existing dims are?
    call h5sget_simple_extent_dims_f(filespace, dims_file(1:rank),&
         maxdims(1:rank),h5error)
    if ( MLSAuxData%name == 'R1A:118.B22D:PT.S0.DACS-4 precision' .and. DEBUG ) then
      print *, 'Before hyperslab'
      write(*, '(a12, 8I8)') 'dims_file: ', dims_file(1:rank)
      write(*, '(a12, 8I8)') 'maxdims    : ', maxdims(1:rank)
    end if
    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, start(1:rank), &
         dims_create(1:rank), h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_HYPERSLAB // trim(MLSAuxData%name) )

    test_type_name: select case (trim(MLSAuxData%type_name))

    case ('real')
       call h5dwrite_f(dset_id, type_id, MLSAuxData%RealField, &
            dims(1:rank), h5error, memspace, filespace)
    case ('double')
       call h5dwrite_f (dset_id, type_id, MLSAuxData%DpField, &
            dims(1:rank), h5error, memspace, filespace)
    case ('integer')
       call h5dwrite_f (dset_id, type_id, MLSAuxData%IntField, &
            dims(1:rank), h5error, memspace, filespace)
    case('character')
       call h5tcopy_f (type_id, s_type_id, h5error)
       IF (h5error /= 0) call MLSMessage (MLSMSG_Error, &
            ModuleName, H5_ERROR_TYPE_COPY // trim(MLSAuxData%name))
       
       call h5tset_size_f (s_type_id, int(string_length, size_t), h5error)
       IF (h5error /= 0) call MLSMessage (MLSMSG_Error, &
            ModuleName, H5_ERROR_TYPE_SET // trim(MLSAuxData%name))

       call h5dwrite_f (dset_id, s_type_id, MLSAuxData%CharField, & 
            dims(1:rank), h5error, memspace, filespace)
  
    end select test_type_name

    if (h5error /= 0) call MLSMessage (MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_WRITE // trim(MLSAuxData%name))

! Do we also write the attributes
    if ( myWrite_attributes ) then
       if (associated(MLSAuxData%FrequencyCoordinates)) & 
            call Write_MLSAuxAttributes (dset_id, &
            trim(MLSAuxData%name), 'FrequencyCoordinates', error, MLSAuxData)
       if (associated(MLSAuxData%HorizontalCoordinates)) & 
            call Write_MLSAuxAttributes (dset_id, &
            trim(MLSAuxData%name), 'HorizontalCoordinates', error, MLSAuxData)
       if (associated(MLSAuxData%VerticalCoordinates)) & 
            call Write_MLSAuxAttributes (dset_id, &
            trim(MLSAuxData%name), 'VerticalCoordinates', error, MLSAuxData)
       if (associated(MLSAuxData%Dimensions)) & 
            call Write_MLSAuxAttributes (dset_id, &
            trim(MLSAuxData%name), 'Dimensions', error, MLSAuxData)
    end if

!
! Close all identifiers.
!
    call h5sclose_f(filespace, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_CLOSE // trim(MLSAuxData%name) )

    call h5sclose_f(memspace, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSPACE_CLOSE // trim(MLSAuxData%name) )

    call h5dclose_f(dset_id, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
         H5_ERROR_DSET_CLOSE // trim(MLSAuxData%name) )

    if ( MLSAuxData%name == 'R1A:118.B22D:PT.S0.DACS-4 precision' .and. DEBUG ) then
      print *, 'After closing all the identifiers'
      write(*, '(a12, 8I8)') 'start      : ', start(1:rank)
      write(*, '(a12, 8I8)') 'dims       : ', dims(1:rank)
      write(*, '(a12, 8I8)') 'dims_create: ', dims_create(1:rank)
      write(*, '(a12, 8I8)') 'maxdims    : ', maxdims(1:rank)
      print *, '(reopening dspace)'
      call h5dopen_f(file_id, trim(MLSAuxData%name), dset_id, h5error)
      call h5dget_space_f(dset_id,dspace_id,h5error)
      call h5sget_simple_extent_ndims_f(dspace_id,rank,h5error)
      call h5sget_simple_extent_dims_f(dspace_id,dims_create(1:rank),&
           maxdims(1:rank),h5error)
      write(*, '(a12, 8I8)') 'dims_create: ', dims_create(1:rank)
      write(*, '(a12, 8I8)') 'maxdims    : ', maxdims(1:rank)
      call h5dclose_f(dset_id, h5error)
    end if
    call revertOutput
  end subroutine Write_MLSAuxData

  subroutine CopyFromDataProducts( dataset, MLSAuxData )
    ! It's unclear to me what the original developer's vision for this
    ! subroutine is. However, since the name of the subroutine is "Copy,"
    ! it is reasonable to assume that the fields Dimensions,
    ! FrequencyCoordinates, VerticalCoordinates, and HorizontalCoordinates
    ! of MLSAuxData will be identical to those of dataset once this
    ! subroutine finishes given that it does not violate the constraints
    ! of MLSAuxData's data type. In addition, since it does not conflict with
    ! the current use of this subroutine, I will mandate that this 
    ! subroutine will allocate those fields for MLSAuxData. This means
    ! that before this subroutine is called, MLSAuxData should not
    ! have those fields allocated, otherwise, the program will stop to
    ! avoid memory leak.
    ! Concerning this subroutine, MLSAuxData_T type has one constraint
    ! that dictates the size of MLSAuxData%Dimensions to be equal to 
    ! the value of MLSAuxData%rank. I'm not aware of constraints regarding
    ! FrequencyCoordinates, VerticalCoordinates, or HorizontalCoordinates.
    ! -haley n-
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type( MLSAuxData_T ), intent(inout) :: MLSAuxData
    type( DataProducts_T ), intent(in) :: dataset
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: i, status, length

    if (associated(MLSAuxData%Dimensions)) &
      call MLSMessage(MLSMSG_Error, ModuleName, & 
         'MLSAuxData%Dimensions is already allocated. ' &
         // 'Quitting to avoid memory leak.')

    if (associated(MLSAuxData%FrequencyCoordinates)) &
      call MLSMessage(MLSMSG_Error, ModuleName, &
         'MLSAuxData%FrequencyCoordinates is already allocated. ' &
         // 'Quitting to avoid memory leak.')

    if (associated(MLSAuxData%VerticalCoordinates)) &
      call MLSMessage(MLSMSG_Error, ModuleName, &
         'MLSAuxData%VerticalCoordinates is already allocated. ' &
         // 'Quitting to avoid memory leak.')

    if (associated(MLSAuxData%HorizontalCoordinates)) &
      call MLSMessage(MLSMSG_Error, ModuleName, &
         'MLSAuxData%HorizontalCoordinates is already allocated. ' &
         // 'Quitting to avoid memory leak.')

    if (associated(dataset%Dimensions)) then 
       allocate(MLSAuxData%Dimensions(MLSAuxData%rank), stat=status)
       addr = 0
       if ( status == 0 ) addr = transfer(c_loc(MLSAuxData%Dimensions(1)), addr)
       call test_allocate ( status, ModuleName, "MLSAuxData%Dimensions", &
         & uBounds = MLSAuxData%rank, &
         & elementSize = storage_size(MLSAuxData%Dimensions ) / 8, address=addr )
       !MLSAuxData%rank should be less than or equal to the size of dataset%Dimensions,
       !but just to be on the safe size, I'll get the min size of the 2. -haley n-
        do i = 1, min(size(dataset%Dimensions), size(MLSAuxData%Dimensions))
           MLSAuxData%Dimensions(i) = trim(dataset%Dimensions(i))
        end do
    end if

    if (associated(dataset%FrequencyCoordinates)) then
       length = size(dataset%FrequencyCoordinates)
       allocate(MLSAuxData%FrequencyCoordinates(length), stat=status)
       addr = 0
       if ( status == 0 ) then
         if ( length>0 ) addr = transfer(c_loc(MLSAuxData%FrequencyCoordinates(1)), addr)
       end if
       call test_allocate ( status, ModuleName, "MLSAuxData%FrequencyCoordinates", &
         & uBounds = length, &
         & elementSize = storage_size(MLSAuxData%FrequencyCoordinates) / 8, &
         & address=addr )
       MLSAuxData%FrequencyCoordinates(1:length) = &
           & dataset%FrequencyCoordinates(1:length)
    end if

    if (associated(dataset%VerticalCoordinates)) then
       length = size(dataset%VerticalCoordinates)
       allocate(MLSAuxData%VerticalCoordinates(length), stat=status)
       addr = 0
       if ( status == 0 ) then
         if ( length>0 ) addr = transfer(c_loc(MLSAuxData%VerticalCoordinates(1)), addr)
       end if
       call test_allocate ( status, ModuleName, "MLSAuxData%VerticalCoordinates", &
         & uBounds = length, &
         & elementSize = storage_size(MLSAuxData%VerticalCoordinates) / 8, &
         & address=addr )
       MLSAuxData%VerticalCoordinates(1:length) = &
           & dataset%VerticalCoordinates(1:length)
    end if

    if (associated(dataset%HorizontalCoordinates)) then
       length = size(dataset%HorizontalCoordinates)
       allocate(MLSAuxData%HorizontalCoordinates(length), stat=status)
       addr = 0
       if ( status == 0 ) then
         if ( length>0 ) addr = transfer(c_loc(MLSAuxData%HorizontalCoordinates(1)), addr)
       end if
       call test_allocate ( status, ModuleName, "MLSAuxData%HorizontalCoordinates", &
         & uBounds = length, &
         & elementSize = storage_size(MLSAuxData%HorizontalCoordinates) / 8 )
       MLSAuxData%HorizontalCoordinates(1:length) = &
           & dataset%HorizontalCoordinates(1:length)
    end if

  end subroutine CopyFromDataProducts

  subroutine CopyToDataProducts( MLSAuxData, dataset )
    type( MLSAuxData_T ), intent(in) :: MLSAuxData
    type( DataProducts_T ), intent(inout) :: dataset
    integer :: i

    if (associated(MLSAuxData%Dimensions) ) then
      do i = 1, size(MLSAuxData%Dimensions)
        dataset%Dimensions(i) = trim(MLSAuxData%Dimensions(i))
      end do
    end if

    if (associated(MLSAuxData%FrequencyCoordinates) ) then
      i = size(MLSAuxData%FrequencyCoordinates)
      dataset%FrequencyCoordinates(1:i) = MLSAuxData%FrequencyCoordinates(1:i)
    end if

    if (associated(MLSAuxData%VerticalCoordinates) ) then
      i = size(MLSAuxData%VerticalCoordinates)
      dataset%VerticalCoordinates(1:i) = MLSAuxData%VerticalCoordinates(1:i)
    end if

    if (associated(MLSAuxData%HorizontalCoordinates) ) then
      i = size(MLSAuxData%HorizontalCoordinates)
      dataset%HorizontalCoordinates(1:i) = MLSAuxData%HorizontalCoordinates(1:i)
    end if

  end subroutine CopyToDataProducts

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSAuxData

! $Log$
! Revision 2.37  2018/04/19 02:00:36  vsnyder
! Compute address for allocate/deallocate tracking.  Remove USE statements for
! unused names.
!
! Revision 2.36  2015/04/29 00:52:16  vsnyder
! Allocate MLSData%IntField before attempting to compute its address
!
! Revision 2.35  2015/03/28 01:12:53  vsnyder
! Some spiffing.
! Added stuff to trace allocate/deallocate addresses -- mostly commented out
! because NAG build 1017 doesn't yet allow arrays as arguments to C_LOC.
!
! Revision 2.34  2015/01/21 19:27:07  pwagner
! Repaired faulty deallocation
!
! Revision 2.33  2014/09/05 00:00:00  vsnyder
! More complete and accurate allocate/deallocate size tracking.
! Convert some local pointer temps to allocatable.
!
! Revision 2.32  2014/03/07 19:19:24  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 2.31  2010/03/25 18:41:02  pwagner
! args to max function now the same integer type
!
! Revision 2.30  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.29  2009/08/21 19:45:05  pwagner
! Reverse effects of hdf5.1.8 that led to consinsistent noMAFs between datasets
!
! Revision 2.28  2009/07/30 20:43:05  honghanh
! Changed CopyFromDataProduct to create arrays before copying to them
! Add array size checking for Dimensions field in CopyFromDataProduct
!
! Revision 2.27  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.26  2007/04/03 20:51:16  pwagner
! Made dims an assumed-shape array in Build_MLSAuxData
!
! Revision 2.25  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.24  2004/10/07 18:51:12  perun
! Add new optional argument to control writing attributes
!
! Revision 2.23  2004/09/02 20:45:32  perun
! Added more fill_value types
!
! Revision 2.22  2003/09/03 19:46:30  perun
! Add fill value for Reals
!
! Revision 2.21  2003/01/02 23:14:50  pwagner
! Repaired misplaced ampersand on line 1812
!
! Revision 2.20  2002/12/20 20:45:57  perun
! Added MaxCharFieldLen for MLSAuxData_T
!
! Revision 2.19  2002/12/18 16:24:43  perun
! Remove memory leaks due to nullify
!
! Revision 2.18  2002/12/06 23:38:12  pwagner
! Helpful pre-nullifying of pointers
!
! Revision 2.17  2002/11/25 05:35:35  jdone
! Added dimension names.
!
! Revision 2.16  2002/11/11 19:20:25  jdone
! name in DataProducts_T now has length of 80
!
! Revision 2.15  2002/10/31 22:12:49  jdone
! Added Build_MLSAuxData.
!
! Revision 2.14  2002/10/24 03:19:00  jdone
! subroutine Allocate_MLSAuxData added.
!
! Revision 2.13  2002/10/21 22:10:32  jdone
! inclusion of string types for i/o
!
! Revision 2.12  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.11  2002/10/05 01:03:08  jdone
! Create_MLSAuxData updated to reduce memory.
!
! Revision 2.10  2002/10/05 00:52:50  jdone
! CreateGroup_MLSAuxData added.
!
! Revision 2.9  2002/10/05 00:30:20  jdone
! Write_MLSAuxData creates & writes dataset if not already created.
!
! Revision 2.8  2002/10/04 22:35:31  jdone
! Added exception handling for h5tequal_f call
!
! Revision 2.7  2002/10/04 22:13:30  jdone
! Replace == for types with h5tequal_f
!
! Revision 2.6  2002/10/03 22:15:37  jdone
! check hdf5 error flags that return rank
!
! Revision 2.5  2002/09/27 23:37:54  pwagner
! More progress toward hdf5-capable l1b files
!
! Revision 2.4  2002/09/26 23:58:04  pwagner
! Moved attributes out of create function; standalone attribute io, too
!
! Revision 2.3  2002/09/09 05:43:39  jdone
! deallocate statements added for read buffers.
!
