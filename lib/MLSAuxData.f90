! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MLSAuxData

  ! Reading and interacting with Level 1B data (HDF5)
  use MLS_DataProducts, only: Deallocate_DataProducts, DataProducts_T
  use HDF5, only: hid_t, hsize_t, H5T_NATIVE_CHARACTER, H5T_NATIVE_DOUBLE, &
       H5T_NATIVE_INTEGER, H5T_IEEE_F32LE, H5S_UNLIMITED_F, H5T_STD_I32LE, & 
       H5T_NATIVE_REAL,H5T_IEEE_F64LE,H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, &
       h5dopen_f, h5dclose_f, h5screate_simple_f, h5pcreate_f, h5dcreate_f, &
       h5pset_chunk_f, h5sclose_f, h5dget_space_f, h5pclose_f, & 
       h5sget_simple_extent_ndims_f, h5sget_simple_extent_dims_f, &
       h5dget_create_plist_f, h5pget_chunk_f, h5dget_type_f, & 
       h5sselect_hyperslab_f, h5dread_f, h5dwrite_f, h5dextend_f, &
       h5acreate_f, h5awrite_f, h5aread_f, h5aclose_f, h5tcopy_f, &
       h5tset_size_f, h5aopen_name_f, h5aget_type_f, h5aget_space_f, &
       h5tequal_f, h5fis_hdf5_f, h5eset_auto_f, h5gcreate_f, h5gclose_f, &
       h5gopen_f
  use MLSCommon, only: r4, r8
  use MLSStrings, only: Lowercase
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_Error, MLSMSG_deallocate, &
       MLSMSG_allocate
  
  implicit NONE

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!     (data types and parameters)
! MLSAuxData_T                   Quantities from an L1B/L2 data file
! NAME_LEN                       Max length of l1b/l2 sds array name
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
  public :: name_len, MLSAuxData_T, Create_MLSAuxData, Read_MLSAuxData, & 
   & Deallocate_MLSAuxData, Write_MLSAuxData, Read_MLSAuxAttributes, &
   & Write_MLSAuxAttributes, CreateGroup_MLSAuxData, Allocate_MLSAuxData, &
   & Build_MLSAuxData, AddMLSAuxDataToDatabase, DestroyMLSAuxDataDatabase
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------
  !
  ! Parameters
  integer, parameter :: name_len = 64  ! Max len of SDS array name
  !
  type MLSAuxData_T
    character (len=name_len) :: name, type_name  
    character(len=20),  dimension(:), pointer :: Dimensions => NULL()
    !
    character(len=27), dimension(:,:,:), pointer :: CharField => NULL()
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

    type (MLSAuxData_T), dimension(:), pointer :: database
    type (MLSAuxData_T), intent(in) :: item

    ! Local variables
    type (MLSAuxData_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddMLSAuxDataToDatabase = newSize

  end function AddMLSAuxDataToDatabase
!-----------------------------------------------------DestroyMLSAuxDataDatabase
  subroutine DestroyMLSAuxDataDatabase ( database )
  ! Deallocates the elements of array of MLSAuxData_T and then the array.
    type (MLSAuxData_T), dimension(:), pointer :: database
    ! Local variables
    integer :: index, status

    if ( associated(database) ) then
      do index = 1, SIZE(database)
          call  Deallocate_MLSAuxData( database(index) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "database" )
    end if

  end subroutine DestroyMLSAuxDataDatabase
!---------------------------------------------------------- Allocate_MLSAuxData
 subroutine Allocate_MLSAuxData( name, data_type, dims, MLSData  )
    ! This should be called when allocating a MLSAuxData structure.
    type( MLSAuxData_T ), intent(inout) :: MLSData
    integer, dimension(3), intent(in) :: dims
    character (len=*), intent(in) :: name, data_type
    ! internal variables
    character(len=480) :: msr
    integer :: status

     MLSData%name        = trim(name) 
     MLSData%type_name   = trim(data_type)

    if ( data_type .eq. 'real') then
     allocate(MLSData%RealField(dims(1),dims(2),dims(3)), stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_allocate // ' MLSAuxData%RealField in ' // name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if ( data_type .eq. 'double') then
     allocate(MLSData%DpField(dims(1),dims(2),dims(3)), stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_allocate // ' MLSAuxData%DpField in ' // name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if ( data_type .eq. 'integer') then
     allocate(MLSData%IntField(dims(1),dims(2),dims(3)), stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_allocate // ' MLSAuxData%IntField in ' // name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if ( data_type .eq. 'character') then
     allocate(MLSData%CharField(dims(1),dims(2),dims(3)), stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_allocate // ' MLSAuxData%CharField in ' // name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

 end subroutine Allocate_MLSAuxData
!---------------------------------------------------------Deallocate_MLSAuxData
 subroutine Deallocate_MLSAuxData( MLSAuxData )
    ! This should be called when deallocating a MLSAuxData structure.
    type( MLSAuxData_T ), intent(inout) :: MLSAuxData
    character(len=480) :: msr
    integer :: status

    if (associated(MLSAuxData%RealField)) then
        deallocate(MLSAuxData%RealField, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' MLSAuxData%RealField in ' // & 
           MLSAuxData%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if (associated(MLSAuxData%IntField)) then 
        deallocate(MLSAuxData%IntField, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' MLSAuxData%IntField in ' // &
           MLSAuxData%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if (associated(MLSAuxData%CharField)) then 
        deallocate(MLSAuxData%CharField, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' MLSAuxData%CharField in ' // &
           MLSAuxData%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if (associated(MLSAuxData%DpField)) then 
        deallocate(MLSAuxData%DpField, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' MLSAuxData%DpField in ' // &
           MLSAuxData%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if (associated(MLSAuxData%FrequencyCoordinates)) then 
        deallocate(MLSAuxData%FrequencyCoordinates, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' MLSAuxData%FrequencyCoordinates in ' // &
           MLSAuxData%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if (associated(MLSAuxData%VerticalCoordinates)) then 
        deallocate(MLSAuxData%VerticalCoordinates, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' MLSAuxData%VerticalCoordinates in ' // &
           MLSAuxData%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if (associated(MLSAuxData%HorizontalCoordinates)) then 
        deallocate(MLSAuxData%HorizontalCoordinates, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' MLSAuxData%HorizontalCoordinates in ' // &
           MLSAuxData%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

    if (associated(MLSAuxData%Dimensions)) then 
        deallocate(MLSAuxData%Dimensions, stat=status)
     if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' MLSAuxData%Dimensions in ' // &
           MLSAuxData%name
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
    endif

 end subroutine Deallocate_MLSAuxData
!--------------------------------------------------------------Build_MLSAuxData
 subroutine Build_MLSAuxData_Character( file_id, dataset, char_data, & 
      char_length, lastIndex)
    type( DataProducts_T ), intent(in) :: dataset
    character (len=*), intent(in) :: char_data
    integer, intent(in) :: char_length
    integer, intent(in), optional ::lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error, status

    dims(1) = 1
    dims(2) = 1
    dims(3) = 1

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),&
          trim(dataset%data_type),dims,MLSData)

      MLSData%CharField(1,1,1) = char_data
      MLSData%rank = 1

      allocate(MLSData%Dimensions(1), stat=status)

      call CopyFromDataProducts(dataset, MLSData)

      if (present(lastIndex)) then 
         if (lastIndex .eq. 1) then 
     call Write_MLSAuxData(file_id, MLSData, error, & 
          write_attributes=.true., & 
          string_length=char_length, index=lastIndex)
         else
     call Write_MLSAuxData(file_id, MLSData, error, & 
          write_attributes=.false., & 
          string_length=char_length, index=lastIndex)
        endif
     else
     call Write_MLSAuxData(file_id, MLSData, error, & 
          write_attributes=.true., & 
          string_length=char_length)
     endif

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Character
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Integer( file_id, dataset, int_data, & 
      lastIndex)
    type( DataProducts_T ), intent(in) :: dataset
    integer, intent(in) :: int_data
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error, status

    dims(1) = 1
    dims(2) = 1
    dims(3) = 1

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dims,MLSData)

     MLSData%IntField(1,1,1) = int_data
     MLSData%rank = 1

      allocate(MLSData%Dimensions(1), stat=status)
      call CopyFromDataProducts(dataset, MLSData)

     if ( present (lastIndex)) then 
        if (lastIndex .eq. 1) then 
        call Write_MLSAuxData(file_id, MLSData, error, &
          write_attributes=.true., &  
             index=lastIndex)
        else
        call Write_MLSAuxData(file_id, MLSData, error, &
          write_attributes=.false., &  
             index=lastIndex)
        endif
     else

        call Write_MLSAuxData(file_id, MLSData, error, &
          write_attributes=.true.)

     endif

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Integer
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Real( file_id, dataset, real_data, lastIndex)
    type( DataProducts_T ), intent(inout) :: dataset
    real, intent(in) :: real_data
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error, status

    dims(1) = 1
    dims(2) = 1
    dims(3) = 1

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dims,MLSData)

     MLSData%RealField(1,1,1) = real_data
     MLSData%rank = 1
      allocate(MLSData%Dimensions(1), stat=status)
      call CopyFromDataProducts(dataset, MLSData)

     if ( present(lastIndex) ) then
        if (lastIndex .eq. 1) then 
        call Write_MLSAuxData(file_id, MLSData, error, &
          write_attributes=.true., index=lastIndex)
        else
        call Write_MLSAuxData(file_id, MLSData, error, &
          write_attributes=.false., index=lastIndex)
        endif
     else
        call Write_MLSAuxData(file_id, MLSData, error)
     endif

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Real
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Double( file_id, dataset, double_data, lastIndex)
    type( DataProducts_T ), intent(in) :: dataset
    real(r8), intent(in) :: double_data
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error, status

    dims(1) = 1
    dims(2) = 1
    dims(3) = 1

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dims,MLSData)
     MLSData%DpField(1,1,1) = double_data
     MLSData%rank = 1


      allocate(MLSData%Dimensions(1), stat=status)
      call CopyFromDataProducts(dataset, MLSData)

     if ( present (lastIndex) ) then 
        if (lastIndex .eq. 1) then 
           call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
                write_attributes=.true.)
        else
           call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
          write_attributes=.false.)
        endif
     else
        call Write_MLSAuxData(file_id, MLSData, error,&           
             write_attributes=.true.)
     endif

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Double
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Real_1d( file_id, dataset, real_data, lastIndex, &
      dims)
    type( DataProducts_T ), intent(in) :: dataset
    real, dimension(:), intent(in) :: real_data
    integer, dimension(3), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error, status

     if (present(dims) ) then        
          do i=1,3
             dim_array(i) = dims(i)
          enddo
     else
          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(real_data))
         dim_array(i) = size(real_data,i)
        enddo
     endif

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

        do i = 1, dim_array(1)
          MLSData%RealField(i,1,1) = real_data(i) 
        enddo

     if ( present(lastIndex)) then 
             MLSData%rank = 2
     else
             MLSData%rank = 1
     endif

          allocate(MLSData%Dimensions(MLSData%rank), stat=status)
          call CopyFromDataProducts(dataset, MLSData)

     if ( present(lastIndex)) then
        if (lastIndex .eq. 1) then 
        call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
          write_attributes=.true.)
        else
        call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
          write_attributes=.false.)
        endif
     else
        call Write_MLSAuxData(file_id, MLSData, error,&           
             write_attributes=.true.) 
     endif

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Real_1d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Double_1d( file_id, dataset, double_data, & 
      lastIndex, dims)
    type( DataProducts_T ), intent(in) :: dataset
    real(r8), dimension(:), intent(in) :: double_data
    integer, dimension(3), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error, status

     if (present(dims) ) then        
          do i=1,3
             dim_array(i) = dims(i)
          enddo
     else
          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(double_data))
         dim_array(i) = size(double_data,i)
        enddo
     endif

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)
    
        do i = 1, dim_array(1)
          MLSData%DpField(i,1,1) = double_data(i) 
       enddo

     if ( present (lastIndex) ) then 
             MLSData%rank = 2
     else
             MLSData%rank = 1
     endif

     allocate(MLSData%Dimensions(MLSData%rank), stat=status)
     call CopyFromDataProducts(dataset, MLSData)

     if ( present (lastIndex) ) then 
        if (lastIndex .eq. 1) then 
           call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
                       write_attributes=.true.)
        else 
        call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
                       write_attributes=.false.)
        endif
     else
        call Write_MLSAuxData(file_id, MLSData, error, &
          write_attributes=.true.)
     endif

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Double_1d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Integer_1d( file_id, dataset, integer_data, &
      lastIndex, dims)
    type( DataProducts_T ), intent(in) :: dataset
    integer, dimension(:), intent(in) :: integer_data
    integer, dimension(3), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,error, status

     if (present(dims) ) then        
          do i=1,3
             dim_array(i) = dims(i)
          enddo
     else
          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(integer_data))
         dim_array(i) = size(integer_data,i)
        enddo
     endif

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

        do i = 1, dim_array(1)
          MLSData%IntField(i,1,1) = integer_data(i) 
       enddo

     if (present (lastIndex) ) then
             MLSData%rank = 2 
     else
             MLSData%rank = 1
     endif

      allocate(MLSData%Dimensions(MLSData%rank), stat=status)
      call CopyFromDataProducts(dataset, MLSData)

     if (present (lastIndex) ) then
        if (lastIndex .eq. 1) then 
        call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
          write_attributes=.true.)
        else
        call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
          write_attributes=.false.)
        endif
     else
        call Write_MLSAuxData(file_id, MLSData, error, &
          write_attributes=.true.)
     endif

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Integer_1d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Real_2d( file_id, dataset, real_data, lastIndex, &
      dims)
    type( DataProducts_T ), intent(in) :: dataset
    real, dimension(:,:), intent(in) :: real_data
    integer, dimension(3), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,k,error, status

     if (present(dims) ) then        
          do i=1,3
             dim_array(i) = dims(i)
          enddo
     else
          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(real_data))
         dim_array(i) = size(real_data,i)
        enddo
     endif

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     do j = 1, dim_array(2)
        do i = 1, dim_array(1)
          MLSData%RealField(i,j,1) = real_data(i,j) 
       enddo
     enddo

     if (present (lastIndex) ) then 
             MLSData%rank = 3
     else
             MLSData%rank = 2
     endif

      allocate(MLSData%Dimensions(MLSData%rank), stat=status)
      call CopyFromDataProducts(dataset, MLSData)

     if (present (lastIndex) ) then 
        if (lastIndex .eq. 1) then 
        call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
                       write_attributes=.true.)
        else
        call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
                       write_attributes=.false.)
        endif
     else
        call Write_MLSAuxData(file_id, MLSData, error, &
          write_attributes=.true.)
     endif

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Real_2d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Double_2d( file_id, dataset, double_data, & 
      lastIndex, dims)
    type( DataProducts_T ), intent(in) :: dataset
    real(r8), dimension(:,:), intent(in) :: double_data
    integer, dimension(3), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,error, status

     if (present(dims) ) then 
          do i=1,3
             dim_array(i) = dims(i)
          enddo
     else
          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(double_data))
         dim_array(i) = size(double_data,i)
        enddo
     endif

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)


     do j = 1, dim_array(2)
        do i = 1, dim_array(1)
          MLSData%DpField(i,j,1) = double_data(i,j) 
       enddo
     enddo

     if (present (lastIndex) ) then 
             MLSData%rank = 3
     else
             MLSData%rank = 2
     endif

      allocate(MLSData%Dimensions(MLSData%rank), stat=status)
      call CopyFromDataProducts(dataset, MLSData)

     if (present (lastIndex) ) then
        if (lastIndex .eq. 1) then 
        call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
          write_attributes=.true.)
        else
        call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
          write_attributes=.false.)
        endif
     else
        call Write_MLSAuxData(file_id, MLSData, error, &
          write_attributes=.true.)
     endif

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Double_2d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Integer_2d( file_id, dataset, integer_data, &
      lastIndex, dims)
    type( DataProducts_T ), intent(in) :: dataset
    integer, dimension(:,:), intent(in) :: integer_data
    integer, dimension(3), intent(in), optional :: dims
    integer, intent(in), optional :: lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,error, status

     if (present(dims) ) then 
          do i=1,3
             dim_array(i) = dims(i)
          enddo
     else
          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(integer_data))
         dim_array(i) = size(integer_data,i)
        enddo
     endif

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

      do j = 1, dim_array(2)
        do i = 1, dim_array(1)
          MLSData%IntField(i,j,1) = integer_data(i,j) 
       enddo
      enddo

      if (present (lastIndex)) then 
             MLSData%rank = 3
      else
             MLSData%rank = 2
      endif

      allocate(MLSData%Dimensions(MLSData%rank), stat=status)
      call CopyFromDataProducts(dataset, MLSData)

      if (present (lastIndex)) then
         if (lastIndex .eq. 1) then 
         call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
          write_attributes=.true.)
         else
         call Write_MLSAuxData(file_id, MLSData, error, index=lastIndex, &
          write_attributes=.false.)
         endif
      else
         call Write_MLSAuxData(file_id, MLSData, error, &
              write_attributes=.true.)
      endif

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Integer_2d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Real_3d( file_id, dataset, real_data, &
      dims)
    type( DataProducts_T ), intent(in) :: dataset
    real, dimension(:,:,:), intent(in) :: real_data
    integer, dimension(3), intent(in), optional :: dims
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,k,error, status

     if (present(dims) ) then 
          do i=1,3
             dim_array(i) = dims(i)
          enddo
     else
          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(real_data))
         dim_array(i) = size(real_data,i)
        enddo
     endif

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     MLSData%rank = 3

     do k = 1, dim_array(3)
      do j = 1, dim_array(2)
        do i = 1, dim_array(1)
          MLSData%RealField(i,j,k) = real_data(i,j,k) 
       enddo
      enddo
     enddo

      allocate(MLSData%Dimensions(3), stat=status)
      call CopyFromDataProducts(dataset, MLSData)

     call Write_MLSAuxData(file_id, MLSData, error,write_attributes=.true.)

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Real_3d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Double_3d( file_id, dataset, double_data, &
      dims)
    type( DataProducts_T ), intent(in) :: dataset
    real(r8), dimension(:,:,:), intent(in) :: double_data
    integer, dimension(3), intent(in), optional :: dims
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,k,error, status

     if (present(dims) ) then 
          do i=1,3
             dim_array(i) = dims(i)
          enddo
     else
          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(double_data))
         dim_array(i) = size(double_data,i)
        enddo
     endif

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     MLSData%rank = 3

     do k = 1, dim_array(3)
      do j = 1, dim_array(2)
        do i = 1, dim_array(1)
          MLSData%DpField(i,j,k) = double_data(i,j,k) 
       enddo
      enddo 
     enddo

      allocate(MLSData%Dimensions(3), stat=status)
      call CopyFromDataProducts(dataset, MLSData)

     call Write_MLSAuxData(file_id, MLSData, error,write_attributes=.true.)
    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData.' ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Double_3d
!------------------------------------------------------------------------------
 subroutine Build_MLSAuxData_Integer_3d(file_id,dataset,integer_data,dims)
    type( DataProducts_T ), intent(in) :: dataset
    integer, dimension(:,:,:), intent(in) :: integer_data
    integer, dimension(3), intent(in), optional :: dims
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,k,error, status

     if (present(dims) ) then 
          do i=1,3
             dim_array(i) = dims(i)
          enddo
     else
          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(integer_data))
         dim_array(i) = size(integer_data,i)
        enddo
     endif

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     MLSData%rank = 3

     do k = 1, dim_array(3)
      do j = 1, dim_array(2)
        do i = 1, dim_array(1)
          MLSData%IntField(i,j,k) = integer_data(i,j,k) 
       enddo
      enddo
     enddo

      allocate(MLSData%Dimensions(3), stat=status)
      call CopyFromDataProducts(dataset, MLSData)

     call Write_MLSAuxData(file_id, MLSData, error,write_attributes=.true.)
    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData.' ) 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Build_MLSAuxData_Integer_3d
  !------------------------------------------------------- Recall_MLSAuxData
 subroutine Recall_MLSAuxData_Character( file_id, dataset, char_data, &
      char_length)
    type( DataProducts_T ), intent(inout) :: dataset
    character (len=*), intent(inout) :: char_data
    integer, intent(inout) :: char_length
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error, status

    dims(1) = 1
    dims(2) = 1
    dims(3) = 1

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),&
          trim(dataset%data_type),dims,MLSData)

     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type), MLSData, error, read_attributes=.true.)

    if (error .eq. 0) then 
      char_data = MLSData%CharField(1,1,1) 

      allocate(dataset%Dimensions(1), stat=status)
      call CopyToDataProducts(MLSData, dataset)

    else 
       call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) )
    endif
 
     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Character
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Integer( file_id, dataset, int_data )
    type( DataProducts_T ), intent(inout) :: dataset
    integer, intent(inout) :: int_data
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error, status

    dims(1) = 1
    dims(2) = 1
    dims(3) = 1

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dims,MLSData)

     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type),MLSData, error, read_attributes=.true.)

    if (error .eq. 0) then 
     int_data = MLSData%IntField(1,1,1)

      allocate(dataset%Dimensions(1), stat=status)
      call CopyToDataProducts(MLSData, dataset)

    else 
       call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) ) 
    endif

     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Integer
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Real(file_id,dataset,real_data)
    type( DataProducts_T ), intent(inout) :: dataset
    real, intent(inout) :: real_data
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error, status

    dims(1) = 1
    dims(2) = 1
    dims(3) = 1

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dims,MLSData)

     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type),MLSData, error, read_attributes=.true.)
    if (error .eq. 0) then 
        real_data = MLSData%RealField(1,1,1)

      allocate(dataset%Dimensions(1), stat=status)
      call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) )
    endif
     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Real
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Double( file_id, dataset, double_data)
    type( DataProducts_T ), intent(inout) :: dataset
    real(r8), intent(inout) :: double_data
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dims
    integer :: error, status

    dims(1) = 1
    dims(2) = 1
    dims(3) = 1

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dims,MLSData)
     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type),MLSData, error, read_attributes=.true.)

    if (error .eq. 0) then 
     double_data = MLSData%DpField(1,1,1)
      allocate(dataset%Dimensions(1), stat=status)
      call CopyToDataProducts(MLSData, dataset)
    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) )
    endif

     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Double
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Real_1d( file_id, dataset, real_data, & 
      firstIndex, lastIndex)
    type( DataProducts_T ), intent(inout) :: dataset
    real, dimension(:), intent(inout) :: real_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error, status, i_first, i_last

        do i=1,3
            dim_array(i) = 1
        enddo

        do i=1,size(shape(real_data))
         dim_array(i) = size(real_data,i)
        enddo

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)
     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type),MLSData, error, read_attributes=.true.)

    if (error .eq. 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          i_first = firstIndex
          i_last  = lastIndex
       else
          i_first = 1
          i_last  = dim_array(1)
       endif

          allocate(dataset%Dimensions(MLSData%rank), stat=status)
          call CopyToDataProducts(MLSData, dataset)

        do i = i_first, i_last
          real_data(i) = MLSData%RealField(i,1,1)  
        enddo

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) )
    endif

     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Real_1d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Double_1d(file_id, dataset, double_data, & 
      firstIndex, lastIndex)
    type( DataProducts_T ), intent(inout) :: dataset
    real(r8), dimension(:), intent(inout) :: double_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,error, status, i_first, i_last

        do i=1,3
            dim_array(i) = 1
        enddo

        do i=1,size(shape(double_data))
         dim_array(i) = size(double_data,i)
        enddo

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type),MLSData, error, read_attributes=.true.)

     if (error /= 0) then 

       if (present (firstIndex) .and. present (lastIndex) ) then
          i_first = firstIndex
          i_last  = lastIndex
       else
          i_first = 1
          i_last  = dim_array(1)
       endif

          allocate(dataset%Dimensions(MLSData%rank), stat=status)
          call CopyToDataProducts(MLSData, dataset)

        do i = i_first, i_last
          double_data(i) = MLSData%DpField(i,1,1) 
        enddo


     else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) )
     endif

     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Double_1d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Integer_1d(file_id,dataset,integer_data, &
      firstIndex, lastIndex)
    type( DataProducts_T ), intent(inout) :: dataset
    integer, dimension(:), intent(inout) :: integer_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,error,status,i_first,i_last

          do i=1,3
             dim_array(i) = 1
          enddo

        do i=1,size(shape(integer_data))
         dim_array(i) = size(integer_data,i)
        enddo

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type),MLSData, error, read_attributes=.true.)

     if (error .eq. 0) then
       if (present (firstIndex) .and. present (lastIndex) ) then
          i_first = firstIndex
          i_last  = lastIndex
       else
          i_first = 1
          i_last  = dim_array(1)
       endif
        do i = i_first, i_last
          integer_data(i) = MLSData%IntField(i,1,1)  
        enddo

      allocate(dataset%Dimensions(MLSData%rank), stat=status)
      call CopyToDataProducts(MLSData, dataset)

     else
        call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) )
     endif

     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Integer_1d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Real_2d( file_id, dataset, real_data, & 
      firstIndex, lastIndex)
    type( DataProducts_T ), intent(inout) :: dataset
    real, dimension(:,:), intent(inout) :: real_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,k,error,status,j_first,j_last

          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(real_data))
         dim_array(i) = size(real_data,i)
        enddo

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type),MLSData, error, read_attributes=.true.)
     if (error .eq. 0) then 

       if (present (firstIndex) .and. present (lastIndex) ) then
          j_first = firstIndex
          j_last  = lastIndex
       else
          j_first = 1
          j_last  = dim_array(2)
       endif

      do j = j_first, j_last
        do i = 1, dim_array(1)
          real_data(i,j) = MLSData%RealField(i,j,1) 
       enddo
      enddo

      allocate(dataset%Dimensions(MLSData%rank), stat=status)
      call CopyToDataProducts(MLSData, dataset)

     else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) )
     endif

     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Real_2d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Double_2d( file_id, dataset, double_data,& 
      firstIndex, lastIndex)
    type( DataProducts_T ), intent(inout) :: dataset
    real(r8), dimension(:,:), intent(inout) :: double_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,error,status,j_first,j_last

          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(double_data))
         dim_array(i) = size(double_data,i)
        enddo

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type), MLSData, error, read_attributes=.true.)
    if (error .eq. 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          j_first = firstIndex
          j_last  = lastIndex
       else
          j_first = 1
          j_last  = dim_array(2)
       endif

     do j = j_first, j_last
        do i = 1, dim_array(1)
          double_data(i,j) = MLSData%DpField(i,j,1) 
       enddo
     enddo

      allocate(dataset%Dimensions(MLSData%rank), stat=status)
      call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) )
    endif
     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Double_2d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Integer_2d( file_id, dataset, integer_data, &
      firstIndex, lastIndex)
    type( DataProducts_T ), intent(inout) :: dataset
    integer, dimension(:,:), intent(inout) :: integer_data
    integer, intent(in), optional :: firstIndex, lastIndex
    integer(hid_t), intent(in) :: file_id

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,error, status,j_first, j_last

          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(integer_data))
         dim_array(i) = size(integer_data,i)
        enddo

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type), MLSData, error, read_attributes=.true.)
    if (error .eq. 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          j_first = firstIndex
          j_last  = lastIndex
       else
          j_first = 1
          j_last  = dim_array(2)
       endif

      do j = j_first, j_last
        do i = 1, dim_array(1)
          integer_data(i,j) = MLSData%IntField(i,j,1)  
       enddo
      enddo

      allocate(dataset%Dimensions(MLSData%rank), stat=status)
      call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) )
    endif
     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Integer_2d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Real_3d( file_id, dataset, real_data, & 
      firstIndex, lastIndex)
    type( DataProducts_T ), intent(inout) :: dataset
    real, dimension(:,:,:), intent(inout) :: real_data
    integer(hid_t), intent(in) :: file_id
    integer, intent(in), optional :: firstIndex, lastIndex

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,k,error, status,k_first, k_last

          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(real_data))
         dim_array(i) = size(real_data,i)
        enddo

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type), MLSData, error, read_attributes=.true.)

    if (error .eq. 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          k_first = firstIndex
          k_last  = lastIndex
       else
          k_first = 1
          k_last  = dim_array(3)
       endif
       
     do k = k_first, k_last
      do j = 1, dim_array(2)
        do i = 1, dim_array(1)
          real_data(i,j,k) = MLSData%RealField(i,j,k)  
       enddo
      enddo
     enddo

      allocate(dataset%Dimensions(3), stat=status)
      call CopyToDataProducts(MLSData, dataset)

    else
       call MLSMessage(MLSMSG_Error, ModuleName, & 
       'Error Writing MLSAuxData for '// trim(dataset%name) )
    endif
     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Real_3d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Double_3d( file_id, dataset, double_data, &
      firstIndex, lastIndex)
    type( DataProducts_T ), intent(inout) :: dataset
    real(r8), dimension(:,:,:), intent(inout) :: double_data
    integer(hid_t), intent(in) :: file_id
    integer, intent(in), optional :: firstIndex, lastIndex

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,k,error, status,k_first, k_last

          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(double_data))
         dim_array(i) = size(double_data,i)
        enddo

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)
     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type), MLSData, error, read_attributes=.true.)
    if (error .eq. 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          k_first = firstIndex
          k_last  = lastIndex
       else
          k_first = 1
          k_last  = dim_array(3)
       endif

     do k = k_first, k_last
      do j = 1, dim_array(2)
        do i = 1, dim_array(1)
          double_data(i,j,k) = MLSData%DpField(i,j,k) 
       enddo
      enddo 
     enddo

      allocate(dataset%Dimensions(3), stat=status)
      call CopyToDataProducts(MLSData, dataset)

    else
     call MLSMessage(MLSMSG_Error,ModuleName,'Error Writing MLSAuxData.')
    endif
     call deallocate_mlsauxdata(MLSData)
 end subroutine Recall_MLSAuxData_Double_3d
!------------------------------------------------------------------------------
 subroutine Recall_MLSAuxData_Integer_3d(file_id,dataset,integer_data, &
      firstIndex, lastIndex)
    type( DataProducts_T ), intent(inout) :: dataset
    integer, dimension(:,:,:), intent(inout) :: integer_data
    integer(hid_t), intent(in) :: file_id
    integer, intent(in), optional :: firstIndex, lastIndex

    type( MLSAuxData_T ) :: MLSData
    integer, dimension(3) :: dim_array
    integer :: i,j,k,error,status,k_first, k_last

          do i=1,3
             dim_array(i) = 1
          enddo
        do i=1,size(shape(integer_data))
         dim_array(i) = size(integer_data,i)
        enddo

     call deallocate_mlsauxdata(MLSData)
     call Allocate_MLSAuxData(trim(dataset%name),& 
          trim(dataset%data_type),dim_array,MLSData)

     call Read_MLSAuxData(file_id, trim(dataset%name), &
          trim(dataset%data_type), MLSData, error, read_attributes=.true.)
    if (error .eq. 0) then

       if (present (firstIndex) .and. present (lastIndex) ) then
          k_first = firstIndex
          k_last  = lastIndex
       else
          k_first = 1
          k_last  = dim_array(3)
       endif
!
! dump the data into the array.
!
     do k = k_first, k_last
      do j = 1, dim_array(2)
        do i = 1, dim_array(1)
          integer_data(i,j,k) = MLSData%IntField(i,j,k)  
       enddo
      enddo
     enddo

      allocate(dataset%Dimensions(3), stat=status)
      call CopyToDataProducts(MLSData, dataset)

    else
     call MLSMessage(MLSMSG_Error,ModuleName,'Error Writing MLSAuxData.')
    endif
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

    endif

  end subroutine CreateGroup_MLSAuxData
! -------------------------------------------------  Create_MLSAuxData ----
  subroutine Create_MLSAuxData(file_id, MLSAuxData, string_length, & 
       write_attributes)
!
! This subroutine creates an entry in the HDF5 file.
!----------------------------------------------------------------------
! External variables
!
    type( MLSAuxData_T ), intent(in) :: MLSAuxData
    integer(hid_t), intent(in)       :: file_id ! From HDF
    integer, intent(in), optional    :: string_length
    logical, intent(in), optional    :: write_attributes
!-----------------------------------------------------------------------
! Internal variables
!
    character(len=480) :: msr
    character(len=name_len) :: aname
    real, dimension(:), pointer :: attr_data
    integer(hsize_t), dimension(7) :: adims
    integer(hsize_t), dimension(3) :: chunk_dims, dims, maxdims
    integer(hsize_t), dimension(1) :: adims_create
    integer(hid_t) :: cparms,dspace_id,dset_id,type_id, &
         attr_id, atype_id, aspace_id, s_type_id
    integer :: i, rank, arank, h5error, status
!-----------------------------------------------------------------------
    nullify(attr_data)
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
       endif
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
       endif
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
       endif
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

    if (trim(MLSAuxData%type_name).ne.'character') then

    call h5dcreate_f(file_id,trim(MLSAuxData%name),type_id,dspace_id, &
          dset_id,h5error,cparms)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSET_CREATE // trim(MLSAuxData%name) ) 

    else
 
     call h5tcopy_f(type_id, s_type_id, h5error)
     if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_TYPE_COPY // trim(MLSAuxData%name))

     call h5tset_size_f(s_type_id, string_length, h5error)
     if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_TYPE_SET // trim(MLSAuxData%name))

     call h5dcreate_f(file_id, trim(MLSAuxData%name), s_type_id, &
       dspace_id, dset_id, h5error, cparms)
     if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_DSET_CREATE // trim(MLSAuxData%name))

    endif

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
       if ( status /= 0 ) then
           msr = MLSMSG_Allocate// ' attr_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif

       do i=1, adims(1)
          attr_data(i) = MLSAuxData%HorizontalCoordinates(i)
       end do

       aname = "HorizontalCoordinates"  

       call h5screate_simple_f(arank, adims_create, aspace_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSPACE_CREATE // trim(aname) // ' in ' // & 
       trim(MLSAuxData%name) )

       call h5acreate_f(dset_id, trim(aname), atype_id, aspace_id, attr_id, &
          h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_ATT_CREATE // trim(aname) // ' in ' // &
       trim(MLSAuxData%name) )

       call h5awrite_f(attr_id, atype_id, attr_data, adims, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_ATT_WRITE // trim(MLSAuxData%name)  )

       call h5aclose_f(attr_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_ATT_CLOSE // trim(MLSAuxData%name) )

       call h5sclose_f(aspace_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSPACE_CLOSE // trim(MLSAuxData%name) )

       deallocate(attr_data,STAT=status)
       if ( status /= 0 ) then
           msr = MLSMSG_deallocate// ' attr_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif

    endif

    if (associated(MLSAuxData%VerticalCoordinates) ) then 
       arank = size(shape(MLSAuxData%VerticalCoordinates))
       atype_id = H5T_IEEE_F32LE
       adims(1) = size(MLSAuxData%VerticalCoordinates)
       adims_create(1) = adims(1)

       allocate(attr_data(adims(1)),STAT=status)
       if ( status /= 0 ) then
           msr = MLSMSG_allocate// ' attr_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif

       do i=1,adims(1)
          attr_data(i) = MLSAuxData%VerticalCoordinates(i)
       end do

       call h5screate_simple_f(arank, adims_create, aspace_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSPACE_CREATE // trim(MLSAuxData%name) )

       aname = "VerticalCoordinates"   
       call h5acreate_f(dset_id, trim(aname), atype_id, aspace_id, attr_id, &
          h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_ATT_CREATE // trim(MLSAuxData%name) )

       call h5awrite_f(attr_id, atype_id, attr_data, adims, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_ATT_WRITE // trim(MLSAuxData%name) )

       call h5aclose_f(attr_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_ATT_CLOSE // trim(MLSAuxData%name) )

       call h5sclose_f(aspace_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSPACE_CLOSE // trim(MLSAuxData%name) )
 
       deallocate(attr_data,STAT=status)
       if ( status /= 0 ) then
           msr = MLSMSG_deallocate// ' attr_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif

    endif

    if (associated(MLSAuxData%FrequencyCoordinates)) then 
       arank = size(shape(MLSAuxData%FrequencyCoordinates))
       atype_id = H5T_IEEE_F32LE
       adims(1) = size(MLSAuxData%FrequencyCoordinates)
       adims_create(1) = adims(1)

       allocate(attr_data(adims(1)),STAT=status)
       if ( status /= 0 ) then
           msr = MLSMSG_allocate// ' attr_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif

       do i=1,adims(1)
          attr_data(i) = MLSAuxData%FrequencyCoordinates(i)
       end do

       call h5screate_simple_f(arank, adims_create, aspace_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSPACE_CREATE // trim(MLSAuxData%name) )

       aname = "FrequencyCoordinates"
       call h5acreate_f(dset_id, trim(aname), atype_id, aspace_id, attr_id, &
          h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_ATT_CREATE // trim(MLSAuxData%name) )

       call h5awrite_f(attr_id, atype_id, attr_data, adims, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_ATT_WRITE // trim(MLSAuxData%name) )

       call h5aclose_f(attr_id, h5error)
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_ATT_CLOSE // trim(MLSAuxData%name) )

       call h5sclose_f(aspace_id, h5error) 
       if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSPACE_CLOSE // trim(MLSAuxData%name) )

       deallocate(attr_data,STAT=status)
       if ( status /= 0 ) then
           msr = MLSMSG_deallocate// ' attr_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif
      endif
     endif
    endif

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
    character(len=*), intent(in)                  :: QuantityName
    character(len=*), intent(in)                  :: AttributeName
    integer(hid_t), intent(in)                    :: dset_id ! From h5dopen_f
    integer, intent(out)                          :: error   ! 0 unless trouble
    type( MLSAuxData_T ), intent(inout), optional :: MLSAuxData
    real, dimension(:), intent(out), optional     :: AttributeData
  ! Private
    character(len=480) :: msr
    character, dimension(:), pointer :: char_data
    real, dimension(:), pointer :: attr_data
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
    nullify(char_data, attr_data)
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
        endif
      case ('HorizontalCoordinates')
        if ( associated ( MLSAuxData%HorizontalCoordinates ) ) then
          arank = 1
          adims(1) = size(MLSAuxData%HorizontalCoordinates)
        else
          error = 1
          return
        endif
      case ('VerticalCoordinates')
        if ( associated ( MLSAuxData%VerticalCoordinates ) ) then
          arank = 1
          adims(1) = size(MLSAuxData%VerticalCoordinates)
        else
          error = 1
          return
        endif
      case ('Dimensions')
        if ( associated (MLSAuxData%Dimensions) ) then
           arank = 1
               is_char  = .true.
               atype_id = H5T_NATIVE_CHARACTER
           do i = 1, arank
              adims(i) = size(MLSAuxData%Dimensions)
           enddo
        else
          error = 1
          return
        endif
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
    if ( status /= 0 ) then
       msr = MLSMSG_allocate// ' char_data.'
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call h5aread_f(attr_id, atype_id, char_data, adims, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
        H5_ERROR_ATT_READ // trim(AttributeName) // ' in ' // &
        QuantityName )

    else

    allocate(attr_data(adims(1)),stat=status)
    if ( status /= 0 ) then
       msr = MLSMSG_allocate// ' attr_data.'
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call h5aread_f(attr_id, atype_id, attr_data, adims, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
        H5_ERROR_ATT_READ // trim(AttributeName) // ' in ' // &
        QuantityName )

    endif

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

    if (associated(attr_data)) then 
    deallocate(attr_data,STAT=status)
     if ( status /= 0 ) then
        msr = MLSMSG_deallocate// ' attr_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
     endif

    if (associated(char_data)) then      
    deallocate(char_data,STAT=status)
     if ( status /= 0 ) then
        msr = MLSMSG_deallocate// ' char_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
     endif
     endif

  end subroutine Read_MLSAuxAttributes
!-------------------------------------------------Write_MLSAuxAttributes ----
  subroutine Write_MLSAuxAttributes(dset_id, QuantityName, AttributeName, & 
       error, MLSAuxData, AttributeData)
  !
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
    character(len=480) :: msr
    character(len=20), dimension(:), pointer :: char_data
    real, dimension(:), pointer :: attr_data
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
    nullify(char_data, attr_data)
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
       if ( status /= 0 ) then
           msr = MLSMSG_allocate// ' attr_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif

        do i=1,adims(1)
           attr_data(i) = MLSAuxData%FrequencyCoordinates(i)
        end do
      case ('HorizontalCoordinates')
        arank = 1
        adims(1) = size(MLSAuxData%HorizontalCoordinates)
        allocate(attr_data(adims(1)),STAT=status)
        if ( status /= 0 ) then
           msr = MLSMSG_allocate// ' attr_data.'
         call MLSMessage(MLSMSG_Error, ModuleName, msr)
        endif

        do i=1,adims(1)
           attr_data(i) = MLSAuxData%HorizontalCoordinates(i)
        end do
      case ('VerticalCoordinates')
        arank = 1
        adims(1) = size(MLSAuxData%VerticalCoordinates)
        allocate(attr_data(adims(1)),STAT=status)
        if ( status /= 0 ) then
           msr = MLSMSG_allocate// ' attr_data.'
         call MLSMessage(MLSMSG_Error, ModuleName, msr)
        endif

        do i=1,adims(1)
           attr_data(i) = MLSAuxData%VerticalCoordinates(i)
        end do
      case ('Dimensions')
        arank = 1
        atype_id = H5T_NATIVE_CHARACTER
        is_char  = .true.
        adims(1) = size(MLSAuxData%Dimensions)
        allocate(char_data(adims(1)),STAT=status)
        if ( status /= 0 ) then
           msr = MLSMSG_allocate// ' char_data.'
         call MLSMessage(MLSMSG_Error, ModuleName, msr)
        endif
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
        if ( status /= 0 ) then
           msr = MLSMSG_allocate // ' attr_data.'
         call MLSMessage(MLSMSG_Error, ModuleName, msr)
        endif

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

     CALL h5tset_size_f(atype_id, attrlen, h5error)
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

    endif

    call h5aclose_f(attr_id, h5error)                                  
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
        H5_ERROR_ATT_CLOSE // trim(AttributeName) // ' in ' // &
        trim(MLSAuxData%name) )

    call h5sclose_f(aspace_id, h5error)                                
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
        H5_ERROR_DSPACE_CLOSE // trim(AttributeName) // ' in ' // & 
    trim(MLSAuxData%name) )

    if (associated(attr_data)) then 
    deallocate(attr_data,STAT=status)
    if ( status /= 0 ) then
        msr = MLSMSG_deallocate// ' attr_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
    endif

    if (associated(char_data)) then 
    deallocate(char_data,STAT=status)
    if ( status /= 0 ) then
        msr = MLSMSG_deallocate// ' char_data.'
        call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
    endif
                                
  end subroutine Write_MLSAuxAttributes
!-------------------------------------------------Read_MLSAuxData ----
  subroutine Read_MLSAuxData(file_id, QuantityName, QuantityType, & 
       MLSAuxData, error, read_attributes, FirstIndex, LastIndex)
!
! This subroutine reads an entry from the HDF5 file.
!
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
    character, dimension(:,:,:), pointer :: character_buffer
    character(len=480) :: msr
    character(len=16) :: myQuantityType
!
    real, dimension(:,:,:), pointer :: real_buffer, double_buffer
!
    integer, dimension(:,:,:), pointer :: integer_buffer
    integer(hsize_t), dimension(7) :: dims
    integer(hsize_t), dimension(3) :: chunk_dims, dims_create, maxdims, start
    integer(hid_t) :: cparms, dset_id, dspace_id, type_id, memspace
    integer        :: rank, h5error, i, j, k, status
!
    logical :: myRead_attributes, is_integer, is_real, is_double, & 
         is_character, is_int32, is_float32, is_float64
!
    error = 0
    myRead_attributes = .false.
    nullify(character_buffer, real_buffer, double_buffer, integer_buffer)
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
      elseif ( is_real .or. is_float32 ) then
        myQuantityType = 'real'
      elseif ( is_double .or. is_float64 ) then
        myQuantityType = 'double'
      elseif ( is_character ) then
        myQuantityType = 'character'
      else
        error = 1
        return
      endif

    endif

    do i = 1, 7
       if (i .le. rank) then 
           dims(i) = dims_create(i) 
       else
           dims(i) = 0
       endif
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
    if ( status /= 0 ) then
      msr = MLSMSG_allocate // ' real_buffer.'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call h5dread_f(dset_id,type_id,& 
         real_buffer,dims,h5error, memspace, dspace_id)

         do k = FirstIndex, LastIndex
            do j = 1, dims(2)
               do i = 1, dims(1)
                  MLSAuxData%RealField(i,j,k) = real_buffer(i,j,k)
               end do 
            end do
         end do

    if (associated(real_buffer)) then 
        nullify(real_buffer)
        deallocate(real_buffer, stat=status)

    if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' real_buffer.'  
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    endif

    case ('double')    

    allocate( double_buffer(dims(1),dims(2),dims(3)),stat=status)
    if ( status /= 0 ) then
      msr = MLSMSG_allocate // ' double_buffer'
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call h5dread_f(dset_id,type_id,& 
         double_buffer,dims,h5error, memspace, dspace_id)
         do k = FirstIndex, LastIndex
            do j = 1, dims(2)
               do i = 1, dims(1)
                  MLSAuxData%DpField(i,j,k) = double_buffer(i,j,k)
               end do 
            end do
         end do    

    if (associated(double_buffer)) then 
        nullify(double_buffer)
        deallocate(double_buffer, stat=status)

    if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' double_buffer.' 
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    endif

    case ('integer')  
    allocate( integer_buffer(dims(1),dims(2),dims(3)),stat=status)
    if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' integer_buffer'  
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call h5dread_f(dset_id,type_id,& 
         integer_buffer,dims,h5error, memspace, dspace_id)
         do k = FirstIndex, LastIndex
            do j = 1, dims(2)
               do i = 1, dims(1)
                  MLSAuxData%IntField(i,j,k) = integer_buffer(i,j,k)
               end do 
            end do
         end do

    if (associated(integer_buffer)) then 
        nullify(integer_buffer)
        deallocate(integer_buffer, stat=status)
    if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' integer_buffer.' 
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
    endif

    case ('character') 
    allocate( character_buffer(dims(1),dims(2),dims(3)),stat=status)

    if ( status /= 0 ) then
      msr = MLSMSG_allocate // ' character_buffer.' 
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    call h5dread_f(dset_id,type_id,& 
         character_buffer,dims,h5error, memspace, dspace_id)
         do k = FirstIndex, LastIndex
            do j = 1, dims(2)
               do i = 1, dims(1)
                  MLSAuxData%CharField(i,j,k) = character_buffer(i,j,k)
               end do 
            end do
         end do

    if (associated(character_buffer)) then 
        nullify(character_buffer)
        deallocate(character_buffer, stat=status)

    if ( status /= 0 ) then
      msr = MLSMSG_deallocate // ' character_buffer.' 
      call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    endif

    end select test_type_dims

    endif

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
   endif
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
  subroutine Write_MLSAuxData(file_id, MLSAuxData, error, &
    & write_attributes, string_length, index)
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
!-------------------------------------------------------------------------
! Internal variables.
!
    character(len=480) :: msr
    character(len=name_len) :: aname
    character(len=1), dimension(:), pointer :: char_data
    real, dimension(:), pointer :: attr_data
    integer(hsize_t), dimension(7) :: adims, dims
    integer(hsize_t), dimension(3) :: chunk_dims, dims_create, maxdims, start
    integer(hsize_t), dimension(1) :: adims_create
    integer(hid_t) :: cparms,dspace_id,dset_id,type_id, &
         attr_id, atype_id, aspace_id, filespace, memspace, s_type_id 
    integer :: i, rank, arank, h5error, status
    logical :: myWrite_attributes
!--------------------------------------------------------------------------

    error = 0
    myWrite_attributes = .false.
    if ( present(write_attributes) ) myWrite_attributes=write_attributes

    nullify(char_data, attr_data)
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
       endif
       endif
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
       endif
       endif
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
       enddo
       endif
       endif
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
       enddo
       endif
       endif
    end select test_type

    call h5eset_auto_f(0, h5error)

    call h5dopen_f(file_id, trim(MLSAuxData%name), dset_id, h5error)

    if (h5error /= 0) then 

    call h5eset_auto_f(1, h5error)

!--------------------------------------------------------------------
    call h5screate_simple_f(rank, dims_create(1:rank),dspace_id,h5error,& 
         maxdims(1:rank))
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSPACE_CREATE // trim(MLSAuxData%name) )

    call h5pcreate_f(H5P_DATASET_CREATE_F,cparms,h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_PROPERTY_CREATE // trim(MLSAuxData%name) )

    call h5pset_chunk_f(cparms,rank,chunk_dims(1:rank),h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_PROPERTY_CHUNK_SET // trim(MLSAuxData%name) )

    if (trim(MLSAuxData%type_name).ne.'character') then
    
    call h5dcreate_f(file_id,trim(MLSAuxData%name),type_id,dspace_id, &
          dset_id,h5error,cparms)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSET_CREATE // trim(MLSAuxData%name) )

    else

     call h5tcopy_f(type_id, s_type_id, h5error)
     if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_TYPE_COPY // trim(MLSAuxData%name))

     call h5tset_size_f(s_type_id, string_length, h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_TYPE_SET // trim(MLSAuxData%name))

     CALL h5dcreate_f(file_id, trim(MLSAuxData%name), s_type_id, &
       dspace_id, dset_id, h5error, cparms)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_DSET_CREATE // trim(MLSAuxData%name))

    endif 
!---------------- close all structures----------------------------------
    call h5sclose_f( dspace_id, h5error )
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSPACE_CLOSE // trim(MLSAuxData%name) )

    call h5pclose_f( cparms,    h5error )
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_PROPERTY_CLOSE // trim(MLSAuxData%name) )

    endif

       if (associated(MLSAuxData%RealField)) then
 
       rank = MLSAuxData%rank
 
       do i=1,rank
          dims_create(i) = size(MLSAuxData%RealField,i)
          dims(i) = dims_create(i)
          chunk_dims(i) = dims(i)
       end do
       endif

       if (associated(MLSAuxData%DpField)) then 

       rank = MLSAuxData%rank

       do i=1,rank
          dims_create(i) = size(MLSAuxData%DpField,i)
          dims(i) = dims_create(i)
          chunk_dims(i) = dims(i)
       end do
       endif

       if (associated(MLSAuxData%CharField)) then 

       rank = MLSAuxData%rank

       do i=1,rank
          dims_create(i) = 1
          dims(i) = 1
          chunk_dims(i) = 1
       end do

       endif
    
       if (associated(MLSAuxData%IntField)) then 

       rank = MLSAuxData%rank

       do i=1,rank
          dims_create(i) = size(MLSAuxData%IntField,i)
          dims(i) = dims_create(i)
          chunk_dims(i) = dims(i)
       end do
       endif

    do i = 1, rank
       start(i) = 0
    end do 

    if (PRESENT(index)) then 
       start(rank) = index-1
       dims_create(rank) = index
       dims(rank) = index
       chunk_dims(rank) = index
    endif 

    call h5dextend_f(dset_id, dims_create, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSET_EXTEND // trim(MLSAuxData%name) )

    call h5dget_space_f(dset_id, filespace, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSET_GETSPACE // trim(MLSAuxData%name) )

    call h5dget_type_f(dset_id, type_id, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSET_GET_TYPE // trim(MLSAuxData%name) )

    dims_create(rank)=1

    call h5screate_simple_f(rank, dims_create(1:rank), memspace, h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSPACE_CREATE // trim(MLSAuxData%name) )

    call h5sselect_hyperslab_f( & 
         filespace, H5S_SELECT_SET_F, start(1:rank), & 
         dims_create(1:rank), h5error)
    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSPACE_HYPERSLAB // trim(MLSAuxData%name) )

    test_type_name: select case (trim(MLSAuxData%type_name))
    case ('real')
     call h5dwrite_f(dset_id, type_id, MLSAuxData%RealField, dims(1:rank), &
          h5error,memspace,filespace)
    case ('double')
     call h5dwrite_f(dset_id, type_id, MLSAuxData%DpField, dims(1:rank), & 
          h5error,memspace,filespace)
    case ('integer')
     call h5dwrite_f(dset_id, type_id, MLSAuxData%IntField, dims(1:rank), &
          h5error,memspace,filespace)
    case('character')

    call h5tcopy_f(type_id, s_type_id, h5error)
     IF (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_TYPE_COPY // trim(MLSAuxData%name))

    call h5tset_size_f(s_type_id, string_length, h5error)
     IF (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_TYPE_SET // trim(MLSAuxData%name))

    call h5dwrite_f(dset_id, s_type_id, MLSAuxData%CharField, & 
         dims(1:rank), h5error,memspace,filespace)

    end select test_type_name

    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
       H5_ERROR_DSET_WRITE // trim(MLSAuxData%name) )

! Do we also write the attributes
   if ( myWrite_attributes ) then
      if ( associated(MLSAuxData%FrequencyCoordinates)) & 
     call Write_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
       & 'FrequencyCoordinates', error, MLSAuxData)
      if ( associated(MLSAuxData%HorizontalCoordinates)) & 
     call Write_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
       & 'HorizontalCoordinates', error, MLSAuxData)
      if ( associated(MLSAuxData%VerticalCoordinates)) & 
     call Write_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
       & 'VerticalCoordinates', error, MLSAuxData)
      if ( associated(MLSAuxData%Dimensions)) & 
     call Write_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
       & 'Dimensions', error, MLSAuxData)
   endif
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

  end subroutine Write_MLSAuxData

  subroutine CopyFromDataProducts( dataset, MLSAuxData )
    type( MLSAuxData_T ), intent(inout) :: MLSAuxData
    type( DataProducts_T ), intent(in) :: dataset
    integer :: i

    if (associated(dataset%Dimensions) ) then 
       do i = 1, size(dataset%Dimensions)
          MLSAuxData%Dimensions(i) = '                    '
          MLSAuxData%Dimensions(i) = trim(dataset%Dimensions(i))
       end do
    endif

    if (associated(dataset%FrequencyCoordinates) ) then 
       do i = 1, size(dataset%FrequencyCoordinates)
          MLSAuxData%FrequencyCoordinates(i) = dataset%FrequencyCoordinates(i)
       enddo 
    endif

    if (associated(dataset%VerticalCoordinates) ) then 
       do i = 1, size(dataset%VerticalCoordinates)
          MLSAuxData%VerticalCoordinates(i) = dataset%VerticalCoordinates(i)
       enddo 
    endif

    if (associated(dataset%HorizontalCoordinates) ) then 
       do i = 1, size(dataset%HorizontalCoordinates)
          MLSAuxData%HorizontalCoordinates(i) =dataset%HorizontalCoordinates(i)
       enddo 
    endif

  end subroutine CopyFromDataProducts

  subroutine CopyToDataProducts( MLSAuxData, dataset )
    type( MLSAuxData_T ), intent(in) :: MLSAuxData
    type( DataProducts_T ), intent(inout) :: dataset
    integer :: i

    if (associated(MLSAuxData%Dimensions) ) then
       do i = 1, size(MLSAuxData%Dimensions)
          dataset%Dimensions(i) = '                    '
          dataset%Dimensions(i) = trim(MLSAuxData%Dimensions(i))
       end do
    endif

    if (associated(MLSAuxData%FrequencyCoordinates) ) then 
       do i = 1, size(MLSAuxData%FrequencyCoordinates)
          dataset%FrequencyCoordinates(i) = MLSAuxData%FrequencyCoordinates(i)
       enddo 
    endif

    if (associated(MLSAuxData%VerticalCoordinates) ) then 
       do i = 1, size(MLSAuxData%VerticalCoordinates)
          dataset%VerticalCoordinates(i) = MLSAuxData%VerticalCoordinates(i)
       enddo 
    endif

    if (associated(MLSAuxData%HorizontalCoordinates) ) then 
       do i = 1, size(MLSAuxData%HorizontalCoordinates)
        dataset%HorizontalCoordinates(i) = MLSAuxData%HorizontalCoordinates(i)
       enddo 
    endif

  end subroutine CopyToDataProducts

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSAuxData

! $Log$
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
