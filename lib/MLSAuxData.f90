! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MLSAuxData

  ! Reading and interacting with Level 1B data (HDF5)

  use HDF5, only: hid_t, hsize_t, H5T_NATIVE_CHARACTER, H5T_NATIVE_DOUBLE, &
       H5T_NATIVE_INTEGER, H5T_IEEE_F32LE, H5S_UNLIMITED_F, &
       H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, &
       h5dopen_f, h5dclose_f, h5screate_simple_f, h5pcreate_f, h5dcreate_f, &
       h5pset_chunk_f, h5sclose_f, h5dget_space_f, h5pclose_f, & 
       h5sget_simple_extent_ndims_f, h5sget_simple_extent_dims_f, &
       h5dget_create_plist_f, h5pget_chunk_f, h5dget_type_f, & 
       h5sselect_hyperslab_f, h5dread_f, h5dwrite_f, h5dextend_f, &
       h5acreate_f, h5awrite_f, h5aread_f, h5aclose_f, h5tcopy_f, &
       h5tset_size_f, h5aopen_name_f, h5aget_type_f, h5aget_space_f

  implicit NONE

!     c o n t e n t s
!     - - - - - - - -

!     (data types and parameters)
! MLSAuxData_T                   Quantities from an L1B/L2 data file
! NAME_LEN                       Max length of l1b/l2 sds array name

!                 (subroutines and functions)
!
! Allocate_MLSAuxData          Allocates a MLSAuxData data structure. 
! Create_MLSAuxData            Creates MLSAuxData in a file. 
! Read_MLSAuxData              Reads all info concerning a MLSAuxData quantity
!                              from a file.
! Write_MLSAuxData             Writes all info concerning a MLSAuxData quantity
!                              to a file.
! Deallocate_MLSAuxData        Called when an MLSAuxData is finished.
!
  private

  public :: MLSAuxData_T, Create_MLSAuxData, Read_MLSAuxData, &
      Write_MLSAuxData, Deallocate_MLSAuxData,&
      NAME_LEN

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Parameters
  integer, parameter :: NAME_LEN = 64  ! Max len of SDS array name

  type MLSAuxData_T
    character (len=name_len) :: name       ! Name of field in file
    character (len=name_len) :: type_name  ! Name of type

    integer :: FrequencyDimension          ! Enumerated type
    integer :: VerticalDimension           ! Enumerated type  
    integer :: HorizontalDimension         ! Enumerated type

    integer, dimension(:), pointer :: FrequencyCoordinates
    integer, dimension(:), pointer :: VerticalCoordinates
    integer, dimension(:), pointer :: HorizontalCoordinates

    character, dimension(:,:,:), pointer :: CharField => NULL()
    real,  dimension(:,:,:), pointer :: DpField => NULL()
    real,  dimension(:,:,:), pointer :: RealField => NULL()
    integer,   dimension(:,:,:), pointer :: IntField => NULL()

    ! all the above dimensioned (noAuxInds,maxMIFs,noMAFs)
  end type MLSAuxData_T

  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_OPEN   = 'HDF5 Error Opening Dataset '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_CREATE = 'HDF5 Error Creating Dataset '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_DSET_WRITE  = 'HDF5 Error Writing Dataset '
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
    H5_ERROR_ATT_READ  = 'HDF5 Error Reading Attribute '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_ATT_CREATE = 'HDF5 Error Creating Attribute '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_ATT_CLOSE  = 'HDF5 Error Closing Attribute '

  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_TYPE_COPY = 'HDF5 Error Copying Datatype '
  CHARACTER(len=*), PUBLIC, PARAMETER :: & 
    H5_ERROR_TYPE_SET  = 'HDF5 Error Setting Datatype '

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
  
contains ! ============================ MODULE PROCEDURES ====================

  !-------------------------------------------  Deallocate_MLSAuxData  -----
 subroutine Deallocate_MLSAuxData( MLSAuxData )
    ! This should be called when deallocating a MLSAuxData structure.
    type( MLSAuxData_T ), intent(inout) :: MLSAuxData
    integer :: status

    if (associated(MLSAuxData%RealField)) then 
        nullify(MLSAuxData%RealField)
        deallocate(MLSAuxData%RealField, stat=status)
! check on status.
    endif

    if (associated(MLSAuxData%IntField)) then 
        nullify(MLSAuxData%IntField)
        deallocate(MLSAuxData%IntField, stat=status)
! check on status.
    endif

    if (associated(MLSAuxData%CharField)) then 
        nullify(MLSAuxData%CharField)
        deallocate(MLSAuxData%CharField, stat=status)
! check on status.
    endif

    if (associated(MLSAuxData%DpField)) then 
        nullify(MLSAuxData%DpField)
        deallocate(MLSAuxData%DpField, stat=status)
! check on status.
    endif

    if (associated(MLSAuxData%FrequencyCoordinates)) then 
        nullify(MLSAuxData%FrequencyCoordinates)
        deallocate(MLSAuxData%FrequencyCoordinates, stat=status)
! check on status.
    endif

    if (associated(MLSAuxData%VerticalCoordinates)) then 
        nullify(MLSAuxData%VerticalCoordinates)
        deallocate(MLSAuxData%VerticalCoordinates, stat=status)
! check on status.
    endif

    if (associated(MLSAuxData%HorizontalCoordinates)) then 
        nullify(MLSAuxData%HorizontalCoordinates)
        deallocate(MLSAuxData%HorizontalCoordinates, stat=status)
! check on status.
    endif

 end subroutine Deallocate_MLSAuxData
! -------------------------------------------------  Create_MLSAuxData ----
  subroutine Create_MLSAuxData(file_id, MLSAuxData)
!
! This subroutine creates an entry in the HDF5 file.
!----------------------------------------------------------------------
! External variables
!
    type( MLSAuxData_T ), intent(in) :: MLSAuxData
    integer(hid_t), intent(in)     :: file_id ! From HDF
!-----------------------------------------------------------------------
! Internal variables
!
    character(len=name_len) :: aname
    real, dimension(3) :: attr_data
    integer(hsize_t), dimension(3) :: chunk_dims, dims, maxdims
    integer(hsize_t), dimension(7) :: adims
    integer(hsize_t), dimension(1) :: adims_create
   integer(hid_t) :: cparms,dspace_id,dset_id,type_id, &
         attr_id, atype_id, aspace_id
    integer :: i, rank_MAF, rank, arank, h5error
!-----------------------------------------------------------------------

    test_type: select case (trim(MLSAuxData%type_name))
    case ('real')
       if (associated(MLSAuxData%RealField)) then 
       rank = size(shape(MLSAuxData%RealField))
       type_id = H5T_IEEE_F32LE
       do i=1,rank
          dims(i) = size(MLSAuxData%RealField, i)
          chunk_dims(i) = dims(i)
          maxdims(i) = H5S_UNLIMITED_F
       end do
       end if
    case ('double')
       if (associated(MLSAuxData%DpField)) then 
       rank = size(shape(MLSAuxData%DpField))
       type_id = H5T_NATIVE_DOUBLE
       do i=1,rank
          dims(i) = size(MLSAuxData%DpField, i)
          chunk_dims(i) = dims(i)
          maxdims(i) = H5S_UNLIMITED_F
       end do
       endif
    case ('integer')
       if (associated(MLSAuxData%IntField)) then 
       rank = size(shape(MLSAuxData%IntField))
       type_id = H5T_NATIVE_INTEGER
       do i=1,rank
          dims(i) = size(MLSAuxData%IntField, i)
          chunk_dims(i) = dims(i)
          maxdims(i) = H5S_UNLIMITED_F
       end do
       endif
    case ('character')
       if (associated(MLSAuxData%CharField)) then 
       rank = size(shape(MLSAuxData%CharField))
       type_id = H5T_NATIVE_CHARACTER
       do i=1, rank
          dims(i) = size(MLSAuxData%CharField, i)
          chunk_dims(i) = dims(i)
          maxdims(i) = H5S_UNLIMITED_F
       end do
       endif
    end select test_type
!--------------------------------------------------------------------
    call h5screate_simple_f(rank, dims(1:rank),dspace_id,h5error,& 
         maxdims(1:rank))
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_CREATE // QuantityName )

    call h5pcreate_f(H5P_DATASET_CREATE_F,cparms,h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_PROPERTY_CREATE // QuantityName )

    call h5pset_chunk_f(cparms,rank,chunk_dims(1:rank),h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_PROPERTY_CHUNK_SET // QuantityName )

    call h5dcreate_f(file_id,trim(MLSAuxData%name),type_id,dspace_id, &
          dset_id,h5error,cparms)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_CREATE // QuantityName ) 

!---------------- write the attributes -----------------------------

    do i = 2, 7
       adims(i) = 0
    end do

    if (associated(MLSAuxData%HorizontalCoordinates)) then 

       arank = size(shape(MLSAuxData%HorizontalCoordinates))
       atype_id = H5T_IEEE_F32LE
       adims(1) = size(MLSAuxData%HorizontalCoordinates)
       adims_create(1) = adims(1)
       do i=1, adims(1)
          attr_data(i) = MLSAuxData%HorizontalCoordinates(i)
       end do

    call h5screate_simple_f(arank, adims_create, aspace_id, h5error)
    aname = "HorizontalCoordinates"  
    call h5acreate_f(dset_id, trim(aname), atype_id, aspace_id, attr_id, &
       h5error)
    call h5awrite_f(attr_id, atype_id, attr_data, adims, h5error)
    call h5aclose_f(attr_id, h5error)
    call h5sclose_f(aspace_id, h5error) 
    endif

    if (associated(MLSAuxData%VerticalCoordinates)) then 
       arank = size(shape(MLSAuxData%VerticalCoordinates))
       atype_id = H5T_IEEE_F32LE
       adims(1) = size(MLSAuxData%VerticalCoordinates)
       adims_create(1) = adims(1)
       do i=1,adims(1)
          attr_data(i) = MLSAuxData%VerticalCoordinates(i)
       end do

    call h5screate_simple_f(arank, adims_create, aspace_id, h5error)
    aname = "VerticalCoordinates"   
    call h5acreate_f(dset_id, trim(aname), atype_id, aspace_id, attr_id, &
       h5error)
    call h5awrite_f(attr_id, atype_id, attr_data, adims, h5error)
    call h5aclose_f(attr_id, h5error)
    call h5sclose_f(aspace_id, h5error) 
    endif

    if (associated(MLSAuxData%FrequencyCoordinates)) then 
       arank = size(shape(MLSAuxData%FrequencyCoordinates))
       atype_id = H5T_IEEE_F32LE
       adims(1) = size(MLSAuxData%FrequencyCoordinates)
       adims_create(1) = adims(1)
       do i=1,adims(1)
          attr_data(i) = MLSAuxData%FrequencyCoordinates(i)
       end do

    call h5screate_simple_f(arank, adims_create, aspace_id, h5error)
    aname = "FrequencyCoordinates"
    call h5acreate_f(dset_id, trim(aname), atype_id, aspace_id, attr_id, &
       h5error)
    call h5awrite_f(attr_id, atype_id, attr_data, adims, h5error)
    call h5aclose_f(attr_id, h5error)
    call h5sclose_f(aspace_id, h5error) 
    endif
!---------------- close all structures----------------------------------
!
    call h5sclose_f( dspace_id, h5error )
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_CLOSE // QuantityName )

    call h5dclose_f( dset_id,   h5error )
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_CLOSE // QuantityName )

    call h5pclose_f( cparms,    h5error )
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_PROPERTY_CLOSE // QuantityName )

  end subroutine Create_MLSAuxData
!-------------------------------------------------Read_MLSAuxData ----
  subroutine Read_MLSAuxData(file_id, QuantityName, QuantityType, & 
       MLSAuxData, FirstMAF, LastMAF)
!
! This subroutine reads an entry from the HDF5 file.
!
    type( MLSAuxData_T ), intent(inout) :: MLSAuxData
    character(len=*), intent(in)   :: QuantityName, QuantityType
    integer(hid_t), intent(in)     :: file_id ! From HDF
    integer, intent(in), optional  :: FirstMAF ! First to read (default 0)
    integer, intent(in), optional  :: LastMAF  ! Last to read 
!
!-------------------------------------------------------------------------
!
    real, dimension(:,:,:), pointer :: real_buffer, double_buffer
    integer, dimension(:,:,:), pointer :: integer_buffer
    character, dimension(:,:,:), pointer :: character_buffer

    integer(hsize_t), dimension(7) :: dims, adims
    integer(hsize_t), dimension(3) :: chunk_dims, dims_create, maxdims, start
    integer(hid_t) :: cparms, dset_id, dspace_id, type_id, memspace, & 
         attr_id, atype_id, aspace_id
    integer        :: rank, h5error, i, j, k, status

!-- assign name and type for dataset

    MLSAuxData%name = trim(QuantityName)
    MLSAuxData%type_name = trim(QuantityType)

!---
    call h5dopen_f(file_id,trim(QuantityName),dset_id,h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_OPEN // trim(QuantityName) )

    call h5dget_space_f(dset_id,dspace_id,h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_CREATE // trim(QuantityName) )

    call h5sget_simple_extent_ndims_f(dspace_id,rank,h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_GET_NDIMS // trim(QuantityName) )

    call h5sget_simple_extent_dims_f(dspace_id,dims_create(1:rank),&
         maxdims(1:rank),h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_GET_DIMS // trim(QuantityName) )
!
    call h5dget_create_plist_f(dset_id,cparms,h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_CREATE // trim(QuantityName) )

    call h5pget_chunk_f(cparms,rank,chunk_dims(1:rank),h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_PROPERTY_GET // trim(QuantityName) )

    call h5screate_simple_f(rank, dims_create(1:rank), memspace, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_CREATE // trim(QuantityName) )

    call h5dget_type_f(dset_id, type_id, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_TYPE // trim(QuantityName) )

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

! based on type chose which h5dread to call.

    if ( .NOT.(PRESENT(FirstMAF) ) ) then 

    test_type: select case (trim(QuantityType))
    case ('real')
    call h5dread_f(dset_id,type_id,MLSAuxData%RealField,dims,& 
         h5error, memspace, dspace_id)
    case ('double')    
    call h5dread_f(dset_id,type_id,MLSAuxData%DpField,dims,&
         h5error,memspace, dspace_id)
    case ('integer')  
    call h5dread_f(dset_id,type_id,MLSAuxData%IntField,dims,& 
         h5error,memspace, dspace_id)
    case ('character') 
    call h5dread_f(dset_id,type_id,MLSAuxData%CharField,dims,&
         h5error,memspace,dspace_id)
    end select test_type
    else
     
    test_type_dims: select case (trim(QuantityType))
    case ('real')
    allocate( real_buffer(dims(1),dims(2),dims(3)),stat=status)
    call h5dread_f(dset_id,type_id,& 
         real_buffer,dims,h5error, memspace, dspace_id)
         do k = FirstMAF, LastMAF
            do j = 1, dims(2)
               do i = 1, dims(1)
                  MLSAuxData%RealField(i,j,k) = real_buffer(i,j,k)
               end do 
            end do
         end do
    case ('double')    

    allocate( double_buffer(dims(1),dims(2),dims(3)),stat=status)
    call h5dread_f(dset_id,type_id,& 
         double_buffer,dims,h5error, memspace, dspace_id)
         do k = FirstMAF, LastMAF
            do j = 1, dims(2)
               do i = 1, dims(1)
                  MLSAuxData%DpField(i,j,k) = double_buffer(i,j,k)
               end do 
            end do
         end do    
    case ('integer')  
    allocate( integer_buffer(dims(1),dims(2),dims(3)),stat=status)
    call h5dread_f(dset_id,type_id,& 
         integer_buffer,dims,h5error, memspace, dspace_id)
         do k = FirstMAF, LastMAF
            do j = 1, dims(2)
               do i = 1, dims(1)
                  MLSAuxData%IntField(i,j,k) = integer_buffer(i,j,k)
               end do 
            end do
         end do
    case ('character') 
    allocate( character_buffer(dims(1),dims(2),dims(3)),stat=status)
    call h5dread_f(dset_id,type_id,& 
         character_buffer,dims,h5error, memspace, dspace_id)
         do k = FirstMAF, LastMAF
            do j = 1, dims(2)
               do i = 1, dims(1)
                  MLSAuxData%CharField(i,j,k) = character_buffer(i,j,k)
               end do 
            end do
         end do
    end select test_type_dims

    endif

!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_READ // trim(QuantityName) )

!
! Close all identifiers.
!
    call h5sclose_f(memspace, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_CLOSE // trim(QuantityName) )
 
    call h5sclose_f(dspace_id, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_CLOSE // trim(QuantityName) )

    call h5dclose_f(dset_id, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_CLOSE // trim(QuantityName) )

    call h5pclose_f(cparms, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_PROPERTY_CLOSE // trim(QuantityName) )

  end subroutine Read_MLSAuxData
  ! -------------------------------------------------  Write_MLSAuxData ----
  subroutine Write_MLSAuxData(file_id, MLSAuxData, startMAF)
!
! This routine assumes that the MLSAuxData dataset is created.
! This subroutine writes an entry to the HDF5 file.
!
    type( MLSAuxData_T ), intent(in) :: MLSAuxData
    integer(hid_t), intent(in)     :: file_id ! From HDF
    integer, intent(in), optional  :: startMAF ! First to write (default 0)
!-------------------------------------------------------------------------
! Internal variables.
!

    integer(hsize_t), dimension(7) :: dims
    integer(hsize_t), dimension(3) :: dims_create, start

    integer(hid_t) :: filespace, memspace, dset_id, dspace_id, type_id

    integer :: h5error, rank, i

!--------------------------------------------------------------------------

    call h5dopen_f(file_id, trim(MLSAuxData%name), dset_id, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_OPEN // trim(QuantityName) )

       if (associated(MLSAuxData%RealField)) then 
       rank = size(shape(MLSAuxData%RealField))
       do i=1,rank
          dims_create(i) = size(MLSAuxData%RealField,i)
          dims(i) = dims_create(i)
       end do
       endif

       if (associated(MLSAuxData%DpField)) then 
       rank = size(shape(MLSAuxData%DpField))
       do i=1,rank
          dims_create(i) = size(MLSAuxData%DpField,i)
          dims(i) = dims_create(i)
       end do
       endif

       if (associated(MLSAuxData%CharField)) then 
       rank = size(shape(MLSAuxData%CharField))
       do i=1,rank
          dims_create(i) = size(MLSAuxData%CharField,i)
          dims(i) = dims_create(i)
       end do
       endif
    
       if (associated(MLSAuxData%IntField)) then 
       rank = size(shape(MLSAuxData%IntField))
       do i=1,rank
          dims_create(i) = size(MLSAuxData%IntField,i)
          dims(i) = dims_create(i)
       end do
       endif

    do i = 1, rank
       start(i) = 0
    end do 
    if (PRESENT(startMAF)) start(rank) = startMAF-1

    call h5dextend_f(dset_id, dims_create(1:rank), h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_EXTEND // trim(QuantityName) )

    call h5dget_space_f(dset_id, filespace, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_GET // trim(QuantityName) )

    call h5dget_type_f(dset_id, type_id, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_GET_TYPE // trim(QuantityName) )

    call h5screate_simple_f(rank, dims_create(1:rank), memspace, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_CREATE // trim(QuantityName) )

    call h5sselect_hyperslab_f( & 
         filespace, H5S_SELECT_SET_F, start(1:rank), & 
         dims_create(1:rank), h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_HYPERSLAB // trim(QuantityName) )

    if ( .NOT.( PRESENT(startMAF) ) ) then 

    test_type: select case (trim(MLSAuxData%type_name))
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
     call h5dwrite_f(dset_id, type_id, MLSAuxData%CharField, dims(1:rank), &
          h5error,memspace,filespace)
    end select test_type

    else

    print *,'Does not yet handle optional MAF I/O.'

    test_type_dims: select case (trim(MLSAuxData%type_name))
    case ('real')
    call h5dwrite_f(dset_id,type_id,& 
         MLSAuxData%RealField(1:dims(1),1:dims(2),startMAF:dims(3)),& 
         dims,h5error,memspace,filespace)
    case ('double')    
    call h5dwrite_f(dset_id,type_id,& 
         MLSAuxData%DpField(1:dims(1),1:dims(2),startMAF:dims(3)),& 
         dims,h5error,memspace,filespace)
    case ('integer')  
    call h5dwrite_f(dset_id,type_id,&
         MLSAuxData%IntField(1:dims(1),1:dims(2),startMAF:dims(3)),&         
         dims,h5error,memspace,filespace)
    case ('character') 
    call h5dwrite_f(dset_id,type_id,&
         MLSAuxData%CharField(1:dims(1),1:dims(2),startMAF:dims(3)),&     
         dims,h5error,memspace,filespace)
    end select test_type_dims    
    endif


!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_WRITE // trim(QuantityName) )
!
! Close all identifiers.
!
    call h5sclose_f(filespace, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_CLOSE // trim(QuantityName) )
 
    call h5sclose_f(memspace, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSPACE_CLOSE // trim(QuantityName) )

    call h5dclose_f(dset_id, h5error)
!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_CLOSE // trim(QuantityName) )

  end subroutine Write_MLSAuxData

end module MLSAuxData

! $Log$
! Revision 2.2  2002/09/09 05:26:15  jdone
! Now reads portion of dataset
!
