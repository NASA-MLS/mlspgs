! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MLSAuxData

  ! Reading and interacting with Level 1B data (HDF5)

  use HDF5, only: hid_t, hsize_t, H5T_NATIVE_CHARACTER, H5T_NATIVE_DOUBLE, &
       H5T_NATIVE_INTEGER, H5T_IEEE_F32LE, H5S_UNLIMITED_F, &
       H5T_STD_I32LE, H5T_NATIVE_REAL, H5T_IEEE_F64LE, &
       H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, &
       h5dopen_f, h5dclose_f, h5screate_simple_f, h5pcreate_f, h5dcreate_f, &
       h5pset_chunk_f, h5sclose_f, h5dget_space_f, h5pclose_f, & 
       h5sget_simple_extent_ndims_f, h5sget_simple_extent_dims_f, &
       h5dget_create_plist_f, h5pget_chunk_f, h5dget_type_f, & 
       h5sselect_hyperslab_f, h5dread_f, h5dwrite_f, h5dextend_f, &
       h5acreate_f, h5awrite_f, h5aread_f, h5aclose_f, h5tcopy_f, &
       h5tset_size_f, h5aopen_name_f, h5aget_type_f, h5aget_space_f
  use MLSCommon, only: r4, r8
  use MLSStrings, only: Lowercase

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
! Read_MLSAuxAttributes        Reads an attribute of a MLSAuxData quantity
!                              from a file.
! Read_MLSAuxData              Reads all info concerning a MLSAuxData quantity
!                              from a file.
! Write_MLSAuxAttributes       Writes an attribute of a MLSAuxData quantity
!                              to a file.
! Write_MLSAuxData             Writes all info concerning a MLSAuxData quantity
!                              to a file.
! Deallocate_MLSAuxData        Called when an MLSAuxData is finished.
!
  private

  public :: MLSAuxData_T, Create_MLSAuxData, Read_MLSAuxData, &
      & Write_MLSAuxData, Deallocate_MLSAuxData, &
      & Read_MLSAuxAttributes, Write_MLSAuxAttributes, &
      & NAME_LEN

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

    real(r8), dimension(:), pointer :: FrequencyCoordinates  => NULL()
    real(r8), dimension(:), pointer :: VerticalCoordinates   => NULL()
    real(r8), dimension(:), pointer :: HorizontalCoordinates => NULL()

    character, dimension(:,:,:), pointer :: CharField        => NULL()
    real(r8),  dimension(:,:,:), pointer :: DpField          => NULL()   
    real(r4),  dimension(:,:,:), pointer :: RealField        => NULL()   
    integer,   dimension(:,:,:), pointer :: IntField         => NULL()

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
    real, dimension(:), pointer :: attr_data
    integer(hsize_t), dimension(3) :: chunk_dims, dims, maxdims
    integer(hsize_t), dimension(7) :: adims
    integer(hsize_t), dimension(1) :: adims_create
   integer(hid_t) :: cparms,dspace_id,dset_id,type_id, &
         attr_id, atype_id, aspace_id
    integer :: i, rank_MAF, rank, arank, h5error
    logical, parameter ::  ATTRIBUTE = .FALSE.
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

!    A question:
! Is it better to do this here in a routine that creates the
! MLSAuxData, or to do it in Write_MLSAuxData ?
! For now we leave it here
! ( Might even break it out into its own separate routines:
!   Read_MLSAuxAttributes and Write_MLSAuxAttributes? )
!---------------- write the attributes -----------------------------

    do i = 2, 7
       adims(i) = 0
    end do

    if (associated(MLSAuxData%HorizontalCoordinates) .and. ATTRIBUTE) then 

       arank = size(shape(MLSAuxData%HorizontalCoordinates))
       atype_id = H5T_IEEE_F32LE
       adims(1) = size(MLSAuxData%HorizontalCoordinates)
       adims_create(1) = adims(1)
       allocate(attr_data(adims(1)))
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
       deallocate(attr_data)
    endif

    if (associated(MLSAuxData%VerticalCoordinates) .and. ATTRIBUTE) then 
       arank = size(shape(MLSAuxData%VerticalCoordinates))
       atype_id = H5T_IEEE_F32LE
       adims(1) = size(MLSAuxData%VerticalCoordinates)
       adims_create(1) = adims(1)
       allocate(attr_data(adims(1)))
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
       deallocate(attr_data)
    endif

    if (associated(MLSAuxData%FrequencyCoordinates) .and. ATTRIBUTE) then 
       arank = size(shape(MLSAuxData%FrequencyCoordinates))
       atype_id = H5T_IEEE_F32LE
       adims(1) = size(MLSAuxData%FrequencyCoordinates)
       adims_create(1) = adims(1)
       allocate(attr_data(adims(1)))
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
       deallocate(attr_data)
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
!-------------------------------------------------Read_MLSAuxAttributes ----
  subroutine Read_MLSAuxAttributes(dset_id, QuantityName, AttributeName, & 
       error, MLSAuxData, AttributeData)
  !
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
    integer, parameter :: OPTARGISMAD = 1
    integer, parameter :: OPTARGISARR = 2
    integer :: which_opt_arg
    real, dimension(:), pointer :: attr_data
    integer(hsize_t), dimension(7) :: adims
    integer(hsize_t), dimension(1) :: adims_create
    integer(hid_t) :: cparms,dspace_id,type_id, &
         attr_id, atype_id, aspace_id
    integer :: i, rank_MAF, rank, arank, h5error
  ! Executable
    error = 0
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
          adims(1) = size(MLSAuxData%FrequencyCoordinates)
        else
          error = 1
          return
        endif
      case ('HorizontalCoordinates')
        if ( associated ( MLSAuxData%HorizontalCoordinates ) ) then
          adims(1) = size(MLSAuxData%HorizontalCoordinates)
        else
          error = 1
          return
        endif
      case ('VerticalCoordinates')
        if ( associated ( MLSAuxData%VerticalCoordinates ) ) then
          adims(1) = size(MLSAuxData%VerticalCoordinates)
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
    call h5aget_space_f(attr_id, aspace_id, h5error)
    arank = 1
    allocate(attr_data(adims(1)))
    call h5aread_f(attr_id, atype_id, attr_data, adims, h5error)
  ! > >     print *, 'reading attribute ' // trim(AttributeName)
  ! > >     print *, 'adims(1) ', adims(1)
  ! > >     print *, 'attr_data ', attr_data
  ! > >     print *, 'h5error ', h5error
    select case (which_opt_arg)
    case (OPTARGISMAD)
      select case (trim(AttributeName))
      case ('FrequencyCoordinates')
        MLSAuxData%FrequencyCoordinates = attr_data(:adims(1))
      case ('HorizontalCoordinates')
        MLSAuxData%HorizontalCoordinates = attr_data(:adims(1))
      case ('VerticalCoordinates')
        MLSAuxData%VerticalCoordinates = attr_data(:adims(1))
      end select
    case (OPTARGISARR)
      AttributeData = attr_data(:adims(1))
    end select
    call h5aclose_f(attr_id, h5error)                                  
    call h5sclose_f(aspace_id, h5error)                                
    deallocate(attr_data)
    error = h5error                                              
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
    integer, parameter :: OPTARGISMAD = 1
    integer, parameter :: OPTARGISARR = 2
    integer :: which_opt_arg
    real, dimension(:), pointer :: attr_data
    integer(hsize_t), dimension(7) :: adims
    integer(hsize_t), dimension(1) :: adims_create
    integer(hid_t) :: cparms,dspace_id,type_id, &
         attr_id, atype_id, aspace_id
    integer :: i, rank_MAF, rank, arank, h5error
  ! Executable
    error = 0
    atype_id = H5T_IEEE_F32LE
  ! Check that exactly 1 of the optional args is supplied
  ! The following will produce 0 if neither, 3 if both; else 1 or 2
    which_opt_arg = 0
    if ( present(MLSAuxData) ) which_opt_arg = which_opt_arg + OPTARGISMAD
    if ( present(AttributeData) ) which_opt_arg = which_opt_arg + OPTARGISARR
    select case (which_opt_arg)
    case (OPTARGISMAD)
      select case (trim(AttributeName))
      case ('FrequencyCoordinates')
        arank = size(shape(MLSAuxData%FrequencyCoordinates))
        adims(1) = size(MLSAuxData%FrequencyCoordinates)
        allocate(attr_data(adims(1)))
        do i=1,adims(1)
           attr_data(i) = MLSAuxData%FrequencyCoordinates(i)
        end do
      case ('HorizontalCoordinates')
        arank = size(shape(MLSAuxData%HorizontalCoordinates))
        adims(1) = size(MLSAuxData%HorizontalCoordinates)
        allocate(attr_data(adims(1)))
        do i=1,adims(1)
           attr_data(i) = MLSAuxData%HorizontalCoordinates(i)
        end do
      case ('VerticalCoordinates')
        arank = size(shape(MLSAuxData%VerticalCoordinates))
        adims(1) = size(MLSAuxData%VerticalCoordinates)
        allocate(attr_data(adims(1)))
        do i=1,adims(1)
           attr_data(i) = MLSAuxData%VerticalCoordinates(i)
        end do
      case default
        error = 1
        return
      end select
    case (OPTARGISARR)
        arank = size(shape(AttributeData))
        adims(1) = size(AttributeData)
        allocate(attr_data(adims(1)))
        do i=1,adims(1)
           attr_data(i) = AttributeData(i)
        end do
    case default
      error = 1
      return
    end select
    adims_create(1) = adims(1)                                         
    call h5screate_simple_f(arank, adims_create, aspace_id, h5error)   
    call h5acreate_f(dset_id, trim(AttributeName), atype_id, &         
      & aspace_id, attr_id, h5error)                                   
    call h5awrite_f(attr_id, atype_id, attr_data, adims, h5error)      
  ! > >     print *, 'writing attribute ' // trim(AttributeName)
  ! > >     print *, 'adims(1) ', adims(1)
  ! > >     print *, 'attr_data ', attr_data
  ! > >     print *, 'h5error ', h5error
    call h5aclose_f(attr_id, h5error)                                  
    call h5sclose_f(aspace_id, h5error)                                
    deallocate(attr_data)
    error = h5error                                              
  end subroutine Write_MLSAuxAttributes
!-------------------------------------------------Read_MLSAuxData ----
  subroutine Read_MLSAuxData(file_id, QuantityName, QuantityType, & 
       MLSAuxData, error, FirstMAF, LastMAF, read_attributes)
!
! This subroutine reads an entry from the HDF5 file.
!
    type( MLSAuxData_T ), intent(inout) :: MLSAuxData
    character(len=*), intent(in)   :: QuantityName
    character(len=*), intent(in)   :: QuantityType ! if 'unknown' will derive
    integer(hid_t), intent(in)     :: file_id ! From HDF
    integer, intent(out)                          :: error   ! 0 unless trouble
    integer, intent(in), optional  :: FirstMAF ! First to read (default 0)
    integer, intent(in), optional  :: LastMAF  ! Last to read 
    logical, intent(in), optional  :: read_attributes  ! read freqCoords, etc 
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
    logical :: myRead_attributes
    character(len=16) :: myQuantityType

    error = 0
    myRead_attributes = .false.
    if ( present(read_attributes) ) myRead_attributes=read_attributes
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

    if ( lowercase(trim(QuantityType)) /= 'unknown' ) then
      myQuantityType = QuantityType
    else
  ! While case select seemed most natural here, these H5T_..
  ! unfortunately are not initialization expressions; therefore
  ! must use clunky old if elseif blocks
  !    select case (type_id)
  !    case (H5T_STD_I32LE, H5T_NATIVE_INTEGER)
      if ( type_id == H5T_STD_I32LE .or. type_id == H5T_STD_I32LE ) then
        myQuantityType = 'integer'
  !    case (H5T_NATIVE_CHARACTER)
      elseif ( type_id == H5T_NATIVE_CHARACTER ) then
        myQuantityType = 'character'
  !    case (H5T_NATIVE_REAL, H5T_IEEE_F32LE)
      elseif ( type_id == H5T_NATIVE_CHARACTER ) then
        myQuantityType = 'real'
  !    case (H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE)
      elseif ( type_id == H5T_NATIVE_CHARACTER ) then
        myQuantityType = 'double'
  !    case default
      else
        error = 1
        return
  !    end select
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

! based on type chose which h5dread to call.

    if ( .NOT.(PRESENT(FirstMAF) ) ) then 

    test_type: select case (trim(myQuantityType))
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
     
    test_type_dims: select case (trim(myQuantityType))
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

    if (associated(real_buffer)) then 
        nullify(real_buffer)
        deallocate(real_buffer, stat=status)
    endif

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

    if (associated(double_buffer)) then 
        nullify(double_buffer)
        deallocate(double_buffer, stat=status)
    endif

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

    if (associated(integer_buffer)) then 
        nullify(integer_buffer)
        deallocate(integer_buffer, stat=status)
    endif

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

    if (associated(character_buffer)) then 
        nullify(character_buffer)
        deallocate(character_buffer, stat=status)
    endif

    end select test_type_dims

    endif

!    if (h5error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, & 
!       H5_ERROR_DSET_READ // trim(QuantityName) )

! Do we also read the attributes
   if ( myRead_attributes ) then
     call Read_MLSAuxAttributes(dset_id, trim(QuantityName), & 
       & 'FrequencyCoordinates', h5error, MLSAuxData)
     call Read_MLSAuxAttributes(dset_id, trim(QuantityName), & 
       & 'HorizontalCoordinates', h5error, MLSAuxData)
     call Read_MLSAuxAttributes(dset_id, trim(QuantityName), & 
       & 'VerticalCoordinates', h5error, MLSAuxData)
   endif

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
  subroutine Write_MLSAuxData(file_id, MLSAuxData, error, &
    & startMAF, write_attributes)
!
! This routine assumes that the MLSAuxData dataset is created.
! This subroutine writes an entry to the HDF5 file.
!
    type( MLSAuxData_T ), intent(in) :: MLSAuxData
    integer(hid_t), intent(in)     :: file_id ! From HDF
    integer, intent(out)           :: error   ! 0 unless trouble
    integer, intent(in), optional  :: startMAF ! First to write (default 0)
    logical, intent(in), optional  :: write_attributes  ! write freqCoords, etc 
!-------------------------------------------------------------------------
! Internal variables.
!

    integer(hsize_t), dimension(7) :: dims
    integer(hsize_t), dimension(3) :: dims_create, start

    integer(hid_t) :: filespace, memspace, dset_id, dspace_id, type_id

    integer :: h5error, rank, i

!--------------------------------------------------------------------------

    logical :: myWrite_attributes

    error = 0
    myWrite_attributes = .false.
    if ( present(write_attributes) ) myWrite_attributes=write_attributes
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
! Do we also write the attributes
   if ( myWrite_attributes ) then
     call Write_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
       & 'FrequencyCoordinates', h5error, MLSAuxData)
     call Write_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
       & 'HorizontalCoordinates', h5error, MLSAuxData)
     call Write_MLSAuxAttributes(dset_id, trim(MLSAuxData%name), & 
       & 'VerticalCoordinates', h5error, MLSAuxData)
   endif

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
! Revision 2.4  2002/09/26 23:58:04  pwagner
! Moved attributes out of create function; standalone attribute io, too
!
! Revision 2.3  2002/09/09 05:43:39  jdone
! deallocate statements added for read buffers.
!
