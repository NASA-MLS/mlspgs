! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

  module MLSDataInfo

  use HDF5, only: hid_t, h5gn_members_f,h5gget_obj_info_idx_f
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
       
  implicit none

!     c o n t e n t s
!     - - - - - - - -
!     (data types and parameters)
! MLSDataInfo_T     Information from a HDF5 data file
! NAME_LEN          Max length of an sds array name
!
!                 (subroutines and functions)
!
! Query_MLSData   Lists all dataset entries in a file.
!
  private

  public :: MLSDataInfo_T, Query_MLSData, NAME_LEN

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Parameters
  integer, parameter :: name_len = 64  ! Max len of SDS array name

  type MLSDataInfo_T
    character(len=name_len), dimension(:), pointer :: name ! Name of field in file
    integer :: number_of_entries
  end type MLSDataInfo_T

 contains ! ======================== MODULE PROCEDURES =======================

! -------------------------------------------------  Query_MLSData ----
  recursive subroutine Query_MLSData(loc_id, loc_name, dataset_info)
!
! This subroutine lists entries in the HDF5 file.
!
! define external variables
!
    type(MLSDataInfo_T), intent(inout) :: dataset_info 
    integer (hid_t), intent(in) :: loc_id
    character(len=*), intent(in) :: loc_name
!
! define internal variables.
!
    integer :: i, nmembers, h5error, count, type_id
    character(len=name_len) :: name_buffer, path_name
    logical :: dataset_found

    call h5gn_members_f(loc_id,loc_name,nmembers,h5error)
    if (h5error /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, 
        & 'Group member access of ' // loc_name // ' failed.')

    if (h5error.eq.0) then 

    do i=0, nmembers-1

    call h5gget_obj_info_idx_f(loc_id,loc_name,i,name_buffer,type_id,h5error)
    if (h5error /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, 
        & 'Retrieval of object information from ' // loc_name // ' failed.')

    if (h5error.eq.0) then 

          path_name = loc_name

         if (type_id .EQ. 2) then  ! H5G_DATASET_F

          dataset_found = .TRUE.
          count = dataset_info%number_of_entries + 1

       if (trim(path_name) /= "/") then 
         dataset_info%name(count) = trim(path_name) // "/" // trim(name_buffer)
       else 
         dataset_info%name(count) = trim(path_name) // trim(name_buffer)
       endif

          dataset_info%number_of_entries = count

         else 

           dataset_found = .FALSE.
 
           if (trim(path_name) /= "/") then 
            path_name = trim(path_name) // "/" // trim(name_buffer)
           else 
            path_name = trim(path_name) // trim(name_buffer)
           endif
 
           call Query_MLSData(loc_id,trim(path_name),dataset_info)

         endif ! H5G_DATASET_F     

       endif ! h5error

      end do ! nmembers

    endif ! h5error

  end subroutine Query_MLSData
end module MLSDataInfo

! $Log$
! Revision 2.1  2002/08/24 00:40:08  jdone
! Creating & adding to repository
!
