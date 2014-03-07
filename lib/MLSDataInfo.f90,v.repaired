! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

  module MLSDataInfo

  use MLSCOMMON, only: NAMELEN
  use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
  use HIGHOUTPUT, only: OUTPUTNAMEDVALUE
       
  implicit none

!     c o n t e n t s
!     - - - - - - - -
!     (data types and parameters)
! MLSDataInfo_T     Information from a HDF5 data file
!
!                 (subroutines and functions)
!
! Query_MLSData   Lists all dataset entries in a file.
!
  private

  public :: MLSDataInfo_T, Query_MLSData

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  type MLSDataInfo_T
    character(len=namelen), dimension(:), pointer :: name => null() ! Name of field in file
    integer :: number_of_entries
  end type MLSDataInfo_T

 contains ! ======================== MODULE PROCEDURES =======================

! -------------------------------------------------  Query_MLSData ----
  recursive subroutine Query_MLSData(loc_id, loc_name, dataset_info)
!
    use HDF5, only: hid_t, H5G_DATASET_F, H5G_LINK_F, &
    ! & H5G_UNKNOWN_F, H5G_GROUP_F, H5G_TYPE_F, &
      & h5gn_members_f,h5gget_obj_info_idx_f
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
    logical, parameter :: OBJINFOBROKEN = .true. ! starting with hdf5 1.8
    integer :: i, nmembers, h5error, count, type_id
    character(len=namelen) :: name_buffer, path_name, new_loc_name
    integer :: nsubmembers
    logical, parameter :: verbose = .false.
    ! call outputNamedValue( 'loc_id', loc_id )
    ! call outputNamedValue( 'loc_name', loc_name )

    call h5gn_members_f(loc_id,loc_name,nmembers,h5error)
    if ( verbose ) call outputNamedValue( 'loc_id', loc_id )
    if ( verbose ) call outputNamedValue( 'loc_name', loc_name )
    if ( verbose ) call outputNamedValue( 'nmembers', nmembers )
    if (h5error /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Group member access of ' // loc_name // ' failed.')

    if (h5error.eq.0) then 

    do i=0, nmembers-1

    call h5gget_obj_info_idx_f(loc_id,loc_name,i,name_buffer,type_id,h5error)
    if (h5error /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Retrieval of object information from ' // loc_name // ' failed.')

    if (h5error.eq.0) then 

          path_name = loc_name
          if (trim(path_name) /= "/") then                                         
            new_loc_name = trim(path_name) // "/" // trim(name_buffer) 
          else                                                                     
            new_loc_name = trim(path_name) // trim(name_buffer)        
          endif                                                                    
          ! call outputNamedValue( 'new loc_name', new_loc_name )

        ! call outputNamedValue( 'H5G_DATASET_F', H5G_DATASET_F )
        ! call outputNamedValue( 'H5G_LINK_F', H5G_LINK_F )
        ! call outputNamedValue( 'H5G_UNKNOWN_F', H5G_UNKNOWN_F )
        ! call outputNamedValue( 'H5G_GROUP_F', H5G_GROUP_F )
        ! call outputNamedValue( 'H5G_TYPE_F', H5G_TYPE_F )
        ! call outputNamedValue( 'type_id', type_id )
         ! if (type_id .EQ. 2) then  ! H5G_DATASET_F
        if ( OBJINFOBROKEN ) then
          ! h5gget_obj_info_idx_f doesn't return a reliable type_id
          ! so we'll call h5gn_members until we get an error
           call h5gn_members_f( loc_id, trim(new_loc_name), nsubmembers, h5error )
           if ( h5error == 0 .and. nsubmembers > 0 ) then
             path_name = new_loc_name
             call Query_MLSData( loc_id, trim(path_name), dataset_info)
           else
             count = dataset_info%number_of_entries + 1
             dataset_info%name(count) = new_loc_name

             dataset_info%number_of_entries = count
           endif
        elseif (type_id .EQ. H5G_DATASET_F) then  ! H5G_DATASET_F

          count = dataset_info%number_of_entries + 1
          dataset_info%name(count) = new_loc_name

          dataset_info%number_of_entries = count

        elseif (type_id .EQ. H5G_LINK_F) then  ! H5G_LINK_F

          count = dataset_info%number_of_entries + 1

          dataset_info%name(count) = new_loc_name

          dataset_info%number_of_entries = count

        else 

          path_name = new_loc_name
 
           call Query_MLSData(loc_id,trim(path_name),dataset_info)

        endif ! type_id     

    endif ! h5error

      end do ! nmembers

    endif ! h5error

  end subroutine Query_MLSData

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSDataInfo

! $Log$
! Revision 2.12  2014/03/07 19:15:51  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 2.11  2013/06/12 02:11:10  vsnyder
! Cruft removal
!
! Revision 2.10  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.9  2009/07/30 00:17:58  pwagner
! Worked around apparent bug in hdf5 1.8
!
! Revision 2.8  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.7  2007/01/12 00:29:28  pwagner
! Renamed routine outputNamedValue
!
! Revision 2.6  2006/06/29 20:36:55  pwagner
! Fixed bug besetting symbolic links
!
! Revision 2.5  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2005/04/29 21:53:29  pwagner
! Nullified name component of data type at declaration
!
! Revision 2.3  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2002/08/25 19:21:09  mjf
! Added a couple of continuation characters.
!
! Revision 2.1  2002/08/24 00:40:08  jdone
! Creating & adding to repository
!
