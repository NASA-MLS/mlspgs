! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module h5_griddeddata

! This module provides functions to read and write griddeddata structures
! as normally read from a specific sort of ascii text file, (which I will
! refer to here as a "climatology file") to a HDF5
! file. This could be used for several things:
! (1) Plotting the contents of a climatology file in a package which 
!     supports HDF5 without having to write a climatology file parser for
!     that package
! (2) provide a drop-in replacement for l3ascii_read_field that supposes 
!     the working climatology for the retrieval is in HDF5 format instead
!     of ASCII.

use GriddedData, only: GriddedData_T
use hdf5, only: H5T_NATIVE_DOUBLE, hid_t, hsize_t, &
  & h5gcreate_f, h5screate_simple_f, h5dcreate_f, &
  & h5dwrite_f, h5dclose_f, h5sclose_f, h5gclose_f
implicit none

private

public::h5_write_griddeddata!,h5_get_uars_clim
integer, public, parameter :: MAX_RANK = 7

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

subroutine h5_write_griddeddata(loc_id,field)
! writes a GriddedData structure to a HDF5 file. 
  type(GriddedData_T), intent(in) :: field
  integer(kind=hid_t),intent(in)::loc_id
!--------------locals-------------------
  integer::error
  integer(kind=hid_t)::group_id,dspace_id,dset_id
  integer(kind=hsize_t),dimension(MAX_RANK)::dims
  integer(kind=hsize_t), dimension(MAX_RANK)::the_dims

  ! create group: all the stuff from this field goes in the group
  call h5gcreate_f(loc_id,field%quantityName, group_id, error)

  ! Write the actual data
  dims(1:6)=(/field%noHeights, field%noLats, field%noLons, field%noLsts, &
       field%noSzas, field%noDates/)
       the_dims = dims
  call h5screate_simple_f(6, dims, dspace_id, error)
  CALL h5dcreate_f(group_id, "field", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field%field, the_dims, error)
  CALL h5dclose_f(dset_id, error)
  CALL h5sclose_f(dspace_id, error)
  ! Write the axes:
  call h5screate_simple_f(1, (/dims(1)/), dspace_id, error)
  CALL h5dcreate_f(group_id, "heights", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field%heights, the_dims, error)
  CALL h5dclose_f(dset_id, error)
  CALL h5sclose_f(dspace_id, error)

  call h5screate_simple_f(1, (/dims(2)/), dspace_id, error)
  CALL h5dcreate_f(group_id, "lats", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  the_dims(:MAX_RANK-1) = dims(2:)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field%lats, the_dims, error)
  CALL h5dclose_f(dset_id, error)
  CALL h5sclose_f(dspace_id, error)

  call h5screate_simple_f(1, (/dims(4)/), dspace_id, error)
  CALL h5dcreate_f(group_id, "lsts", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  the_dims(:MAX_RANK-3) = dims(4:)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field%lsts, the_dims, error)
  CALL h5dclose_f(dset_id, error)
  CALL h5sclose_f(dspace_id, error)

  call h5screate_simple_f(1, (/dims(6)/), dspace_id, error)
  CALL h5dcreate_f(group_id, "dateStarts", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  the_dims(:MAX_RANK-5) = dims(6:)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field%dateStarts, the_dims, error)
  CALL h5dcreate_f(group_id, "dateEnds", H5T_NATIVE_DOUBLE, dspace_id, &
       dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field%dateEnds, the_dims, error)
  CALL h5dclose_f(dset_id, error)
  CALL h5sclose_f(dspace_id, error)

  call h5gclose_f(group_id,error)

end subroutine h5_write_griddeddata

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module h5_griddeddata


! $Log$
! Revision 2.10  2013/06/12 02:15:22  vsnyder
! Cruft removal
!
! Revision 2.9  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.8  2009/06/16 17:20:26  pwagner
! Restrict USEd items to only the ones needed
!
! Revision 2.7  2005/07/12 17:12:50  pwagner
! New hdf5 library will drop integer dimension interfaces
!
! Revision 2.6  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.5  2002/10/08 00:09:09  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/10/01 22:03:55  pwagner
! Fixed RCS Ident Block
!
! Revision 2.3  2002/10/01 20:04:16  bwknosp
! Added Log Info
!
