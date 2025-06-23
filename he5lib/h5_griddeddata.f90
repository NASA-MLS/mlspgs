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
use l3ascii
use hdf5
implicit none

private

public::h5_write_griddeddata!,h5_get_uars_clim
integer, public, parameter :: MAX_RANK = 7

contains

subroutine h5_write_griddeddata(loc_id,field)
! writes a GriddedData structure to a HDF5 file. 
  type(GriddedData_T), intent(in) :: field
  integer(kind=hid_t),intent(in)::loc_id
!--------------locals-------------------
  integer::noHeights,error
  integer(kind=hid_t)::group_id,dspace_id,dset_id
  integer(kind=hsize_t),dimension(MAX_RANK)::dims
  integer, dimension(MAX_RANK)::the_dims
  noHeights=field%noHeights
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

end module h5_griddeddata
