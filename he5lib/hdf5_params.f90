module hdf5_params

! This module is a QnD way to provide some of the HDF5 parameters used by 
! HDF-EOS5 programs. These are not needed by the Fortran 90 interface to
! straight HDF5 -- it provides its own parameters via the module HDF5.

integer,public,parameter::H5F_ACC_RDWR=1, H5F_ACC_TRUNC=2, H5F_ACC_RDONLY=0
integer,public,parameter::H5S_UNLIMITED=-1

! is this right......
integer,public,parameter::H5T_NATIVE_DOUBLE=2, H5T_NATIVE_FLOAT=1
integer,public,parameter::H5T_NATIVE_INT=0,H5T_NATIVE_CHAR=3

! These kind types depend on (a) your f9x compiler and (b) how
! you built the hdf5 library. In particular, you may have disabled hsize_t
! because your C compiler doesn't do "long long" right in which case both 
! will be ordinary integers. (They are only in here for an experiment that 
! proved unnecessary -- The Fortran interface to HDF-EOS doesn't know abuot 
! HSIZE_T. )
integer,public,parameter::HSIZE_T=4,HID_T=3


end module hdf5_params
