! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!
! JPD:: This program iterates through the groups and subgroups of an HDF5 file 
!       and lists the full path of each dataset in that HDF5 file.
!
     PROGRAM MLS_h5ls

     use HDF5, only: HID_T, H5F_ACC_RDONLY_F, h5fopen_f, h5fclose_f, h5fis_hdf5_f   
     use MLSDataInfo, only: MLSDataInfo_T, Query_MLSData

     IMPLICIT NONE

     type(MLSDataInfo_T) :: dataset_info

     CHARACTER(LEN=255) :: filename          ! filename
     INTEGER(HID_T) :: file_id               ! File identifier
     INTEGER     ::  i, count, status, error ! Counting indices & Error flags
     INTEGER, PARAMETER ::  max_nsds = 1000  ! Maximum number of datasets in file.
     LOGICAL     :: is_hdf5

!---------------------------- RCS Ident Info -------------------------------
  character(len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
     CALL h5open_f(error)

     write(*,21) 
 21  format("Enter the name of the HDF5 file. Datasets in the file will be listed shortly.")

     read(*,*) filename

     call h5fis_hdf5_f(trim(filename), is_hdf5, error)

     if (is_hdf5) then 

     CALL h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

     IF (error .eq. 0) then  

!
! The structure, dataset_info, is initialized below.
!
             allocate( dataset_info%name(1:max_nsds), stat=status)             
             dataset_info%name = ''
             dataset_info%number_of_entries = 0
!
! The real work is performed by MLSDataInfo::Query_MLSData.
!
      call Query_MLSData(file_id,'/', dataset_info)

     count = dataset_info%number_of_entries
!
! List the datasets.
!
        write(*,101) count, filename
 101    format(1x,"There are ",1i5," dataset entries in the HDF5 file, ",1a64)

     do i = 1, count
        write(*,102) dataset_info%name(i) 
 102    format(1x,1a64)
     end do
!
! Prevent memory leaks.
!
             if (associated(dataset_info%name)) then 
                 nullify(dataset_info%name)
                 deallocate(dataset_info%name, stat=status)
             endif

     ENDIF ! if h5fopen_f
                  
     CALL h5fclose_f(file_id, error)

     endif ! is_hdf5

     CALL h5close_f(error)
!-------------------------------------------------------------
     STOP
     END PROGRAM MLS_h5ls

! $Log$
