! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

     PROGRAM MLS_ReadWriteData

     use HDF5, only: HID_T, H5F_ACC_RDONLY_F, H5F_ACC_TRUNC_F, &
       & h5fopen_f,h5fclose_f,h5fis_hdf5_f, h5fcreate_f   
     use MLSAuxData, only: MLSAuxData_T, Read_MLSAuxData, Write_MLSAuxData, & 
          Create_MLSAuxData, Deallocate_MLSAuxData

     IMPLICIT NONE

     type(MLSAuxData_T) :: write_data
     type(MLSAuxData_T) :: read_data

     CHARACTER(LEN=9), PARAMETER :: filename = 'myhdf5.h5' ! filename
     INTEGER(HID_T) :: file_id               ! File identifier
     INTEGER     ::  i, j, k, status, error ! Counting indices & Error flags
     INTEGER, PARAMETER ::  max    = 2 ! Maximum number of #frequencies.
     INTEGER, PARAMETER ::  maxMIF = 2 ! Maximum number of #MIFs.
     INTEGER, PARAMETER ::  maxMAF = 3 ! Maximum number of #MAFs.
     INTEGER, PARAMETER ::  maxVC  = 72 ! Maximum number of #vert. coords
     INTEGER, PARAMETER ::  maxHC  = 144 ! Maximum number of #horz. coords
     INTEGER, DIMENSION(maxMAF) :: counterMAF
     LOGICAL     :: is_hdf5
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

!-----------------------------------------------------------------------------
     CALL h5open_f(error)
       if ( error /= 0 ) then
         print *, 'Error in h5open_f'
         stop
       endif

     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
       if ( error /= 0 ) then
         print *, 'Error in creating file'
         stop
       endif
!
! The structure, dataset_info, is initialized below.
!
       allocate( write_data%RealField(max,maxMIF,maxMAF), stat=status)
       if ( status /= 0 ) then
         print *, 'Error in allocating write_data%RealField'
         stop
       endif
       write_data%name = 'myDataset'
       write_data%type_name = 'real'
!
! assign test data.
!       
       do k = 1, maxMAF
          do j = 1, maxMIF
             do i = 1, max
                write_data%RealField(i,j,k) = i*100 + j*10 + k
             end do
          end do
       end do

       allocate( write_data%FrequencyCoordinates(max), stat=status)
       if ( status /= 0 ) then
         print *, 'Error in allocating write_data%FrequencyCoordinates'
         stop
       endif
       do i=1, max
         write_data%FrequencyCoordinates(i) = i*44.1
       enddo
       allocate( write_data%VerticalCoordinates(maxVC), stat=status)
       if ( status /= 0 ) then
         print *, 'Error in allocating write_data%VerticalCoordinates'
         stop
       endif
       do i=1, maxVC
         write_data%VerticalCoordinates(i) = i*1000./maxVC
       enddo
       allocate( write_data%HorizontalCoordinates(maxHC), stat=status)
       if ( status /= 0 ) then
         print *, 'Error in allocating write_data%HorizontalCoordinates'
         stop
       endif
       do i=1, maxHC
         write_data%HorizontalCoordinates(i) = i*360./maxHC
       enddo

!
! Create structure within file.
!
       call Create_MLSAuxData(file_id, write_data)
!
! Write data to the file.
!
       call Write_MLSAuxData(file_id, write_data, write_attributes=.true.)
       print *,'write_data%RealField = ',write_data%RealField

!
! Prevent memory leaks.
!

       call Deallocate_MLSAuxData(write_data)                  

       CALL h5fclose_f(file_id, error)

       CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
!
! Read data
!
       allocate( read_data%RealField(max,maxMIF,maxMAF), stat=status)         
       allocate( read_data%FrequencyCoordinates(max), stat=status)
       allocate( read_data%VerticalCoordinates(maxVC), stat=status)
       allocate( read_data%HorizontalCoordinates(maxHC), stat=status)
       call Read_MLSAuxData(file_id,'myDataset','real',read_data, &
         & read_attributes=.true.)
       print *,'read_data%RealField = ',read_data%RealField
       print *,'read_data%FrequencyCoordinates = ',read_data%FrequencyCoordinates
       print *,'read_data%VerticalCoordinates = ',read_data%VerticalCoordinates
       print *,'read_data%HorizontalCoordinates = ',read_data%HorizontalCoordinates
!
! Prevent memory leaks.
!

             call Deallocate_MLSAuxData(read_data)                  

     CALL h5fclose_f(file_id, error)

     CALL h5close_f(error)
!-------------------------------------------------------------
     STOP
     END PROGRAM MLS_ReadWritedata

! $Log$
