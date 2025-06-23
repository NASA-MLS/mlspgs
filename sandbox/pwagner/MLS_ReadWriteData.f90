! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

     PROGRAM MLS_ReadWriteData

     use MLSAuxData, only: MLSAuxData_T, Read_MLSAuxData, Write_MLSAuxData, & 
       & Create_MLSAuxData, Deallocate_MLSAuxData
     use MLSHDF5, only: IsHDF5AttributeInFile, IsHDF5DSInFile, &
       & IsHDF5AttributePresent, IsHDF5DSPresent, &
       & GetHDF5DSRank, GetHDF5DSDims, GetHDF5DSQType, LoadFromHDF5DS
     use L1BData, only: dump, l1bdata_t, findl1bdata, ReadL1BData, &
       & AssembleL1BQtyName
     use MLSCommon, only: FileNameLen, i4
     use MLSFiles, only: mls_io_gen_closeF, mls_io_gen_openF, HDFVERSION_5
     use MLSStrings, only: CompressString
     use Hdf, only: DFACC_RDONLY
     use HDF5, only: HID_T, H5F_ACC_RDONLY_F, H5F_ACC_TRUNC_F, &
       & h5fopen_f,h5fclose_f,h5fis_hdf5_f, h5fcreate_f, hSize_t, &
       & h5dOpen_f, h5dget_type_f, h5dClose_f, H5TEqual_f, &
       & H5S_SCALAR_F, H5T_NATIVE_INTEGER, H5T_NATIVE_CHARACTER, &
       & H5T_NATIVE_DOUBLE, H5T_NATIVE_REAL, HID_T, HSIZE_T, &
       & H5T_STD_I32LE, H5T_IEEE_F32LE, H5T_IEEE_F64LE

     IMPLICIT NONE

     type(MLSAuxData_T) :: write_data
     type(MLSAuxData_T) :: read_data
     type(l1bdata_t) :: l1b_data

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
     integer     :: rank
     integer(kind=hSize_t), dimension(:), pointer :: dims
     character(len=FileNameLen) :: the_file_name
    integer :: type_id                 ! typeID for DS
    integer :: SETID                    ! ID for DS
     character(len=*), parameter :: l1b_path = &
       & '/users/jdone/mlspgs/l1/outputs-15-September-2002/'
     character(len=*), parameter :: l1boa_file = &
       & 'MLS-Aura_L1BOA_V0-5-C01_2002-032.dat'
     character(len=*), parameter :: l1bradd_file = &
       & 'MLS-Aura_L1BRADD_V0-5-C01_2002-032.dat'
     character(len=*), parameter :: l1bradf_file = &
       & 'MLS-Aura_L1BRADF_V0-5-C01_2002-032.dat'
     character(len=16) :: QType
     character(len=32) :: QName
     integer :: l1boa_id, l1brad(2), noMAFs
    integer(i4)  :: ErrType
    integer(i4)  :: record_length
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
       call Write_MLSAuxData(file_id, write_data, error, &
         & write_attributes=.true.)
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
       call Read_MLSAuxData(file_id, 'myDataset', 'real', read_data, &
         & error, read_attributes=.true.)
       print *,'read_data%RealField = ',read_data%RealField
       print *,'read_data%FrequencyCoordinates = ',read_data%FrequencyCoordinates
       print *,'read_data%VerticalCoordinates = ',read_data%VerticalCoordinates
       print *,'read_data%HorizontalCoordinates = ',read_data%HorizontalCoordinates
!
! Prevent memory leaks.
!

             call Deallocate_MLSAuxData(read_data)                  

     CALL h5fclose_f(file_id, error)

!-------------------------------------------------------------

!-------------- Now test parts of MLSHDF5 -----------------
     is_hdf5 = IsHDF5DSInFile(filename, 'myDataset')
     print *, 'Dataset in file? ', is_hdf5
     is_hdf5 = IsHDF5DSInFile(filename, 'notmyDataset')
     print *, 'Dataset in file (wrong DS name)? ', is_hdf5
     is_hdf5 = IsHDF5AttributeInFile(filename, 'myDataset', 'FrequencyCoordinates')
     print *, 'Frequency coordinates attribute in file? ', is_hdf5
     is_hdf5 = IsHDF5AttributeInFile(filename, 'myDataset', 'HorizontalCoordinates')
     print *, 'Horizontal coordinates attribute in file? ', is_hdf5
     is_hdf5 = IsHDF5AttributeInFile(filename, 'myDataset', 'VerticalCoordinates')
     print *, 'Vertical coordinates attribute in file? ', is_hdf5
     is_hdf5 = IsHDF5AttributeInFile(filename, 'notmyDataset', 'VerticalCoordinates')
     print *, 'Vertical coordinates attribute in file (wrong DS name)? ', is_hdf5
     is_hdf5 = IsHDF5AttributeInFile(filename, 'myDataset', 'DiagonalCoordinates')
     print *, 'Diagonal coordinates attribute in file? ', is_hdf5
       CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
     is_hdf5 = IsHDF5DSPresent(file_id, 'myDataset')
     print *, 'Dataset present? ', is_hdf5
     is_hdf5 = IsHDF5AttributePresent(file_id, 'myDataset', 'FrequencyCoordinates')
     print *, 'Frequency coordinates attribute present? ', is_hdf5
     is_hdf5 = IsHDF5AttributePresent(file_id, 'myDataset', 'HorizontalCoordinates')
     print *, 'Horizontal coordinates attribute present? ', is_hdf5
     is_hdf5 = IsHDF5AttributePresent(file_id, 'myDataset', 'VerticalCoordinates')
     print *, 'Vertical coordinates attribute present? ', is_hdf5
     
     call GetHDF5DSRank(file_id, 'myDataset', rank)
     print *, 'DS rank ', rank
     allocate(dims(rank))
     call GetHDF5DSDims(file_id, 'myDataset', dims)
     print *, 'DS dimensions ', dims
     call GetHDF5DSQType(file_id, 'myDataset', QType)
     print *, 'DS quantity type ', trim(QType)
    call h5dOpen_f ( File_ID, 'myDataset', setID, status )
    print *, file_id, setid 
    call h5dget_type_f(setID,type_id,status)
    print *, type_id
    print *, 'H5T_NATIVE_INTEGER, H5T_NATIVE_CHARACTER, ' // &
       & 'H5T_NATIVE_DOUBLE, H5T_NATIVE_REAL, ' // &
       & 'H5T_STD_I32LE, H5T_IEEE_F32LE, H5T_IEEE_F64LE'
    print *, H5T_NATIVE_INTEGER, H5T_NATIVE_CHARACTER,  &
       & H5T_NATIVE_DOUBLE, H5T_NATIVE_REAL, &
       & H5T_STD_I32LE, H5T_IEEE_F32LE, H5T_IEEE_F64LE
    call h5tequal_f(type_id, H5T_IEEE_F32LE, is_hdf5, status)
    print *, 'are the 2 types the same? ', is_hdf5 
    call h5tequal_f(type_id, H5T_IEEE_F64LE, is_hdf5, status)
    print *, 'are the 2 dissimilar types the same? ', is_hdf5 
    call h5dClose_f ( setID, status )
     
       ! allocate( read_data%FrequencyCoordinates(max), stat=status)
       ! allocate( read_data%VerticalCoordinates(maxVC), stat=status)
       ! allocate( read_data%HorizontalCoordinates(maxHC), stat=status)
       allocate( read_data%RealField(dims(1),dims(2),dims(3)), stat=status) 
       call LoadFromHDF5DS(file_id, 'myDataset', read_data%RealField)
       print *,'read_data%RealField = ',read_data%RealField
       deallocate( read_data%RealField)

       print *, 'read just the 1st "MAF" '
!       allocate( read_data%RealField(1,1,1), stat=status) 
       allocate( read_data%RealField(dims(1),dims(2),1), stat=status) 
       call LoadFromHDF5DS(file_id, 'myDataset', read_data%RealField, &
         &  (/0,0,0/), (/int(dims(1:2)),1/) )
       print *,'read_data%RealField = ',read_data%RealField
       deallocate( read_data%RealField)

       print *, 'read just the last "MAF" '
!       allocate( read_data%RealField(1,1,1), stat=status) 
       allocate( read_data%RealField(dims(1),dims(2),1), stat=status) 
       call LoadFromHDF5DS(file_id, 'myDataset', read_data%RealField, &
         &  (/0,0,int(dims(3))-1/), (/int(dims(1:2)),1/) )
       print *,'read_data%RealField = ',read_data%RealField
       deallocate( read_data%RealField)
     CALL h5fclose_f(file_id, error)
!------------------------------------------------------------------

! Now try to read some data from l1b files
    print *, 'Now trying to read data from hdf5 l1b files'
    the_file_name = trim(l1b_path) // trim(l1boa_file)
    l1boa_id = mls_io_gen_openF('hg', .true., ErrType, &
    & record_length, DFACC_RDONLY, &
    & trim(the_file_name), hdfVersion=HDFVERSION_5)
    the_file_name = trim(l1b_path) // trim(l1bradd_file)
    l1brad(1) = mls_io_gen_openF('hg', .true., ErrType, &
    & record_length, DFACC_RDONLY, &
    & trim(the_file_name), hdfVersion=HDFVERSION_5)
    the_file_name = trim(l1b_path) // trim(l1bradf_file)
    l1brad(2) = mls_io_gen_openF('hg', .true., ErrType, &
    & record_length, DFACC_RDONLY, &
    & trim(the_file_name), hdfVersion=HDFVERSION_5)
    print *, 'ids of l1boa, l1bradd, l1bradf: ', l1boa_id, l1brad
    print *, '1st--l1boa'
    call ReadL1BData ( l1boa_id, '/sc/ECI', l1b_data, &
      & noMAFs, status, &  
      & firstMAF=21, lastMAF=21, &
      & NeverFail= .true., hdfVersion=HDFVERSION_5 )                                          
    if (status /= 0) print*, 'trouble reading l1bdata: ', status
    call dump ( l1b_data, 1)  
    print *, '2nd--l1brad'
    type_id = findl1bdata(l1brad, '/R2:190.B2F:H2O.S0.FB25-2', HDFVERSION_5)
    print *, 'file id: ', type_id
    call ReadL1BData ( type_id, '/R2:190.B2F:H2O.S0.FB25-2', l1b_data, &
      & noMAFs, status, &
      & firstMAF=21, lastMAF=21, &
      & NeverFail= .true., hdfVersion=HDFVERSION_5 )
    if (status /= 0) print*, 'trouble reading l1bdata: ', status
    call dump ( l1b_data, 1)
    print *, 'See if the following are in l1bradf'
    is_hdf5 = IsHDF5DSPresent(l1brad(2), '/counterMAF')
    print *, '/counterMAF ', is_hdf5
    is_hdf5 = IsHDF5DSPresent(l1brad(2), '/HDFEOS INFORMATION/coremetadata.0')
    print *, '/HDFEOS INFORMATION/coremetadata.0', is_hdf5
    is_hdf5 = IsHDF5DSPresent(l1brad(2), '/R1A:118.B1F:PT.S0.FB25-1')
    print *, '/R1A:118.B1F:PT.S0.FB25-1', is_hdf5
    status = mls_io_gen_closeF('hg', l1boa_id, hdfVersion=HDFVERSION_5 )
    status = mls_io_gen_closeF('hg', l1brad(1), hdfVersion=HDFVERSION_5 )
    status = mls_io_gen_closeF('hg', l1brad(2), hdfVersion=HDFVERSION_5 )
!------------------------------------------------------------------
     CALL h5close_f(error)

! Some tests of MLSStrings
     QName = ' GHz.tp ECI'
     print *, 'Before commpression: ', trim(QName)
     QName = CompressString(QName)
     print *, 'After commpression: ', trim(QName)
     do i=4, 5
       print *, 'HDF', i, ' l1b quantity names)'
       QName = AssembleL1BQtyName('R1A:118.B1F:PT.S0.FB25-1', i, .false.)
       print *, trim(QName)
       QName = AssembleL1BQtyName('scAngle', i, .false., 'GHz')
       print *, trim(QName)
       QName = AssembleL1BQtyName('scanRate', i, .false., 'GHz')
       print *, trim(QName)
       QName = AssembleL1BQtyName('ECI', i, .TRUE., 'GHz')
       print *, trim(QName)
       QName = AssembleL1BQtyName('GeocAlt', i, .TRUE., 'GHz')
       print *, trim(QName)
       QName = AssembleL1BQtyName('SolarZenith', i, .TRUE., 'GHz')
       print *, trim(QName)
       QName = AssembleL1BQtyName('scAngle', i, .false., 'THz')
       print *, trim(QName)
       QName = AssembleL1BQtyName('scanRate', i, .false., 'THz')
       print *, trim(QName)
       QName = AssembleL1BQtyName('ECI', i, .TRUE., 'THz')
       print *, trim(QName)
       QName = AssembleL1BQtyName('GeocAlt', i, .TRUE., 'THz')
       print *, trim(QName)
       QName = AssembleL1BQtyName('SolarZenith', i, .TRUE., 'THz')
       print *, trim(QName)
       QName = AssembleL1BQtyName('ECI', i, .false., 'sc')
       print *, trim(QName)
       QName = AssembleL1BQtyName('Vel', i, .false., 'sc')
       print *, trim(QName)
       QName = AssembleL1BQtyName('VelECI', i, .false., 'sc')
       print *, trim(QName)
     enddo
     STOP
     END PROGRAM MLS_ReadWritedata

! $Log$
