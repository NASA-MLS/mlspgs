! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!
! Converts each l1b dataset in a (list of) hdf4 file(s) to its
! counterpart in a (list of) HDF5 file(s).
! Usage: l1b_h4toh5 -h4 hdf4_file_name -h5 hdf5_file_name [..]
!
     PROGRAM l1b_h4toh5

     use MLSDataInfo, only: MLSDataInfo_T, Query_MLSData
     use MLSHDF5, only: SaveAsHDF5DS, mls_h5close, mls_h5open
     use HDF5, only: HID_T, H5F_ACC_TRUNC_F, h5fopen_f, h5fclose_f, &
       & h5fis_hdf5_f, h5gcreate_f, h5gclose_f
     use HDF, only: DFACC_RDONLY, &
       & DFNT_CHAR8, DFNT_INT32, DFNT_FLOAT64, &
       & sfselect, sfn2index, sfginfo, sfendacc, sffinfo, sfstart, sfend
     use L1BData, only: AssembleL1BQtyName, ReadL1BData, &
       & l1bdata_t
     use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
     use MLSFiles, only: split_path_name
     use MLSMessageModule, ONLY : MLSMessage, MLSMSG_Error, MLSMSG_Warning
     use MLSStrings, only: GetIntHashElement

     IMPLICIT NONE

     type(MLSDataInfo_T) :: dataset_info

     CHARACTER(LEN=255) :: hdf4filename, hdf5filename          ! filenames
     integer            :: hdf4_id, which_l1b
     integer            :: n_datasets, n_globattrs
     INTEGER(HID_T)     :: hdf5_id               ! File identifier
     INTEGER     ::  i, count, status, error ! Counting indices & Error flags
     INTEGER, PARAMETER ::  max_nsds = 1000  ! Maximum number of datasets in file.
     LOGICAL            :: is_hdf5
     CHARACTER(LEN=255) :: hdf4_path
     CHARACTER(LEN=128) :: hdf4_barefn, dummy
     integer            :: cm_index
     integer            :: cm_id
     integer            :: sds_id
     integer            :: rank
     integer            :: data_type
     integer            :: n_attrs
     integer            :: noMAF
     integer, dimension(3)            :: dim_sizes
     character(len=*), parameter :: file_types = 'radd, radf, boa'
!---------------------------- RCS Ident Info -------------------------------
  character(len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
     
  CALL mls_h5open(error)
  do      ! Loop over files to convert
     call get_filename(hdf4filename, hdf5filename)
     if ( hdf5filename == ' ' ) exit
     hdf4_id = sfstart(hdf4filename, DFACC_RDONLY)
     if ( hdf4_id == -1 ) then
       print *, 'Error starting hdf4 file: ', trim(hdf4filename)
       stop
     endif
     ! Check to see which of the l1b files it is: radd, radf, or boa
     call split_path_name(hdf4filename, hdf4_path, hdf4_barefn)
     which_l1b = GetIntHashElement(file_types, (/1,2,3/), &
       & hdf4_barefn, error, countEmpty=.true., part_match=.true.)
     if ( error /= 0 .or. which_l1b < 1 .or. which_l1b > 3) then
       print *, 'Unrecognized l1b file type: ', trim(hdf4_barefn)
       stop
     endif
     status = sffinfo(hdf4_id, n_datasets, n_globattrs)
     if ( status /= 0 ) then
       print *, 'Error querying hdf4 file: ', trim(hdf4_barefn)
       stop
     endif
     ! Get some important info about file so that we can allocate
     ! arrays with proper sizes
     cm_index = sfn2index(hdf4_id, 'counterMAF')
     if ( cm_index == -1 ) then
       print *, 'counterMAF field not found in hdf4 file: ', trim(hdf4_barefn)
       stop
     endif
     cm_id = sfselect(hdf4_id, cm_index)
     status = sfginfo ( cm_id, dummy, rank, dim_sizes, data_type, &
      n_attrs )
     noMAF = dim_sizes(1)
     print *, noMAF, ' MAFs in file ', trim( hdf4_barefn)
     status = sfendacc(cm_id)
     
     CALL h5fopen_f(trim(hdf5filename), H5F_ACC_TRUNC_F, hdf5_id, error)
     
     ! Now create some groups on the hdf5 file depending on the l1b type
     select case(which_l1b)
     case(1)
       call radd_group(hdf5_id)
     case(2)
       call radf_group(hdf5_id)
     case(3)
       call boa_group(hdf5_id)
     case default
       print*, 'Should not have gotten here'
       stop
     end select

     if (error .eq. 0) then
       do i=1, n_datasets
         sds_id = sfselect(hdf4_id, i)
         call convert_sds(sds_id, hdf4_id, hdf5_id, noMAF)
         status = sfendacc(sds_id)
       enddo
     endif ! if h5fopen_f
     error = sfend(hdf4_id)
     CALL h5fclose_f(hdf5_id, error)

  enddo        ! Loop over filenames

     CALL mls_h5close(error)
!-------------------------------------------------------------
     STOP
  contains

!------------------------- convert_sds ---------------------
    subroutine convert_sds(sds_id, hdf4_id, hdf5_id, NoMAF)
    ! read sds, convert to new format, write to hdf5
    integer, intent(in)             :: sds_id     
    integer, intent(in)             :: hdf5_id    
    integer, intent(in)             :: hdf4_id    
    integer, intent(in)             :: NoMAF      
    !
    type(l1bdata_t)                 :: L1BDATA    
    integer                         :: NoMAFsRead      
    integer                         :: status     
    integer                         :: rank       
    integer                         :: data_type  
    integer                         :: n_attrs    
    character(len=64)               :: sds_name   
    character(len=64)               :: sds5_name   
    integer, dimension(3)           :: dims       
    integer, dimension(3)           :: start      
    integer, dimension(3)           :: stride     
    integer, dimension(3)           :: edge
    integer, external               :: sfrcdata
    ! Executable
    edge = 1
    dims = 1
    status = sfginfo(sds_id, sds_name, rank, dims, data_type, n_attrs)
    if ( status /= 0 ) then
      print *, 'Trouble converting sds_id ', sds_id
      return
    endif    
    ! Set "slab" dimensions
    edge = dims(1:rank)
    start = 0
    stride = 1
    ! Any special cases go here
    select case(trim(sds_name))
    case('MAFStartTimeUTC')
      allocate(L1BDATA%CharField(1,1,1))
      status = sfrcdata(sds_id, start(1:1), stride(1:1), edge(1:1), &
        & L1BDATA%CharField)
      call SaveAsHDF5DS( hdf5_id, '/MAFStartTimeUTC', &
        & trim(L1BDATA%CharField(1,1,1)) )
      deallocate(L1BDATA%CharField)
    case default
      ! These are the standard cases:
      ! subdivided according to data type
      call ReadL1BData( hdf4_id, trim(sds_name), L1bData, NoMAFsRead, &  
        &  status, hdfVersion=4 )                                        
      if ( NoMAFsRead /= NoMAF ) then
        CALL MLSMessage(MLSMSG_Warning, ModuleName, & 
          & 'Number of mafs read different from expected number in sd ' // &
          & trim(sds_name))
      endif
      sds5_name = AssembleL1BQtyName(trim(sds_name), 5, .FALSE.)
      select case ( data_type )

      case ( DFNT_CHAR8 )
        print *, 'Sorry--unable to convert char array yet'
        stop
      case ( DFNT_INT32 )
      !  allocate(L1BDATA%IntField(1:dims(1),1:dims(2),1:dims(3)))
      !  status = sfrdata(sds_id, start(1:rank), stride(1:rank), edge(1:rank), &
      !    & L1BDATA%IntField)
        call SaveAsHDF5DS( hdf5_id, trim(sds5_name), &
          & L1BDATA%IntField )
        deallocate(L1BDATA%IntField)
      case ( DFNT_FLOAT64 )
      !  allocate(L1BDATA%IntField(1:dims(1),1:dims(2),1:dims(3)))
      !  status = sfrdata(sds_id, start(1:rank), stride(1:rank), edge(1:rank), &
      !    & L1BDATA%IntField)
        call SaveAsHDF5DS( hdf5_id, trim(sds5_name), &
          & L1BDATA%DpField )
        deallocate(L1BDATA%DpField)
      case default
      end select
    end select
  end subroutine convert_sds

!------------------------- get_filename ---------------------
    subroutine get_filename(hdf4filename, hdf5filename)
    ! Added for command-line processing
     CHARACTER(LEN=*), intent(out)    :: hdf4filename, hdf5filename ! filenames
     CHARACTER(LEN=len(hdf4filename)) :: filename
     integer ::                          error = 0
     integer, save ::                    i = 1
     integer, parameter               :: ihdf4 = 1
     integer, parameter               :: ihdf5 = ihdf4+1
     logical ::                          got(ihdf4:ihdf5)
  ! Get inputfile name, process command-line args
  ! (which always start with -)
  ! Proceed until you have a matched pair of hdf4 and hdf5 files
    got = .false.
    do while(.not. all(got))
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      if ( filename(1:1) /= '-' ) exit
      error = 1
      if ( filename(1:3) == '-h ' ) then
        call print_help
      else if ( filename(1:4) == '-h4 ' ) then
        call getarg ( i+1+hp, hdf4filename )
        error = 0
        i = i + 1
        got(ihdf4) = .true.
      else if ( filename(1:4) == '-h5 ' ) then
        call getarg ( i+1+hp, hdf5filename )
        error = 0
        i = i + 1
        got(ihdf5) = .true.
      else
        call print_help
      end if
      i = i + 1
    end do
    if ( error /= 0 ) then
      call print_help
    endif
    i = i + 1
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage: MLS_h5ls [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename => use filename'
      write (*,*) '          -h          => print brief help'
      stop
  end subroutine print_help
!------------------------- radd_group ---------------------
    subroutine radd_group(hdf5_id)
    ! Dummy arguments
    integer, intent(in)    :: hdf5_id
    
  end subroutine radd_group
!------------------------- radf_group ---------------------
    subroutine radf_group(hdf5_id)
    ! Dummy arguments
    integer, intent(in)    :: hdf5_id
    
  end subroutine radf_group
!------------------------- boa_group ---------------------
    subroutine boa_group(hdf5_id)
    ! Dummy arguments
    integer, intent(in)    :: hdf5_id
    
    ! Local parameters
    integer(hid_t)         :: group_id
    integer(hid_t)         :: subgroup_id
    integer                :: h5error
!
! Names of HDF5 Groups.
!
  CHARACTER(len=*), PARAMETER :: SC_GROUP_NAME    = 'sc'
  CHARACTER(len=*), PARAMETER :: GHZ_GROUP_NAME   = 'GHz'
  CHARACTER(len=*), PARAMETER :: THZ_GROUP_NAME   = 'THz'
  CHARACTER(len=*), PARAMETER :: TP_SUBGROUP_NAME = 'tp'
  CHARACTER(len=*), PARAMETER :: &
    H5_ERROR_GROUP_OPEN   = 'HDF5 Error Opening Group '
  CHARACTER(len=*), PARAMETER :: & 
    H5_ERROR_GROUP_CREATE = 'HDF5 Error Creating Group '
  CHARACTER(len=*), PARAMETER :: & 
    H5_ERROR_GROUP_CLOSE  = 'HDF5 Error Closing Group '

! Create the spacecraft group.
! 
     CALL h5gcreate_f(hdf5_id, SC_GROUP_NAME, group_id, h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_GROUP_CREATE // SC_GROUP_NAME)
!
! Close the spacecraft group.
!
     CALL h5gclose_f(group_id, h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_GROUP_CLOSE // SC_GROUP_NAME)
!
!-------------------------------------------------------------------
! 
! Open the GHz group.
!
     CALL h5gcreate_f(hdf5_id, GHZ_GROUP_NAME, group_id, h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_GROUP_CREATE // GHZ_GROUP_NAME)
!
! Create tangent point subgroup and the necessary interior dataspaces.
!
     CALL h5gcreate_f(group_id, TP_SUBGROUP_NAME, subgroup_id, & 
       h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_GROUP_CREATE // TP_SUBGROUP_NAME)
!
! Close the tangent point subgroup.
!
     CALL h5gclose_f(subgroup_id, h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_GROUP_CLOSE // TP_SUBGROUP_NAME)
!
! Close the GHz group.
!
     CALL h5gclose_f(group_id, h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_GROUP_CLOSE // GHZ_GROUP_NAME)
!
!-------------------------------------------------------------------
!
! Open and create the THz group.
!
     CALL h5gcreate_f(hdf5_id, THZ_GROUP_NAME, group_id, h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_GROUP_CREATE // THZ_GROUP_NAME)

!
! Create the THz tangent point subgroup.
!
     CALL h5gcreate_f(group_id, TP_SUBGROUP_NAME, subgroup_id, & 
       h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_GROUP_CREATE // TP_SUBGROUP_NAME)

!
! Close the THz/tp subgroup and the THz group.
!
     CALL h5gclose_f(subgroup_id, h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_GROUP_CLOSE // TP_SUBGROUP_NAME)

     CALL h5gclose_f(group_id, h5error)
     IF (h5error /= 0) CALL MLSMessage(MLSMSG_Error, ModuleName, & 
      H5_ERROR_GROUP_CLOSE // THZ_GROUP_NAME) 
  end subroutine boa_group
END PROGRAM l1b_h4toh5

! $Log$
! Revision 1.1  2002/10/29 00:56:57  pwagner
! First commit
!
