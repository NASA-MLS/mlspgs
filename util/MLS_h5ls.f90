! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!
! JPD:: This program iterates through the groups and subgroups of an HDF5 file 
!       and lists the full path of each dataset in that HDF5 file.
!
     PROGRAM MLS_h5ls

     use Hdf, only: DFACC_RDONLY
     use L1BData, only: L1BData_T, ReadL1BData, deallocateL1BData
     use MLSDataInfo, only: MLSDataInfo_T, Query_MLSData
     use MLSFiles, only: HDFVERSION_5, mls_sfend, mls_sfstart
     use MLSFillValues, only: isFinite
     use MLSHDF5, only: mls_h5close, mls_h5open
     use HDF5, only: HID_T, H5F_ACC_RDONLY_F, h5fopen_f, h5fclose_f, &
       & h5fis_hdf5_f   
     use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG

     IMPLICIT NONE
     logical :: checkForNaNs

     type(MLSDataInfo_T) :: dataset_info

     CHARACTER(LEN=255) :: filename          ! filename
     integer            :: n_filenames
     INTEGER(HID_T) :: file_id               ! File identifier
     INTEGER     ::  i, count, status, error ! Counting indices & Error flags
     INTEGER, PARAMETER ::  max_nsds = 10000  ! Maximum number of datasets in file.
     LOGICAL     :: is_hdf5
     type(l1bdata_t) :: HDFDATA
     integer :: L1FileHandle
     integer :: noMAFs
     logical :: verbose

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------
     
     ! print *, 'Your file name: ', filename

     CALL mls_h5open(error)
  n_filenames = 0
  checkForNaNs = .false.
  verbose = .false.
  do      ! Loop over filenames
     call get_filename(filename, n_filenames)
     if ( filename == ' ' ) exit
     n_filenames = n_filenames + 1
     call h5fis_hdf5_f(trim(filename), is_hdf5, error)

     if (is_hdf5) then 

     CALL h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

     IF (error .eq. 0) then  
     if ( verbose ) then
       write (*, *) ' '
       write (*, *) '------------------------------------------------------'
       write (*, *) 'file: ', trim(filename)
       write (*, *) '------------------------------------------------------'
     endif
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
!        write(*,101) count, filename
! 101    format(1x,"There are ",1i5," dataset entries in the HDF5 file, ",1a64)
        write(*,*) ' There are ', count, &
          & ' dataset entries in the HDF5 file, ', trim(filename)

     do i = 1, count
        write(*,102) dataset_info%name(i) 
 102    format(1x,1a64)
        if ( .not. checkForNaNs ) cycle
        ! L1FileHandle = mls_sfstart( trim(filename), &
        !   & DFACC_RDONLY, HDFVERSION_5 )
        call ReadL1BData ( file_id, trim(dataset_info%name(i)), HDFDATA, &
          & noMAFs, status, NEVERFAIL=.true., L2AUX=.true. )
        if ( associated(HDFDATA%DpField) ) then
          if ( .not. all( isFinite(HDFDATA%DpField) ) ) then
            write(*, *) '**** NaNs or Infs found in dataset'
          endif
        else
            write(*, *) '           (Unable to check this dataset)'
        endif
        call deallocateL1BData( HDFDATA )
        ! status = mls_sfend( L1FileHandle, HDFVersion_5 )
        
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

  enddo        ! Loop over filenames

     CALL mls_h5close(error)
!-------------------------------------------------------------
     STOP
  contains

!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames)
    ! Added for command-line processing
     CHARACTER(LEN=255), intent(out) :: filename          ! filename
     integer, intent(in) ::             n_filenames
     integer ::                         error = 0
     integer, save ::                   i = 1
  ! Get inputfile name, process command-line args
  ! (which always start with -)
!    i = 1
!    error = 0
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      if ( filename(1:1) /= '-' ) exit
      error = 1
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:6) == '-check' ) then
        error = 0
        checkForNaNs = .true.
      elseif ( filename(1:2) == '-v' ) then
        error = 0
        verbose = .true.
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        error = 0
        i = i + 1
        exit
      else
        call print_help
      end if
      i = i + 1
    end do
    ! print *, error, ' ', trim(filename)
    if ( error /= 0 ) then
      call print_help
    endif
    i = i + 1
    if (trim(filename) == ' ' .and. n_filenames == 0) then

    ! Last chance to enter filename
      print *,  "Enter the name of the HDF5 file. " // &
       &  "Datasets in the file will be listed shortly."
      read(*,'(a)') filename
    endif
    
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
END PROGRAM MLS_h5ls

! $Log$
! Revision 1.5  2005/06/22 19:27:33  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.4  2004/03/25 18:39:28  pwagner
! Uses mls_h5open/close
!
! Revision 1.3  2002/10/03 23:25:41  pwagner
! Less likely to snip part of path/filename in printing
!
! Revision 1.2  2002/08/29 17:53:55  pwagner
! Now takes filenames on command line; Lahey has another internal error while compiling
!
! Revision 1.1  2002/08/28 22:49:21  pwagner
! First commit
!
