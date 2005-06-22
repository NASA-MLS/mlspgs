! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=================================
program l1bdiff ! diffs two l1b or L2AUX files
!=================================

   use Dump_0, only: DUMP
   use Hdf, only: DFACC_CREATE, DFACC_READ
   use HDF5, only: h5fis_hdf5_f, &
     & H5GCLOSE_F, H5GOPEN_F, H5DOPEN_F, H5DCLOSE_F, h5gcreate_f
   use L1BData, only: L1BData_T, NAME_LEN, &
     & DeallocateL1BData, Diff, ReadL1BData
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: FILENOTFOUND, WILDCARDHDFVERSION, &
     & mls_exists, mls_hdf_version, mls_sfstart, mls_sfend, &
     & MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
     & HDFVERSION_4, HDFVERSION_5
   use MLSHDF5, only: GetAllHDF5DSNames, saveAsHDF5DS, &
     & IsHDF5AttributePresent, mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSStringLists, only: GetStringElement, NumStringElements
   use output_m, only: output
   use Time_M, only: Time_Now, USE_WALL_CLOCK
   
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! Brief description of program
! diffs two l1b or l2aux files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"
! Then run it
! LF95.Linux/test [options] [input files]

  type options_T
    logical     :: verbose = .false.
    logical     :: list = .false.
    logical     :: stats = .false.
    logical     :: rms = .false.
    character(len=255) :: referenceFileName= 'default.h5'  ! reference filename
  end type options_T
  
  type ( options_T ) :: options


  integer, parameter ::          MAXDS = 50
  integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*NAME_LEN
  integer, parameter ::          MAXFILES = 100
  logical ::          columnsOnly
  character(len=255) :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  i, count, status, error ! Counting indices & Error flags
  logical     :: is_hdf5
  character (len=MAXSDNAMESBUFSIZE) :: mySdList
  real        :: t1
  real        :: t2
  real        :: tFile
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  USE_WALL_CLOCK = .true.
  CALL mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename(filename, n_filenames, options)
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     call h5fis_hdf5_f(trim(filename), is_hdf5, error)
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no input files supplied'
    stop
  elseif ( options%referenceFileName == 'default.h5' ) then
    options%referenceFileName = filenames(n_filenames)
    n_filenames = n_filenames - 1
  endif
  call time_now ( t1 )
  if ( options%verbose .and. .not. options%list) &
    & print *, 'Compare l1b data to: ', trim(options%referenceFileName)
  do i=1, n_filenames
    call time_now ( tFile )
    if ( options%list ) then
      print *, 'DS Names in: ', trim(filenames(i))
      call GetAllHDF5DSNames (trim(filenames(i)), '/', mysdList)
      call dump(mysdList, 'DS names')
    else 
      if ( options%verbose ) then
        print *, 'diffing from: ', trim(filenames(i))
      endif
      call myDiff(trim(filenames(i)), &
      & trim(options%referenceFileName), (i==1), &
      & HDFVERSION_5, options)
      call sayTime('diffing this file', tFile)
    endif
  enddo
  if ( .not. options%list) call sayTime('diffing all files')
  call mls_h5close(error)
contains
!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(LEN=255), intent(out)       :: filename          ! filename
     integer, intent(in)                   :: n_filenames
     type ( options_T ), intent(inout) :: options
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
     character(LEN=16)       :: number
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:3) == '-r ' ) then
        call getarg ( i+1+hp, options%referenceFileName )
        i = i + 1
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:3) == '-l ' ) then
        options%list = .true.
        exit
      else if ( filename(1:5) == '-rms ' ) then
        options%rms = .true.
        exit
      else if ( filename(1:3) == '-s ' ) then
        options%stats = .true.
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else
        call print_help
      end if
      i = i + 1
    end do
    if ( error /= 0 ) then
      call print_help
    endif
    i = i + 1
    if (trim(filename) == ' ' .and. n_filenames == 0) then

    ! Last chance to enter filename
      print *,  "Enter the name of the HDF5 l1b or l2aux file. "
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:l1bdiff [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename     => add filename to list of filenames'
      write (*,*) '                           (can do the same w/o the -f)'
      ! write (*,*) '          -r reffile      => compare sds to reffile'
      write (*,*) '          -v              => switch on verbose mode'
      write (*,*) '          -l              => just list sd names in files'
      write (*,*) '          -rms            => just print mean, rms'
      write (*,*) '          -s              => just show statistics'
      write (*,*) '          -h              => print brief help'
      stop
  end subroutine print_help
!------------------------- SayTime ---------------------
  subroutine SayTime ( What, startTime )
    character(len=*), intent(in) :: What
    real, intent(in), optional :: startTime
    real :: myt1
    if ( present(startTime) ) then
      myt1 = startTime
    else
      myt1 = t1
    endif
    call time_now ( t2 )
    call output ( "Timing for " // what // " = " )
    call output ( dble(t2 - myt1), advance = 'yes' )
  end subroutine SayTime

  ! ---------------------- myDiff  ---------------------------
  subroutine myDiff(file1, file2, create2, hdfVersion, options)
  !------------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine prints the diffs between them

    ! Arguments

    character (len=*), intent(in) :: file1 ! Name of file 1
    character (len=*), intent(in) :: file2 ! Name of file 2
    logical, intent(in) :: create2
    integer, intent(in) :: hdfVersion
    type ( options_T ), intent(in) :: options

    ! Local
    integer :: sdfid1
    integer :: sdfid2
    integer :: grpid
    integer :: status
    integer :: the_hdfVersion
    logical :: file_exists
    integer :: file_access
    integer :: noSds
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    logical, parameter            :: countEmpty = .true.
    logical :: isl1boa
    ! type (L2AUXData_T) :: l2aux
    type(l1bdata_t) :: L1BDATA  ! Result
    type(l1bdata_t) :: L1BDATA2 ! for diff
    integer :: i
    integer :: NoMAFs
    character (len=80) :: sdName
    
    ! Executable code
    the_hdfVersion = HDFVERSION_5
    the_hdfVersion = hdfVersion
    file_exists = ( mls_exists(trim(File1)) == 0 )
    if ( .not. file_exists ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File 1 not found; make sure the name and path are correct' &
        & // trim(file1) )
    endif
    if ( the_hdfVersion == WILDCARDHDFVERSION ) then
      the_hdfVersion = mls_hdf_version(File1, hdfVersion)
      if ( the_hdfVersion == FILENOTFOUND ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'File 1 not found; make sure the name and path are correct' &
          & // trim(file1) )
    endif
    call GetAllHDF5DSNames (trim(File1), '/', mysdList)
    if ( options%verbose) then
      call output ( '============ DS names in ', advance='no' )
      call output ( trim(file1) //' ============', advance='yes' )
    endif
    if ( mysdList == ' ' ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'No way yet to find sdList in ' // trim(File1) )
      return
    else
      if ( options%verbose ) call dump(mysdList, 'DS names')
    endif

    isl1boa = (index(trim(mysdList), '/GHz') > 0)
    file_access = DFACC_READ
    sdfid1 = mls_sfstart(File1, DFACC_READ, hdfVersion=hdfVersion)
    if (sdfid1 == -1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      &  'Failed to open l1b file ' // trim(File1) )
    end if
	 call h5gOpen_f (sdfid1,'/', grpID, status)
    if ( status /= 0 ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to open group to read attribute in l2aux file' )
    endif
    sdfId2 = mls_sfstart(trim(file2), file_access, &
              & hdfVersion=hdfVersion)
    if (sdfid2 == -1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Failed to open l1b ' // trim(File2) )
    end if
    noSds = NumStringElements(trim(mysdList), countEmpty)
    if ( noSds < 1 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'No sdNames cp to file--unable to count sdNames in ' // trim(mysdList) )
    endif
    ! May need to create groups '/sc', '/GHz', and '/THz'
    if ( isL1boa .and. file_access == DFACC_CREATE ) then
      call h5gcreate_f(sdfId2, '/sc', grpID, status)
      if ( status /= 0 ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to create group /sc in ' // trim(File2) )
      endif
      call h5gcreate_f(sdfId2, '/GHz', grpID, status)
      if ( status /= 0 ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to create group /GHz in ' // trim(File2) )
      endif
      call h5gcreate_f(sdfId2, '/THz', grpID, status)
      if ( status /= 0 ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to create group /THz in ' // trim(File2) )
      endif
    endif
    ! Loop over sdNames in file 1
    do i = 1, noSds
      call GetStringElement (trim(mysdList), sdName, i, countEmpty )
      ! Allocate and fill l2aux
      if ( options%verbose ) print *, 'About to read ', trim(sdName)
        call ReadL1BData ( sdfid1, trim(sdName), L1bData, NoMAFs, status, &
          & hdfVersion=the_hdfVersion, NEVERFAIL=.true. )
        if ( status /= 0 ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
          cycle
        endif
        call ReadL1BData ( sdfid2, trim(sdName), L1bData2, NoMAFs, status, &
          & hdfVersion=the_hdfVersion, NEVERFAIL=.true. )
        if ( status /= 0 ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
          call DeallocateL1BData ( l1bData )
          cycle
        endif
        call diff(L1bData, L1bData2, stats=options%stats, rms=options%rms)
        call DeallocateL1BData ( l1bData )
        call DeallocateL1BData ( l1bData2 )
    enddo
	 call h5gClose_f (grpID, status)
    if ( status /= 0 ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
       & 'Unable to close group in l2aux file: ' // trim(File1) // ' after diffing' )
    endif
	 status = mls_sfend(sdfid1, hdfVersion=the_hdfVersion)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2aux file: " // trim(File1) // ' after diffing')
	 status = mls_sfend(sdfid2, hdfVersion=the_hdfVersion)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2aux file: " // trim(File2) // ' after diffing')
  end subroutine myDiff

!==================
end program l1bdiff
!==================

! $Log$
