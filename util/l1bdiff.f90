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

   use Dump_0, only: DIFFRMSMEANSRMS, rmsFormat, DIFF, DUMP, SELFDIFF
   use Hdf, only: DFACC_CREATE, DFACC_READ
   use HDF5, only: h5fis_hdf5_f, &
     & H5GCLOSE_F, H5GOPEN_F, H5DOPEN_F, H5DCLOSE_F, h5gcreate_f
   use L1BData, only: L1BData_T, NAME_LEN, &
     & DeallocateL1BData, Diff, ReadL1BData
   use MACHINE, only: FILSEP, HP, IO_ERROR, GETARG
   use MLSCommon, only: R8
   use MLSFiles, only: FILENOTFOUND, WILDCARDHDFVERSION, &
     & mls_exists, mls_hdf_version, mls_sfstart, mls_sfend, &
     & HDFVERSION_4, HDFVERSION_5
   use MLSHDF5, only: GetAllHDF5DSNames, saveAsHDF5DS, &
     & IsHDF5AttributePresent, mls_h5open, mls_h5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSStringLists, only: GetStringElement, NumStringElements
   use MLSStrings, only: WriteIntsToChars
   use output_m, only: resumeOutput, suspendOutput, output, outputNamedValue
   use Time_M, only: Time_Now, time_config
   
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
    logical     :: halfwaves = .false.
    logical     :: self = .false.
    logical     :: silent = .false.
    logical     :: timing = .false.
    logical     :: unique = .false.
    logical     :: verbose = .false.
    logical     :: list = .false.
    logical     :: stats = .false.
    logical     :: rms = .false.
    logical     :: direct = .false.
    logical     :: l2aux = .false.
    logical     :: ascii = .false. ! If true, diff even character fields
    integer     :: maf1 = 0
    integer     :: maf2 = 0
    integer     :: moff = 0
    integer     :: numDiffs = 0
    character(len=255) :: sdList= '*'  ! what SDs to diff
    character(len=255) :: skipList= ''  ! what SDs to skip
    character(len=255) :: referenceFileName= 'default.h5'  ! reference filename
  end type options_T
  
  type ( options_T ) ::          options
  integer, parameter ::          MAXDS = 300
  integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*NAME_LEN
  integer, parameter ::          MAXFILES = 100
  logical ::                     columnsOnly
  character(len=8)   ::          dumpOptions
  character(len=255) ::          filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer     ::                 i, status, error ! Counting indices & Error flags
  logical     ::                 is_hdf5
  integer ::                     maf1, maf2
  character (len=MAXSDNAMESBUFSIZE) :: mySdList
  integer            ::          n_filenames
  character(len=16) ::           string
  real        ::                 t1
  real        ::                 t2
  real        ::                 tFile
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  time_config%use_wall_clock = .true.
  DIFFRMSMEANSRMS = .true.
  CALL mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_options(filename, n_filenames, options)
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
  elseif ( options%verbose .and. options%silent ) then
    print *, 'Sorry--either verbose or silent; cant be both'
    stop
  elseif ( options%self ) then
    ! We'll do diffs within each file
  elseif ( options%referenceFileName == 'default.h5' ) then
    options%referenceFileName = filenames(n_filenames)
    n_filenames = n_filenames - 1
  endif
  if ( options%rms ) rmsFormat = '(1pe9.2)'
  if ( options%silent ) call suspendOutput
  dumpOptions = '-'
  if ( options%rms ) dumpOptions = trim(dumpOptions) // 'r'
  if ( options%stats ) dumpOptions = trim(dumpOptions) // 's'
  if ( options%unique ) dumpOptions = trim(dumpOptions) // 'u'
  call time_now ( t1 )
  if ( options%verbose .and. .not. options%list ) &
    & print *, 'Compare l1b data to: ', trim(options%referenceFileName)
  do i=1, n_filenames
    call time_now ( tFile )
    if ( options%list ) then
      print *, 'DS Names in: ', trim(filenames(i))
      call GetAllHDF5DSNames (trim(filenames(i)), '/', mysdList)
      call dump(mysdList, 'DS names')
    elseif ( options%self )then
      if ( options%verbose ) then
        print *, 'diffing from: ', trim(filenames(i))
      endif
      call mySelfDiff(trim(filenames(i)), HDFVERSION_5, options)
      if ( options%timing ) call sayTime('diffing this file', tFile)
    else 
      if ( options%verbose ) then
        print *, 'diffing from: ', trim(filenames(i))
      endif
      call myDiff(trim(filenames(i)), &
      & trim(options%referenceFileName), (i==1), &
      & HDFVERSION_5, options)
      if ( options%timing ) call sayTime('diffing this file', tFile)
    endif
  enddo
  if ( options%timing ) call sayTime('diffing all files')
  call resumeOutput
  if ( options%silent .and. options%numDiffs > 0 ) then
    ! write(string, '(i)') options%numDiffs
    call WriteIntsToChars ( options%numDiffs, string )
    call print_string(string)
  endif
  call mls_h5close(error)
contains
!------------------------- get_options ---------------------
    subroutine get_options( filename, n_filenames, options )
    ! Added for command-line processing of options, filenames
     character(LEN=255), intent(out)       :: filename          ! filename
     integer, intent(in)                   :: n_filenames
     type ( options_T ), intent(inout)     :: options
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
      elseif ( filename(1:5) == '-maf ' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%maf1, options%maf2
        i = i + 1
        exit
      elseif ( filename(1:6) == '-moff ' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%moff
        i = i + 1
        exit
      elseif ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, options%sdList )
        i = i + 1
        exit
      elseif ( filename(1:3) == '-r ' ) then
        call getarg ( i+1+hp, options%referenceFileName )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-half' ) then
        options%halfWaves = .true.
        exit
      elseif ( filename(1:6) == '-self ' ) then
        options%self = .true.
        exit
      elseif ( filename(1:8) == '-silent ' ) then
        options%silent = .true.
        exit
      elseif ( filename(1:8) == '-unique ' ) then
        options%unique = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:7) == '-ascii ' ) then
        options%ascii = .true.
        exit
      elseif ( filename(1:7) == '-l2aux ' ) then
        options%l2aux = .true.
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
      elseif ( filename(1:6) == '-skip ' ) then
        call getarg ( i+1+hp, options%skipList )
        i = i + 1
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
    
  end subroutine get_options
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:l1bdiff [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename     => add filename to list of filenames'
      write (*,*) '                           (can do the same w/o the -f)'
      write (*,*) '          -d list         => just diff the SDs in list'
      write (*,*) '          -ascii          => diff even character-valued fields'
      write (*,*) '          -l2aux          => the files are l2aux, not l1b'
      write (*,*) '          -v              => switch on verbose mode'
      write (*,*) '          -self           => dump successive differences'
      write (*,*) '                             between values in same file'
      write (*,*) '          -half           => show no. of 1/2 waves'
      write (*,*) '          -silent         => switch on silent mode'
      write (*,*) '                            (printing only if diffs found)'
      write (*,*) '          -unique         => dump only unique elements'
      write (*,*) '          -l              => just list sd names in files'
      write (*,*) '          -maf m1,m2      => just diff in the range [m1,m2]'
      write (*,*) '          -moff offset    => 2nd data set starts after 1st'
      write (*,*) '          -rms            => just print mean, rms'
      write (*,*) '          -s              => just show statistics'
      write (*,*) '          -skip list      => skip diffing the SDs in list'
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
  subroutine myDiff( file1, file2, create2, hdfVersion, options )
  !------------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine prints the diffs between them

    ! Arguments

    character (len=*), intent(in) :: file1 ! Name of file 1
    character (len=*), intent(in) :: file2 ! Name of file 2
    logical, intent(in)           :: create2
    integer, intent(in)           :: hdfVersion
    type ( options_T )            :: options

    ! Local
    logical, parameter            :: countEmpty = .true.
    logical :: file_exists
    integer :: file_access
    integer :: grpid
    integer :: i
    logical :: isl1boa
    type(l1bdata_t) :: L1BDATA  ! Result
    type(l1bdata_t) :: L1BDATA2 ! for diff
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    integer :: NoMAFs
    integer :: noSds
    integer :: numDiffs
    integer :: sdfid1
    integer :: sdfid2
    character (len=80) :: sdName
    integer :: status
    integer :: the_hdfVersion
    
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
    if ( options%sdList /= '*' ) then
      call output( ' Selected SDs to diff', advance='yes' )
      mySDList = options%sdList
      call dump(mysdList, 'DS names')
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
    ! (But skip PCF and HDFEOS INFORMATION/coremetadata.0)
    do i = 1, noSds
      call GetStringElement (trim(mysdList), sdName, i, countEmpty )
      if ( sdName == 'PCF' .or. &
        &  sdName == 'HDFEOS INFORMATION/coremetadata.0' .or. &
        &  sdName == 'l2cf' .or. &
        &  index(options%skipList, trim(sdName)) > 0 ) cycle
      ! Allocate and fill l2aux
      ! if ( options%verbose ) print *, 'About to read ', trim(sdName)
        call ReadL1BData ( sdfid1, trim(sdName), L1bData, NoMAFs, status, &
          & hdfVersion=the_hdfVersion, NEVERFAIL=.true., l2aux=options%l2aux )
        if ( status /= 0 ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
          call DeallocateL1BData ( l1bData )
          cycle
        endif
      ! if ( options%verbose ) print *, 'About to read ', trim(sdName), ' (2nd)'
        call ReadL1BData ( sdfid2, trim(sdName), L1bData2, NoMAFs, status, &
          & hdfVersion=the_hdfVersion, NEVERFAIL=.true., l2aux=options%l2aux )
        if ( status /= 0 ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
          call DeallocateL1BData ( l1bData )
          call DeallocateL1BData ( l1bData2 )
          cycle
        endif
        if ( associated(L1bData%charField) .and. .not. options%ascii ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	& 'Skipping diff of char-valued ' // trim(sdName) )
          call DeallocateL1BData ( l1bData )
          call DeallocateL1BData ( l1bData2 )
          cycle
        endif
        maf1 = 1
        if ( options%maf1 > 0 ) maf1 = options%maf1
        maf2 = NoMAFs
        if ( options%maf2 > 0 ) maf2 = options%maf2
        ! if ( options%verbose ) print *, 'About to diff'
        if ( .not. options%silent ) print *, 'About to diff ', trim(sdName)
        if ( options%direct .or. .not. associated(L1bData%dpField) ) then
          ! print *, 'About to do it direct'
          call diff(L1bData, L1bData2, numDiffs=numDiffs, options=dumpOptions )
        else
          ! print *, 'details=0'
          call diff(L1bData, L1bData2, details=0, &
            & numDiffs=numDiffs, options=dumpOptions )
          ! print *, 'done diffing'
          
          numDiffs = numDiffs + count( L1bData%dpField /= L1bData2%dpField )
          if ( .true. .and. associated(L1bData%dpField) .and. &
            & associated(L1bData2%dpField)) then
!             call outputNamedValue( 'maf1', maf1 )
!             call outputNamedValue( 'maf2', maf2 )
!             call outputNamedValue( 'options%moff', options%moff )
!             print *, shape(L1bData%dpField(:,:,maf1:maf2))
!             print *, shape(L1bData2%dpField(:,:,maf1+options%moff:maf2+options%moff))
!             print *, L1bData%dpField(1,1,maf1)
!             print *, L1bData%dpField(1,1,maf2)
!             print *, L1bData2%dpField(1,1,maf1+options%moff)
!             print *, L1bData2%dpField(1,1,maf2+options%moff)
            if ( options%silent ) then
            elseif ( .false. .and. all(L1bData%dpField(:,:,maf1:maf2) == &
              & L1bData2%dpField(:,:,maf1+options%moff:maf2+options%moff)) ) then
              print *, '(The two fields are exactly equal)'
            elseif( options%moff /= 0 .or. options%maf1 > 0 ) then
              ! print *, 'About to diff dpFields'
              call diff( L1bData%dpField(:,:,maf1:maf2), &
                & '(1)', &
                & L1bData2%dpField(:,:,maf1+options%moff:maf2+options%moff), &
                & '(2)', &
                & options=dumpOptions )
            else
              call diff( L1bData%dpField, &
                & '(1)', &
                & L1bData2%dpField, &
                & '(2)', &
                & options=dumpOptions )
            endif
          else
            ! print *, 'About to form dpField = dpField1 - dpField2'
            L1bData%dpField = L1bData%dpField - L1bData2%dpField(:,:,1+options%moff:)
            if ( options%silent ) then
            elseif ( .false. .and. all(L1bData%dpField == 0.d0) ) then
              print *, '(The two fields are exactly equal)'
            elseif ( options%maf1 /= 0 .and. options%maf2 /= 0 ) then
              ! print *, shape( L1bData%dpField(:,:,options%maf1:options%maf2) )
              call dump( L1bData%dpField(:,:,options%maf1:options%maf2), &
                & 'l1bData%dpField', options=dumpOptions )
            else
              ! print *, shape( L1bData%dpField )
              call dump( L1bData%dpField, &
                & 'l1bData%dpField', options=dumpOptions )
            endif
          endif
        endif
        call DeallocateL1BData ( l1bData )
        call DeallocateL1BData ( l1bData2 )
        options%numDiffs = options%numDiffs + numDiffs
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

  ! ---------------------- mySelfDiff  ---------------------------
  subroutine mySelfDiff( file1, hdfVersion, options )
  !------------------------------------------------------------------------

    ! Given file name file1
    ! This routine prints its selfdiffs

    ! Arguments

    character (len=*), intent(in) :: file1 ! Name of file 1
    integer, intent(in)           :: hdfVersion
    type ( options_T )            :: options

    ! Local
    logical, parameter            :: countEmpty = .true.
    logical :: file_exists
    integer :: file_access
    integer :: grpid
    integer :: i
    logical :: isl1boa
    type(l1bdata_t) :: L1BDATA  ! Result
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    integer :: NoMAFs
    integer :: noSds
    integer :: numDiffs
    integer :: sdfid1
    character (len=80) :: sdName
    integer :: status
    integer :: the_hdfVersion
    
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
    if ( options%sdList /= '*' ) then
      call output( ' Selected SDs to diff', advance='yes' )
      mySDList = options%sdList
      call dump(mysdList, 'DS names')
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
    noSds = NumStringElements(trim(mysdList), countEmpty)
    if ( noSds < 1 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'No sdNames cp to file--unable to count sdNames in ' // trim(mysdList) )
    endif
    ! Loop over sdNames in file 1
    ! (But skip PCF and HDFEOS INFORMATION/coremetadata.0)
    do i = 1, noSds
      call GetStringElement (trim(mysdList), sdName, i, countEmpty )
      if ( sdName == 'PCF' .or. &
        &  sdName == 'HDFEOS INFORMATION/coremetadata.0' .or. &
        &  sdName == 'l2cf' .or. &
        &  index(options%skipList, trim(sdName)) > 0 ) cycle
      ! Allocate and fill l2aux
      ! if ( options%verbose ) print *, 'About to read ', trim(sdName)
      call ReadL1BData ( sdfid1, trim(sdName), L1bData, NoMAFs, status, &
        & hdfVersion=the_hdfVersion, NEVERFAIL=.true., l2aux=options%l2aux )
      if ( status /= 0 ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
        call DeallocateL1BData ( l1bData )
        cycle
      endif
      if ( associated(L1bData%charField) ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Skipping diff of char-valued ' // trim(sdName) )
      elseif ( associated(L1bData%intField) ) then
        if ( size(L1bData%intField, 1) > 1 .or. &
          &  size(L1bData%intField, 2) > 1 ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Skipping diff of char-valued ' // trim(sdName) )
        else
          call selfdiff( L1bData%intField(1, 1, :), trim(sdName), &
            & waves=options%halfWaves, options=dumpOptions )
        endif
      elseif ( associated(L1bData%dpField) ) then
        if ( size(L1bData%dpField, 1) > 1 .or. &
          &  size(L1bData%dpField, 2) > 1 ) then
	       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Skipping diff of char-valued ' // trim(sdName) )
        else
          call selfdiff( L1bData%dpField(1, 1, :), trim(sdName), &
            & waves=options%halfWaves, options=dumpOptions )
        endif
      endif
      call DeallocateL1BData ( l1bData )
      cycle
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
  end subroutine mySelfDiff

!------------------------- print_string ---------------------
  subroutine print_string(string)
    character(len=*), intent(in) :: string
    write(*,'(a)') trim(string)
  end subroutine print_string

!==================
end program l1bdiff
!==================

! $Log$
! Revision 1.11  2008/07/10 00:14:43  pwagner
! Added -self, -half, -unique options
!
! Revision 1.10  2008/04/10 20:22:37  pwagner
! Less voluminous output
!
! Revision 1.9  2008/02/28 01:34:07  pwagner
! Normally should skip diffing character-valued fields
!
! Revision 1.8  2007/11/28 19:35:21  pwagner
! May specify that files are l2aux
!
! Revision 1.7  2007/07/18 00:15:25  pwagner
! -rms outputs fractional diffs more clearly
!
! Revision 1.6  2006/11/22 18:33:10  pwagner
! New optional args to diff l1b files with different number MAFs
!
! Revision 1.5  2006/06/14 16:42:38  pwagner
! Should not run out of memory unless direct reset to TRUE
!
! Revision 1.4  2006/01/14 00:58:31  pwagner
! Added -silent option
!
! Revision 1.3  2005/10/29 00:13:56  pwagner
! Removed unused procedures from use statements
!
! Revision 1.2  2005/09/23 21:01:13  pwagner
! use_wall_clock now a component of time_config
!
! Revision 1.1  2005/06/22 19:27:32  pwagner
! Reworded Copyright statement, moved rcs id
!
