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
program l1bdump ! dumps an l1b or L2AUX file
!=================================

   use Dump_0, only: Dump
   use Dump_Options, only: Rmsformat
   use HDF, only: Dfacc_Read
   use HDF5, only: H5fis_HDF5_F, &
     & H5gclose_F, H5gopen_F
   use L1bData, only: L1bData_T, Namelen, &
     & DeallocateL1BData, Dump, ReadL1BData
   use Machine, only: Hp, Getarg
   use MLSFiles, only: Filenotfound, WildcardHDFversion, &
     & MLS_Exists, MLS_HDF_Version, MLS_Sfstart, MLS_Sfend, &
     & HDFversion_5
   use MLSHDF5, only: GetallHDF5dsnames, MLS_H5open, MLS_H5close
   use MLSMessagemodule, only: MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSStringlists, only: GetstringElement, NumstringElements
   use MLSStrings, only: Replace, Streq
   use Output_M, only: ResumeOutput, SuspendOutput, Output
   use Printit_M, only: Set_Config
   use Time_M, only: Time_Now, Time_Config
   
   implicit none

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

! Brief description of program
! dumps an l1b or l2aux file

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"
! Then run it
! LF95.Linux/test [options] [input files]

  type options_T
    logical     :: silent = .false.
    logical     :: timing = .false.
    logical     :: unique = .false.
    logical     :: debug = .false.
    logical     :: verbose = .false.
    logical     :: list = .false.
    logical     :: stats = .false.
    logical     :: rms = .false.
    logical     :: table = .false.
    logical     :: direct = .true.
    logical     :: l2aux = .false.
    logical     :: ascii = .false. ! If true, dump even character fields
    integer     :: hdfVersion = HDFVERSION_5
    integer     :: maf1 = 0
    integer     :: maf2 = 0
    integer     :: moff = 0
    character(len=255) :: group= ''     ! if SDs are within a group, group path
    character(len=255) :: sdList= '*'   ! what SDs to dump
    character(len=255) :: skipList= ''  ! what SDs to skip
    character(len=80)  :: dumpOptions       = ' '
  end type options_T
  
  type ( options_T ) ::          options
  integer, parameter ::          MAXDS = 300
  integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*namelen
  integer, parameter ::          MAXFILES = 100
  ! character(len=8)   ::          options%dumpOptions
  character(len=255) ::          filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer     ::                 i, status, error ! Counting indices & Error flags
  logical     ::                 is_hdf5
  integer ::                     maf1, maf2
  character (len=MAXSDNAMESBUFSIZE) :: mySdList
  integer            ::          n_filenames
  real        ::                 t1
  real        ::                 t2
  real        ::                 tFile
  ! 
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  time_config%use_wall_clock = .true.
  CALL mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_options(filename, n_filenames, options)
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     call h5fis_hdf5_f(trim(filename), is_hdf5, error)
     if ( options%hdfVersion==HDFVERSION_5 .and. .not. is_hdf5 ) then
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
  endif
  if ( options%rms ) rmsFormat = '(1pe9.2)'
  if ( options%silent ) call suspendOutput
  options%dumpOptions = '-'
  if ( options%rms ) options%dumpOptions = trim(options%dumpOptions) // 'r'
  if ( options%stats ) options%dumpOptions = trim(options%dumpOptions) // 's'
  if ( options%unique ) options%dumpOptions = trim(options%dumpOptions) // 'u'
  if ( options%silent ) options%dumpOptions = trim(options%dumpOptions) // 'h'
  if ( options%direct ) options%dumpOptions = trim(options%dumpOptions) // 'd'
  if ( options%table ) options%dumpOptions = trim(options%dumpOptions) // 'b'
  call time_now ( t1 )
  do i=1, n_filenames
    call time_now ( tFile )
    if ( options%list ) then
      print *, 'DS Names in: ', trim(filenames(i))
      if ( len_trim(options%group) < 1 ) then
        call GetAllHDF5DSNames (trim(filenames(i)), '/', mysdList)
      else
        call GetAllHDF5DSNames (trim(filenames(i)), trim(options%group), mysdList)
      endif
      call dump(mysdList, 'DS names')
    else 
      if ( options%verbose ) then
        print *, 'dumping from: ', trim(filenames(i))
      endif
      call mydump(trim(filenames(i)), &
      & options%hdfVersion, options)
      if ( options%timing ) call sayTime('dumping this file', tFile)
    endif
  enddo
  if ( options%timing ) call sayTime('dumping all files')
  call resumeOutput
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
      else if ( filename(1:4) == '-hdf' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%hdfVersion
        i = i + 1
        exit
      elseif ( filename(1:5) == '-opt ' ) then
        call getarg ( i+1+hp, options%dumpOptions )
        i = i + 1
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
      else if ( filename(1:2) == '-t' ) then
        options%table = .true.
        exit
      elseif ( filename(1:6) == '-skip ' ) then
        call getarg ( i+1+hp, options%skipList )
        i = i + 1
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-g ' ) then
        call getarg ( i+1+hp, options%group )
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
      & 'Usage:l1bdump [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options: -f filename => add filename to list of filenames'
      write (*,*) '                  (can do the same w/o the -f)'
      write (*,*) '   -d list     => just dump the SDs in list'
      write (*,*) '   -g group    => find SDs under group path'
      write (*,*) '   -ascii      => dump even character-valued fields'
      write (*,*) '   -l2aux      => the files are l2aux, not l1b'
      write (*,*) '   -v          => switch on verbose mode'
      write (*,*) '   -hdf version=> hdf version (default is 5)'
      write (*,*) '   -opt opts   => pass opts to dump routines'
      write (*,*) '   -silent     => switch on silent mode'
      write (*,*) '                 (printing only if dumps found)'
      write (*,*) '   -unique     => dump only unique elements'
      write (*,*) '   -l          => just list sd names in files'
      write (*,*) '   -maf m1,m2  => just dump in the range [m1,m2]'
      write (*,*) '   -moff offset=> 2nd data set starts after 1st'
      write (*,*) '   -rms        => just print mean, rms'
      write (*,*) '   -s          => just show number of dumperences'
      write (*,*) '   -skip list  => skip dumping the SDs in list'
      write (*,*) '   -t[able]    => table of % vs. magnitudes (pdf)'
      write (*,*) '   -h          => print brief help'
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

  ! ---------------------- mydump  ---------------------------
  subroutine mydump( file1, hdfVersion, options )
  !------------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine prints the dumps between them

    ! Arguments

    character (len=*), intent(in) :: file1 ! Name of file 1
    integer, intent(in)           :: hdfVersion
    type ( options_T )            :: options

    ! Local
    logical, parameter            :: countEmpty = .true.
    integer :: data_type
    integer, dimension(7) :: dimsizes
    logical :: file_exists
    integer :: grpid
    integer :: i
    type(l1bdata_t) :: L1BDATA  ! Result
    logical :: mustdump
    character(len=8) :: myOptions
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    integer :: NoMAFs
    integer :: noSds
    integer :: nsize
    integer :: num_Attrs
    integer :: rank
    integer :: sdfid1
    character (len=80) :: sdName
    integer :: sds_id
    real :: stime
    integer :: the_hdfVersion
    
    ! hdf4 externals
    integer :: sfstart, sffinfo, sfselect, sfginfo, sfend
    external :: sfstart, sffinfo, sfselect, sfginfo, sfend
    
    ! Executable code
    stime = t2
    ! the_hdfVersion = HDFVERSION_5
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
    if ( options%verbose ) print *, 'the hdf version: ', the_hdfVersion
    if ( the_hdfVersion == HDFVERSION_5 ) then
      if ( len_trim(options%group) < 1 ) then
        call GetAllHDF5DSNames (trim(File1), '/', mysdList)
      else
        call GetAllHDF5DSNames (trim(File1), trim(options%group), mysdList)
      endif
      if ( options%verbose ) then
        call output ( '============ DS names in ', advance='no' )
        call output ( trim(file1) //' ============', advance='yes' )
      endif
      if ( mysdList == ' ' .and. .not. options%silent ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'No way yet to find sdList in ' // trim(File1) )
        return
      else
        if ( options%verbose ) call dump(mysdList, 'DS names')
      endif
      if ( options%sdList /= '*' ) then
        call output( ' Selected SDs to dump', advance='yes' )
        mySDList = options%sdList
        call dump(mysdList, 'DS names')
      endif

      sdfid1 = mls_sfstart(File1, DFACC_READ, hdfVersion=hdfVersion)
      if (sdfid1 == -1 ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        &  'Failed to open l1b file ' // trim(File1) )
      end if
	   call h5gOpen_f (sdfid1,'/', grpID, status)
      if ( status /= 0 .and. .not. options%silent ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          	  & 'Unable to open group to read attribute in l2aux file' )
      endif
      noSds = NumStringElements(trim(mysdList), countEmpty)
      if ( noSds < 1 .and. .not. options%silent ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'No sdNames cp to file--unable to count sdNames in ' // trim(mysdList) )
      endif
    else
      sdfid1 = sfstart( trim(File1), DFACC_READ )
      status = sffinfo( sdfid1, noSds, nsize ) 
      ! print *, 'status ', status
      ! print *, 'noSds ', noSds
      mysdList = '*'
    endif
    ! Loop over sdNames in file 1
    ! (But skip PCF and HDFEOS INFORMATION/coremetadata.0)
    do i = 1, noSds
      if ( the_hdfversion == HDFVERSION_5 ) then
        call GetStringElement (trim(mysdList), sdName, i, countEmpty )
        if ( len_trim(options%group) > 1 ) sdName = trim(options%group) // sdName
      else
        sds_id = sfselect( sdfid1, i-1 ) ! "c" arrays begin with index 0
        status = sfginfo( sds_id, sdName, rank, dimsizes, data_type, num_attrs ) 
        ! print *, 'sdName: ', sdName
      endif
!      if ( sdName == 'PCF' .or. &
!        &  sdName == 'HDFEOS INFORMATION/coremetadata.0' .or. &
!        &  sdName == 'l2cf' .or. &
      if ( any( &
        & streq( &
        & (/ 'PCF ', 'meta', 'l2cf', 'utcp', 'leap', 'LCF ' /), &
        & sdname, options='-Pw' ) ) .or. &
        &  index(options%skipList, trim(sdName)) > 0 ) cycle
   !    print *, trim(sdName)
   !    print *, streq( '*PCF*', sdname, '-w' )
   !    print *, streq( '*meta*', sdname, '-w' )
   !    print *, streq( '*l2cf*', sdname, '-w' )
        ! Allocate and fill l2aux
      if ( options%debug ) print *, 'About to read ', trim(sdName)
      if ( options%maf2 > options%maf1 ) then
        call ReadL1BData ( sdfid1, trim(sdName), L1bData, NoMAFs, status, &
          & hdfVersion=the_hdfVersion, &
          & firstMAF=options%maf1, LastMAF=options%maf2, &
          & NEVERFAIL=.true., l2aux=options%l2aux )
      else
        call ReadL1BData ( sdfid1, trim(sdName), L1bData, NoMAFs, status, &
          & hdfVersion=the_hdfVersion, NEVERFAIL=.true., l2aux=options%l2aux )
      endif
      if ( status /= 0 .and. .not. options%silent ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
        call DeallocateL1BData ( l1bData )
        cycle
      endif
      if ( options%timing ) call SayTime( 'Reading l1bdata 1', stime )
      stime = t2

      if ( associated(L1bData%charField) .and. .not. options%ascii .and. &
        & .not. options%silent ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Skipping dump of char-valued ' // trim(sdName) )
        call DeallocateL1BData ( l1bData )
        cycle
      endif
      maf1 = 1
      if ( options%maf1 > 0 ) maf1 = options%maf1
      maf2 = NoMAFs
      if ( options%maf2 > 0 ) maf2 = options%maf2
      if ( options%verbose ) then
        print *, 'About to dump :: ' // trim(sdName) // ' ::'
      elseif ( .not. options%silent ) then
        print *, ':: ' // trim(sdName) // ' ::'
      endif
      mustdump = .true.
      myOptions = options%dumpOptions
      if ( mustdump ) myOptions = Replace ( myOptions, 'h', ' ' )
        call dump(L1bData, options=myOptions )
      if ( options%timing ) call SayTime( 'Doing the dump', stime )
      stime = t2
    enddo ! Loop of datasets
    if ( the_hdfVersion == HDFVERSION_5) call h5gClose_f (grpID, status)
    if ( status /= 0 .and. .not. options%silent ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
       & 'Unable to close group in l2aux file: ' // trim(File1) // ' after dumping' )
    endif
    if ( the_hdfVersion == HDFVERSION_5 ) then
	   status = mls_sfend(sdfid1, hdfVersion=the_hdfVersion)
      if ( status /= 0 ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "Unable to close L2aux file: " // trim(File1) // ' after dumping')
    else
      ! print *, 'Ending Mydump'
      status = sfend( sdfid1 )
    endif
  end subroutine mydump

!==================
end program l1bdump
!==================

! $Log$
! Revision 1.4  2016/10/04 22:22:38  pwagner
! Builds properly with some Dumps moved to Dump_1
!
! Revision 1.3  2014/03/07 21:43:51  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 1.2  2013/08/23 02:51:47  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.1  2013/06/01 00:40:26  pwagner
! First commit
!
