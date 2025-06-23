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
program l2auxdump ! dumps datasets, attributes from L2AUX files
!=================================

   use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
   use Dates_Module, only: Tai93s2hid
   use Dump_1, only: Dump
   use Dump_Options, only: DefaultMaxLon, DefaultWidth, DumpDumpOptions, &
     & IntPlaces
   use HDF, only: Dfacc_Read
   use HDF5, only: H5fis_HDF5_F, H5gclose_F, H5gopen_F
   use HighOutput, only: Dump
   use L1BData, only: L1BData_T, NameLen, PrecisionSuffix, &
     & DeallocateL1BData, ReadL1BData
   use Machine, only: Hp, Getarg
   use MLSCommon, only: R8
   use MLSFiles, only: FileNotFound, &
     & MLS_Exists, MLS_SFStart, MLS_SFEnd, &
     & HDFVersion_5, MLS_HDF_Version, WildCardHDFVersion
   use MLSFillValues, only: IsNaN
   use MLSHDF5, only: MaxNDSNames, DumpHDF5Attributes, DumpHDF5DS, &
     & GetAllHDF5AttrNames, GetAllHDF5DSNames, &
     & MLS_H5Open, MLS_H5Close
   use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning, &
     & MLSMessage
   use MLSStats1, only: FillValueRelation, Stat_T, StatsOnOneLine, &
     & Dump, Statistics
   use MLSStringLists, only: Catlists, GetStringElement, Intersection, &
     & NumStringElements, StringElementNum
   use MLSStrings, only: Indexes, Lowercase, Streq, Trim_Safe
   use Output_M, only: Output, SwitchOutput
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
! dumps datasets or attributes from l2aux files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"
! Then run it
! LF95.Linux/test [options] [input files]

  type options_T
    logical             :: verbose            = .false. ! Print (lots) extra
    logical             :: la                 = .false.
    logical             :: ls                 = .false.
    logical             :: anyNaNs            = .false. ! Just say if any NaNs
    logical             :: printMaxVal        = .false. ! Just print maxval
    logical             :: timereads          = .false. ! Just time how long to read
    logical             :: isL1BOA            = .false.
    logical             :: radiances          = .false.
    logical             :: TAI                = .false.
    logical             :: useFillValue       = .false.
    character(len=16)   :: dumpOptions        = ' '
    character(len=128)  :: DSName      = '' ! Extra dataset if attributes under one
    character(len=128)  :: root        = '/'
    character(len=1024) :: attributes = ''
    character(len=1024) :: datasets   = '*'
    character(len=255)  :: skipList= ''  ! what SDs to skip
    character(len=1)    :: fillValueRelation = '='
    real                :: fillValue  = 0.e0
    integer             :: firstMAF = -1
    integer             :: lastMAF = -1
    integer             :: width = 10
  end type options_T
  
  type ( options_T ) :: options
  integer, parameter ::          MAXDS = MAXNDSNAMES
  integer, parameter ::          MAXSDNAMESBUFSIZE = MAXDS*namelen
  integer, parameter ::          MAXFILES = 4000
  integer, parameter ::          hdfVersion = HDFVERSION_5
  character(len=255) ::          filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer     ::                 i, status, error ! Counting indices & Error flags
  logical     ::                 is_hdf5
  character (len=MAXSDNAMESBUFSIZE) :: mySdList
  integer            ::          n_filenames
  real(r8), dimension(:), pointer :: precisions => null()
  real(r8), dimension(:), pointer :: radiances => null()
  integer     ::                 sdfid1
  real        ::                 t1
  real        ::                 t2
  real        ::                 tFile
  ! 
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  call switchOutput( 'stdout' )
  time_config%use_wall_clock = .true.
  INTPLACES = '8'
  DEFAULTMAXLON = 32
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
  endif
  if ( options%verbose ) call dumpSettings ( options, n_filenames, filenames ) 
  if ( options%useFillValue ) then
    fillValueRelation = options%fillValueRelation
  endif
  call time_now ( t1 )
  do i=1, n_filenames
    call time_now ( tFile )
    if ( options%ls ) then
      print *, 'DS Names in: ', trim(filenames(i)), ' ', trim(options%root)
      call GetAllHDF5DSNames (trim(filenames(i)), trim(options%root), mysdList)
      call dump(mysdList, 'DS names')
    endif 
    if ( options%la ) then
      sdfid1 = mls_sfstart(filenames(i), DFACC_READ, hdfVersion=hdfVersion)
      print *, 'Attribute Names in: ', trim(filenames(i)), ' ', &
        & trim(options%root) // '/' // trim(options%DSName)
      if ( options%DSName /= ' ' ) then
        call GetAllHDF5AttrNames ( sdfid1, mysdList, &
          & DSName=trim(options%root) // '/' // trim(options%DSName) )
      else
        call GetAllHDF5AttrNames ( sdfid1, mysdList, &
          & groupName=trim(options%root) )
      endif
      call dump(mysdList, 'Attribute names')
      status = mls_sfend(sdfid1, hdfVersion=hdfVersion)
    endif
    if ( (options%attributes // options%datasets) == ' ' ) cycle
    if ( options%verbose .or. n_filenames > 1 ) then
      print *, 'Reading from: ', trim(filenames(i))
    endif
    if ( .not. ( options%radiances .or. options%TAI ) ) then
      sdfid1 = mls_sfstart( filenames(i), DFACC_READ, hdfVersion=hdfVersion )
      if ( sdfid1 == -1 ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        &  'Failed to open l2aux file ' // trim(filenames(i)) )
      end if
    end if
    if ( options%datasets /= ' ' ) then
      if ( options%radiances .or.  options%TAI .or. options%timereads &
        & .or. options%anyNaNs .or. options%printMaxVal &
        & .or. options%firstMAF > 0 ) then
        call dumpRadiances ( filenames(i), hdfVersion, options )
        sdfid1 = mls_sfstart( filenames(i), DFACC_READ, hdfVersion=hdfVersion )
      elseif ( options%useFillValue ) then
        call DumpHDF5DS ( sdfid1, trim(options%root), trim(options%datasets), &
          & fillValue=options%fillValue, options=options%dumpOptions )
      else
        call DumpHDF5DS ( sdfid1, trim(options%root), trim(options%datasets), &
          & options=options%dumpOptions )
      endif
    endif
    if ( options%attributes /= ' ' ) then
      if ( options%DSName /= ' ' ) then
        call DumpHDF5Attributes ( sdfid1, trim(options%attributes), &
          & DSName=trim(options%root) // '/' // trim(options%DSName), &
          & options=options%dumpOptions )
      else
        ! print *, 'Trying to dump ', trim(options%attributes), ' from ', trim(options%root)
        call DumpHDF5Attributes ( sdfid1, trim(options%attributes), &
          & groupName=trim(options%root), options=options%dumpOptions )
      endif
    endif
    ! print *, 'About to mls_sfend on ', sdfid1
	 status = mls_sfend(sdfid1, hdfVersion=hdfVersion)
    if ( options%verbose ) call sayTime('reading this file', tFile)
  enddo
  if ( options%verbose .and. .not. ( options%la .or. options%ls ) ) &
    &  call sayTime('reading all files')
  call mls_h5close(error)
contains
!------------------------- dumpSettings ---------------------
    subroutine dumpSettings( options, n_filenames, filenames )
    ! Added for command-line processing
     integer, intent(in)              :: n_filenames
     character(len=255), dimension(:) :: filenames
     type ( options_T ), intent(in)   :: options
     ! Local variables
     integer :: i
     ! print *, 'laconic?            ', options%laconic
     print *, 'verbose?            ', options%verbose
     print *, 'list attributes  ?  ', options%la   
     print *, 'list datasets  ?    ', options%ls
     print *, 'say if any NaNs?    ', options%anyNaNs
     print *, 'print max value?    ', options%printMaxVal
     print *, 'just time reads?    ', options%timereads
     print *, 'radiances only    ? ', options%radiances
     print *, 'TAI only          ? ', options%TAI
     print *, 'useFillValue  ?     ', options%useFillValue
     print *, 'root                ', trim_safe(options%root)
     print *, 'fillValue           ', options%fillValue
     print *, 'fillValueRelation   ', options%fillValueRelation
     print *, 'first maf           ', options%firstMAF
     print *, 'last maf            ', options%lastMAF
     print *, 'width               ', options%width
     print *, 'DSName              ', trim_safe(options%DSName)
     print *, 'attributes          ', trim_safe(options%attributes)
     print *, 'datasets            ', trim_safe(options%datasets)
     print *, 'num files           ', n_filenames
     do i=1, n_filenames
       print *, i, trim(filenames(i))
     enddo
    end subroutine dumpSettings

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
      elseif ( filename(1:3) == '-A ' ) then
        options%attributes = '*'
        options%datasets = ''
        exit
      elseif ( filename(1:3) == '-D ' ) then
        options%datasets = '*'
        exit
      elseif ( filename(1:3) == '-nA ' ) then
        options%attributes = ''
        exit
      elseif ( filename(1:3) == '-nD ' ) then
        options%datasets = ''
        exit
      else if ( filename(1:3) == '-a ' ) then
        call getarg ( i+1+hp, options%attributes )
        i = i + 1
        options%datasets = ''
        exit
      else if ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, options%datasets )
        i = i + 1
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:4) == '-la ' ) then
        options%la = .true.
        options%datasets = ''
        exit
      elseif ( filename(1:4) == '-ls ' ) then
        options%ls = .true.
        options%datasets = ''
        exit
      else if ( filename(1:3) == '-r ' ) then
        call getarg ( i+1+hp, options%root )
        i = i + 1
        exit
      else if ( filename(1:7) == '-first ' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%firstMAF
        i = i + 1
        exit
      else if ( filename(1:6) == '-last ' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%lastMAF
        i = i + 1
        exit
      else if ( filename(1:4) == '-w ' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%width
        defaultWidth = options%width
        i = i + 1
        exit
      elseif ( filename(1:4) == '-one' ) then
        statsOnOneLine = .true.
        options%dumpOptions = trim(options%dumpOptions) // 'v'
        exit
      else if ( filename(1:4) == '-fv ' ) then
        call getarg ( i+1+hp, number )
        read(number, *) options%fillValue
        options%useFillValue = .true.
        i = i + 1
        exit
      else if ( filename(1:5) == '-fvr ' ) then
        call getarg ( i+1+hp, options%fillValueRelation )
        options%useFillValue = .true.
        i = i + 1
        exit
      else if ( filename(1:3) == '-o ' ) then
        call getarg ( i+1+hp, options%dumpOptions )
        if ( index( options%dumpOptions, '?' ) > 0 ) then
          call DumpDumpOptions( "?" )
          stop
        endif
        i = i + 1
        exit
      else if ( filename(1:5) == '-radi' ) then
        options%radiances = .true.
        exit
      else if ( filename(1:4) == '-TAI' ) then
        options%TAI = .true.
        exit
      else if ( filename(1:3) == '-rd ' ) then
        call getarg ( i+1+hp, options%DSName )
        i = i + 1
        exit
      else if ( filename(1:3) == '-t ' ) then
        options%timereads = .true.
        exit
      else if ( lowercase(filename(1:4)) == '-nan' ) then
        options%anyNaNs = .true.
        exit
      else if ( lowercase(filename(1:4)) == '-max' ) then
        options%printMaxVal = .true.
        options%isL1BOA = .true. ! So we don't try to read precisions
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
      print *,  "Enter the name of the HDF5 l2aux file. "
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Usage:l2auxdump [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options:'
      write (*,*) '  -f filename     => add filename to list of filenames'
      write (*,*) '                   (can do the same w/o the -f)'
      write (*,*) '  -v              => switch on verbose mode'
      write (*,*) '  -w width        => use width when dumping char arrays'
      write (*,*) '  -la             => just list attribute names in files'
      write (*,*) '  -ls             => just list sd names in files'
      write (*,*) '  -skip list      => skip dumping the SDs in list'
      write (*,*) '  -max            => print just max value'
      write (*,*) '  -NaN            => just say if there are any NaNs'
      write (*,*) '  -t              => just time reads'
      write (*,*) '  -radiances      => show radiances only'
      write (*,*) '  -TAI            => show TAI times only'
      write (*,*) '  -o opts         => pass opts to dump routines'
      write (*,*) '                  e.g., "-rs" to dump only rms, stats'
      write (*,*) '                  e.g., "?" to list available ones'
      write (*,*) '  -one            => print statistics on one line (dont)'
      write (*,*) '  -r root         => limit to group based at root'
      write (*,*) '                     (default is "/")'
      write (*,*) '  -rd DSName      => limit attributes to root/DSName'
      write (*,*) '                     (default is group attributes at root)'
      write (*,*) '  -first maf1     => read l1b starting with maf1'
      write (*,*) '  -last maf1      => read l1b ending with maf1'
      write (*,*) '  -fv value       => filter rms, % around value'
      write (*,*) '  -fvr relation   => one of {"=","<",">"}'
      write (*,*) '                      we filter values standing in'
      write (*,*) '                      this relation with fillValue'
      write (*,*) '  -A              => dump all attributes'
      write (*,*) '  -D              => dump all datasets (default)'
      write (*,*) '  -nA             => do not dump attributes (default)'
      write (*,*) '  -nD             => do not dump datasets'
      write (*,*) '  -a a1,a2,..     => dump just attributes named a1,a2,..'
      write (*,*) '  -d d1,d2,..     => dump just datasets named a1,a2,..'
      write (*,*) '                 (may use \* as wild card)'
      write (*,*) '  -h              => print brief help'
      stop
  end subroutine print_help

  ! ---------------------- dumpradiances  ---------------------------
  subroutine dumpradiances( file1, hdfVersion, options )
  !------------------------------------------------------------------------

    ! Given file names file1,
    ! This routine prints the radiances based on the pattern
    ! that if a DS named 'x' is a radiance, then
    ! the file must also contain 'x precision'

    ! Arguments

    character (len=*), intent(in)     :: file1 ! Name of file
    integer, intent(in)               :: hdfVersion
    type ( options_T )                :: options

    ! Local
    logical, parameter                :: countEmpty = .true.
    logical                           :: file_exists
    integer                           :: grpid
    integer                           :: i
    integer                           :: iPrec
    logical                           :: isl1boa
    type(l1bdata_t)                   :: L1BPrecision  ! Result
    type(l1bdata_t)                   :: L1BRadiance   ! Result
    type(Stat_T)                      :: L1BStat
    character (len=MAXSDNAMESBUFSIZE) :: matches
    character (len=MAXSDNAMESBUFSIZE) :: mySdList
    integer                           :: NoMAFs
    integer                           :: noSds
    integer                           :: sdfid1
    character (len=80)                :: sdName
    integer, dimension(3)             :: shp
    integer                           :: status
    integer                           :: the_hdfVersion
    character (len=80)                :: which
    
    ! Executable code
    if ( .not. any( indexes(options%dumpOptions, (/ 'r', 's' /) ) > 0 ) )  then
      which = '*'
    else
      which = ' '
      if ( index(options%dumpOptions, 'r') > 0 ) which = 'rms'
      if ( index(options%dumpOptions, 's') > 0 ) which = catLists( which, 'max,min,mean,stddev' )
    endif 
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
    if ( options%verbose ) then
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

    options%isl1boa = options%isl1boa .or. (index(trim(mysdList), 'GHz/') > 0)
    if ( options%isl1boa .and. .not. &
      & ( &
      & options%timereads .or. options%anyNaNs .or. options%printMaxVal .or. &
      & options%TAI .or. options%firstMAF > 0 ) &
      & ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'l1boa file contains no radiances ' // trim(File1) )
      return
    endif
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
    ! print *, 'Looping over ', noSds
    do i = 1, noSds
      call GetStringElement (trim(mysdList), sdName, i, countEmpty )
      ! print *, trim(sdName)
      if ( index( lowercase(trim(sdName)), PRECISIONSUFFIX ) > 0 ) cycle
      if ( options%TAI .and. index( lowercase(trim(sdName)), 'tai' ) < 1 ) cycle
      iPrec = StringElementNum( mysdList, trim(sdName) // PRECISIONSUFFIX, &
        & countEmpty )
      ! print *, iPrec, options%isl1boa
      if ( iPrec < 1 .and. .not. options%isl1boa ) cycle
      ! We won't try to dump metadata or obvious non-radiance data
      if ( any( &
        & streq( &
        & (/ 'PCF ', 'meta', 'l2cf', 'utcp', 'leap', 'LCF ' /), &
        & sdname, options='-Pw' ) ) .or. &
        &  index(options%skipList, trim(sdName)) > 0 ) cycle
      ! Did we ask for any datasets by name? With wildcard?
      if ( len_trim(options%datasets) > 1 .and. &
        & index( options%datasets, '*' ) > 0 ) then
        matches = Intersection( trim(sdName), trim(options%datasets), options='-w' )
        if ( len_trim(matches) < 1 ) cycle
      elseif ( options%datasets /= '*' ) then
        status = stringElementNum(options%datasets, trim(sdName), countEmpty)
        if ( status < 1 ) cycle
      endif
      ! Allocate and fill l2aux
      if ( options%verbose ) print *, 'About to read ', trim(sdName)
      if ( options%firstMAF > -1 ) then
        call ReadL1BData ( sdfid1, trim(sdName), L1bRadiance, &
          & NoMAFs, status, firstMAF=options%firstMAF, lastMAF=options%lastMAF, &
          & hdfVersion=the_hdfVersion, NEVERFAIL=.true., L2AUX=.true. )
      else
        call ReadL1BData ( sdfid1, trim(sdName), L1bRadiance, NoMAFs, status, &
          & hdfVersion=the_hdfVersion, NEVERFAIL=.true., L2AUX=.true. )
      endif
      if ( status /= 0 ) then
	     call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to find ' // trim(sdName) // ' in ' // trim(File1) )
        cycle
      endif
      if ( options%isl1boa ) then
        status = 0 ! No radiance precisions in l1boa
      else
        if ( options%verbose ) print *, 'About to read ', trim(sdName) // PRECISIONSUFFIX
        if ( options%firstMAF > -1 ) then
          call ReadL1BData ( sdfid1, trim(sdName)  // PRECISIONSUFFIX, L1bPrecision, &
            & NoMAFs, status, firstMAF=options%firstMAF, lastMAF=options%lastMAF, &
            & hdfVersion=the_hdfVersion, NEVERFAIL=.true., L2AUX=.true. )
        else
          call ReadL1BData ( sdfid1, trim(sdName)  // PRECISIONSUFFIX, L1bPrecision, &
            & NoMAFs, status, &
            & hdfVersion=the_hdfVersion, NEVERFAIL=.true., L2AUX=.true. )
        endif
      endif
      if ( status /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Unable to find ' // trim(sdName)  // PRECISIONSUFFIX // &
          & ' in ' // trim(File1) )
        cycle
      elseif ( .not. associated(L1bRadiance%DpField) ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'no f.p. values of ' // trim(sdName)  // &
          & ' in ' // trim(File1) )
        cycle
      endif
      if ( options%timereads ) cycle
      if ( options%anyNaNs ) then
        if ( .not. associated(L1bRadiance%DpField) ) then
          print *, trim(sdName), ' not d.p. and so not checked'
        elseif ( any(isNan(L1bRadiance%DpField)) ) then
          print *, 'NaNs found in ', trim(sdName)
        elseif ( options%verbose ) then
          print *, 'no NaNs found in ', trim(sdName)
        endif
        cycle
      elseif ( options%printMaxVal ) then
        if ( .not. associated(L1bRadiance%DpField) ) then
          print *, trim(sdName), ' not d.p. and so not checked'
        elseif ( options%verbose ) then
          print *, 'max(' // trim(sdName) // ')', &
            & maxval(L1bRadiance%DpField)
        else
          print *, maxval(L1bRadiance%DpField)
        endif
        cycle
      endif
      ! Convert tai93s to hours-in-day
      if ( options%TAI ) L1bRadiance%DpField = tai93s2hid( L1bRadiance%DpField )
      shp = shape(L1bRadiance%DpField)
      if ( any( indexes(options%dumpOptions, (/ 'r', 's' /) ) > 0 ) )  then
      elseif ( options%useFillValue ) then
        call dump( L1bRadiance%DpField, name=trim(sdName), &
          & options='s', FillValue=real(options%fillValue, r8) )
      else
        call dump( L1bRadiance%DpField, name=trim(sdName), &
          & options=options%dumpOptions )
        if ( associated(L1bPrecision%DpField) ) &
          &  call dump( L1bPrecision%DpField, name=trim(sdName) // precisionSuffix, &
          & options=options%dumpOptions )
      endif
      L1BStat%count = 0
      if ( options%verbose ) print *, 'About to reshape radiances'
      call allocate_test ( radiances, product(shp), 'radiances', ModuleName )
      radiances = reshape( L1bRadiance%DpField, (/ product(shp) /) )
      call DeallocateL1BData ( l1bRadiance )
      if ( associated(L1bPrecision%DpField) ) then
        if ( options%verbose ) print *, 'About to reshape precisions'
        call allocate_test ( precisions, product(shp), 'precisions', ModuleName )
        precisions = reshape( L1bPrecision%DpField, (/ product(shp) /) )
        call DeallocateL1BData ( l1bPrecision )
        if ( options%verbose ) print *, 'About to call statistics'
        call statistics(&
          & radiances, &
          & L1BStat, &
          & precision=precisions )
        call dump ( L1BStat, which )
      endif
      call deallocate_test ( radiances, 'radiances', ModuleName )
      call deallocate_test ( precisions, 'precisions', ModuleName )
    enddo
	 call h5gClose_f ( grpID, status )
    if ( status /= 0 ) then
	   call MLSMessage ( MLSMSG_Warning, ModuleName, &
       & 'Unable to close group in l2aux file: ' // trim(File1) // ' after diffing' )
    endif
	 status = mls_sfend( sdfid1, hdfVersion=the_hdfVersion )
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2aux file: " // trim(File1) // ' after diffing' )
  end subroutine dumpradiances

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

!==================
end program l2auxdump
!==================

! $Log$
! Revision 1.24  2017/12/01 00:28:06  pwagner
! Now obeys -maf1, -maf2 cmdline options
!
! Revision 1.23  2016/08/09 22:45:26  pwagner
! Consistent with splitting of Dunp_0
!
! Revision 1.22  2016/04/05 23:57:38  pwagner
! -one cmdline option added; prints name on each line
!
! Revision 1.21  2016/01/22 00:38:04  pwagner
! May use wildcard as part of dataset names
!
! Revision 1.20  2016/01/12 00:52:12  pwagner
! May override DEFAULTWIDTH when dumping char array
!
! Revision 1.19  2015/04/17 23:03:03  pwagner
! Added -TAI option to convert tai93s to hours-in-day
!
! Revision 1.18  2015/01/24 00:06:03  pwagner
! Added commandline option to detect NaNs lin l1b files
!
! Revision 1.17  2014/03/07 21:47:21  pwagner
! Name_Len changed to nameLen; got from MLSCommon
!
! Revision 1.16  2014/01/09 00:31:26  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 1.15  2013/08/23 02:51:47  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.14  2012/06/14 00:02:44  pwagner
! Changed how dump options input; now use -o 'opts'
!
! Revision 1.13  2011/01/04 00:50:45  pwagner
! Shortened field width of MAFStartTimeUTC dumps to 32
!
! Revision 1.12  2010/07/23 17:53:27  pwagner
! Now able to limit dump to a range of mafs
!
! Revision 1.11  2010/03/11 23:35:55  pwagner
! verbose output looks better
!
! Revision 1.10  2009/11/20 23:03:38  pwagner
! -shape option just dumps array rank, shape
!
! Revision 1.9  2009/08/18 20:42:15  pwagner
! Dumping counterMAF array looks better
!
! Revision 1.8  2009/06/16 22:38:28  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 1.7  2009/05/08 17:33:09  pwagner
! New command line option -lac to prevent printing non-data
!
! Revision 1.6  2008/09/09 16:53:34  pwagner
! Many changes, hoping some are correct
!
! Revision 1.5  2007/08/17 00:42:38  pwagner
! Needed to increase MAXDS
!
! Revision 1.4  2007/02/06 23:20:19  pwagner
! -radiances options dumps only radiances
!
! Revision 1.3  2006/07/13 18:11:01  pwagner
! May set fillValue and -relation for rms, pctage
!
! Revision 1.2  2006/06/29 20:39:45  pwagner
! Repaired a few bugs
!
! Revision 1.1  2006/06/28 00:06:54  pwagner
! First commit
!
