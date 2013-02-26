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
PROGRAM L2GPDump ! dumps L2GPData files
!=================================

   use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
   use BITSTUFF, only: ISBITSET
   use DUMP_0, only: DUMP
   use HDF, only: DFACC_READ
   use HDF5, only: H5FCLOSE_F, H5GOPEN_F, H5GCLOSE_F, H5FIS_HDF5_F   
   use INTRINSIC, only: L_SWATH
   use L2GPDATA, only: L2GPDATA_T, L2GPNAMELEN, MAXSWATHNAMESBUFSIZE, RGP, &
     & CONTRACTL2GPRECORD, DUMP, DUMPRANGE, READL2GPDATA, DESTROYL2GPCONTENTS, &
     & SETUPNEWL2GPRECORD
   use MACHINE, only: HP, GETARG
   use MLSCOMMON, only: MLSFILE_T
   use MLSFILES, only: HDFVERSION_5, INITIALIZEMLSFILE, MLS_INQSWATH, &
     & CLOSE_MLSFILE, OPEN_MLSFILE, SPLIT_PATH_NAME
   use MLSHDF5, only: MLS_H5OPEN, MLS_H5CLOSE
   use MLSMESSAGEMODULE, only: MLSMESSAGECONFIG, MLSMSG_ERROR, MLSMSG_WARNING, &
     & MLSMESSAGE
   use MLSSTRINGLISTS, only: CATLISTS, EXPANDSTRINGRANGE, &
     & GETSTRINGELEMENT, NUMSTRINGELEMENTS, READINTSFROMLIST, &
     & STRINGELEMENT, STRINGELEMENTNUM
   use MLSSTRINGS, only: READNUMSFROMCHARS
   use OUTPUT_M, only: BLANKS, NEWLINE, OUTPUT, OUTPUTNAMEDVALUE, &
     & RESUMEOUTPUT, SUSPENDOUTPUT
   
   implicit none

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program dumps L2GPData files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [filenames]

  integer, parameter ::  max_nsds = 1000  ! Maximum number of datasets in file.
  integer, parameter :: MAXNCHUNKS = 50

  type options_T
     character(len=255) ::  chunks = '*' ! wild card means 'all'
     integer, dimension(2)::chunkRange = 0
     real, dimension(2)  :: hoursRange = 0.
     logical ::             debug = .true.
     logical ::             verbose = .false.
     logical ::             senseInquiry = .true.
     integer ::             details = 1
     integer ::             width = 5
     logical ::             columnsOnly = .false.
     logical ::             attributesToo = .false.
     character(len=16)  ::  dumpOptions = ''
     character(len=255) ::  dsInquiry = ''
     character(len=255) ::  attrInquiry = ''
     character(len=255) ::  fields = ''
     character(len=255) ::  swaths = '*' ! wildcard, meaning all swaths
     character(len=255) ::  geoBoxNames = '' ! which geolocation names to box
     integer            ::  nGeoBoxDims = 0
     real, dimension(4) ::  geoBoxLowBound
     real, dimension(4) ::  geoBoxHiBound
     real    ::             ConvergenceCutOff = -1. ! Show % above, below this
     real    ::             PrecisionCutOff = -1. ! Show % above, below this
     real    ::             QualityCutOff = -1. ! Show % above, below this
     logical ::             StatusBits = .false. ! Show % with various status bits set
     logical ::             merge = .false. ! Show % after merging input files
  end type options_T

  type ( options_T ) :: options
  character(LEN=255) :: filename          ! filename
  integer            :: n_filenames
  integer     ::  error ! Counting indices & Error flags
  logical     :: is_hdf5
  logical     :: is_present
  integer, save                   :: numGood = 0
  integer, save                   :: numGoodPrec = 0
  integer, save                   :: numNotUseable = 0
  integer, save                   :: numOddStatus = 0
  real, dimension(3), save        :: numTest = 0.
  integer, parameter              :: MAXNUMBITSUSED = 9
  integer, dimension(MAXNUMBITSUSED), parameter :: bitNumber = &
    & (/ 0, 1, 2, 4, 5, 6, 7, 8, 9 /)
  integer, dimension(MAXNUMBITSUSED, 2), save :: bitCounts = 0
  character(len=*), parameter     :: bitNames = &
    & '  dontuse,   bewary,     info,    hicld,    locld,   nogmao,abandoned,   toofew,    crash'
  ! 
  MLSMessageConfig%useToolkit = .false.   
  MLSMessageConfig%logFileUnit = -1       
  CALL mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename(filename, n_filenames, options)
     if ( filename == ' ' ) exit
     n_filenames = n_filenames + 1
     call h5fis_hdf5_f(trim(filename), is_hdf5, error)
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(filename)
     endif
     if ( options%dsInquiry /= ' ' ) then
       call suspendOutput
       is_present = IsDSInFile( trim(filename), trim(options%dsInquiry) )
       call Respond( options%senseInquiry, is_present, &
         & trim(filename), trim(options%dsInquiry) )
     elseif ( options%attrInquiry /= ' ' ) then
       call suspendOutput
       is_present = IsAttributeInFile( trim(filename), trim(options%attrInquiry) )
       call Respond( options%senseInquiry, is_present, &
         & trim(filename), trim(options%attrInquiry) )
     else
       if ( options%verbose ) print *, 'Dumping swaths in ', trim(filename)
       call dump_one_file(trim(filename), options)
     endif
     call resumeOutput
  enddo
  if ( options%merge ) then
    if ( options%ConvergenceCutOff > -1. ) &
      & call showPercentages( numTest(1), numGood, 'convergence', &
      & options%ConvergenceCutOff, 'above' )
    if ( options%QualityCutOff > -1. ) &
      & call showPercentages( numTest(2), numGood, 'Quality', &
      & options%QualityCutOff, 'below' )
    if ( options%PrecisionCutOff > -1. ) &
      & call showPercentages( numTest(3), numGoodPrec, 'precision', &
      & options%PrecisionCutOff, 'above' )
    if ( options%statusBits ) call showStatusPct
    call ShowSummary
  endif
  call mls_h5close(error)
contains
!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(len=255), intent(out) :: filename          ! filename
     integer, intent(in) ::             n_filenames
     type ( options_T ) :: options
     integer ::                         error = 1
     integer, save ::                   i = 1
     character(len=255) :: argstr
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
      elseif ( filename(1:3) == '-0 ' ) then
        options%details = 0
      elseif ( filename(1:3) == '-1 ' ) then
        options%details = -1
      elseif ( filename(1:3) == '-2 ' ) then
        options%details = -2
      elseif ( filename(1:3) == '-a ' ) then
        options%attributesToo = .true.
      elseif ( filename(1:3) == '-c ' ) then
        options%columnsOnly = .true.
      else if ( filename(1:6) == '-chunk' ) then
        call getarg ( i+1+hp, options%chunks )
        i = i + 1
        if ( index(options%chunks, '-') > 0 ) then
          ! call ReadIntsFromList( options%chunks, options%chunks, &
          !  & inseparator='-' )
          print *, 'Sorry--ReadIntsFromList not coded for this format'
          stop
        else
          call ReadIntsFromList( options%chunks, options%chunkRange )
          if ( options%chunkRange(2) < options%chunkRange(1) ) &
            & options%chunkRange(2) = options%chunkRange(1)
        endif
      else if ( filename(1:5) == '-hour' ) then
        call getarg ( i+1+hp, argstr )
        call ReadNumsFromChars( StringElement (argstr, 1, .true. ), &
          & options%hoursRange(1) )
        call ReadNumsFromChars( StringElement (argstr, 2, .true. ), &
          & options%hoursRange(2) )
        if ( options%hoursRange(2) < options%hoursRange(1) ) &
          & options%hoursRange(2) = options%hoursRange(1)
        print *, 'argstr: ', argstr
        print *, 'hours: ', options%hoursRange
        i = i + 1
      else if ( filename(1:5) == '-conv' ) then
        call getarg ( i+1+hp, argstr )
        read( argstr, * ) options%convergenceCutOff
        i = i + 1
      else if ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, options%dumpOptions )
        i = i + 1
      else if ( filename(1:4) == '-geo' ) then
        call getarg ( i+1+hp, filename )
        options%geoBoxNames = catLists( options%geoBoxNames, filename )
        i = i + 1
        options%nGeoBoxDims = min( options%nGeoBoxDims + 1, 4 )
        call getarg ( i+1+hp, filename )
        read( filename, * ) options%geoBoxLowBound(options%nGeoBoxDims), &
          & options%geoBoxHiBound(options%nGeoBoxDims)
        i = i + 1
      else if ( filename(1:6) == '-inqat' ) then
        call getarg ( i+1+hp, options%attrInquiry )
        i = i + 1
      elseif ( filename(1:2) == '-m' ) then
        options%merge = .true.
      else if ( filename(1:6) == '-inqds' ) then
        call getarg ( i+1+hp, options%dsInquiry )
        i = i + 1
      else if ( filename(1:7) == '-ninqat' ) then
        options%senseInquiry = .false.
        call getarg ( i+1+hp, options%attrInquiry )
        i = i + 1
      else if ( filename(1:7) == '-ninqds' ) then
        options%senseInquiry = .false.
        call getarg ( i+1+hp, options%dsInquiry )
        i = i + 1
      else if ( filename(1:3) == '-l ' ) then
        call getarg ( i+1+hp, options%fields )
        i = i + 1
      else if ( filename(1:5) == '-prec' ) then
        call getarg ( i+1+hp, argstr )
        read( argstr, * ) options%precisionCutOff
        i = i + 1
      else if ( filename(1:2) == '-w' ) then
        call getarg ( i+1+hp, argstr )
        read( argstr, * ) options%width
        i = i + 1
      else if ( filename(1:5) == '-qual' ) then
        call getarg ( i+1+hp, argstr )
        read( argstr, * ) options%qualityCutOff
        i = i + 1
      else if ( filename(1:3) == '-s ' ) then
        call getarg ( i+1+hp, options%swaths )
        i = i + 1
      elseif ( filename(1:5) == '-stat' ) then
        options%statusBits = .true.
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
      & 'Usage: l2gpdump [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) &
      & ' optionally restrict dumps to certain fields, chunks, etc.'
      write (*,*) ' Options: -f filename => use filename'
      write (*,*) '          -h          => print brief help'
      write (*,*) '          -chunks cl  => dump only chunks named in cl'
      write (*,*) '          -geo name lo,hi  '
      write (*,*) '                      => dump only geobox low <= geo <= hi'
      write (*,*) '                      where geo in {latitude, longitude, time, pressure}'
      write (*,*) '                      (may be repeated)'
      write (*,*) '                      if hi < lo then dump is outside geobox'
      write (*,*) '          -d opts     => pass opts to dump routines'
      write (*,*) '                          e.g., "-rs" to dump only rms, stats'
      write (*,*) '          -[n]inqattr attr'
      write (*,*) '                      => print only if attribute attr [not] present'
      write (*,*) '          -[n]inqds ds'
      write (*,*) '                      => print only if dataset ds [not] present'
      write (*,*) '          -l list     => dump only fields named in list'
      write (*,*) '          -s slist    => dump only swaths named in slist'
      write (*,*) '          -c          => dump only column abundances'
      write (*,*) '          -a          => dump attributes, too'
      write (*,*) '          -v          => verbose'
      write (*,*) '          -w width    => dump with width items on each line'
      write (*,*) '     (details level)'
      write (*,*) '          -0          => dump only scalars, 1-d array'
      write (*,*) '          -1          => dump only scalars'
      write (*,*) '          -2          => dump only swath names'

      write (*,*) '    (The following options print only summaries)'
      write (*,*) '          -conv x     => show % nonconverged by x cutoff'
      write (*,*) '          -prec x     => show % with precision > x cutoff'
      write (*,*) '          -qual x     => show % with quality < x cutoff'
      write (*,*) '          -status     => show % with various status bits set'
      write (*,*) '          -m[erge]    => merge data from all files (dont)'

      write (*,*) '    (Notes)'
      write (*,*) ' (1) by default, dumps all fields in allswaths, but not attributes'
      write (*,*) ' (2) by default, detail level is -1'
      write (*,*) ' (3) details levels, -l options are all mutually exclusive'
      write (*,*) ' (4) the list of chunks may include the range operator "-"'
      write (*,*) ' (5) -conv, -qual, -prec, and -status all turn off detailed dumps'
      stop
  end subroutine print_help
  
  function IsAttributeInFile( file, attribute ) result(sooDesu)
    use MLSHDF5, only: IsHDF5ItemPresent
    use HDF5, only: h5fopen_f, H5F_ACC_RDONLY_F
    ! Dummy args
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: attribute
    logical :: sooDesu
    ! Local variables
    integer :: fileID
    integer :: grpID
    integer :: status
    character(len=len(attribute)) :: path, name
    ! TRUE if attribute in file
    call h5fopen_f ( trim(file), H5F_ACC_RDONLY_F, fileID, status )
    if ( status /= 0 ) call defeat('Opening file')
    call split_path_name ( attribute, path, name )
    call h5gopen_f( fileID, trim(path), grpID, status )
    if ( status /= 0 ) call defeat('Opening group')
    sooDesu = IsHDF5ItemPresent ( grpID, name, '-a' )
    call h5gclose_f(grpID, status)
    if ( status /= 0 ) call defeat('Closing group')
    call h5fclose_f(fileID, status)
    if ( status /= 0 ) call defeat('Closing file')
  end function IsAttributeInFile

  function IsDSInFile( file, DS ) result(sooDesu)
    use MLSHDF5, only: IsHDF5ItemPresent
    use HDF5, only: h5fopen_f, H5F_ACC_RDONLY_F
    ! Dummy args
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: DS
    logical :: sooDesu
    ! Local variables
    integer :: fileID
    integer :: status
    integer :: grpID
    character(len=len(DS)) :: path, name
    ! TRUE if DS in file
    call h5fopen_f ( trim(file), H5F_ACC_RDONLY_F, fileID, status )
    call split_path_name ( DS, path, name )
    call h5gopen_f( fileID, trim(path), grpID, status )
    sooDesu = IsHDF5ItemPresent ( grpID, name, '-d' )
    call h5gclose_f(grpID, status)
    call h5fclose_f(fileID, status)
  end function IsDSInFile
  
  subroutine Defeat(msg)
    character(len=*), intent(in) :: msg
    call resumeOutput
    call output('Serious error: ' // msg, advance='yes')
    call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'stopping' )
  end subroutine Defeat
  
  subroutine Respond( sense, test, file, name )
    ! print only if sense matches test
    ! Dummy args
    logical, intent(in)          :: sense
    logical, intent(in)          :: test
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: name
    character(len=*), parameter  :: Found = 'found'
    character(len=*), parameter  :: notFound = 'not found'
    character(len=16)            :: answer
    ! Executable
    if ( sense .neqv. test ) return
    if ( sense ) then
      answer = found
    else
      answer = notfound
    endif
    call resumeOutput
    call output(trim(name) // ' ' // trim(answer) &
      & // ' in ' // trim(file), advance='yes' )
  end subroutine Respond

   subroutine dump_one_file(filename, options)
    ! Dummy args
    character(len=*), intent(in) :: filename          ! filename
    type ( options_T ) :: options
    ! Local variables
    logical, parameter            :: countEmpty = .true.
    integer :: File1
    integer :: i
    integer :: listsize
    type (L2GPData_T) :: l2gp
    type(MLSFile_T)                :: MLSFile
    integer :: noSwaths
    integer :: record_length
    integer :: status
    character (len=L2GPNameLen) :: swath
    character (len=MAXSWATHNAMESBUFSIZE) :: SwathList
    ! Get swath list
    noSwaths = mls_InqSwath ( filename, SwathList, listSize, &
           & hdfVersion=HDFVERSION_5)
    if ( options%details < -1 ) then
      call output('swaths in ' // trim(filename), advance='yes')
      call dump(trim(swathList))
      return
    endif
    ! Executable code
    noSwaths = NumStringElements(trim(swathList), countEmpty)
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths to dump--unable to count swaths in ' // trim(swathList) )
    endif
    status = InitializeMLSFile ( MLSFile, type=l_swath, access=DFACC_READ, &
     & name=filename, HDFVersion=HDFVERSION_5 )
    ! Loop over swaths in file 1
    do i = 1, noSwaths
      call GetStringElement (trim(swathList), swath, i, countEmpty )
      ! Is this one of the swaths we wished to dump?
      if ( options%swaths /= '*' ) then
        status = stringElementNum(options%swaths, trim(swath), countEmpty)
        if ( status < 1 ) cycle
      endif
      ! Allocate and fill l2gp
      ! print *, 'Reading swath from file: ', trim(swath)
      call ReadL2GPData ( trim(filename), trim(swath), l2gp, &
           & hdfVersion=HDFVERSION_5 )
      ! Dump the swath- and file-level attributes
      if ( options%attributesToo ) then
        call open_MLSFile( MLSFile )
        file1 = MLSFile%FileID%f_id
        call dump(file1, l2gp)
        call close_MLSFile ( MLSFile )
      endif
      if ( options%ConvergenceCutOff > -1. .or. options%QualityCutOff > -1. .or. &
        & options%PrecisionCutOff > -1. .or. options%statusBits ) then
        call myPcts( options, l2gp, swath )
      else
        call myDump( options, l2gp, swath )
      endif
      call DestroyL2GPContents ( l2gp )
    enddo
   end subroutine dump_one_file

   subroutine myPcts( options, inl2gp, swath )
     ! Args
     type ( options_T ), intent(in)  :: options
     type (L2GPData_T), intent(in)   :: inl2gp
     character(len=*), intent(in)    :: swath
     ! Internal variables
     integer                         :: bitindex, i, n
     type (L2GPData_T)               :: l2gp
     logical, dimension(:), pointer  :: negativePrec => null() ! true if all prec < 0
     logical, dimension(:), pointer  :: oddStatus => null() ! true if all status odd
     ! Executable
     if ( options%verbose ) print *, 'swath: ', trim(swath)
     ! Contract to just the levels and instances requested
     n = options%nGeoBoxDims
     if ( any(options%chunkRange /= 0) ) then
       call ContractL2GPRecord ( inl2gp, l2gp, chunks=options%chunkRange )
     elseif ( any(options%hoursRange /= 0.) ) then
       print *, 'Contracting according to hours in day: ', options%hoursRange
       call ContractL2GPRecord ( inl2gp, l2gp, hoursInDay=options%hoursRange )
     elseif ( n > 0 ) then
       call ContractL2GPRecord ( inl2gp, l2gp, &
         & options%geoBoxNames, &
         & options%geoBoxLowBound(1:n), options%geoBoxHiBound(1:n) )
     else
       call SetupNewL2GPRecord ( l2gp, proto=inl2gp )
     endif
     ! print *, 'chunkRange ', options%chunkRange
     ! print *, 'Contracted l2gp'
     ! call dump( l2gp )
     if ( .not. options%merge ) then
       bitCounts = 0
       numGood = 0
       numGoodPrec = 0
       numNotUseable = 0
       numOddStatus = 0
       numTest = 0.
     endif
     call allocate_test( negativePrec, l2gp%nTimes, 'negativePrec', ModuleName )
     call allocate_test( oddStatus, l2gp%nTimes, 'oddStatus', ModuleName )
     do i=1, l2gp%nTimes
       negativePrec(i) = all( l2gp%l2GPPrecision(:,:,i) < 0._rgp )
     enddo
     do i=1, l2gp%nTimes
       oddStatus(i) = mod(l2gp%status(i), 2) > 0
     enddo
     numGood = numGood + count( .not. ( negativePrec .or. &
       & (mod(l2gp%status, 2) > 0) ) )
     if ( options%ConvergenceCutOff > -1. ) then
       numTest(1) = numTest(1) + count( .not. ( negativePrec .or. &
         & (mod(l2gp%status, 2) > 0) .or. &
         & (l2gp%Convergence < options%ConvergenceCutOff) ) )
       if ( .not. options%merge ) &
         & call showPercentages( numTest(1), numGood, 'convergence', &
         & options%ConvergenceCutOff, 'above' )
     endif
     if ( options%QualityCutOff > -1. ) then
       numTest(2) = numtest(2) + count( .not. ( negativePrec .or. &
         & (mod(l2gp%status, 2) > 0) .or. &
         & (l2gp%Quality > options%QualityCutOff) ) )
       print *, 'numTest(2) ', numTest(2)
       print *, 'numGood ', numGood
       if ( .not. options%merge ) &
         & call showPercentages( numTest(2), numGood, 'Quality', &
         & options%QualityCutOff, 'below' )
     endif
     if ( options%PrecisionCutOff > -1. ) then
       ! numGood = 0
       ! numTest = 0
       do i=1, l2gp%nTimes
         if ( negativePrec(i) .or. &
           & mod(l2gp%status(i), 2) > 0 ) cycle
         numGoodPrec = numGoodPrec + &
           & l2gp%nLevels*max(1, l2gp%nFreqs)
         numTest(3) = numTest(3) + &
           & count( l2gp%l2gpPrecision(:,:,i) > options%PrecisionCutOff )
       enddo
       if ( .not. options%merge ) &
         & call showPercentages( numTest(3), numGoodPrec, 'precision', &
         & options%PrecisionCutOff, 'above' )
     endif
     ! numGood = count( .not. ( negativePrec .or. &
     !   & (mod(l2gp%status, 2) > 0) ) )
     numNotUseable = numNotUseable + count ( negativePrec .or. oddStatus )
     numOddStatus = numOddStatus + count( oddStatus )
     if ( options%statusBits ) then
       ! First, and last 3 bits are special
       ! For bit 0 we filter out only points with precision < 0
       ! For the last 3
       ! we want % of crashed chunks, so we don't filter out at all
       bitCounts(1, 2) = bitCounts(1, 2) + count( .not. negativePrec )
       bitCounts(1, 1) = bitCounts(1, 1) + count( .not. ( negativePrec .or. &
         & (mod(l2gp%status, 2) == 0) ) )
       do bitindex=MAXNUMBITSUSED-2, MAXNUMBITSUSED
         bitCounts(bitindex, 2) = bitCounts(bitindex, 2) + l2gp%nTimes
         bitCounts(bitindex, 1) = bitCounts(bitindex, 1) + &
           & count(isBitSet( l2gp%status, bitNumber(bitindex) ) )
       enddo
       ! call outputNamedValue ( 'max status', maxval(l2gp%status) )
       ! call outputNamedValue ( 'min status', minval(l2gp%status) )
       do bitindex=2, MAXNUMBITSUSED-3
         ! Bit 
         bitCounts(bitindex, 2) = numGood
         bitCounts(bitindex, 1) = bitCounts(bitindex, 1) + &
           & count( .not. ( negativePrec .or. &
           & (mod(l2gp%status, 2) > 0) ) .and. &
           & isBitSet( l2gp%status, bitNumber(bitindex) ) )
       enddo
       if ( .not. options%merge ) then
         call showStatusPct
         call showSummary
       endif
     endif
     call deallocate_test( negativePrec, 'negativePrec', ModuleName )
     call deallocate_test( oddStatus, 'oddStatus', ModuleName )
     call DestroyL2GPContents( l2gp )
   end subroutine myPcts

   subroutine myDump( options, l2gp, swath )
     ! Args
     type ( options_T ), intent(in)  :: options
     type (L2GPData_T), intent(in)   :: l2gp
     character(len=*), intent(in)    :: swath
     ! Internal variables
     integer, dimension(MAXNCHUNKS)  :: chunks
     type (L2GPData_T)               :: contractedl2gp
     integer                         :: nChunks
     ! Dump the actual swath
     if ( options%verbose ) print *, 'swath: ', trim(swath)
     if ( options%nGeoBoxDims > 0 ) then
       call dumpRange( l2gp, &
         & options%geoBoxNames, options%geoBoxLowBound, options%geoBoxHiBound, &
         & columnsOnly=options%columnsOnly, details=options%details, &
         & fields=options%fields, options=options%dumpOptions )
     elseif ( any(options%hoursRange /= 0.) ) then
       print *, 'Contracting according to hours in day: ', options%hoursRange
       call ContractL2GPRecord ( l2gp, contractedl2gp, &
         & hoursInDay=options%hoursRange )
       call dump( contractedl2gp, &
         & columnsOnly=options%columnsOnly, details=options%details, &
         & fields=options%fields, &
         & width=options%width, options=options%dumpOptions )
       call showSummary
       call DestroyL2GPContents ( contractedl2gp )
     elseif ( options%chunks == '*' ) then
       call dump(l2gp, options%columnsOnly, options%details, options%fields, &
         & width=options%width, options=options%dumpOptions)
       call showSummary
     else
       call ExpandStringRange(options%chunks, chunks, nchunks)
       if ( nchunks < 1 ) return
       call dump(l2gp, chunks(1:nChunks), &
         & options%columnsOnly, options%details, options%fields, &
         &width=options%width, options=options%dumpOptions)
       call showSummary
     endif
   end subroutine myDump

   subroutine showPercentages( numTest, numGood, name, CutOff, aorb )
     ! output % figures derived from numTest/numGood
     ! Args
     real, intent(in)    ::          numTest
     integer, intent(in) ::          numGood
     real, intent(in) ::             cutOff
     character(len=*), intent(in) :: name
     character(len=*), intent(in) :: aorb ! 'above' or 'below'
     ! Internal variables
     ! Executable
     call output( '% ' )
     call output( trim(name) )
     call blanks(2)
     call output( aorb )
     call blanks(2)
     call output( cutOff )
     call blanks(1)
     call output( (100.*numTest)/max(1, numGood) )
     if ( options%debug ) then
       call blanks(3)
       call output( int(numTest) )
       call blanks(3)
       call output( numGood )
     endif
     call newline
   end subroutine showPercentages

   subroutine showStatusPct
     ! output % figures derived from numTest/numGood
     ! Internal variables
     integer :: bitindex
     logical, parameter :: DEEBUG = .false.
     ! Executable
     if ( DEEBUG ) then
       do bitindex=2, MAXNUMBITSUSED
         call outputNamedValue ( 'bitNumber(bitindex)', bitNumber(bitindex) )
       enddo
     endif
     if ( options%debug ) then
       call output( 'valid data counts' )
       call blanks(4)
       call output( bitCounts(1, 2), advance='yes' )
     endif
     call output( '% valid data with bits set', advance='yes' )
     call output( 'bit' )
     call blanks(4)
     call output( 'desc' )
     call blanks(4)
     call output( '%', advance='yes' )
     do bitindex=1, MAXNUMBITSUSED
       call output( bitNumber(bitIndex) )
       call blanks(2)
       call output( trim( stringElement( bitNames, bitIndex, &
         & countEmpty=.true. ) ) )
       call blanks(4)
       call output( (100.*bitCounts(bitindex, 1)) / max(1, bitCounts(bitindex, 2) ) )
       if ( options%debug ) then
         call blanks(4)
         call output( bitCounts(bitindex, 1) )
         call blanks(4)
         call output( bitCounts(bitindex, 2) )
       endif
       call newline
     enddo
   end subroutine showStatusPct

   subroutine showSummary
     ! output number good, unuseable profiles
     call outputNamedValue( 'Number of good profiles', numGood )
     call outputNamedValue( 'Number of unuseable profiles', numNotUseable )
     call outputNamedValue( 'Number with odd status set', numOddStatus )
   end subroutine showSummary

!==================
end program L2GPDump
!==================

! $Log$
! Revision 1.11  2011/05/26 20:47:55  pwagner
! Repaired use of Cutoffs for Quality and Precision (we hope)
!
! Revision 1.10  2011/02/18 23:10:27  pwagner
! Passes -d opts to dump routines
!
! Revision 1.9  2010/03/31 18:15:30  pwagner
! Removed commented-out, outdated stuff
!
! Revision 1.8  2009/05/14 22:04:27  pwagner
! New commandline arg -w width
!
! Revision 1.7  2009/04/13 20:43:17  pwagner
! Fixed a bug preventing macros file from using its own macros properly
!
! Revision 1.6  2008/12/03 00:15:04  pwagner
! Must use MLSFile_T interfaces instead of mls_io_gen_..
!
! Revision 1.5  2008/09/09 16:51:38  pwagner
! Added geolocation box to dump subsetted data
!
! Revision 1.4  2007/10/12 23:38:56  pwagner
! Shows num profiles good, unuseable, with odd status
!
! Revision 1.3  2007/06/14 21:45:42  pwagner
! Many bugs corrected regarding percentages
!
! Revision 1.2  2007/02/06 23:19:06  pwagner
! Can show percentages with thresholds in convergence, status, prec.
!
! Revision 1.1  2006/08/10 23:06:13  pwagner
! First commit
!
