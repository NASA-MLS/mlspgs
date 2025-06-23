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
program l2gpcat ! catenates split L2GPData files, e.g. dgg
!=================================

   use Dump_0, only: Dump
   use HDF, only: DFACC_Create, DFACC_RDonly, DFACC_RdWr
   use HDF5, only: H5FIs_HDF5_F
   use HighOutput, only: OutputNamedValue
   use Intrinsic, only: L_Swath
   use Io_Stuff, only: Read_TextFile
   use L2GPData, only: L2GPData_T, L2GPNameLen, MaxSwathNamesBufSize, RGP, &
     & AppendL2GPData, CpHE5GlobalAttrs, CpL2GPData, DestroyL2GPContents, &
     & ExtractL2GPRecord, ReadL2GPData, WriteL2GPData
   use Machine, only: Hp, Getarg
   use MLSCommon, only: MLSFile_T, L2MetaData_T
   use MLSFiles, only: MLS_Exists, MLS_CloseFile, MLS_OpenFile, &
     & HDFVersion_4, HDFVersion_5, MLS_InqSwath, InitializeMLSFile
   use MLSHDF5, only: MLS_H5Open, MLS_H5Close
   use MLSFillValues, only: Monotonize
   use MLSFinds, only: FindFirst, FindLast
   use MLSStringLists, only: CatLists, GetStringElement, &
     & Intersection, NumStringElements, RemoveElemFromList, &
     & StringElement, StringElementNum
   use MLSStrings, only: Asciify
   use Output_M, only: Output
   use PrintIt_M, only: Set_Config
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
! Catenate L2GPData from list of input files to a single output file
! May thus "unsplit" the dgg files

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files] -o [output file]
  type options_T
    logical     ::           timing = .false.           ! Show detailed timings
    logical     ::           verbose = .false.          ! print extra
    character(len=255) ::    glAttrFile= ''             ! file with global attrs
    character(len=255) ::    inputFile= ''              ! file with list of inputs       
    character(len=255) ::    outputFile= 'default.he5'  ! output filename       
    logical ::               append   = .false.         ! append swaths with same name
    logical ::               catenate = .false.         ! catenate swaths with same name
    logical ::               columnsOnly = .false.
    logical ::               ignoreFills = .false.      ! catenate swaths with same name
    logical ::               noDupSwaths = .false.      ! cp 1st, ignore rest   
    logical ::               monotonize = .true.        ! apply to geod. ang.
    character(len=3) ::      convert= ' '               ! e.g., '425'
    character(len=255) ::    swathNames = ' '           ! which swaths to copy
    character(len=255) ::    rename = ' '               ! how to rename them
    integer            ::    overlap = 0                ! max profiles in overlap
    integer, dimension(2) :: freqs = 0                  ! Keep range of freqs   
    integer, dimension(2) :: levels = 0                 ! Keep range of levels   
    integer, dimension(2) :: nProfiles = 0              ! Discard if nProfiles outside range
    integer, dimension(2) :: profiles = 0               ! Keep range of profiles   
  end type options_T
  
  type ( options_T ) :: options

  integer, parameter ::          MAXFILES = 750
  logical, parameter ::          countEmpty = .true.
  logical     :: createdYet
  logical, parameter ::          DEEBUG = .false.
  integer                                 :: file1Handle
  integer                                 :: file2Handle
  character(len=255) :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: fileNames
  type(MLSFile_T)                         :: L2GPFile1
  type(MLSFile_T)                         :: L2GPFile2
  type(MLSFile_T), dimension(MAXFILES)    :: L2GPFiles
  integer            :: n_filenames
  integer     ::  i, j, status, error ! Counting indices & Error flags
  integer     :: elem
  integer     ::  hdfversion1
  integer     ::  hdfversion2
  logical     :: is_hdf5
  integer, parameter :: S2US  = 1000000 ! How many microseconds in a s
  integer, parameter :: DELAY = 1*S2US  ! How long to sleep in microseconds
  real        :: t1
  real        :: t2
  real        :: tLast
  real        :: tFile
  character(len=255) ::    rename = ' '               ! how to rename them
  character(len=L2GPNameLen)          :: swath
  character(len=MaxSwathNamesBufSize) :: swathList
  character(len=MaxSwathNamesBufSize) :: swathList1
  character(len=MaxSwathNamesBufSize) :: swathListAll
  character(len=MaxSwathNamesBufSize) :: swathListOut
  type (L2Metadata_T) :: l2metaData
  integer :: listSize
  integer :: NUMSWATHSPERFILE
  integer :: NUMSWATHSSOFAR
  ! 
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  time_config%use_wall_clock = .true.
  CALL mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename(filename, n_filenames, options)
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     if ( mls_exists(trim(filename)) /= 0 ) then
       print *, 'Sorry--file not found: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  createdYet = .false.
  if ( len_trim(options%inputFile) > 0 ) then
    call read_textfile ( options%inputFile, filenames, nLines=n_filenames )
    ! Must replace all nulls with spaces
    do i=1, n_filenames
      filenames(i) = asciify( filenames(i), how='snip' )
    enddo
    if ( options%verbose ) then
      call dump( filenames(1:n_filenames), width=1, options='-t' )
    endif
  endif
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no input files to copy'
    stop
  endif
  ! Check that the hdfversions of the input files accord with convert mode
  status = 0
  do i=1, n_filenames
   call h5fis_hdf5_f(filenames(i), is_hdf5, error)
   select case (options%convert)
   case ('425')
     if ( is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf4 file: ', trim(filenames(i))
       status = 1
       cycle
     else
       hdfVersion1 = HDFVERSION_4
       hdfVersion2 = HDFVERSION_5
     endif
   case ('524')
     if ( .not. is_hdf5 ) then
       print *, 'Sorry--not recognized as hdf5 file: ', trim(filenames(i))
       status = 1
       cycle
     else
       hdfVersion1 = HDFVERSION_5
       hdfVersion2 = HDFVERSION_4
     endif
   case default
     if ( .not. is_hdf5 ) then
       hdfVersion1 = HDFVERSION_4
     else
       hdfVersion1 = HDFVERSION_5
     endif
     hdfVersion2 = hdfVersion1
   end select
  enddo
  call time_now ( t1 )

  swathListAll = ''
  numswathssofar = 0
  if ( options%catenate ) then
    call catenate_swaths
  elseif ( options%append ) then
    call append_swaths
  else
    call copy_swaths
  endif
  ! Do we copy global attributes from a particular file
  if ( len_trim(options%glAttrFile) > 0 ) then
    status = InitializeMLSFile( l2gpFile1, type=l_swath, access=DFACC_RDONLY, &
      & content='l2gp', name=options%glAttrFile, hdfVersion=HDFVERSION_5 )
    status = InitializeMLSFile( l2gpFile2, type=l_swath, access=DFACC_RDWR, &
      & content='l2gp', name=options%outputFile, hdfVersion=HDFVERSION_5 )
    call MLS_OpenFile( l2gpFile1 )
    call MLS_OpenFile( l2gpFile2 )
    File1Handle = l2gpFile1%FileID%f_id
    File2Handle = l2gpFile2%FileID%f_id
    call cpHE5GlobalAttrs ( File1Handle, File2Handle, status )
    if ( status /= 0 ) &
      & call output ( '(Global Attributes missing) ' // &
      & trim(options%glAttrFile), advance='yes')
    call MLS_CloseFile( l2gpFile1 )
    call MLS_CloseFile( l2gpFile2 )
  endif
  call mls_h5close(error)
contains
!------------------------- append_swaths ---------------------
!  Add to or overwrite existing swaths
  subroutine append_swaths
    ! Internal variables
    integer :: jj
    type( MLSFile_T ) :: l2gpFile
    type(L2GPData_T) :: ol2gp
    ! Executable
    status = InitializeMLSFile( l2gpFile, type=l_swath, access=DFACC_CREATE, &
      & content='l2gp', name='unknown', hdfVersion=HDFVERSION_5 )
    l2gpFile%name = options%outputFile
    l2gpFile%stillOpen = .false.
    if ( options%verbose ) print *, 'Append l2gp data to: ', &
      & trim(options%outputFile)
    numswathssofar = mls_InqSwath ( filenames(1), SwathList, listSize, &
           & hdfVersion=HDFVERSION_5)
    if ( len_trim(options%swathNames) > 1 ) swathList = options%swathNames
    swathListOut = swathList
    if ( len_trim(options%rename) > 1 ) swathListOut = options%rename
    if ( DEEBUG ) then
      print *, 'swaths in file'
      print *, trim(swathList)
    endif
    call GetStringElement( swathList, swath, 1, countEmpty )
    call GetStringElement( swathListOut, rename, 1, countEmpty )
    do jj=1, NumStringElements( swathList, countEmpty )
      call GetStringElement( swathList, swath, jj, countEmpty )
      call GetStringElement( swathListOut, rename, jj, countEmpty )
      call time_now ( tFile )
      if ( options%verbose ) print *, 'Appending swath: ', trim(swath)
      do i=1, n_filenames
        if ( options%verbose ) print *, 'Reading from: ', trim(filenames(i))
        call ReadL2GPData( trim(filenames(i)), swath, ol2gp )
        ! if ( i == 1 ) then
        ! What do we call the swath in the output file?
        call AppendL2GPData ( ol2gp, options%outputFile, rename )
      enddo

    enddo
  end subroutine append_swaths

!------------------------- catenate_swaths ---------------------
  subroutine catenate_swaths
    if ( options%ignoreFills ) then
      call catenate_non_Fills
    else
      call catenate_trimming_overlaps
    endif
  end subroutine catenate_swaths

!------------------------- catenate_non_Fills ---------------------
! catenate the profiles not containing Fills
! In terms of the profile numbers, we will 
! write starting with N+1 and continue up to M
! So, if a represents values already in the file, and 
! b represents value we will overwrite, then after overwriting we will have

! a[1]   ..  a[N] b[k] b[k+1] .. b[M] a[M+1] ..
! <---   N   --->  <---   M-k+1  --->

  subroutine catenate_non_Fills
    ! Internal variables
    integer :: i1
    integer :: i2
    integer :: iBloc
    integer :: jj
    type(L2GPData_T) :: l2gp
    type( MLSFile_T ) :: l2gpFile
    integer :: nBlocs
    integer, parameter :: numFilesPerBloc = 40
    integer :: numProfs
    integer :: numTotProfs
    type(L2GPData_T) :: ol2gp
    integer :: l2FileHandle
    integer :: N ! last profile number of a
    integer :: M ! size of b
    integer :: status
    logical :: wroteAlready
    ! Executable
    if ( options%verbose ) print *, 'Catenate l2gp data to: ', &
      & trim(options%outputFile)
    numswathssofar = mls_InqSwath ( filenames(1), SwathList, listSize, &
           & hdfVersion=HDFVERSION_5)
    if ( DEEBUG ) then
      print *, 'swaths in file'
      print *, trim(swathList)
    endif
    if ( options%swathNames /= ' ' ) then
      swathList1 = swathList
      swathList = Intersection( options%swathNames, swathList1 )
    endif
    status = InitializeMLSFile( l2gpFile, type=l_swath, access=DFACC_CREATE, &
      & content='l2gp', name='unknown', hdfVersion=HDFVERSION_5 )
    l2gpFile%name = options%outputFile
    l2gpFile%stillOpen = .false.
    call mls_openFile( L2GPFile, Status )
    l2FileHandle = l2gpFile%FileID%f_id
    nBlocs = (n_filenames - 1)/numFilesPerBloc + 1

    call GetStringElement( swathList, swath, 1, countEmpty )
    do jj=1, NumStringElements( swathList, countEmpty ) ! Loop of swaths
      call GetStringElement( swathList, swath, jj, countEmpty )
      call time_now ( tFile )
      if ( options%verbose ) print *, 'Catenating swath: ', trim(swath)
      N = 0
      M = 0
      wroteAlready = .false.
      numTotProfs = 0
      i2 = 0
      do iBloc = 1, nBlocs ! Loop of blocs
        i1 = i2 + 1
        i2 = min( i2 + numFilesPerBloc, n_filenames )
        ! 1st open the files
        do i=i1, i2
          status = InitializeMLSFile( L2GPFiles(i), type=l_swath, access=DFACC_RDONLY, &
            & content='l2gp', name=filenames(i), hdfVersion=HDFVERSION_5 )
          L2GPFiles(i)%name = filenames(i)
          L2GPFiles(i)%stillOpen = .false.
          call mls_openFile( L2GPFiles(i), Status )
        enddo
        ! Read each file, putting its values into ol2gp
        do i=i1, i2
          if ( options%verbose ) print *, 'Reading from: ', trim(filenames(i))
          if ( i == 1 ) then
            call ReadL2GPData( L2GPFiles(i)%FileID%f_id, swath, ol2gp, numProfs )
            cycle
          endif
          call ReadL2GPData( L2GPFiles(i)%FileID%f_id, swath, l2gp, numProfs )
          ! We will use ChunkNumbers to determine N and M
          N = FindFirst( l2gp%chunkNumber > 0 ) - 1
          M = FindLast( l2gp%chunkNumber > 0 )
          ol2gp%latitude       (N+1:M)     = l2gp%latitude       (N+1:M)
          ol2gp%longitude      (N+1:M)     = l2gp%longitude      (N+1:M)
          ol2gp%solarTime      (N+1:M)     = l2gp%solarTime      (N+1:M)
          ol2gp%solarZenith    (N+1:M)     = l2gp%solarZenith    (N+1:M)
          ol2gp%losAngle       (N+1:M)     = l2gp%losAngle       (N+1:M)
          ol2gp%geodAngle      (N+1:M)     = l2gp%geodAngle      (N+1:M)
          ol2gp%time           (N+1:M)     = l2gp%time           (N+1:M)

          ol2gp%chunkNumber    (N+1:M)     = l2gp%chunkNumber      (N+1:M)
          ol2gp%status         (N+1:M)     = l2gp%status           (N+1:M)
          ol2gp%quality        (N+1:M)     = l2gp%quality          (N+1:M)
          ol2gp%convergence    (N+1:M)     = l2gp%convergence      (N+1:M)  
          ol2gp%AscDescMode    (N+1:M)     = l2gp%AscDescMode      (N+1:M)  

          ol2gp%l2gpValue      (:,:,N+1:M) = l2gp%l2gpValue        (:,:,N+1:M)    
          ol2gp%l2gpPrecision  (:,:,N+1:M) = l2gp%l2gpPrecision    (:,:,N+1:M)    

          call DestroyL2GPContents( l2gp )
        enddo
        ! Last, close each file
        do i=i1, i2
          call mls_closeFile( L2GPFiles(i), Status )
        enddo
      enddo  ! End Loop of blocs
      if ( options%timing ) call sayTime('Reading this swath', tFile)
      call WriteL2GPData ( ol2gp, l2FileHandle, swath )
      if ( options%timing ) call sayTime('Writing this swath', tFile)
      call DestroyL2GPContents( ol2gp )
    enddo ! End Loop of swaths
    call mls_closeFile( L2GPFile, Status )
    if ( options%timing ) call sayTime('catenating all swaths')
  end subroutine catenate_non_Fills

!------------------------- catenate_trimming_overlaps ---------------------
! Identify and trim overlaps as follows:
! let {a1 a2 .. aN} be an array of geodetic angles in one
! swath while {b1 b2 .. bM} is an array in the next swath
! where a and b overlap by a certain number of profiles which we will discover:
! let k be the 1st b s.t. (b - a[N]) > eps
! Then keeping the a array values in the overlap region will result in
!
! a[1]   ..  a[N] b[k] b[k+1] .. b[M] a[M+1] ..
! <---   N   --->  <---   M-k+1  --->

! A later wrinkle:
! We may choose not to keep the full N profiles from the previous
! swath, instead lopping off a few (kLopOff in number)
! a[1]   ..  a[N-kL] b[k-kL] b[kL+1] .. b[M]
! <---   N-kL   --->  <---   M-k+kL+1  --->

  subroutine catenate_trimming_overlaps
    ! Internal variables
    integer :: j
    integer :: jj
    integer :: k
    integer :: kLopOff
    type(L2GPData_T) :: l2gp
    type( MLSFile_T ) :: l2gpFile
    integer, parameter :: MAXNUMPROFS = 3500
    integer :: numProfs
    integer :: numTotProfs
    type(L2GPData_T) :: ol2gp
    real(rgp), dimension(MAXNUMPROFS) :: a
    real(rgp), dimension(MAXNUMPROFS) :: b
    integer :: N ! size of a
    integer :: M ! size of b
    logical :: wroteAlready
    real(rgp), parameter :: eps = 0.01_rgp
    ! Executable
    a = 0.
    b = 0.
    if ( options%verbose ) print *, 'Catenate l2gp data to: ', &
      & trim(options%outputFile)
    numswathssofar = mls_InqSwath ( filenames(1), SwathList, listSize, &
           & hdfVersion=HDFVERSION_5)
    if ( DEEBUG ) then
      print *, 'swaths in file'
      print *, trim(swathList)
    endif
    if ( options%swathNames /= ' ' ) then
      swathList1 = swathList
      swathList = Intersection( options%swathNames, swathList1 )
    endif
    call GetStringElement( swathList, swath, 1, countEmpty )
    do jj=1, NumStringElements( swathList, countEmpty )
      ! For an unexplained reason, in our tests, 
      ! everything slows down after 48 swaths
      ! Here is something we aree trying
      if ( mod( jj, 10 ) < -1 ) then
        call mls_h5close( error )
        call mls_h5open( error )
      endif
      call GetStringElement( swathList, swath, jj, countEmpty )
      call time_now ( tFile )
      tLast = tFile
      if ( options%verbose ) print *, 'Catenating swath: ', trim(swath)
      N = 0
      M = 0
      wroteAlready = .false.
      numTotProfs = 0
      do i=1, n_filenames
        if ( options%verbose ) print *, 'Reading from: ', trim(filenames(i))
        call ReadL2GPData( trim(filenames(i)), swath, ol2gp, numProfs )
        if ( options%timing ) call sayTime( 'Reading swath ' // trim(swath), tLast )
        tlast = t2
        ! Check if nProfiles is in range
        if ( any(options%nProfiles /= 0) ) then
          if ( ol2gp%nTimes < options%nProfiles(1) .or. &
            & ol2gp%nTimes > options%nProfiles(2) ) then
            print *, 'Warning--skipping ', ol2gp%nTimes, ' profiles in ', trim(filenames(i))
            call DestroyL2GPContents( ol2gp )
            cycle
          endif
        endif
        if ( options%monotonize ) &
          & call Monotonize( ol2gp%GeodAngle )
          ! & call Monotonize( ol2gp%GeodAngle, FillValue=-999.99 )
        ! if ( i == 1 ) then
        if ( .not. wroteAlready ) then
          status = InitializeMLSFile( l2gpFile, type=l_swath, access=DFACC_CREATE, &
            & content='l2gp', name='unknown', hdfVersion=HDFVERSION_5 )
          l2gpFile%name = options%outputFile
          l2gpFile%stillOpen = .false.
          if ( DEEBUG ) then
            print *, 'About to write ', trim(swath)
          endif
          call WriteL2GPData ( ol2gp, l2gpFile, swath )
          if ( options%timing ) call sayTime( 'Writing swath', tLast )
          tlast = t2
          N = ol2gp%nTimes
          a(1:N) = ol2gp%GeodAngle
          if ( options%verbose ) call dump( a(1:N), '1st a' )
          call DestroyL2GPContents( ol2gp )
          numTotProfs = numTotProfs + N
          wroteAlready = .true.
          cycle
        endif
        M = ol2gp%nTimes
        b(1:M) = ol2gp%GeodAngle
        ! 1st, is there nearly 360 deg offset between a, and b?
        ! (If so, we'll assume a periodic wrap occurred)
        if ( abs(b(1)-a(N)) > 300._rgp ) then
        ! if ( abs(b(1)-maxval(a)) > 300._rgp ) then
          ! add or subtract the requsite 360 deg from a
          j = int( sign(1.01_rgp, b(1)-a(N)) )
          if ( DEEBUG ) call outputNamedvalue ( 'j', j )
          a(1:N) = a(1:N) + j*360._rgp
        endif
        ! Find k
        k = FindFirst( (b(1:M) - a(N)) > eps )
        kLopOff = 0
        ! print *, 'k: ', k
        if ( options%overlap > 0 ) then
          kLopOff = k - min( k, options%overlap )
          k = k - kLopOff
        endif
        ! Are our overlaps uniform?
        if ( k < 2 ) then
          if ( options%verbose ) print *, 'Warning--overlaps not found', k, N, M
          if ( options%verbose ) call dump( a(1:N), 'a' )
          if ( options%verbose ) call dump( b(1:M), 'b' )
          if ( DEEBUG ) print *, 'About to append ', trim(swath)
          call AppendL2GPData ( ol2gp, options%outputFile, &
          & swath, offset=max(0,numTotProfs) )
          if ( options%timing ) call sayTime( 'Appending swath (k < 2 )', tLast )
          tlast = t2
          ! Could this be a bug in the HDFEOS library?
          !if ( i > 1 .and. ol2gp%nTimes > numTotProfs ) &
          !  & call AppendL2GPData ( ol2gp, options%outputFile, &
          !  & swath, offset=max(0,numTotProfs) )
          numTotProfs = numTotProfs + M
        else
          ! print *, 'k (before extraction): ', k
          call ExtractL2GPRecord ( ol2gp, l2gp, rTimes=(/ k, M /) )
          if ( options%timing ) &
            & call sayTime( 'Extracting record (k > 1 )', tLast )
          tlast = t2
          ! print *, 'k (after extraction): ', k
          ! print *, 'should be writing'
          if ( options%verbose ) &
            & call dump( l2gp%GeodAngle, 'appended Geod. angle' )
          if ( DEEBUG ) print *, 'About to append ', trim(swath)
          ! call usleep ( delay ) ! Should we make this parallel%delay?
          call AppendL2GPData ( l2gp, options%outputFile, &
          & swath, offset=max(0,numTotProfs-kLopOff) )
          if ( options%timing ) call sayTime( 'Appending swath (k > 1 )', tLast )
          tlast = t2
          ! Could this be a bug in the HDFEOS library?
          !if ( i > 1 .and. l2gp%nTimes > numTotProfs .and. .false. ) &
          !  & call AppendL2GPData ( l2gp, options%outputFile, &
          !  & swath, offset=max(0,numTotProfs) )
          call DestroyL2GPContents( l2gp )
          numTotProfs = numTotProfs + M - k - kLopOff + 1
          ! print *, 'k,N,M,total', k, N, M, numTotProfs
        endif
        call DestroyL2GPContents( ol2gp )
        ! Next time what was "b" will serve as "a"
        a = b
        N = M
        if ( DEEBUG ) then
          call ReadL2GPData( options%outputFile, swath, ol2gp )
          ! print *, 'nTimes(total): ', ol2gp%nTimes
          if ( options%verbose ) &
            & call dump( ol2gp%GeodAngle, 'stored Geod. angle' )
          call DestroyL2GPContents( ol2gp )
        endif
      enddo
      if ( options%timing ) call sayTime('catenating this swath', tFile)
    enddo
    if ( options%timing ) call sayTime('catenating all swaths')
  end subroutine catenate_trimming_overlaps

!------------------------- copy_swaths ---------------------
  subroutine copy_swaths
    ! logical, parameter :: DEEBUG = .true.
    if ( options%verbose ) &
      & print *, 'Copy l2gp data to: ', trim(options%outputFile)
    do i=1, n_filenames
      call time_now ( tFile )
      if ( options%verbose ) print *, 'Copying from: ', trim(filenames(i))
      if ( options%verbose .and. len_trim(options%swathNames) > 0 ) &
        & print *, 'swath names: ', trim(options%swathNames)
      if ( options%noDupSwaths .or. options%swathNames /= ' ' ) then
        numswathsperfile = mls_InqSwath ( trim(filenames(i)), &
          & swathList, listSize, hdfVersion=hdfVersion1)
        if ( DEEBUG ) then
          print *, 'swaths in file'
          print *, trim(swathList)
        endif
        if ( DEEBUG ) then
          print *, 'all swaths'
          print *, trim(swathListAll)
        endif
        if ( options%swathNames /= ' ' ) then
          swathList1 = swathList
          swathList = Intersection( options%swathNames, swathList1 )
          if ( swathList == ' ' ) cycle
          rename = ' '
          do j=1, NumStringElements( swathList, countEmpty )
            call GetStringElement( swathList, swath, j, countEmpty )
            elem = StringElementNum( options%swathNames, swath, countEmpty )
            rename = catLists( rename, &
              & StringElement( options%rename, elem, countempty ) )
          enddo
        elseif ( numswathssofar > 0 ) then
          ! Remove any duplicates
          do j=1, numswathssofar
            call GetStringElement(swathListAll, swath, j, countEmpty)
            swathList1 = swathList
            if ( DEEBUG ) then
              print *, 'Removing ', trim(swath)
              print *, trim(swathListAll)
            endif
            call RemoveElemFromList (swathList1, swathList, trim(swath))
            ! Crude hAck--really should fix removeElem procedure
            if ( swathList(1:1) == ',' ) then
              swathList1 = swathList(2:)
              swathList = swathList1
            endif
          enddo
          if ( DEEBUG ) then
            print *, 'swaths to cp'
            print *, trim(swathList)
          endif
        endif
        if (len_trim(swathList) < 1 ) return ! Because there are none to copy
        if ( any( (/options%freqs(2), options%levels(2), &
          & options%profiles(2)/) > 0 ) &
          & ) then
          call cpL2GPData( l2metaData, trim(filenames(i)), &
          & trim(options%outputFile), create2=.not. createdYet, &
          & hdfVersion1=hdfVersion1, hdfVersion2=hdfVersion2, &
          & swathList=trim(swathList), rename=rename, &
          & notUnlimited=.true., andGlAttributes=.true., &
          & rFreqs=options%freqs, rLevels=options%levels, rTimes=options%profiles)
        else
          call cpL2GPData( l2metaData, trim(filenames(i)), &
          & trim(options%outputFile), create2=.not. createdYet, &
          & hdfVersion1=hdfVersion1, hdfVersion2=hdfVersion2, &
          & swathList=trim(swathList), rename=rename, &
          & notUnlimited=.true., andGlAttributes=.true.)
        endif
        swathList1 = swathListAll
        swathListAll = catlists(swathList1, swathList)
        numswathssofar = NumStringElements(swathListAll, countEmpty)
      else
        if ( any( (/options%freqs(2), options%levels(2), options%profiles(2)/) &
          & > 0 ) ) then
          call cpL2GPData( l2metaData, trim(filenames(i)), &
          & trim(options%outputFile), create2=.not. createdYet, &
          & hdfVersion1=hdfVersion1, hdfVersion2=hdfVersion2, &
          & notUnlimited=.true., andGlAttributes=.true., &
          & rFreqs=options%freqs, rLevels=options%levels, rTimes=options%profiles)
        else
          call cpL2GPData( l2metaData, trim(filenames(i)), &
          & trim(options%outputFile), create2=.not. createdYet, &
          & hdfVersion1=hdfVersion1, hdfVersion2=hdfVersion2, &
          & notUnlimited=.true., andGlAttributes=.true.)
        endif
      endif
      if ( options%timing ) call sayTime('copying this file', tFile)
      createdYet = .true.
    enddo
    if ( options%timing ) call sayTime('copying all files')
  end subroutine copy_swaths
!------------------------- get_filename ---------------------
    subroutine get_filename( filename, n_filenames, options )
    ! Added for command-line processing
     character(LEN=255), intent(out) :: filename          ! filename
     integer, intent(in)             :: n_filenames
     type ( options_T ), intent(inout) :: options
     ! character(LEN=*), intent(inout) :: outputFile        ! output filename
     ! Local variables
     integer ::                         error = 1
     integer, save ::                   i = 1
  ! Get inputfile name, process command-line args
  ! (which always start with -)
    do
      call getarg ( i+hp, filename )
      ! print *, i, ' th Arg: ', trim(filename)
      error = 0
      if ( filename(1:1) /= '-' ) exit
      if ( filename(1:3) == '-h ' ) then
        call print_help
      elseif ( filename(1:3) == '-o ' ) then
        call getarg ( i+1+hp, options%outputFile )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-425 ' ) then
        options%convert = '425'
        exit
      elseif ( filename(1:5) == '-524 ' ) then
        options%convert = '524'
        exit
      elseif ( filename(1:3) == '-t ' ) then
        options%timing = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      elseif ( filename(1:4) == '-app' ) then
        options%append = .true.
        exit
      elseif ( filename(1:4) == '-cat' ) then
        options%catenate = .true.
        exit
      elseif ( filename(1:4) == '-ign' ) then
        options%ignoreFills = .true.
        exit
      elseif ( filename(1:3) == '-no' ) then
        options%noDupSwaths = .true.
        exit
      else if ( filename(1:3) == '-F ' ) then
        call getarg ( i+1+hp, options%inputFile )
        i = i + 1
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:3) == '-g ' ) then
        call getarg ( i+1+hp, options%glAttrFile )
        i = i + 1
      else if ( filename(1:3) == '-s ' ) then
        call getarg ( i+1+hp, options%swathNames )
        i = i + 1
        exit
      else if ( filename(1:3) == '-r ' ) then
        call getarg ( i+1+hp, options%rename )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-freq' ) then
        call igetarg ( i+1+hp, options%freqs(1) )
        i = i + 1
        call igetarg ( i+1+hp, options%freqs(2) )
        i = i + 1
        exit
      elseif ( filename(1:6) == '-level' ) then
        call igetarg ( i+1+hp, options%levels(1) )
        i = i + 1
        call igetarg ( i+1+hp, options%levels(2) )
        i = i + 1
        exit
      elseif ( filename(1:6) == '-nprof' ) then
        call igetarg ( i+1+hp, options%nProfiles(1) )
        i = i + 1
        call igetarg ( i+1+hp, options%nProfiles(2) )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-over' ) then
        call igetarg ( i+1+hp, options%overlap )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-prof' ) then
        call igetarg ( i+1+hp, options%profiles(1) )
        i = i + 1
        call igetarg ( i+1+hp, options%profiles(2) )
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
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
    write (*,*) &
    & 'Usage:l2gpcat [options] [filenames]'
    write (*,*) 'Options: '
    write (*,*) '-f filename   => add filename to list of filenames'
    write (*,*) '                 (can do the same w/o the -f)'
    write (*,*) '-F infile     => read list of filenames from infile'
    write (*,*) '-g glattrfile => read global attrs from glattrfile'
    write (*,*) '-o ofile      => copy swaths to ofile'
    write (*,*) '-425          => convert from hdf4 to hdf5'
    write (*,*) '-524          => convert from hdf5 to hdf4'
    write (*,*) '-t            => show detailed timings'
    write (*,*) '-v            => switch on verbose mode'
    write (*,*) '-append       => append or overwrite swaths with same name'
    write (*,*) '-cat          => catenate swaths with same name'
    write (*,*) '-ignoreFills  => ignore profiles with Fill Values'
    write (*,*) '                 useful for combining independently-run chunks'
    write (*,*) '-nodup        => if dup swath names, cp 1st only'
    write (*,*) '-freqs m n    => keep only freqs in range m n'
    write (*,*) '-levels m n   => keep only levels in range m n'
    write (*,*) '-profiles m n => keep only profiles in range m n'
    write (*,*) '-nprofiles m n => keep only if nProfs is in range m n'
    write (*,*) '-s name1,name2,..'
    write (*,*) '              => copy only swaths so named; otherwise all'
    write (*,*) '-r rename1,rename2,..'
    write (*,*) '              => if and how to rename the copied swaths'
    write (*,*) '-overlap n    => max num profiles in overlap'
    write (*,*) '-h            => print brief help'
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
!------------------------- igetarg ---------------------
  subroutine igetarg ( pos, iarg )
   integer, intent(in) :: pos
   integer, intent(out) :: iarg
   character(len=16) :: arg
   call getarg ( pos, arg )
   read(arg, *) iarg
  end subroutine igetarg

!==================
end program L2GPcat
!==================

! $Log$
! Revision 1.26  2017/08/25 00:20:21  pwagner
! Fixed bug when rolling up new dual-phase nrt
!
! Revision 1.25  2017/05/17 22:21:59  pwagner
! Works properly with dual-l2 nrt scripts
!
! Revision 1.24  2016/08/25 22:52:52  pwagner
! Can now successfully cat all 351 chunks of Pleiades-style dgg
!
! Revision 1.23  2016/06/13 23:27:46  pwagner
! Added -g commandline option
!
! Revision 1.22  2016/06/10 16:13:18  pwagner
! Added commandline option -F; upped max num of input files tto 750
!
! Revision 1.21  2016/02/11 19:53:20  pwagner
! options to turn timing, verbose mode on
!
! Revision 1.20  2015/08/05 20:37:03  pwagner
! Option -ign can catenate non-Fills
!
! Revision 1.19  2014/09/12 22:22:14  pwagner
! Fixed sense errors in len_trim tests
!
! Revision 1.18  2014/09/12 00:04:08  pwagner
! Added -append commandline option to overwrite swath values in target file
!
! Revision 1.17  2014/01/09 00:31:26  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 1.16  2013/08/23 02:51:47  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.15  2013/05/30 20:41:09  pwagner
! Reduce amount of debug printing
!
! Revision 1.14  2009/10/27 21:10:03  pwagner
! Brought into compliance with cpL2GPData api
!
! Revision 1.13  2008/10/13 23:32:41  pwagner
! Changed meaning of -overlap option; useful now in rolling up nrts
!
! Revision 1.12  2008/09/25 23:13:20  pwagner
! May exclude swaths when nProfiles outside range
!
! Revision 1.11  2008/02/28 01:36:24  pwagner
! -cat option catenates different pieces of same swath
!
! Revision 1.10  2006/05/19 22:47:40  pwagner
! Fixed bug stopping us from creating -o file if first input file not copied
!
! Revision 1.9  2006/05/19 20:55:36  pwagner
! May rename copied swaths
!
! Revision 1.8  2006/04/06 23:04:21  pwagner
! Optionally cp only ranges of freq, level, profile
!
! Revision 1.7  2005/10/29 00:13:56  pwagner
! Removed unused procedures from use statements
!
! Revision 1.6  2005/09/23 21:01:13  pwagner
! use_wall_clock now a component of time_config
!
! Revision 1.5  2005/06/22 19:27:33  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.4  2004/12/06 19:13:12  pwagner
! With -nodup option ignores dup swath names after 1st
!
! Revision 1.3  2004/08/07 00:15:55  pwagner
! All stringlist stuff was moved from mlsstrings to mlsstringlists
!
! Revision 1.2  2004/05/06 21:50:48  pwagner
! Uses mls_h5open/close
!
! Revision 1.1  2004/04/30 18:54:22  pwagner
! First commit
!
