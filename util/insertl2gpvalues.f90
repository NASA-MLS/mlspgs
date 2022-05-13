! Copyright 2022, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=================================
program insertL2GPValues ! inserts values of L2GPData files, e.g. nrt
!=================================

   use Dump_0, only: Dump
   use HDF, only: Dfacc_RdOnly, Dfacc_Rdwr
   use HDF5, only: HSize_T
   use HighOutput, only: OutputNamedValue
   use Intrinsic, only: L_HDF, L_Swath
   use L2GPData, only: L2GPData_T, L2GPnamelen, MaxSwathNamesBufSize, Rgp, &
     & AppendL2GPData, DestroyL2GPContents, Dump, &
     & ReadL2GPData
   use Machine, only: Hp, Getarg
   use MLSCommon, only: MLSFile_T
   use MLSFiles, only: HDFVersion_5, Dump, MLS_CloseFile, MLS_Exists, &
     & MLS_Inqswath, InitializeMLSFile
   use MLSFinds, only: FindFirst
   use MLSHDF5, only: GetAllHDF5DSNames, GetHDF5DSDims, LoadFromHDF5DS, &
     & MLS_H5open, MLS_H5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMessage, MLSMSG_Warning
   use MLSStringLists, only: GetStringElement, &
     & Intersection, NumStringElements, StringElementNum
   use Output_M, only: Blanks, Output
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
! Insert new values into L2GPData from list of files
! using values read from a "new values" hdf file

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"

! Then run it something like this
!    /users/pwagner/mlspgs/tests/lib/NAG.Linux-CemtOS7/test \ 
!       -s H2O-StdProd -d ANN_Prediction -p ANN_Precision \
!       -DGM MLS-Aura_L2GP-DGM_*.h5 \
!       -Vf h2o_prediction.h5 \
!       -Wf /users/fwerner/Documents/database/neural_network_weights/weights/v05-0x/MLS-Aura_ANN-H2O_v05-0x-01_20220415.h5 \
!       MLS-Aura_L2GP-DGG_*.he5 

! We'll use it to modify existing swaths, giving them new values
! and precisions predicted by a different retrieval, e.g.
! a trained neural network.

! To create a file of predicted values, you must run's Frank's python script.
! A typical way is the following
! /testing/workspace/pwagner/t2--t-negovlaps 179: echo /usr/local/anaconda3/bin/python3 $H2ONNSCRIPT $L1BRADG $L1BRADD $L1BOA $Weights_File $Prediction_File
! /usr/local/anaconda3/bin/python3 /users/fwerner/Documents/software/python/ann/h2o_prediction.py /data/emls/l1b/v05.01/2022/009/MLS-Aura_L1BRADG_v05-01-c01_2022d009.h5 /data/emls/l1b/v05.01/2022/009/MLS-Aura_L1BRADD_v05-01-c01_2022d009.h5 /data/emls/l1b/v05.01/2022/009/MLS-Aura_L1BOA_v05-01-c01_2022d009.h5 /users/fwerner/Documents/database/neural_network_weights/weights/v05-0x/MLS-Aura_ANN-H2O_v05-0x-01_20220415.h5 prediction.h5

  type options_T
    logical               :: debug   = .false.
    logical               :: verbose = .true.
    logical               :: dryrun  = .false.
    character(len=255)    :: swathNames = ' '        ! which swath to modify
    character(len=255)    :: DSNames    = ' '        ! which new DS holds values
    character(len=255)    :: PrecNames  = ' '        ! which new DS holds precs
    character(len=255)    :: newValues  = ' '        ! file with new values
    character(len=255)    :: L1BOAFile  = ' '        ! L1BOA file
    character(len=255)    :: DGMFile    = ' '        ! DGM file
    character(len=255)    :: Weights    = ' '        ! Weights file
  end type options_T
  
  type ( options_T ) :: options

  integer, parameter ::          MAXFILES = 10 ! Usually 1
  logical, parameter ::          countEmpty = .true.
  character(len=255) :: filename          ! input filename
  character(len=255), dimension(MAXFILES) :: filenames
  integer            :: n_filenames
  integer     ::  status, error ! Counting indices & Error flags
  real        :: t1
  real        :: t2
  real        :: tFile
  character(len=L2GPNameLen)          :: DSName
  character(len=MAXSWATHNAMESBUFSIZE) :: DSList
  character(len=L2GPNameLen)          :: swath
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList1
  integer :: listSize
  integer :: NUMSWATHSSOFAR
  integer, allocatable, dimension(:)  :: nearestProfiles
  integer, allocatable, dimension(:)  :: Output_Pressure_Levels_Indices
  real (rgp), dimension(:,:), allocatable :: GeodAngle
  real (rgp), dimension(:,:), allocatable :: GeodLat
  ! 
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  time_config%use_wall_clock = .true.
  CALL mls_h5open(error)
  n_filenames = 0
  do      ! Loop over filenames
     call get_filename( filename, n_filenames, options )
     if ( filename(1:1) == '-' ) cycle
     if ( filename == ' ' ) exit
     if ( mls_exists(trim(filename)) /= 0 ) then
       print *, 'Sorry--file not found: ', trim(filename)
       cycle
     endif
     n_filenames = n_filenames + 1
     filenames(n_filenames) = filename
  enddo
  if ( n_filenames == 0 ) then
    if ( options%verbose ) print *, 'Sorry no files to insert'
  else
    call time_now ( t1 )

    numswathssofar = 0
    if ( options%verbose ) &
      & call Dump ( filenames(1:n_filenames), 'l2gp files', width=1 )
    call insert_swaths
  endif
  call mls_h5close(error)
contains
!------------------------- insert_swaths ---------------------
! Identify and insert values

  subroutine insert_swaths
    ! Internal variables
    integer :: i
    integer :: iSwathNum
    integer :: jj
    type(L2GPData_T) :: L2GP
    type( MLSFile_T ) :: L2GPFile
    type( MLSFile_T ) :: DGMFile
    type( MLSFile_T ) :: L1BOAFile
    type( MLSFile_T ) :: valuesFile
    type( MLSFile_T ) :: WeightsFile
    ! type (L1BData_T)  :: l1bField ! L1B data
    integer :: nPr
    integer :: numProfs
    integer :: numTotProfs
    integer :: N ! size of a
    integer, dimension(:,:,:), allocatable  :: neartemp
    real (rgp), dimension(:,:), allocatable :: precisions
    real (rgp), dimension(:,:), allocatable :: values
    integer(kind=hSize_t), dimension(3)     :: DIMS
    ! --------------------------------------------------------------------
    ! This is the number of instances in the new values file
    integer :: nvTimes
    ! --------------------------------------------------------------------
    ! These are the number of pressure levels each has been given
    ! integer, dimension(nDS), parameter :: NPrs = (/ 41, 14, 12, 13, 37, 41 /)
    
    ! Note that in the latest netCDF file, all species are on 37 pressure
    ! surfaces except for the _HR species. This choice happens to match
    ! the resolution of our l2gp products. Will future netCDF files from
    ! Joaquim adhere to this convention, too?
    integer, dimension(1), parameter :: NPrs = (/ 42 /)
    ! --------------------------------------------------------------------

    ! --------------------------------------------------------------------
    ! Executable
    ! ------------------------------------------------------------------------
    ! 1st data sets to read:    GHz/GeodLat and GHz/GeodAngle
    ! read from                        L1BOA File
    status = InitializeMLSFile( L1BOAFile, type=l_hdf, access=DFACC_RDOnly, &
      & content='hdf', name=trim(options%L1BOAFile), hdfVersion=HDFVERSION_5 )
    call GetAllHDF5DSNames ( L1BOAFile, DSList )
    if ( options%verbose ) then
      print *, 'DS names in l1boa file ' // trim(DSList)
    endif
    DSName = 'GHz/GeodLat'
    call GetHDF5DSDims ( L1BOAFile, DSName, DIMS )
    NvTimes = DIMS(2)
    if ( options%verbose ) print *, 'DIMS: ', DIMS
    allocate(GeodLat(125,nvTimes))
    call LoadFromHDF5DS ( L1BOAFile, DSName, GeodLat )
    if ( options%verbose ) call dump( GeodLat(1,:), 'GeodLat' )
    DSName = 'GHz/GeodAngle'
    allocate(GeodAngle(125,nvTimes))
    call LoadFromHDF5DS ( L1BOAFile, DSName, GeodAngle )
    if ( options%verbose ) call dump( GeodAngle(1,:), 'GeodAngle' )
    
    ! 1st data set to read:          nearest profile
    ! read from                        DGM File
    ! Alas, this array fails to do what it promised
    if ( .false. ) then
      status = InitializeMLSFile( DGMFile, type=l_hdf, access=DFACC_RDOnly, &
        & content='hdf', name=trim(options%DGMFile), hdfVersion=HDFVERSION_5 )
      call GetAllHDF5DSNames ( DGMFile, DSList )
      if ( options%verbose ) then
        print *, 'DS names in hdf file ' // trim(DSList)
      endif
      DSName = 'nearest profile'
      call GetHDF5DSDims ( DGMFile, DSName, DIMS )
      NvTimes = DIMS(3)
      if ( options%verbose ) print *, 'DIMS: ', DIMS
      allocate(neartemp(1,1,nvTimes))
      allocate(nearestProfiles(nvTimes))
      call LoadFromHDF5DS ( DGMFile, DSName, neartemp )
      nearestProfiles = neartemp(1,1,:)
      if ( options%verbose ) call dump( nearestProfiles, 'nearestProfiles' )
      deallocate(neartemp)
    endif
    
    ! ------------------------------------------------------------------------
    ! 2nd data set to read:          Output_Pressure_Levels_Indices
    ! read from                        Weights File
    
    ! We don't need these if the prediction and precision arrays
    ! already fill in the vertical levels missing from the Weights file.
    if ( .false. ) then
      status = InitializeMLSFile( WeightsFile, type=l_hdf, access=DFACC_RDOnly, &
        & content='hdf', name=trim(options%Weights), hdfVersion=HDFVERSION_5 )
      call GetAllHDF5DSNames ( WeightsFile, DSList )
      if ( options%verbose ) then
        print *, 'DS names in hdf file ' // trim(DSList)
      endif
      DSName = 'Output_Pressure_Levels_Indices'
      call GetHDF5DSDims ( WeightsFile, DSName, DIMS )
      NPr = DIMS(1)
      if ( options%verbose ) print *, 'DIMS: ', DIMS
      allocate(Output_Pressure_Levels_Indices(Npr))
      call LoadFromHDF5DS ( WeightsFile, DSName, Output_Pressure_Levels_Indices )
    endif
    
    ! ------------------------------------------------------------------------
    ! 3rd and 4th data sets to read:    values and precisions
    ! read from                                Values File
    status = InitializeMLSFile( valuesFile, type=l_hdf, access=DFACC_RDOnly, &
      & content='hdf', name=trim(options%newValues), hdfVersion=HDFVERSION_5 )
    call GetAllHDF5DSNames ( valuesFile, DSList )
    if ( options%verbose ) then
      print *, 'DS names in hdf file ' // trim(DSList)
    endif
    DSName = options%DSNames
    call GetHDF5DSDims ( valuesFile, DSName, DIMS )
    NPr     = DIMS(1)
    NvTimes = DIMS(2)
    if ( options%verbose ) print *, 'DIMS: ', DIMS
    if ( options%verbose ) print *, 'Num swaths: ', NumStringElements( DSList, countEmpty )
    numTotProfs = 0
    ! do jj=1, NumStringElements( DSList, countEmpty )
    jj = 1
      call GetStringElement( DSList, DSName, jj, countEmpty )
      call time_now ( tFile )
      swath = options%swathNames
      if ( options%verbose ) print *, 'inserting values into swath: ', trim(swath)
      if ( options%verbose ) print *, 'using values from DS: ', trim(swath)
      call Blanks( 72, fillChar='-', advance='yes' )
      if ( options%verbose ) then
        call OutputNamedValue ( 'DSName', swath )
      else
        call OutputNamedValue ( 'Looking to write DSName', swath )
      endif
      DSName = swath
      allocate(values(nPr, nvTimes))
      allocate(precisions(nPr, nvTimes))
      if ( options%verbose ) &
        & call OutputNamedValue ( 'shp(values)', (/ shape(values) /) )
      call LoadFromHDF5DS ( valuesFile, options%DSNames, values )
      if ( options%verbose ) &
        & call Dump( Values(:,1), 'values used in overwriting', width=5 )
      call LoadFromHDF5DS ( valuesFile, options%PrecNames, precisions )
      if ( options%verbose ) &
        & call Dump( precisions(:,1), 'precisions used in overwriting', width=5 )
      N = 0
!       do i=1, n_filenames
      i = 1
        if ( options%verbose ) &
          & print *, 'Beginning with L2GP file ' // trim(filenames(i))
        numswathssofar = mls_InqSwath ( filenames(i), SwathList, listSize, &
             & hdfVersion=HDFVERSION_5)
!         if ( options%swathNames /= ' ' ) then
!           swathList1 = swathList
!           swathList = Intersection( options%swathNames, swathList1 )
!         endif
        swathList = options%swathNames
        !
        iSwathNum = StringElementNum( swathList, trim(swath), countEmpty )
        call dump ( swathList, 'List of swaths' )
        print *, 'swath to read: ', trim(swath)
        if ( iSwathNum < 1 ) then
          if ( options%debug ) then
            call OutputNamedValue ( 'swathName', swath )
            call OutputNamedValue ( 'swathList', swathList )
            call MLSMessage( MLSMSG_Warning, ModuleName, &
              & "swathName not found in list; so not replaced" )
          endif
          return ! was cycle
        else
          call OutputNamedValue ( 'L2GP File', trim(filenames(i)) )
          call Output( 'Replacing ' // trim(swath) // ' with ' // &
            & trim(DSName) // ' values', advance='yes' )
        endif
        status = InitializeMLSFile( L2GPFile, type=l_swath, access=DFACC_RDWR, &
            & content='L2GP', name=trim(filenames(i)), hdfVersion=HDFVERSION_5 )
        if ( options%verbose ) print *, 'Reading from: ', trim(filenames(i))
        call ReadL2GPData( L2GPFile, swath, L2GP, numProfs )
        print *, 'Num profiles: ', numProfs
        if ( options%debug ) call Dump( L2GP )
        if ( options%verbose ) &
          & call Dump( L2GP%L2GPValue(1,1,:), 'L2gp before overwriting', width=5 )
        call insertL2GPDataByOverwrite ( L2GP, values, precisions )
        if ( options%verbose ) &
          & call Dump( L2GP%L2GPValue(1,1,:), 'L2gp after overwriting', width=5 )

!    stop
        if ( .not. options%dryrun ) &
          & call AppendL2GPData ( L2GP, L2GPFile, swath )
        N = L2GP%nTimes
        ! a(1:N) = L2GP%L2GPValue(1,1,1:N)
        if ( options%verbose ) &
          & call dump( L2GP%L2GPValue(1,1,1:N), 'new values', width=5 )
        if ( options%debug ) call Dump( L2GP )
        call DestroyL2GPContents( L2GP )
        if ( options%debug ) then
          call ReadL2GPData( L2GPFile, swath, L2GP, numProfs )
          call output ( 'insert L2GP as reread', advance='yes' )
          call Dump( L2GP )
          call DestroyL2GPContents( L2GP )
        endif
        numTotProfs = numTotProfs + N
        call Dump( L2GPFile )
        call MLS_CloseFile( L2GPFile )
!       enddo
      deallocate(precisions)
      deallocate(values)
      call sayTime('insertting this swath', tFile)
      call Blanks( 72, advance='yes' )
    ! enddo
    call MLS_CloseFile( ValuesFile )
    call MLS_CloseFile( WeightsFile )
    call MLS_CloseFile( DGMFile )
    call sayTime('insertting all swaths')
  end subroutine insert_swaths

  ! (1) overwriting the current l2gpvalue using values
  ! and precisions. We'll attempt to match profiles to the MAFs
  ! and also match pressure levels
  subroutine insertL2GPDataByOverwrite ( L2GP, values, precisions )
    ! Args
    ! Note that we'll also use nearestProfiles and
    ! Output_Pressure_Levels_Indices
    ! Should they be passed as args, too?
    type(L2GPData_T)                            :: L2GP
    real(rgp), dimension(:,:), intent(in)       :: precisions   
    real(rgp), dimension(:,:), intent(in)       :: values  
    ! internal variables
    integer                                     :: k
    integer                                     :: MAF
    integer                                     :: time ! Like chunk number
    integer                                     :: ntimes ! num of values
    ! Executable
    ! ntimes = size(values, 2)
    ntimes = L2GP%nTimesTotal
    print *, 'Shape (values): ', shape(values)
    print *, 'Shape (l2gpvalues): ', shape(L2GP%l2gpValue)
    print *, 'ntimes: ', ntimes
    do time = 1, ntimes
      ! MAF = nearestMAF ( time, nearestProfiles )
      MAF = matchingMAF ( time, GeodLat, GeodAngle, L2GP )
      if ( MAF > 0 ) then
!         do k = 1, size(Output_Pressure_Levels_Indices)
!           L2GP%l2gpValue(1, Output_Pressure_Levels_Indices(k), time) = &
!             & values(k, MAF)
!           L2GP%l2gpPrecision(1, Output_Pressure_Levels_Indices(k), time) = &
!             & precisions(k, MAF)
!         enddo
        do k = 1, size(values, 1)
          if ( values(k, MAF) > -999.00_rgp ) then
            L2GP%l2gpValue(1, k, time) = &
              & values(k, MAF)
            L2GP%l2gpPrecision(1, k, time) = &
              & precisions(k, MAF)
          else
            L2GP%l2gpPrecision(1, k, time) = -1._rgp
          endif
        enddo
      else
        print *, 'No MAF close to profile ', time
      endif
    enddo
  end subroutine insertL2GPDataByOverwrite
  
  function matchingMAF ( profile, GeodLat, GeodAngle, L2GP ) result ( MAF )
    ! Return the nearest MAF to the input profile number
    ! Do we care that the MAF index numbers conventionally start at 0?
    ! Well .. since we'll use the returned integer value as an index in an array
    ! whose lower bound is 1 instead of 0,  then "No".
    
    ! However, once again we see how ill-advised we were when we ever
    ! sought to make MAFs a 0-based array. From that point on, every use
    ! of MAF causes us to question whether we're making an off-by-one
    ! error. 
    
    ! Moreover, if we ever find we've made an off-by-one error, we won't
    ! know whether it's due to the MAFs being 0-based or due to some other
    ! indexing error.
     integer, intent(in)                 :: profile
     real (rgp), dimension(:,:)          :: GeodLat
     real (rgp), dimension(:,:)          :: GeodAngle
     type(L2GPData_T)                    :: L2GP
     integer                             :: MAF
     !
     integer, dimension(1)               :: intarray
     ! Executable
     ! This is the method we use in l2/NeuralNet_m.f90 to find the
     ! MAF matching our profile. See lines 465-469
     intarray = minloc( &
          & abs(GeodAngle(36,:)-L2GP%GeodAngle(profile)) &
          & + &
          & abs(GeodLat(36,:)-L2GP%latitude(profile)) &
          & ) - 1
     MAF = intarray(1) + 1 ! Because we'll use this MAF as an index, start at 1
  end function matchingMAF
  
  function nearestMAF ( profile, nearestProfiles ) result ( MAF )
    ! Return the nearest MAF to the input profile number
    ! Do we care that the MAF index numbers conventionally start at 0?
    ! Well .. since we'll use the returned integer value as an index in an array
    ! whose lower bound is 1 instead of 0,  then "No".
    
    ! However, once again we see how ill-advised we were when we ever
    ! sought to make MAFs a 0-based array. From that point on, every use
    ! of MAF causes us to question whether we're making an off-by-one
    ! error. 
    
    ! Moreover, if we ever find we've made an off-by-one error, we won't
    ! know whether it's due to the MAFs being 0-based or due to some other
    ! indexing error.
     integer, intent(in)                 :: profile
     integer, dimension(:), intent(in)   :: nearestProfiles
     integer                             :: MAF
     ! Executable
     MAF = FindFirst ( nearestProfiles, profile )
  end function nearestMAF

!------------------------- get_filename ---------------------
    subroutine get_filename(filename, n_filenames, options)
    ! Added for command-line processing
     character(len=255), intent(out) :: filename          ! filename
     integer, intent(in)             :: n_filenames
     type ( options_T ), intent(inout) :: options
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
      elseif ( filename(1:4) == '-deb' ) then
        options%debug = .true.
        exit
      elseif ( filename(1:4) == '-dry' ) then
        options%dryrun = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:4) == '-DGM' ) then
        call getarg ( i+1+hp, options%DGMFile )
        i = i + 1
        exit
      else if ( filename(1:4) == '-Lf' ) then
        call getarg ( i+1+hp, options%L1BOAFile )
        i = i + 1
        exit
      else if ( filename(1:4) == '-Vf ' ) then
        call getarg ( i+1+hp, options%newValues )
        i = i + 1
        exit
      else if ( filename(1:4) == '-Wf ' ) then
        call getarg ( i+1+hp, options%Weights )
        i = i + 1
        exit
      else if ( filename(1:3) == '-s ' ) then
        call getarg ( i+1+hp, options%swathNames )
        i = i + 1
        exit
      else if ( filename(1:3) == '-d ' ) then
        call getarg ( i+1+hp, options%DSNames )
        i = i + 1
        exit
      else if ( filename(1:3) == '-p ' ) then
        call getarg ( i+1+hp, options%PrecNames )
        i = i + 1
        exit
!       elseif ( filename(1:3) == '-sp' ) then
!         call igetarg ( i+1+hp, options%spread )
!         i = i + 1
!         exit
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
      print *,  "Enter the name of the L2GP file. " // &
       &  "The default output file name will be used."
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Overwrite the swath values in l2gp files, e.g. nrt dgg files,'
      write (*,*) &
      & 'using values stored in a different file (made by a n-n script)'
      write (*,*) &
      & 'Usage: insertl2gpvalues [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options:'
      write (*,*) ' -dryrun       => dont execute, just describe'
      write (*,*) ' -f filename   => add filename to list of l2gp filenames'
      write (*,*) '                  (can do the same w/o the -f)'
      write (*,*) ' -d ds_name    => get values to insert from ds_name'
      write (*,*) ' -p ds_name    => get precisions to insert from ds_name'
      write (*,*) ' -s swaths     => insert new values into swaths'
!      write (*,*) ' -DGM dgmname  => read nearest profile dgmname'
      write (*,*) ' -Lf name      => name of L1BOA file'
      write (*,*) ' -Vf vname     => read new values from vname'
!      write (*,*) ' -Wf wname     => read Output_Pressure_Levels_Indices from wname'
      write (*,*) ' -debug        => switch on debug mode'
      write (*,*) ' -v            => switch on verbose mode'
      write (*,*) '                   (we must then assume the surfaces match)'
      write (*,*) ' -h            => print brief help'
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
end program insertL2GPValues
!==================

! $Log$
