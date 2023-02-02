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
   use HDF, only: DFACC_Create, Dfacc_RdOnly, Dfacc_Rdwr
   use HDF5, only: HSize_T
   use HighOutput, only: OutputNamedValue
   use Intrinsic, only: L_HDF, L_Swath
   use L2GPData, only: L2GPData_T, L2GPnamelen, MaxSwathNamesBufSize, Rgp, &
     & AppendL2GPData, DestroyL2GPContents, Dump, &
     & ReadL2GPData, SetupNewL2GPRecord
   use Machine, only: Hp, Getarg
   use MLSCommon, only: MLSFile_T
   use MLSFiles, only: HDFVersion_5, &
     & Dump, MLS_CloseFile, MLS_OpenFile, MLS_Exists, &
     & MLS_Inqswath, InitializeMLSFile
   use MLSFinds, only: FindFirst, FindLast
   use MLSHDF5, only: GetAllHDF5DSNames, GetHDF5DSDims, LoadFromHDF5DS, &
     & MLS_H5open, MLS_H5close
   use MLSMessageModule, only: MLSMessageConfig, MLSMessage, MLSMSG_Warning
   use MLSStringLists, only: GetStringElement, &
     & Intersection, NumStringElements, StringElementNum
   use Output_M, only: Blanks, Output
   use PCFHdr, only: GlobalAttributes, GlobalAttributes_T, &
     & DumpGlobalAttributes, He5_ReadGlobalAttr, He5_WriteGlobalAttr
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
! Usage (1): nrt
!   insertl2gpvalues -s H2O -d ANN_Prediction -p ANN_Precision \
!       -Lf L1BOA_File \
!       -Vf New_Values_File \
!       -conv 1 \
!       -qual 0 \
!       -stat 68 \
!       L2GP_File`

! We'll use it to modify existing swaths, giving them new values
! and precisions predicted by a different retrieval, e.g.
! a trained neural network.

! To create a file of predicted values, you must run's Frank's python script.
! A typical way is the following
! /testing/workspace/pwagner/t2--t-negovlaps 179: echo /usr/local/anaconda3/bin/python3 $H2ONNSCRIPT $L1BRADG $L1BRADD $L1BOA $Weights_File $Prediction_File
! /usr/local/anaconda3/bin/python3 /users/fwerner/Documents/software/python/ann/h2o_prediction.py /data/emls/l1b/v05.01/2022/009/MLS-Aura_L1BRADG_v05-01-c01_2022d009.h5 /data/emls/l1b/v05.01/2022/009/MLS-Aura_L1BRADD_v05-01-c01_2022d009.h5 /data/emls/l1b/v05.01/2022/009/MLS-Aura_L1BOA_v05-01-c01_2022d009.h5 /users/fwerner/Documents/database/neural_network_weights/weights/v05-0x/MLS-Aura_ANN-H2O_v05-0x-01_20220415.h5 prediction.h5


! Usage (2): dmp (derived meteorological products
!    insertl2gpvalues -dmp \ 
!       -s H2O-StdProd \
!       -Vf DMP_File \
!       L2GP_File

! rm DMP.he5 DMP.out ;
! /users/pwagner/mlspgs/bin/IFC.Linux.ifc19-OL8/insertl2gpvalues -v \
!  -Af MLS-Aura_L2GP-CO_v05-01-c01_2022d325.he5 \
!  -Vf MLS-Aura_L2EDMP-GEOS5124-v201_v05-01-c01_2023d008.he5 -d Altitude \
!  -dmp DMP.he5 > DMP.out

! We'll use it to convert the older DMP files to genuine hdfeos files
! that can be read by users with, e.g., python-based tools.
  type options_T
    logical               :: debug   = .false.
    logical               :: silent  = .false.
    logical               :: verbose = .false.
    logical               :: dryrun  = .false.
    logical               :: DMP     = .false.
    character(len=255)    :: swathNames = ' '        ! which swath to modify
    character(len=255)    :: DSNames    = ' '        ! which new DS holds values
    character(len=255)    :: PrecNames  = ' '        ! which new DS holds precs
    character(len=255)    :: newValues  = ' '        ! file with new values
    character(len=255)    :: L1BOAFile  = ' '        ! L1BOA file
    character(len=255)    :: L2AttrFile = ' '        ! L2GP file with attributes
    character(len=255)    :: Weights    = ' '        ! Weights file
    integer               :: NewStatus  = -999       ! Replace value for status
    double precision      :: NewConverg = -999.99d0  ! for convergence
    double precision      :: NewQuality = -999.99d0  ! for Quality
    double precision      :: PressBottom= -999.99d0  ! Bottommost Pr overwritten
    double precision      :: PressTop   = -999.99d0  ! Topmost Pr overwritten
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
  ! DSName and swath must be long enough to hold the full path
  ! set by everything beginning with "/HDFEOS/.."s
  character(len=255)                  :: DSName
  character(len=MAXSWATHNAMESBUFSIZE) :: DSList
  character(len=255)                  :: swath
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList
  character(len=MAXSWATHNAMESBUFSIZE) :: swathList1
  integer :: ACCESS
  integer :: listSize
  integer :: NUMSWATHSSOFAR
  integer, allocatable, dimension(:)  :: nearestProfiles
  integer, allocatable, dimension(:)  :: Output_Pressure_Levels_Indices
  real (rgp), dimension(:,:), allocatable :: GeodAngle
  real (rgp), dimension(:,:), allocatable :: GeodLat
  real (rgp), dimension(:,:), allocatable :: Pressure
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
     if ( options%DMP )then
       ! n_filenames = n_filenames + 1
       ! filenames(n_filenames) = filename
       ! cycle
       ! We need a Fortran no op here
       listSize = 1 ! This ought to harmless enough
     elseif ( mls_exists(trim(filename)) /= 0 ) then
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
    if ( options%silent ) then
      options%debug   = .false.
      options%verbose = .false.
    endif

    numswathssofar = 0
    if ( options%verbose ) &
      & call Dump ( filenames(1:n_filenames), 'l2gp files', width=1 )
    if ( options%DMP ) then
      if ( mls_exists(trim(filenames(1))) /= 0 ) then
        if ( .not. options%silent ) print *, 'Must create: ', trim(filenames(1))
        ACCESS = DFACC_CREATE
      else
        if ( .not. options%silent ) print *, 'Must add to: ', trim(filenames(1))
        ACCESS = DFACC_RDWR
      endif
      call DMP_swaths ( ACCESS )
    else
      call insert_swaths
    endif
  endif
  call mls_h5close(error)
contains

!------------------------- DMP_swaths ---------------------
! Identify and insert values from DMP products

  subroutine DMP_swaths ( ACCESS )
    ! Dummy args
    integer, intent (in)            :: ACCESS
    ! Internal variables
    integer :: i
    integer :: iSwathNum
    ! integer :: jj
    type(L2GPData_T)  :: L2GP
    type( MLSFile_T ) :: L2GPFile
    type( MLSFile_T ) :: valuesFile
    type( MLSFile_T ) :: AttributesFile  ! If needed to copy global attributes from
    integer :: nPr
    integer :: nPrs
    integer :: numProfs
    integer :: numTotProfs
    integer :: NvTimes
    real (rgp), dimension(:), allocatable   :: Pressure
    real (rgp), dimension(:,:), allocatable :: precisions
    real (rgp), dimension(:,:), allocatable :: values
    integer(kind=hSize_t), dimension(3)     :: DIMS
    type (GlobalAttributes_T)               :: gAttributes
    character(len=*), parameter             :: identifier_product_doi = &
      & "10.5067/AURA/MLS/DATA2521"
    character(len=16)                       :: ProcessLevel
    integer                                 :: DayofYear
    double precision                        :: TAI93At0zOfGranule
    character(len=64)                       :: HostName
    character(len=256)                      :: MiscNotes
    
    ! --------------------------------------------------------------------
    ! This is the number of instances in the new values file
    ! integer :: nvTimes
    ! --------------------------------------------------------------------
    ! These are the number of pressure levels each has been given
    ! integer, dimension(nDS), parameter :: NPrs = (/ 41, 14, 12, 13, 37, 41 /)
    
    ! Note that in the latest netCDF file, all species are on 37 pressure
    ! surfaces except for the _HR species. This choice happens to match
    ! the resolution of our l2gp products. Will future netCDF files from
    ! Joaquim adhere to this convention, too?
    ! integer, dimension(1), parameter :: NPrs = (/ 42 /)
    ! --------------------------------------------------------------------
    ! The DMP use different pressure levels from any other standard species
! HDF5 "/data/emls/dmp/l2edmp/v05.01/2023/008/MLS-Aura_L2EDMP-GEOS5124-v201_v05-01-c01_2023d008.he5" {
! DATASET "/HDFEOS/GEOLOCATION/Pressure" {
!    DATATYPE  H5T_IEEE_F32LE
!    DATASPACE  SIMPLE { ( 61 ) / ( 61 ) }
!    DATA {
!    (0): 1000, 825.404, 681.292, 562.341, 464.159, 383.119, 316.228, 261.016,
!    (8): 215.443, 177.828, 146.78, 121.153, 100, 82.5404, 68.1292, 56.2341,
!    (16): 46.4159, 38.3119, 31.6228, 26.1016, 21.5443, 17.7828, 14.678,
!    (23): 12.1153, 10, 8.25404, 6.81292, 5.62341, 4.64159, 3.83119, 3.16228,
!    (31): 2.61016, 2.15443, 1.77828, 1.4678, 1.21153, 1, 0.8254, 0.68129,
!    (39): 0.56234, 0.46416, 0.38312, 0.31623, 0.26102, 0.21544, 0.17783,
!    (46): 0.14678, 0.12115, 0.1, 0.08254, 0.06813, 0.05623, 0.04642, 0.03831,
!    (54): 0.03162, 0.0261, 0.02154, 0.01778, 0.01468, 0.01212, 0.01
!    }
!    ATTRIBUTE "MissingValue" {
!       DATATYPE  H5T_IEEE_F32LE
!       DATASPACE  SCALAR
!       DATA {
!       (0): 1e+15
!       }
!    }
!    ATTRIBUTE "Title" {
!       DATATYPE  H5T_STRING {
!          STRSIZE 9;
!          STRPAD H5T_STR_NULLTERM;
!          CSET H5T_CSET_ASCII;
!          CTYPE H5T_C_S1;
!       }
!       DATASPACE  SCALAR
!       DATA {
!       (0): "Pressure"
!       }
!    }
!    ATTRIBUTE "UniqueFieldDefinition" {
!       DATATYPE  H5T_STRING {
!          STRSIZE 13;
!          STRPAD H5T_STR_NULLTERM;
!          CSET H5T_CSET_ASCII;
!          CTYPE H5T_C_S1;
!       }
!       DATASPACE  SCALAR
!       DATA {
!       (0): "MLS-Specific"
!       }
!    }
!    ATTRIBUTE "Units" {
!       DATATYPE  H5T_STRING {
!          STRSIZE 4;
!          STRPAD H5T_STR_NULLTERM;
!          CSET H5T_CSET_ASCII;
!          CTYPE H5T_C_S1;
!       }
!       DATASPACE  SCALAR
!       DATA {
!       (0): "hPa"
!       }
!    }
!    ATTRIBUTE "_FillValue" {
!       DATATYPE  H5T_IEEE_F32LE
!       DATASPACE  SCALAR
!       DATA {
!       (0): 1e+15
!       }
!    }
! }
! }

    ! --------------------------------------------------------------------
    ! Executable
    ! ------------------------------------------------------------------------
    ! 1st data sets to read:    /HDFEOS/GEOLOCATION/Pressure
    ! read from                        Values File
    status = InitializeMLSFile( ValuesFile, type=l_hdf, access=DFACC_RDOnly, &
      & content='hdf', name=trim(options%newValues), hdfVersion=HDFVERSION_5 )
    ! call GetAllHDF5DSNames ( ValuesFile, DSList )
    ! if ( options%verbose ) then
    !   print *, 'DS names in values file ' // trim(DSList)
    ! endif
    DSName = '/HDFEOS/GEOLOCATION/Pressure'
    call GetHDF5DSDims ( ValuesFile, DSName, DIMS )
    NPrs = DIMS(1)
    if ( options%verbose ) print *, 'NPrs: ', NPrs
    allocate(Pressure(NPrs))
    call LoadFromHDF5DS ( ValuesFile, DSName, Pressure )
    if ( options%verbose ) call dump( Pressure, 'Pressure' )

    ! The typicaal l2gp values are stored in a roundabout way something like this:
!   ~/mlspgs/util 212: h5dump -g /HDFEOS/SWATHS/Altitude/"Data Fields" /data/emls/dmp/l2edmp/v05.01/2023/008/* | grep Altitude
! GROUP "/HDFEOS/SWATHS/Altitude/Data Fields" {
!    DATASET "Altitude" {
!          (0): "Altitude"
!       HARDLINK "/HDFEOS/SWATHS/Altitude/Data Fields/Altitude"
  
    DSName = '/HDFEOS/SWATHS/' // trim(options%DSNames) // &
      & '/Data Fields/' // trim(options%DSNames)
    call GetHDF5DSDims ( valuesFile, DSName, DIMS )
    swath = options%DSNames
    NPr     = DIMS(1)
    NvTimes = DIMS(2)
    if ( options%verbose ) print *, 'NPr, NvTimes: ', DIMS
    if ( options%verbose ) print *, 'Num swaths: ', NumStringElements( DSList, countEmpty )
    if ( NvTimes < 1 ) then
      NPr = 1
      NvTimes = Dims(1)
    endif
    numTotProfs = 0
    ! jj = 1
    ! call GetStringElement( DSList, DSName, jj, countEmpty )
    call time_now ( tFile )
    ! swath = options%swathNames
    if ( options%verbose ) print *, 'inserting values into swath: ', trim(swath)
    if ( options%verbose ) print *, 'using values from DS: ', trim(swath)
    if ( .not. options%silent ) call Blanks( 72, fillChar='-', advance='yes' )
    if ( options%verbose ) then
      call OutputNamedValue ( 'DSName', DSName )
    elseif ( .not. options%silent ) then
      call OutputNamedValue ( 'Looking to write DSName', DSName )
    endif
    ! DSName = swath
    allocate(values(nPr, nvTimes))
    allocate(precisions(nPr, nvTimes))
    if ( options%verbose ) &
      & call OutputNamedValue ( 'shp(values)', (/ shape(values) /) )
    if ( .not. options%silent ) print *, 'Trying to read dataset: ', trim(DSName)
    if ( NPr > 1 ) then
      call LoadFromHDF5DS ( valuesFile, DSName, values )
    else
      call LoadFromHDF5DS ( valuesFile, DSName, values(1,:) )
    endif
    if ( options%verbose ) &
      & call Dump( Values(:,1), 'values used in overwriting', width=5 )
    ! call LoadFromHDF5DS ( valuesFile, options%PrecNames, precisions )
    precisions = 0.
    if ( options%verbose ) &
      & call Dump( precisions(:,1), 'precisions used in overwriting', width=5 )
    i = 1
    if ( options%verbose ) &
      & print *, 'Beginning with L2GP file ' // trim(filenames(i))

    status = InitializeMLSFile( L2GPFile, type=l_swath, access=ACCESS, &
        & content='L2GP', name=trim(filenames(i)), hdfVersion=HDFVERSION_5 )
    ! if ( options%verbose ) print *, 'Reading from: ', trim(filenames(i))
    ! call ReadL2GPData( L2GPFile, swath, L2GP, numProfs )
    call SetupNewL2GPRecord ( L2GP, nFreqs=1, NLevels=NPrs, NTimes=NvTimes )
    call LoadFromHDF5DS ( ValuesFile, '/HDFEOS/GEOLOCATION/Latitude          ', L2GP%Latitude          )
    call LoadFromHDF5DS ( ValuesFile, '/HDFEOS/GEOLOCATION/LineOfSightAngle  ', L2GP%losAngle  )
    call LoadFromHDF5DS ( ValuesFile, '/HDFEOS/GEOLOCATION/LocalSolarTime    ', L2GP%SolarTime    )
    call LoadFromHDF5DS ( ValuesFile, '/HDFEOS/GEOLOCATION/Longitude         ', L2GP%Longitude         )
    call LoadFromHDF5DS ( ValuesFile, '/HDFEOS/GEOLOCATION/OrbitGeodeticAngle', L2GP%geodAngle)
    if ( NPr > 1 ) &
  & call LoadFromHDF5DS ( ValuesFile, '/HDFEOS/GEOLOCATION/Pressure          ', L2GP%Pressures         )
    call LoadFromHDF5DS ( ValuesFile, '/HDFEOS/GEOLOCATION/SolarZenithAngle  ', L2GP%SolarZenith  )
    call LoadFromHDF5DS ( ValuesFile, '/HDFEOS/GEOLOCATION/Time              ', L2GP%Time              )
    ! call LoadFromHDF5DS ( ValuesFile, '/HDFEOS/GEOLOCATION/seconds_into_day  ', L2GP%seconds_into_day  )
    L2GP%NFreqs       = 1
    L2GP%NLevels      = NPrs
    L2GP%NTimes       = NvTimes
    L2GP%Name         = options%DSNames
    L2GP%MissingL2GP  = real(1.E15, kind=rgp)
    ! L2GP%MissingValue = 1.E15
    
    ! Quality flagging not applicable for DMP
    L2GP%Status       = 0
    L2GP%Quality      = 0
    L2GP%Convergence  = 0
    L2GP%ChunkNumber  = 0
    
    if ( NPr < 2 ) then
      deallocate ( L2GP%L2GPValue )
      deallocate ( L2GP%L2GPPrecision )
      L2GP%NLevels      = 1
      allocate ( L2GP%L2GPValue(1, 1, NvTimes) )
      allocate ( L2GP%L2GPPrecision(1, 1, NvTimes) )
    endif

    if ( options%debug ) call Dump( L2GP )
    ! call insertL2GPDataByOverwrite ( L2GP, values, precisions )
    L2GP%L2GPValue(1,:,:) = values
    if ( options%debug ) print *, 'shape of l2gp values: ', shape(L2GP%L2GPValue)
    L2GP%L2GPPrecision(1,:,:) = precisions
    L2GP%Pressures = Pressure
    if ( options%verbose ) &
      & call Dump( L2GP%L2GPValue(1,1,:), 'L2gp after overwriting', width=5 )

!    stop
    if ( len_trim(options%L2AttrFile) > 1 .and. .not. options%dryrun ) then
      status = InitializeMLSFile( AttributesFile, type=l_swath, access=DFACC_RDOnly, &
        & content='L2GP', name=trim(options%L2AttrFile), hdfVersion=HDFVERSION_5 )
      call MLS_OpenFile( AttributesFile )
      if ( options%verbose ) print *, 'Reading from: ', trim(options%L2AttrFile)
      if ( .not. options%silent ) call Dump ( AttributesFile )
      call he5_readglobalattr ( AttributesFile%fileID%f_ID, gAttributes, &
        & ProcessLevel, DayofYear, TAI93At0zOfGranule, &
        & HostName, MiscNotes )
      gAttributes%doi                = identifier_product_doi
      gAttributes%ProcessLevel       = 'L2'
      gAttributes%DayofYear          = DayofYear
      gAttributes%TAI93At0zOfGranule = TAI93At0zOfGranule
      gAttributes%HostName           = HostName
      gAttributes%MiscNotes          = '(Not relevant)'
      gAttributes%ProductionLoc      = 'SCF'
      GlobalAttributes               = gAttributes
      if ( .not. options%silent ) call DumpGlobalAttributes
      call he5_writeglobalattr ( L2GPFile, DOI=.true., &
        & skip_if_already_there=.false. )
      call MLS_CloseFile( AttributesFile )
    endif
    if ( .not. options%dryrun ) &
      & call AppendL2GPData ( L2GP, L2GPFile, L2GP%Name )
    if ( options%verbose ) &
      & call dump( L2GP%L2GPValue(1,1,1:NvTimes), 'new values', width=5 )
    if ( options%debug ) call Dump( L2GP )
    call DestroyL2GPContents( L2GP )
    if ( options%debug ) then
      call ReadL2GPData( L2GPFile, swath, L2GP, numProfs )
      call output ( 'insert L2GP as reread', advance='yes' )
      call Dump( L2GP )
      call DestroyL2GPContents( L2GP )
    endif
    numTotProfs = numTotProfs + NvTimes
    if ( .not. options%silent ) call Dump( L2GPFile )
    if ( L2GPFile%StillOpen ) call MLS_CloseFile( L2GPFile )
    deallocate(precisions)
    deallocate(values)
    deallocate(Pressure)
    call sayTime('inserting this swath', tFile)
    if ( .not. options%silent ) call Blanks( 72, advance='yes' )
    call MLS_CloseFile( ValuesFile )
    call sayTime('inserting all DMP swaths')
  end subroutine DMP_swaths

!------------------------- insert_swaths ---------------------
! Identify and insert values derived from nrt python scripts

  subroutine insert_swaths
    ! Internal variables
    integer :: i
    integer :: iSwathNum
    integer :: jj
    type(L2GPData_T) :: L2GP
    type( MLSFile_T ) :: L2GPFile
    type( MLSFile_T ) :: L1BOAFile
    type( MLSFile_T ) :: valuesFile
    type( MLSFile_T ) :: WeightsFile
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
    ! What we'll do instead is match MAFs and profiles using a method
    ! similar to what we did in mlsl2

    ! ------------------------------------------------------------------------
    ! 2nd data set to read:          Output_Pressure_Levels_Indices
    ! read from                        Weights File
    
    ! We don't need these if the prediction and precision arrays
    ! already fill in the vertical levels missing from the Weights file.

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
    jj = 1
    call GetStringElement( DSList, DSName, jj, countEmpty )
    call time_now ( tFile )
    swath = options%swathNames
    if ( options%verbose ) print *, 'inserting values into swath: ', trim(swath)
    if ( options%verbose ) print *, 'using values from DS: ', trim(swath)
    if ( .not. options%silent ) call Blanks( 72, fillChar='-', advance='yes' )
    if ( options%verbose ) then
      call OutputNamedValue ( 'DSName', swath )
    elseif ( .not. options%silent ) then
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
    i = 1
    if ( options%verbose ) &
      & print *, 'Beginning with L2GP file ' // trim(filenames(i))
    numswathssofar = mls_InqSwath ( filenames(i), SwathList, listSize, &
         & hdfVersion=HDFVERSION_5)
    swathList = options%swathNames
    !
    iSwathNum = StringElementNum( swathList, trim(swath), countEmpty )
    if ( .not. options%silent ) call dump ( swathList, 'List of swaths' )
    if ( .not. options%silent ) print *, 'swath to read: ', trim(swath)
    if ( iSwathNum < 1 ) then
      if ( options%debug ) then
        call OutputNamedValue ( 'swathName', swath )
        call OutputNamedValue ( 'swathList', swathList )
        call MLSMessage( MLSMSG_Warning, ModuleName, &
          & "swathName not found in list; so not replaced" )
      endif
      return ! was cycle
    elseif ( .not. options%silent ) then
      call OutputNamedValue ( 'L2GP File', trim(filenames(i)) )
      call Output( 'Replacing ' // trim(swath) // ' with ' // &
        & trim(DSName) // ' values', advance='yes' )
    endif
    status = InitializeMLSFile( L2GPFile, type=l_swath, access=DFACC_RDWR, &
        & content='L2GP', name=trim(filenames(i)), hdfVersion=HDFVERSION_5 )
    if ( options%verbose ) print *, 'Reading from: ', trim(filenames(i))
    call ReadL2GPData( L2GPFile, swath, L2GP, numProfs )
    if ( .not. options%silent ) print *, 'Num profiles: ', numProfs
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
    if ( .not. options%silent ) call Dump( L2GPFile )
    if ( L2GPFile%StillOpen ) call MLS_CloseFile( L2GPFile )
    deallocate(precisions)
    deallocate(values)
    call sayTime('inserting this swath', tFile)
    if ( .not. options%silent ) call Blanks( 72, advance='yes' )
    call MLS_CloseFile( ValuesFile )
    call MLS_CloseFile( WeightsFile )
    call sayTime('inserting all swaths')
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
    integer                                     :: k1 ! PressBottom level
    integer                                     :: k2 ! PressTop level
    integer                                     :: MAF
    integer                                     :: time ! Like chunk number
    integer                                     :: ntimes ! num of values
    ! Executable
    ! ntimes = size(values, 2)
    ntimes = L2GP%nTimesTotal
    if ( .not. options%silent ) print *, 'Shape (values): ', shape(values)
    if ( .not. options%silent ) print *, 'Shape (l2gpvalues): ', shape(L2GP%l2gpValue)
    if ( .not. options%silent ) print *, 'ntimes: ', ntimes
    k1 = 1
    k2 = size(values, 1)
    if ( options%PressBottom > -999.99d0 ) &
      & k1 = FindFirst( (options%PressBottom - L2GP%Pressures) >= 0.d0 )
    if ( options%PressTop > -999.99d0 ) &
      & k2 = FindLast( (options%PressTop - L2GP%Pressures) <= 0.d0 )
    ! FindFirst and FindLast may return 0
    if ( k1 < 1 ) k1 = 1
    if ( k2 < 1 ) k2 = size(values, 1)
    if ( options%verbose ) then
      if ( k1 > 1 .or. k2 < size(values, 1) ) &
        call outputNamedValue ( 'k1, k2', (/ k1, k2 /) )
    endif
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
!        do k = 1, size(values, 1)
        do k = k1, k2
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
        if ( .not. options%silent ) print *, 'No MAF close to profile ', time
      endif
    enddo
    ! Have we been asked to overwrite Convergence, Quality, Status fields?
    if ( options%NewConverg > 0. ) then
      L2GP%Convergence = options%NewConverg
    endif
    if ( options%NewStatus > 0. ) then
      L2GP%Status = options%NewStatus
    endif
    if ( options%NewQuality > 0. ) then
      L2GP%Quality = options%NewQuality
    endif
  end subroutine insertL2GPDataByOverwrite
  
  ! This function attempts to replicate how we matched MAFs to profiles in mlsl2
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
      elseif ( filename(1:4) == '-dmp' ) then
        options%DMP = .true.
        exit
      elseif ( filename(1:4) == '-dry' ) then
        options%dryrun = .true.
        exit
      elseif ( filename(1:7) == '-silent' ) then
        options%silent = .true.
        exit
      elseif ( filename(1:3) == '-v ' ) then
        options%verbose = .true.
        exit
      else if ( filename(1:3) == '-f ' ) then
        call getarg ( i+1+hp, filename )
        i = i + 1
        exit
      else if ( filename(1:4) == '-Af' ) then
        call getarg ( i+1+hp, options%L2AttrFile )
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
      elseif ( filename(1:5) == '-conv' ) then
        call dgetarg ( i+1+hp, options%NewConverg )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-qual' ) then
        call dgetarg ( i+1+hp, options%NewQuality )
        i = i + 1
        exit
      elseif ( filename(1:5) == '-stat' ) then
        call igetarg ( i+1+hp, options%NewStatus )
        i = i + 1
        exit
      elseif ( filename(1:6) == '-prbot' ) then
        call dgetarg ( i+1+hp, options%PressBottom )
        i = i + 1
        exit
      elseif ( filename(1:6) == '-prtop' ) then
        call dgetarg ( i+1+hp, options%PressTop )
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
      print *,  "Enter the name of the L2GP file. " // &
       &  "The default output file name will be used."
      read(*,'(a)') filename
    endif
    
  end subroutine get_filename
!------------------------- print_help ---------------------
  subroutine print_help
  ! Print brief but helpful message
      write (*,*) &
      & 'Either:'
      write (*,*) &
      & '(1) Overwrite the swath values in l2gp files, e.g. nrt dgg files,'
      write (*,*) &
      & 'using values stored in a different file (made by a n-n script)'
      write (*,*) '                   (we must then assume the surfaces match)'
      write (*,*) &
      & '(2) Convert an older DMP file to genuine hdfeos format'
      write (*,*) &
      & '    Defaults to usage (1)'
      write (*,*) &
      & 'Usage: insertl2gpvalues [options] [filenames]'
      write (*,*) &
      & ' If no filenames supplied, you will be prompted to supply one'
      write (*,*) ' Options:'
      write (*,*) ' in the following (1) means relevant for usage (1)'
      write (*,*) ' in the following (2) means relevant for usage (2)'
      write (*,*) ' -dryrun       => dont execute, just describe'
      write (*,*) ' -f filename   => add filename to list of l2gp filenames'
      write (*,*) '                  (can do the same w/o the -f)'
      write (*,*) ' -d ds_name    => get values to insert from ds_name'
      write (*,*) ' -p ds_name    => get precisions to insert from ds_name'
      write (*,*) ' -s swaths     => insert new values into swaths'
      write (*,*) ' -Af name      => name of file to copy glob. attr. from (2)'
      write (*,*) ' -Lf name      => name of L1BOA file (1)'
      write (*,*) ' -Vf vname     => read new values from file vname'
      write (*,*) ' -debug        => switch on debug mode'
      write (*,*) ' -v            => switch on verbose mode'
      write (*,*) ' -silent       => switch on silent mode'
      write (*,*) ' -dmp          => convert older dmp file to genuine hefeos (2)'
      write (*,*) ' -conv val     => set values of Convergence field to val'
      write (*,*) ' -qual val     => set values of Quality field to val'
      write (*,*) ' -stat val     => set values of Status field to val'
      write (*,*) ' -prbot val    => overwrite only heights at val and above'
      write (*,*) ' -prtop val    => overwrite only heights at val and below'
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
    if ( options%silent ) return
    call output ( "Timing for " // what // " = " )
    call output ( dble(t2 - myt1), advance = 'yes' )
  end subroutine SayTime
!------------------------- dgetarg ---------------------
  subroutine dgetarg ( pos, darg )
   integer, intent(in) :: pos
   double precision, intent(out) :: darg
   character(len=16) :: arg
   call getarg ( pos, arg )
   read(arg, *) darg
  end subroutine dgetarg
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
! Revision 1.4  2022/07/13 20:46:40  pwagner
! Added new cmdline options to restrict height range
!
! Revision 1.3  2022/07/08 20:34:10  pwagner
! 3 new cmdline args to set conv, qual, status
!
! Revision 1.2  2022/05/13 18:02:27  pwagner
! Numerous bugfixes; removed DGM file
!
! Revision 1.1  2022/05/13 17:59:20  pwagner
! First commit
!
