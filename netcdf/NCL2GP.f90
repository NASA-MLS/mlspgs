! Copyright 2020, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module NCL2GP

  ! This module contains routines to
  ! output and input L2GP std products in NetCDF4 format.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dump_1, only: Dump
  use NetCDF ! Everything
  use HDF, only: Dfacc_Rdonly, Dfacc_Read, Dfacc_Rdwr
  use HighOutput, only: BeVerbose, OutputNamedValue
  use Intrinsic ! "units" Type Literals, Beginning With L
  use L2GPData, only: CharAttrLen, Col_Species_Keys, Col_Species_Hash, &
    & Data_Field1, Data_Field2, Dim_Name1, Dim_Name12, Dim_Name123, &
    & L2GPData_T, L2GPNameLen, MaxChunkTimes, MaxNLevels, &
    & Max_DimL, Max_DimL1, Max_DimL12, Max_DimL123, NumGeolocFields, &
    & RGP, &
    & DestroyL2GPContents, Dump, SetupNewL2GPRecord, &
    & WriteMastersFileAttributes
  use MLSCommon, only: L2MetaData_T, &
    & R4, R8, MLSFile_T, UndefinedIntegerValue
  use MLSFiles, only: Dump, InitializeMLSFile, MLS_CloseFile, MLS_OpenFile
  use MLSMessagemodule, only: MLSMSG_Error, MLSMSG_Warning, MLSMessage
  use MLSStats1, only: MLSMin
  use MLSStringlists, only: ExtractSubstring, GetHashElement, GetStringElement, &
    & List2Array, NumStringElements, ReplaceSubstring
  use MLSStrings, only: Lowercase
  use Output_M, only: Output
  use PCfhdr, only: GlobalAttributes_T, GlobalAttributes, &
    & DumpGlobalAttributes
  use SDPToolKit, only: Max_Orbits
  use Time_M, only: SayTime, ConfigureSayTime, Time_Now
  use Trace_M, only: Trace_Begin, Trace_End

  implicit none
  private

  public :: AppendNCL2GPData, CpNCL2GPData, DumpNCL2GP, IsNCL2GPInFile, &
    & ReadNCL2GPData, ReadNCGlobalAttr, &
    & WriteNCFileAttr, WriteNCGlobalAttr, WriteNCL2GPData

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface cpNCL2GPData
    module procedure CpNCL2GPData_fileID
    ! module procedure CpNCL2GPData_fileName
    ! module procedure CpNCL2GPData_MLSFile
  end interface
  
  interface DumpNCL2GP
    ! module procedure DumpNCL2GPData_fileID
    ! module procedure DumpNCL2GPData_fileName
    module procedure DumpNCL2GPData_MLSFile
  end interface
  
  interface ReadNCGlobalAttr
    module procedure ReadNCGlobalAttr_FileID
  end interface

  interface ReadNCL2GPData
    module procedure ReadNCL2GPData_fileID
    module procedure ReadNCL2GPData_fileName
    module procedure ReadNCL2GPData_MLSFile
  end interface
  
  interface WriteNCGlobalAttr
    module procedure WriteNCGlobalAttr_FileID
  end interface

  interface WriteNCL2GPData
    module procedure WriteNCL2GPData_fileID
    module procedure WriteNCL2GPData_MLSFile
  end interface
! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! AppendNCL2GPData  Appends or inserts new L2GPData into a file
! CpNCL2GPData      Copies L2GPData from one file to another
! IsNCGlAttrInFile  Is the global attribute already in the file?
! IsNCL2GPInFile    Is the L2GP already in the file?
! ReadNCL2GPData    Reads an L2GPData from a file
! WriteNCFileAttr   Writes the file attributes at file level
! WriteNCGlobalAttr Writes the global attributes at file level
! WriteNCL2GPData   Writes a new L2GPData to a file
! WriteNCMastrAttr  Rewrites the global attributes based on the master
! === (end of toc) ===

! === (start of api) ===
! AppendNCL2GPData ( L2GPData_T l2gp, MLSFile_T L2GPFile, char* swathName, &
!     [int offset], [int LastProfile], [int TotalNumProfs], [log CreateSwath], &
!     [int MaxChunkSize] )
! === (end of api) ===

  ! Print debugging stuff?
  logical, parameter          :: DEEBUG = .false.  
  real                        :: t2
  integer, parameter          :: MaxFlds = 24
  integer, parameter          :: MaxRank = 7

contains ! ======================= Public Procedures =========================

  ! ---------------------- AppendNCL2GPData  ---------------------------

  subroutine AppendNCL2GPData( l2gp, l2gpFile, &
    & swathName, offset, lastProfile, TotNumProfs, &
    & createSwath, maxchunksize )
    ! sticks l2gp into the swath swathName in the file pointed at by
    ! l2FileHandle,starting at the profile number "offset" (First profile
    ! in the file has offset==0). If this runs off the end of the swath, 
    ! it is lengthened automagically. 
    ! This call has been altered recently, so that it can be used to create
    ! a swath as well as adding to one. 

    ! Arguments

    type(MLSFile_T)                :: L2GPFile

    ! This is a L2GPData_T structure containing all the data to be written
    type (L2GPData_T), intent(INOUT) :: l2gp
    ! This is the name the swath is given in the file. By default it is
    ! the name contained in l2gp
    character (len=*), optional, intent(in) ::swathName!default->l2gp%swathName
    ! This (offset) is the point in the swath at which the data is written. 
    ! First profile in the file has offset==0. If the swath in the file is 
    ! shorter than offset + ( num of profiles in l2gp) then it grows by magic
    integer, intent(in), optional::offset
    ! TotNumProfs is a new argument. It seems only to be used if we are 
    ! creating a swath, rather than adding to one. In that case I guess
    ! it is the total number of profiles in the swath created. I also 
    ! guess that this is done so that we can avoid growing and re-growing 
    ! the swath.
    integer, intent(in), optional ::TotNumProfs
    integer, intent(in), optional ::lastProfile
    logical, intent(in), optional :: createSwath
    integer, optional, intent(in) :: maxchunksize
    ! Local
    logical :: alreadyOpen
    type (L2GPData_T) :: largerl2gp
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: myLastProfile
    character (len=L2GPNameLen) :: myswathName
    logical :: notUnlimited
    integer :: status
    logical :: swath_exists
    logical :: timing
    real :: tFile ! How long have we been fooling with this file?
    ! logical, parameter :: DEEBUG = .false.

    ! Executable code
    call trace_begin ( me, 'AppendNCL2GPData', cond=.false. )
    ! call Dump ( l2gp%chunkNumber, 'Appending l2gp%chunkNumber' )
    call time_now ( tFile )
    call configureSayTime ( tFile )
    status = 0
    timing = DEEBUG .or. BeVerbose( 'l2gp', 2 )

    if (present(lastProfile)) then
      myLastProfile = lastProfile 
    elseif (present(TotNumProfs)) then
      myLastProfile = TotNumProfs 
    else
      myLastProfile = L2GP%nTimesTotal
    endif
    alreadyOpen = L2GPFile%stillOpen
    if ( .not. alreadyOpen ) then
      if ( DEEBUG ) print *, 'Needed to open file ', trim(L2GPFile%name)
      call mls_openFile(L2GPFile, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l2gp file', MLSFile=L2GPFile)
    endif
    if ( timing ) then
      call sayTime( 'Opening file ' // trim(L2GPFile%name), tFile, t2 )
      tFile = t2
    endif
    if ( L2GPFile%access == DFACC_RDONLY )  &
      & call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'l2gp file is rdonly', MLSFile=L2GPFile)
    myswathName = l2gp%name
    if ( present(swathName) ) myswathName = swathName
    
    swath_exists = IsNCL2GPInFile( L2GPFile, myswathName )

    if ( timing ) then
      call sayTime( 'Checking that swath exists', tFile, t2 )
      tFile = t2
    endif
    if ( swath_exists ) then
      if(DEEBUG) print *, 'OK, swath already exists, but HDFEOS bug no longer scares us'
    else
      ! Must create swath in file w/o disturbing other swaths
      if(DEEBUG) print *, 'Must create swath'
      if(DEEBUG) print *, 'Will have ', myLastProfile, ' profiles'
      if(DEEBUG) print *, 'instead of ', l2gp%nTimes, ' profiles'
      ! By default allow limited; 
      ! may force unlimited by setting avoidUnlimitedDims to FALSE
      notUnlimited = .false. ! ( avoidUnlimitedDims .and. present(totNumProfs) )
      if ( present(maxchunksize) ) then
        if ( maxchunksize > l2gp%nTimes ) then
          call SetupNewL2GPRecord ( largerl2gp, proto=l2gp, &
            & nTimes=maxchunksize )
          call OutputNCL2GP_createFile_MF ( largerl2gp, L2GPFile, &
            & myswathName, notUnlimited=notUnlimited )
          call DestroyL2GPContents ( largerl2gp )
        else
          call OutputNCL2GP_createFile_MF ( l2gp, L2GPFile, &
            & myswathName, notUnlimited=notUnlimited )
        endif
      else
        call OutputNCL2GP_createFile_MF ( l2gp, L2GPFile, &
          & myswathName, notUnlimited=notUnlimited )
      endif
    endif

    if(DEEBUG) then
      if ( present(offset) ) then
        print*,"offset=",offset,"myLastProfile=",myLastProfile,&
        "size(l2gp%l2gpValue,3)=",size(l2gp%l2gpValue,3)
      else
        print*,"no offset; myLastProfile=",myLastProfile,&
        "size(l2gp%l2gpValue,3)=",size(l2gp%l2gpValue,3)
      endif
    endif
    if ( size(l2gp%l2gpValue,3) == 0 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & "No profiles in this chunk", MLSFile=L2GPFile )

    else
      ! actual_ntimes = l2gp%nTimes
      ! l2gp%nTimes = max(myLastProfile - offset + 1, 1)
      if ( all(l2gp%chunkNumber == -999) ) &
        & call output( 'all chunk numbers are -999', advance='yes' )
      call OutputNCL2GP_writeGeo_MF (l2gp, l2GPFile, &
        & myswathName, offset)
      if ( timing ) then
        call sayTime( 'Writing geolocations', tFile, t2 )
        tFile = t2
      endif
      call OutputNCL2GP_writeData_MF (l2gp, l2GPFile, &
        & myswathName, offset)
      if ( timing ) then
        call sayTime( 'Writing data', tFile, t2 )
        tFile = t2
      endif
      if ( .not. swath_exists ) then
        call OutputNCL2GP_attributes_MF ( l2gp, l2GPFile, swathName )
        ! NetCDF does not support soft links or aliases
        ! call SetL2GP_aliases_MF (l2gp, l2GPFile, swathName)
        if ( timing ) then
          call sayTime( 'Writing attributes', tFile, t2 )
          tFile = t2
        endif
      end if
      ! l2gp%nTimes = actual_ntimes
    end if

    if ( .not. alreadyOpen )  call mls_closeFile(L2GPFile, Status)
    L2GPFile%errorCode = status
    L2GPFile%lastOperation = 'append'
    call trace_end ( 'AppendNCL2GPData', cond=.false. )
  end subroutine AppendNCL2GPData

  ! ---------------------- CpNCL2GPData_fileID  ---------------------------
  subroutine CpNCL2GPData_fileID( l2metaData, file1, file2, swathList, &
    & notUnlimited, rename, ReadData, &
    & HGrid, rFreqs, rLevels, rTimes, options )
    !------------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine copies swathList from 1 to 2
    ! and, depending on options, makes some repairs or applies a sanity filter
    ! which guarantees that any profiles with Fill values among the
    ! geolocations are marked with DANGERWILLROBINSON status

    use HGridsDatabase, only: HGrid_T
    use L2GPData, only: DestroyL2GPContents, ExtractL2GPRecord, &
      & FilterL2GP, RepairL2GP

    ! Arguments

    type (L2Metadata_T) :: l2metaData
    integer, intent(in)           :: file1 ! handle of file 1
    integer, intent(in)           :: file2 ! handle of file 2
    character (len=*), intent(in) :: swathList ! copy only these; no wildcard
    logical, optional, intent(in) :: notUnlimited
    logical, optional, intent(in) :: ReadData
    character (len=*), optional, intent(in) :: rename
    type (HGrid_T), optional, intent(in)    :: HGrid
    integer, dimension(2), intent(in), optional :: rFreqs  ! subscript range
    integer, dimension(2), intent(in), optional :: rLevels ! subscript range
    integer, dimension(2), intent(in), optional :: rTimes  ! subscript range
    character (len=*), optional, intent(in) :: options ! E.g., '-v'

    ! Local variables
!     integer :: chunk
    logical, parameter            :: countEmpty = .true.
    logical :: filter
    integer :: i
    type (L2GPData_T) :: l2gp
    type (L2GPData_T) :: reducedl2gp
    integer, dimension(2) :: myFreqs  ! subscript range
    integer, dimension(2) :: myLevels ! subscript range
    integer, dimension(2) :: myTimes  ! subscript range
    character (len=8) :: myOptions
    integer :: noSwaths
    logical :: renameSwaths
    logical :: repair
    character (len=L2GPNameLen) :: swath
    character (len=L2GPNameLen) :: swath2
    logical :: verbose
    
    ! Executable code
    myOptions = ' '
    if ( present(options) ) myOptions = options
    verbose = ( index(myOptions, 'v') > 0 )
    repair  = ( index(myOptions, 'r') > 0 )
    filter  = ( index(myOptions, 'f') > 0 )
    renameSwaths  = present(rename)
    if ( renameSwaths ) renameSwaths = ( rename /= ' ' )
    if ( repair .and. .not. present(HGrid) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'cpL2GPFile must be given HGrid to repair L2GPData' )
    if ( .not. repair .and. present(HGrid) ) &
      & call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'cpL2GPFile was given HGrid but no order to repair L2GPData' )
    noSwaths = NumStringElements( trim(swathList), countEmpty )
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths to cp to file--unable to count swaths in ' // trim(swathList) )
    endif
    if ( verbose ) call dump(swathlist, 'swath names')
    if ( renameSwaths ) then
      noSwaths = min ( noSwaths, NumStringElements(trim(rename), countEmpty) )
      if ( verbose ) call dump(rename, 'swath names (copied)')
    endif
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths cp to file--unable to count swaths in ' // trim(rename) )
    endif
    ! Loop over swaths in file 1
    do i = 1, noSwaths
      call GetStringElement (trim(swathList), swath, i, countEmpty )
      if ( renameSwaths ) then
        call GetStringElement (trim(rename), swath2, i, countEmpty )
      else
        swath2 = swath
      endif
      ! Allocate and fill l2gp
      if ( DEEBUG ) print *, 'Reading swath from file: ', trim(swath)
      call ReadNCL2GPData ( file1, trim(swath), l2gp, &
           & ReadData=ReadData )
      if ( all(l2gp%chunkNumber == -999) ) &
        & call output( 'all chunk numbers are -999', advance='yes' )
      if ( DEEBUG ) then
        print *, 'Writing swath to file: ', trim(swath)
        print *, 'l2gp%nFreqs:  ', l2gp%nFreqs
        print *, 'l2gp%nLevels: ', l2gp%nLevels
        print *, 'l2gp%nTimes:  ', l2gp%nTimes
        print *, 'shape(l2gp%l2gpvalue):  ', shape(l2gp%l2gpvalue)
        print *, 'l2gp%MissingL2GP:  ', l2gp%MissingL2GP
        print *, 'l2gp%Missingvalue:  ', l2gp%Missingvalue
      endif
      ! call dump ( l2gp%chunkNumber, 'chunkNumber as Read' )
      ! chunk = FindFirst( l2gp%chunkNumber /= -999 )
      ! call outputNamedValue ( 'index of non-Fill chunkNumber', chunk )
      ! Possibly repair or filter l2gp
      if ( repair ) call RepairL2GP( l2gp, HGrid, options=options )
      if ( filter ) call FilterL2GP( l2gp )
      if ( all(l2gp%chunkNumber == -999) ) &
        & call output( 'After repair and filtering, all chunk numbers are -999', advance='yes' )
      ! Possibly extract a reduced l2pg
      if ( present(rFreqs) .or. present(rLevels) .or. present(rTimes) ) then
        if ( present(rFreqs) ) then
          myFreqs = rFreqs
        else
          myFreqs(1) = 1
          myFreqs(2) = l2gp%nFreqs
        endif
        if ( present(rLevels) ) then
          myLevels = rLevels
        else
          myLevels(1) = 1
          myLevels(2) = l2gp%nLevels
        endif
        if ( present(rTimes) ) then
          myTimes = rTimes
        else
          myTimes(1) = 1
          myTimes(2) = l2gp%nTimes
        endif
        call ExtractL2GPRecord ( l2gp, reducedl2gp, myFreqs, myLevels, myTimes )
        call WriteNCL2GPData(reducedl2gp, file2, trim(swath2), &
          & notUnlimited=notUnlimited )
        call DestroyL2GPContents ( reducedl2gp )
      else
      ! Write the filled l2gp to file2
        ! call dump ( l2gp%chunkNumber, 'chunkNumber as Written' )
        ! chunk = FindFirst( l2gp%chunkNumber /= -999 )
        ! call outputNamedValue ( 'index of non-Fill chunkNumber', chunk )
        call WriteNCL2GPData ( l2gp, file2, trim(swath2), &
          & notUnlimited=notUnlimited )
      endif

      l2metaData%minLat = MLSMin( l2gp%latitude, l2gp%MissingValue )
      l2metaData%maxLat = maxval(l2gp%latitude)
      l2metaData%minLon = MLSMin( l2gp%longitude, l2gp%MissingValue )
      l2metaData%maxLon = maxval(l2gp%longitude)

      call DestroyL2GPContents ( l2gp )
    enddo
       
  end subroutine CpNCL2GPData_fileID
  
  ! ---------------------- CpNCGlobalAttr  ---------------------------
  subroutine CpNCGlobalAttr( File1Handle, File2Handle, status )
  ! Copy global attributes from file 1 to file 2
    ! Args
    integer, intent(in)            :: File1Handle                      
    integer, intent(in)            :: File2Handle                      
    integer, intent(out)           :: status                           
    ! Internal variables
    type(GlobalAttributes_T)       :: gAttributes                      
    ! Executable
    if ( DEEBUG ) then
      call output( 'Before reading global attributes', advance='yes' )
      call outputNamedValue( 'StartUTC', trim(GlobalAttributes%StartUTC) )
    endif
    call ReadNCGlobalAttr ( File1Handle, gAttributes, status )
    if ( status == 0 ) then
      GlobalAttributes = gAttributes
      call WriteNCGlobalAttr( File2Handle )
    endif
  end subroutine CpNCGlobalAttr

  ! ---------------------- DumpNCL2GPData_MLSFile  -----------------------------

  subroutine DumpNCL2GPData_MLSFile( L2GPFile, swathname, &
       firstProf, lastProf, DumpData, DumpAttributes )
    !------------------------------------------------------------------------

    ! Given a NetCDF formatted file,
    ! This routine Reads and then Dumps an L2GP structure
    ! returning a filled data structure and the !
    ! number of profiles Dump.

    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    type(MLSFile_T)               :: L2GPFile
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    logical, optional, intent(in) :: DumpAttributes
    logical, optional, intent(in) :: DumpData

    ! Local
    logical                       :: myDumpData
    logical                       :: myDumpAttributes
    integer                       :: numProfs ! Number actually Dump
    integer                       :: status
    type( l2GPData_T )            :: l2gp ! Result
    
    ! Executable code
    status = 0
    myDumpData = .true.
    if ( present(DumpData) ) myDumpData = DumpData
    if ( present(DumpAttributes) ) myDumpAttributes = DumpAttributes
    call ReadNCL2GPData_MLSFile( L2GPFile, swathname, l2gp, numProfs=numProfs, &
       & firstProf=firstProf, lastProf=lastProf, &
       & ReadData=.true. )
    if ( myDumpData ) &
      & call Dump ( l2gp )
    if ( myDumpAttributes ) &
      & call DumpL2GP_attributes_NC_MF( L2GPFile, l2gp, swathName )
    L2GPFile%errorCode = status
    L2GPFile%lastOperation = 'Dump'
    call DestroyL2GPContents ( l2gp )
  end subroutine DumpNCL2GPData_MLSFile

  !------------------------------------------  IsNCL2GPInFile  -----
  function IsNCL2GPInFile ( File, swath ) result( ItIs )
    use MLSNetCDF4, only: MLS_Swath_In_File
    ! return TRUE only if the swath is already in the file

    ! Dummy arguments
    type (MLSFile_T), intent(in)   :: File
    character(len=*), intent(in)   :: swath
    logical                        :: ItIs
    ! Executable
    call MLS_OpenFile( File )
    Itis = MLS_Swath_In_File ( File%FileID%f_id, swath )
  end function IsNCL2GPInFile

  ! --------------------------------------  OutputNCL2GP_createFile_MF  -----
  ! Create the swath group
  ! Populate the new groug with its
  ! (a) Dimensions
  ! (2) Geolocation fields
  ! (3) Data fields
  subroutine OutputNCL2GP_createFile_MF (l2gp, L2GPFile, &
    & swathName, nLevels, notUnlimited, compressTimes)

    use MLSNetCDF4, only: MLS_DFldsetup, MLS_GFldsetup, MLS_Swcreate, &
      & MLS_SWDetach, MLS_SWDefdim
    ! Brief description of subroutine
    ! This subroutine sets up the structural definitions in an empty L2GP file.

    ! Arguments
    type(MLSFile_T)                :: L2GPFile
    type( L2GPData_T ), intent(inout) :: l2gp
    character (len=*), optional, intent(in) :: swathName ! Defaults to l2gp%swathName
    integer, optional, intent(in) :: nLevels
    logical, optional, intent(in) :: notUnlimited   !               as nTimes
    logical, optional, intent(in) :: compressTimes  ! don't store nTimesTotal

    ! Variables

    character (len=132) :: NAME   ! From l2gp%name
    character (len=32) :: MYDIM1, MYDIM12, MYDIM123

    ! These are hdf5 chunks, _not_ mls along-track processing chunks 
    integer, dimension(7) :: Chunk_Dims
    integer, dimension(7) :: DimIDs
    integer :: CHUNK_RANK
    integer :: CHUNKTIMES, CHUNKFREQS, CHUNKLEVELS

    integer :: SWID, STATUS
    logical :: myNotUnlimited
    logical :: mycompressTimes
    logical :: deebughere
    integer :: f_dimId       ! NetCDF dimension identifier for nFreqs
    integer :: x_dimId       ! NetCDF dimension identifier for nTimes
    integer :: y_dimId       ! NetCDF dimension identifier for nLevels
    ! Executable
    deebughere = DEEBUG !   .or. .TRUE.
    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    myNotUnlimited = .false.
    if ( present ( notUnlimited ) ) myNotUnlimited = notUnlimited
    mycompressTimes = .false.
    if ( present (compressTimes ) ) mycompressTimes = compressTimes
    ! swid = mls_SWcreate(L2GPFile%FileID%f_id, trim(name), &
    swid = mls_SWcreate( L2GPFile, trim(name) )
    L2GPFile%fileID%sd_id = swid
    
    ! Work out the chunking
    if ( myNotUnlimited ) then
      chunkTimes = max ( min ( MAXCHUNKTIMES, l2gp%nTimes ), 1 )
    else
      chunktimes = MAXCHUNKTIMES      ! was 1
    endif
    chunkfreqs = max ( l2gp%nFreqs, 1)
    if(present(nLevels))then
       chunklevels = nLevels
    else
      chunklevels = min(l2gp%nLevels, 500)     ! was .., 5)
      chunklevels = max(chunklevels, 1)
    endif
    
    ! The Dimensions
    ! Define the dimensions. NetCDF will hand back an ID for each. 
    if ( mycompressTimes ) then
      if ( deebughere ) call outputNamedValue ('nTimes',  max(l2gp%nTimes,1) )
      x_dimid = mls_swdefdim( swid, 'nTimes', max(l2gp%nTimes,1) )
    else
      if ( deebughere ) call outputNamedValue ('nTimes',  max(l2gp%nTimesTotal,1) )
      x_dimid = mls_swdefdim( swid, 'nTimes', max(l2gp%nTimesTotal,1) )
    endif
    y_dimid = mls_swdefdim( swid, 'nLevels', l2gp%nLevels )
    ! dimids =  (/ y_dimid, x_dimid, z_dimid /)
    if ( l2gp%nFreqs > 0 ) then
      f_dimid = mls_swdefdim( swid, 'nFreqs', l2gp%nFreqs )
    endif
    
    ! Create the swath within the file

    if ( deebughere )  print *, 'About to sw_create ', TRIM(name)
    if ( deebughere )  print *, 'myNotUnlimited ', myNotUnlimited
    if ( deebughere )  print *, 'l2gp%MissingL2GP ', l2gp%MissingL2GP
    if ( deebughere )  print *, 'l2gp%MissingValue ', l2gp%MissingValue
    if ( swid == -1 ) then
       call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Failed to create swath ' // TRIM(name) &
        & // ' (maybe has the same name as another swath in this file?)', &
            & MLSFile=L2GPFile )
    end if

    if ( myNotUnlimited ) then
      myDim1 = DIM_NAME1
      myDim12 = DIM_NAME12
      myDim123 = DIM_NAME123
    else
      myDim1 = MAX_DIML1
      myDim12 = MAX_DIML12
      myDim123 = MAX_DIML123
    endif
    if ( deebughere ) then
      print *, 'myDim1 ', myDim1
      print *, 'myDim12 ', myDim12
      print *, 'myDim123 ', myDim123
      print *, 'nTimes ', l2gp%nTimes
      print *, 'nTimesTotal ', l2gp%nTimesTotal
      print *, 'nLevels ', l2gp%nLevels
      print *, 'nFreqs ', l2gp%nFreqs
    endif
    ! Define horizontal geolocation fields using x_dimid
    chunk_rank=1
    chunk_dims=1
    chunk_dims(1)=CHUNKTIMES
    if ( deebughere )  print *, 'chunk_dims ', chunk_dims(1:chunk_rank)
    status = mls_gfldsetup(swid, 'Latitude', 'nTimes', MYDIM1, &
      & nf90_real, chunk_rank, chunk_dims, (/x_DimID /), &
      & rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'Longitude', 'nTimes', MYDIM1, &
      & nf90_real, chunk_rank, chunk_dims, (/x_DimID /), &
      & rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'Time', 'nTimes', MYDIM1, &
      & nf90_double, chunk_rank, chunk_dims, (/x_DimID /), &
      & dFill=real(l2gp%MissingValue, r8))

    status = mls_gfldsetup(swid, 'LocalSolarTime', 'nTimes', MYDIM1, &
      & nf90_real, chunk_rank, chunk_dims, (/x_DimID /), &
      & rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'SolarZenithAngle', 'nTimes', MYDIM1, &
      & nf90_real, chunk_rank, chunk_dims, (/x_DimID /), &
      & rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'LineOfSightAngle', 'nTimes', MYDIM1, &
      & nf90_real, chunk_rank, chunk_dims, (/x_DimID /), &
      & rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'OrbitGeodeticAngle', 'nTimes', MYDIM1, &
      & nf90_real, chunk_rank, chunk_dims, (/x_DimID /), &
      & rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'ChunkNumber', 'nTimes', MYDIM1, &
      & nf90_int, chunk_rank, chunk_dims, (/x_DimID /), &
      & iFill=UndefinedIntegerValue)

    if ( l2gp%nLevels > 0 ) then

      chunk_rank=1
      chunk_dims(1)=CHUNKLEVELS
      status = mls_gfldsetup(swid, 'Pressure', 'nLevels', MAX_DIML, &
      & nf90_real, chunk_rank, chunk_dims, (/y_DimID /), &
      & rFill=l2gp%MissingValue)
    end if

    if ( l2gp%nFreqs > 0 ) then

      chunk_rank=1
      chunk_dims(1)=CHUNKFREQS
      status = mls_gfldsetup(swid, 'Frequency', 'nFreqs', MAX_DIML, &
      & nf90_real, chunk_rank, chunk_dims, (/f_DimID /), &
      & rFill=l2gp%MissingValue)
    end if

    ! Define data fields using above dimensions

    if ( (l2gp%nFreqs > 0) .and. (l2gp%nLevels > 0) ) then
       chunk_rank=3
       chunk_dims(1:3)=(/ CHUNKFREQS,CHUNKLEVELS,CHUNKTIMES /)
       dimids(1:chunk_rank) =  (/ f_dimid, y_dimid, x_dimid /)

      status = mls_dfldsetup(swid, 'L2gpValue', 'nFreqs,nLevels,nTimes', &
      & MYDIM123, &
      & nf90_real, chunk_rank, chunk_dims, dimids(1:chunk_rank), &
      & rFill=l2gp%MissingL2GP)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nFreqs,nLevels,nTimes', &
      & MYDIM123, &
      & nf90_real, chunk_rank, chunk_dims, dimids(1:chunk_rank), &
      & rFill=l2gp%MissingValue)

    else if ( l2gp%nLevels > 0 ) then
       chunk_rank=2
       chunk_dims(1:7)=(/ CHUNKLEVELS,CHUNKTIMES,37,38,39,47,49/)
       dimids(1:chunk_rank) =  (/ y_dimid, x_dimid /)

      status = mls_dfldsetup(swid, 'L2gpValue', 'nLevels,nTimes', &
      & MYDIM12, &
      & nf90_real, chunk_rank, chunk_dims, dimids(1:chunk_rank), &
      & rFill=l2gp%MissingL2GP)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nLevels,nTimes', &
      & MYDIM12, &
      & nf90_real, chunk_rank, chunk_dims, dimids(1:chunk_rank), &
      & rFill=l2gp%MissingValue)

    else
       chunk_rank=1
       chunk_dims(1)=CHUNKTIMES
       dimids(1:chunk_rank) =  (/ x_dimid /)

      status = mls_dfldsetup(swid, 'L2gpValue', 'nTimes', &
      & MYDIM1, &
      & nf90_real, chunk_rank, chunk_dims, dimids(1:chunk_rank), &
      & rFill=l2gp%MissingL2GP)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nTimes', &
      & MYDIM1, &
      & nf90_real, chunk_rank, chunk_dims, dimids(1:chunk_rank), &
      & rFill=l2gp%MissingValue)

    end if

    if ( deebughere )  print *, 'chunk_dims ', chunk_dims(1:chunk_rank)
    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    dimids(1:chunk_rank) =  (/ x_dimid /)
    if ( deebughere )  print *, 'chunk_dims ', chunk_dims(1:chunk_rank)
    status = mls_dfldsetup(swid, 'Status', 'nTimes', &
    & MYDIM1, &
    & nf90_int, chunk_rank, chunk_dims, dimids(1:chunk_rank), &
    & iFill=l2gp%MissingStatus)

    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    if ( deebughere )  print *, 'chunk_dims ', chunk_dims(1:chunk_rank)
    status = mls_dfldsetup(swid, 'Quality', 'nTimes', &
    & MYDIM1, &
    & nf90_real, chunk_rank, chunk_dims, dimids(1:chunk_rank), &
    & rFill=l2gp%MissingValue)

    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    status = mls_dfldsetup(swid, 'Convergence', 'nTimes', &
    & MYDIM1, &
    & nf90_real, chunk_rank, chunk_dims, dimids(1:chunk_rank), &
    & rFill=l2gp%MissingValue)

    ! Detach from the HE5_SWath interface.This stores the swath info within the
    ! file and must be done before writing or reading data to or from the
    ! swath. (May be un-necessary for HDF5 -- test program works OK without.)
    ! 
    status = mls_SWdetach(swid)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Failed to detach from swath interface after definition.', &
            & MLSFile=L2GPFile )
    end if

  !--------------------------------------
  end subroutine OutputNCL2GP_createFile_MF
  !--------------------------------------

  ! ---------------------- ReadNCL2GPData_fileID  -----------------------------

  subroutine ReadNCL2GPData_fileID ( L2FileHandle, swathname, l2gp, numProfs, &
       firstProf, lastProf, ReadData )
    !------------------------------------------------------------------------

    ! Given a file handle,
    ! This routine reads an L2GP file, in either hdfVersion,
    ! returning a filled data structure and the !
    ! number of profiles read.
    ! if present, hdfVersion must be one of HDFVERSION_4, HDFVERSION_5

    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    integer, intent(in) :: L2FileHandle ! Returned by swopen
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( l2GPData_T ), intent(out) :: l2gp ! Result
    integer, intent(out), optional :: numProfs ! Number actually read
    logical, optional, intent(in) :: ReadData

    ! Local
    integer :: status
    type( MLSFile_T ) :: l2gpFile
    
    ! Executable code
    status = InitializeMLSFile ( l2gpFile, type=l_NetCDF4, access=DFACC_RDONLY, &
      & content='l2gp', name='unknown' )
    l2gpFile%FileID%f_id = l2FileHandle
    l2gpFile%stillOpen = .true.
    call ReadNCL2GPData ( l2gpFile, swathname, l2gp, numProfs, &
       firstProf, lastProf, ReadData )
  end subroutine ReadNCL2GPData_fileID

  ! ---------------------- ReadNCL2GPData_fileName  -----------------------------

  subroutine ReadNCL2GPData_fileName ( fileName, swathname, l2gp, numProfs, &
       firstProf, lastProf, ReadData )
    !------------------------------------------------------------------------

    ! Given a file name,
    ! This routine reads an L2GP file, in either hdfVersion,
    ! returning a filled data structure and the !
    ! number of profiles read.
    ! hdfVersion may be WILDCARDHDFVERSION
    use MLSFiles, only: HDFVersion_5
    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    character (len=*), intent(in) :: fileName ! Name of swath
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( l2GPData_T ), intent(out) :: l2gp ! Result
    integer, intent(out), optional :: numProfs ! Number actually read
    logical, optional, intent(in) :: ReadData

    ! Local
    type(MLSFile_T)                :: MLSFile
    integer :: status
    
    ! Executable code
!     status = InitializeMLSFile ( MLSFile, type=l_hdf, access=DFACC_READ, &
    status = InitializeMLSFile ( MLSFile, type=l_NetCDF4, access=DFACC_READ, &
     & hdfVersion=HDFVERSION_5, name=trim(fileName) )
    call MLS_OpenFile( MLSFile )
!     L2FileHandle = MLSFile%FileID%f_id
!     call ReadNCL2GPData_fileID ( L2FileHandle, swathname, l2gp, numProfs=numProfs, &
!        & firstProf=firstProf, lastProf=lastProf, &
!        & ReadData=ReadData )
    call ReadNCL2GPData_MLSFile ( MLSFile, swathname, l2gp, numProfs=numProfs, &
       & firstProf=firstProf, lastProf=lastProf, &
       & ReadData=ReadData )
    call MLS_CloseFile( MLSFile )
  end subroutine ReadNCL2GPData_fileName

  ! ---------------------- ReadNCL2GPData_MLSFile  -----------------------------

  subroutine ReadNCL2GPData_MLSFile(L2GPFile, swathname, l2gp, numProfs, &
       firstProf, lastProf, ReadData)
    !------------------------------------------------------------------------

    ! Given a file,
    ! This routine reads an L2GP structure, in either hdfVersion,
    ! returning a filled data structure and the !
    ! number of profiles read.
    ! if present, hdfVersion must be one of HDFVERSION_4, HDFVERSION_5

    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    type(MLSFile_T)                :: L2GPFile
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( l2GPData_T ), intent(out) :: l2gp ! Result
    integer, intent(out), optional :: numProfs ! Number actually read
    logical, optional, intent(in) :: ReadData

    ! Local
    logical :: alreadyOpen
    integer :: status
    
    ! Executable code
    status = 0
    alreadyOpen = L2GPFile%stillOpen
    if ( .not. alreadyOpen ) then
      call mls_openFile(L2GPFile, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l2gp file', MLSFile=L2GPFile)
    endif
    call ReadNCL2GPData_MF_NC4 ( L2GPFile, swathname, l2gp,&
      & numProfs, firstProf, lastProf, ReadData )
    ! if ( .not. alreadyOpen )  call mls_closeFile ( L2GPFile, Status )
    L2GPFile%errorCode = status
    L2GPFile%lastOperation = 'read'
  end subroutine ReadNCL2GPData_MLSFile

  ! ------------------- ReadNCL2GPData_MF_NC4 ----------------

  subroutine ReadNCL2GPData_MF_NC4( L2GPFile, swathname, l2gp, &
    & numProfs, firstProf, lastProf, ReadData )
  use MLSNetCDF4, only: MLS_SWAttach, MLS_SWDetach, &
    & MLS_SWRdLAttr, MLS_SWRdfld
    !------------------------------------------------------------------------

    ! This routine reads an L2GP file that has been written in NetCDF format.
    ! Returns a filled data structure and the !
    ! number of profiles read.

    ! All the ReadConvergence harrumphing is because convergence is newly added
    ! Arguments

    character (len=*), intent(in)   :: swathname ! Name of swath
    type(MLSFile_T)                 :: L2GPFile
    integer, intent(in), optional   :: firstProf, lastProf ! Defaults to first and last
    type( L2GPData_T ), intent(out) :: l2gp ! Result
    integer, intent(out),optional   :: numProfs ! Number actually read
    logical, optional, intent(in)   :: ReadData

    ! Local Parameters
    character (len=*), parameter :: MLSMSG_INPUT = 'Error in input argument '
    integer, parameter           :: MAXDIMSIZE = 4
    ! Local Variables
    character (len=80) :: DF_Name
    character (len=80) :: DF_Precision
    character (len=80) :: list
    character (len=8)  :: dimName
    character (len=480) :: msr

    integer, dimension(MAXDIMSIZE) :: dims
    integer :: first
    integer :: freq
    integer :: grpId
    integer :: i
    integer :: lev
    integer :: nDims
    integer :: nFlds
    integer :: swid
    integer :: status
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: nFreqs
    integer :: nFreqsOr1
    integer :: nLevelsOr1
    integer :: nLevels
    integer :: nTimes
    integer :: nTimesTotal
    integer :: myNumProfs
    logical :: firstCheck, lastCheck
    integer, dimension(MaxRank) :: dimIds
    integer :: varId

    real(rgp), dimension(1)         :: MissingValue
    real(r4), pointer, dimension(:) :: REALFREQ
    real(r4), pointer, dimension(:) :: REALSURF
    real(r4), pointer, dimension(:) :: REALPROF
    real(r4), pointer, dimension(:,:,:) :: REAL3
    logical :: ReadingData
    logical :: deeBugHere
    integer, dimension(MaxRank) :: start
    integer, dimension(MaxRank) :: edge
    integer, dimension(MaxRank) :: stride
    integer, dimension(MaxFlds) :: varIDs
    ! Executable code
    call trace_begin ( me,  'ReadNCL2GPData_MF_NC4', cond=.false. )
    deeBugHere = DEEBUG ! .or. .true.
    nullify ( realFreq, realSurf, realProf, real3 )
    ! Don't fail when trying to read an mls-specific field 
    ! if the file is from another Aura instrument
    ReadingData = .true.
    if ( present(ReadData) ) ReadingData = ReadData
    if ( deeBugHere ) call Dump ( L2GPFile, Details=2 )
    ! Attach to the swath for reading
    l2gp%Name = swathname
    ! We have suffered surprises when hefeos character fields 
    ! were not initialized
    list = ' '
    
    if ( deeBugHere ) print *, 'Trying to attach ', trim(l2gp%Name)
    swid = mls_SWattach( L2GPFile, l2gp%Name )
    if ( deeBugHere ) print *, 'swid: ',swid
    DF_Name = DATA_FIELD1
    DF_Precision = DATA_FIELD2
    if ( deeBugHere ) print *, 'DF_NAME: ',DF_NAME
    if ( deeBugHere ) print *, 'DF_Precision: ',DF_Precision
    ! Here we read the Missing Value attribute
    ! print *, 'Trying to read missing data'
    status = MLS_SWRdLAttr( swid, 'L2gpValue', 'MissingValue', &
      & nf90_real, RealBuffer=MissingValue )
    if ( status /= 0 ) then
      call MLSMessage( MLSMSG_Warning, ModuleName, &
       &'Failed to read Missing Value attribute for ' &
       & // trim(swathname), MLSFile=L2GPFile )
    else
      l2gp%MissingL2GP = MissingValue(1)
    endif
    if (swid == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
         &'Failed to attach to NetCDF4 swath interface for reading' &
         & // trim(swathname), MLSFile=L2GPFile)

    ! Get dimension information
    ! Look inside Data Fields
!     print *, 'seek group name and id'
    call check( &
      & nf90_inq_ncid( swId, 'Data Fields', grpId ), &
      & 'ReadNCL2GPData_MF_NC4 (varId):' // trim(l2gp%name) )
    lev = 0
    freq = 0
    nTimesTotal = 0
    nFreqs = 0
!     print *, 'Trying to read L2gpValue'
    call check( nf90_inq_varid(grpId, 'L2gpValue', varid), &
      & 'ReadNCL2GPData_MF_NC4 (varId):' // trim(l2gp%name) )
    call check( nf90_inquire_variable(grpId, varid, dimids=dimids, nDims=nDims), &
      & 'ReadNCL2GPData_MF_NC4 (dimIds):' // trim(l2gp%name) )
    do i=1, NDims
      call check( nf90_inquire_dimension ( &
        & grpId, dimids(i), name=dimname, len=edge(i) &
        & ), &
        & 'ReadNCL2GPData_MF_NC4 (dim(i)):' // trim(l2gp%name) )
      select case (trim(dimname))
      case ('nLevels')
        nlevels = edge(i)
      case ('nFreqs')
        nFreqs = edge(i)
      case ('nTimes')
        nTimes = edge(i)
      case ('nTimesTotal')
        nTimesTotal = edge(i)
      end select
    enddo
    ! print *, 'nDims as read from NetCDF ', nDims
    ! print *, 'edge as read from NetCDF ', edge
    nDims = min( nDims, MAXDIMSIZE )

    dims(1:nDims) = edge(1:nDims)                                  
    if ( nDims < MAXDIMSIZE ) dims(nDims+1:) = 0 ! Just to make sure they're defined, not junk

    ! Get data and geolocation field information
    call check( nf90_inq_varIds(swid, nFlds, varIds), &
      & 'ReadNCL2GPData_MF_NC4 (varIds):' // trim(l2gp%name) )

    if ( deeBugHere ) print *, 'swathName: ', l2gp%name
    if ( deeBugHere ) print *, 'ndims: ', ndims
    if ( deeBugHere ) print *, 'dims: ', dims
    if (nDims == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'Failed to get dimension information on NetCDF4 swath ' // &
      & trim(swathname), MLSFile=L2GPFile)
    if ( index(list,'nLevels') /= 0 ) lev = 1
    if ( index(list,'Freq') /= 0 ) freq = 1
    if ( nTimesTotal == 0 ) nTimesTotal = nTimes

    l2gp%nTimes      = nTimes     
    l2gp%nTimesTotal = nTimesTotal
    l2gp%nLevels     = nLevels    
    l2gp%nFreqs      = nFreqs     


    ! Check optional input arguments

    firstCheck = present(firstProf)
    lastCheck = present(lastProf)

    if (firstCheck) then
      ! Note that if time is an umlimited dimension, HDF-EOS won't 
      ! nTimes is wrong.
       if ( (firstProf >= l2gp%nTimes) &
         .or. (firstProf < 0) ) then
          msr = MLSMSG_INPUT // 'firstProf'
          call MLSMessage(MLSMSG_Error, ModuleName, msr, MLSFile=L2GPFile)
       else
          first = firstProf
       endif

    else

       first = 0

    endif

    if (lastCheck) then

       if (lastProf < first) then
          msr = MLSMSG_INPUT // 'lastProf'
          call MLSMessage(MLSMSG_Error, ModuleName, msr, MLSFile=L2GPFile)
       endif

       ! If user has supplied "last" _and_ time is unlimited, we have
       ! to believe the user about how many profiles there are.
       ! This is _crap_ and is a temporary workaround. 
!       if(timeIsUnlim) then
!         myNumProfs = lastProf - first + 1
!         nTimes=lastprof-first+1
!         l2gp%nTimes=nTimes
!       endif
       if (lastProf >= nTimes) then
          myNumProfs = nTimes - first
       else
          myNumProfs = lastProf - first + 1
       endif

    else

       myNumProfs = l2gp%nTimes - first

    endif

    ! Allocate result
    if ( deeBugHere ) then
      print *, 'nFreqs: ', nFreqs
      print *, 'nLevels: ', nLevels
      print *, 'mynumProfs: ', mynumProfs
    endif
    call SetupNewL2GPRecord (l2gp, nFreqs=nFreqs, nLevels=nLevels, &
      &  nTimes=mynumProfs, FillIn = .true.)

    ! Skip reading any of the data?
    if ( .not. readingData ) then
      status = mls_SWdetach(swid)
      if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
           &detach from swath interface after reading.', MLSFile=L2GPFile)
      if (present(numProfs)) numProfs=myNumProfs
      call trace_end (  'ReadNCL2GPData_MF_NC4', cond=.false. )
      if(DEEBUG) call outputNamedValue( 'no data read; nTimes', myNumProfs )
      return
    endif
    ! Allocate temporary arrays

    nFreqsOr1=max(nFreqs,1)
    nLevelsOr1=max(nLevels, 1)

    call Allocate_test ( realProf, myNumProfs, 'realProf', ModuleName )
    call Allocate_test ( realSurf, l2gp%nLevels, 'realSurf', ModuleName )
    call Allocate_test ( realFreq, l2gp%nFreqs, 'realFreq', ModuleName )
    call Allocate_test ( real3, nFreqsOr1, nLevelsOr1, myNumProfs, 'real3', ModuleName )

    ! Read the horizontal geolocation fields

    start(1) = 1 ! 0
    start(2) = 1 ! 0
    start(3) = 1 ! first
    stride = 1
    edge(1) = nFreqsOr1
    edge(2) = nLevelsOr1
    edge(3) = myNumProfs
    status=0

    if ( DeeBugHere ) then
    print *, 'start ', start
    print *, 'stride ', stride
    print *, 'edge ', edge
    print *, 'shape(realProf) ', shape(realProf)
    endif
    status = mls_SWrdfld(swid, 'Latitude', start(3:3), stride(3:3), &
      edge(3:3), realProf)
    l2gp%latitude = realProf
    print *, 'min  max (latitude) ', minval(realProf), maxval(realProf)

    status = mls_SWrdfld(swid, 'Longitude', start(3:3), stride(3:3), edge(3:3),&
      &    realProf)
    l2gp%longitude = realProf

    status = mls_SWrdfld(swid, 'Time', start(3:3), stride(3:3), edge(3:3),&
      &    l2gp%time)

    status = mls_SWrdfld(swid, 'LocalSolarTime', start(3:3), stride(3:3), edge(3:3),&
      &    realProf)
    l2gp%solarTime = realProf

    status = mls_SWrdfld(swid, 'SolarZenithAngle', start(3:3), stride(3:3), edge(3:3),&
      &    realProf)
    l2gp%solarZenith = realProf

    ! These next 3 are MLS-specific
    status = mls_SWrdfld(swid, 'LineOfSightAngle', start(3:3), stride(3:3), edge(3:3),&
      &    realProf)
    l2gp%losAngle = realProf

    status = mls_SWrdfld(swid, 'OrbitGeodeticAngle', start(3:3), stride(3:3), edge(3:3),&
      &   realProf)
    l2gp%geodAngle = realProf

    status = mls_SWrdfld(swid, 'ChunkNumber', start(3:3), stride(3:3), edge(3:3),&
      &    l2gp%chunkNumber)

    ! Read the pressures vertical geolocation field, if it exists

    if (NLevels /= 0) then
       status = mls_SWrdfld(swid,'Pressure',start(2:2),stride(2:2), edge(2:2),&
         & realSurf)
       l2gp%pressures = realSurf
    endif

    ! Read the frequency geolocation field, if it exists

    if (freq == 1) then

       edge(1) = l2gp%nFreqs

       status = mls_SWrdfld(swid,'Frequency',start(1:1),stride(1:1),edge(1:1),&
         & realFreq)
       l2gp%frequency = realFreq

    endif

    ! Read the data fields that may have 1-3 dimensions
    ! Why would start now be set back to 0 instead of 1?
    ! start = 0

    if ( any(freq == (/0,1/)) ) then

!       status = mls_SWrdfld(swid, trim(DF_Name), start, stride, edge, real3)
!       status = mls_SWrdfld(swid, 'L2gpValue', start, stride, edge, real3)
       status = mls_SWrdfld(swid, 'L2gpValue', start(2:3), stride(2:3), edge(2:3), real3(1,:,:))
       l2gp%l2gpValue = real3
       if ( DeeBugHere ) then
       print *, 'Min, max (real3) ', minval(real3), maxval(real3)
       print *, 'Min, max (l2gpval) ', minval(l2gp%l2gpValue), maxval(l2gp%l2gpValue)
       print *, 'shape (l2gpval) ', shape(l2gp%l2gpValue)
       print *, 'status ', status       
!        print *, '3rd index'    
!        call Dump( l2gp%l2gpValue(1,1,:), trim(swathname) )
       print *, '2nd index'    
       call Dump( l2gp%l2gpValue(1,:,1), trim(swathname) )
       endif

       status = mls_SWrdfld(swid, trim(DF_Precision), start(2:3), stride(2:3), edge(2:3), real3(1,:,:))
       l2gp%l2gpPrecision = real3
       if ( DeeBugHere ) then
       print *, 'Min, max (precision) ', minval(l2gp%l2gpPrecision), maxval(l2gp%l2gpPrecision)
       print *, 'shape (precision) ', shape(l2gp%l2gpPrecision)
       endif

    else if ( lev == 1) then

       status = mls_SWrdfld( swid, trim(DF_Name), start(2:3), stride(2:3), &
            edge(2:3), real3(1,:,:) )
       l2gp%l2gpValue = real3

       status = mls_SWrdfld( swid, trim(DF_Precision), start(2:3), stride(2:3), &
            edge(2:3), real3(1,:,:) )
       l2gp%l2gpPrecision = real3

    else

       print *, start(3:3),stride(3:3),edge(3:3)
       print *, shape(real3(1,1,:))
       status = mls_SWrdfld(swid,trim(DF_Name),start(3:3),stride(3:3),edge(3:3),&
         &   real3(1,1,:) )
       l2gp%l2gpValue = real3

       status = mls_SWrdfld( swid, trim(DF_Precision), start(3:3), stride(3:3), edge(3:3), &
            real3(1,1,:) )
       l2gp%l2gpPrecision = real3

    endif
    
    ! Read the data fields that are 1-dimensional

    l2gp%status = l2gp%MissingStatus ! l2gp%MissingValue ! So it has a value.
    status = mls_swrdfld( swid, 'Status',start(3:3),stride(3:3),edge(3:3),&
      & l2gp%status )
     print *, 'Min, max (status) ', minval(l2gp%status), maxval(l2gp%status)

    status = mls_SWrdfld(swid, 'Quality', start(3:3), stride(3:3),&
      edge(3:3),realProf)
    l2gp%quality = realProf

    l2gp%Convergence = l2gp%MissingValue ! l2gp%MissingValue ! So it has a value.
    status = mls_swrdfld( swid, 'Convergence',start(3:3),stride(3:3),edge(3:3),&
      & l2gp%convergence )


    ! Deallocate local variables
    call Deallocate_test ( realProf, 'realProf', ModuleName )
    call Deallocate_test ( realSurf, 'realSurf', ModuleName )
    call Deallocate_test ( realFreq, 'realFreq', ModuleName )
    call Deallocate_test ( real3, 'real3', ModuleName )

    !  After reading, detach from swath
    status = mls_SWdetach ( swid )
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &detach from NC4 swath interface after reading.', MLSFile=L2GPFile)
    !print*," leaving ReadNCL2GPData_hdf: first/last/read=",&
    !  firstprof,lastprof,myNumProfs
    ! Set numProfs if wanted
    if (present(numProfs)) numProfs=myNumProfs
    ! Must we close the MLS file somehow?
    if ( deeBugHere ) call Dump ( L2GPFile, Details=2 )
!    call h5fclose_f ( L2GPFile%fileID%f_id, status )
!     if (status == -1) call MLSMessage( MLSMSG_Error, ModuleName, 'Failed to &
!          & close NC4 file after reading.', MLSFile=L2GPFile )
    call trace_end (  'ReadNCL2GPData_MF_NC4', cond=.false. )

  end subroutine ReadNCL2GPData_MF_NC4

!------------------------------------------------------------
   subroutine ReadNCGlobalAttr_FileID ( fileID, gAttributes, status )
!------------------------------------------------------------
! Should eventually be moved to PCFHdr module
    use MLSNetCDF4, only: MLS_SwRdattr
! Brief description of subroutine
! This subroutine reads the global attributes for a NetCDF4 file

! Arguments

      integer, intent(in)            :: fileID
      type(GlobalAttributes_T)       :: gAttributes        
      integer                        :: status
! Internal variables
      logical, parameter :: DeeBug = .false.
      integer, dimension(1) :: ibuf
! Executable
      if ( DeeBug ) print *, 'Reading NetCDF4 Global attributes'
      status = MLS_SwRdattr(fileID, &
       & 'identifier_product_doi', nf90_char, 1, &
       &  gAttributes%DOI)
      status = MLS_SwRdattr(fileID, &
       & 'ProductionLocation', nf90_char, 1, &
       &  gAttributes%productionLoc)
      status = MLS_SwRdattr(fileID, &
            & 'OrbitNumber', nf90_int, max_orbits, &
            &  intbuffer=gAttributes%OrbNum)
      status = MLS_SwRdattr(fileID, &
            & 'OrbitPeriod', nf90_double, max_orbits, &
            &  Doublebuffer=gAttributes%OrbPeriod)
      status = MLS_SwRdattr(fileID, &
       & 'InstrumentName', nf90_char, 1, &
       &  gAttributes%InstrumentName)
      status = MLS_SwRdattr(fileID, &
         & 'HostName', nf90_char, 1, &
         &  gAttributes%HostName)
      status = MLS_SwRdattr(fileID, &
       & 'ProcessLevel', nf90_char, 1, &
       &  gAttributes%ProcessLevel)
      status = MLS_SwRdattr(fileID, &
       & 'PGEVersion', nf90_char, 1, &
       &  gAttributes%PGEVersion)
      status = MLS_SwRdattr(fileID, &
       & 'StartUTC', nf90_char, 1, &
       &  gAttributes%StartUTC)
      status = MLS_SwRdattr(fileID, &
       & 'EndUTC', nf90_char, 1, &
       &  gAttributes%EndUTC)

      status = MLS_SwRdattr(fileID, &
       & 'GranuleDayOfYear', nf90_int, 1, &
       &  intBuffer=ibuf )
      gAttributes%DayOfYear=ibuf(1)
      status = MLS_SwRdattr(fileID, &
       & 'GranuleDay', nf90_int, 1, &
       &  int=gAttributes%GranuleDay )
      status = MLS_SwRdattr(fileID, &
       & 'GranuleMonth', nf90_int, 1, &
       &  int=gAttributes%GranuleMonth )
      status = MLS_SwRdattr(fileID, &
       & 'GranuleYear', nf90_int, 1, &
       &  int=gAttributes%GranuleYear )

      status = MLS_SwRdattr(fileID, &
       & 'TAI93At0zOfGranule', nf90_double, 1, &
       &  Doublesca= gAttributes%TAI93At0zOfGranule )
      if ( lowercase(gAttributes%ProcessLevel(1:2)) == 'l2' ) then
        status = MLS_SwRdattr(fileID, &
         & 'FirstMAF', nf90_int, 1, &
         &  int=gAttributes%FirstMAFCtr )
        status = MLS_SwRdattr(fileID, &
         & 'LastMAF', nf90_int, 1, &
         &  int=gAttributes%LastMAFCtr )
      endif

      status = MLS_SwRdattr(fileID, &
       & 'geos5 type', nf90_char, 1, &
       &  gAttributes%geos5_type)
      status = MLS_SwRdattr(fileID, &
       & 'identifier_product_doi', nf90_char, 1, &
       &  gAttributes%doi)
      status = MLS_SwRdattr(fileID, &
        & 'MiscNotes', nf90_char, 1, &
        &  gAttributes%MiscNotes)
      status = MLS_SwRdattr(fileID, &
        & 'HostName', nf90_char, 1, &
        &  gAttributes%HostName)
      status = MLS_SwRdattr(fileID, &
        & 'ProductionLocation', nf90_char, 1, &
        &  gAttributes%ProductionLoc)
      if ( deebug ) then
        call output( 'Reading global attributes', advance='yes' )
        call dumpGlobalAttributes ( gAttributes )
      endif
!------------------------------------------------------------
   end subroutine ReadNCGlobalAttr_FileID
!------------------------------------------------------------

!------------------------------------------------------------
   subroutine WriteNCFileAttr ( MLSFile )
!------------------------------------------------------------

   use MLSNetCDF4, only: MLS_SwWrattr
! Brief description of subroutine
! This subroutine writes the components of an MLSFile_t 
! as attributes for an NetCDF4-formatted file

! Arguments

      type(MLSFile_T)       :: MLSFile
! Local variables
      integer :: fileID
      integer :: status

      ! Executable code
      if ( .not. MLSFile%stillOpen ) then
        call MLS_OpenFile( MLSFile )
      endif
      fileID = MLSFile%FileID%f_id
      status = MLS_SwWrattr ( fileID, &
       & 'content', nf90_char, 1, &
       &  charBuffer=MLSFile%content )
      status = MLS_SwWrattr ( fileID, &
       & 'lastOperation', nf90_char, 1, &
       &  charBuffer=MLSFile%lastOperation )
      status = MLS_SwWrattr ( fileID, &
       & 'name', nf90_char, 1, &
       &  charBuffer=MLSFile%name )
      status = MLS_SwWrattr ( fileID, &
       & 'ShortName', nf90_char, 1, &
       &  charBuffer=MLSFile%ShortName )
      status = MLS_SwWrattr ( fileID, &
       & 'typeStr', nf90_char, 1, &
       &  charBuffer=MLSFile%typeStr )
      status = MLS_SwWrattr (fileID, &
       & 'type', nf90_int, 1, &
       &  intbuffer=(/ MLSFile%type /) )
      status = MLS_SwWrattr (fileID, &
       & 'access', nf90_int, 1, &
       &  int= MLSFile%access  )
      status = MLS_SwWrattr ( fileID, &
       & 'HDFVersion', nf90_int, 1, &
       &  intbuffer=(/ MLSFile%HDFVersion /) )
      status = MLS_SwWrattr ( fileID, &
       & 'PCFID', nf90_int, 1, &
       &  int= MLSFile%PCFId  )
      status = MLS_SwWrattr ( fileID, &
       & 'recordlength', nf90_int, 1, &
       &  int= MLSFile%recordlength  )
      status = MLS_SwWrattr ( fileID, &
       & 'errorCode', nf90_int, 1, &
       &  int=MLSFile%errorCode )
      status = MLS_SwWrattr ( fileID, &
       & 'PCFIDRange', nf90_int, 2, &
       &  intBuffer=(/ MLSFile%PCFIDRange%Bottom, MLSFile%PCFIDRange%Top /) )
      status = MLS_SwWrattr ( fileID, &
       & 'FileID', nf90_int, 3, &
       &  intBuffer=(/MLSFile%FileID%f_id, MLSFile%FileID%grp_id, MLSFile%FileID%sd_id/) )

!------------------------------------------------------------
   end subroutine WriteNCFileAttr
!------------------------------------------------------------

  !----------------------------------------  WriteNCL2GPData_fileID  -----
  ! This subroutine is an amalgamation of the last three
  ! Should be renamed CreateAndWriteNCL2GPData
  subroutine WriteNCL2GPData_fileID ( l2gp, l2FileHandle, swathName, &
    & notUnlimited )

    ! Arguments

    integer, intent(in) :: l2FileHandle ! From swopen
    type (L2GPData_T), intent(INOUT) :: l2gp
    character (len=*), optional, intent(in) ::swathName!default->l2gp%swathName
    logical, optional, intent(in) :: notUnlimited
    ! Exectuable code

    ! Local
    integer :: status
    type( MLSFile_T ) :: l2gpFile

    ! Executable code
    status = InitializeMLSFile(l2gpFile, type=l_NetCDF4, access=DFACC_RDWR, &
      & content='l2gp', name='unknown')
    l2gpFile%FileID%f_id = l2FileHandle
    l2gpFile%stillOpen = .true.
    call WriteNCL2GPData( l2gp, l2gpFile, swathName, notUnlimited )

  end subroutine WriteNCL2GPData_fileID

  !----------------------------------------  WriteNCL2GPData_MLSFile  -----
  ! This subroutine is an amalgamation of the last three
  ! Should be renamed CreateAndWriteNCL2GPData
  subroutine WriteNCL2GPData_MLSFile(l2gp, L2GPFile, swathName, &
    & notUnlimited)

    ! Arguments

    type(MLSFile_T)                :: L2GPFile
    type (L2GPData_T), intent(INOUT) :: l2gp
    character (len=*), optional, intent(in) ::swathName!default->l2gp%swathName
    logical, optional, intent(in) :: notUnlimited
    ! Local
    logical :: alreadyOpen
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: status
    logical, parameter :: DEEBUG = .false.
    ! Executable code

    call trace_begin ( me, 'WriteNCL2GPData_MLSFile', cond=.false. )
    status = 0
    alreadyOpen = L2GPFile%stillOpen
    if ( .not. alreadyOpen ) then
      print *, 'Opening ' // trim(L2GPFile%Name)
      call mls_openFile(L2GPFile, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l2gp file', MLSFile=L2GPFile)
    endif
!     call Dump ( L2GPFile )
    if ( L2GPFile%access == DFACC_RDONLY )  &
      & call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'l2gp file is rdonly', MLSFile=L2GPFile)
    if ( DEEBUG ) print *, 'About to OutputNCL2GP_createFile_MF'
    call OutputNCL2GP_createFile_MF (l2gp, L2GPFile, &
      & swathName, notUnlimited=notUnlimited)
    call OutputNCL2GP_writeGeo_MF (l2gp, L2GPFile, &
      & swathName)
    call OutputNCL2GP_writeData_MF (l2gp, L2GPFile, &
      & swathName)
    if ( DEEBUG ) print *, 'Outputting attributes'
    call OutputNCL2GP_attributes_MF (l2gp, L2GPFile, swathName)
    if ( DEEBUG ) print *, 'NetCDF unable to create aliases'
    ! call SetL2GP_aliases_MF (l2gp, L2GPFile, swathName)

    if ( .not. alreadyOpen )  call mls_closeFile(L2GPFile, Status)
    L2GPFile%errorCode = status
    L2GPFile%lastOperation = 'write'
    call trace_end ( 'WriteNCL2GPData_MLSFile', cond=.false. )
  end subroutine WriteNCL2GPData_MLSFile
  
!------------------------------------------------------------
   subroutine WriteNCGlobalAttr_FileID ( fileID, dayNum, DOI, MiscNotes, &
     & skip_if_already_there )
!------------------------------------------------------------
! Should eventually be moved to PCFHdr module
    use MLSNetCDF4, only: MLS_SwwrAttr, MLS_IsGlAtt
    use PCFHdr, only: GranuleDayOfYear, GranuleDay, GranuleMonth, &
      & ProcessLevelFun
! Brief description of subroutine
! This subroutine writes the global attributes for a NetCDF4 file

! Arguments

      integer, intent(in) :: fileID
      integer, intent(in), optional :: DayNum
      logical, intent(in), optional :: Doi          ! Write DOI
      logical, intent(in), optional :: MiscNotes    ! Write MiscNotes
      logical, intent(in), optional :: Skip_if_already_there
! Internal variables
      integer :: status
      logical :: myDOI
      logical, parameter :: DeeBug = .false.
      logical :: my_skip
! Executable
      if ( deebug ) print *, 'Writing NetCDF4 Global attributes'
      myDOI = .false.
      if ( present(DOI) ) myDOI=DOI
      my_skip = .false.
      if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
      if ( deebug ) then
        call output( 'Writing global attributes', advance='yes' )
        call OutputNamedValue ( 'myDOI', myDOI )
        call OutputNamedValue ( 'DOI', GlobalAttributes%DOI )
        call OutputNamedValue ( 'ProductionLocation', &
          & GlobalAttributes%productionLoc )
        call dumpGlobalAttributes
      endif
      if ( len_trim(GlobalAttributes%DOI) > 0 ) &
       & status = MLS_Swwrattr(fileID, &
       & 'identifier_product_doi', nf90_char, 1, &
       &  charBuffer=GlobalAttributes%DOI)
      if ( len_trim(GlobalAttributes%productionLoc) > 0 ) &
       & status = MLS_Swwrattr(fileID, &
       & 'ProductionLocation', nf90_char, 1, &
       &  charBuffer=GlobalAttributes%productionLoc)
      if ( my_skip ) then
        if ( MLS_IsGlAtt ( fileID, 'OrbitNumber' ) ) &
          & return
      endif
      if (present(dayNum)) then
         status = MLS_Swwrattr(fileID, &
            & 'OrbitNumber', nf90_int, max_orbits, &
            &  intBuffer=GlobalAttributes%OrbNumDays(:,dayNum))
         status = MLS_Swwrattr(fileID, &
            & 'OrbitPeriod', nf90_double, max_orbits, &
            &  DoubleBuffer=GlobalAttributes%OrbPeriodDays(:,dayNum))
      else
         status = MLS_Swwrattr(fileID, &
            & 'OrbitNumber', nf90_int, max_orbits, &
            &  intBuffer=GlobalAttributes%OrbNum)
         status = MLS_Swwrattr(fileID, &
            & 'OrbitPeriod', nf90_double, max_orbits, &
            &  DoubleBuffer=GlobalAttributes%OrbPeriod)
      end if
      status = MLS_Swwrattr(fileID, &
       & 'InstrumentName', nf90_char, 1, &
       &  charBuffer=GlobalAttributes%InstrumentName)
      if ( len_trim(GlobalAttributes%HostName) > 0 ) then
        status = MLS_Swwrattr(fileID, &
         & 'HostName', nf90_char, 1, &
         &  charBuffer=GlobalAttributes%HostName)
      else
        status = MLS_Swwrattr(fileID, &
         & 'HostName', nf90_char, 1, &
         &  charBuffer=GlobalAttributes%ProductionLoc)
      endif
      status = MLS_Swwrattr(fileID, &
       & 'ProcessLevel', nf90_char, 1, &
       &  charBuffer=ProcessLevelFun())
      status = MLS_Swwrattr(fileID, &
       & 'PGEVersion', nf90_char, 1, &
       &  charBuffer=GlobalAttributes%PGEVersion)
      status = MLS_Swwrattr(fileID, &
       & 'StartUTC', nf90_char, 1, &
       &  charBuffer=GlobalAttributes%StartUTC)
      status = MLS_Swwrattr(fileID, &
       & 'EndUTC', nf90_char, 1, &
       &  charBuffer=GlobalAttributes%EndUTC)
      if ( GlobalAttributes%GranuleDay < 1 ) return
      status = MLS_Swwrattr(fileID, &
       & 'GranuleMonth', nf90_int, 1, &
       &  intBuffer=(/ GranuleMonth() /) )
      status = MLS_Swwrattr(fileID, &
       & 'GranuleDay', nf90_int, 1, &
       &  intBuffer=(/ GranuleDay() /) )
      status = MLS_Swwrattr(fileID, &
       & 'GranuleDayOfYear', nf90_int, 1, &
       &  intBuffer=(/ GranuleDayOfYear() /) )
      status = MLS_Swwrattr(fileID, &
       & 'GranuleYear', nf90_int, 1, &
       &  int= GlobalAttributes%GranuleYear )
      status = MLS_Swwrattr(fileID, &
       & 'TAI93At0zOfGranule', nf90_double, 1, &
       &  Doublesca= GlobalAttributes%TAI93At0zOfGranule )
      if ( index(lowercase(ProcessLevelFun()), 'l2') > 0 ) then
        status = MLS_Swwrattr(fileID, &
         & 'FirstMAF', nf90_int, 1, &
         &  int= GlobalAttributes%FirstMAFCtr  )
        status = MLS_Swwrattr(fileID, &
         & 'LastMAF', nf90_int, 1, &
         &  int= GlobalAttributes%LastMAFCtr )
      endif
      status = MLS_Swwrattr(fileID, &
      & 'MiscNotes', nf90_char, 1, &
      &  GlobalAttributes%MiscNotes)
      status = MLS_Swwrattr(fileID, &
      & 'eos5 type', nf90_char, 1, &
      &  GlobalAttributes%geos5_type)
!------------------------------------------------------------
   end subroutine WriteNCGlobalAttr_FileID
!------------------------------------------------------------

! ======================= Private Procedures ===========================
  ! --------------------------------------  check  -----
  ! Check that the NetCDF operation was successful; i.e., nf90_noerr
  ! If not, print mesg (unless silent)
  ! and halt with error (unless FailureOK)
  subroutine check( status, mesg, FailureOK, silent )
    use MLSNetCDF4, only: NCError
    integer, intent ( in )                          :: status ! how was it?
    character(len=*), intent ( in )                 :: mesg
    logical, optional, intent ( in )                :: FailureOK ! Don't halt
    logical, optional, intent ( in )                :: silent ! Don't print
    !
    logical                                         :: MyFailOK
    logical                                         :: MySilent
    ! Executable
    NCError = status ! So the caller knows in case we don't exit
    MyFailOK = .false.
    if ( present(FailureOK) ) MyFailOK = FailureOK
    MySilent = .false.
    if ( present(Silent) ) MySilent = Silent
!     print *, 'status   ', status
!     print *, 'MySilent ', MySilent
    if( status /= nf90_noerr ) then 
      if ( .not. MySilent ) &
        & call output( trim(nf90_strerror(status)), advance='yes' )
      if ( MyFailOK ) return
      print *, trim(nf90_strerror(status))
      call MLSMessage ( MLSMSG_Error, moduleName,  &
          & trim(mesg) )
    end if
  end subroutine check  
  
  !----------------------------------------  OutputNCL2GP_attributes_MF  -----
  subroutine OutputNCL2GP_attributes_MF( l2gp, L2GPFile, swathName )

  use MLSNetCDF4, only: MLS_ISGlatt, &
    & MLS_SWAttach, MLS_SWDetach, MLS_SWWrattr, MLS_SWWrlattr
    ! Brief description of subroutine
    ! This subroutine writes the attributes for an l2gp
    ! These include
    ! Global attributes which are file level
    ! Swath level
    ! Geolocation field attributes
    ! Data field attributes
    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    type(MLSFile_T)                :: L2GPFile
    character (len=*), intent(IN), optional :: swathName ! Defaults->l2gp%name

    ! Parameters
    character (len=*), parameter :: NOUNITS = 'NoUnits'

    ! The following pair of string list encode the Units attribute
    ! corresponding to each Title attribute; e.g., the Units for Latitude is deg
    character (len=*), parameter :: GeolocationTitles = &
      & 'Latitude,Longitude,Time,LocalSolarTime,SolarZenithAngle,' // &
      & 'LineOfSightAngle,OrbitGeodeticAngle,ChunkNumber,Pressure,Frequency'
    character (len=*), parameter :: GeolocationUnits = &
      & 'deg,deg,s,h,deg,' // &
      & 'deg,deg,NoUnits,hPa,GHz'
    character (len=*), parameter :: GeoUniqueFieldDefinition = &
      & 'HMT,HMT,AS,HMT,HMT,' // &
      & 'M,M,M,AS,M'   ! These are abbreviated values
    character (len=*), parameter :: UniqueFieldDefKeys = &
      & 'HM,HMT,MT,AS,M'
    character (len=*), parameter :: UniqueFieldDefValues = &
      & 'HIRDLS-MLS-Shared,HIRDLS-MLS-TES-Shared,MLS-TES-Shared,' // &
      & 'Aura-Shared,MLS-Specific'  ! Expanded values
    ! The following associate UniqueFieldDefs with species names
    character (len=*), parameter :: Species = &
      & 'Temperature,BrO,CH3CN,CO,ClO,GPH,HCl,HCN,H2O,H2O2,' // &
      & 'HNO3,HOCl,HO2,N2,N2O,OH,O2,O3,RHI,SO2'
    character (len=*), parameter :: SpUniqueFieldDefinition = &
      & 'HMT,M,M,MT,M,M,M,M,HMT,M,' // &
      & 'HMT,M,M,M,HM,M,M,HMT,M,M'   ! These are abbreviated values
    ! logical, parameter :: DEEBUG = .true.

    ! Variables
    character (len=132) :: name     ! Either swathName or l2gp%name
    character(len=CHARATTRLEN) :: abbr_uniq_fdef
    character(len=CHARATTRLEN) :: expnd_uniq_fdef
    integer :: field
    character(len=CHARATTRLEN) :: field_name
    logical :: isColumnAmt
    logical :: isTPPressure
    integer :: rgp_type
    character(len=CHARATTRLEN) :: species_name ! Always lower case
    integer :: status
    integer :: swid
    character(len=CHARATTRLEN) :: temp_name
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theTitles
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theUnits
    character(len=CHARATTRLEN) :: units_name
    ! Begin
    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    if ( DEEBUG ) print *, 'About to wr attrs to: ', trim(name)
    if ( rgp == r4 ) then
      rgp_type = nf90_real
    elseif ( rgp == r8 ) then
      rgp_type = nf90_double
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & 'Attributes have unrecognized numeric data type; should be r4 or r8', &
        & MLSFile=L2GPFile)
    endif
    call List2Array(GeolocationTitles, theTitles, .true.)
    call List2Array(GeolocationUnits, theUnits, .true.)

    ! - -   G l o b a l   A t t r i b u t e s   - -
    if ( WRITEMASTERSFILEATTRIBUTES .and. &
      & .not. MLS_ISGLATT(L2GPFile%fileID%f_id, 'StartUTC') ) then
      if ( DEEBUG ) print *, 'Writing Global attributes; e.g., StartUTC'
      call WriteNCGlobalAttr( L2GPFile%fileID%f_id )
    endif

    swid = mls_SWattach (L2GPFile, name)
    
    !   - -   S w a t h   A t t r i b u t e s   - -
    status = mls_swwrattr( swid, trim(l2gp%verticalCoordinate), &
      & rgp_type, size(l2gp%pressures), &
      & realbuffer=l2gp%pressures )
    field_name = l2gp%verticalCoordinate ! 'Pressure'
    status = mls_swwrattr(swid, 'VerticalCoordinate', nf90_char, 1, &
      & charBuffer=field_name)
    if ( DeeBug ) print *, 'Missing L2GP: ', real(l2gp%MissingL2GP)
    if ( DeeBug ) print *, 'Missing value: ', real(l2gp%MissingValue)
    
    !   - -   G e o l o c a t i o n   A t t r i b u t e s   - -
    if ( DEEBUG ) print *, 'About to wr loc attrs to: ', trim(name)
    do field=1, NumGeolocFields
      ! Take care not to write attributes to "missing fields"
      if ( trim(theTitles(field)) == 'Frequency' &
        & .and. l2gp%nFreqs < 1 ) then
        field_name = ''
      elseif ( trim(theTitles(field)) == 'Pressure' &
        & .and. l2gp%nLevels < 1 ) then
        field_name = ''
      else
        field_name = theTitles(field)
        if ( trim(theTitles(field)) == 'Pressure' ) &
          & field_name = l2gp%verticalCoordinate
        call GetHashElement (GeolocationTitles, &
          & GeoUniqueFieldDefinition, trim(theTitles(field)), &
          & abbr_uniq_fdef, .false.)
        call GetHashElement (UniqueFieldDefKeys, &
          & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
          & expnd_uniq_fdef, .false.)
        if ( DEEBUG ) print *, 'Field Title ', trim(theTitles(field))
        status = mls_swwrlattr(swid, trim(theTitles(field)), 'Title', &
          & nf90_char, 1, charBuffer=theTitles(field))
        if ( DEEBUG ) print *, 'Units ', trim(theUnits(field))
        status = mls_swwrlattr(swid, trim(theTitles(field)), 'Units', &
          & nf90_char, 1, charBuffer=theUnits(field))

        if ( trim(theTitles(field)) == 'Time' ) then
          status = mls_swwrlattr(swid, trim(theTitles(field)), &
            & 'MissingValue', &
            & nf90_double, 1, DoubleSca= real(l2gp%MissingValue, r8) )
        elseif ( trim(theTitles(field)) == 'ChunkNumber' ) then
          status = mls_swwrlattr(swid, trim(theTitles(field)), &
            & 'MissingValue', &
            & nf90_int, 1, int= UndefinedIntegerValue )
        else
          status = mls_swwrlattr(swid, trim(theTitles(field)), &
            & 'MissingValue', &
            & rgp_type, 1, realsca= real(l2gp%MissingValue, rgp) )
        endif
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "Unable to write local attribute to " // trim(theTitles(field)) )
        if ( DEEBUG ) print *, 'Uniquefielddef ', trim(expnd_uniq_fdef)
        status = mls_swwrlattr(swid, trim(theTitles(field)), &
          & 'UniqueFieldDefinition', &
          & nf90_char, 1, charBuffer=expnd_uniq_fdef)
        ! print *, 'status : ', status
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "Unable to write local attribute to " // trim(theTitles(field)) )
      endif
    enddo
    if ( DEEBUG ) print *, 'Data'
    !   - -   D a t a   A t t r i b u t e s   - -
    ! call GetQuantityAttributes ( l2gp%quantityType, &
    !  & units_name, expnd_uniq_fdef)
    field_name = Name
    species_name = lowercase(name)
    ! The following special cases are handled by crude, despicable hacks
    ! It would be better to check on the quantity type, but paw hasn't
    ! succeeded in getting that to work properly and reliably
    isColumnAmt = ( index(species_name, 'column') > 0 )
    isTPPressure = ( index(species_name, 'tpp') > 0 )
    if ( isColumnAmt ) then
      ! This next failed sometimes
      temp_name = species_name
      call ExtractSubString(Name, species_name, 'column', 'wmo')
      if ( species_name == ' ' ) species_name=temp_name
      ! So we'll try another tack:
      temp_name = species_name
      call ReplaceSubString ( temp_name, species_name, 'column', '' )
      if ( species_name == ' ' ) species_name=temp_name
      temp_name = species_name
      call GetStringElement( temp_name, species_name, &
        & 1, countEmpty=.true., inseparator='-' )
      if ( species_name == ' ' ) species_name=temp_name
    else
      temp_name = species_name
      call GetStringElement( temp_name, species_name, &
        & 1, countEmpty=.true., inseparator='-' )
      if ( species_name == ' ' ) species_name=temp_name
    endif
    ! Hopefully by now we've turned species_name into one of the species
    if ( DEEBUG ) &
      & print *, 'full name: ', trim(name), ' Species: ', trim(species_name)
    call GetHashElement (lowercase(Species), &
      & SpUniqueFieldDefinition, trim(species_name), &
      & abbr_uniq_fdef, .false.)
    call GetHashElement (UniqueFieldDefKeys, &
      & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
      & expnd_uniq_fdef, .false.)
    if ( expnd_uniq_fdef == '' .or. expnd_uniq_fdef == ',' ) &
      & expnd_uniq_fdef = 'MLS-Specific'
    select case (trim(lowercase(species_name)))
    case ('temperature')
      units_name = 'K'
    case ('gph')
      units_name = 'm'
    case ('geoheight')
      units_name = 'km'
    case ('rhi')
      units_name = '%rhi'
    case ('iwc', 'iwp')
      units_name = 'g/m^3'
    case default
      units_name = 'vmr'
    end select
    if ( isColumnAmt ) then
      ! units_name = 'DU'
      if ( DEEBUG ) &
        & call dump( .true., col_species_keys, col_species_hash, &
        & 'column species units' )
      call GetHashElement (col_species_keys, &
      & col_species_hash, trim(lowercase(species_name)), &
      & units_name, .true.)
      if ( DEEBUG ) &
      & print *, 'units_name: ', trim(units_name)
      if ( units_name == ',' ) units_name = 'molcm2'
    endif
    if ( isTPPressure ) units_name = 'hPa'
    if ( DEEBUG ) print *, 'Title ', trim(field_name)
    status = mls_swwrlattr(swid, 'L2gpValue', 'Title', &
      & nf90_char, 1, charBuffer=field_name)
    if ( DEEBUG ) print *, 'Units ', trim(units_name)
    status = mls_swwrlattr(swid, 'L2gpValue', 'Units', &
      & nf90_char, 1, charBuffer=units_name)
    status = mls_swwrlattr(swid, 'L2gpValue', 'MissingValue', &
      & rgp_type, 1, realsca= real(l2gp%MissingL2GP, rgp) )
    if ( DEEBUG ) print *, 'Title ', trim(expnd_uniq_fdef)
    status = mls_swwrlattr(swid, 'L2gpValue', &
      & 'UniqueFieldDefinition', &
      & nf90_char, 1, charBuffer=expnd_uniq_fdef)
    status = mls_swwrlattr(swid, 'L2gpPrecision', 'Title', &
      & nf90_char, 1, charBuffer=trim(field_name)//'Precision')
    status = mls_swwrlattr(swid, 'L2gpPrecision', 'Units', &
      & nf90_char, 1, charBuffer=units_name)
    status = mls_swwrlattr(swid, 'L2gpPrecision', 'MissingValue', &
      & rgp_type, 1, realsca= real(l2gp%MissingValue, rgp) )
    status = mls_swwrlattr(swid, 'L2gpPrecision', &
      & 'UniqueFieldDefinition', &
      & nf90_char, 1, charBuffer=expnd_uniq_fdef)

    ! ('Status' data field newly written)
    status = mls_swwrlattr(swid, 'Status', 'Title', &
      & nf90_char, 1, charBuffer=trim(field_name)//'Status')
    status = mls_swwrlattr(swid, 'Status', 'Units', &
      & nf90_char, 1, NOUNITS)
    status = mls_swwrlattr(swid, 'Status', 'MissingValue', &
      & nf90_int, 1, int=l2gp%MissingStatus )
    status = mls_swwrlattr(swid, 'Status', &
      & 'UniqueFieldDefinition', &
      & nf90_char, 1, charBuffer='MLS-Specific')
    
    status = mls_swwrlattr(swid, 'Quality', 'Title', &
      & nf90_char, 1, charBuffer=trim(field_name)//'Quality')
    status = mls_swwrlattr(swid, 'Quality', 'Units', &
      & nf90_char, 1, charBuffer=NOUNITS)
    status = mls_swwrlattr(swid, 'Quality', 'MissingValue', &
      & rgp_type, 1, realsca=real(l2gp%MissingValue, rgp) )
    status = mls_swwrlattr(swid, 'Quality', &
      & 'UniqueFieldDefinition', &
      & nf90_char, 1, charBuffer='MLS-Specific')
    
    status = mls_swwrlattr(swid, 'Convergence', 'Title', &
      & nf90_char, 1, charBuffer=trim(field_name)//'Convergence')
    status = mls_swwrlattr(swid, 'Convergence', 'Units', &
      & nf90_char, 1, charBuffer=NOUNITS)
    status = mls_swwrlattr(swid, 'Convergence', 'MissingValue', &
      & rgp_type, 1, realsca=real(l2gp%MissingValue, rgp) )
    status = mls_swwrlattr(swid, 'Convergence', &
      & 'UniqueFieldDefinition', &
      & nf90_char, 1, charBuffer='MLS-Specific')
    
    status = mls_SWdetach( swid )
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from NetCDF4 interface', MLSFile=L2GPFile )
    end if


  !-------------------------------------
  end subroutine OutputNCL2GP_attributes_MF
  !-------------------------------------

  !-----------------------------------------  OutputNCL2GP_writeGeo_MF  -----
  subroutine OutputNCL2GP_writeGeo_MF (l2gp, L2GPFile, &
    & swathName,offset)

  use MLSNetCDF4, only: MLS_SWAttach, MLS_SWDetach, MLS_SWWrfld
    ! Brief description of subroutine
    ! This subroutine writes the geolocation fields to an L2GP output file.

    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    type(MLSFile_T)                :: L2GPFile
    character (len=*), intent(in), optional :: swathName ! Defaults->l2gp%name
    integer,intent(in),optional::offset

    ! Variables

    character (len=132) :: name ! Either swathName or l2gp%name
!     integer :: chunk
    integer :: status, swid,myOffset
    integer :: start(2), stride(2), edge(2)
    logical, parameter:: DEEBUG = .false.

    ! Begin
    if (present(offset)) then
       myOffset=offset
    else
       myOffset=0
    endif

    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif

    swid = mls_SWattach (L2GPFile, name)

    ! Write data to the fields
    edge = 0
    stride = 1
    start = 1 + myOffset ! Please do not set to zero
    edge(1) = l2gp%nTimes

    if (DEEBUG) then
      call output( 'Writing geolocations', advance='yes' )
      call outputNamedValue( 'stride, start, edge', &
        & (/ stride(1), start(1), edge(1) /) )
      call outputNamedValue( 'shape(Geod. Angle', shape(l2gp%geodAngle) )
    endif

    status = mls_SWwrfld(swid, 'Latitude', start, stride, edge, &
         real(l2gp%latitude))

    status = mls_SWwrfld(swid, 'Longitude', start, stride, edge, &
         real(l2gp%longitude))

    status = mls_SWwrfld(swid, 'Time', start, stride, edge, &
         l2gp%time)

    status = mls_SWwrfld(swid, 'LocalSolarTime', start, stride, edge, &
        real(l2gp%solarTime))

    status = mls_SWwrfld(swid, 'SolarZenithAngle', start, stride, edge, &
         real(l2gp%solarZenith))

    status = mls_SWwrfld(swid, 'LineOfSightAngle', start, stride, edge, &
         real(l2gp%losAngle))

    status = mls_SWwrfld(swid, 'OrbitGeodeticAngle', start, stride, edge, &
         real(l2gp%geodAngle))

    status = mls_SWwrfld(swid, 'ChunkNumber', start, stride, edge, &
         l2gp%chunkNumber)

    if ( l2gp%nLevels > 0 ) then
       edge(1) = l2gp%nLevels
       start(1) = 1 ! needed because offset may have made this /=0
       status = mls_SWwrfld(swid, 'Pressure', start(1:1), stride(1:1), edge(1:1), &
            real(l2gp%pressures))
    end if

    if ( l2gp%nFreqs > 0 ) then
       edge(1) = l2gp%nFreqs
       start(1) = 1 ! needed because offset may have made this /=0
       status = mls_SWwrfld(swid, 'Frequency', start, stride, edge, &
            real(l2gp%frequency))
    end if

    ! Detach from the swath interface.  

    status = mls_SWdetach( swid )
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach from NC4 swath interface', MLSFile=L2GPFile )
    end if

  !------------------------------------
  end subroutine OutputNCL2GP_writeGeo_MF
  !------------------------------------

  !----------------------------------------  OutputNCL2GP_writeData_MF  -----
  subroutine OutputNCL2GP_writeData_MF(l2gp, L2GPFile, &
    & swathName,offset)

  use MLSNetCDF4, only: MLS_SWAttach, MLS_SWDetach, MLS_SWWrfld
    ! Brief description of subroutine
    ! This subroutine writes the data fields to an L2GP output file.
    ! For now, you have to write all of l2gp, but you can choose to write
    ! it at some offset into the file
    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    type(MLSFile_T)                :: L2GPFile
    character (len=*), intent(in), optional :: swathName ! Defaults->l2gp%name
    integer,intent(in),optional::offset
    ! Parameters
    ! logical, parameter :: DEEBUG = .false.

    ! Variables

    character (len=132) :: name     ! Either swathName or l2gp%name

    integer :: status,myOffset
    integer :: start(3), stride(3), edge(3)
    integer :: swid
    real :: tFile ! How long have we been fooling with this file?
    

    ! Begin
    if (present(offset)) then
       myOffset=offset
    else
       myOffset=0
    endif

    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif

    start = 1
    stride = 1
    start(3) = 1 + myOffset ! Please do not set to zero
    edge(1) = l2gp%nFreqs
    edge(2) = l2gp%nLevels
    edge(3) = l2gp%nTimes
    call time_now ( tFile )
    if (DEEBUG) then
      call output( 'Writing data now ' // trim(name), advance='yes' )
      call outputNamedValue( 'stride, start, edge', (/ stride(3), start(3), edge(3) /) )
      call outputNamedValue( 'shape(Geod. Angle', shape(l2gp%geodAngle) )
    endif
    swid = mls_SWattach (L2GPFile, name)
    if ( DEEBUG ) then
      call sayTime( 'Attaching swath ' // trim(name), tFile, t2 )
      tFile = t2
    endif
    if ( l2gp%nFreqs > 0 ) then
       ! Value and Precision are 3-D fields
       if (DEEBUG) print *, 'start, stride, edge ', start, stride, edge
       ! What the devil???
       ! Why aren't we writing a 3-d array here??
       status = mls_SWwrfld(swid, 'L2gpValue', start, stride, edge, &
            & reshape(l2gp%l2gpValue, edge) )
       status = mls_SWwrfld(swid, 'L2gpPrecision', start, stride, edge, &
            & reshape(l2gp%l2gpPrecision, edge) )

       if ( DEEBUG ) then
         call sayTime( 'Witing 3d values and precision', tFile, t2 )
         tFile = t2
       endif
    else if ( l2gp%nLevels > 0 ) then
       ! Value and Precision are 2-D fields
      
       if (DEEBUG) print *, 'start, stride, edge: ', start(2:3), stride(2:3), edge(2:3)
       status = mls_SWwrfld( swid, 'L2gpValue', start(2:3), stride(2:3), &
            edge(2:3), l2gp%l2gpValue(1,:,:))
       if ( DEEBUG ) then
         call sayTime( 'Witing 2d values', tFile, t2 )
         tFile = t2
       endif
       status = mls_SWwrfld( swid, 'L2gpPrecision', start(2:3), stride(2:3), &
            edge(2:3), l2gp%l2gpPrecision(1,:,:))
       if ( DEEBUG ) then
         call sayTime( 'Witing 2d precision', tFile, t2 )
         tFile = t2
       endif
    else

       ! Value and Precision are 1-D fields
       if (DEEBUG) print *, 'start, stride, edge: ', start(3), stride(3), edge(3)
       status = mls_SWwrfld( swid, 'L2gpValue', start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpValue(1,1,:) ))
       if ( DEEBUG ) then
         call sayTime( 'Witing 1d values', tFile, t2 )
         tFile = t2
       endif
       status = mls_SWwrfld( swid, 'L2gpPrecision', start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpPrecision(1,1,:) ))
       if ( DEEBUG ) then
         call sayTime( 'Witing 1d precision', tFile, t2 )
         tFile = t2
       endif
    end if
    if ( status /= 0 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to write l2gp values/precicions NC4 swath interface', MLSFile=L2GPFile )
    end if

    ! 1-D status & quality fields

    status = mls_swwrfld(swid, 'Status', start(3:3), stride(3:3), edge(3:3), &
       &   l2gp%status)

    status = mls_SWwrfld(swid, 'Quality', start(3:3), stride(3:3), edge(3:3), &
         real(l2gp%quality))

    status = mls_SWwrfld(swid, 'Convergence', start(3:3), stride(3:3), edge(3:3), &
         real(l2gp%convergence))


    if ( DEEBUG ) then
      call sayTime( 'Witing 1d status, Quality, Asc/Desc', tFile, t2 )
      tFile = t2
    endif
    !     Detach from the swath interface.
    status = mls_SWdetach(swid)
    if(DEEBUG) print *, 'Detached from swid ', swid
    if(DEEBUG) print *, 'file handle ', L2GPFile%FileID%f_id
    if(DEEBUG) print *, 'status ', status
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach from NC4 swath interface', MLSFile=L2GPFile )
    end if


  !-------------------------------------
  end subroutine OutputNCL2GP_writeData_MF
  !-------------------------------------

  !----------------------------------------  DumpL2GP_attributes_NC_MF  -----
  subroutine DumpL2GP_attributes_NC_MF( L2GPFile, l2gp, swathName )
  use HDFEOS5, only: MLS_Chartype
  use MLSHDFEOS, only: MLS_SWAttach, MLS_SWDetach
  use MLSNetCDF4, only: MLS_SwRdLattr
  use PCFHdr, only:  GlobalAttributes_T, &
    & DumpGlobalAttributes
    ! Brief description of subroutine
    ! This subroutine reads and then dumps the attributes for an l2gp
    ! that has been written to a NetCDF-formatted file.
    ! These include
    ! Global attributes which are file level
    ! Swath level
    ! Geolocation field attributes
    ! Data field attributes
    
    ! Arguments

    type(MLSFile_T)                         :: L2GPFile
    type( L2GPData_T ), intent(inout)       :: l2gp
    character (len=*), optional, intent(in) :: swathName ! Defaults->l2gp%name

    ! Variables
    logical                           :: alreadyOpen
    type (GlobalAttributes_T)         :: gAttributes
    integer                           :: DayofYear
    double precision                  :: TAI93At0zOfGranule
    real(rgp), dimension(MAXNLEVELS)  :: pressures

    character (len=132) :: name     ! Either swathName or l2gp%name
    ! The following pair of string list encode the Units attribute
    ! corresponding to each Title attribute; e.g., the Units for Latitude is deg
    character (len=*), parameter :: GeolocationTitles = &
      & 'Latitude,Longitude,Time,LocalSolarTime,SolarZenithAngle,' // &
      & 'LineOfSightAngle,OrbitGeodeticAngle,ChunkNumber,Pressure,Frequency'
    character (len=*), parameter :: GeolocationUnits = &
      & 'deg,deg,s,h,deg,' // &
      & 'deg,deg,NoUnits,hPa,GHz'
    character (len=*), parameter :: GeoUniqueFieldDefinition = &
      & 'HMT,HMT,AS,HMT,HMT,' // &
      & 'M,M,M,AS,M'   ! These are abbreviated values
    character (len=*), parameter :: UniqueFieldDefKeys = &
      & 'HM,HMT,MT,AS,M'
    character (len=*), parameter :: UniqueFieldDefValues = &
      & 'HIRDLS-MLS-Shared,HIRDLS-MLS-TES-Shared,MLS-TES-Shared,' // &
      & 'Aura-Shared,MLS-Specific'  ! Expanded values
    ! The following associate UniqueFieldDefs with species names
    character (len=*), parameter :: Species = &
      & 'Temperature,BrO,CH3CN,CO,ClO,GPH,HCl,HCN,H2O,H2O2,' // &
      & 'HNO3,HOCl,HO2,N2,N2O,OH,O2,O3,RHI,SO2'
    character (len=*), parameter :: SpUniqueFieldDefinition = &
      & 'HMT,M,M,MT,M,M,M,M,HMT,M,' // &
      & 'HMT,M,M,M,HM,M,M,HMT,M,M'   ! These are abbreviated values

    integer                    :: field
    logical                    :: isColumnAmt
    integer                    :: status
    integer                    :: swid
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theTitles
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theUnits
    character(len=CHARATTRLEN) :: field_name
    character(len=CHARATTRLEN) :: units_name
    character(len=CHARATTRLEN) :: uniqueness
    character(len=CHARATTRLEN) :: species_name
    character(len=CHARATTRLEN) :: abbr_uniq_fdef
    character(len=CHARATTRLEN) :: expnd_uniq_fdef
    real(rgp), dimension(1)    :: MissingValue
    ! Begin
    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    call outputNamedValue ( 'L2GP Attributes: (swath name) ', trim(name), &
      & options = '-H' )
    if ( .not. any( rgp /= (/ r4, r8 /) ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & 'Attributes have unrecognized numeric data type; should be r4 or r8' )
    endif
    call List2Array(GeolocationTitles, theTitles, .true.)
    call List2Array(GeolocationUnits, theUnits, .true.)

    ! - -   G l o b a l   A t t r i b u t e s   - -
    alreadyOpen = L2GPFile%stillOpen
    if ( .not. alreadyOpen ) then
      if ( DEEBUG ) print *, 'Needed to open file ', trim(L2GPFile%name)
      call mls_openFile(L2GPFile, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l2gp file', MLSFile=L2GPFile)
    endif
    call output ( '(Global Attributes) ', advance='yes' )
    call ReadNCGlobalAttr_FileID ( L2GPFile%fileID%f_id, gAttributes, status )
    if ( status /= 0 ) then
      call output ('No global attributes found in file', advance='yes')
    else
      ! We stopped reading MiscNotes in the PCFHdr module to avoid
      ! letting the master task step on the value each slave had carefully
      ! written there.
      ! Therefore, we must read it here, now, along with other miscellania.
      ! (Isn't Miscellania a state bordering Missouri and Pennsylvania?)
!       status = MLS_SwRdattr( L2GPFile%fileID%f_id, &
!       & 'MiscNotes', nf90_char, 1, &
!       &  CharBuffer=gAttributes%MiscNotes )
!       status = MLS_SwRdattr( L2GPFile%fileID%f_id, &
!       & 'ProductionLocation', nf90_char, 1, &
!       &  CharBuffer=gAttributes%ProductionLoc )
!       status = MLS_SwRdattr( L2GPFile%fileID%f_id, &
!       & 'HostName', nf90_char, 1, &
!       &  CharBuffer=gAttributes%HostName )
!       status = MLS_SwRdattr( L2GPFile%fileID%f_id, &
!       & 'identifier_product_doi', nf90_char, 1, &
!       &  CharBuffer=gAttributes%DOI )
      call dumpGlobalAttributes ( gAttributes ) 
      call dump_int ( DayOfYear, 'Granule day of year:' )
      call dump_r8 ( TAI93At0zOfGranule, 'Equator crossing time (tai93):' )
    endif
    swid = mls_SWattach ( L2GPFile, name )
    if ( swid == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to attach swath ' // trim(name))
    end if
    
    !   - -   S w a t h   A t t r i b u t e s   - -
    call output ( '(Swath Attributes) ', advance='yes')
    status = MLS_SwRdLattr( swid, 'L2gpValue', 'Pressure', nf90_real, &
      & RealBuffer=l2gp%pressures)
    call dump ( pressures, 'Vertical coordinates:' )
    status = MLS_SwRdLattr( swid, 'L2gpValue', 'VerticalCoordinate', MLS_CHARTYPE, &
      & CharBuffer=field_name)
    call dump_chars ( field_name, 'Vertical coordinates type:' )
    status = MLS_SWRdLAttr( swid, 'L2gpValue', 'MissingValue', &
      & nf90_real, RealBuffer=MissingValue )
    call dump ( MissingValue, 'MissingValues (L2GPValues only):' )
    status = MLS_SwRdLattr( swid, 'L2gpPrecision', 'MissingValue', &
      & nf90_real, RealBuffer=MissingValue )
    call dump ( MissingValue, 'MissingValues (other fields):' )
    
    !   - -   G e o l o c a t i o n   A t t r i b u t e s   - -
    call output ( '(Geolocation Attributes) ', advance='yes')
    do field=1, NumGeolocFields
      ! Take care not to write attributes to "missing fields"
      if ( trim(theTitles(field)) == 'Frequency' &
        & .and. l2gp%nFreqs < 1 ) then
        field_name = ''
      elseif ( trim(theTitles(field)) == 'Pressure' &
        & .and. l2gp%nLevels < 1 ) then
        field_name = ''
      else
        call GetHashElement (GeolocationTitles, &
          & GeoUniqueFieldDefinition, trim(theTitles(field)), &
          & abbr_uniq_fdef, .false.)
        call GetHashElement (UniqueFieldDefKeys, &
          & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
          & expnd_uniq_fdef, .false.)
        status = MLS_swrdlattr(swid, trim(theTitles(field)), 'Title', MLS_CHARTYPE, &
      & CharBuffer=field_name)
        status = MLS_swrdlattr(swid, trim(theTitles(field)), 'Units', MLS_CHARTYPE, &
      & CharBuffer=units_name)
        call dump_chars ( field_name, 'Field title:' )
        call dump_chars ( units_name, 'Units:' )

        status = MLS_swrdlattr(swid, trim(theTitles(field)), &
          & 'UniqueFieldDefinition', MLS_CHARTYPE, &
      & CharBuffer=uniqueness)
        call dump_chars ( uniqueness, 'Unique field definition:' )
      endif
    enddo
    !   - -   D a t a   A t t r i b u t e s   - -
    call output ( '(Data Attributes) ', advance='yes')
    field_name = Name
    species_name = name
    isColumnAmt = ( index(species_name, 'Column') > 0 )
    if ( isColumnAmt ) then
      call ExtractSubString(Name, species_name, 'Column', 'wmo')
    endif
    call GetHashElement (lowercase(Species), &
      & SpUniqueFieldDefinition, trim(lowercase(species_name)), &
      & abbr_uniq_fdef, .false.)
    call GetHashElement (UniqueFieldDefKeys, &
      & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
      & expnd_uniq_fdef, .false.)
    select case (trim(lowercase(species_name)))
    case ('temperature')
      units_name = 'K'
    case ('gph')
      units_name = 'm'
    case ('rhi')
      units_name = '%rhi'
    case default
      units_name = 'vmr'
    end select
    if ( isColumnAmt ) units_name = 'DU'
    status = MLS_swrdlattr(swid, 'L2gpValue', 'Title', MLS_CHARTYPE, &
      & CharBuffer=field_name)
    status = MLS_swrdlattr(swid, 'L2gpValue', 'Units', MLS_CHARTYPE, &
      & CharBuffer=units_name)
    status = MLS_swrdlattr(swid, 'L2gpValue', &
      & 'UniqueFieldDefinition', MLS_CHARTYPE, &
      & CharBuffer=uniqueness)
    call dump_chars ( field_name, 'Field title:' )
    call dump_chars ( units_name, 'Units:' )
    call dump_chars ( uniqueness, 'Unique field definition:' )
    status = MLS_swrdlattr(swid, 'L2gpPrecision', 'Title', MLS_CHARTYPE, &
      & CharBuffer=field_name)
    status = MLS_swrdlattr(swid, 'L2gpPrecision', 'Units', MLS_CHARTYPE, &
      & CharBuffer=units_name)
    status = MLS_swrdlattr(swid, 'L2gpPrecision', &
      & 'UniqueFieldDefinition', MLS_CHARTYPE, &
      & CharBuffer=uniqueness)
    call dump_chars ( field_name, 'Field title:' )
    call dump_chars ( units_name, 'Units:' )
    call dump_chars ( uniqueness, 'Unique field definition:' )

    ! ('Status' data field newly written)
    status = MLS_swrdlattr(swid, 'Status', 'Title', MLS_CHARTYPE, &
      & CharBuffer=field_name)
    status = MLS_swrdlattr(swid, 'Status', 'Units', MLS_CHARTYPE, &
      & CharBuffer=units_name)
    status = MLS_swrdlattr(swid, 'Status', &
      & 'UniqueFieldDefinition', MLS_CHARTYPE, &
      & CharBuffer=uniqueness)
    call dump_chars ( field_name, 'Field title:' )
    call dump_chars ( units_name, 'Units:' )
    call dump_chars ( uniqueness, 'Unique field definition:' )
    
    status = MLS_swrdlattr(swid, 'Quality', 'Title', MLS_CHARTYPE, &
      & CharBuffer=field_name)
    status = MLS_swrdlattr(swid, 'Quality', 'Units', MLS_CHARTYPE, &
      & CharBuffer=units_name)
    status = MLS_swrdlattr(swid, 'Quality', &
      & 'UniqueFieldDefinition', MLS_CHARTYPE, &
      & CharBuffer=uniqueness)
    call dump_chars ( field_name, 'Field title:' )
    call dump_chars ( units_name, 'Units:' )
    call dump_chars ( uniqueness, 'Unique field definition:' )
    
    status = mls_SWdetach( swid )
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from swath interface' )
    end if
    contains
    subroutine dump_chars(value, name)
      ! arguments
      character(len=*), intent(in) :: value
      character(len=*), intent(in) :: name
      ! Executable
      call output(trim(name) // ' ', advance='no')
      call output(trim(value), advance='yes')
    end subroutine dump_chars
    subroutine dump_int(value, name)
      ! arguments
      integer, intent(in) :: value
      character(len=*), intent(in) :: name
      ! Executable
      call output(trim(name) // ' ', advance='no')
      call output(value, advance='yes')
    end subroutine dump_int
    subroutine dump_r8(value, name)
      ! arguments
      real(r8), intent(in) :: value
      character(len=*), intent(in) :: name
      ! Executable
      call output(trim(name) // ' ', advance='no')
      call output(value, advance='yes')
    end subroutine dump_r8

  !-------------------------------------
  end subroutine DumpL2GP_attributes_NC_MF
  !-------------------------------------

  !----------------------------------------  OutputL2GP_attributes_MF  -----
  ! Oops, we already coded a subroutine that does this.
  ! What we really need is a subroutine to read l2gp attributes
  subroutine OutputL2GP_attributes_MF(l2gp, L2GPFile, swathName)

  use HDFEOS5, only: HE5T_Native_Real, HE5T_Native_Double
  use MLSNetCDF4, only: MLS_ISGlatt, &
    & MLS_SWAttach, MLS_SWDetach, MLS_SWWrlattr
    ! Brief description of subroutine
    ! This subroutine writes the attributes for an l2gp
    ! These include
    ! Global attributes which are file level
    ! Swath level
    ! Geolocation field attributes
    ! Data field attributes
    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    type(MLSFile_T)                :: L2GPFile
    character (len=*), intent(in), optional :: swathName ! Defaults->l2gp%name

    ! Parameters

    ! The following pair of string list encode the Units attribute
    ! corresponding to each Title attribute; e.g., the Units for Latitude is deg
    character (len=*), parameter :: GeolocationTitles = &
      & 'Latitude,Longitude,Time,LocalSolarTime,SolarZenithAngle,' // &
      & 'LineOfSightAngle,OrbitGeodeticAngle,ChunkNumber,Pressure,Frequency'
    character (len=*), parameter :: GeolocationUnits = &
      & 'deg,deg,s,h,deg,' // &
      & 'deg,deg,NoUnits,hPa,GHz'

    ! Variables
    character (len=132) :: name     ! Either swathName or l2gp%name
    integer :: rgp_type
    ! character(len=CHARATTRLEN) :: species_name ! Always lower case
    integer :: status
    integer :: swid
    ! character(len=CHARATTRLEN) :: temp_name
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theTitles
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theUnits
    ! character(len=CHARATTRLEN) :: units_name
    ! Begin
    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    if ( DEEBUG ) print *, 'About to wr attrs to: ', trim(name)

    if ( rgp == r4 ) then
      rgp_type = HE5T_NATIVE_REAL
    elseif ( rgp == r8 ) then
      rgp_type = HE5T_NATIVE_DOUBLE
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & 'Attributes have unrecognized numeric data type; should be r4 or r8', &
        & MLSFile=L2GPFile)
    endif
    call List2Array(GeolocationTitles, theTitles, .true.)
    call List2Array(GeolocationUnits, theUnits, .true.)

    ! - -   G l o b a l   A t t r i b u t e s   - -
    if ( WRITEMASTERSFILEATTRIBUTES .and. &
      & .not. MLS_ISGLATT(L2GPFile%fileID%f_id, 'StartUTC') ) &
      call WriteNCGlobalAttr_FileID ( L2GPFile%fileID%f_id )
! 
    swid = mls_SWattach (L2GPFile, name)
!     
    !   - -   S w a t h   A t t r i b u t e s   - -
    ! Under NteCDF, every local attribute must be attached to a dataset
    ! or what they call a 'variable'
    ! So we'll attach swath-level attributes to the L2gpValue.
    
!     status = he5_swwrattr(swid, trim(l2gp%verticalCoordinate), &
!       & rgp_type, hsize(size(l2gp%pressures)), &
!       & l2gp%pressures)
    status = MLS_Swwrlattr ( swid, &
    & 'L2gpValue', trim(l2gp%verticalCoordinate), rgp_type, size(l2gp%pressures), &
    & realbuffer=l2gp%pressures )
!     field_name = l2gp%verticalCoordinate ! 'Pressure'
!     status = mls_swwrattr(swid, 'VerticalCoordinate', MLS_CHARTYPE, 1, &
!       & field_name)
!     if ( DeeBug ) print *, 'Missing L2GP: ', real(l2gp%MissingL2GP)
!     if ( DeeBug ) print *, 'Missing value: ', real(l2gp%MissingValue)
!     if ( SWATHLEVELMISSINGVALUE ) &
!       & status = he5_swwrattr(swid, 'MissingValue', rgp_type, hsize(1), &
!       & (/ real(l2gp%MissingL2GP, rgp) /) )
!     
!     !   - -   G e o l o c a t i o n   A t t r i b u t e s   - -
!     if ( DEEBUG ) print *, 'About to wr loc attrs to: ', trim(name)
!     do field=1, NumGeolocFields
!       ! Take care not to write attributes to "missing fields"
!       if ( trim(theTitles(field)) == 'Frequency' &
!         & .and. l2gp%nFreqs < 1 ) then
!         field_name = ''
!       elseif ( trim(theTitles(field)) == 'Pressure' &
!         & .and. l2gp%nLevels < 1 ) then
!         field_name = ''
!       else
!         field_name = theTitles(field)
!         if ( trim(theTitles(field)) == 'Pressure' ) &
!           & field_name = l2gp%verticalCoordinate
!         call GetHashElement (GeolocationTitles, &
!           & GeoUniqueFieldDefinition, trim(theTitles(field)), &
!           & abbr_uniq_fdef, .false.)
!         call GetHashElement (UniqueFieldDefKeys, &
!           & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
!           & expnd_uniq_fdef, .false.)
!         if ( DEEBUG ) print *, 'Field Title ', trim(theTitles(field))
!         status = mls_swwrlattr(swid, trim(theTitles(field)), 'Title', &
!           & MLS_CHARTYPE, 1, theTitles(field))
!         if ( DEEBUG ) print *, 'Units ', trim(theUnits(field))
!         status = mls_swwrlattr(swid, trim(theTitles(field)), 'Units', &
!           & MLS_CHARTYPE, 1, theUnits(field))
! 
!         if ( trim(theTitles(field)) == 'Time' ) then
!           status = he5_swwrlattr(swid, trim(theTitles(field)), &
!             & 'MissingValue', &
!             & HE5T_NATIVE_DOUBLE, hsize(1), (/ real(l2gp%MissingValue, r8) /) )
!         elseif ( trim(theTitles(field)) == 'ChunkNumber' ) then
!           status = he5_swwrlattr(swid, trim(theTitles(field)), &
!             & 'MissingValue', &
!             & HE5T_NATIVE_INT, hsize(1), (/ UndefinedIntegerValue /) )
!         else
!           status = he5_swwrlattr(swid, trim(theTitles(field)), &
!             & 'MissingValue', &
!             & rgp_type, hsize(1), (/ real(l2gp%MissingValue, rgp) /) )
!         endif
!         if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & "Unable to write local attribute to " // trim(theTitles(field)) )
!         if ( DEEBUG ) print *, 'Uniquefielddef ', trim(expnd_uniq_fdef)
!         status = mls_swwrlattr(swid, trim(theTitles(field)), &
!           & 'UniqueFieldDefinition', &
!           & MLS_CHARTYPE, 1, expnd_uniq_fdef)
!         ! print *, 'status : ', status
!         if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & "Unable to write local attribute to " // trim(theTitles(field)) )
!       endif
!     enddo
!     if ( DEEBUG ) print *, 'Data'
!     !   - -   D a t a   A t t r i b u t e s   - -
!     ! call GetQuantityAttributes ( l2gp%quantityType, &
!     !  & units_name, expnd_uniq_fdef)
!     field_name = Name
!     species_name = lowercase(name)
!     ! The following special cases are handled by crude, despicable hacks
!     ! It would be better to check on the quantity type, but paw hasn't
!     ! succeeded in getting that to work properly and reliably
!     isColumnAmt = ( index(species_name, 'column') > 0 )
!     isTPPressure = ( index(species_name, 'tpp') > 0 )
!     if ( isColumnAmt ) then
!       ! This next failed sometimes
!       temp_name = species_name
!       call ExtractSubString(Name, species_name, 'column', 'wmo')
!       if ( species_name == ' ' ) species_name=temp_name
!       ! So we'll try another tack:
!       temp_name = species_name
!       call ReplaceSubString ( temp_name, species_name, 'column', '' )
!       if ( species_name == ' ' ) species_name=temp_name
!       temp_name = species_name
!       call GetStringElement( temp_name, species_name, &
!         & 1, countEmpty=.true., inseparator='-' )
!       if ( species_name == ' ' ) species_name=temp_name
!     else
!       temp_name = species_name
!       call GetStringElement( temp_name, species_name, &
!         & 1, countEmpty=.true., inseparator='-' )
!       if ( species_name == ' ' ) species_name=temp_name
!     endif
!     ! Hopefully by now we've turned species_name into one of the species
!     if ( DEEBUG ) &
!       & print *, 'full name: ', trim(name), ' Species: ', trim(species_name)
!     call GetHashElement (lowercase(Species), &
!       & SpUniqueFieldDefinition, trim(species_name), &
!       & abbr_uniq_fdef, .false.)
!     call GetHashElement (UniqueFieldDefKeys, &
!       & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
!       & expnd_uniq_fdef, .false.)
!     if ( expnd_uniq_fdef == '' .or. expnd_uniq_fdef == ',' ) &
!       & expnd_uniq_fdef = 'MLS-Specific'
!     select case (trim(lowercase(species_name)))
!     case ('temperature')
!       units_name = 'K'
!     case ('gph')
!       units_name = 'm'
!     case ('geoheight')
!       units_name = 'km'
!     case ('rhi')
!       units_name = '%rhi'
!     case ('iwc', 'iwp')
!       units_name = 'g/m^3'
!     case default
!       units_name = 'vmr'
!     end select
!     if ( isColumnAmt ) then
!       ! units_name = 'DU'
!       if ( DEEBUG ) &
!         & call dump( .true., col_species_keys, col_species_hash, &
!         & 'column species units' )
!       call GetHashElement (col_species_keys, &
!       & col_species_hash, trim(lowercase(species_name)), &
!       & units_name, .true.)
!       if ( DEEBUG ) &
!       & print *, 'units_name: ', trim(units_name)
!       if ( units_name == ',' ) units_name = 'molcm2'
!     endif
!     if ( isTPPressure ) units_name = 'hPa'
!     if ( DEEBUG ) print *, 'Title ', trim(field_name)
!     status = mls_swwrlattr(swid, 'L2gpValue', 'Title', &
!       & MLS_CHARTYPE, 1, field_name)
!     if ( DEEBUG ) print *, 'Units ', trim(units_name)
!     status = mls_swwrlattr(swid, 'L2gpValue', 'Units', &
!       & MLS_CHARTYPE, 1, units_name)
!     status = he5_swwrlattr(swid, 'L2gpValue', 'MissingValue', &
!       & rgp_type, hsize(1), (/ real(l2gp%MissingL2GP, rgp) /) )
!     if ( DEEBUG ) print *, 'Title ', trim(expnd_uniq_fdef)
!     status = mls_swwrlattr(swid, 'L2gpValue', &
!       & 'UniqueFieldDefinition', &
!       & MLS_CHARTYPE, 1, expnd_uniq_fdef)
!     status = mls_swwrlattr(swid, 'L2gpPrecision', 'Title', &
!       & MLS_CHARTYPE, 1, trim(field_name)//'Precision')
!     status = mls_swwrlattr(swid, 'L2gpPrecision', 'Units', &
!       & MLS_CHARTYPE, 1, units_name)
!     status = he5_swwrlattr(swid, 'L2gpPrecision', 'MissingValue', &
!       & rgp_type, hsize(1), (/ real(l2gp%MissingValue, rgp) /) )
!     status = mls_swwrlattr(swid, 'L2gpPrecision', &
!       & 'UniqueFieldDefinition', &
!       & MLS_CHARTYPE, 1, expnd_uniq_fdef)
! 
!     ! ('Status' data field newly written)
!     status = mls_swwrlattr(swid, 'Status', 'Title', &
!       & MLS_CHARTYPE, 1, trim(field_name)//'Status')
!     status = mls_swwrlattr(swid, 'Status', 'Units', &
!       & MLS_CHARTYPE, 1, NOUNITS)
!     status = he5_swwrlattr(swid, 'Status', 'MissingValue', &
!       & HE5T_NATIVE_INT, hsize(1), (/ l2gp%MissingStatus /) )
!     status = mls_swwrlattr(swid, 'Status', &
!       & 'UniqueFieldDefinition', &
!       & MLS_CHARTYPE, 1, 'MLS-Specific')
!     
!     status = mls_swwrlattr(swid, 'Quality', 'Title', &
!       & MLS_CHARTYPE, 1, trim(field_name)//'Quality')
!     status = mls_swwrlattr(swid, 'Quality', 'Units', &
!       & MLS_CHARTYPE, 1, NOUNITS)
!     status = he5_swwrlattr(swid, 'Quality', 'MissingValue', &
!       & rgp_type, hsize(1), (/ real(l2gp%MissingValue, rgp) /) )
!     status = mls_swwrlattr(swid, 'Quality', &
!       & 'UniqueFieldDefinition', &
!       & MLS_CHARTYPE, 1, 'MLS-Specific')
!     
!     status = mls_swwrlattr(swid, 'Convergence', 'Title', &
!       & MLS_CHARTYPE, 1, trim(field_name)//'Convergence')
!     status = mls_swwrlattr(swid, 'Convergence', 'Units', &
!       & MLS_CHARTYPE, 1, NOUNITS)
!     status = he5_swwrlattr(swid, 'Convergence', 'MissingValue', &
!       & rgp_type, hsize(1), (/ real(l2gp%MissingValue, rgp) /) )
!     status = mls_swwrlattr(swid, 'Convergence', &
!       & 'UniqueFieldDefinition', &
!       & MLS_CHARTYPE, 1, 'MLS-Specific')
!     if ( AscDescModeIsField ) then
!       status = mls_swwrlattr(swid, 'AscDescMode', 'Title', &
!         & MLS_CHARTYPE, 1, trim(field_name)//'AscDescMode ')
!       status = mls_swwrlattr(swid, 'AscDescMode', 'Units', &
!         & MLS_CHARTYPE, 1, NOUNITS)
!       status = he5_swwrlattr(swid, 'AscDescMode', 'MissingValue', &
!         & HE5T_NATIVE_INT, hsize(1), (/ 0 /) )
!       status = mls_swwrlattr(swid, 'AscDescMode', &
!         & 'UniqueFieldDefinition', &
!         & MLS_CHARTYPE, 1, 'MLS-Specific')
!       endif
!     
    status = mls_SWdetach( swid )
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from swath interface', MLSFile=L2GPFile )
    end if


  !-------------------------------------
  end subroutine OutputL2GP_attributes_MF
  !-------------------------------------

  ! ----------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module NCL2GP

! $Log$
! Revision 1.8  2023/09/15 16:39:33  pwagner
! ReadNCGlobalAttr now public
!
! Revision 1.7  2022/12/22 22:56:41  pwagner
! Added two new Dumping subroutines; fixed bugs
!
! Revision 1.6  2022/12/08 18:56:25  pwagner
! Fixed typo; bypass ReadNCL2GPData_fileID
!
! Revision 1.5  2022/11/16 23:26:37  pwagner
! Fixed bugs in ReadNCL2GPData_MF_NC4
!
! Revision 1.4  2022/09/02 21:25:20  pwagner
! Fixed some but not all bugs in reading
!
! Revision 1.3  2022/02/03 18:43:15  pwagner
! Reduce the amount debug printing
!
! Revision 1.2  2020/03/19 22:35:23  pwagner
! Repaired more bugs than you can shake a stick at.
!
! Revision 1.1  2020/03/06 00:24:19  pwagner
! First commit
!
