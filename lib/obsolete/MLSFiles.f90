! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
module MLSFiles               ! Utility file routines
!===============================================================================
  use HDF, only: sfstart, sfend
  use HDFEOS, only: gdclose, gdopen, swclose, swopen, swinqswath
  use machine, only: io_error
  use MLSCommon, only: i4, NameLen, BareFNLen
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSStrings, only: Capitalize, LowerCase, Reverse, SortArray
  use output_m, only: blanks, output
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS, Pgs_smf_getMsg, &
    & PGSd_IO_Gen_RSeqFrm, PGSd_IO_Gen_RSeqUnf, & 
    & PGSd_IO_Gen_RDirFrm, PGSd_IO_Gen_RDirUnf, & 
    & PGSd_IO_Gen_WSeqFrm, PGSd_IO_Gen_WSeqUnf, & 
    & PGSd_IO_Gen_WDirFrm, PGSd_IO_Gen_WDirUnf, & 
    & PGSd_IO_Gen_USeqFrm, PGSd_IO_Gen_USeqUnf, & 
    & PGSd_IO_Gen_UDirFrm, PGSd_IO_Gen_UDirUnf, & 
    & PGSd_IO_Gen_ASeqFrm, PGSd_IO_Gen_ASeqUnf, &
    & PGS_IO_GEN_CloseF, PGS_IO_GEN_OpenF, PGSd_PC_FILE_PATH_MAX, &
    & UseSDPToolkit
  implicit none

  private 

  PUBLIC :: GetPCFromRef, get_free_lun, mls_io_gen_openF, &
  & mls_io_gen_closeF, split_path_name, &
  & mls_hdf_version, mls_inqswath, mls_sfstart, mls_sfend

  !------------------- RCS Ident Info -----------------------
  character(LEN=130) :: Id = &
    "$Id$"
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !----------------------------------------------------------

!     c o n t e n t s
!     - - - - - - - -

! GetPCFromRef       Turns a FileName into the corresponding PC
! get_free_lun       Gets a free logical unit number
! mls_hdf_version    Returns one of 'hdf4', 'hdf5', or '????'
! mls_inqswath       A wrapper for doing swingswath for versions 4 and 5
! mls_io_gen_openF   Opens a generic file using either the toolbox or else a Fortran OPEN statement
! mls_io_gen_closeF  Closes a generic file using either the toolbox or else a Fortran OPEN statement
! mls_sfstart        Opens an hdf file for writing metadata
! mls_sfend          Closes a file opened by mls_sfstart
! split_path_name    splits the input full_file_name into its components path and name
! hdf2hdf5_fileaccess
!                    Translates version 4 hdf access codes to version 5

  ! ((( Not yet ready for hdfeos5 versions of files )))
  !       The latter will incorporate the
  !       routines in he5lib

   ! Assume hdf files w/o explicit hdfVersion field are this
   ! 4 corresponds to hdf4, 5 to hdf5 in L2GP, L2AUX, etc.
   integer, parameter, public :: HDFVERSION_4 = 4
   integer, parameter, public :: HDFVERSION_5 = 5
   integer, parameter         :: DEFAULT_HDFVERSION = HDFVERSION_4
   
  ! Given this hdfVersion, try to autorecognize hdfversion
  ! then perform appropriate version of open/close; i.e., forgiving
  ! (must *not* be 4 or 5)
  integer, parameter, public :: WILDCARDHDFVERSION=HDFVERSION_4+HDFVERSION_5

  ! This isn't NameLen because it may have a path prefixed
  integer, parameter :: MAXFILENAMELENGTH=PGSd_PC_FILE_PATH_MAX

  ! These are error codes that may be returned by GetPCFromRef

  integer, parameter, public :: NAMENOTFOUND=-1
  integer, parameter, public :: INVALIDPCRANGE=NAMENOTFOUND-1
  integer, parameter, public :: CANTALLOCATENAMEARRAY=INVALIDPCRANGE-1

  ! These are error codes that may be returned by mls_io_gen_openF

  integer, parameter, public :: UNKNOWNFILEACCESSTYPE=-999
  integer, parameter, public :: UNKNOWNTOOLBOXMODE=UNKNOWNFILEACCESSTYPE+1
  integer, parameter, public :: NOFREEUNITS=UNKNOWNTOOLBOXMODE+1
  integer, parameter, public :: MUSTSUPPLYFILENAMEORPC=NOFREEUNITS+1
  integer, parameter, public :: NOPCIFNOTOOLKIT=MUSTSUPPLYFILENAMEORPC+1

  ! Now we have the legal unit numbers that files may be assigned

  integer, parameter :: bottom_unit_num=1
  integer, parameter :: top_unit_num=99

contains

  ! ---------------------------------------------  GetPCFromRef  -----

  ! This function takes a FileName as an arg and a range of PC numbers
  ! [PCBottom, PCTop] which are integers
  ! It returns thePC corresponding to the FileName

  ! FileName may be a fragment such as l2gp_temp of a longer name
  ! such as mls_l2_temp_v0.5_01-01-2004.dat

  ! If no matching file found, it sets ErrType=NAMENOTFOUND
  ! otherwise ErrType=0

  ! This is useful because all the Toolbox routines refer to files
  ! by their PC numbers, not their names

  ! Optionally you may require the match to be case-sensitive
  !   (by default it is not: l2_temp will match MLS_L2_TEMP_...)

  ! If you pass in a path, it will require that the paths also match

  ! Optionally you may pass in a version number and a debug flag

  ! optionally returns the exact name of the matching file

  function GetPCFromRef(FileName, PCBottom, PCTop, &
    & caseSensitive, ErrType, versionNum, debugOption, path, ExactName) &
    & result (thePC)

    ! Dummy arguments
    character (LEN=*), intent(IN)            :: FileName
    integer(i4),  intent(IN)                 :: PCBottom, PCTop
    integer(i4)                              :: thePC, notThePC
    integer(i4),  intent(OUT)                :: ErrType
    logical,  intent(IN)                     :: caseSensitive
    integer(i4),  optional                   :: versionNum
    logical,  optional, intent(IN)           :: debugOption
    character (LEN=*),  optional, intent(IN) :: path
    character (LEN=*), optional, intent(out) :: ExactName

    ! Local variables
    character (LEN=BareFNLen), dimension(:), allocatable &
     &                                :: nameArray
    integer, dimension(:), allocatable &
     &                                :: intArray
    character (LEN=*), parameter      :: UNASSIGNEDFILENAME = '*'
    character (LEN=MAXFILENAMELENGTH) :: MatchName, TryName, NameOnly
    character (LEN=MAXFILENAMELENGTH) :: PhysicalName, MatchPath
	 integer                       ::     version, returnStatus
    integer                       ::     numberPCs
    logical ::                            debug

   if(.not. UseSDPToolkit) then
      ErrType = NOPCIFNOTOOLKIT
      thePC = 0
      return
   endif

    if(present(debugOption)) then
      debug = debugOption
   else
      debug = .false.
   endif
   
    if(debug) then
      call output('get pc from ref', advance='yes')
      call output('FileName: ' // trim(FileName), advance='yes')
      call output('lower PCF limit: ' )
      call output(PCBottom, advance='yes')
      call output('upper PCF limit: ' )
      call output(PCTop, advance='yes')
      call output('case sensitive?: ' )
      call output(caseSensitive, advance='yes')
      call output('version number: ' )
      call output(versionNum, advance='yes')
    endif

    thePC = 0
    if(PCTop < PCBottom) then
      ErrType = INVALIDPCRANGE
      return
    endif

    if(caseSensitive) then
      MatchName = FileName
    else
      MatchName = Capitalize(FileName)
    endif

    if(debug) then
      call output('getting ref from pc:', advance='no')
    endif

    numberPCs = PCTop - PCBottom + 1
    Allocate(nameArray(numberPCs), intArray(numberPCs), STAT=ErrType)
    if ( ErrType /= 0 ) then
      ErrType = CANTALLOCATENAMEARRAY
      return
    endif

    ErrType = NAMENOTFOUND
    nameArray = UNASSIGNEDFILENAME
    do thePC = PCBottom, PCTop

      if(present(versionNum)) then
        version = versionNum
      else
        version = 1
      endif

      returnStatus = Pgs_pc_getReference(thePC, version, &
        & PhysicalName)

      if ( returnStatus == PGS_S_SUCCESS ) then

        if(.not. caseSensitive) then
          TryName = Capitalize(PhysicalName)
        else
          TryName = PhysicalName
        endif
        
        call split_path_name(TryName, MatchPath, NameOnly)
        nameArray(thePC-PCBottom+1) = NameOnly
      endif

    enddo
    ! Sort the file names from short to long
    ! to prevent unwanted matches between "O3" and "HNO3"
    call SortArray(nameArray, intArray, caseSensitive, &
     & sortedArray=nameArray, shorterFirst=.true., leftRight='r')
    do notThePC = 1, numberPCs
      thePC = intArray(notThePC) + PCBottom - 1         
      NameOnly = nameArray(notThePC)            
      if ( index(NameOnly, trim(MatchName)) /= 0 )then  
        ErrType = 0                                     
        exit                                            
      endif                                             
    enddo

    if(present(versionNum)) then  
      version = versionNum        
    else                          
      version = 1                 
    endif                         
    returnStatus = Pgs_pc_getReference(thePC, version, &
        & PhysicalName)
    if ( returnStatus == PGS_S_SUCCESS ) then             

      if(.not. caseSensitive) then                        
        TryName = Capitalize(PhysicalName)                
      else                                                
        TryName = PhysicalName                            
      endif                                               
                                                          
      call split_path_name(TryName, MatchPath, NameOnly)
    else
      ErrType = NAMENOTFOUND
    endif  
    if(present(path) .and. ErrType == 0) then
        if ( index(MatchPath, trim(path)) == 0 )then
          ErrType = NAMENOTFOUND
        endif
    endif

    if(present(ExactName) .and. ErrType == 0) then
      ExactName = PhysicalName
    endif

    Deallocate(nameArray, intArray)
  end function GetPCFromRef

! ---------------------------------------------- get_free_lun ------

! This function returns a free logical unit number

  INTEGER(i4) FUNCTION get_free_lun()
  LOGICAL :: exist                    ! Flag from inquire
  LOGICAL :: opened                   ! Flag from inquire
  DO get_free_lun = bottom_unit_num, top_unit_num
    INQUIRE(UNIT=get_free_lun, EXIST=exist, OPENED=opened)
    IF(exist .and. .not. opened) EXIT
  END DO
  IF (opened .or. .not. exist) CALL MLSMessage ( MLSMSG_Error, moduleName,  &
     "No logical unit numbers available" )
  END FUNCTION get_free_lun

  ! ---------------------------------------------  mls_io_gen_openF  -----

  ! This function opens a generic file using either the toolbox
  ! or else a Fortran OPEN statement
  ! according to toolbox_mode:
  ! 'sw' for swath files opened with swopen
  ! 'gd' for grid files opened with gdopen
  ! 'pg' for generic files opened with pgs_io_gen_openF
  ! 'op' for l3ascii files opened with simple fortran 'open'

  ! It returns theFileHandle corresponding to the FileName or the PC

  ! If given a FileName as an arg and a range of PC numbers
  ! [PCBottom, PCTop] which are integers
  ! it will attempt to find a corresponding PC

  ! If given a PC it will attempt to find the corresponding FileName

  ! toolbox_mode                  meaning
  ! PGS_IO_GEN_OpenF              use PGS_IO_Gen_OpenF or fail
  !      swopen                   use swopen or fail
  !      gdopen                   use gdopen or fail
  !       open                    use Fortran or fail

  ! If the FileName is not found, it sets ErrType=NAMENOTFOUND
  ! otherwise ErrType=0

  ! This is useful because all the Toolbox routines refer to files
  ! by their PC numbers, not their names

  function mls_io_gen_openF(toolbox_mode, caseSensitive, ErrType, &
    & record_length, FileAccessType, &
    & FileName, PCBottom, PCTop, versionNum, unknown, thePC, &
    & hdfVersion, debugOption) &
    &  result (theFileHandle)

    ! Dummy arguments
    integer(i4),  intent(OUT)  :: ErrType
    integer(i4),  intent(OUT)  :: record_length
    logical,  intent(IN)       :: caseSensitive
    character (LEN=*), intent(IN) :: toolbox_mode
    integer(i4), intent(IN)       :: FileAccessType
    integer(i4)  :: theFileHandle
    character (LEN=*), optional, intent(IN)   :: FileName
    integer(i4),  optional, intent(IN)   :: PCBottom, PCTop
    integer(i4), optional, intent(IN)                :: thePC
    integer(i4),  optional     :: versionNum
    logical, optional, intent(in) :: unknown
    logical, optional, intent(in) :: debugOption

    integer, optional, intent(in) :: hdfVersion

    ! Local
    integer :: myhdfVersion

    logical, parameter :: DEFAULT_PRINT_EVERY_OPEN=.false.
    integer, parameter :: FH_ON_ERROR=-99
    integer, parameter :: DEFAULTRECLEN=0
    integer, parameter :: KEYWORDLEN=12			! Max length of keywords in OPEN(...)
    character (LEN=MAXFILENAMELENGTH) :: myName
    integer(i4) :: myPC
    integer                       :: version, returnStatus, your_version
    logical       :: tiedup
    character (LEN=KEYWORDLEN) :: access, action, form, position, status
    character (LEN=2) :: the_eff_mode
    integer                       :: unit
    logical ::                            debug
    logical ::                            PRINT_EVERY_OPEN

    ! begin

    if(present(debugOption)) then
      debug = debugOption
   else
      debug = .false.
   endif

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = DEFAULT_HDFVERSION
    endif
   PRINT_EVERY_OPEN = DEFAULT_PRINT_EVERY_OPEN .or. debug
    ! In case of premature return
    theFileHandle = FH_ON_ERROR
    record_length = DEFAULTRECLEN
    myPC = FH_ON_ERROR
    returnStatus = 0           ! In case using Toolkit but supplied FileName

    if(present(versionNum)) then
      version = versionNum
    else
      version = 1
    endif

   the_eff_mode = LowerCase(toolbox_mode(1:2))
   
   if(UseSDPToolkit) then
    your_version = version
   ! Using Toolkit
    if(present(thePC)) then
      myPC = thePC
      returnStatus = Pgs_pc_getReference(thePC, your_version, &
        & myName)
      if ( debug ) then
          call output('Call to Pgs_pc_getReference', &
          & advance='yes')
          call output('returnStatus: ', advance='no')
          call blanks(2)
          call output(returnStatus, advance='yes')
          call output('thePC: ', advance='no')
          call blanks(2)
          call output(thePC, advance='yes')
          call output('your_version: ', advance='no')
          call blanks(2)
          call output(your_version, advance='yes')
      endif
    elseif(present(FileName)) then
      myName = FileName
      if(LowerCase(toolbox_mode(1:2)) == 'pg') then
        myPC = GetPCFromRef(trim(FileName), PCBottom, PCTop, &
          &	 caseSensitive, returnStatus, your_version, &
          & debugOption=debugOption)
      endif

    else
      ErrType = MUSTSUPPLYFILENAMEORPC
      if ( debug ) & 
       & call output('Must supply file name or pc to mls_io_gen_openF', &
       & advance='yes')
      return
    endif

   ! Not Using Toolkit
   ! Must supply FileName, use generic Fortran open
   elseif(.not. present(FileName)) then
      ErrType = NOPCIFNOTOOLKIT
      if ( debug ) & 
       & call output('No toolkit: must supply file name to mls_io_gen_openF', &
       & advance='yes')
      return
      
   else
      myName = FileName
      if(the_eff_mode == 'pg') the_eff_mode = 'op'
   endif

   if ( debug ) then
       call output('Arguments and options in call to mls_io_gen_openF', &
       & advance='yes')
       call output('Mode: ', advance='no')
       call blanks(2)
       call output(the_eff_mode, advance='yes')
       call output('File Name: ', advance='no')
       call blanks(2)
       call output(trim(myName), advance='yes')
       if(UseSDPToolkit) then
          call output('PCF-supplied number: ', advance='no')
          call blanks(2)
          call output(myPC, advance='yes')
       endif
       call output('Case sensitive? ', advance='no')
       call blanks(2)
       call output(caseSensitive, advance='yes')
       call output('File access type ', advance='no')
       call blanks(2)
       call output(FileAccessType, advance='yes')
   endif

    your_version = version
    select case (the_eff_mode)

    case('pg')
      if(returnStatus == 0) then
         ErrType = PGS_IO_Gen_OpenF(myPC, FileAccessType, record_length, &
          & theFileHandle, your_version)
      else
        ErrType = returnStatus
       endif

    case('sw')
      if(returnStatus == 0) then
        theFileHandle = swopen(trim(myName), FileAccessType)
        if(theFileHandle <= 0) then
          ErrType = min(theFileHandle, -1)
        else
          ErrType = 0
        endif
      else
        ErrType = returnStatus
      endif

    case('gd')
      if(returnStatus == 0) then
        theFileHandle = gdopen(trim(myName), FileAccessType)
        if(theFileHandle <= 0) then
          ErrType = min(theFileHandle, -1)
        else
          ErrType = 0
        endif
      else
        ErrType = returnStatus
      endif

    case('op')
      theFileHandle = FH_ON_ERROR
      if(FileAccessType == PGSd_IO_Gen_RSeqFrm) then
        status = 'old'
        access = 'sequential'
        form = 'formatted'
        action = 'read'
        position = 'rewind'
      elseif(FileAccessType == PGSd_IO_Gen_RSeqUnf) then
        status = 'old'
        access = 'sequential'
        form = 'unformatted'
        action = 'read'
        position = 'rewind'
      elseif(FileAccessType == PGSd_IO_Gen_RDirFrm) then
        status = 'old'
        access = 'direct'
        form = 'formatted'
        action = 'read'
        position = 'rewind'
      elseif(FileAccessType == PGSd_IO_Gen_RDirUnf) then
        status = 'old'
        access = 'direct'
        form = 'unformatted'
        action = 'read'
        position = 'rewind'

      elseif(FileAccessType == PGSd_IO_Gen_WSeqFrm) then
        status = 'new'
        access = 'sequential'
        form = 'formatted'
        action = 'write'
        position = 'rewind'
      elseif(FileAccessType == PGSd_IO_Gen_WSeqUnf) then
        status = 'new'
        access = 'sequential'
        form = 'unformatted'
        action = 'write'
        position = 'rewind'
      elseif(FileAccessType == PGSd_IO_Gen_WDirFrm) then
        status = 'new'
        access = 'direct'
        form = 'formatted'
        action = 'write'
        position = 'rewind'
      elseif(FileAccessType == PGSd_IO_Gen_WDirUnf) then
        status = 'new'
        access = 'direct'
        form = 'unformatted'
        action = 'write'
        position = 'rewind'

      elseif(FileAccessType == PGSd_IO_Gen_USeqFrm) then
        status = 'old'
        access = 'sequential'
        form = 'formatted'
        action = 'readwrite'
        position = 'rewind'
      elseif(FileAccessType == PGSd_IO_Gen_USeqUnf) then
        status = 'old'
        access = 'sequential'
        form = 'unformatted'
        action = 'readwrite'
        position = 'rewind'
      elseif(FileAccessType == PGSd_IO_Gen_UDirFrm) then
        status = 'old'
        access = 'direct'
        form = 'formatted'
        action = 'readwrite'
        position = 'rewind'
      elseif(FileAccessType == PGSd_IO_Gen_UDirUnf) then
        status = 'old'
        access = 'direct'
        form = 'unformatted'
        action = 'readwrite'
        position = 'rewind'

      elseif(FileAccessType == PGSd_IO_Gen_ASeqFrm) then
        status = 'old'
        access = 'sequential'
        form = 'formatted'
        action = 'readwrite'
        position = 'append'
      elseif(FileAccessType == PGSd_IO_Gen_ASeqUnf) then
        status = 'old'
        access = 'sequential'
        form = 'unformatted'
        action = 'readwrite'
        position = 'append'

      else
        ErrType = UNKNOWNFILEACCESSTYPE
        if ( debug ) &
         & call output('Unknown file access type', advance='yes')
        return
      endif

      tiedup = .true.

      do unit = bottom_unit_num, top_unit_num
      	inquire ( unit=unit, opened=tiedup )
        if (.not. tiedup) then
          exit
        endif
      enddo

      if(tiedup) then
        ErrType = NOFREEUNITS
        if ( debug ) &
         & call output('No free io units available', advance='yes')
        return
      endif

      if (present(unknown)) then
        if (unknown) status = 'unknown'
      end if

      if(access /= 'direct') then
        open(unit=unit, access=access, action=action, form=form, &
          & position=position, status=status, file=trim(myName), iostat=ErrType)
      else
        open(unit=unit, access=access, action=action, form=form, &
          & status=status, file=trim(myName), iostat=ErrType)
      endif

      if(ErrType /= 0 .or. PRINT_EVERY_OPEN) then
        call output( 'Fortran opening unit ', advance='no')
        call output(  unit, advance='yes')
        call output( 'access ' // access, advance='yes')
        call output( 'action ' // action, advance='yes')
        call output( 'form ' // form, advance='yes')
        call output( 'position ' // position, advance='yes')
        call output( 'status ' // status, advance='yes')
        call output( 'file ' // trim(myName), advance='yes')
      endif

      if(ErrType /= 0) then
        call output( 'iostat ', advance='no')
        call output(  ErrType, advance='yes')
        call io_error('io error in MLSFiles: mls_io_gen_openF' // &
          & ' Fortran open', ErrType, trim(myName))
      else
        theFileHandle = unit
      endif

    case default
      ErrType = UNKNOWNTOOLBOXMODE

    end select

    if(ErrType /= 0) then
      theFileHandle = FH_ON_ERROR
    endif
    
    if( debug ) then
       call output('Error Type (0 means none): ', advance='no')
       call blanks(2)
       call output(ErrType, advance='yes')
       call output('record_length: ', advance='no')
       call blanks(2)
       call output(record_length, advance='yes')
       call output('theFileHandle: ', advance='no')
       call blanks(2)
       call output(theFileHandle, advance='yes')
    endif

  end function mls_io_gen_openF

  ! ---------------------------------------------  mls_io_gen_closeF  -----

  ! This function closes a generic file using either the toolbox
  ! or else a Fortran CLOSE statement
  ! according to toolbox_mode:
  ! 'sw' for swath files opened with swopen
  ! 'gd' for grid files opened with gdopen
  ! 'pg' for generic files opened with pgs_io_gen_openF
  ! 'cl' for l3ascii files opened with simple fortran 'open'

  ! It returns a non-zero error status only if unsuccessful

  ! If must be given a FileHandle as an arg
  ! (A later version may allow choice between file handle and file name)
  ! (This version only uses a file name in autodetecting its hdf version)
  function mls_io_gen_closeF(toolbox_mode, theFileHandle, &
    & FileName, hdfVersion) &
    &  result (ErrType)

    ! Dummy arguments
    integer(i4)  :: ErrType
    integer(i4), intent(IN)  :: theFileHandle
    character (LEN=*), intent(IN)   :: toolbox_mode
    character (LEN=*), optional, intent(IN)   :: FileName

    integer, optional, intent(in) :: hdfVersion

    ! Local
    integer :: myhdfVersion

    logical, parameter :: PRINT_EVERY_CLOSE=.false.
    character (LEN=2) :: the_eff_mode

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = DEFAULT_HDFVERSION
    endif

   if(UseSDPToolkit) then
   ! Using Toolkit
    the_eff_mode = LowerCase(toolbox_mode(1:2))
   ! Not Using Toolkit
   else
      the_eff_mode = LowerCase(toolbox_mode(1:2))
      if(the_eff_mode == 'pg') the_eff_mode = 'cl'
   endif

    select case (the_eff_mode)

    case('pg')
      ErrType = PGS_IO_Gen_CLoseF(theFileHandle)

    case('sw')
      ErrType = swclose(theFileHandle)

    case('gd')
      ErrType = gdclose(theFileHandle)

    case('cl')
      close(unit=theFileHandle, iostat=ErrType)		

      if(ErrType /= 0 .or. PRINT_EVERY_CLOSE) then
        call output( 'Fortran closing unit ', advance='no')
        call output(  theFileHandle, advance='yes')
      endif

      if(ErrType /= 0) then
        call output( 'iostat ', advance='no')
        call output(  ErrType, advance='yes')
        call io_error('io error in MLSFiles: mls_io_gen_closeF' // &
          & ' Fortran close', ErrType, 'unknown')
      endif

    case default
      ErrType = UNKNOWNTOOLBOXMODE

    end select

  end function mls_io_gen_closeF

  ! ---------------------------------------------  split_path_name  -----

  ! This routine splits the input full_file_name
  ! into its components path and name
  ! where path may include one or more "/" or slash elements
  ! (but one must be the terminating one; e.g., 'System/')
  ! while name must have none
  ! special cases by example: full_file_name -> (path, name)
  ! look.ma.no.slash -> (' ', 'look.ma.no.slash')
  ! Luke/I/am/your/father/ -> ('Luke/I/am/your/father/', ' ')

  ! optionally you may supply the slash divider
  ! which must be a single character

  subroutine split_path_name(full_file_name, path, name, slash)

    ! Arguments

    character (len=*), intent(in) :: full_file_name
    character (len=*), intent(out) :: path
    character (len=*), intent(out) :: name
    character (len=1), optional, intent(in) :: slash

    ! Local

    character (len=1) :: mySlash
    character (len=MAXFILENAMELENGTH) :: mirrored_ffn
    integer :: loc
!   logical, parameter :: DEBUG = .false.

    ! Begin

    if(present(slash)) then
      mySlash = slash
    else
      mySlash = '/'
    endif

    if(len(full_file_name) <= 0) then
      path = ' '
      name = ' '
      return
    endif

    mirrored_ffn = Reverse(full_file_name)
    loc = index(mirrored_ffn, mySlash)
    

    if(loc <= 0) then
      path = ' '
      name = adjustl(full_file_name)
    elseif(loc == 1) then
      path = adjustl(full_file_name)
      name = ' '
    else
      path = adjustl(Reverse(mirrored_ffn(loc:)))
      name = adjustl(Reverse(mirrored_ffn(:loc-1)))
    endif

  end subroutine split_path_name

  ! ---------------------------------------------  mls_inqswath  -----

  ! This function acts as a wrapper to allow hdf5 or hdf4 routines to be called

  function mls_inqswath(FileName, swathList, strBufSize, hdfVersion)

    ! Arguments

      character (len=*), intent(in) :: FILENAME
      character (len=*), intent(out) :: SWATHLIST
      integer, intent(out):: STRBUFSIZE
      integer :: mls_inqswath

    integer, optional, intent(in) :: hdfVersion

    ! Local
    integer :: myhdfVersion

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = DEFAULT_HDFVERSION
    endif
    ! begin
    mls_inqswath = swinqswath(FileName, swathList, strBufSize)

  end function mls_inqswath

  ! ---------------------------------------------  mls_sfstart  -----

  ! This function acts as a wrapper to allow hdf5 or hdf4 routines to be called

  function mls_sfstart(FileName, FileAccess, hdfVersion)

    ! Arguments

      character (len=*), intent(in) :: FILENAME
    integer(i4), intent(IN)       :: FileAccess
    integer                       :: mls_sfstart

    integer, optional, intent(in) :: hdfVersion

    ! Local
    integer :: myhdfVersion

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = DEFAULT_HDFVERSION
    endif
    ! begin
    mls_sfstart = sfstart(FileName, FileAccess)

  end function mls_sfstart

  ! ---------------------------------------------  mls_sfend  -----

  ! This function acts as a wrapper to allow hdf5 or hdf4 routines to be called

  function mls_sfend(sdid, hdfVersion)

    ! Arguments

      integer, intent(IN)       :: sdid
      integer :: mls_sfend

    integer, optional, intent(in) :: hdfVersion

    ! Local
    integer :: myhdfVersion

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = DEFAULT_HDFVERSION
    endif
    ! begin
    mls_sfend = sfend(sdid)

  end function mls_sfend

  ! ---------------------------------------------  mls_hdf_version  -----

  ! This function returns a character string depending on hdf version:
  ! hdf version         returned value
  !    hdf4                 hdf4
  !    hdf5                 hdf5
  !  unknown                ????

  function mls_hdf_version(FileName, preferred_version, AccessType) &
   & result (hdf_version)
!  function mls_hdf_version(FileName)  result (hdf_version)

    ! Arguments

      character (len=*), intent(in) :: FILENAME
   !   character (len=4)             :: hdf_version
      integer(i4), optional, intent(in)  :: PREFERRED_VERSION
      integer(i4), optional, intent(in)  :: ACCESSTYPE
      integer(i4)                        :: HDF_VERSION

      integer :: returnStatus
      logical :: is_hdf5
    ! begin
!      hdf_version = 'hdf4'
      hdf_version = 4

  end function mls_hdf_version

!====================
end module MLSFiles
!====================

!
! $Log$
! Revision 2.32  2002/06/04 18:11:39  bill
! added a get free lun module--wgr
!
! Revision 2.31  2002/03/14 23:31:45  pwagner
! HDFVERSION_4 and 5 now public
!
! Revision 2.30  2002/02/19 23:12:40  pwagner
! Eliminated unwanted matches between o3 and hno3
!
! Revision 2.29  2002/01/31 00:35:56  pwagner
! Brought mls_io_gen_closeF into compilance with he5lib version
!
! Revision 2.28  2002/01/29 23:45:26  pwagner
! Added WILDCARDHDFVERSION as public param
!
! Revision 2.27  2002/01/29 00:47:41  pwagner
! Converted mls_hdf_version to integer function
!
! Revision 2.26  2002/01/23 21:47:31  pwagner
! Begun to make hdf5-capable; not yet, though
!
! Revision 2.25  2002/01/18 18:51:22  pwagner
! Fixed bug when calling mls_open w/o toolkit
!
! Revision 2.24  2002/01/18 00:53:37  pwagner
! Added inqswath, sfend, sfstart wrappers; not ready yet for hdf5
!
! Revision 2.23  2002/01/11 00:44:37  pwagner
! Fixed bug where ErrType was left unset
!
! Revision 2.22  2002/01/09 23:43:26  pwagner
! Added toc; removed debugging stuff
!
! Revision 2.21  2001/07/18 23:59:16  pwagner
! Fixed bug with version variable in mls_io_gen_openf
!
! Revision 2.20  2001/06/19 22:22:28  pwagner
! Set record_length to DEFAULTRECLEN unless otherwise
!
! Revision 2.19  2001/05/08 23:35:44  pwagner
! Fixed bug in the_eff_mode
!
! Revision 2.18  2001/05/08 21:36:38  livesey
! Changed PRINT_EVERY_OPEN to .false.
!
! Revision 2.17  2001/05/07 23:25:02  pwagner
! Detachable from toolkit
!
! Revision 2.16  2001/04/25 21:52:50  livesey
! Added optional `unknown' argument
!
! Revision 2.15  2001/04/25 20:32:29  livesey
! Bug fix, trim filename on output
!
! Revision 2.14  2001/04/17 23:44:54  pwagner
! Fixed bug with ExactName in getpc..
!
! Revision 2.13  2001/04/17 22:10:55  pwagner
! getpcfromref takes optional path; otherwise matches only filename
!
! Revision 2.12  2001/04/13 23:50:55  pwagner
! split_path_name now works properly
!
! Revision 2.11  2001/04/13 00:18:10  pwagner
! Prints only if debug present AND .TRUE.
!
! Revision 2.10  2001/04/10 20:05:07  livesey
! Tidied up
!
! Revision 2.9  2001/04/09 23:42:38  pwagner
! Resets version properly each time
!
! Revision 2.8  2001/04/07 00:16:28  pwagner
! fixed MAXFILENAMELENGTH to not core dump in toolbox
!
! Revision 2.7  2001/04/03 23:52:29  pwagner
! Added split_path_name
!
! Revision 2.6  2001/03/27 17:30:59  pwagner
! Added mls_io_gen_closeF
!
! Revision 2.5  2001/03/24 00:30:14  pwagner
! Now complains only if an error, and then via output
!
! Revision 2.4  2001/03/22 01:09:31  pwagner
! Added file name to Fortran open statement
!
! Revision 2.3  2001/03/21 00:48:43  pwagner
! Corrected mls_io_gen_openF
!
! Revision 2.2  2001/03/20 00:41:28  pwagner
! Added mls_io_gen_openF
!
! Revision 2.1  2001/03/07 01:02:37  pwagner
! First commit
!
