! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
module MLSFiles               ! Utility file routines
  !===============================================================================
  use Hdf, only: DFACC_CREATE, DFACC_READ, sfstart, sfend
  use Hdf5, only: h5fopen_f, h5fclose_f, h5fis_hdf5_f
  use Hdf5_params, only: H5F_ACC_RDONLY, H5F_ACC_RDWR, H5F_ACC_TRUNC
  use HDFEOS, only: gdclose, gdopen, swclose, swopen, swinqswath
  use HDFEOS5, only: he5_swclose, he5_swopen, he5_swinqswath, &
    & he5_gdopen, he5_gdclose
  use machine, only: io_error
  use MLSCommon, only: i4
  use MLSStrings, only: Capitalize, LowerCase, Reverse
  use output_m, only: blanks, output
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS, &
    & PGSd_IO_Gen_RSeqFrm, PGSd_IO_Gen_RSeqUnf, & 
    & PGSd_IO_Gen_RDirFrm, PGSd_IO_Gen_RDirUnf, & 
    & PGSd_IO_Gen_WSeqFrm, PGSd_IO_Gen_WSeqUnf, & 
    & PGSd_IO_Gen_WDirFrm, PGSd_IO_Gen_WDirUnf, & 
    & PGSd_IO_Gen_USeqFrm, PGSd_IO_Gen_USeqUnf, & 
    & PGSd_IO_Gen_UDirFrm, PGSd_IO_Gen_UDirUnf, & 
    & PGSd_IO_Gen_ASeqFrm, PGSd_IO_Gen_ASeqUnf, &
    & PGS_IO_GEN_CloseF, PGS_IO_GEN_OpenF, PGSd_PC_FILE_PATH_MAX, &
    & UseSDPToolkit
!   In the long run, we'll try putting interfaces to these in SDPToolkit.f90
!   Until then, just declare them as external
!    & PGS_MET_SFstart, PGS_MET_SFend, &
  implicit none

  private 

  public :: GetPCFromRef, mls_io_gen_openF, &
  & mls_io_gen_closeF, split_path_name, &
  & mls_hdf_version, mls_inqswath, mls_sfstart, mls_sfend

  !------------------- RCS Ident Info -----------------------
  character(LEN=130) :: Id = &
    "$Id$"
  !----------------------------------------------------------

!     c o n t e n t s
!     - - - - - - - -

! GetPCFromRef       Turns a FileName into the corresponding PC
! mls_hdf_version    Returns one of 'hdf4', 'hdf5', or '????'
! mls_inqswath       A wrapper for doing swingswath for versions 4 and 5
! mls_io_gen_openF   Opens a generic file using either the toolbox or else a Fortran OPEN statement
! mls_io_gen_closeF  Closes a generic file using either the toolbox or else a Fortran OPEN statement
! mls_sfstart        Opens an hdf file for writing metadata
! mls_sfend          Closes a file opened by mls_sfstart
! split_path_name    splits the input full_file_name into its components path and name
! hdf2hdf5_fileaccess
!                    Translates version 4 hdf access codes to version 5

   ! Assume hdf files w/o explicit hdfVersion field are this
   ! 4 corresponds to hdf4, 5 to hdf5 in L2GP, L2AUX, etc.
   integer, parameter :: HDFVERSION_4 = 4
   integer, parameter :: HDFVERSION_5 = 5
   integer, parameter :: DEFAULT_HDFVERSION = HDFVERSION_5
   
  ! Given this hdfVersion, try to autorecognize hdfversion
  ! then perform appropriate version of open/close; i.e., forgiving
  ! (must *not* be 4 or 5)
  integer, parameter, public :: WILDCARDHDFVERSION=HDFVERSION_4+HDFVERSION_5

  ! This isn't NameLen because it may have a path prefixed
  integer, parameter :: MAXFILENAMELENGTH=PGSd_PC_FILE_PATH_MAX

  ! These are error codes that may be returned by GetPCFromRef

  integer, parameter, public :: NAMENOTFOUND=-1
  integer, parameter, public :: INVALIDPCRANGE=NAMENOTFOUND-1

  ! These are error codes that may be returned by mls_io_gen_openF

  integer, parameter, public :: UNKNOWNFILEACCESSTYPE=-999
  integer, parameter, public :: UNKNOWNTOOLBOXMODE=UNKNOWNFILEACCESSTYPE+1
  integer, parameter, public :: NOFREEUNITS=UNKNOWNTOOLBOXMODE+1
  integer, parameter, public :: MUSTSUPPLYFILENAMEORPC=NOFREEUNITS+1
  integer, parameter, public :: NOPCIFNOTOOLKIT=MUSTSUPPLYFILENAMEORPC+1
  integer, parameter, public :: NOSUCHHDFVERSION=NOPCIFNOTOOLKIT+1

  ! These are error codes that may be returned by mls_io_gen_openF
  ! or by mls_hdf_version if you call it directly

  integer, parameter, public :: MUSTSUPPLYFILENAME=NOSUCHHDFVERSION+1
  integer, parameter, public :: ERRORINH5FFUNCTION=MUSTSUPPLYFILENAME+1
  integer, parameter, public :: WRONGHDFVERSION=ERRORINH5FFUNCTION+1
  
  ! Whether to use PGS_MET commands in mls_sf(start)(end)
  ! (The alternative is to use hdf5 calls directly)
  logical, parameter :: PGS_MET4MLS_SF = .true.
  
  ! The only legal unit numbers that files may be assigned
  ! for use by Fortran opens, closes, reads and writes
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
    integer(i4)                              :: thePC
    integer(i4),  intent(OUT)                :: ErrType
    logical,  intent(IN)                     :: caseSensitive
    integer(i4),  optional                   :: versionNum
    logical,  optional, intent(IN)           :: debugOption
    character (LEN=*),  optional, intent(IN)       :: path
    character (LEN=*), optional, intent(out) :: ExactName

    ! Local variables

    character (LEN=MAXFILENAMELENGTH) :: MatchName, TryName, NameOnly
    character (LEN=MAXFILENAMELENGTH) :: PhysicalName, MatchPath
    integer                       ::     version, returnStatus
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

    ErrType = NAMENOTFOUND

    if(debug) then
      call output('getting ref from pc:', advance='no')
    endif

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

        if ( index(NameOnly, trim(MatchName)) /= 0 )then
          ErrType = 0
          exit
        endif

      endif

    enddo

    if(present(path) .and. ErrType == 0) then
        if ( index(MatchPath, trim(path)) == 0 )then
          ErrType = NAMENOTFOUND
        endif
    endif

    if(present(ExactName) .and. ErrType == 0) then
      ExactName = PhysicalName
    endif

  end function GetPCFromRef

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

   PRINT_EVERY_OPEN = DEFAULT_PRINT_EVERY_OPEN .or. debug
   returnStatus = 0           ! In case using Toolkit but supplied FileName
    ! In case of premature return
    theFileHandle = FH_ON_ERROR
    record_length = DEFAULTRECLEN

    if(present(versionNum)) then
      version = versionNum
    else
      version = 1
    endif

   if(UseSDPToolkit) then
    the_eff_mode = LowerCase(toolbox_mode(1:2))
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
      the_eff_mode = LowerCase(toolbox_mode(1:2))
      if(the_eff_mode == 'pg') the_eff_mode = 'op'
   endif

   ! hdfVersion unimportant for 'pg' or 'op' operations
   if (the_eff_mode == 'pg' .or. the_eff_mode == 'op' ) then
     myhdfVersion = WILDCARDHDFVERSION
   else
     myhdfVersion = mls_hdf_version(trim(myName), hdfVersion, FileAccessType)             
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
       call output('PCF-supplied number: ', advance='no')
       call blanks(2)
       call output(myPC, advance='yes')
       call output('Case sensitive? ', advance='no')
       call blanks(2)
       call output(caseSensitive, advance='yes')
       call output('File access type ', advance='no')
       call blanks(2)
       call output(FileAccessType, advance='yes')
       call output('Return status ', advance='no')
       call blanks(2)
       call output(returnStatus, advance='yes')
   endif

   if ( myhdfVersion < 0 ) then
     ErrType = myhdfVersion
     return
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
      if(returnStatus /= 0) then
        ErrType = returnStatus
        return
      elseif(myhdfVersion == HDFVERSION_5) then
        theFileHandle = he5_swopen(trim(myName), &
          & hdf2hdf5_fileaccess(FileAccessType))
      elseif(myhdfVersion == HDFVERSION_4) then
        theFileHandle = swopen(trim(myName), FileAccessType)
      else
        ErrType = NOSUCHHDFVERSION
        return
      endif
      if(theFileHandle <= 0) then                                          
        ErrType = min(theFileHandle, -1)                                   
      else                                                                 
        ErrType = 0                                                        
      endif                                                                
      if ( debug ) then                                                    
          call output('Args and results from he5_swopen', &                
          & advance='yes')                                                 
          call output('File Name: ', advance='no')                         
          call blanks(2)                                                   
          call output(trim(myName), advance='yes')                         
          call output('File Access: ', advance='no')                       
          call blanks(2)                                                   
          call output(hdf2hdf5_fileaccess(FileAccessType), advance='yes')  
          call output('theFileHandle: ', advance='no')                     
          call blanks(2)                                                   
          call output(theFileHandle, advance='yes')                        
      endif                                                                

    case('gd')
      if(returnStatus /= 0) then
        ErrType = returnStatus
        return
      elseif(myhdfVersion == HDFVERSION_5) then
        theFileHandle = he5_gdopen(trim(myName), &
          & hdf2hdf5_fileaccess(FileAccessType))
      elseif(myhdfVersion == HDFVERSION_4) then
        theFileHandle = gdopen(trim(myName), FileAccessType)
      else
        ErrType = NOSUCHHDFVERSION
        return
      endif
      if(theFileHandle <= 0) then         
        ErrType = min(theFileHandle, -1)  
      else                                
        ErrType = 0                       
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

  ! It must be given a FileHandle as an arg
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

    ! begin

   if(UseSDPToolkit) then
   ! Using Toolkit
    the_eff_mode = LowerCase(toolbox_mode(1:2))
   ! Not Using Toolkit
   else
      the_eff_mode = LowerCase(toolbox_mode(1:2))
      if(the_eff_mode == 'pg') the_eff_mode = 'cl'
   endif

   ! hdfVersion unimportant for 'pg' or 'cl' operations
   if (the_eff_mode == 'pg' .or. the_eff_mode == 'cl' ) then
     ErrType = 0
   elseif ( present(FileName) ) then
     myhdfVersion = mls_hdf_version(trim(FileName), hdfVersion, DFACC_READ)
     if ( myhdfVersion < 0 ) then
        ErrType = myhdfVersion
        return
     endif
   elseif ( present(hdfVersion) ) then
     myhdfVersion = hdfVersion
   else
     myhdfVersion = DEFAULT_HDFVERSION
   endif
    select case (the_eff_mode)

    case('pg')
      ErrType = PGS_IO_Gen_CLoseF(theFileHandle)

    case('sw')
      if(myhdfVersion == HDFVERSION_5) then
        ErrType = he5_swclose(theFileHandle)
      elseif(myhdfVersion == HDFVERSION_4) then
        ErrType = swclose(theFileHandle)
      else
        ErrType = NOSUCHHDFVERSION
      endif

    case('gd')
      if(myhdfVersion == HDFVERSION_5) then
        ErrType = he5_gdclose(theFileHandle)
      elseif(myhdfVersion == HDFVERSION_4) then
        ErrType = gdclose(theFileHandle)
      else
        ErrType = NOSUCHHDFVERSION
      endif

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

  ! ---------------------------------------------  hdf2hdf5_fileaccess  -----

  ! This function converts hdf4 file access types to
  ! corresponding hdf5 numbers

  function hdf2hdf5_fileaccess(FileAccesshdf4) result (FileAccesshdf5)

    ! Arguments

    integer(i4), intent(IN)       :: FileAccesshdf4
    integer(i4)                   :: FileAccesshdf5
    
!    integer(i4), parameter        :: H5F_ACC_RDONLY = 0
!    integer(i4), parameter        :: H5F_ACC_RDWR   = 1
!    integer(i4), parameter        :: H5F_ACC_TRUNC  = 2

    ! begin
    select case (FileAccesshdf4)

    case(DFACC_CREATE)
      FileAccesshdf5 = H5F_ACC_TRUNC

    case(DFACC_READ)
      FileAccesshdf5 = H5F_ACC_RDONLY

    case default
      FileAccesshdf5 = H5F_ACC_RDWR

    end select

  end function hdf2hdf5_fileaccess

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
    if(myhdfVersion == HDFVERSION_5) then
      mls_inqswath = he5_swinqswath(trim(FileName), swathList, strBufSize)
    elseif(myhdfVersion == HDFVERSION_4) then
      mls_inqswath = swinqswath(FileName, swathList, strBufSize)
    else                          
      mls_inqswath = NOSUCHHDFVERSION  
    endif                         

  end function mls_inqswath

  ! ---------------------------------------------  mls_sfstart  -----

  ! This function acts as a wrapper to allow hdf5 or hdf4 routines to be called

  function mls_sfstart(FileName, FileAccess, hdfVersion)

    ! Arguments

      character (len=*), intent(in) :: FILENAME
    integer(i4), intent(IN)       :: FileAccess
    integer, optional, intent(in) :: hdfVersion
    integer                       :: mls_sfstart
    integer                       :: returnStatus
    integer                       :: nameLength
    integer                       :: access_prp_default
    integer                       :: myhdfVersion
    
    integer, parameter :: h5p_default_f = 0
    integer, external :: PGS_MET_SFstart

    ! begin
   if ( present(hdfVersion) ) then
     myhdfVersion = hdfVersion
   else
     myhdfVersion = DEFAULT_HDFVERSION
   endif
   if ( myhdfVersion == HDFVERSION_4) then
     mls_sfstart = sfstart (FileName, FileAccess)
     return
   endif
!    mls_sfstart = PGS_MET_SFstart(FileName, hdf2hdf5_fileaccess(FileAccess))
   if ( PGS_MET4MLS_SF ) then
     returnStatus = PGS_MET_SFstart(trim(FileName), &
      & hdf2hdf5_fileaccess(FileAccess), mls_sfstart)
   else
     access_prp_default = h5p_default_f
     call h5fopen_f(trim(FileName), hdf2hdf5_fileaccess(FileAccess), &
      & mls_sfstart, returnStatus, access_prp_default)
   endif
!     nameLength = LEN(FileName)
!     returnStatus = h5fopen_c(FileName, nameLength, &
!      & hdf2hdf5_fileaccess(FileAccess), access_prp_default, &
!      & mls_sfstart)
!     mls_sfstart he5_swopen(trim(
     if (returnStatus /= 0 ) then
       call output ('Try again--PGS_MET_SFstart still unhappy; returns ')
       call output (returnStatus, advance='yes')
       mls_sfstart = -1
      endif

  end function mls_sfstart

  ! ---------------------------------------------  mls_sfend  -----

  ! This function acts as a wrapper to allow hdf5 or hdf4 routines to be called

  function mls_sfend(sdid, hdfVersion)

    ! Arguments

    integer, intent(IN)       :: sdid  
    integer :: mls_sfend            
    integer, optional, intent(in) :: hdfVersion

    integer                       :: myhdfVersion
    integer, external :: PGS_MET_SFend

    ! begin
   if ( present(hdfVersion) ) then
     myhdfVersion = hdfVersion
   else
     myhdfVersion = DEFAULT_HDFVERSION
   endif
   if ( myhdfVersion == HDFVERSION_4) then
     mls_sfend = sfend (sdid)
     return
   endif
!    mls_sfend = h5fclose_c(sdid)
   if ( PGS_MET4MLS_SF ) then
     mls_sfend = PGS_MET_SFend(sdid)
   else
     call h5fclose_f(sdid, mls_sfend)
   endif
!    mls_sfend = 0

  end function mls_sfend

  ! ---------------------------------------------  mls_hdf_version  -----

  ! This function returns the hdf version of the file:
  ! hdf version         returned value   integer
  !    hdf4                 hdf4           4
  !    hdf5                 hdf5           5
  !  unknown                ????      NOSUCHHDFVERSION

  function mls_hdf_version(FileName, preferred_version, AccessType) &
   & result (hdf_version)

    ! Arguments

      character (len=*), intent(in)  :: FILENAME
   !   character (len=4)             :: hdf_version
      integer(i4), optional, intent(in)  :: PREFERRED_VERSION
      integer(i4), optional, intent(in)  :: ACCESSTYPE
      integer(i4)                        :: HDF_VERSION

      integer :: returnStatus
      integer :: myPreferred_Version
      integer :: myAccessType
      logical :: is_hdf5

    ! begin

    if (FileName == '') then
      HDF_VERSION = MUSTSUPPLYFILENAME
      return
    endif
    
    if ( present(preferred_version) ) then
      myPreferred_Version = preferred_version
    else
      myPreferred_Version = WILDCARDHDFVERSION
    endif
    if ( present(accesstype) ) then
      myAccessType = accesstype
    else
      myAccessType = DFACC_READ
    endif

   ! If the file is being newly created, no way to find its hdf version yet
    if ( myAccessType == DFACC_CREATE ) then
      hdf_version = myPreferred_Version
      return
    endif
    
    returnStatus = 0
    is_hdf5 = (DEFAULT_HDFVERSION == 5)
    call h5fis_hdf5_f(trim(FileName), is_hdf5, returnStatus)
    if ( returnStatus /= 0 ) then
!      hdf_version = '????'
      hdf_version = ERRORINH5FFUNCTION
    elseif ( is_hdf5 ) then
!      hdf_version = 'hdf5'
      hdf_version = 5
    else
!      hdf_version = 'hdf4'
      hdf_version = 4
    endif

   ! hdf_version ==  myPreferred_Version ?
   if ( myPreferred_Version /= WILDCARDHDFVERSION &
     & .and. &
     & hdf_version /= myPreferred_Version ) &
     & hdf_version = WRONGHDFVERSION

  end function mls_hdf_version

  !====================
end module MLSFiles
!====================

!
! $Log$
! Revision 1.7  2002/01/29 23:46:30  pwagner
! Added WILDCARDHDFVERSION as public param
!
! Revision 1.6  2002/01/29 00:48:43  pwagner
! Now should handle both hdfVersions; not tested yet
!
! Revision 1.5  2002/01/26 00:14:47  pwagner
! Accepts hdfVersion optional arg; restored hdf5 module use
!
! Revision 1.4  2002/01/23 22:41:21  pwagner
! Handles optional hdfVersion parameter
!
! Revision 1.3  2002/01/23 00:54:10  pwagner
! Added mls_hdf_version function
!
! Revision 1.2  2002/01/18 23:55:49  pwagner
! Various strategems to try writing metadata; all fail
!
! Revision 1.1  2002/01/18 00:51:22  pwagner
! First commit
!
