! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
module MLSFiles               ! Utility file routines
  !===============================================================================
  use HDFEOS, only: gdclose, gdopen, swclose, swopen
  use machine, only: io_error
  use MLSCommon, only: i4, NameLen
  use MLSStrings, only: Capitalize, LowerCase, Reverse
  use output_m, only: output
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS, Pgs_smf_getMsg, &
    & PGSd_IO_Gen_RSeqFrm, PGSd_IO_Gen_RSeqUnf, & 
    & PGSd_IO_Gen_RDirFrm, PGSd_IO_Gen_RDirUnf, & 
    & PGSd_IO_Gen_WSeqFrm, PGSd_IO_Gen_WSeqUnf, & 
    & PGSd_IO_Gen_WDirFrm, PGSd_IO_Gen_WDirUnf, & 
    & PGSd_IO_Gen_USeqFrm, PGSd_IO_Gen_USeqUnf, & 
    & PGSd_IO_Gen_UDirFrm, PGSd_IO_Gen_UDirUnf, & 
    & PGSd_IO_Gen_ASeqFrm, PGSd_IO_Gen_ASeqUnf, &
    & PGS_IO_GEN_CloseF, PGS_IO_GEN_OpenF, PGSd_PC_FILE_PATH_MAX
  implicit none
  public

  private :: ID

  !------------------- RCS Ident Info -----------------------
  character(LEN=130) :: Id = &
    "$Id$"
  !----------------------------------------------------------

  ! This isn't NameLen because it may have a path prefixed
  integer, parameter :: MAXFILENAMELENGTH=PGSd_PC_FILE_PATH_MAX

  ! These are error codes that may be returned by GetPCFromRef

  integer, parameter :: NAMENOTFOUND=-1
  integer, parameter :: INVALIDPCRANGE=NAMENOTFOUND-1

  ! These are error codes that may be returned by mls_io_gen_openF

  integer, parameter :: UNKNOWNFILEACCESSTYPE=-999
  integer, parameter :: UNKNOWNTOOLBOXMODE=UNKNOWNFILEACCESSTYPE+1
  integer, parameter :: NOFREEUNITS=UNKNOWNTOOLBOXMODE+1
  integer, parameter :: MUSTSUPPLYFILENAMEORPC=NOFREEUNITS+1

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
    integer(i4)                              :: thePC
    integer(i4),  intent(OUT)                :: ErrType
    logical,  intent(IN)                     :: caseSensitive
    integer(i4),  optional                   :: versionNum
    logical,  optional, intent(IN)           :: debugOption
    character (LEN=*),  optional, intent(IN)       :: path
    character (LEN=*), optional, intent(out) :: ExactName

    ! Local variables

    character (LEN=MAXFILENAMELENGTH) :: MatchName, TryName
    character (LEN=MAXFILENAMELENGTH) :: PhysicalName, MatchPath
    integer                       ::     version, returnStatus
   logical ::                            debug

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

      !		if(debug) then
      !			call output('pc: ' )
      !			call output(thePC, advance='no')
      !			call output('              version: ' )
      !			call output(version, advance='yes')
      !		endif

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
        
        call split_path_name(TryName, MatchPath, PhysicalName)

        if ( index(PhysicalName, trim(MatchName)) /= 0 )then
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
    & FileName, PCBottom, PCTop, versionNum, thePC) &
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

    ! Local variables

    logical, parameter :: PRINT_EVERY_OPEN=.true.
    integer, parameter :: FH_ON_ERROR=-99
    integer, parameter :: KEYWORDLEN=12			! Max length of keywords in OPEN(...)
    character (LEN=MAXFILENAMELENGTH) :: myName
    integer(i4) :: myPC
    integer                       :: version, returnStatus
    logical       :: tiedup
    character (LEN=KEYWORDLEN) :: access, action, form, position, status
    integer                       :: unit

    ! begin

    ! In case of premature return
    theFileHandle = FH_ON_ERROR

    if(present(versionNum)) then
      version = versionNum
    else
      version = 1
    endif

    if(present(thePC)) then
      myPC = thePC
      returnStatus = Pgs_pc_getReference(thePC, version, &
        & myName)
    elseif(present(FileName)) then
      myName = FileName
      if(LowerCase(toolbox_mode(1:2)) == 'pg') then
        myPC = GetPCFromRef(FileName, PCBottom, PCTop, &
          &	 caseSensitive, returnStatus, versionNum)
      endif
    else
      ErrType = MUSTSUPPLYFILENAMEORPC
      return
    endif

    select case (LowerCase(toolbox_mode(1:2)))

    case('pg')
      if(returnStatus == 0) then
        ErrType = PGS_IO_Gen_OpenF(myPC, FileAccessType, record_length, &
          & theFileHandle, version)
      endif

    case('sw')
      if(returnStatus == 0) then
        theFileHandle = swopen(myName, FileAccessType)
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
        theFileHandle = gdopen(myName, FileAccessType)
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
        return
      endif


      if(access /= 'direct') then
        open(unit=unit, access=access, action=action, form=form, &
          & position=position, status=status, file=myName, iostat=ErrType)
      else
        open(unit=unit, access=access, action=action, form=form, &
          & status=status, file=myName, iostat=ErrType)
      endif

      if(ErrType /= 0 .or. PRINT_EVERY_OPEN) then
        call output( 'Fortran opening unit ', advance='no')
        call output(  unit, advance='yes')
        call output( 'access ' // access, advance='yes')
        call output( 'action ' // action, advance='yes')
        call output( 'form ' // form, advance='yes')
        call output( 'position ' // position, advance='yes')
        call output( 'status ' // status, advance='yes')
        call output( 'file ' // myName, advance='yes')
      endif

      if(ErrType /= 0) then
        call output( 'iostat ', advance='no')
        call output(  ErrType, advance='yes')
        call io_error('io error in MLSFiles: mls_io_gen_openF' // &
          & ' Fortran open', ErrType, myName)
      else
        theFileHandle = unit
      endif

    case default
      ErrType = UNKNOWNTOOLBOXMODE

    end select

    if(ErrType /= 0) then
      theFileHandle = FH_ON_ERROR
    endif

  end function mls_io_gen_openF

  ! ---------------------------------------------  mls_io_gen_closeF  -----

  ! This function closes a generic file using either the toolbox
  ! or else a Fortran OPEN statement
  ! according to toolbox_mode:
  ! 'sw' for swath files opened with swopen
  ! 'gd' for grid files opened with gdopen
  ! 'pg' for generic files opened with pgs_io_gen_openF
  ! 'cl' for l3ascii files opened with simple fortran 'open'

  ! It returns a non-zero error status only if unsuccessful

  ! If must be given a FileHandle as an arg
  ! (A later version may allow choice between file handle and file name)
  function mls_io_gen_closeF(toolbox_mode, theFileHandle) &
    &  result (ErrType)

    ! Dummy arguments
    integer(i4)  :: ErrType
    integer(i4), intent(IN)  :: theFileHandle
    character (LEN=*), intent(IN)   :: toolbox_mode

    ! Local variables

    logical, parameter :: PRINT_EVERY_CLOSE=.true.

    select case (LowerCase(toolbox_mode(1:2)))

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
   logical, parameter :: DEBUG = .false.

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
    
    if(DEBUG) then
      print*, trim(mirrored_ffn), ' ', loc
   endif

    if(loc <= 0) then
      path = ' '
      name = adjustl(full_file_name)
    elseif(loc == 1) then
      path = adjustl(full_file_name)
      name = ' '
    else
    if(DEBUG) then
      print*, trim(mirrored_ffn(loc:))
      print*, trim(mirrored_ffn(:loc-1))
   endif
      path = adjustl(Reverse(mirrored_ffn(loc:)))
      name = adjustl(Reverse(mirrored_ffn(:loc-1)))
    endif

  end subroutine split_path_name

  !====================
end module MLSFiles
!====================

!
! $Log$
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
!
