! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
module MLSFiles               ! Utility file routines
  !===============================================================================
  use Hdf, only: DFACC_CREATE, DFACC_RDONLY, DFACC_READ, DFACC_RDWR, &
    & sfstart, sfend
  use HDF5, only: hid_t, h5fopen_f, h5fcreate_f, h5fclose_f, h5fis_hdf5_f, &
      H5F_ACC_RDONLY_F, H5F_ACC_RDWR_F, H5F_ACC_TRUNC_F, H5F_ACC_EXCL_F
  use HDFEOS, only: gdclose, gdopen, swclose, swopen, swinqswath
  use HDFEOS5, only: he5_swclose, he5_swopen, he5_swinqswath, &
    & he5_gdopen, he5_gdclose, &
    & HE5F_ACC_TRUNC, HE5F_ACC_RDONLY, HE5F_ACC_RDWR
  use machine, only: io_error
  use MLSCommon, only: i4, BareFNLen, FileNameLen, MLSFile_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_Warning
  use MLSStrings, only: Capitalize, LowerCase, Reverse, SortArray
  use output_m, only: blanks, output
  use SDPToolkit, only: &
    & HDF5_ACC_CREATE, HDF5_ACC_RDONLY, HDF5_ACC_RDWR,  &
    & Pgs_pc_getReference, PGS_S_SUCCESS, &
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

  public :: AddFileToDataBase, GetPCFromRef, get_free_lun, mls_io_gen_openF, &
  & mls_io_gen_closeF, split_path_name, RmFileFromDataBase, &
  & mls_hdf_version, mls_inqswath, mls_sfstart, mls_sfend, &
  & mls_openFile, mls_closeFile, MLSFile_T, Deallocate_filedatabase, &
  & open_MLSFile, close_MLSFile, Dump

  !------------------- RCS Ident Info -----------------------
  character(LEN=130) :: Id = &
    "$Id$"
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !----------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (parameters)
! HDFVERSION_4       integer corresponding to hdf4
! HDFVERSION_5       integer corresponding to hdf5
! WILDCARDHDFVERSION integer corresponding to autorecognize hdf4/hdf5
! NAMENOTFOUND       GetPCFromRef unable to find named file
! INVALIDPCRANGE     GetPCFromRef given invalid pc range
! CANTALLOCATENAMEARRAY
!                    GetPCFromRef unable to allocate memory needed
! UNKNOWNFILEACCESSTYPE
!                    mls_io_gen_openF given invalid file access type
! UNKNOWNTOOLBOXMODE mls_io_gen_openF given invalid toolbox mode
! NOFREEUNITS        mls_io_gen_openF ran out of free units
! MUSTSUPPLYFILENAMEORPC
!                    mls_io_gen_openF given neither required argument
! NOPCIFNOTOOLKIT    mls_io_gen_openF must be given filename if toolkitless
! NOSUCHHDFVERSION   mls_io_gen_openF given invalid hdf version
! MUSTSUPPLYFILENAME Required filename arg missing
! ERRORINH5FFUNCTION Internal hdf5 library function error
! WRONGHDFVERSION    Supplied hdf version differs from actual file version

!     (subroutines and functions)
! AddFileToDataBase  Enters a FileName, id, etc. into the database
! close_MLSFile      Closes an mls file (of any type)
! Deallocate_filedatabase
!                    Deallocates file database, closing any still-open files
! Dump               Dumps file info: type, access, name, etc.
! GetPCFromRef       Turns a FileName into the corresponding PC
! get_free_lun       Gets a free logical unit number
! mls_closeFile      Closes a file opened by mls_openFile
! mls_hdf_version    Returns one of 'hdf4', 'hdf5', or '????'
! mls_inqswath       A wrapper for doing swingswath for versions 4 and 5
! mls_io_gen_openF   Opens a generic file using the toolbox or Fortran OPEN
! mls_io_gen_closeF  Closes a generic file using the toolbox or Fortran CLOSE
! mls_openFile       Opens an hdf5 file
! mls_sfend          Closes a file opened by mls_sfstart
! mls_sfstart        Opens an hdf file for writing metadata
! open_MLSFile       Opens an mls file (of any type)
! RmFileFromDataBase Removes a FileName, id, etc. from the database
! split_path_name    splits the input path/name into path and name
! === (end of toc) ===

! (The following 2 are currently private, but could be made public if needed)
! hdf2hdf5_fileaccess
!                    Translates version 4 hdf access codes to version 5
! he2he5_fileaccess
!                    Translates version 4 hdfeos access codes to version 5

! === (start of api) ===
! i4 GetPCFromRef (char* FileName, i4 PCBottom, i4 PCTop,
!     log caseSensitive, i4 ErrType, [i4 versionNum], [log debugOption],
!     [char* path], [char* ExactName])
! i4 mls_io_gen_openF (char* toolbox_mode, log caseSensitive, i4 ErrType, 
!     i4 record_length, i4 FileAccessType, 
!     [char* FileName], [i4 PCBottom], [i4 PCTop], [i4 versionNum],
!     [log unknown], [i4 thePC], [int hdfVersion], [log debugOption])
! i4 mls_io_gen_closeF (char* toolbox_mode, i4 theFileHandle, 
!     [char* FileName], [int hdfVersion])
! split_path_name (char* full_file_name, char* path, char* name, [char slash])
! int mls_inqswath (char* FileName, char* swathList, int strBufSize,
!     [int hdfVersion])
! int mls_sfstart (char* FileName, i4 FileAccess, [int hdfVersion])
! int mls_sfend (int sdid, [int hdfVersion])
! i4 mls_hdf_version (char* FileName, [i4 preferred_version], [i4 AccessType])
! mls_openFile (char* filename, char* access, hid_t file_id)
! mls_closeFile (hid_t file_id)
! === (end of api) ===

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
  integer, parameter, public :: NOSUCHHDFVERSION=NOPCIFNOTOOLKIT+1

  ! These are error codes that may be returned by mls_io_gen_openF
  ! or by mls_hdf_version if you call it directly

  integer, parameter, public :: MUSTSUPPLYFILENAME=NOSUCHHDFVERSION+1
  integer, parameter, public :: ERRORINH5FFUNCTION=MUSTSUPPLYFILENAME+1
  integer, parameter, public :: WRONGHDFVERSION=ERRORINH5FFUNCTION+1
  
  ! Whether to use PGS_MET commands in mls_sf(start)(end)
  ! (The alternative is to use hdf5 calls directly)
  logical, parameter :: PGS_MET4MLS_SF = .true.
  ! Whether to pass hdf5_acc types to PGS_MET 
  ! (The alternative is to pass h5f_acc directly)
  ! Contradicting what the documents say,
  ! currently (TK 5.2.7.4) the toolkit routine
  ! PGS_MET_HDFSDStart.c assumes the File Access is of
  ! one of the following types:
  ! (1) hdf4_acc: DFACC_RDWR, DFACC_RDONLY, or DFACC_CREATE
  !      (equiv. to HDF4_ACC_RDWR, etc.)
  ! (2) hdf5_acc: HDF5_ACC_RDWR, HDF5_ACC_RDONLY, or HDF5_ACC_CREATE
  !          * * *    w a r n i n g   * * *
  ! According to mls/hirdls telecom of 03-12-2002 this will change
  ! in a future release of the toolkit
  ! That means some programmer (paw?) will have to revisit
  ! this and see what the toolkit code actually does,
  ! not just what the docs say
  logical, parameter :: HDF5_ACC_TYPES_TO_MET = .true.
  integer, parameter :: HDF5_ACC_DEFAULT = HE5F_ACC_RDWR ! HDF5_ACC_RDWR
  
  ! The only legal unit numbers that files may be assigned
  ! for use by Fortran opens, closes, reads and writes
  integer, parameter :: bottom_unit_num=1
  integer, parameter :: top_unit_num=99

  interface DUMP
    module procedure Dump_MLSFile
    module procedure Dump_FileDatabase
  end interface

! The following data type was moved to MLSCommon in an attempt to code around
! Lahey's compiler that can cause compile times to exceed mission lifetime
! >   ! Information describing the files used by the mls software
! >   ! Stop passing file handles back & forth bewteen routines
! >   ! -- pass one of these instead
! >   TYPE MLSFile_T
! >     CHARACTER (LEN=8) :: Type=""  ! e.g., 'ascii', 'hdf', 'swath', 'binary'
! >     CHARACTER (LEN=8) :: access=""  ! e.g., 'rdonly', 'write', 'rdwrite'
! >     CHARACTER (LEN=8) :: content=""  ! e.g., 'l1brad', 'l2gp', 'l2aux'
! >     CHARACTER (LEN=FileNameLen) :: Name=""  ! its name (usu. w/path)
! >     INTEGER :: File_Id=0     ! The HDF ID (handle) or io unit for the file
! >     INTEGER :: PCF_Id=0      ! The PCF ID (ref), if any,  for the file
! >     INTEGER :: HDFVersion=0  ! Which hdf version is the file if hdf(eos)
! >     LOGICAL :: StillOpen=.false.
! >   END TYPE MLSFile_T
contains

  !-----------------------------------------  AddFileToDatabase  -----
  integer function AddFileToDatabase ( DATABASE, ITEM )

  ! This routine adds a vector to a database of such vectors, 
  ! creating the database if necessary.

    ! Dummy arguments
    type (MLSFile_T), dimension(:), pointer :: DATABASE
    type (MLSFile_T), intent(in) ::            ITEM

    ! Local variables
    type (MLSFile_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddFileToDatabase = newSize
  end function AddFileToDatabase

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
  ! 'gd' for grid files opened with gdopen
  ! 'hg' for hdf4/5 files opened with sfstart/h5fopen_f; e.g. l1brad
  ! 'pg' for generic files opened with pgs_io_gen_openF
  ! 'op' for l3ascii files opened with simple fortran 'open'
  ! 'sw' for swath files opened with swopen

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
    & hdfVersion, debugOption, inp_rec_length) &
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
    integer, optional, intent(in) :: inp_rec_length

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
          & he2he5_fileaccess(FileAccessType))
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
          call output(he2he5_fileaccess(FileAccessType), advance='yes')  
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
          & he2he5_fileaccess(FileAccessType))
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

    case('hg')
      theFileHandle = -1   ! This is the error value expected by hdf4/5
      if(returnStatus /= 0) then
        ErrType = returnStatus
        return
      elseif(myhdfVersion == HDFVERSION_5) then
      
        call h5fopen_f(trim(myName), &
          & hdf2hdf5_fileaccess(FileAccessType), theFileHandle, ErrType)
      elseif(myhdfVersion == HDFVERSION_4) then
        theFileHandle = sfstart(trim(myName), FileAccessType)
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

      if ( present(inp_rec_length) ) then
        open(unit=unit, recl=inp_rec_length, form=form, &
          & status=status, file=trim(myName), iostat=ErrType)
      elseif(access /= 'direct') then
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

  ! ---------------------------------------------  open_MLSFile  -----

  ! This routine opens an mls file using either the toolbox
  ! or else a Fortran OPEN statement
  ! as appropriate, filling the fields of the MLSFile as needed
  
  ! It must be supplied an MLSFile_t for its arg

  ! Naturally, at least one of the arg's fields must be non-default values:
  ! MLSFile%{Name, PC_Id}

  ! It generates an error if unsuccessful

  subroutine open_MLSFile(MLSFile, PCBottom, PCTop, inp_rec_length, unknown)

    ! Dummy arguments
    type (MLSFile_T)          :: MLSFile
    integer(i4),  optional, intent(IN)   :: PCBottom, PCTop
    logical, optional, intent(in) :: unknown
    integer, optional, intent(in) :: inp_rec_length

    ! Local
    integer(i4)  :: ErrType
    logical,  parameter      :: caseSensitive = .true.
    integer(i4)       :: FileAccessType

    logical, parameter :: DEFAULT_PRINT_EVERY_OPEN=.false.
    integer, parameter :: FH_ON_ERROR=-99
    integer, parameter :: DEFAULTRECLEN=0
    integer, parameter :: KEYWORDLEN=12			! Max length of keywords in OPEN(...)
    integer(i4) :: myPC
    integer                       :: version, returnStatus
    logical       :: tiedup
    character (LEN=KEYWORDLEN) :: access, action, form, position, status
    integer                       :: record_length
    integer                       :: unit

    ! begin
    version = 1
    if ( present(PCBottom) .and. present(PCTop) ) then
      if ( MLSFile%PCF_Id > 0 ) then
        returnStatus = Pgs_pc_getReference(MLSFile%PCF_Id, version, &
          & MLSFile%Name)
      elseif ( MLSFile%Name /= '' ) then
          MLSFile%PCF_Id = GetPCFromRef(trim(MLSFile%Name), PCBottom, PCTop, &
            &  caseSensitive, returnStatus, version, &
            & debugOption=.false.)

      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
         & 'You must supply either a name or a PCFid to open a file' )
      endif
    endif
    select case (MLSFile%access)
    case('rdonly')
      FileAccessType = DFACC_RDONLY
    case('write')
      FileAccessType = DFACC_CREATE
    case('rdwrite')
      FileAccessType = DFACC_RDWR
    case default
      FileAccessType = DFACC_RDWR

    end select

    select case (MLSFile%Type)

    case('asc-tk', 'bin-tk')
     ! Using Toolkit
      if ( MLSFile%PCF_Id > 0 ) then
        returnStatus = Pgs_pc_getReference(MLSFile%PCF_Id, version, &
          & MLSFile%Name)
      elseif ( MLSFile%Name /= '' ) then
          MLSFile%PCF_Id = GetPCFromRef(trim(MLSFile%Name), PCBottom, PCTop, &
            &  caseSensitive, returnStatus, version, &
            & debugOption=.false.)

      else
        CALL MLSMessage ( MLSMSG_Error, moduleName,  &
         & 'You must supply either a name or a PCFid to open a file' )
      endif
      ErrType = PGS_IO_Gen_OpenF(MLSFile%PCF_Id, FileAccessType, record_length, &
          & MLSFile%File_id, version)

    case('swath')
      if(MLSFile%HDFVersion == HDFVERSION_5) then
        MLSFile%File_Id = he5_swopen(trim(MLSFile%Name), &
          & he2he5_fileaccess(FileAccessType))
      elseif(MLSFile%HDFVersion == HDFVERSION_4) then
         MLSFile%File_Id = swopen(trim(MLSFile%Name), FileAccessType)
      else
        ErrType = NOSUCHHDFVERSION
      endif

    case('grid')
      if(MLSFile%HDFVersion == HDFVERSION_5) then
        MLSFile%File_Id = he5_gdopen(trim(MLSFile%Name), &
          & he2he5_fileaccess(FileAccessType))
      elseif(MLSFile%HDFVersion == HDFVERSION_4) then
        MLSFile%File_Id = gdopen(trim(MLSFile%Name), FileAccessType)
      else
        ErrType = NOSUCHHDFVERSION
      endif

    case('hdf')
      if(MLSFile%HDFVersion == HDFVERSION_5) then
        call h5fopen_f(trim(MLSFile%Name), &
          & hdf2hdf5_fileaccess(FileAccessType), MLSFile%File_Id, ErrType)
      elseif(MLSFile%HDFVersion == HDFVERSION_4) then
        MLSFile%File_Id = sfstart(trim(MLSFile%Name), FileAccessType)
      else
        ErrType = NOSUCHHDFVERSION
      endif

    case('ascii', 'binary')
      if(FileAccessType == DFACC_RDONLY) then
        status = 'old'
        access = 'sequential'
        form = 'formatted'
        action = 'read'
        position = 'rewind'
      elseif(FileAccessType == DFACC_CREATE) then
        status = 'new'
        access = 'sequential'
        form = 'formatted'
        action = 'write'
        position = 'rewind'
      elseif(FileAccessType == DFACC_RDWR) then
        status = 'old'
        access = 'sequential'
        form = 'formatted'
        action = 'readwrite'
        position = 'rewind'
      else
        ErrType = UNKNOWNFILEACCESSTYPE
      endif
      
      if ( MLSFile%Type == 'binary' ) then
        access = 'direct'
        form = 'unformatted'
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
      endif

      if (present(unknown)) then
        if (unknown) status = 'unknown'
      end if

      if ( present(inp_rec_length) ) then
        open(unit=unit, recl=inp_rec_length, form=form, &
          & status=status, file=trim(MLSFile%Name), iostat=ErrType)
      elseif(access /= 'direct') then
        open(unit=unit, access=access, action=action, form=form, &
          & position=position, status=status, file=trim(MLSFile%Name), iostat=ErrType)
      else
        open(unit=unit, access=access, action=action, form=form, &
          & status=status, file=trim(MLSFile%Name), iostat=ErrType)
      endif

      if(ErrType /= 0) then
        call output( 'Fortran opening unit ', advance='no')
        call output(  unit, advance='yes')
        call output( 'access ' // access, advance='yes')
        call output( 'action ' // action, advance='yes')
        call output( 'form ' // form, advance='yes')
        call output( 'position ' // position, advance='yes')
        call output( 'status ' // status, advance='yes')
        call output( 'file ' // trim(MLSFile%Name), advance='yes')
        call output( 'iostat ', advance='no')
        call output(  ErrType, advance='yes')
        call io_error('io error in MLSFiles: mls_io_gen_openF' // &
          & ' Fortran open', ErrType, trim(MLSFile%Name))
      else
        MLSFile%File_Id = unit
      endif

    case default
      ErrType = UNKNOWNTOOLBOXMODE

    end select

    if(ErrType /= 0) then
      call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot open file " // trim(MLSFile%Name) )
    endif
    
  end subroutine open_MLSFile

  ! ---------------------------------------------  mls_io_gen_closeF  -----

  ! This function closes a generic file using either the toolbox
  ! or else a Fortran CLOSE statement
  ! according to toolbox_mode:
  ! 'cl' for l3ascii files opened with simple fortran 'open'
  ! 'gd' for grid files opened with gdopen
  ! 'hg' for hdf4/5 files opened with sfstart/h5fopen_f; e.g. l1brad
  ! 'pg' for generic files opened with pgs_io_gen_openF
  ! 'sw' for swath files opened with swopen

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

    case('hg')
      if(myhdfVersion == HDFVERSION_5) then
        call h5fclose_f(theFileHandle, ErrType)
      elseif(myhdfVersion == HDFVERSION_4) then
        ErrType = sfend(theFileHandle)
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

  ! ---------------------------------------------  close_MLSFile  -----

  ! This subroutine closes an mls file using either the toolbox
  ! or else a Fortran CLOSE statement
  ! as appropriate

  ! It generates an error if unsuccessful

  ! It must be supplied an MLSFile_t for its arg
  subroutine close_MLSFile(MLSFile)

    ! Dummy arguments
    type (MLSFile_T) ::            MLSFile

    ! Local
    integer :: ErrType

    ! begin

    select case (MLSFile%Type)

    case('asc-tk', 'bin-tk')
      ErrType = PGS_IO_Gen_CLoseF(MLSFile%File_Id)

    case('swath')
      if(MLSFile%HDFVersion == HDFVERSION_5) then
        ErrType = he5_swclose(MLSFile%File_Id)
      elseif(MLSFile%HDFVersion == HDFVERSION_4) then
        ErrType = swclose(MLSFile%File_Id)
      else
        ErrType = NOSUCHHDFVERSION
      endif

    case('grid')
      if(MLSFile%HDFVersion == HDFVERSION_5) then
        ErrType = he5_gdclose(MLSFile%File_Id)
      elseif(MLSFile%HDFVersion == HDFVERSION_4) then
        ErrType = gdclose(MLSFile%File_Id)
      else
        ErrType = NOSUCHHDFVERSION
      endif

    case('hdf')
      if(MLSFile%HDFVersion == HDFVERSION_5) then
        call h5fclose_f(MLSFile%File_Id, ErrType)
      elseif(MLSFile%HDFVersion == HDFVERSION_4) then
        ErrType = sfend(MLSFile%File_Id)
      else
        ErrType = NOSUCHHDFVERSION
      endif

    case('ascii', 'binary')
      close(unit=MLSFile%File_Id, iostat=ErrType)		

      if(ErrType /= 0) then
        call output( 'Fortran closing unit ', advance='no')
        call output(  MLSFile%File_Id, advance='yes')
        call output( 'iostat ', advance='no')
        call output(  ErrType, advance='yes')
        call io_error('io error in MLSFiles: mls_io_gen_closeF' // &
          & ' Fortran close', ErrType, 'unknown')
      endif

    case default
      ErrType = UNKNOWNTOOLBOXMODE

    end select
    if ( ErrType /= 0 ) &
      &  CALL MLSMessage ( MLSMSG_Error, moduleName,  &
      & "Unable to close file " // trim(MLSFile%Name) )

  end subroutine close_MLSFile

  !-----------------------------------------  RmFileFromDatabase  -----
  integer function RmFileFromDatabase ( DATABASE, ITEM )

  ! This routine removes a vector from a database of such vectors, 
  ! deallocating the database if necessary.
  ! Alas, doesn't work--we need to know how to undecorate character tree
  ! first before we will be able to make it work; sorry (P. Wagner)

    ! Dummy arguments
    type (MLSFile_T), dimension(:), pointer :: DATABASE
    type (MLSFile_T), intent(in) ::            ITEM

    ! Local variables
    type (MLSFile_T), dimension(:), pointer :: tempDatabase
    logical, parameter                     :: okToDeallocEmptyDB = .FALSE.
    include "rmItemFromDatabase.f9h"
    call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot yet (ever?) rm File from database" ) 

    RmFileFromDatabase = newSize
  end function RmFileFromDatabase

  ! ----------------------  Deallocate_filedatabase  -----

  ! This routine deallocates the file database, closing any files
  ! that may still be open

  subroutine Deallocate_filedatabase(database)

    ! Arguments
    type (MLSFile_T), dimension(:), pointer :: DATABASE
    ! Local variables
    integer :: i, error
    ! Executable
    if ( .not. associated(database) ) return
    do i=1, size(database)
      if ( database(i)%StillOpen ) call close_MLSFile(database(i))
    enddo
    Deallocate ( database, stat=error )
    if ( error /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, &
        & ModuleName, "Cannot deallocate file database" )
  end subroutine Deallocate_filedatabase

  ! ------------------------------------------ Dump_FileDataBase ------------

  subroutine Dump_FileDataBase ( database, Name )

    ! Dummy arguments
    type (MLSFile_T), intent(in) ::          database(:)
    character(len=*), intent(in), optional :: Name

    ! Local variables
    integer :: i, dim
    
    call output ( '============ MLS File Data Base ============', advance='yes' )
    call output ( ' ', advance='yes' )
    if ( present(name) ) then
      call output ( 'MLS File Database name: ', advance='no' )
      call output ( name, advance='yes' )
    endif
    if ( size(database) < 1 ) then
      call output ( '**** MLS File Database empty ****', advance='yes' )
      return
    endif
    do i = 1, size(database)
      call dump(database(i))
    end do
      
  end subroutine Dump_FileDataBase

  ! ------------------------------------------ Dump_MLSFile ------------

  subroutine Dump_MLSFile ( MLSFile )

    ! Dummy arguments
    type (MLSFile_T), intent(in) ::          MLSFile

    ! Executable code
    call output ( 'MLS File Info: ')                                  
    call output ( '(name) ')                                  
    call output ( trim(MLSFile%Name), advance='yes')                                  
    call output ( '    Type         : ')                                
    call output ( trim(MLSFile%Type), advance='yes')                                  
    call output ( '    Access       : ')                                
    call output ( trim(MLSFile%access), advance='yes')                                  
    call output ( '    content      : ')                                
    call output ( trim(MLSFile%content), advance='yes')                                  
    call output ( '    File ID      : ')                                
    call output ( MLSFile%File_Id, advance='yes')                                  
    call output ( '    PCF ID       : ')                                
    call output ( MLSFile%PCF_Id, advance='yes')                                  
    call output ( '    hdf version  : ')                                
    call output ( MLSFile%HDFVersion, advance='yes')                                  
    call output ( '    Open?        : ')                                
    call output ( MLSFile%StillOpen, advance='yes')                                  
  end subroutine Dump_MLSFile

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

  ! ---------------------------------------------  he2he5_fileaccess  -----

  ! This function converts hdfeos2 file access types to
  ! corresponding hdfeos5 numbers

  function he2he5_fileaccess(FileAccesshdf4) result (FileAccesshdf5)

    ! Arguments

    integer(i4), intent(IN)       :: FileAccesshdf4
    integer(i4)                   :: FileAccesshdf5
    
!    integer(i4), parameter        :: H5F_ACC_RDONLY = 0
!    integer(i4), parameter        :: H5F_ACC_RDWR   = 1
!    integer(i4), parameter        :: H5F_ACC_TRUNC  = 2

    ! begin
    select case (FileAccesshdf4)

    case(DFACC_CREATE)
      FileAccesshdf5 = HE5F_ACC_TRUNC   ! H5F_ACC_TRUNC

    case(DFACC_READ)                      ! also , DFACC_RDONLY
      FileAccesshdf5 = HE5F_ACC_RDONLY   ! H5F_ACC_RDONLY

    case default
      FileAccesshdf5 = HE5F_ACC_RDWR   ! H5F_ACC_RDWR

    end select

  end function he2he5_fileaccess

  ! ---------------------------------------------  hdf2hdf5_fileaccess  -----

  ! This function converts hdf4 file access types to
  ! corresponding hdf5 numbers

  function hdf2hdf5_fileaccess(FileAccesshdf4) result (FileAccesshdf5)

    ! Arguments

    integer(i4), intent(IN)       :: FileAccesshdf4
    integer(kind(H5F_ACC_EXCL_F)) :: FileAccesshdf5
    
    ! begin
    select case (FileAccesshdf4)

    case(DFACC_CREATE)
      FileAccesshdf5 = H5F_ACC_EXCL_F  ! H5F_ACC_TRUNC_F

    case(DFACC_READ)    ! also , DFACC_RDONLY
      FileAccesshdf5 = H5F_ACC_RDONLY_F

    case default
      FileAccesshdf5 = H5F_ACC_RDWR_F

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
  ! Right now, it works for hdf4 files in general, but only for adding
  ! metadata to hdf5 files
  
  ! Therefore, when the grand unified hdf4/hdf5 interfaces are
  ! implemented this will probably need to take an added arg:
  ! the logical addingMetadata

  function mls_sfstart(FileName, FileAccess, hdfVersion, addingmetadata)

    ! Arguments

    character (len=*), intent(in) :: FILENAME
    integer(i4), intent(IN)       :: FileAccess ! (one of the hdf4 types)
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: addingmetadata
    ! Local variables
    integer                       :: mls_sfstart
    integer                       :: returnStatus
    integer                       :: access_prp_default
    integer                       :: myhdfVersion
    integer                       :: myAccess
    logical                       :: myaddingmetadata
    
    integer, parameter :: h5p_default_f = 0
    integer, external :: PGS_MET_SFstart
    logical, parameter :: DEBUG = .false.

    ! begin
   returnStatus = 0
   myaddingmetadata = .false.
   if ( present (addingmetadata) ) myaddingmetadata = addingmetadata
   if ( DEBUG ) then
     call output ('Entering mls_sfstart with args ', advance='no')
     call output ('File name: ', advance='no')
     call output (trim(Filename), advance='no')
     call blanks (2)
     call output ('FileAccess: ', advance='no')
     call output (FileAccess, advance='no')
     if ( present(hdfVersion) ) then
       call blanks (2)
       call output ('hdfVersion: ', advance='no')
       call output (hdfVersion, advance='no')
     endif
     call blanks (2)
     call output ('adding meta data?: ', advance='no')
     call output (myaddingmetadata, advance='yes')
   endif
   if ( present(hdfVersion) ) then
     myhdfVersion = hdfVersion
   else
     myhdfVersion = DEFAULT_HDFVERSION
   endif
   if ( myhdfVersion == HDFVERSION_4) then
     mls_sfstart = sfstart (FileName, FileAccess)
     return
   elseif ( myhdfVersion /= HDFVERSION_5) then
     mls_sfstart = WRONGHDFVERSION
     return
   endif
!   if ( PGS_MET4MLS_SF ) then
   if ( myaddingmetadata ) then
     if ( .not. HDF5_ACC_TYPES_TO_MET ) then
       myAccess = he2he5_fileaccess(FileAccess)
     elseif ( FileAccess == DFACC_RDWR ) then
       myAccess = HE5F_ACC_RDWR  ! HDF5_ACC_RDWR
     elseif ( FileAccess == DFACC_CREATE ) then
       myAccess = HE5F_ACC_TRUNC ! HDF5_ACC_CREATE
     elseif ( FileAccess == DFACC_RDONLY ) then
       myAccess = HE5F_ACC_RDONLY ! HDF5_ACC_RDONLY
     else
       myAccess = HDF5_ACC_DEFAULT
     endif
     returnStatus = PGS_MET_SFstart(trim(FileName), myAccess, mls_sfstart)
   else
     access_prp_default = h5p_default_f    ! Can't figure out what this means
! >      print *, 'About to call h5fopen_f'
! >      print *, 'FileAccess: ', FileAccess
! >      print *, 'DFACC_CREATE: ', DFACC_CREATE
! >      print *, 'DFACC_RDWR: ', DFACC_RDWR
! >      print *, 'DFACC_RDONLY: ', DFACC_RDONLY
! >      print *, 'hdf2hdf5_fileaccess(FileAccess): ', hdf2hdf5_fileaccess(FileAccess)
! >      print *, 'H5F_ACC_RDWR_F: ', H5F_ACC_RDWR_F
! >      print *, 'H5F_ACC_EXCL_F: ', H5F_ACC_EXCL_F
! >      print *, 'H5F_ACC_RDONLY_F: ', H5F_ACC_RDONLY_F
! >      print *, 'H5F_ACC_TRUNC_F: ', H5F_ACC_TRUNC_F
     ! call h5fopen_f(trim(FileName), hdf2hdf5_fileaccess(FileAccess), &
     ! & mls_sfstart, returnStatus)
!      & mls_sfstart, returnStatus, access_prp_default)  ! so abandoning it
     select case (FileAccess)
     case (DFACC_CREATE)
       call h5fcreate_f(trim(filename), H5F_ACC_EXCL_F, mls_sfstart, returnStatus)
       ! call mls_openFile(filename, 'create', mls_sfstart, HDFVERSION_5)
     case (DFACC_RDWR)
       call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, mls_sfstart, returnStatus)
       ! call mls_openFile(filename, 'update', mls_sfstart, HDFVERSION_5)
     case (DFACC_RDONLY)
       call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, mls_sfstart, returnStatus)
       ! call mls_openFile(filename, 'readonly', mls_sfstart, HDFVERSION_5)
     case default
       call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, mls_sfstart, returnStatus)
       ! call mls_openFile(filename, 'readonly', mls_sfstart, HDFVERSION_5)
     end select
   endif
   if ( returnStatus /= 0 .and. myAddingMetaData) then                                            
     call output ('Try again--PGS_MET_SFstart still unhappy; returns ')    
     call output (returnStatus, advance='yes')                             
     mls_sfstart = -1                                                      
   elseif ( returnStatus /= 0) then                                            
     call output ('Try again--h5fopen_f still unhappy; returns ')    
     call output (returnStatus, advance='yes')                             
     mls_sfstart = -1                                                      
   endif                                                                  
   if ( DEBUG ) then
     call output ('Returning from mls_sfstart an sdid: ', advance='no')
     call output (mls_sfstart, advance='yes')
   endif

  end function mls_sfstart

  ! ---------------------------------------------  mls_sfend  -----

  ! This function acts as a wrapper to allow hdf5 or hdf4 routines to be called
  ! Right now, it works for hdf4 files in general, but only for adding
  ! metadata to hdf5 files
  
  ! Therefore, when the grand unified hdf4/hdf5 interfaces are
  ! implemented this will probably need to take an added arg:
  ! the logical addingMetadata

  function mls_sfend(sdid, hdfVersion, addingMetadata)

    ! Arguments

    integer, intent(IN)       :: sdid  
    integer :: mls_sfend            
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: addingmetadata
    ! Local variables
    logical                       :: myaddingmetadata
    integer                       :: myhdfVersion
    integer, external :: PGS_MET_SFend
    logical, parameter :: DEBUG = .false.

    ! begin
   myaddingmetadata = .false.
   if ( present (addingmetadata) ) myaddingmetadata = addingmetadata
   if ( DEBUG ) then
     call output ('Entering mls_sfend with args ', advance='no')
     call output ('sdid: ', advance='no')
     call output (sdid, advance='no')
     call blanks (2)
     call output ('adding meta data?: ', advance='no')
     call output (myaddingmetadata, advance='yes')
   endif
   if ( present(hdfVersion) ) then
     myhdfVersion = hdfVersion
   else
     myhdfVersion = DEFAULT_HDFVERSION
   endif
   if ( myhdfVersion == HDFVERSION_4) then
     mls_sfend = sfend (sdid)
     return
   elseif ( myhdfVersion /= HDFVERSION_5) then
     mls_sfend = WRONGHDFVERSION
     return
   endif
!    mls_sfend = h5fclose_c(sdid)
!   if ( PGS_MET4MLS_SF ) then
   if ( myaddingmetadata ) then
     mls_sfend = PGS_MET_SFend(sdid)
   else
     call h5fclose_f(sdid, mls_sfend)
   endif
!    mls_sfend = 0
   if ( DEBUG ) then
     call output ('Returning from mls_sfend a status: ', advance='no')
     call output (mls_sfend, advance='yes')
   endif

  end function mls_sfend

  ! ---------------------------------------------  mls_hdf_version  -----

  ! This function returns the hdf version of the file:
  ! hdf version         returned value   integer
  !    hdf4                 hdf4           4
  !    hdf5                 hdf5           5
  !  unknown                ????      ERRORINH5FFUNCTION

  function mls_hdf_version(FileName, preferred_version, AccessType) &
   & result (hdf_version)

   ! Arguments

    character (len=*), intent(in)  :: FILENAME                                  
   ! character (len=4)             :: hdf_version                               
    integer(i4), optional, intent(in)  :: PREFERRED_VERSION                     
    integer(i4), optional, intent(in)  :: ACCESSTYPE   ! one of the hdf4 types  
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
    is_hdf5 = (DEFAULT_HDFVERSION == HDFVERSION_5) ! was ... == 5
    call h5fis_hdf5_f(trim(FileName), is_hdf5, returnStatus)
    if ( returnStatus /= 0 ) then
!      hdf_version = '????'
      hdf_version = ERRORINH5FFUNCTION
    elseif ( is_hdf5 ) then
!      hdf_version = 'hdf5'
      hdf_version = HDFVERSION_5       ! 5
    else
!      hdf_version = 'hdf4'
      hdf_version = HDFVERSION_4       ! 4
    endif

   ! hdf_version ==  myPreferred_Version ?
   if ( myPreferred_Version /= WILDCARDHDFVERSION &
     & .and. &
     & hdf_version /= myPreferred_Version ) &
     & hdf_version = WRONGHDFVERSION

  end function mls_hdf_version
!-----------------------------------------------

  subroutine mls_openFile(filename, access, file_id, hdfVersion)
!
! External Variables
!
    character(len=*), intent(in) :: filename, access
    integer, intent(in) :: hdfVersion
    integer(i4), intent(out) :: file_id
!
! Internal Variables
!
    integer :: error
    logical :: is_hdf5

    error = 0
   select case (hdfVersion)
    case (HDFVERSION_4) 
     if (lowercase(trim(access)) == 'create') then
       file_id = mls_sfstart(trim(filename),DFACC_CREATE,hdfVersion)
     elseif (lowercase(trim(access)) == 'update') then
       file_id = mls_sfstart(trim(filename),DFACC_RDWR,hdfVersion)
     elseif (lowercase(trim(access)) == 'readonly') then
       file_id = mls_sfstart(trim(filename),DFACC_RDONLY,hdfVersion)
     endif
     if (file_id .eq. -1) error = file_id 
    case (HDFVERSION_5)
      if (trim(access) == 'create') then 
          call h5fcreate_f(trim(filename), H5F_ACC_EXCL_F, file_id, error)
      elseif (trim(access) == 'update') then
          is_hdf5 = (DEFAULT_HDFVERSION == HDFVERSION_5) 
          call h5fis_hdf5_f(trim(filename), is_hdf5, error)
       if (is_hdf5) then 
          call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
       endif
      elseif (trim(access) == 'readonly') then 
          is_hdf5 = (DEFAULT_HDFVERSION == HDFVERSION_5) 
          call h5fis_hdf5_f(trim(filename), is_hdf5, error)
       if (is_hdf5) then 
          call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
       endif
      else           ! same as 'update'
          is_hdf5 = (DEFAULT_HDFVERSION == HDFVERSION_5) 
          call h5fis_hdf5_f(trim(filename), is_hdf5, error)
       if (is_hdf5) then 
          call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
       endif
      endif
    case default
! >      if (lowercase(trim(access)) == 'create') then 
! >        file_id = mls_sfstart(trim(filename),DFACC_CREATE,hdfVersion)
! >      elseif (lowercase(trim(access)) == 'update') then
! >        file_id = mls_sfstart(trim(filename),DFACC_RDWR,hdfVersion)
! >      elseif (lowercase(trim(access)) == 'readonly') then
! >        file_id = mls_sfstart(trim(filename),DFACC_RDONLY,hdfVersion)
! >      endif
! >      if (file_id .eq. -1) error = file_id 
    end select 

   if (error /= 0) then 
      call MLSMessage (MLSMSG_Error, ModuleName, & 
           "Error: cannot "//lowercase(trim(access))//" file: "//trim(filename))
   endif
 
  end subroutine mls_openFile

!-----------------------------------------------
  subroutine mls_closeFile(file_id,hdfVersion)
!
! External Variables
!
    integer(i4), intent(in) :: file_id
    integer, intent(in) :: hdfVersion
!
! Internal Variables
!
    integer :: error

    select case (hdfVersion)
    case (HDFVERSION_4) 
       error = sfend(file_id)
    case (HDFVERSION_5) 
       call h5fclose_f(file_id, error)
    case default
       error = sfend(file_id)       
    end select 

    if (error /= 0) then 
       call MLSMessage (MLSMSG_Error, ModuleName, "Error closing file.")
    endif

  end subroutine mls_closeFile
!-----------------------------------------------

!====================
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSFiles
!====================

!
! $Log$
! Revision 2.45  2002/12/06 23:38:57  pwagner
! mls_sfstart improved with addingMetaData optional arg
!
! Revision 2.44  2002/12/05 19:44:24  pwagner
! Moved MLSFile_T from MLSFiles to MLSCommon
!
! Revision 2.43  2002/12/04 01:16:25  pwagner
! Open_ Close_ Dump_ Deallocate_ on new MLSFile_T added
!
! Revision 2.42  2002/12/02 23:37:27  pwagner
! First halting steps toward reorg of how mls treats files--will be done via MLSFile_T
!
! Revision 2.41  2002/11/07 21:22:40  jdone
! HDF4/HDF5 capabilities integrated
!
! Revision 2.40  2002/10/29 01:02:26  pwagner
! Reverted DEFAULT_HDFVERSION to 4; added api
!
! Revision 2.39  2002/10/08 23:46:03  pwagner
! Separate he2he5_fileaccess and hdf2hdf5_fileaccess functions
!
! Revision 2.38  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.37  2002/09/27 23:38:38  pwagner
! Added hg type to mls_gen_close
!
! Revision 2.36  2002/09/26 23:59:39  pwagner
! added 'hg' toolbox_mode to mls_io_gen_openF
!
! Revision 2.35  2002/08/08 22:40:42  jdone
! New routines to Open/Close HDF5 Files
!
! Revision 2.34  2002/07/23 00:06:55  pwagner
! debug set to false in sfstart and sfend
!
! Revision 2.33  2002/07/11 22:21:01  pwagner
! These hdf5-savvy versions transferred from he5lib
!
! Revision 1.13  2002/06/12 17:59:37  livesey
! Brought get_free_lun in from lib version
!
! Revision 1.12  2002/05/28 23:11:23  pwagner
! Changed to comply with hdf5.1.4.3/hdfeos5.1.2
!
! Revision 1.11  2002/03/14 23:31:57  pwagner
! HDFVERSION_4 and 5 now public
!
! Revision 1.10  2002/03/13 18:32:05  pwagner
! No longer dumps core after writing metadata to hdf5 files
!
! Revision 1.9  2002/02/19 23:39:29  pwagner
! Eliminated unwanted match between o3 and hno3
!
! Revision 1.8  2002/01/31 00:36:41  pwagner
! Repaired comment statements
!
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
