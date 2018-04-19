! Copyright 2016, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module psivm

! An alternative to pvm

! Based on exchanging Messages via text files left 
! in uniquely-named subdirectories of /tmp
  use HDF, only: DFACC_RdOnly
  use intrinsic, only: l_ascii
  use io_stuff, only: read_textfile, write_textfile
  use machine, only: create_script, execute, execute_command_line, getenv, getids
  use MLSCommon, only: MLSFile_T
  use MLSFiles, only: initializeMLSFile
  use MLSStringLists, only: GetStringElement, NumStringElements

  implicit none
  private

! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (datatypes)
! psiprocess_t  

!     (subroutines and functions)
! InitPsiVM         Set up the global variables to be used
! InitPsiVMD        We are the demon; set up the subdirectory /tmp/PID
! PsiVMFJoinGroup   Set global variables associated with this group
! === (end of toc) ===

! === (start of api) ===
! int InitPsiVM ( )
! int InitPsiVMD ( )
! PsiVMFJoinGroup( char* Name, int status)
! Notes:
! We will require use of files in /tmp
! === (end of api) ===

  public :: psiprocess_t
  public :: initPSIVM, initPSIVMD, PsiVMFJoinGroup
  
  character(len=1024), parameter :: PSIVMTMP = ' '
  
  type psiprocess_t
    character(len=1024)        :: psiPath  = ' '
    character(len=1024)        :: hostname = ' '
    integer                    :: pid      = 0
  end type

  !
  interface discoverMessage
    module procedure discoverMessage_byName, discoverMessage_byTid
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  integer, public                           :: GID
  integer, public                           :: PID
  integer, public                           :: demonsPID
  character(len=16), public                 :: groupName
  character(len=16), public                 :: HostName
  character(len=128)                        :: psivmroot = ' '
  character(len=16), public                 :: USERNAME
  type(psiprocess_t), pointer, dimension(:) :: psiprocessDB => null()
  character(len=1024)                       :: DBFile

contains ! =====     Public Procedures     =============================

  ! Initpsivm
  ! Set global variables associated with this process
  function initPSIVM() result(status)
    integer         :: status
    call getenv( 'HOSTNAME', HOSTNAME )
    call getenv( 'USERNAME', USERNAME )
    call getids( PID, GID )
    status = 0
    DBFile = trim(psivmroot) // '/' // 'PIDDB.txt'
  end function initPSIVM

  ! Initpsivmd
  ! Set global variables associated with the demon
  function initPSIVMD() result(status)
    integer              :: status
    character(len=512)   :: cmd_out
    ! Local variables
    character(len=32)    :: header
    ! Executable
    call getenv( 'USERNAME', USERNAME )
    call getids( PID, GID )
    DBFile = trim(psivmroot) // '/' // 'PIDDB.txt'
    status = 0
    ! 1st--check whether demon is running already
    call execute( 'ps aux | grep psivmd | grep -v grep', status, cmd_out )
    if ( len_trim(cmd_out) > 0 ) then
      status = 1
      return
    endif
    ! mkdir the psivmroot
    ! write( psivmroot, * ) PID
    ! psivmroot = PSIVMTMP // '/' // adjustl(psivmroot)
    call execute( 'mkdir ' // trim(PSIVMTMP), status )
    ! Create the PID database
    write( header, * ) PID, 1, trim(HOSTNAME)
    call write_textfile( DBFile, header )
  end function initPSIVMD

  ! psivmfgettid
  ! Get master tid associated with this group
  subroutine psivmfgettid( Name, tag, tid )
    character(len=*), intent(in) :: Name
    integer, intent(in)          :: tag
    integer, intent(out)         :: tid
    ! Local variables
    character(len=128), dimension(13) :: lines
    character(len=1024)               :: scriptname
    ! Executable
  end subroutine psivmfgettid

  ! PsiVMFJoinGroup
  ! Set global variables associated with this group
  ! Create directory
  subroutine PsiVMFJoinGroup( Name, status)
    character(len=*), intent(in) :: Name
    integer, intent(out)         :: status
    ! Local variables
    character(len=128), dimension(13) :: lines
    character(len=1024)               :: scriptname
    ! Executable
    call getenv( 'USERNAME', USERNAME )
    call getids( PID, GID )
    ! call getDemonPID( USERNAME, demonPID )
    ! write( demonChars, * ) demonPID
    status = 0
    groupName = Name
    ! Create script to create group (if it doesn't exist aleady)
    lines( 1)   = '#!/bin/sh'
    lines( 2)   = '#creategroup.sh'
    lines( 3)   = 'DEMONROOT=' // trim(PSIVMTMP)
    lines( 4)   = 'if [ ! -d "$DEMONROOT" ]'
    lines( 5)   = 'then'
    lines( 6)   = '  mkdir "$DEMONROOT"'
    lines( 7)   = 'fi'
    lines( 8)   = 'GROUPDDIR=' // trim(PSIVMTMP) // '/'  // trim(groupName)
    lines( 9)   = 'if [ ! -d "$GROUPDDIR" ]'
    lines(10)   = 'then'
    lines(11)   = '  mkdir "$GROUPDDIR"'
    lines(12)   = 'fi'
    lines(13)   = 'exit 0'
    ! Now run this script
    scriptname = trim(PSIVMTMP) // '/' // 'create' // trim(groupName) // '.sh'
    call create_script ( scriptname, lines, thenRun=.true. )
  end subroutine PsiVMFJoinGroup

  ! psivmfMyTid
  ! Get my own tid associated with my username and hostname
  subroutine psivmfMyTid( tid )
    integer, intent(out)         :: tid
    ! Local variables
    character(len=128)           :: cmd_out
    character(len=16)            :: PIDChars
    character(len=1024)          :: scriptname     
    integer                      :: status         
    ! Executable
    call getenv( 'USERNAME', USERNAME )
    call getenv( 'HOSTNAME', HOSTNAME )
    write( PIDChars, * ) PID
    call Execute ( 'grep ' // trim(USERNAME) // ' ' // trim(DBFile) // &
      & ' | grep ' // trim(HOSTNAME) // ' |  grep ' // PIDChars,  status, cmd_out )
    if ( len_trim(cmd_out) < 1 ) then
      ! Ask demon to add us
      call askDemonToAddUs ( PID, USERNAME, HOSTNAME, tid )
    else
      read( cmd_out, * ) tid
    endif
  end subroutine psivmfMyTid

! =====     Private Procedures     =====================================

  !-------------------------------------------  AddpsiprocessToDatabase  -----
  integer function AddpsiprocessToDatabase( DATABASE, ITEM )

    ! This function adds an psiprocess data type to a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where psiprocess is put.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (psiprocess_t), dimension(:), pointer :: DATABASE
    type (psiprocess_t), intent(in) :: ITEM

    ! Local variables
    type (psiprocess_T), dimension(:), pointer :: tempDatabase
    include "addItemToDatabase.f9h" 

    AddpsiprocessToDatabase = newSize
  end function AddpsiprocessToDatabase

  subroutine askDemonToAddUs ( PID, USERNAME, HOSTNAME, tid )
    ! To FileDVB
    integer, intent(in)              :: PID
    character(len=*), intent(in)     :: USERNAME
    character(len=*), intent(in)     :: HOSTNAME
    integer, intent(out)             :: tid
    ! Internal variables
    integer, parameter               :: AddUsTag = 300
    logical                          :: exist
    character(len=32)                :: header
    character(len=1024)              :: messageFile
    type(MLSFile_T)                  :: MLSFile
    character(len=16)                :: PIDChars
    ! We will send a message to the demon
    ! The message will consist of a file named
    ! add.PID.USERNAME.HOSTNAME.txt in PSIVMTMP
    write( PIDChars, * ) PID
    messageFile = trim(PSIVMTMP) // '/' // trim(adjustl(PIDChars)) // '.' // &
      & trim(USERNAME) // '.' // trim(HOSTNAME) // '.txt'
    write( header, * ) PID, 1, trim(HOSTNAME)
    call write_textfile ( messageFile, header )
    messageFile = trim(PSIVMTMP) // '/' // trim(adjustl(PIDChars)) // '.' // &
      & trim(USERNAME) // '.' // trim(HOSTNAME) // '.reply'
    do
      ! See if the demon replied
      call discoverMessage ( messageFile, exist, MLSFile )
      if ( exist ) exit
    enddo
  end subroutine askDemonToAddUs
  
!--------------------------- decodeMessage ------------------------------
  ! Decodes the Message passed in the form of a Message string
  subroutine decodeMessage ( Message, tag, id )
    ! Dummy args
    character(len=*), intent(in)    :: Message
    integer, intent(out)            :: tag
    integer, intent(out)            :: id
    ! Internal variables
    character(len=128)  :: string
    ! Executable
    call GetStringElement ( message, string, &
      & 1, countEmpty=.false., inseparator=' ' )
    call ReadIntFromChars( string, tag )
    call GetStringElement ( message, string, &
      & 2, countEmpty=.false., inseparator=' ' )
    call ReadIntFromChars( string, id )
  end subroutine decodeMessage

!--------------------------- discoverMessage ------------------------------
  ! Returns a handle the first file matching Name
  subroutine discoverMessage_byName ( Name, exist, MLSFile )
    ! Dummy args
    character(len=*), intent(in)     :: Name
    logical, intent(out)             :: exist
    type(MLSFile_T), intent(out)     :: MLSFile
    ! Internal variables
    integer                          :: status
    ! Executable
    if ( len_trim(Name) < 1 ) return
    inquire( file=trim(Name), exist=exist )
    if ( .not. exist ) return
    ! We have a match! Ding! Ding! Ding!
    status = InitializeMLSFile ( MLSFile, name=Name, &
      & content = 'message', type=l_ascii, access=DFACC_RDONLY )
  end subroutine discoverMessage_byName

  ! Returns a handle the first file found in the subdirectory
  ! owned by mytid matching tid and tag (if > 0)
  subroutine discoverMessage_byTID ( mytid, tid, tag, MLSFile )
    ! Dummy args
    character(len=*), intent(in)     :: mytid
    integer, intent(in)              :: tid ! Filter available messages, if > 0
    integer, intent(in)              :: tag ! Filter available messages, if > 0
    type(MLSFile_T), intent(out)     :: MLSFile
    ! Internal variables
    character(len=1024)              :: dirName
    character(len=1024)              :: tmpName
    character(len=4096)              :: fileList
    character(len=4096)              :: fileName
    integer                          :: i
    integer                          :: itsId
    integer                          :: itsTag
    integer                          :: n
    integer                          :: status
    ! Executable
    if ( len_trim(mytid) < 1 ) return
    dirName = PSIVMTMP // '/' // mytid
    tmpName = trim(dirName) // '/temp.txt'
    call execute_command_line ( '/bin/rm -f ' // trim(tmpName), exitstat=status )
    call execute_command_line ( 'ls ' // trim(dirName) &
      & // ' > ' // trim(tmpName), &
      & exitstat=status )
    call read_textfile( tmpName, fileList )
    n = NumStringElements( fileList, countEmpty=.false., &
      & inseparator=' ' )
    do i=1, n
      call GetStringElement ( fileList, fileName, &
        & i, countEmpty=.false., inseparator=' ' )
      if ( len_trim(fileName) < 1 ) cycle
      call read_header( trim(dirName) // '/' // fileName, itsId, itsTag )
      if ( tid > 0 .and. tid /= itsId ) cycle
      if ( tag > 0 .and. tag /= itsTag ) cycle
      ! We have a match! Ding! Ding! Ding!
      status = InitializeMLSFile ( MLSFile, name=trim(dirName) // '/' // fileName, &
        & content = 'message', type=l_ascii, access=DFACC_RDONLY )
      return
    enddo
  end subroutine discoverMessage_byTID

  ! Encodes the Message passed in the form of a Message string
  subroutine EncodeMessage ( tag, id, Message )
    ! Dummy args
    integer, intent(in)              :: tag
    integer, intent(in)              :: id
    character(len=1024), intent(out) :: Message
    ! Internal variables
    character(len=128)  :: string
    ! Executable
    write( Message, * ) tag, id, HOSTNAME
  end subroutine EncodeMessage
  
  ! Get PID of psivm demon with your username
  ! That only works if we and the demon are running on the same machine
  ! 
  subroutine getDemonPID( USERNAME, demonPID )
    character(len=*), intent(in) :: USERNAME
    integer, intent(out)         :: demonPID
    ! Local variables
    integer                      :: status
    character(len=128)           :: cmd_out
    ! Executable
    call execute( &
      & "ps aux | grep pvmd | grep" // trim(USERNAME) // &
      & " | grep -v grep | awk '{print $2}'", &
      & status, cmd_out )
    read( cmd_out, * ) demonPID
  end subroutine getDemonPID

  ! Read the header string, assuming it's tag, id
  subroutine read_header( Name, itsId, itsTag )
    ! Dummy args
    character(len=*), intent(in)    :: Name
    integer, intent(out)            :: itsID
    integer, intent(out)            :: itsTAG
    ! local variables
    integer :: unit
    call get_lun ( unit )
    open( unit, file=name, form = 'formatted' )
    read ( unit, * ) itsTag, itsID
    close( unit )
  end subroutine read_header

  ! Read the lastPID from PID file, assuming it's the first line and
  ! contains PID, tid
  subroutine read_lastPID( Name, lastPID, lastTID )
    ! Dummy args
    character(len=*), intent(in)    :: Name
    integer, intent(out)            :: lastPID
    integer, intent(out)            :: lastTID
    ! local variables
    integer :: unit
    call get_lun ( unit )
    open( unit, file=name, form = 'formatted' )
    read ( unit, * ) lastPID, lastTID
    close( unit )
  end subroutine read_lastPID

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module psivm

! $Log$
! Revision 2.2  2018/04/19 02:00:36  vsnyder
! Compute address for allocate/deallocate tracking.  Remove USE statements for
! unused names.
!
! Revision 2.1  2016/02/26 21:02:31  pwagner
! First commit
!
